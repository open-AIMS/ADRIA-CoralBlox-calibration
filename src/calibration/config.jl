struct CalibConfig
    sample_bounds::Vector{Tuple{Float64,Float64}}
    biogroups_ordering::Vector{Int64}
    param_idxs::Vector{Int64}
    coral_param_names::Vector{Symbol}
    growth_accel_names::Vector{String}
    depth_atten_names::Vector{String}
    dist_std_names::Vector{String}
    dist_mean_names::Vector{String}
end

function CalibConfig(dom)::CalibConfig
    # Snapshot of `ADRIA.linear_extensions()` (m/yr) on 2026-10-26, the values the current
    # size class widths were derived from. Hard-coded so the search box does not drift with
    # ADRIA's defaults, which are themselves a calibration output.
    linear_extension_mean = [
        0.0266831 0.0456663   0.055042   0.0699477 0.0631659 0.0726078 0.0
        0.0267012 0.0314565   0.0305842  0.0312019 0.0343453 0.0392649 0.0
        0.0192487 0.0183389   0.016142   0.0181973 0.0178266 0.0160711 0.0
        0.0127533 0.00815496  0.00798983 0.0100552 0.0143346 0.0153482 0.0
        0.0127479 0.00800844  0.00814327 0.0100825 0.0146049 0.0154561 0.0
    ]
    # Size classes 2 onwards cannot be skipped within a timestep, so their upper bound is
    # also limited to the base linear extension that keeps the worst-case effective rate
    # (base x scale x (1 + growth acceleration height)) within the class width. Class 1 and
    # the terminal class have no such limit.
    worst_case_growth = LIN_EXT_SCALE_UB * (1.0 + last(GROWTH_ACCEL_HEIGHT_BOUNDS))
    growing_classes = 2:(N_SIZE_CLASSES - 1)
    lin_ext_ub_matrix = linear_extension_mean .* LIN_EXT_UB_FACTOR
    lin_ext_ub_matrix[:, growing_classes] .= min.(
        lin_ext_ub_matrix[:, growing_classes],
        ADRIA.bin_widths().data[:, growing_classes] ./ worst_case_growth
    )

    lin_ext_lb = flatten_group_size(linear_extension_mean .* collect(LIN_EXT_LB_FACTORS))
    lin_ext_ub = flatten_group_size(lin_ext_ub_matrix)

    mb_rate_mean = [
        0.312661  0.194444  0.211039   0.192857   0.157895   0.142857   0.142857
        0.223847  0.130748  0.0915385  0.123348   0.110294   0.110294   0.110294
        0.218824  0.128571  0.078534   0.0833333  0.0833333  0.0833333  0.0833333
        0.238342  0.0799508 0.0446043  0.026387   0.0135135  0.016      0.0272109
        0.282609  0.0799508 0.0446043  0.026387   0.0135135  0.016      0.0272109
    ]
    mb_rate_lb = flatten_group_size(clamp.(mb_rate_mean .* MB_RATE_LB_FACTOR, 0.0, 1.0))
    mb_rate_ub = flatten_group_size(clamp.(mb_rate_mean .* MB_RATE_UB_FACTOR, 0.0, 1.0))

    all_coral_params = ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.Coral)

    # `dist_mean` and `dist_std` are deliberately excluded from the coral block: both are
    # calibrated per functional group (5 values each) rather than per group and size class
    # (35), so neither maps 1:1 onto ADRIA field names like the rest of the block. They get
    # their own appended blocks below, expanded back to all 35 ADRIA columns in `setup_run`.
    lin_ext_idx, mbrate_idx =
        extract_param_group_idx.([all_coral_params], target_param_names()[1:2])

    # Group-major, size-fastest - the order `repeat(...; inner=N_SIZE_CLASSES)` expands into
    _fieldnames(needle) = string.(
        all_coral_params[sc_fg_param_idxs(needle, all_coral_params), :fieldname]
    )
    dist_mean_names = _fieldnames("dist_mean")
    dist_std_names = _fieldnames("dist_std")

    coral_param_idx = vcat(lin_ext_idx, mbrate_idx)
    coral_params = all_coral_params[sort(coral_param_idx), :]

    lin_ext_idx = sc_fg_param_idxs("linear_extension", coral_params)
    mbrate_idx = sc_fg_param_idxs("mb_rate", coral_params)

    mbrate_scale_idx = findall(occursin.(Ref("mb_rate_scale"), string.(coral_params.fieldname)))
    lin_ext_scale_idx = findall(occursin.(Ref("linear_extension_scale"), string.(coral_params.fieldname)))

    sample_bounds = collect(zip(coral_params.lower_bound, coral_params.upper_bound))
    set_bounds!(sample_bounds, lin_ext_idx, lin_ext_lb, lin_ext_ub)
    set_bounds!(sample_bounds, mbrate_idx, mb_rate_lb, mb_rate_ub)

    n_mbrate_scale = length(mbrate_scale_idx)
    set_bounds!(sample_bounds, mbrate_scale_idx,
        fill(MB_RATE_SCALE_LB, n_mbrate_scale), fill(MB_RATE_SCALE_UB, n_mbrate_scale))

    n_lin_ext_scale = length(lin_ext_scale_idx)
    set_bounds!(sample_bounds, lin_ext_scale_idx,
        fill(LIN_EXT_SCALE_LB, n_lin_ext_scale), fill(LIN_EXT_SCALE_UB, n_lin_ext_scale))

    coral_start_idx = 1
    coral_end_idx = length(sample_bounds)

    biogroups_ordering = sort(unique(dom.loc_data.CB_CALIB_GROUPS))
    n_biogroups = length(biogroups_ordering)

    growth_acc_start_idx = coral_end_idx + 1
    biogroup_accel_bounds = Matrix{Tuple{Float64,Float64}}(undef, n_biogroups, N_GROWTH_ACCEL_PARAMS)
    biogroup_accel_bounds[:, 1] .= [GROWTH_ACCEL_STEEPNESS_BOUNDS]
    biogroup_accel_bounds[:, 2] .= [GROWTH_ACCEL_HEIGHT_BOUNDS]
    biogroup_accel_bounds[:, 3] .= [GROWTH_ACCEL_MIDPOINT_BOUNDS]
    append!(sample_bounds, accel_params_array_to_vec(biogroup_accel_bounds))
    growth_acc_end_idx = length(sample_bounds)

    sc_dist_start_idx = growth_acc_end_idx + 1
    append!(sample_bounds, fill((SC_DIST_LB, SC_DIST_UB), n_biogroups))
    sc_dist_end_idx = length(sample_bounds)

    # Depth attenuation of surface DHW: two GBR-wide scalars, appended last so the existing
    # blocks keep their indices. Order must match DEPTH_ATTEN_PARAM_NAMES.
    depth_atten_start_idx = sc_dist_end_idx + 1
    append!(sample_bounds, [EFF_DHW_BASE_BOUNDS, EFF_DHW_MIX_BOUNDS])
    depth_atten_end_idx = length(sample_bounds)

    # DHW tolerance mean and std: one value per functional group each. Both ADRIA defaults
    # are already constant within a group, so taking every N_SIZE_CLASSES-th element recovers
    # the 5 distinct group values.
    #
    # Anchor `dist_mean` to the :legacy literature values, not the spec default: the :calib
    # defaults are themselves a calibration output, so anchoring to them would let the search
    # box drift with the ADRIA pin and re-centre on each round's own answer. `dist_std` has no
    # :legacy variant and needs none - its default is literature-derived.
    dist_std_anchor = ADRIA.dist_std(; n_sizes=N_SIZE_CLASSES)[1:N_SIZE_CLASSES:end]
    dist_mean_anchor =
        ADRIA.dist_mean(; version=:legacy, n_sizes=N_SIZE_CLASSES)[1:N_SIZE_CLASSES:end]
    @assert length(dist_std_anchor) == N_DIST_STD_PARAMS &&
        length(dist_mean_anchor) == N_DIST_MEAN_PARAMS (
        "Expected $(N_DIST_MEAN_PARAMS) group-level dist_mean/dist_std anchors, got " *
        "$(length(dist_mean_anchor))/$(length(dist_std_anchor))"
    )

    dist_std_start_idx = depth_atten_end_idx + 1
    append!(sample_bounds, collect(zip(
        dist_std_anchor .* DHW_TOL_STD_LB_FACTOR,
        dist_std_anchor .* DHW_TOL_STD_UB_FACTOR
    )))
    dist_std_end_idx = length(sample_bounds)

    dist_mean_start_idx = dist_std_end_idx + 1
    append!(sample_bounds, collect(zip(
        dist_mean_anchor .* DHW_TOL_MEAN_LB_FACTOR,
        dist_mean_anchor .* DHW_TOL_MEAN_UB_FACTOR
    )))
    dist_mean_end_idx = length(sample_bounds)

    param_idxs = [
        coral_start_idx, coral_end_idx,
        growth_acc_start_idx, growth_acc_end_idx,
        sc_dist_start_idx, sc_dist_end_idx,
        depth_atten_start_idx, depth_atten_end_idx,
        dist_std_start_idx, dist_std_end_idx,
        dist_mean_start_idx, dist_mean_end_idx
    ]

    coral_param_names = coral_params.fieldname
    growth_accel_names = accel_params_array_to_vec(
        generate_growth_accel_names(collect(1:n_biogroups))
    )

    # Fail loudly if ADRIA's DepthAttenuation factor names drift from the ones this package
    # writes into the parameter vector and the exported NetCDF - a silent mismatch would
    # leave the calibrated values unused, with the ADRIA defaults quietly standing in.
    adria_depth_names = string.(
        ADRIA.component_params(ADRIA.model_spec(dom), ADRIA.DepthAttenuation).fieldname
    )
    @assert sort(adria_depth_names) == sort(DEPTH_ATTEN_PARAM_NAMES) (
        "ADRIA's DepthAttenuation factors $(adria_depth_names) do not match " *
        "DEPTH_ATTEN_PARAM_NAMES $(DEPTH_ATTEN_PARAM_NAMES)"
    )

    return CalibConfig(
        sample_bounds, biogroups_ordering, param_idxs, coral_param_names,
        growth_accel_names, copy(DEPTH_ATTEN_PARAM_NAMES), dist_std_names, dist_mean_names
    )
end
