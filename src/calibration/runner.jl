_to_group_size(flat_vec::Vector{Float64})::Matrix{Float64} =
    permutedims(reshape(flat_vec, (N_SIZE_CLASSES, N_TAXA)), (2, 1))

argmin_missing(x) = argmin(ifelse.(ismissing.(x), Inf, x))
argmax_missing(x) = argmax(ifelse.(ismissing.(x), -Inf, x))

function _relative_deviation_penalty(calibrated::Float64, default::Float64)::Float64
    denom = default == 0.0 ? 1.0 : abs(default)
    rel_dev = abs(calibrated - default) / denom
    return 2.0 / (1.0 + exp(-rel_dev)) - 1.0
end

"""
    validate_gbr_wide_scalar_mean(linear_ext_scalar, survival_scalar)::Float64

Penalise mean scalar values outside [0.95, 1.05] for linear extension and [-0.05, 0.05]
for background mortality. Returns 0.0 when all scalars are valid.
"""
function validate_gbr_wide_scalar_mean(
    linear_ext_scalar::Matrix{Float64},
    survival_scalar::Matrix{Float64}
)::Float64
    mean_lin_ext_scalar = dropdims(mean(linear_ext_scalar, dims=2), dims=2)
    mean_survival_scalar = dropdims(mean(survival_scalar, dims=2), dims=2)
    lin_ext_outside = (!).(0.95 .<= mean_lin_ext_scalar .<= 1.05)
    survival_outside = (!).(-0.05 .<= mean_survival_scalar .<= 0.05)
    return sum(abs.(mean_lin_ext_scalar .- 1.0)[lin_ext_outside]) +
           sum(abs.(mean_survival_scalar)[survival_outside])
end

"""
    validate_linear_extension_coefficients(linear_ext_vals, linear_ext_coefs)::Float64

Return total excess above size-class bin widths across all locations. Returns 0.0 when valid.
"""
function validate_linear_extension_coefficients(
    linear_ext_vals::Matrix{Float64},
    linear_ext_coefs::Matrix{Float64}
)::Float64
    size_class_bins = ADRIA.bin_widths()
    dist_from_valid = 0.0
    for j in 1:size(linear_ext_coefs, 2)
        tmp = (linear_ext_vals .* linear_ext_coefs[:, j]) .- size_class_bins
        tmp[tmp .< 0] .= 0.0
        dist_from_valid += sum(tmp)
    end
    return dist_from_valid
end

function _obj_func(
    init_values::Vector{Float64},
    dom,
    cfg::CalibConfig,
    observations::LocationDataStore,
    location_classification::AbstractVector{Int64}
)::Float64
    local_dom, scen = setup_run(
        dom, init_values;
        param_names=cfg.coral_param_names,
        growth_accel_names=cfg.growth_accel_names,
        depth_atten_names=cfg.depth_atten_names,
        dist_std_names=cfg.dist_std_names,
        dist_mean_names=cfg.dist_mean_names,
        param_idxs=cfg.param_idxs,
        observations=observations,
        biogroup_ord=cfg.biogroups_ordering,
    )

    scale_factors = get_scale_factors(local_dom, scen)
    corals = ADRIA.to_coral_spec(scen[1, :])
    loc_k_areas = ADRIA.loc_k_area(local_dom)
    loc_areas = ADRIA.loc_area(local_dom)
    linear_ext = _to_group_size(corals.linear_extension)

    validity_violation = validate_linear_extension_coefficients(linear_ext, scale_factors[:, 1, :]) +
                         validate_gbr_wide_scalar_mean(scale_factors[:, 1, :], scale_factors[:, 2, :])
    if validity_violation > 0.0
        return 1e4 + 1e4 * validity_violation
    end

    res = nothing
    try
        res = ADRIA.run_model(local_dom, scen[1, :]; apply_allee_effect=false)
    catch err
        _param_assert_patterns = (
            "!any(recruits_scale_factor .<= 0)",
            "Survival rate should be <= 1",
        )
        err_msg = sprint(showerror, err)
        is_param_assert = err isa AssertionError &&
            any(contains(err_msg, p) for p in _param_assert_patterns)
        is_param_assert || rethrow(err)

        # Match only the 35 per-(group, size class) params: a bare "linear_ext"/"mb_rate"
        # substring also catches the 60 biogroup scale factors, so the broadcast below would
        # throw DimensionMismatch from inside this recovery path.
        comp_params = ADRIA.component_params(local_dom.model, :Coral)
        lin_ext_overage = _relative_deviation_penalty.(
            corals.linear_extension,
            comp_params[sc_fg_param_idxs("linear_extension", comp_params), :val]
        )
        mb_rate_overage = _relative_deviation_penalty.(
            corals.mb_rate,
            comp_params[sc_fg_param_idxs("mb_rate", comp_params), :val]
        )
        return 5e5 + sum(lin_ext_overage) + sum(mb_rate_overage) + sum(scale_factors .> 1.0)
    end

    loc_cover = dropdims(sum(res.raw; dims=(2, 3)); dims=(2, 3)) .* loc_k_areas' ./ loc_areas'
    reef_error_series = reef_error(loc_cover; observations=observations)
    reef_perf = year_weighted_error(reef_error_series, reef_observation_counts(observations))

    taxa_cover = dropdims(sum(res.raw; dims=3); dims=3)
    fg_corr = reef_taxa_error(taxa_cover; observations=observations)

    has_obs = vec(any(.!ismissing.(observations.ltmp_coral_cover), dims=2))
    valid_survey_rows = findall(has_obs)
    valid_domain_rows = observations.ltmp_cover_to_domain[valid_survey_rows]

    each_obs_reef = eachrow(observations.ltmp_coral_cover[valid_survey_rows, :])
    obs_peak_col = argmax_missing.(each_obs_reef)
    obs_trough_col = argmin_missing.(each_obs_reef)
    obs_peaks = Float64.(observations.ltmp_coral_cover[CartesianIndex.(valid_survey_rows, obs_peak_col)])
    obs_troughs = Float64.(observations.ltmp_coral_cover[CartesianIndex.(valid_survey_rows, obs_trough_col)])
    modelled_peaks = loc_cover[CartesianIndex.(obs_peak_col, valid_domain_rows)]
    modelled_troughs = loc_cover[CartesianIndex.(obs_trough_col, valid_domain_rows)]

    return reef_perf +
            mean(abs.(modelled_peaks .- obs_peaks)) * 0.5 +
            mean(abs.(modelled_troughs .- obs_troughs)) * 0.5 +
            fg_corr * 2.0
end

mutable struct _StagnationState
    best_fitness::Float64
    last_improved_step::Int
    reason::String
end

"""
    _stagnation_shutdown_check!(
        oc,
        state::_StagnationState;
        patience::Int
    )::Nothing

Stop the optimiser (via `BlackBoxOptim.shutdown!`) once `patience` steps have passed without a
new best fitness. BlackBoxOptim's own `MaxStepsWithoutProgress` parameter is parsed but never
enforced for the single-objective `TopListArchive` path this project uses (verified against its
source — it's only wired up for the Borg multi-objective archive), so this replaces it.

# Arguments
- `oc` : BlackBoxOptim `OptRunController` passed into the callback
- `state` : mutable stagnation-tracking state, updated in place

# Keyword arguments
- `patience` : number of steps without a new best fitness before shutting down
"""
function _stagnation_shutdown_check!(
    oc,
    state::_StagnationState;
    patience::Int
)::Nothing
    current_best = best_fitness(oc)
    current_step = oc.num_steps

    if current_best < state.best_fitness
        state.best_fitness = current_best
        state.last_improved_step = current_step
        return nothing
    end

    if current_step - state.last_improved_step >= patience
        state.reason = "No fitness improvement in $(patience) steps " *
                        "(last improved at step $(state.last_improved_step))."
        @info "Stopping: $(state.reason)"
        BlackBoxOptim.shutdown!(oc)
    end

    return nothing
end

function _save_results_callback(
    oc,
    dom,
    cfg::CalibConfig,
    observations::LocationDataStore,
    out_dir::String,
    last_save::Ref{Float64};
    time_interv::Real=1800.0,
    step_interv::Int=1000
)::Nothing
    start_time = replace(string(unix2datetime(oc.start_time)), "T" => "_", ":" => "")
    elapsed = round(oc.last_report_time - oc.start_time; digits=2)

    if last_save[] == 0.0
        is_save_point = true
    else
        is_save_point = oc.num_steps % step_interv == 0
        if !is_save_point && (oc.last_report_time - last_save[]) > time_interv
            is_save_point = true
        end
    end

    is_save_point || return nothing

    last_save[] = datetime2unix(now(UTC))
    best_state = best_candidate(oc)
    intermediate_path = joinpath(out_dir, "intermediate_results.dat")
    serialize(intermediate_path, best_state)
    write_params_metadata(intermediate_path, length(best_state))

    interim_res = progress_run(
        best_state, dom,
        cfg.coral_param_names, cfg.growth_accel_names, cfg.depth_atten_names,
        cfg.dist_std_names, cfg.dist_mean_names, cfg.param_idxs, observations,
        cfg.biogroups_ordering
    )
    CairoMakie.save(
        joinpath(out_dir, "taxa_cover", "calib_progress_taxa_cover_$(start_time)_$(elapsed).png"),
        taxa_cover_proportions(interim_res.raw)
    )
    CairoMakie.save(
        joinpath(out_dir, "taxa_pop", "calib_progress_pop_cover_$(start_time)_$(elapsed).png"),
        taxa_population_proportions(interim_res.raw)
    )
    @info "Saved intermediate progress"

    return nothing
end

"""
    run_calibration(
        dom,
        cfg::CalibConfig,
        calib_data::CalibrationData,
        location_classification::AbstractVector{Int64},
        stagnation_patience::Int;
        init_cover_path::String,
        out_dir::String,
        init_guess_path::String="",
        config::Union{CalibrationConfig,Nothing}=nothing,
        rng_seed::Int=isnothing(config) ? 42 : config.rng_seed,
        n_threads::Union{Int,Nothing}=nothing,
        time_interv::Real=1800.0,
        step_interv::Int=1000
    )::Nothing

Run BlackBoxOptim to calibrate coral parameters. Calls `construct_cover!` on `dom` in place.

# Arguments
- `dom` : ADRIA domain
- `cfg` : `CalibConfig` describing the parameter search space
- `calib_data` : `CalibrationData` holding the calibration/test observation split
- `location_classification` : per-location classification vector used to construct initial cover
- `stagnation_patience` : stop once this many steps pass without a new best fitness (fully
  replaces the old fixed-`MaxSteps` budget; see `_stagnation_shutdown_check!`)

# Keyword arguments
- `init_cover_path` : path to the serialised initial-cover vector
- `out_dir` : directory for result and progress files. Final result is written to
  `results.dat`; periodic progress saves go to `intermediate_results.dat`, both under `out_dir`.
  The labelled NetCDF products (`params/calibrated_params.nc`, `params/historic_init_cover.nc`)
  are written here too, via `export_calibration_products`
- `init_guess_path` : filename (relative to `out_dir`) for a warm-start vector; `""` disables
- `config` : optional `CalibrationConfig`; when given, defaults `rng_seed` to `config.rng_seed`
  so it can't silently drift from `config.toml`'s `operation.rng_seed`
- `rng_seed` : RNG seed passed to BlackBoxOptim
- `n_threads` : worker threads (`Threads.nthreads()` by default)
- `time_interv` : minimum seconds between progress saves
- `step_interv` : minimum steps between progress saves
"""
function run_calibration(
    dom,
    cfg::CalibConfig,
    calib_data::CalibrationData,
    location_classification::AbstractVector{Int64},
    stagnation_patience::Int;
    init_cover_path::String,
    out_dir::String,
    init_guess_path::String="",
    config::Union{CalibrationConfig,Nothing}=nothing,
    rng_seed::Int=isnothing(config) ? 42 : config.rng_seed,
    n_threads::Union{Int,Nothing}=nothing,
    time_interv::Real=1800.0,
    step_interv::Int=1000
)::Nothing
    init_cover = deserialize(init_cover_path)
    construct_cover!(dom, init_cover, location_classification)

    best_score_file = isempty(init_guess_path) ? "" : joinpath(out_dir, init_guess_path)
    best_init_state = if !isempty(best_score_file) && isfile(best_score_file)
        state = load_calibrated_params(best_score_file, length(cfg.sample_bounds))
        @assert all(first.(cfg.sample_bounds) .<= state .<= last.(cfg.sample_bounds)) "Initial state is out of bounds"
        state
    else
        nothing
    end

    available_threads = something(n_threads, Threads.nthreads())
    threads_display = available_threads == 1 ? available_threads : available_threads - 1
    @info "Using $(threads_display) threads."
    @info "Using RNG seed $(rng_seed) for this calibration run (config.toml: operation.rng_seed)."

    last_save = Ref(0.0)
    stagnation_state = _StagnationState(Inf, 0, "")
    obj = (init_values) -> _obj_func(
        init_values, dom, cfg, calib_data.calibration_store, location_classification
    )
    callback = (oc) -> begin
        _save_results_callback(
            oc, dom, cfg, calib_data.calibration_store, out_dir, last_save;
            time_interv=time_interv, step_interv=step_interv
        )
        _stagnation_shutdown_check!(oc, stagnation_state; patience=stagnation_patience)
    end

    mkpath(joinpath(out_dir, "taxa_cover"))
    mkpath(joinpath(out_dir, "taxa_pop"))

    opts = (
        SearchRange=cfg.sample_bounds,
        MaxSteps=0, # disabled — stagnation_patience (via _stagnation_shutdown_check!) fully
                    # replaces the fixed step budget. 0 must be explicit: BlackBoxOptim's own
                    # default is 10000, so omitting this would silently cap every run at 10000
                    # steps.
        NThreads=available_threads - 1,
        CallbackFunction=callback,
        CallbackInterval=0,
        RngSeed=rng_seed,
        RandomizeRngSeed=false,
    )

    res = if isnothing(best_init_state)
        bboptimize(obj; opts...)
    else
        @info "Using initial guess."
        bboptimize(obj, best_init_state; opts...)
    end

    stop_reason_str = isempty(stagnation_state.reason) ?
        BlackBoxOptim.stop_reason(res) : stagnation_state.reason
    @info "Best fitness: $(best_fitness(res))"
    @info "Stop reason: $(stop_reason_str)"
    final_path = joinpath(out_dir, "results.dat")
    best_params = best_candidate(res)
    serialize(final_path, best_params)
    write_params_metadata(final_path, length(best_params))

    export_calibration_products(
        dom, init_cover, location_classification, best_params, cfg.param_idxs,
        cfg.coral_param_names, cfg.growth_accel_names, cfg.depth_atten_names,
        calib_data.combined_store, cfg.biogroups_ordering;
        out_dir=out_dir
    )

    return nothing
end
