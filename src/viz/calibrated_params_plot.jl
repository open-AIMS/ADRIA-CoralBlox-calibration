const _PARAM_EDGE_TOL = 0.02

function _strip_fg_sc_label(name::AbstractString, suffix::AbstractString)::String
    return replace(replace(name, "_$(suffix)" => ""), "_" => " ")
end

"""
    params_vector_from_dataset(ds, n_params, param_idxs, coral_param_names)::Vector{Float64}

Invert `build_params_dataset`: repack a `calibrated_params.nc` dataset into the flat
sample-vector layout. The initial-cover lambda block is not stored in the NetCDF and stays `NaN`.
"""
function params_vector_from_dataset(
    ds, n_params::Int, param_idxs::Vector{Int64}, coral_param_names
)::Vector{Float64}
    params = fill(NaN, n_params)
    names = string.(coral_param_names)
    fg_sc(v) = Array(permutedims(v, (:functional_group, :size_class)))
    fg_bg(v) = Array(permutedims(v, (:functional_group, :cb_calib_group)))

    les = occursin.("linear_extension_scale", names)
    mbs = occursin.("mb_rate_scale", names)
    le = occursin.("linear_extension", names) .& .!les
    mb = occursin.("mb_rate", names) .& .!mbs

    coral = zeros(length(names))
    coral[le] .= vec(permutedims(fg_sc(ds.linear_extension), (2, 1)))
    coral[mb] .= vec(permutedims(fg_sc(ds.mb_rate), (2, 1)))
    coral[les] .= vec(permutedims(fg_bg(ds.linear_extension_scale), (2, 1)))
    coral[mbs] .= vec(permutedims(fg_bg(ds.mb_rate_scale), (2, 1)))
    params[param_idxs[1]:param_idxs[2]] .= coral

    ga = Array(permutedims(ds.growth_acceleration, (:cb_calib_group, :accel_param)))
    params[param_idxs[3]:param_idxs[4]] .= vec(ga)
    params[param_idxs[7]:param_idxs[8]] .= Array(ds.depth_attenuation)
    params[param_idxs[9]:param_idxs[10]] .= Array(ds.dist_std)
    # `dist_mean` holds the end-of-history tolerance per (group, size class) when
    # `dist_mean_initial` exists; the search-space value is the latter
    dist_mean = :dist_mean_initial in keys(ds.cubes) ? ds.dist_mean_initial : ds.dist_mean
    params[param_idxs[11]:param_idxs[12]] .= Array(dist_mean)
    return params
end

"""
    plot_calibrated_params(params, sample_bounds, param_idxs, coral_param_names,
        growth_accel_names, depth_atten_names, biogroups_ordering; dist_mean_end, fig_size)

Plot every calibrated parameter against its search bounds in a three-column figure. Each
block of related parameters has its own x-axis, since magnitudes differ by orders. Bars span
the search bounds; points are the calibrated values, drawn in red when within
$(round(Int, 100 * _PARAM_EDGE_TOL))% of the bound range of either edge (the search ran into its box).

`params` holds the calibrated initial `dist_mean` (the searched value). If `dist_mean_end`
(`(n_groups, n_sizes)`, the end-of-history tolerance stored in `calibrated_params.nc`) is
given, it is drawn as green diamonds, one per size class, next to the initial value. The
diamond size grows with the size class index while its shade gets lighter, and larger classes
are drawn first so the smaller ones stay in front.
"""
function plot_calibrated_params(
    params::Vector{Float64},
    sample_bounds::Vector{<:Tuple{Float64,Float64}},
    param_idxs::Vector{Int64},
    coral_param_names,
    growth_accel_names::Vector{String},
    depth_atten_names::Vector{String},
    biogroups_ordering::Vector{Int64};
    dist_mean_end::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    fig_size::Tuple{Int,Int}=(1800, 1500)
)::Figure
    fg_names = string.(ADRIA.functional_group_names())
    coral_names = string.(coral_param_names)
    coral_offset = param_idxs[1] - 1

    les_mask = occursin.("linear_extension_scale", coral_names)
    mbs_mask = occursin.("mb_rate_scale", coral_names)
    le_mask = occursin.("linear_extension", coral_names) .& .!les_mask
    mb_mask = occursin.("mb_rate", coral_names) .& .!mbs_mask

    coral_block(mask, suffix) = (
        coral_offset .+ findall(mask),
        [_strip_fg_sc_label(n, suffix) for n in coral_names[mask]]
    )
    scale_block(mask, prefix) = (
        coral_offset .+ findall(mask),
        [replace(replace(n, "$(prefix)_cb_group_" => "grp "), "_" => " ") for n in coral_names[mask]]
    )

    accel_idxs = collect(param_idxs[3]:param_idxs[4])
    accel_labels = [
        replace(replace(n, "growth_acceleration_cb_group_" => "grp "), "_" => " ")
        for n in growth_accel_names
    ]

    blocks = [
        [
            ("Linear extension (m)", coral_block(le_mask, "linear_extension")...),
            ("DHW tolerance mean (● initial, ◆ end, by size class)", collect(param_idxs[11]:param_idxs[12]), fg_names),
            ("DHW tolerance std", collect(param_idxs[9]:param_idxs[10]), fg_names),
            ("Depth attenuation", collect(param_idxs[7]:param_idxs[8]), depth_atten_names),
        ],
        [
            ("Background mortality rate", coral_block(mb_mask, "mb_rate")...),
            ("Growth acceleration", accel_idxs, accel_labels),
        ],
        [
            ("Linear extension scale", scale_block(les_mask, "linear_extension_scale")...),
            ("Mortality rate scale", scale_block(mbs_mask, "mb_rate_scale")...),
        ],
    ]

    fig = Figure(; size=fig_size)
    for (col, col_blocks) in enumerate(blocks)
        gl = fig[1, col] = GridLayout()
        for (row, (title, idxs, labels)) in enumerate(col_blocks)
            ax = Axis(
                gl[row, 1];
                title=title,
                titlealign=:left,
                yticks=(1:length(idxs), labels),
                yreversed=true,
                yticklabelsize=11,
                xticklabelsize=11,
                yticksvisible=false,
                yaxisposition=:left,
                limits=(nothing, (0.3, length(idxs) + 0.7)),
            )
            lb = [sample_bounds[i][1] for i in idxs]
            ub = [sample_bounds[i][2] for i in idxs]
            vals = params[idxs]
            ys = 1:length(idxs)

            rangebars!(ax, collect(ys), lb, ub; direction=:x, color=("#1E88E5", 0.45), linewidth=6)
            at_edge = (vals .- lb .<= _PARAM_EDGE_TOL .* (ub .- lb)) .|
                      (ub .- vals .<= _PARAM_EDGE_TOL .* (ub .- lb))
            scatter!(
                ax, vals, collect(ys);
                color=ifelse.(at_edge, "#D81B60", "#222222"), markersize=7
            )
            if startswith(title, "DHW tolerance mean") && !isnothing(dist_mean_end)
                n_sizes = size(dist_mean_end, 2)
                size_idx = repeat(1:n_sizes; inner=size(dist_mean_end, 1))
                # largest size class first, so the smaller (darker) diamonds sit on top
                draw_order = sortperm(size_idx; rev=true)
                scatter!(
                    ax,
                    vec(Float64.(dist_mean_end))[draw_order],
                    repeat(collect(ys); outer=n_sizes)[draw_order];
                    marker=:diamond, color=size_idx[draw_order], colormap=Reverse(:Greens),
                    colorrange=(1, n_sizes + 3), markersize=5 .+ 1.5 .* size_idx[draw_order],
                    strokewidth=0.5, strokecolor=("#1B5E20", 0.6)
                )
            end
            rowsize!(gl, row, Auto(length(idxs) + 2))
        end
    end

    return fig
end
