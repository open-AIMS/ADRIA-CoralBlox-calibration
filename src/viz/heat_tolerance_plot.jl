"""
    mean_heat_tolerance_by_group_size(dhw_log, cover_log, init_cover, init_tol, k_area)

Cover-weighted mean heat tolerance (DHW-weeks) across locations, as a
(timesteps, groups * size classes) matrix. `dhw_log` and `cover_log` are the
`coral_dhw_log` / `coral_cover_log` outputs of `ADRIA.run_model` (timesteps x
groups*sizes x locations); the model only writes them from the 2nd timestep, so the
first row is filled from `init_cover` and `init_tol` (both groups*sizes x locations).
Combinations with no cover at any location are NaN.
"""
function mean_heat_tolerance_by_group_size(
    dhw_log, cover_log, init_cover, init_tol, k_area
)::Matrix{Float64}
    tol = Float64.(Array(dhw_log))
    cover = Float64.(Array(cover_log))
    tol[1, :, :] .= init_tol
    cover[1, :, :] .= init_cover

    weights = cover .* reshape(k_area, 1, 1, :)
    w_sum = dropdims(sum(weights; dims=3); dims=3)
    mean_tol = dropdims(sum(tol .* weights; dims=3); dims=3) ./ w_sum
    mean_tol[w_sum .== 0.0] .= NaN
    return mean_tol
end

"""
    plot_heat_tolerance_timeseries(mean_tol; years, n_sizes, exclude_last=1)::Figure

One panel per functional group (rows) and size class (columns) showing the mean heat
tolerance time series. `mean_tol` is the output of `mean_heat_tolerance_by_group_size`.
The last `exclude_last` size classes are not drawn (class 7 is the terminal "end of line"
class).
"""
function plot_heat_tolerance_timeseries(
    mean_tol::AbstractMatrix; years, n_sizes::Int, exclude_last::Int=1, fig_size=(1800, 1100)
)::Figure
    group_names = ADRIA.human_readable_name(ADRIA.functional_group_names(); title_case=true)
    n_groups = size(mean_tol, 2) ÷ n_sizes
    group_colors = [COLORS[Symbol(n)] for n in ADRIA.functional_group_names()]
    n_shown = n_sizes - exclude_last
    ylims = extrema(filter(!isnan, vec(mean_tol)))
    ylims = (ylims[1] - 0.5, ylims[2] + 0.5)

    f = Figure(; size=fig_size)
    for g in 1:n_groups, s in 1:n_shown
        ax = Axis(
            f[g, s];
            title=g == 1 ? "Size class $s" : "",
            ylabel=s == 1 ? "$(group_names[g])\nHeat tolerance (DHW-weeks)" : "",
            xlabel=g == n_groups ? "Year" : "",
            limits=(nothing, ylims),
            xticklabelsvisible=g == n_groups,
            yticklabelsvisible=s == 1
        )
        lines!(ax, years, mean_tol[:, (g - 1) * n_sizes + s]; color=group_colors[g], linewidth=2.5)
    end
    Label(f[0, 1:n_shown], "Mean heat tolerance by functional group and size class"; fontsize=22)
    return f
end
