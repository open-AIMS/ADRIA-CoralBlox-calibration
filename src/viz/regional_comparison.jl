BASE_OPTS_REGIONAL_COMPARISON::Dict{Symbol,Any} = Dict(
    :model_line_color => COLORS[:model_line_color],
    :model_band_color => COLORS[:model_band_color],
    :historic_line_color => COLORS[:historic_line_color],
    :historic_band_color => COLORS[:historic_band_color],
    :markersize => 8
)

"""
# Example
```
management_areas = unique(dom.loc_data.management_area_short)
management_area_masks = [dom.loc_data.management_area_short .== area for area in management_areas]

management_area_stats = region_stats(Symbol.(management_area_masks), management_areas)

plot_regional_comparison(management_area_stats)
```
"""
function plot_regional_comparison(
    regional_stats::NamedTuple,
    regional_stats_sorted_keys::Tuple;
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    legend_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    opts::Dict{Symbol,Any} = merge(BASE_OPTS_REGIONAL_COMPARISON, opts)
    f = Figure(; fig_opts...)
    plot_regional_comparison!(
        f, regional_stats, regional_stats_sorted_keys;
        axis_opts=axis_opts, opts=opts
    )

    if get(legend_opts, :show_legend, true)
        legend_regional_comparison!(f[end+1, :]; legend_opts=legend_opts)
    end

    return f
end
function plot_regional_comparison!(
    grid::Union{GridPosition,GridLayout,Figure},
    regional_stats::NamedTuple,
    regional_stats_sorted_keys::Tuple;
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    n_plots = length(regional_stats_sorted_keys)
    invert_positions = get(opts, :invert_positions, true)
    positions = fig_coord.(1:n_plots, Ref(n_plots); invert=invert_positions)

    last_row_idx = maximum(getindex.(positions, 1))
    for (idx, key) in enumerate(regional_stats_sorted_keys)
        n_locations = regional_stats[key].n_locations
        _rmse = trunc(regional_stats[key].rmse, digits=2)
        _mae = trunc(regional_stats[key].mae, digits=2)
        _srcc = trunc(regional_stats[key].srcc, digits=2)

        # position = fig_coord(idx, n_plots; invert=true)
        label_text = "$((string(key))) ($n_locations reefs)"
        group_title = "$label_text\nRMSE: $_rmse | SRCC : $_srcc"
        ax = Axis(
            grid[positions[idx]...];
            title=group_title,
            xlabel="Year",
            ylabel="Relative cover",
            xticklabelsvisible=positions[idx][1] == last_row_idx ? true : false,
            xlabelvisible=positions[idx][1] == last_row_idx ? true : false,
            ylabelvisible=positions[idx][2] == 1 ? true : false,
            yticklabelsvisible=positions[idx][2] == 1 ? true : false,
            limits=(nothing, (0, 1)),
            # limits=axis_limits,
            xticks=get(axis_opts, :xticks, START_YEAR:3:END_YEAR),
            xticklabelrotation=get(axis_opts, :xticklabelrotation, π / 4),
            axis_opts...
        )

        plot_regional_comparison!(
            ax,
            regional_stats[key].model_stats,
            replace(regional_stats[key].historic_stats, missing => NaN),
            opts,
        )

        if get(opts, :showtextlabel, false)
            tlsize = get(opts, :textlabelsize, 16)
            tlbackground = get(opts, :textlabelbackground, :white)
            tlposition = get(opts, :textlabelposition, Point2f(2008, 0.9))
            textlabel!(
                ax,
                tlposition;
                text_align=(:left, :top),
                text=label_text,
                fontsize=tlsize,
                strokewidth=0,
                background_color=tlbackground,
                cornerradius=0.0
            )
        end
    end

    if get(opts, :show_title, true)
        titlesize = get(opts, :titlesize, 18)
        Label(grid[0, :], get(opts, :title, ""); fontsize=titlesize)
    end

    return nothing
end
function plot_regional_comparison!(
    ax::Axis,
    model_stats::Matrix{Float64},
    historic_stats::Matrix{Float64},
    opts::Dict{Symbol,Any};
    xdata=START_YEAR:END_YEAR
)::Nothing
    band!(
        ax, xdata, model_stats[1, :], model_stats[3, :];
        color=(get(opts, :model_band_color, COLORS[:model_band_color]))
    )
    band!(
        ax, xdata, historic_stats[1, :], historic_stats[3, :];
        color=get(opts, :historic_band_color, COLORS[:historic_band_color])
    )

    scatterlines!(
        ax, xdata, collect(model_stats[2, :]);
        color=get(opts, :model_line_color, COLORS[:model_line_color]),
        markersize=get(opts, :markersize, 20)
    )
    scatterlines!(
        ax, xdata, historic_stats[2, :];
        color=get(opts, :historic_line_color, COLORS[:historic_line_color]),
        markersize=get(opts, :markersize, 20)
    )

    return nothing
end

function legend_regional_comparison!(
    fig::Union{GridPosition,GridLayout,Figure};
    legend_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    model_leg_element = [
        LineElement(color=pop!(legend_opts, :model_line_color, COLORS[:model_line_color])),
        PolyElement(color=pop!(legend_opts, :model_band_color, COLORS[:model_band_color]))
    ]
    observed_leg_element = [
        LineElement(color=pop!(legend_opts, :historic_line_color, COLORS[:historic_line_color])),
        PolyElement(color=pop!(legend_opts, :historic_band_color, COLORS[:historic_band_color]))
    ]
    leg_labels = ["Observed", "Model"]

    Legend(
        fig,
        [observed_leg_element, model_leg_element],
        leg_labels;
        orientation=:horizontal,
        legend_opts...
    )
    return nothing
end
function legend_regional_comparison!(
    fig::Union{GridPosition,GridLayout,Figure},
    opts_test::Dict{Symbol,Any},
    opts_calibration::Dict{Symbol,Any})::Nothing
    model_leg_element_test = [
        LineElement(color=opts_test[:model_line_color]),
        PolyElement(color=opts_test[:model_band_color])
    ]
    historic_leg_element_test = [
        LineElement(color=opts_test[:historic_line_color]),
        PolyElement(color=opts_test[:historic_band_color])
    ]
    model_leg_element_calib = [
        LineElement(color=opts_calibration[:model_line_color]),
        PolyElement(color=opts_calibration[:model_band_color])
    ]
    historic_leg_element_calib = [
        LineElement(color=opts_calibration[:historic_line_color]),
        PolyElement(color=opts_calibration[:historic_band_color])
    ]

    leg_labels = [
        "CoralBlox test Reefs",
        "LTMP test Reefs",
        "CoralBlox calibration Reefs",
        "LTMP calibration Reefs"]

    leg_elements = [
        model_leg_element_test,
        historic_leg_element_test,
        model_leg_element_calib,
        historic_leg_element_calib
    ]

    Legend(
        fig,
        leg_elements,
        leg_labels;
        orientation=:horizontal
    )

    return nothing
end

"""
    plot_all_regions(dom, model_results; region_masks=[NORTH_MASK, CENTRAL_MASK, SOUTH_MASK], ref_years=START_YEAR:END_YEAR)::Figure
    plot_all_regions(north_results, north_obs, central_results, central_obs, south_results, south_obs)::Figure

Plot the modelled results against observations/modelled observation results.
"""
function plot_all_regions(
    dom,
    raw_data::Array{Float64,4},
    ltmp_north, ltmp_central, ltmp_south;
    region_masks::Vector{BitVector},
    ref_years=START_YEAR:END_YEAR,
    fig_title="",
    fig_size=(800, 800)
)::Figure
    s_rac = (dropdims(sum(raw_data, dims=(2, 3)), dims=(2, 3)) .*
             ADRIA.loc_k_area(dom)') ./
            ADRIA.loc_area(dom)'

    north_res = ADRIA.DataCube(s_rac[:, region_masks[1]];
        timesteps=ref_years, locations=1:count(region_masks[1]))
    central_res = ADRIA.DataCube(s_rac[:, region_masks[2]];
        timesteps=ref_years, locations=1:count(region_masks[2]))
    south_res = ADRIA.DataCube(s_rac[:, region_masks[3]];
        timesteps=ref_years, locations=1:count(region_masks[3]))

    f = plot_all_regions(
        north_res, ltmp_north, central_res, ltmp_central, south_res, ltmp_south;
        fig_title=fig_title, fig_size=fig_size
    )

    return f
end
function plot_all_regions(
    north_results, north_obs, central_results, central_obs, south_results, south_obs;
    fig_title="", fig_size=(800, 500)
)::Figure
    f = Figure(; size=fig_size)

    plot_region!(
        f[1, 1],
        "North GBR\n$(size(north_results, 2)) test locations",
        north_obs,
        north_results
    )
    plot_region!(
        f[1, 2],
        "Central GBR\n$(size(central_results, 2)) test locations",
        central_obs,
        central_results
    )
    plot_region!(
        f[1, 3],
        "South GBR\n$(size(south_results, 2)) test locations",
        south_obs,
        south_results
    )

    obs_el = [
        PolyElement(color=(:black, 0.4)),
        LineElement(color=:black)
    ]
    sim_el = [
        PolyElement(color=(:red, 0.4)),
        LineElement(color=:red)
    ]

    Legend(
        f[2, :],
        [obs_el, sim_el],
        ["LTMP", "CoralBlox"],
        orientation=:horizontal
    )

    # linkyaxes!(north_ax, central_ax, south_ax)
    # Label(f[0, :], fig_title; fontsize=20)

    # resize_to_layout!(f)
    return f
end

"""
    plot_region(f::Figure, row::Int64, col::Int64, title::String, obs::DataFrame, sim::YAXArray; showlegend::Bool = false, legend_row::Int64 = 2, legend_col::Int64 = 2)::Axis

"""
function plot_region!(
    g::GridPosition,
    title::String,
    obs::DataFrame,
    sim::YAXArray;
)::Nothing
    ax = Axis(
        g,
        title=title,
        xlabel="year",
        ylabel="Relative Absolute Cover",
        limits=(nothing, (0.0, 1.0))
    )

    # mean_agg = dropdims(mean(sim, dims=:locations), dims=:locations)
    confints = ADRIA.analysis.series_confint(read(sim))
    xs::Vector{Float64} = collect(sim.timesteps)

    band!(ax, obs.Year, obs.lower, obs.upper, color=(:black, 0.4))
    band!(ax, xs, confints[:, 1], confints[:, 3], color=(:red, 0.4))

    lines!(ax, obs.Year, obs.response, color=:black)
    lines!(ax, xs, confints[:, 2], color=:red, linewidth=2)

    return nothing
end
