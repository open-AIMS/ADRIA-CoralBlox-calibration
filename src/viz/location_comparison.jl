"""
    plot_location_comparison(dom, raw_data::Array{Float64,4}, reef_id::String, cyc_scens, dhw_scens, disturbances; observations::LocationDataStore)::Figure
    plot_location_comparison(dom, raw_data::Array{Float64,4}, ltmp_loc_index::Int64, cyc_scens, dhw_scens, disturbances; observations::LocationDataStore)::Figure

Plot LTMP Manta Tow Coral Cover against the modelled cover output given the LTMP location
index.
The returned figure has four subplots by default: model versus observed series
(`model_vs_obs`), model disturbances (`model_dist`), ltmp disturbances (`ltmp_dist`) and
benthic composition (`benthic`). These four keywords are used in some `opts` to show/hide
these plots and to set axis and figure opts. See `filter_opts` dosctring for more info.

# Arguments
- `opts` :
    - `:short_title` : If true, display only RMSE and SCC. Otherwise shows RMSE, μ bnch,
    PCC, SCC, MAE and BIAS
    - `:show_ltmp_dist` : defaults to true.
    - `:show_model_dist` : defaults to true.
    - `:show_benthic_composition` : defaults to true.
"""
function plot_location_comparison(
    raw_data::Array{Float64,4},
    reef_id::String,
    dhw_scens::YAXArray{T},
    cyc_scens::YAXArray{Float64,3},
    disturbances;
    dom,
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    observations::LocationDataStore,
)::Figure where {T<:Real}
    ltmp_loc_index = findfirst(observations.ltmp_unique_ids .== reef_id)
    return plot_location_comparison(
        raw_data,
        ltmp_loc_index,
        dhw_scens,
        cyc_scens,
        disturbances;
        dom=dom,
        opts=opts,
        fig_opts=fig_opts,
        observations=observations,
    )
end
function plot_location_comparison!(
    grid::Union{GridPosition,GridLayout},
    raw_data::Array{Float64,4},
    reef_id::String,
    dhw_scens::YAXArray{T},
    cyc_scens::YAXArray{Float64,3},
    disturbances;
    dom,
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    observations::LocationDataStore,
)::Nothing where {T<:Real}
    ltmp_loc_index = findfirst(observations.ltmp_unique_ids .== reef_id)
    plot_location_comparison!(
        grid,
        raw_data,
        ltmp_loc_index,
        dhw_scens,
        cyc_scens,
        disturbances;
        dom=dom,
        opts=opts,
        axis_opts=axis_opts,
        fig_opts=fig_opts,
        observations=observations,
    )
    return nothing
end
function plot_location_comparison(
    raw_data::Array{Float64,4},
    ltmp_loc_idx::Int64,
    dhw_scens::YAXArray{T},
    cyc_scens::YAXArray{Float64,3},
    disturbances;
    dom,
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    observations::LocationDataStore,
)::Figure where {T<:Real}
    fig = Figure(; fig_opts...)
    g = fig[1, 1] = GridLayout()
    plot_location_comparison!(
        g,
        raw_data,
        ltmp_loc_idx,
        dhw_scens,
        cyc_scens,
        disturbances;
        dom=dom,
        opts=opts,
        axis_opts=axis_opts,
        fig_opts=fig_opts,
        observations=observations,
    )
    return fig
end
function plot_location_comparison!(
    grid::Union{GridPosition,GridLayout},
    raw_data::Array{Float64,4},
    ltmp_loc_idx::Int64,
    dhw_scens::YAXArray{T},
    cyc_scens::YAXArray{Float64,3},
    disturbances;
    dom,
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    fig_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    observations::LocationDataStore,
)::Nothing where {T<:Real}
    base_title = _location_err_title(raw_data, ltmp_loc_idx; dom=dom, opts, observations=observations)
    title_text = pop!(fig_opts, :title, base_title)
    show_legends = get(opts, :show_legends, true)
    current_row = 1

    subplot_keys = ["model_vs_obs", "model_dist", "ltmp_dist", "benthic"]

    # Plot LTMP and CoralBlox covers
    _except_patterns = setdiff(subplot_keys, ["model_vs_obs"])
    plot_modelled_v_ltmp!(
        grid, current_row, raw_data, ltmp_loc_idx;
        dom=dom,
        observations=observations,
        axis_opts=filter_opts(axis_opts, "model_vs_obs"; except_patterns=_except_patterns),
        opts=filter_opts(opts, "model_vs_obs"; except_patterns=_except_patterns)
    )
    show_legends && legend_modelled_v_ltmp!(grid[current_row, 2])
    current_row += 1

    # These are RME_UNIQUE_ID. Need to select locations from LocationDataStore
    reef_id = observations.ltmp_unique_ids[ltmp_loc_idx]

    if get(opts, :show_model_dist, true)
        # Need to use RME_GBRMPA_ID to select locations from dhw_scens and cyc_scens
        spatial_data = observations.domain_gpkg[:, [:RME_UNIQUE_ID, :RME_GBRMPA_ID]]
        reef_gbrmpa_id = spatial_data[spatial_data.RME_UNIQUE_ID.==reef_id, :RME_GBRMPA_ID][1]

        loc_dhw_scens = dhw_scens[locations=At(reef_gbrmpa_id)]
        loc_cyc_scens = cyc_scens[locs=At(reef_gbrmpa_id)]

        _except_patterns = setdiff(subplot_keys, ["model_dist"])
        plot_coralblox_disturbances!(
            grid[current_row, 1], loc_dhw_scens, loc_cyc_scens;
            axis_opts=filter_opts(axis_opts, "model_dist"; except_patterns=_except_patterns),
            opts=filter_opts(opts, "model_dist"; except_patterns=_except_patterns)
        )
        show_legends && legend_coralblox_disturbances!(grid[current_row, 2])

        current_row += 1
    end

    # Plot CoralBlox Observed DHW and Cyclone Categories
    if get(opts, :show_ltmp_dist, true)
        if reef_id ∈ disturbances.locations
            loc_ltmp_disturbances = disturbances[locations=At(reef_id)]
            _except_patterns = setdiff(subplot_keys, ["ltmp_dist"])
            plot_ltmp_disturbances!(
                grid[current_row, 1], loc_ltmp_disturbances;
                axis_opts=filter_opts(axis_opts, "ltmp_dist"; except_patterns=_except_patterns)
            )

            show_legends && legend_ltmp_disturbances!(grid, current_row, loc_ltmp_disturbances)
            current_row += 1
        else
            @warn("Reef $reef_id / $ltmp_loc_idx not found in LTMP disturbances DataFrame")
        end
    end

    if get(opts, :show_benthic_composition, true)
        cb_loc_id = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
        cover = raw_data[:, :, :, cb_loc_id]
        _except_patterns = setdiff(subplot_keys, ["benthic"])
        plot_taxa_props!(
            grid[current_row, 1], cover;
            axis_opts=filter_opts(axis_opts, "benthic"; except_patterns=_except_patterns),
            opts=filter_opts(opts, "benthic"; except_patterns=_except_patterns)
        )
        show_legends && legend_taxa_props!(grid[current_row, 2])
    end

    titlesize = pop!(fig_opts, :titlesize, 20)
    titlehalign = pop!(fig_opts, :titlehalign, :center)
    titlevalign = pop!(fig_opts, :titlevalign, :center)
    titlefont = pop!(fig_opts, :titlefont, :regular)
    titlejustification = pop!(fig_opts, :titlejustification, :left)
    Label(grid[0, :], title_text, fontsize=titlesize, tellwidth=false, halign=titlehalign,
        valign=titlevalign, font=titlefont, justification=titlejustification, padding=(0, 0, -15, 0)
    )

    return nothing
end

function plot_modelled_v_ltmp!(
    fig::Union{Figure,GridLayout,GridPosition},
    ax_row::Int64,
    raw_data,
    ltmp_loc_idx;
    dom,
    observations::LocationDataStore,
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    _axis_opts = merge(
        Dict(
            :titlesize => FONT_SIZES.title,
            :ylabel => "Relative coral cover",
            BASE_AXIS_OPTS...
        ),
        axis_opts
    )
    ax = Axis(fig[ax_row, 1]; axis_opts...)
    return plot_modelled_v_ltmp!(ax, raw_data, ltmp_loc_idx; dom=dom, observations=observations, opts=opts)
end
function plot_modelled_v_ltmp!(
    ax::Axis,
    raw_data,
    ltmp_loc_idx;
    dom,
    observations::LocationDataStore,
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    domain_loc_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    loc_k_area = ADRIA.loc_k_area(dom)[domain_loc_idx]
    loc_area = ADRIA.loc_area(dom)[domain_loc_idx]

    loc_cover = dropdims(
        sum(raw_data[:, :, :, domain_loc_idx], dims=(2, 3)),
        dims=(2, 3)) .* loc_k_area ./ loc_area

    obs_loc_data = observations.ltmp_coral_cover[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))
    obs_tf = (START_YEAR:END_YEAR)[not_missing_obs]

    markersize = pop!(opts, :markersize, 15)
    obs = scatter!(
        ax, obs_tf, obs_loc_data[not_missing_obs];
        color=(COLORS[:model_vs_obs_color_obs], 0.8), markersize=markersize
    )

    linewidth = pop!(opts, :linewidth, 3)
    sim = lines!(
        ax, 2008:2022, loc_cover;
        color=(COLORS[:model_vs_obs_color_model], 0.9), linewidth=linewidth
    )
    return obs, sim
end

function legend_modelled_v_ltmp!(
    f::Union{Figure,GridLayout,GridPosition,GridSubposition};
    legend_opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    obs_el = MarkerElement(color=COLORS[:model_vs_obs_color_obs], marker=:circle, markersize=15)
    model_el = LineElement(color=COLORS[:model_vs_obs_color_model], linewidth=3)
    Legend(
        f,
        [obs_el, model_el],
        ["Observed", "Model"],
        "Coral Cover",
        valign=:top;
        legend_opts...
    )
end

function plot_coralblox_disturbances!(
    fig::Union{Figure,GridLayout,GridPosition,GridSubposition},
    loc_dhw_scens,
    loc_cyclone_scens;
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Nothing
    axis_opts = merge(
        Dict(
            BASE_AXIS_OPTS...,
            :ylabel => "DHW",),
        axis_opts
    )
    ax = Axis(fig; axis_opts...)

    plot_dhw_scens!(ax, loc_dhw_scens; opts=opts)
    plot_cyclone_scens!(ax, loc_cyclone_scens; opts=opts)

    return nothing
end

function plot_dhw_scens!(
    ax::Axis, dhw_scens;
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Nothing
    dhw_data::Vector{Float64} = Float64.(collect(dhw_scens))
    timesteps::Vector{Int64} = collect(dhw_scens.timesteps.val.data)
    linewidth = get(opts, :linewidth, 4)
    lines!(
        ax, timesteps, dhw_data;
        color=COLORS[:model_dist_color_dhw], linestyle=:dot, linewidth=linewidth
    )
    return nothing
end

function plot_cyclone_scens!(
    ax::Axis, cyc_scens;
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Nothing
    sum_cyc_scens = dropdims(sum(cyc_scens; dims=2), dims=2)
    any(sum_cyc_scens .> 0) || return nothing
    timesteps_vec = collect(sum_cyc_scens.timesteps.val.data)
    target_years = timesteps_vec[collect(sum_cyc_scens) .> 0]
    linewidth = get(opts, :linewidth, 4)
    vlines!(
        ax, target_years;
        ymin=0, color=COLORS[:model_dist_color_dist], linestyle=:dash, linewidth=linewidth
    )
    return nothing
end

function legend_coralblox_disturbances!(
    f::Union{Figure,GridLayout,GridPosition,GridSubposition};
    legend_opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    dhw_el = LineElement(linewidth=4, color=COLORS[:model_dist_color_dhw], linestyle=:dot)
    dist_el = LineElement(linewidth=4, color=COLORS[:model_dist_color_dist], linestyle=:dash)
    return Legend(
        f, [dhw_el, dist_el], ["DHW", "Cyclone/COTS"], "Disturbances", valign=:top;
        legend_opts...
    )
end

function plot_ltmp_disturbances!(
    fig::Union{Figure,GridLayout,GridPosition,GridSubposition}, loc_disturbances; axis_opts::Dict=BASE_AXIS_OPTS
)::Nothing
    limits = get(axis_opts, :limits, nothing)
    xlimits = !isnothing(limits) ? limits[1] : nothing
    axis_opts[:limits] = (xlimits, (0, 1))

    axis_opts = Dict(
        axis_opts...,
        :yticks => [0, 1],
        :yticklabelsize => 0,
        :ylabel => "Disturbances",
    )

    ax_ltmp = Axis(fig; axis_opts...)


    disturbance_vlines, disturbance_types = plot_ltmp_disturbances!(ax_ltmp, loc_disturbances)

    return nothing
end
function plot_ltmp_disturbances!(ax::Axis, loc_disturbances)
    disturbances_types = collect(loc_disturbances.disturbances)
    n_disturbances_types = length(disturbances_types)
    disturbances_years = collect(loc_disturbances.timesteps)
    plotted_vlines = []
    plotted_types = []
    for (disturbance_type_idx, disturbance_type) in enumerate(disturbances_types)
        disturbance_years_mask = loc_disturbances[disturbances=At(disturbance_type)] .== 1
        if any(disturbance_years_mask)
            vline = vlines!(
                ax,
                disturbances_years[disturbance_years_mask],
                color=disturbance_type_idx,
                linewidth=3,
                colormap=COLORMAPS.ltmp_disturbances,
                colorrange=(1, n_disturbances_types)
            )
            push!(plotted_vlines, vline)
            push!(plotted_types, disturbance_type)
        end
    end

    return plotted_vlines, plotted_types
end

function legend_ltmp_disturbances!(
    f::Figure, row::Int64, loc_disturbances::YAXArray{Int16,2}; col::Int64=2
)
    disturbance_mask = read(dropdims(sum(loc_disturbances, dims=1), dims=1) .> 0)
    colors = ColorSchemes.seaborn_bright6[disturbance_mask]
    labels = collect(loc_disturbances.disturbances.val)[disturbance_mask]
    return legend_ltmp_disturbances!(f, row, colors, labels; col=col)
end
function legend_ltmp_disturbances!(
    f::Figure, row::Int64, disturbances::YAXArray{Int16,3}; col::Int64=2
)
    colors = ColorSchemes.seaborn_bright6
    labels = collect(disturbances.disturbance)
    return legend_ltmp_disturbances!(f, row, colors, labels; col=col)
end
function legend_ltmp_disturbances!(
    f::Figure, row::Int64, colors, labels::Vector{String};
    col::Int64=2
)
    n_colors = length(colors)
    ltmp_dist_ele = [LineElement(color=colors[i], linewidth=2) for i in 1:n_colors]
    return Legend(f[row, col], ltmp_dist_ele, labels, "LTMP Disturbances", valign=:top)
end

"""
    plot_taxa_props(loc::String, cover)::Figure

Plot the modelled coral composition for a given location.

# Example
```julia
# ADRIA single run results
rs = ADRIA.run_model(...)
# Index of location of interest
location_idx = ...

f = plot_taxa_props(<Location Name or ID>, rs.raw[:, :, location_idx])
```
"""
function plot_taxa_props(
    loc::String,
    cover;
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Figure
    f = Figure(; size=(1300, 900))
    _axis_opts = merge(
        Dict(
            :xlabel => "year",
            :ylabel => "Cover Proportion",
            :title => "$(loc): Functional Group Cover Proportions",
            :limits => (nothing, nothing, 0, 1)
        ),
        axis_opts
    )

    plot_taxa_props!(f[1, 1], cover; axis_opts=axis_opts, opts=opts)
    return f
end
function plot_taxa_props!(
    fig::Union{Figure,GridLayout,GridPosition,GridSubposition},
    cover;
    axis_opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Nothing
    axis_opts = merge(Dict(:xlabel => "Year", :ylabel => "Proportional Cover",), axis_opts)
    taxa_axis = Axis(fig; axis_opts...)
    taxa_plot = plot_taxa_props!(taxa_axis, cover; opts=opts)
    return nothing
end
function plot_taxa_props!(
    ax::Axis,
    cover;
    opts::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::Axis
    color = COLORMAPS.taxa
    cover = permutedims(cover, (1, 3, 2))
    cover = dropdims(sum(cover, dims=2), dims=2) ./ dropdims(sum(cover, dims=(2, 3)), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022
    linewidth = pop!(opts, :linewidth, 3)
    series!(
        xs,
        cover;
        color=color,
        labels=String.(ADRIA.functional_group_names()),
        linewidth=linewidth,
    )
    return ax
end

function taxa_props_legend_elements()
    taxa = String.(ADRIA.functional_group_names())
    taxa_colors = COLORS[Symbol.(taxa)]
    n_taxa = length(taxa)
    line_el_attrs = (linewidth=2, linepoints=[Point2f(0, 0.5), Point2f(1.0, 0.5)], colgap=0)
    return [
        LineElement(color=taxa_colors[i], ; line_el_attrs...) for i in 1:n_taxa
    ]
end

function legend_taxa_props!(
    f::Union{Figure,GridPosition,GridLayout};
    legend_opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Nothing
    taxa = String.(ADRIA.functional_group_names())
    taxa_labels = titlecase.(join.(split.(taxa, "_"), " "))
    Legend(f, taxa_props_legend_elements(), taxa_labels, "Functional groups"; legend_opts...)
    return nothing
end

"""
    plot_observation_locs(
        calibration_store::LocationDataStore, test_store::LocationDataStore
    )::Figure

Map showing calibration and test locations with different colors.
"""
function plot_observation_locs(
    calibration_store::LocationDataStore, test_store::LocationDataStore
)::Figure
    domain_gpkg = calibration_store.domain_gpkg
    calib_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈Ref(calibration_store.ltmp_unique_ids), :]
    test_gpkg = domain_gpkg[domain_gpkg.UNIQUE_ID.∈Ref(test_store.ltmp_unique_ids), :]

    f = Figure(size=FIG_SIZE[:map])
    ax = GeoAxis(
        f[1, 1],
        dest="+proj=latlong +datum=WGS84",
        title="Observation Locations",
        titlesize=FONT_SIZES[:title]
    )
    try
        poly!(ax, domain_gpkg.geom, color=:black)
    catch
        poly!(ax, domain_gpkg.geometry, color=:black)
    end

    obs = (Calibration=:red, Test=:blue)

    try
        scatter!(ax, calib_gpkg.X_COORD, calib_gpkg.Y_COORD; markersize=10, color=obs.Calibration, alpha=0.5)
        scatter!(ax, test_gpkg.X_COORD, test_gpkg.Y_COORD; markersize=10, color=obs.Test, alpha=0.5)
    catch
        scatter!(ax, calib_gpkg.LON, calib_gpkg.LAT; markersize=10, color=obs.Calibration, alpha=0.5)
        scatter!(ax, test_gpkg.LON, test_gpkg.LAT; markersize=10, color=obs.Test, alpha=0.5)
    end


    els = [MarkerElement(color=c, marker=:circle, markersize=10) for c in values(obs)]
    labels = collect(string.(keys(obs)))

    Legend(
        f[2, 1],
        els,
        labels,
        valign=:top,
        orientation=:horizontal
    )

    return f
end
