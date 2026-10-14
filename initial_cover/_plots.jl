using WGLMakie, GeoMakie, GraphMakie
function plot_class_size_props(
    init_cover,
    class_idx
)::Figure
    class_state = init_cover[(class_idx-1)*11+1:class_idx*11]
    taxa_names = ADRIA.functional_group_names()
    taxa_prop = class_state[7:11]
    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="taxonomy", xticks=(1:5, String.(taxa_names)), ylabel="taxa size lambda", title="class: $(class_idx)")
    barplot!(
        ax,
        1:5,
        taxa_prop
    )
    return f
end

function plot_class_properties(
    init_cover,
    class_idx
)::Figure
    class_state = init_cover[(class_idx-1)*11+1:class_idx*11]
    taxa_names = ADRIA.functional_group_names()
    taxa_prop = class_state[2:6] ./ sum(class_state[2:6])
    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="taxonomy", xticks=(1:5, String.(taxa_names)), ylabel="taxa proportions", title="class: $(class_idx)")
    barplot!(
        ax,
        1:5,
        taxa_prop
    )
    return f
end



function plot_ltmp(ltmp_n, ltmp_c, ltmp_s)::Nothing
    f = Figure(; size=(900, 1600))
    Axis(
        f[1, 1],
        title="North GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )

    band!(ltmp_n.Year, ltmp_n.lower, ltmp_n.upper, color=(:black, 0.4))
    lines!(ltmp_n.Year, ltmp_n.response, color=:black)
    Axis(
        f[2, 1],
        title="Central GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )


    band!(ltmp_c.Year, ltmp_c.lower, ltmp_c.upper, color=(:black, 0.4))
    lines!(ltmp_c.Year, ltmp_c.response, color=:black)
    Axis(
        f[3, 1],
        title="South GBR",
        xlabel="Year",
        ylabel="Relative Absolute Cover"
    )

    band!(ltmp_s.Year, ltmp_s.lower, ltmp_s.upper, color=(:black, 0.4))
    lines!(ltmp_s.Year, ltmp_s.response, color=:black)
    save("LTMP_data.png", f)
    return nothing
end

function location_comparison(
    raw_data,
    ltmp_loc_idx;
    obs_data=raw_ltmp_reef_data,
    obs_idxs=ltmp_reefmod_idxs,
    obs_loc_labels=ltmp_reef_data.RME_UNIQUE_ID,
    loc_k_areas=loc_k_area(dom),
    loc_areas=site_area(dom)
)::Figure
    loc_cover = dropdims(sum(raw_data, dims=2), dims=2) .* loc_k_areas' ./ loc_areas'
    if obs_idxs[ltmp_loc_idx] == -1
        return Figure()
    end

    obs_loc_data = obs_data[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))
    obs_tf = (2008:2022)[not_missing_obs]

    if !any(not_missing_obs)
        return Figure()
    end

    sim_data = loc_cover[:, obs_idxs[ltmp_loc_idx]]
    reef_id = obs_loc_labels[ltmp_loc_idx]

    f = Figure()
    Axis(f[1, 1], xlabel="Year", ylabel="relative total area", title="Location $(reef_id)")
    obs = scatter!(obs_tf, obs_loc_data[not_missing_obs], color=:transparent, strokewidth=2, strokecolor=:black, markersize=15)
    sim = lines!(2008:2022, sim_data, color=:red)
    Legend(
        f[1, 2],
        [obs, sim],
        ["LTMP", "CoralBlox"]
    )
    save("Outputs/loc_plots/loc_$(reef_id).png", f)
    return f
end

function taxa_cover_proportions(raw_data)::Figure
    cover = reshape(raw_data, (15, 7, 5, 3806))
    cover = dropdims(sum(cover, dims=4), dims=4)
    cover ./= sum(cover, dims=(2, 3))
    cover = dropdims(sum(cover, dims=2), dims=2)
    cover = permutedims(cover, (2, 1))
    xs = 2008:2022

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Cover Proportion",
        title="Functional Group Cover Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, cover, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)
    return f
end


function taxa_population_proportions(raw_data)::Figure
    sc_mean_area = reshape(
        permutedims(ADRIA.colony_areas()[2], (2, 1)),
        (1, 35, 1)
    )
    population = raw_data ./ sc_mean_area
    population = reshape(population, (15, 7, 5, 3806))
    population = dropdims(sum(population, dims=4), dims=4)
    population ./= sum(population, dims=(2, 3))
    population = dropdims(sum(population, dims=2), dims=2)
    population = permutedims(population, (2, 1))
    xs = 2008:2022

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title="Functional Group Population Proportions",
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, population, color=:Paired_5, labels=String.(ADRIA.functional_group_names()))
    Legend(f[1, 2], ax, framevisible=false)
    return f
end


"""
    temporal_size_class_proportions(raw_data)

Calculate the percentage of coral population occupied by each size class.
"""
function temporal_size_class_proportions(raw_data)::Figure
    sc_mean_area = reshape(
        permutedims(ADRIA.colony_areas()[2], (2, 1)),
        (1, 35, 1)
    )
    population = raw_data ./ sc_mean_area
    population = reshape(population, (15, 7, 5, 3806))
    population = dropdims(sum(population, dims=4), dims=4)
    population ./= sum(population, dims=2)
    population = permutedims(population, (3, 2, 1))

    fg_names = ADRIA.functional_group_names()
    xs = 2008:2022
    col = :PuBuGn_9

    f = Figure(; size=(1600, 900))
    ax = Axis(
        f[1, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[1]),
        limits=(nothing, nothing, 0, 1)
    )
    sr = series!(xs, population[1, :, :], color=col, labels="Size Class: " .* string.(1:7))
    Axis(
        f[1, 2];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[2]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[2, :, :], color=col)
    Axis(
        f[2, 1];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[3]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[3, :, :], color=col)
    Axis(
        f[2, 2];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[4]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[4, :, :], color=col)
    Axis(
        f[1, 3];
        xlabel="year",
        ylabel="Population Proportion",
        title=String(fg_names[5]),
        limits=(nothing, nothing, 0, 1)
    )
    series!(xs, population[5, :, :], color=col)
    Legend(f[2, 3], ax, framevisible=false)
    colsize!(f.layout, 3, Relative(1 / 3))
    resize_to_layout!(f)
    return f
end


function plot_region(
    f::Figure,
    row::Int64,
    col::Int64,
    title::String,
    obs::DataFrame,
    sim::YAXArray;
    showlegend::Bool=false,
    legend_row::Int64=2,
    legend_col::Int64=2
)::Axis

    ax = Axis(
        f[row, col],
        title=title,
        xlabel="year",
        ylabel="Relative Absolute Cover",
        width=400,
        height=400
    )

    mean_agg = dropdims(mean(sim, dims=:sites), dims=:sites)
    confints = ADRIA.analysis.series_confint(mean_agg.data)

    xs::Vector{Float64} = collect(sim.timesteps)
    ADRIA_series = series!(xs, sim[:, :, 1].data[:, :]'; solid_color=(:red, 0.01), linewidth=1, labels=nothing)
    ADRIA_line = lines!(xs, confints[:, 2], color=:red, linewidth=5)

    xs = obs.Year
    lower::Vector{Float64} = obs.lower
    upper::Vector{Float64} = obs.upper
    response::Vector{Float64} = obs.response

    LTMP_band = band!(xs, lower, upper, color=(:black, 0.4))
    LTMP_line = lines!(xs, response, color=:black)

    if showlegend
        Legend(
            f[legend_row, legend_col],
            [[LTMP_line, LTMP_band], [ADRIA_line, ADRIA_series]],
            ["LTMP", "CoralBlox"]
        )
    end

    return ax
end


function plot_residual(
    f::Figure,
    ax_row::Int64,
    ax_col::Int64,
    ltmp_data::DataFrame,
    ADRIA_data::YAXArray,
    title::String
)::Nothing

    Axis(
        f[ax_row, ax_col],
        title=title,
        xlabel="year",
        ylabel="Residuals"
    )

    ltmp_ind::Int64 = findfirst(x -> x >= start_year, ltmp_data.Year)
    ltmp_end::Int64 = findfirst(x -> x >= end_year, ltmp_data.Year) - 1

    ltmp_xs = ltmp_data.Year[ltmp_ind:ltmp_end]
    ltmp_resp = ltmp_data.response[ltmp_ind:ltmp_end]

    adria_mean = dropdims(mean(
            ADRIA_data[timesteps=At(ltmp_xs)], dims=(:sites, :scenarios)
        ), dims=(:sites, :scenarios))

    residuals::Vector{Float64} = ltmp_resp .- adria_mean
    scatter!(ltmp_xs, residuals, color=:black)

    return nothing
end

function plot_residuals(
    filename::String; fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:size => (1600, 1600))
)::Figure

    f = Figure(; fig_opts...)

    plot_residual(
        f,
        1,
        1,
        ltmp_north,
        north_res,
        "North GBR"
    )

    plot_residual(
        f,
        2,
        1,
        ltmp_central,
        central_res,
        "Central GBR"
    )

    plot_residual(
        f,
        3,
        1,
        ltmp_south,
        south_res,
        "South GBR"
    )

    save(filename, f)

    return f
end

function create_compare_plot(
    filename::String; fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:size => (1600, 1600))
)::Figure

    f = Figure(; fig_opts...)

    plot_region(
        f,
        1,
        1,
        "North GBR",
        ltmp_north,
        north_res
    )
    plot_region(
        f,
        1,
        2,
        "Central GBR",
        ltmp_central,
        central_res
    )
    plot_region(
        f,
        2,
        1,
        "South GBR",
        ltmp_south,
        south_res;
        showlegend=true,
        legend_row=2,
        legend_col=2
    )
    resize_to_layout!(f)

    save(filename, f)
    return f
end
