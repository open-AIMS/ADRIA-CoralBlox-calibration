# Result analysis and plot generation script.
#
# Loads the calibration products written by scripts/run_loc_calib.jl and writes plots into
# config.out_dir. If config.out_dir/params/{calibrated_params.nc,historic_init_cover.nc} are
# missing but config.out_dir/results.dat is present, they are regenerated from results.dat
# (via export_calibration_products) without re-running the BlackBoxOptim search.
#
# Usage (from the repository root):
#   julia --project=. scripts/plot/viz_results.jl                 plot only the figures that are missing
#   julia --project=. scripts/plot/viz_results.jl --force         redo every figure
#   julia --project=. scripts/plot/viz_results.jl NAME [NAME...]  redo exactly these figures
#   julia --project=. scripts/plot/viz_results.jl --list          show every figure and whether it exists
#
# A figure counts as present when all its output files exist (for figure groups written to a
# folder, when that folder is not empty). Anything expensive (domain, model run) is only
# loaded when a selected figure needs it.

using ADRIA
using Statistics
using YAXArrays
using DataFrames
using Serialization: deserialize

using ADRIACoralBloxCalibration
import ADRIACoralBloxCalibration.viz
using ADRIACoralBloxCalibration.common
using ADRIACoralBloxCalibration.calibration
using CairoMakie: save

config = load_config(joinpath(@__DIR__, "..", "..", "config.toml"))

OUT_DIR = config.out_dir
PARAMS_DIR = joinpath(OUT_DIR, "params")
CALIB_PARAMS_FN = joinpath(PARAMS_DIR, "calibrated_params.nc")
INIT_COVER_FN = joinpath(PARAMS_DIR, "historic_init_cover.nc")

# ----- Lazily loaded inputs ---------------------------------------------------

_cache = Dict{Symbol,Any}()
_once(f, key::Symbol) = get!(f, _cache, key)

function ensure_products!()
    _once(:products) do
        if !isfile(CALIB_PARAMS_FN) || !isfile(INIT_COVER_FN)
            results_fn = joinpath(OUT_DIR, "results.dat")
            isfile(results_fn) || error(
                "Missing calibration product(s) under $PARAMS_DIR and no $results_fn to rebuild " *
                "them from - run scripts/run_loc_calib.jl first."
            )

            @info "Calibration products missing; rebuilding from $results_fn"
            rebuild_dom = load_domain(config)
            location_classification = load_location_classification(config.loc_class_path)
            cfg = CalibConfig(rebuild_dom)
            rebuild_calib_data = build_calibration_data(
                rebuild_dom, config.ltmp_reef_data_path, config.composition_path;
                out_dir=OUT_DIR
            )
            init_cover = deserialize(config.init_cover_path)
            best_params = load_calibrated_params(results_fn, length(cfg.sample_bounds))

            export_calibration_products(
                rebuild_dom, init_cover, location_classification.consecutive_classification,
                best_params, cfg.param_idxs, cfg.coral_param_names, cfg.growth_accel_names,
                cfg.depth_atten_names, rebuild_calib_data.combined_store, cfg.biogroups_ordering;
                dist_std_names=cfg.dist_std_names, dist_mean_names=cfg.dist_mean_names,
                out_dir=OUT_DIR
            )
        end
        true
    end
end

# ADRIA folds calibrated_params.nc into the domain's model spec, so param_table returns the
# calibrated scenario directly - no results.dat and no setup_run needed here.
get_dom() = _once(:dom) do
    ensure_products!()
    dom = load_domain(config; calib_params_fn=CALIB_PARAMS_FN)

    init_cover_da = open_dataset(INIT_COVER_FN).init_coral_cover
    @assert size(init_cover_da) == size(dom.init_coral_cover) (
        "$INIT_COVER_FN is $(size(init_cover_da)) but the domain expects " *
        "$(size(dom.init_coral_cover)) - it was written for a different domain."
    )
    dom.init_coral_cover .= Array(init_cover_da)
    dom
end

get_regional_data() = _once(:regional_data) do
    viz.load_regional_analysis_data(get_dom(), config.ltmp_modelled_obs_path, config.ltmp_shp_path)
end

get_calib_data() = _once(:calib_data) do
    build_calibration_data(get_dom(), config.ltmp_reef_data_path, config.composition_path)
end

get_scen() = _once(:scen) do
    dom = get_dom()
    scen = ADRIA.param_table(dom)
    # The file's `dist_mean` is the end-of-history tolerance (forecast start); the historical
    # run must start from the calibrated initial value
    ds = open_dataset(CALIB_PARAMS_FN)
    if :dist_mean_initial in keys(ds.cubes)
        scen[!, CalibConfig(dom).dist_mean_names] .=
            repeat(Array(ds.dist_mean_initial); inner=N_SIZE_CLASSES)'
    end
    scen
end

get_run() = _once(:run) do
    ADRIA.setup()  # run_model calls setup(), which would otherwise overwrite the ENV flags below
    ENV["ADRIA_LOG_DHW_TOLS"] = "true"  # needed for the heat tolerance time series
    ENV["ADRIA_LOG_COVER"] = "true"
    ADRIA.run_model(get_dom(), get_scen()[1, :]; apply_allee_effect=false)
end

# ----- Figures ----------------------------------------------------------------

struct PlotJob
    name::String
    outputs::Vector{String}  # relative to OUT_DIR; a trailing "/" marks a folder
    default::Bool            # whether it is part of the default (missing-only) run
    run::Function
end

function plot_obs_loc_map()
    cd = get_calib_data()
    save(joinpath(OUT_DIR, "obs_loc_map.png"), viz.plot_observation_locs(cd.calibration_store, cd.test_store))
end

# Calibrated parameters (from calibrated_params.nc) against their search bounds. The
# initial size-class lambdas are not stored in the NetCDF, so they are not plotted.
function plot_calibrated_params()
    plot_cfg = CalibConfig(get_dom())
    params_ds = open_dataset(CALIB_PARAMS_FN)
    f = viz.plot_calibrated_params(
        viz.params_vector_from_dataset(
            params_ds, length(plot_cfg.sample_bounds),
            plot_cfg.param_idxs, plot_cfg.coral_param_names
        ),
        plot_cfg.sample_bounds, plot_cfg.param_idxs, plot_cfg.coral_param_names,
        plot_cfg.growth_accel_names, plot_cfg.depth_atten_names, plot_cfg.biogroups_ordering;
        dist_mean_end=:dist_mean_initial in keys(params_ds.cubes) ?
                      Array(permutedims(params_ds.dist_mean, (:functional_group, :size_class))) :
                      nothing
    )
    save(joinpath(OUT_DIR, "calibrated_params.png"), f)
end

plot_taxa_cover() = save(joinpath(OUT_DIR, "locs_taxa_cov.png"), viz.taxa_cover_proportions(get_run().raw))
plot_taxa_pop() = save(joinpath(OUT_DIR, "locs_taxa_pop.png"), viz.taxa_population_proportions(get_run().raw))
plot_size_class() = save(
    joinpath(OUT_DIR, "locs_size.png"), viz.temporal_size_class_proportions(get_run().raw; fig_size=(900, 600))
)

# Mean heat tolerance by functional group and size class (last size class excluded)
function plot_heat_tolerance()
    dom, scen, rs = get_dom(), get_scen(), get_run()
    mean_tol = viz.mean_heat_tolerance_by_group_size(
        rs.coral_dhw_log, rs.coral_cover_log, Array(dom.init_coral_cover),
        reshape(ADRIA.to_coral_spec(scen[1, :]).dist_mean, :, 1) .* ones(1, size(dom.init_coral_cover, 2)),
        ADRIA.loc_k_area(dom)
    )
    f = viz.plot_heat_tolerance_timeseries(
        mean_tol; years=collect(ADRIA.timesteps(dom)), n_sizes=dom.coral_growth.n_sizes, exclude_last=1
    )
    save(joinpath(OUT_DIR, "heat_tolerance_by_group_size.png"), f)
end

# Regional comparison plots (03_b)
function plot_regional()
    cd, rd = get_calib_data(), get_regional_data()
    viz.save_regional_analysis_plots(
        get_run().raw, get_dom(), cd.calibration_store, cd.test_store, OUT_DIR,
        rd.ltmp_north, rd.ltmp_central, rd.ltmp_south,
        rd.north_mask, rd.central_mask, rd.south_mask
    )
end

# Metric analysis plots (03_c)
function plot_metrics()
    cd = get_calib_data()
    @info "Plotting metrics"
    viz.save_metric_analysis_plots(get_run().raw, get_dom(), cd.calibration_store, cd.test_store, OUT_DIR)
    @info "Finished plotting metrics"
end

# Per-location time-series plots (03_d)
function plot_location_timeseries()
    dom, cd = get_dom(), get_calib_data()
    root_path = dirname(dirname(dirname(@__FILE__)))
    disturbances = open_dataset(joinpath(root_path, "datasets", "ltmp_data", "disturbances.nc")).layer
    viz.save_location_timeseries_plots(
        get_run().raw, dom, cd.calibration_store, cd.test_store, OUT_DIR,
        dom.dhw_scens[scenarios=1], dom.cyclone_mortality_scens[scenarios=1], disturbances
    )
end

# Extreme-score reef tables (RMSE diff, PCC, SRCC; test & calibration)
function _extremes_table(ids, names, scores; n::Int=2)
    order = sortperm(scores)
    lo, hi = order[1:n], order[(end - n + 1):end]
    idx = vcat(lo, hi)
    return DataFrame(
        Group=vcat(fill("Lowest", n), fill("Highest", n)),
        ID=ids[idx],
        Name=names[idx],
        Score=round.(scores[idx]; digits=3),
    )
end

# Console summary: error statistics and extreme-score reefs. Writes no files.
function print_summary()
    dom, cd, raw = get_dom(), get_calib_data(), get_run().raw

    n_test_locs = length(cd.test_store.ltmp_unique_ids)
    test_sccs = [
        viz.collect_error_stats(raw, id, dom; observations=cd.test_store).srcc
        for id in 1:n_test_locs
    ]
    test_sccs_sortperm = sortperm(test_sccs)
    @info "Three highest SCC test reefs: $(cd.test_store.ltmp_unique_ids[test_sccs_sortperm][end-2:end])"
    @info "Three lowest  SCC test reefs: $(cd.test_store.ltmp_unique_ids[test_sccs_sortperm][1:3])"

    stats_calib = viz.collect_error_stats(raw, dom; observations=cd.calibration_store)
    stats_test = viz.collect_error_stats(raw, dom; observations=cd.test_store)
    n_calib = length(cd.calibration_store.ltmp_unique_ids)
    n_test = length(cd.test_store.ltmp_unique_ids)
    @info "Calibration locations where model outperforms benchmark: $(sum(stats_calib.rmse_model .< stats_calib.rmse_benchmark)) / $n_calib"
    @info "Test  locations where model outperforms benchmark: $(sum(stats_test.rmse_model .< stats_test.rmse_benchmark)) / $n_test"
    @info "Mean model calib. RMSE: $(mean(stats_calib.rmse_model))"
    @info "Mean model test.  RMSE: $(mean(stats_test.rmse_model))"
    @info "Mean model calib. SRCC: $(mean(stats_calib.srcc))"
    @info "Mean model test.  SRCC: $(mean(stats_test.srcc))"

    for (store, store_label) in ((cd.calibration_store, "Calibration"), (cd.test_store, "Test"))
        ids = store.ltmp_unique_ids
        names = dom.loc_data.GBRMPA_ID[store.ltmp_cover_to_domain]

        rmse_stats = rmse_diff_stats(raw, store, dom)
        pcc_stats = correlation_stats(raw, store, dom; correlation_metric=:pearson)
        srcc_stats = correlation_stats(raw, store, dom; correlation_metric=:spearman)

        println("\n$store_label RMSE diff (Benchmark - Model):")
        println(_extremes_table(ids, names, rmse_stats.diff))
        println("\n$store_label PCC:")
        println(_extremes_table(ids, names, pcc_stats.corr))
        println("\n$store_label SRCC:")
        println(_extremes_table(ids, names, srcc_stats.corr))
    end
end

JOBS = [
    PlotJob("obs_loc_map", ["obs_loc_map.png"], true, plot_obs_loc_map),
    PlotJob("calibrated_params", ["calibrated_params.png"], true, plot_calibrated_params),
    PlotJob("taxa_cover", ["locs_taxa_cov.png"], true, plot_taxa_cover),
    PlotJob("taxa_pop", ["locs_taxa_pop.png"], true, plot_taxa_pop),
    PlotJob("size_class", ["locs_size.png"], true, plot_size_class),
    PlotJob("heat_tolerance", ["heat_tolerance_by_group_size.png"], true, plot_heat_tolerance),
    PlotJob("regional", ["regional_analysis/"], true, plot_regional),
    PlotJob("metrics", ["metrics/"], true, plot_metrics),
    PlotJob("location_timeseries", ["calibration_locations/", "test_locations/"], true, plot_location_timeseries),
    PlotJob("summary", String[], false, print_summary),
]

function is_present(rel::String)::Bool
    path = joinpath(OUT_DIR, rel)
    return endswith(rel, "/") ? (isdir(path) && !isempty(readdir(path))) : isfile(path)
end
is_missing(job::PlotJob)::Bool = !all(is_present, job.outputs)

# ----- Selection and run ------------------------------------------------------

flags = filter(startswith("--"), ARGS)
requested = filter(!startswith("--"), ARGS)
valid_flags = ["--force", "--list", "--help"]
bad_flags = setdiff(flags, valid_flags)
isempty(bad_flags) || error("Unknown option(s) $(bad_flags); valid: $(join(valid_flags, ", "))")

job_names = [job.name for job in JOBS]
unknown = setdiff(requested, job_names)
isempty(unknown) || error("Unknown figure(s) $(unknown); valid: $(join(job_names, ", "))")

if "--help" in flags
    println("Usage: viz_results.jl [--force] [--list] [NAME ...]\nFigures: $(join(job_names, ", "))")
elseif "--list" in flags
    println("Figures in $OUT_DIR:")
    for job in JOBS
        status = isempty(job.outputs) ? "console only (run by name)" : (is_missing(job) ? "MISSING" : "present")
        println("  ", rpad(job.name, 22), status)
    end
else
    selected = if !isempty(requested)
        filter(job -> job.name in requested, JOBS)
    elseif "--force" in flags
        filter(job -> job.default, JOBS)
    else
        filter(job -> job.default && is_missing(job), JOBS)
    end

    if isempty(selected)
        @info "All figures already exist in $OUT_DIR; name a figure or use --force to redo it."
    else
        mkpath(OUT_DIR)
        @info "Plotting: $(join([job.name for job in selected], ", "))"
        for job in selected
            @info "Figure: $(job.name)"
            job.run()
        end
    end
end
nothing
