# Result analysis and plot generation script.
#
# Loads the calibration products written by scripts/run_loc_calib.jl and writes plots. If
# config.out_dir/params/{calibrated_params.nc,historic_init_cover.nc} are missing but
# config.out_dir/results.dat is present, they are regenerated from results.dat (via
# export_calibration_products) without re-running the BlackBoxOptim search.

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
    calib_data = build_calibration_data(
        rebuild_dom, config.ltmp_reef_data_path, config.composition_path;
        out_dir=OUT_DIR
    )
    init_cover = deserialize(config.init_cover_path)
    best_params = load_calibrated_params(results_fn, length(cfg.sample_bounds))

    export_calibration_products(
        rebuild_dom, init_cover, location_classification.consecutive_classification,
        best_params, cfg.param_idxs, cfg.coral_param_names, cfg.growth_accel_names,
        cfg.depth_atten_names, calib_data.combined_store, cfg.biogroups_ordering;
        out_dir=OUT_DIR
    )
end

# ADRIA folds calibrated_params.nc into the domain's model spec, so param_table returns the
# calibrated scenario directly - no results.dat and no setup_run needed here.
dom = load_domain(config; calib_params_fn=CALIB_PARAMS_FN)

init_cover_da = open_dataset(INIT_COVER_FN).init_coral_cover
@assert size(init_cover_da) == size(dom.init_coral_cover) (
    "$INIT_COVER_FN is $(size(init_cover_da)) but the domain expects " *
    "$(size(dom.init_coral_cover)) - it was written for a different domain."
)
dom.init_coral_cover .= Array(init_cover_da)

regional_data = viz.load_regional_analysis_data(
    dom, config.ltmp_modelled_obs_path, config.ltmp_shp_path
)

calib_data = build_calibration_data(
    dom,
    config.ltmp_reef_data_path, config.composition_path
)

# ----- Run model --------------------------------------------------------------

scen = ADRIA.param_table(dom)
rs_raw = ADRIA.run_model(dom, scen[1, :]; apply_allee_effect=false)
mkpath(OUT_DIR)

# ----- Generate plots ---------------------------------------------------------

# Observation locations map
f_obs_loc_map = viz.plot_observation_locs(calib_data.calibration_store, calib_data.test_store)
save(joinpath(OUT_DIR, "obs_loc_map.png"), f_obs_loc_map)

# Calibrated parameters (from calibrated_params.nc) against their search bounds. The
# initial size-class lambdas are not stored in the NetCDF, so they are not plotted.
plot_cfg = CalibConfig(dom)
f_params = viz.plot_calibrated_params(
    viz.params_vector_from_dataset(
        open_dataset(CALIB_PARAMS_FN), length(plot_cfg.sample_bounds),
        plot_cfg.param_idxs, plot_cfg.coral_param_names
    ),
    plot_cfg.sample_bounds, plot_cfg.param_idxs, plot_cfg.coral_param_names,
    plot_cfg.growth_accel_names, plot_cfg.depth_atten_names, plot_cfg.biogroups_ordering
)
save(joinpath(OUT_DIR, "calibrated_params.png"), f_params)

# Functional group cover proportions
f_taxa_cover = viz.taxa_cover_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_cov.png"), f_taxa_cover)

f_taxa_pop = viz.taxa_population_proportions(rs_raw.raw)
save(joinpath(OUT_DIR, "locs_taxa_pop.png"), f_taxa_pop)

f_size_class = viz.temporal_size_class_proportions(rs_raw.raw; fig_size=(900, 600))
save(joinpath(OUT_DIR, "locs_size.png"), f_size_class)

# Summary stats
n_test_locs = length(calib_data.test_store.ltmp_unique_ids)
test_sccs = [
    viz.collect_error_stats(rs_raw.raw, id, dom; observations=calib_data.test_store).srcc
    for id in 1:n_test_locs
]
test_sccs_sortperm = sortperm(test_sccs)
@info "Three highest SCC test reefs: $(calib_data.test_store.ltmp_unique_ids[test_sccs_sortperm][end-2:end])"
@info "Three lowest  SCC test reefs: $(calib_data.test_store.ltmp_unique_ids[test_sccs_sortperm][1:3])"

stats_calib = viz.collect_error_stats(rs_raw.raw, dom; observations=calib_data.calibration_store)
stats_test = viz.collect_error_stats(rs_raw.raw, dom; observations=calib_data.test_store)
n_calib = length(calib_data.calibration_store.ltmp_unique_ids)
n_test = length(calib_data.test_store.ltmp_unique_ids)
@info "Calibration locations where model outperforms benchmark: $(sum(stats_calib.rmse_model .< stats_calib.rmse_benchmark)) / $n_calib"
@info "Test  locations where model outperforms benchmark: $(sum(stats_test.rmse_model .< stats_test.rmse_benchmark)) / $n_test"
@info "Mean model calib. RMSE: $(mean(stats_calib.rmse_model))"
@info "Mean model test.  RMSE: $(mean(stats_test.rmse_model))"
@info "Mean model calib. SRCC: $(mean(stats_calib.srcc))"
@info "Mean model test.  SRCC: $(mean(stats_test.srcc))"

# Regional comparison plots (03_b)
viz.save_regional_analysis_plots(
    rs_raw.raw, dom, calib_data.calibration_store, calib_data.test_store, OUT_DIR,
    regional_data.ltmp_north, regional_data.ltmp_central, regional_data.ltmp_south,
    regional_data.north_mask, regional_data.central_mask, regional_data.south_mask
)

@info "Plotting metrics"
# Metric analysis plots (03_c)
viz.save_metric_analysis_plots(
    rs_raw.raw, dom, calib_data.calibration_store, calib_data.test_store, OUT_DIR
)
@info "Finished plotting metrics"

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

for (store, store_label) in (
    (calib_data.calibration_store, "Calibration"), (calib_data.test_store, "Test")
)
    ids = store.ltmp_unique_ids
    names = dom.loc_data.GBRMPA_ID[store.ltmp_cover_to_domain]

    rmse_stats = rmse_diff_stats(rs_raw.raw, store, dom)
    pcc_stats = correlation_stats(rs_raw.raw, store, dom; correlation_metric=:pearson)
    srcc_stats = correlation_stats(rs_raw.raw, store, dom; correlation_metric=:spearman)

    println("\n$store_label RMSE diff (Benchmark - Model):")
    println(_extremes_table(ids, names, rmse_stats.diff))
    println("\n$store_label PCC:")
    println(_extremes_table(ids, names, pcc_stats.corr))
    println("\n$store_label SRCC:")
    println(_extremes_table(ids, names, srcc_stats.corr))
end

# Per-location time-series plots (03_d)
root_path = dirname(dirname(dirname(@__FILE__)))
cyc_scens = dom.cyclone_mortality_scens[scenarios=1]
dhw_scens  = dom.dhw_scens[scenarios=1]
disturbances = open_dataset(
    joinpath(root_path, "datasets", "ltmp_data", "disturbances.nc")
).layer

viz.save_location_timeseries_plots(
    rs_raw.raw, dom, calib_data.calibration_store, calib_data.test_store, OUT_DIR,
    dhw_scens, cyc_scens, disturbances
)
nothing
