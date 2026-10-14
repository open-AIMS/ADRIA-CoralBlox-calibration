"""
Diagnostic script: checks whether `reef_error`'s per-year error grows with year (i.e. whether
error compounds over a free-running ADRIA simulation, since it integrates forward from
`init_cover` with no re-anchoring to observations mid-run).

Originally run 2026-07-23 against a real calibrated parameter set, giving
`Spearman(year, error) ≈ 0.89` — a strong, close-to-monotonic upward trend (2009: 0.038 →
2022: 0.150). This result motivated switching `_obj_func`'s `reef_perf` term from
`temporal_variability` (flat mean/std blend) to `year_weighted_error` (linearly increasing
per-year weight) — see `year_weighted_error`'s docstring in
[src/common/perf_metrics.jl](../src/common/perf_metrics.jl) for the full rationale.

Re-run this whenever the loss function or model changes materially, to confirm the drift
assumption still holds.
"""

using ADRIACoralBloxCalibration
using ADRIACoralBloxCalibration.common
using ADRIACoralBloxCalibration.calibration
using ADRIACoralBloxCalibration.viz
using ADRIA
using Serialization
using StatsBase

config = load_config(joinpath(@__DIR__, "..", "config.toml"))
dom = load_domain(config)

cfg = CalibConfig(dom)
calib_data = build_calibration_data(
    dom, config.ltmp_reef_data_path, config.composition_path;
    out_dir=config.out_dir
)

# Pick a parameter vector to inspect. Swap this for any serialised result you want to check
# (e.g. an in-progress run's intermediate save), not necessarily the final "best" one.
params = load_calibrated_params(
    joinpath(config.out_dir, "results.dat"), length(cfg.sample_bounds)
)

res = progress_run(
    params, dom,
    cfg.coral_param_names, cfg.growth_accel_names, cfg.depth_atten_names,
    cfg.dist_std_names, cfg.dist_mean_names, cfg.param_idxs,
    calib_data.calibration_store, cfg.biogroups_ordering
)

loc_k_areas = ADRIA.loc_k_area(dom)
loc_areas = ADRIA.loc_area(dom)
loc_cover = dropdims(sum(res.raw; dims=(2, 3)); dims=(2, 3)) .* loc_k_areas' ./ loc_areas'

err_series = reef_error(loc_cover; observations=calib_data.calibration_store)

years = 2008:2022
for (yr, err) in zip(years, err_series)
    println("$(yr): $(round(err, digits=4))")
end

# Spearman correlation between year index and error.
# Close to +1 => error grows ~monotonically with year (supports drift/lead-time weighting).
# Close to 0 => no clear trend (a couple of bad years likely dominate instead).
println("\nSpearman(year, error) = $(corspearman(collect(1:length(err_series)), err_series))")
