"""
Central named constants for model dimensions and parameter bounds.

This file must be included before any file that references these constants
(see `common.jl` / `plot.jl` include order). Values mirror the literals
previously hardcoded across `cover_construction.jl`, `param_bounds.jl`,
`plotting.jl`, and the numbered pipeline scripts — no numeric value should be
changed here without re-deriving every dependent reshape/index.
"""

# --- Core domain dimensions -------------------------------------------------

# Calibration/simulation period (inclusive). Earliest possible start for
# ReefModDomain is 2008.
const START_YEAR = 2008
const END_YEAR = 2022

# Number of coral functional groups / taxa (ADRIA.functional_group_names())
const N_TAXA = 5

# Number of coral size classes (ADRIA.bin_edges() / bin_widths())
const N_SIZE_CLASSES = 7

# Number of reef locations in the RME domain
const N_LOCATIONS = 3806

# Number of timesteps in a single model run output (raw_data first dimension)
const N_TIMESTEPS = 15

# --- Derived / parameter-vector dimensions ----------------------------------

# Number of flattened coral parameters per location sample (N_TAXA * N_SIZE_CLASSES)
const N_PARAMS = N_TAXA * N_SIZE_CLASSES

# Stride of a single location's sample vector in `vec_sample`:
# 1 relative cover + N_TAXA taxonomy weights + N_TAXA size-distribution lambdas
const STRIDE = 1 + 2 * N_TAXA

# Number of growth-acceleration parameters per biogroup (steepness, height, midpoint)
const N_GROWTH_ACCEL_PARAMS = 3

# --- Parameter bound constants (param_bounds.jl) ----------------------------

# Linear extension / background mortality bounds: +/-20% of the mean. From size class 2
# onwards the linear extension upper bound is also limited to the largest base value that
# keeps worst-case growth within the class width (see CalibConfig).
const LIN_EXT_LB_FACTOR = 0.8
const LIN_EXT_UB_FACTOR = 1.2
const MB_RATE_LB_FACTOR = 0.8
const MB_RATE_UB_FACTOR = 1.2

# Linear extension scale factor bounds
const LIN_EXT_SCALE_LB = 0.7
const LIN_EXT_SCALE_UB = 1.5

# Background mortality scale factor bounds
const MB_RATE_SCALE_LB = -1.0
const MB_RATE_SCALE_UB = 1.0

# Growth acceleration parameter bounds (steepness, height, midpoint)
const GROWTH_ACCEL_STEEPNESS_BOUNDS = (-20.0, -15.0)
const GROWTH_ACCEL_HEIGHT_BOUNDS = (0.0, 2.0)
const GROWTH_ACCEL_MIDPOINT_BOUNDS = (0.0, 0.3)

# Size-class distribution (lambda) bounds
const SC_DIST_LB = 0.25
const SC_DIST_UB = 2.0

# Depth attenuation bounds (ADRIA `effective_dhw_at_depth`). Two GBR-wide scalars, in the
# order they occupy in the parameter vector.
#
# Both are free parameters with no empirical anchor (see ADRIA's `effective_dhw_at_depth`
# docstring), so these boxes are deliberately wide rather than centred on the ADRIA defaults
# of 0.04 / 12.0. Note that RME assigns depth_med = 7 m to every location, so there is no
# spatial depth gradient to identify them against - they act on the whole domain uniformly
# and are only weakly constrained by the cover trajectories. Treat calibrated values as
# "what best reproduces cover", not as measured optical properties.
const DEPTH_ATTEN_PARAM_NAMES = ["eff_dhw_base", "eff_dhw_mix"]
const EFF_DHW_BASE_BOUNDS = (0.0, 1.0)
const EFF_DHW_MIX_BOUNDS = (0.0, 40.0)
const N_DEPTH_ATTEN_PARAMS = length(DEPTH_ATTEN_PARAM_NAMES)

# DHW-tolerance distribution bounds, as factors of the ADRIA literature default.
#
# Both `dist_mean` and `dist_std` are calibrated per functional group only (5 values each),
# broadcast across that group's 7 size classes when written into the scenario. That matches
# ADRIA's own priors, which are `repeat(<5 group values>; inner=n_sizes)` for both - already
# constant within a group - so the 30 extra degrees of freedom per parameter were
# unsupported by the literature they came from.
const DHW_TOL_MEAN_LB_FACTOR = 0.5
const DHW_TOL_MEAN_UB_FACTOR = 1.5

# `dist_std` is deliberately tighter than `dist_mean`. The two trade off along a ridge - many
# (mean, std) pairs produce near-identical bleaching at a given DHW - and letting both hit
# opposite corners of a +/-50% box spans a coefficient of variation of 0.26 to 2.32 (the
# ADRIA defaults all sit at CV 0.774). The high end is degenerate: a tolerance distribution
# that wide, truncated to [HEAT_LB, mean+HEAT_UB], is nearly flat, making bleaching almost
# insensitive to DHW. Holding std to +/-25% narrows the span to 0.39-1.94, removing the
# CV > 2 corner. Std is genetic variance and far less likely to be 2x wrong than the mean.
const DHW_TOL_STD_LB_FACTOR = 0.75
const DHW_TOL_STD_UB_FACTOR = 1.25

# Number of calibrated `dist_mean` / `dist_std` values: one per functional group each
const N_DIST_MEAN_PARAMS = N_TAXA
const N_DIST_STD_PARAMS = N_TAXA

# Serialized-parameter-vector schema version, recorded in each result's `.meta.toml`
# sidecar (see `write_params_metadata`/`load_calibrated_params`). Bump whenever
# CalibConfig's param_idxs/sample_bounds shape changes, or ADRIA's size class edges change
# (per-size-class values are meaningless under different edges), so older serialized
# results are rejected rather than silently reinterpreted under the current schema.
const PARAM_SCHEMA_VERSION = 7
