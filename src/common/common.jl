module common

using ADRIA
using ADRIA: RMEDomain, GDF
using DataFrames
using StatsBase
using Statistics
using CSV
using TOML
using YAXArrays
using Dates
using Serialization
using Random
using Bootstrap

import Base: copy

import ..ADRIACoralBloxCalibration:
    LocationDataStore,
    ltmp_cover_idx_to_domain,
    composition_idx_to_domain,
    get_ltmp_loc_unique_id,
    extract_param_group_idx,
    extract_param_group

include("constants.jl")
include("calib_setup.jl")
include("cover_construction.jl")
include("bootstrap_stats.jl")
include("perf_metrics.jl")
include("dist_mean_penalty.jl")
include("param_bounds.jl")
include("params_extraction.jl")
include("params_metadata.jl")

export LocationDataStore,
       ltmp_cover_idx_to_domain,
       composition_idx_to_domain,
       get_ltmp_loc_unique_id,
       extract_param_group_idx,
       extract_param_group,
       squared_expo_cdf,
       size_class_proportion,
       size_class_distribution,
       construct_location_cover!,
       compute_cover,
       construct_cover!,
       temporal_correlation,
       temporal_correlation_penalty,
       rmse,
       bias,
       MAE,
       MAE_series,
       MAEE,
       MAEE_series,
       temporal_variability,
       year_weighted_error,
       reef_observation_counts,
       calib_func,
       collect_error_stats,
       rmse_diff,
       rmse_diff_stats,
       nse_stats,
       correlation_stats,
       bias_stats,
       block_bootstrap_ci,
       iid_bootstrap_ci,
       bootstrap_median_ci,
       bootstrap_mean_ci,
       average_cc,
       reef_taxa_error,
       bray_curtis,
       reef_error,
       START_YEAR,
       END_YEAR,
       N_TAXA,
       N_SIZE_CLASSES,
       N_LOCATIONS,
       N_TIMESTEPS,
       N_PARAMS,
       STRIDE,
       N_GROWTH_ACCEL_PARAMS,
       LIN_EXT_LB_FACTORS,
       LIN_EXT_UB_FACTOR,
       MB_RATE_LB_FACTOR,
       MB_RATE_UB_FACTOR,
       LIN_EXT_SCALE_LB,
       LIN_EXT_SCALE_UB,
       MB_RATE_SCALE_LB,
       MB_RATE_SCALE_UB,
       GROWTH_ACCEL_STEEPNESS_BOUNDS,
       GROWTH_ACCEL_HEIGHT_BOUNDS,
       GROWTH_ACCEL_MIDPOINT_BOUNDS,
       SC_DIST_LB,
       SC_DIST_UB,
       DHW_TOL_MEAN_LB_FACTOR,
       DHW_TOL_MEAN_UB_FACTOR,
       DHW_TOL_STD_LB_FACTOR,
       DHW_TOL_STD_UB_FACTOR,
       DIST_MEAN_ORDER_SCALE,
       DIST_MEAN_ORDER_WEIGHT,
       dist_mean_order_penalty,
       N_DIST_MEAN_PARAMS,
       N_DIST_STD_PARAMS,
       DEPTH_ATTEN_PARAM_NAMES,
       N_DEPTH_ATTEN_PARAMS,
       EFF_DHW_BASE_BOUNDS,
       EFF_DHW_MIX_BOUNDS,
       scale_factor_vec_to_array,
       scale_factor_array_to_vec,
       generate_scale_factor_names,
       accel_params_array_to_vec,
       accel_params_vec_to_array,
       generate_growth_accel_names,
       target_param_names,
       sc_fg_param_idxs,
       set_bounds!,
       flatten_group_size,
       insert_init_loc_cover!,
       get_scale_factors,
       setup_run,
       build_params_dataset,
       build_init_cover_dataset,
       export_calibration_products,
       CalibrationConfig,
       load_config,
       load_domain,
       load_location_classification,
       PARAM_SCHEMA_VERSION,
       write_params_metadata,
       load_calibrated_params

end
