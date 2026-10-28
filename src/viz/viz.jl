module viz

using ADRIA
using ADRIA: GDF, AG
using ADRIA: GDF.GeoInterface as GI

using CairoMakie
using GeometryBasics: Polygon, Point2f
using ColorSchemes
using GeoMakie
using DataFrames
using Dates
using Statistics
using StatsBase
using YAXArrays
using ProgressMeter
using CSV

import ..ADRIACoralBloxCalibration:
    LocationDataStore,
    ltmp_cover_idx_to_domain,
    composition_idx_to_domain,
    get_ltmp_loc_unique_id

using ..common

# ----- Visual constants -------------------------------------------------------

const GBRMPA_MAINLAND_PATH = joinpath(
    dirname(dirname(@__DIR__)),
    "datasets", "spatial_data", "Great_Barrier_Reef_Features.geojson"
)
const GBRMPA_MAINLAND_GPKG = GDF.read(GBRMPA_MAINLAND_PATH)
# Convert mainland geometries to pure Julia at load time while the GDAL layer is open.
# ArchGDAL geometry ptrs from GDF.read are borrowed GDAL references that become stale
# after the layer closes — preconversion avoids NULL ptr errors at render time.
# Handles both wkbPolygon and wkbMultiPolygon features present in the GeoJSON.
_gdal_ring_to_poly(ring) =
    Polygon([Point2f(AG.getx(ring, j), AG.gety(ring, j)) for j in 0:AG.ngeom(ring)-1])
const GBRMPA_MAINLAND_POLYS = mapreduce(vcat, GBRMPA_MAINLAND_GPKG.geometry) do g
    if AG.getgeomtype(g) == AG.wkbMultiPolygon
        [_gdal_ring_to_poly(AG.getgeom(AG.getgeom(g, i), 0)) for i in 0:AG.ngeom(g)-1]
    else
        [_gdal_ring_to_poly(AG.getgeom(g, 0))]
    end
end

const inch = 96
const pt = 4 / 3
const cm = inch / 2.54

const COLORMAPS = (
    ltmp_disturbances=:seaborn_bright6,
    taxa=:seaborn_colorblind6
)

const COLORS = (
    model_vs_obs_color_model="#1E88E5",
    model_vs_obs_color_obs="#D81B60",
    model_dist_color_dhw="#eb4034",
    model_dist_color_dist="#31004D",
    model_line_color="#1E88E5",
    model_band_color=("#1E88E5", 0.4),
    historic_line_color="#D81B60",
    historic_band_color=("#D81B60", 0.3),
    NamedTuple(zip(
        Symbol.(ADRIA.functional_group_names()),
        ColorSchemes.seaborn_colorblind6[1:5],
    )
    )...
)

const FONT_SIZES = (
    xlabel=20,
    ylabel=20,
    title=20
)

const BASE_AXIS_OPTS::Dict = Dict(
    :xlabelsize => FONT_SIZES.xlabel,
    :ylabelsize => FONT_SIZES.ylabel,
    :xticks => collect(START_YEAR:END_YEAR),
    :xticklabelrotation => pi / 6,
    :limits => ((START_YEAR - 0.5, END_YEAR + 0.5), nothing),
)

const FIG_SIZE::Dict = Dict(
    :map => (500, 600),
)

# ----- Plot helpers and submodule files ---------------------------------------

include("plot_helpers.jl")
include("regional_data_loader.jl")

include("location_comparison.jl")
include("regional_comparison.jl")
include("metrics.jl")
include("progress.jl")
include("initial_cover_plots.jl")
include("regional_analysis.jl")

# ----- Analysis script wrappers (from 03_b, 03_c, 03_d) ----------------------

include("regional_analysis_plots.jl")
include("metric_analysis_plots.jl")
include("location_timeseries_plots.jl")
include("calibrated_params_plot.jl")
include("heat_tolerance_plot.jl")

# ----- Exports ----------------------------------------------------------------

export LocationDataStore,
       ltmp_cover_idx_to_domain,
       collect_error_stats,
       rmse_diff,
       correlation_stats,
       region_stats,
       plot_regional_comparison,
       plot_all_regions,
       plot_location_comparison,
       plot_location_comparison!,
       plot_modelled_v_ltmp!,
       plot_rmse_diff_map,
       plot_correlation_map,
       plot_metrics_heatmap,
       plot_rmse_scatter,
       plot_nse_scatter,
       plot_correlation_scatter,
       plot_observation_locs,
       plot_taxa_props,
       taxa_cover_proportions,
       taxa_population_proportions,
       progress_run,
       save_regional_analysis_plots,
       save_metric_analysis_plots,
       save_location_timeseries_plots,
       plot_calibrated_params,
       mean_heat_tolerance_by_group_size,
       plot_heat_tolerance_timeseries,
       RegionalAnalysisData,
       load_regional_analysis_data

end
