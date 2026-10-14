module calibration

using ADRIA
using ADRIA: GDF
using DataFrames
using Random
using Statistics
using YAXArrays
using Serialization
using Dates
using CSV
using BlackBoxOptim
import CairoMakie

import ..ADRIACoralBloxCalibration: LocationDataStore

using ..common
using ..viz

include("config.jl")
include("data_split.jl")
include("runner.jl")

export CalibConfig,
       CalibrationData,
       build_calibration_data,
       validate_gbr_wide_scalar_mean,
       validate_linear_extension_coefficients,
       run_calibration

end
