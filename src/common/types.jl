"""
Shared data types for the calibration pipeline.

Extracted from common.jl so that ADRIACoralBloxCalibration (and its submodules) can
define these types once and all scripts can obtain them via `using ADRIACoralBloxCalibration`.
"""

"""
Immutable struct providing a consistent interface for indexing between the arrays
that refer to the same locations across different data sources.
"""
struct LocationDataStore
    domain_gpkg::DataFrame                          # ADRIA domain location table (one row per reef site)
    ltmp_unique_ids::Vector{String}                 # RME_UNIQUE_ID for each LTMP survey site in this store
    ltmp_coral_cover::Array{Union{Missing,Float64},2}  # observed total coral cover [location × timestep]
    coral_composition::Array{Union{Missing,Float64},3} # observed species-group cover [timestep × group × location]
    ltmp_cover_to_domain::Vector{Int64}             # maps each LTMP cover location index → domain row index
    composition_to_domain::Vector{Int64}            # maps each composition location index → domain row index
end

function ltmp_cover_idx_to_domain(loc_data_store::LocationDataStore, ltmp_cover_idx::Int64)
    return loc_data_store.ltmp_cover_to_domain[ltmp_cover_idx]
end

function composition_idx_to_domain(loc_data_store::LocationDataStore, composition_idx::Int64)
    return loc_data_store.composition_to_domain[composition_idx]
end

function get_ltmp_loc_unique_id(loc_data_store::LocationDataStore, ltmp_idx::Int64)::String
    return loc_data_store.ltmp_unique_ids[ltmp_idx]
end

function get_composition_loc_unique_id(
    loc_data_store::LocationDataStore, composition_idx::Int64
)::String
    return loc_data_store.coral_composition.location[composition_idx]
end

function extract_param_group_idx(model_spec::DataFrame, needle::String)::Vector{Int64}
    needle_pos = contains.(string.(model_spec.fieldname), needle)
    return findall(needle_pos)
end
function extract_param_group_idx(
    model::ADRIA.Model, component, needle::String
)::Vector{Int64}
    comp_params = ADRIA.component_params(model, component)
    return extract_param_group_idx(comp_params, needle)
end

function extract_param_group(model::ADRIA.Model, component, needle::String)::DataFrame
    group_pos = extract_param_group_idx(model, component, needle)
    comp_params = ADRIA.component_params(model, component)
    return comp_params[group_pos, :]
end
function extract_param_group(model_spec::DataFrame, needle::String)::DataFrame
    group_pos = extract_param_group_idx(model_spec, needle)
    return model_spec[group_pos, :]
end
