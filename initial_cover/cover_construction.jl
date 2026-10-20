using ADRIA

"""
CDF of squared exponential distribution with X ~ π/4 * E^2. Where E is an exponential distributed
"""
function squared_expo_cdf(lambda::Float64, x::Float64)::Float64
    return 1 - exp(-4 * lambda * sqrt(x) / π)
end

function size_class_proportion(
    lambda::Float64,
    lower_bound::Float64,
    upper_bound::Float64
)::Float64
    return squared_expo_cdf(lambda, upper_bound) - squared_expo_cdf(lambda, lower_bound)
end

function size_class_distribution(
    lambda::Float64, bin_edges::AbstractVector{Float64}
)::Vector{Float64}
    area_bounds::Vector{Float64} = (x -> π / 4 * x^2).(bin_edges)
    area_props::Vector{Float64} = size_class_proportion.(
        lambda,
        area_bounds[1:end-1],
        area_bounds[2:end]
    )
    area_props = area_props ./ sum(area_props)
    return area_props
end

"""
    construct_location_cover!(preallocated, location_sample, bin_edges)::Nothing

Location Sample element description

1        - Location relative habitable cover
[2 - 6]  - Taxonomy cover weights
[7 - 11] - Taxonomy size exponential parametrisation

# Arguments
- `preallocated` : An abstract Matrix of size [n_sizes ⋅ n_taxa]
- `location_sample` : An abstract vector containing the embedded from the sampler
- `bin_edges` : An abstract vector defining the diameter bounds of the taxa
"""
function construct_location_cover!(
    preallocated::AbstractMatrix{Float64},
    location_sample::AbstractVector{Float64},
    bin_edges::AbstractMatrix{Float64}
)::Nothing
    n_taxa::Int64 = size(bin_edges, 1)
    taxonomy_covers::Vector{Float64} = location_sample[2:(1+n_taxa)] .* (
        location_sample[1] ./ sum(location_sample[2:(1+n_taxa)])
    )

    for taxa in 1:n_taxa
        preallocated[:, taxa] .= size_class_distribution(location_sample[6+taxa], bin_edges[taxa, :])
    end

    preallocated .*= taxonomy_covers'
    return nothing
end

"""
    compute_cover(vec_sample, location_types, n_locations)::Matrix{Float64}

Compute per-location initial coral cover from a flat location-type sample vector.
Returns a plain Matrix{Float64} of shape [n_sizes*n_taxa, n_locations].
"""
function compute_cover(
    vec_sample::Vector{Float64},
    location_types::AbstractVector{Int64},
    n_locations::Int64
)::Matrix{Float64}
    n_location_types = maximum(location_types)
    bin_edges::Matrix{Float64} = ADRIA.bin_edges().data
    n_taxa::Int64 = size(bin_edges, 1)
    n_size_classes::Int64 = size(bin_edges, 2) - 1
    temporary_cover::Matrix{Float64} = zeros(Float64, n_size_classes, n_taxa)
    cover::Matrix{Float64} = zeros(Float64, n_size_classes * n_taxa, n_locations)
    location_mask::BitVector = BitVector(undef, n_locations)

    stride::Int64 = 1 + 2 * n_taxa
    for loc_type in 1:n_location_types
        @views construct_location_cover!(
            temporary_cover,
            vec_sample[(1+stride*(loc_type-1)):(stride*loc_type)],
            bin_edges
        )
        location_mask .= location_types .== loc_type
        cover[:, location_mask] .= vec(temporary_cover)
    end
    return cover
end

function construct_cover!(
    dom::Domain,
    vec_sample::Vector{Float64},
    location_types::AbstractVector{Int64}
)::Nothing
    dom.init_coral_cover .= compute_cover(vec_sample, location_types, ADRIA.n_locations(dom))
    return nothing
end
