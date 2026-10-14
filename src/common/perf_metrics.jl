using StatsBase

function temporal_correlation(err_series)::Float64
    return cor(1:length(err_series), err_series)
end

function temporal_correlation_penalty(err_series; threshold::Float64=0.3)::Float64
    corr::Float64 = temporal_correlation(err_series)
    return abs(corr) < threshold ? 1.0 : 2 * abs(corr) + 1.0
end

"""
    rmse(modelled, observed)

Calculate Root Mean Squared Error
"""
function rmse(modelled, observed)
    return sqrt(mean((modelled .- observed) .^ 2.0))
end

"""
    bias(modelled, observed)
"""
function bias(modelled, observed)
    return mean(modelled .- observed)
end

"""
Mean Absolute Error
"""
MAE(sim, obs) = mean(abs.(sim .- obs))

"""
Element-wise absolute error, i.e. `MAE` without the final `mean` reduction.
"""
MAE_series(sim, obs) = abs.(sim .- obs)

"""
Mean Absolute Exponential Error.

Assign error that increases exponentially with distance to observed/"true" data.
"""
function MAEE(sim, obs)
    abs_err = abs.(sim .- obs)
    mean(ℯ .^ ((abs_err) .* (1.0 .+ abs_err ./ 1.0)) .- 1.0)
end

function MAEE_series(sim, obs)
    abs_err = abs.(sim .- obs)
    return ℯ .^ ((abs_err) .* (1.0 .+ abs_err / 1.0)) .- 1.0
end

function temporal_variability(x::AbstractVector{<:Real}; w=[0.9, 0.1])
    return mean([mean(x), std(x)], weights(w))
end

"""
    reef_observation_counts(observations::LocationDataStore)::Vector{Int}

Number of LTMP reefs with a non-missing cover observation in each year.
"""
function reef_observation_counts(observations::LocationDataStore)::Vector{Int}
    return vec(count(!ismissing, observations.ltmp_coral_cover; dims=1))
end

"""
    year_weighted_error(err_series::AbstractVector{<:Real}, counts::AbstractVector{<:Real})::Float64

Weighted mean of a per-year error series (as returned by `reef_error`). Weight = year index
(later years matter more) × `sqrt(counts)` (years with more reporting reefs are less noisy,
so trusted more; `sqrt` rather than raw count since reefs surveyed in the same year plausibly
don't fail independently — a bad year hits several reefs via the same disturbance forcing).
Years with no observation data (`err_series[i] == 0`) are excluded entirely, not zero-weighted.

# Why year-weighting
ADRIA runs forward from `init_cover` with no re-anchoring to observations mid-run, so
parameter error compounds year over year — a flat/endpoint-only score lets the optimizer
nail the easy early years while drifting badly by the end. Verified empirically:
`Spearman(year, error) ≈ 0.89` on a real calibrated run (2009: 0.038 → 2022: 0.150).
Reproduce via [scripts/check_error_by_year.jl](../../scripts/check_error_by_year.jl).

Replaced `temporal_variability`'s `0.9·mean + 0.1·std` blend plus a `first_non_zero`/
`last_non_zero` double-count that had no justification for weighting the first year (it's
the easiest to fit — zero time to drift). No dispersion (`std`) term: it would fight the
year-weighting itself, which *wants* error to grow with year. Full discussion:
`sandbox/feature_docs/13_break-up-obj-func.md` Q1/Q2 (not version-controlled — this
docstring is the durable record).
"""
function year_weighted_error(
    err_series::AbstractVector{<:Real}, counts::AbstractVector{<:Real}
)::Float64
    n = length(err_series)
    weights = Float64.(1:n) .* sqrt.(counts)
    valid = err_series .!= 0
    return sum(err_series[valid] .* weights[valid]) / sum(weights[valid])
end

function calib_func(sim, obs)
    return temporal_variability(MAEE_series(sim, obs))
end

function _loc_cover(raw_data, dom)
    loc_k_areas = ADRIA.loc_k_area(dom)
    loc_areas = ADRIA.loc_area(dom)
    dropdims(sum(raw_data, dims=(2, 3)), dims=(2, 3)) .* loc_k_areas' ./ loc_areas'
end

"""
    collect_error_stats(raw_data, dom; observations)
    collect_error_stats(raw_data, ltmp_loc_idx, dom; observations)

Error stats for all locations or a single location.
"""
function collect_error_stats(
    raw_data::Array{Float64,4},
    dom;
    observations::LocationDataStore,
)::NamedTuple
    n_locations = size(observations.ltmp_coral_cover, 1)

    # Vectors to hold error stats for each location
    rmse_model, rmse_benchmark, mae, pcc, srcc, bias =
        collect.(eachcol(repeat(zeros(Float64, n_locations), 1, n_locations)))

    for loc in 1:n_locations
        error_stats = collect_error_stats(raw_data, loc, dom; observations=observations)

        rmse_model[loc] = error_stats.rmse_model
        rmse_benchmark[loc] = error_stats.rmse_benchmark
        mae[loc] = error_stats.mae
        pcc[loc] = error_stats.pcc
        srcc[loc] = error_stats.srcc
        bias[loc] = error_stats.bias
    end

    error_names = (:rmse_model, :rmse_benchmark, :mae, :pcc, :srcc, :bias)
    return NamedTuple{error_names}((rmse_model, rmse_benchmark, mae, pcc, srcc, bias))
end
function collect_error_stats(
    raw_data::Array{Float64,4},
    ltmp_loc_idx,
    dom;
    observations::LocationDataStore,
)
    loc_cover = _loc_cover(raw_data, dom)

    obs_loc_data = observations.ltmp_coral_cover[ltmp_loc_idx, :]
    not_missing_obs = (!).(ismissing.(obs_loc_data))

    if !any(not_missing_obs)
        return NaN, NaN, NaN, NaN, NaN
    end

    obs_tf = (START_YEAR:END_YEAR)[not_missing_obs]

    domain_idx::Int64 = ltmp_cover_idx_to_domain(observations, ltmp_loc_idx)
    sim_data = loc_cover[:, domain_idx]

    rmse_::Float64 = rmse(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    pcc_::Float64 = cor(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    mae_::Float64 = MAE(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    bias_::Float64 = bias(sim_data[not_missing_obs], obs_loc_data[not_missing_obs])
    srcc_::Float64 = corspearman(
        sim_data[not_missing_obs], Vector{Float64}(obs_loc_data[not_missing_obs])
    )

    μ_obs = mean(obs_loc_data[not_missing_obs])
    s = length(sim_data[not_missing_obs])
    benchmark_::Float64 = rmse(fill(μ_obs, s), obs_loc_data[not_missing_obs])

    error_names = (:rmse_model, :rmse_benchmark, :mae, :pcc, :srcc, :bias)
    return NamedTuple{error_names}((rmse_, benchmark_, mae_, pcc_, srcc_, bias_))
end

"""
Bray-Curtis dissimilarity between two composition vectors, in [0, 1] (0 = identical, 1 =
fully disjoint).
"""
bray_curtis(sim, obs) = sum(abs.(sim .- obs)) / sum(sim .+ obs)

"""
    reef_taxa_error(cover; observations)::Float64

Mean compositional mismatch between simulated and observed functional group cover across
LTMP reefs and years, via Bray-Curtis dissimilarity. Returns a value in [0, 1], 0 = perfect
match. Reefs with no observations get a neutral 0.5.
"""
function reef_taxa_error(
    cover;
    observations::LocationDataStore
)
    n_locs = length(observations.composition_to_domain)
    loc_scores = Vector{Float64}(undef, n_locs)

    for (j, idx) in enumerate(observations.composition_to_domain)
        non_missing_mask = vec(all(.!ismissing.(observations.coral_composition[:, :, j]); dims=2))
        observed_ids = eachindex(2008:2022)[non_missing_mask]

        if isempty(observed_ids)
            loc_scores[j] = 0.5
            continue
        end

        bc_vals = [
            bray_curtis(cover[id, :, idx], Float64.(observations.coral_composition[id, :, j]))
            for id in observed_ids
        ]

        loc_scores[j] = mean(bc_vals)
    end

    return mean(loc_scores)
end

"""
    reef_error(cover; observations)::Vector{Float64}

Mean reef-level MAE per year, averaged across all LTMP reefs with observations in that year.

For each LTMP reef location, computes the element-wise MAE between simulated `cover`
and observed coral cover at each observed timestep, then accumulates into a `n_years`
error series. Missing observations are excluded from both the error sum and the per-year
reef count, so each element of the returned vector is the unweighted mean over whichever
reefs had valid data that year.

Returns a `Vector{Float64}` of length `n_years` (one entry per modelled timestep).
"""
function reef_error(
    cover;
    observations::LocationDataStore,
)::Vector{Float64}
    n_years::Int64 = size(cover, 1)
    err_series::Vector{Float64} = zeros(Float64, n_years)
    tmp_err::Vector{Float64} = zeros(Float64, n_years)
    err_counts::Vector{Int64} = zeros(Int64, n_years)
    not_missing::BitVector = BitVector(fill(true, n_years))

    domain_idx::Int64 = -1
    each_obs_ltmp_cover_series = eachrow(observations.ltmp_coral_cover)
    for (row_idx, loc_obs) in enumerate(each_obs_ltmp_cover_series)
        domain_idx = ltmp_cover_idx_to_domain(observations, row_idx)

        not_missing .= (!).(ismissing.(loc_obs))
        tmp_err[not_missing] .= MAE_series(
            cover[not_missing, domain_idx], loc_obs[not_missing]
        )

        err_series[not_missing] .+= tmp_err[not_missing]
        err_counts[not_missing] .+= 1.0
    end

    if any(err_counts .== 0)
        @debug "No reef level observation data for some years."
        err_counts[err_counts .== 0] .= 1
    end

    return err_series ./ err_counts
end

function rmse_diff(
    rs_raw::Array{Float64,4}, observations::LocationDataStore, dom
)::Vector{Float64}
    n_test_locs = length(observations.ltmp_cover_to_domain)
    model_rmse = zeros(Float64, n_test_locs)
    benchmark_rmse = zeros(Float64, n_test_locs)
    for i in 1:n_test_locs
        error_stats = collect_error_stats(rs_raw, i, dom; observations=observations)
        model_rmse[i] = error_stats.rmse_model
        benchmark_rmse[i] = error_stats.rmse_benchmark
    end
    return benchmark_rmse .- model_rmse
end

function _rmse_diff_statistic(rows::AbstractMatrix)
    s, o = rows[:, 1], rows[:, 2]
    return rmse(fill(mean(o), length(o)), o) - rmse(s, o)
end

function _bias_statistic(rows::AbstractMatrix)
    s, o = rows[:, 1], rows[:, 2]
    return bias(s, o)
end

function _nse_statistic(rows::AbstractMatrix)
    s, o = rows[:, 1], rows[:, 2]
    rmse_benchmark = rmse(fill(mean(o), length(o)), o)
    return 1 - (rmse(s, o) / rmse_benchmark)^2
end

"""
    _per_reef_bootstrap_stats(rs_raw, observations, dom, statistic; B=2000, rng_seed=1, ci_level=0.95)::NamedTuple

Per-reef bootstrap CI for `statistic(rows)` (a function of the `n_years x 2` matrix of
paired `[sim obs]` rows) — block bootstrap over observed years for reefs with
`n_years >= 5`, iid bootstrap for `n_years == 4`, plus a median + second-level bootstrap
CI aggregated over the block-bootstrap-eligible (`n_years >= 5`) reefs only. Returns a
`NamedTuple` with `estimate`, `ci_lo`, `ci_hi`, `se`, `n_years`, `block_eligible` (one
entry per reef), and `median`/`median_lo`/`median_hi`/`median_se` (aggregate, over
eligible reefs only). Shared by [`rmse_diff_stats`](@ref) and [`nse_stats`](@ref).
"""
function _per_reef_bootstrap_stats(
    rs_raw::Array{Float64,4},
    observations::LocationDataStore,
    dom,
    statistic::Function;
    B::Int=2000,
    rng_seed::Int=1,
    ci_level::Float64=0.95,
)::NamedTuple
    n_locs = length(observations.ltmp_cover_to_domain)
    estimate = zeros(Float64, n_locs)
    ci_lo = zeros(Float64, n_locs)
    ci_hi = zeros(Float64, n_locs)
    se = zeros(Float64, n_locs)
    n_years = zeros(Int, n_locs)
    block_eligible = falses(n_locs)

    loc_cover = _loc_cover(rs_raw, dom)
    Random.seed!(rng_seed)

    for i in 1:n_locs
        obs_loc_data = observations.ltmp_coral_cover[i, :]
        not_missing = (!).(ismissing.(obs_loc_data))

        domain_idx::Int64 = ltmp_cover_idx_to_domain(observations, i)
        sim = loc_cover[not_missing, domain_idx]
        obs = Vector{Float64}(obs_loc_data[not_missing])
        n = length(obs)
        n_years[i] = n
        block = n >= 5
        block_eligible[i] = block

        # Rows kept paired so sim/obs resample together
        data = hcat(sim, obs)
        stats = block ?
                block_bootstrap_ci(statistic, data; B=B, ci_level=ci_level) :
                iid_bootstrap_ci(statistic, data; B=B, ci_level=ci_level)
        estimate[i] = stats.estimate
        ci_lo[i] = stats.lo
        ci_hi[i] = stats.hi
        se[i] = stats.se
    end

    agg = bootstrap_median_ci(estimate[block_eligible]; B=B, ci_level=ci_level)

    return (
        estimate=estimate, ci_lo=ci_lo, ci_hi=ci_hi, se=se, n_years=n_years,
        block_eligible=block_eligible,
        median=agg.median, median_lo=agg.lo, median_hi=agg.hi, median_se=agg.se,
    )
end

"""
    rmse_diff_stats(rs_raw, observations, dom; B=2000, rng_seed=1, ci_level=0.95)::NamedTuple

Per-reef `benchmark_RMSE - model_RMSE` with bootstrap CIs — see
[`_per_reef_bootstrap_stats`](@ref) for the resampling/aggregation method. Returns a
`NamedTuple` with `diff`, `ci_lo`, `ci_hi`, `se`, `n_years`, `block_eligible`, and
`median`/`median_lo`/`median_hi`/`median_se`.
"""
function rmse_diff_stats(
    rs_raw::Array{Float64,4},
    observations::LocationDataStore,
    dom;
    B::Int=2000,
    rng_seed::Int=1,
    ci_level::Float64=0.95,
)::NamedTuple
    stats = _per_reef_bootstrap_stats(
        rs_raw, observations, dom, _rmse_diff_statistic; B=B, rng_seed=rng_seed, ci_level=ci_level
    )
    return (
        diff=stats.estimate, ci_lo=stats.ci_lo, ci_hi=stats.ci_hi, se=stats.se,
        n_years=stats.n_years, block_eligible=stats.block_eligible,
        median=stats.median, median_lo=stats.median_lo, median_hi=stats.median_hi,
        median_se=stats.median_se,
    )
end

"""
    bias_stats(rs_raw, observations, dom; B=2000, rng_seed=1, ci_level=0.95)::NamedTuple

Per-reef bias `mean(modelled - observed)` with bootstrap CIs — see
[`_per_reef_bootstrap_stats`](@ref) for the resampling/aggregation method. Positive values
indicate the model overpredicts coral cover on average, negative values indicate
underprediction. Bootstrap CIs computed identically to [`rmse_diff_stats`](@ref).

In addition to the signed median (shared with the other `*_stats` functions), this also
returns the median *absolute* per-reef bias (`median_abs_bias`, with its own bootstrap CI
via [`bootstrap_median_ci`](@ref)). The signed median can be near zero even when every reef
has a substantial bias, if roughly as many reefs over- as underpredict — the signed
aggregate answers "is there a net directional error across reefs" (relevant for
regional-scale use), while `median_abs_bias` answers "how far off is a typical individual
reef" (relevant since the model is not intended as a precise per-reef predictor). Both are
needed together; the signed median alone can misleadingly read as "the model has almost no
bias." Median (not mean) is used for the absolute-value aggregate too, consistent with every
other aggregate in this module — `|bias|` is right-skewed (a few poorly-fit reefs create a
long tail), so the mean would be inflated by those outliers relative to what's typical;
empirically the mean ran ~20-40% above the median on this project's calibration/test sets.

Returns a `NamedTuple` with `bias`, `ci_lo`, `ci_hi`, `se`, `n_years`, `block_eligible`,
`median`/`median_lo`/`median_hi`/`median_se`, and
`median_abs_bias`/`median_abs_bias_lo`/`median_abs_bias_hi`/`median_abs_bias_se`.
"""
function bias_stats(
    rs_raw::Array{Float64,4},
    observations::LocationDataStore,
    dom;
    B::Int=2000,
    rng_seed::Int=1,
    ci_level::Float64=0.95,
)::NamedTuple
    stats = _per_reef_bootstrap_stats(
        rs_raw, observations, dom, _bias_statistic; B=B, rng_seed=rng_seed, ci_level=ci_level
    )
    abs_agg = bootstrap_median_ci(
        abs.(stats.estimate[stats.block_eligible]); B=B, ci_level=ci_level
    )
    return (
        bias=stats.estimate, ci_lo=stats.ci_lo, ci_hi=stats.ci_hi, se=stats.se,
        n_years=stats.n_years, block_eligible=stats.block_eligible,
        median=stats.median, median_lo=stats.median_lo, median_hi=stats.median_hi,
        median_se=stats.median_se,
        median_abs_bias=abs_agg.median, median_abs_bias_lo=abs_agg.lo,
        median_abs_bias_hi=abs_agg.hi, median_abs_bias_se=abs_agg.se,
    )
end

"""
    nse_stats(rs_raw, observations, dom; B=2000, rng_seed=1, ci_level=0.95)::NamedTuple

Per-reef Nash-Sutcliffe efficiency `1 - (RMSE_model / RMSE_benchmark)^2` — a skill score
normalized by each reef's own observed variance, so reefs with little natural variability
(and thus little room to beat the benchmark in absolute terms) aren't penalized relative
to more volatile reefs. Bootstrap CIs computed identically to
[`rmse_diff_stats`](@ref) — see [`_per_reef_bootstrap_stats`](@ref). Returns a
`NamedTuple` with `nse`, `ci_lo`, `ci_hi`, `se`, `n_years`, `block_eligible`, and
`median`/`median_lo`/`median_hi`/`median_se`.
"""
function nse_stats(
    rs_raw::Array{Float64,4},
    observations::LocationDataStore,
    dom;
    B::Int=2000,
    rng_seed::Int=1,
    ci_level::Float64=0.95,
)::NamedTuple
    stats = _per_reef_bootstrap_stats(
        rs_raw, observations, dom, _nse_statistic; B=B, rng_seed=rng_seed, ci_level=ci_level
    )
    return (
        nse=stats.estimate, ci_lo=stats.ci_lo, ci_hi=stats.ci_hi, se=stats.se,
        n_years=stats.n_years, block_eligible=stats.block_eligible,
        median=stats.median, median_lo=stats.median_lo, median_hi=stats.median_hi,
        median_se=stats.median_se,
    )
end

"""
    correlation_stats(rs_raw, observations, dom; correlation_metric=:spearman, B=2000,
                       rng_seed=1, ci_level=0.95)::NamedTuple

Per-reef correlation coefficient (`:spearman` or `:pearson`) between simulated and
observed cover, with bootstrap CIs — see [`_per_reef_bootstrap_stats`](@ref) for the
resampling/aggregation method (shared with [`rmse_diff_stats`](@ref)/[`nse_stats`](@ref)).
A degenerate resample (e.g. all-identical resampled years) makes the correlation
undefined (zero-variance denominator); rejected and redrawn like any other non-finite
bootstrap replicate, same as NSE's 0/0 case. Returns a `NamedTuple` with `corr`, `ci_lo`,
`ci_hi`, `se`, `n_years`, `block_eligible`, and `median`/`median_lo`/`median_hi`/`median_se`.
"""
function correlation_stats(
    rs_raw::Array{Float64,4},
    observations::LocationDataStore,
    dom;
    correlation_metric::Symbol=:spearman,
    B::Int=2000,
    rng_seed::Int=1,
    ci_level::Float64=0.95,
)::NamedTuple
    correlation_function = if correlation_metric == :spearman
        corspearman
    elseif correlation_metric == :pearson
        cor
    else
        error("Invalid value for `correlation_metric`")
    end

    statistic = rows -> correlation_function(rows[:, 1], rows[:, 2])
    stats = _per_reef_bootstrap_stats(
        rs_raw, observations, dom, statistic; B=B, rng_seed=rng_seed, ci_level=ci_level
    )
    return (
        corr=stats.estimate, ci_lo=stats.ci_lo, ci_hi=stats.ci_hi, se=stats.se,
        n_years=stats.n_years, block_eligible=stats.block_eligible,
        median=stats.median, median_lo=stats.median_lo, median_hi=stats.median_hi,
        median_se=stats.median_se,
    )
end

"""
Use the Fisher transformation to calculate the average correlation coefficient.

Refs:
- https://stats.stackexchange.com/questions/8019/averaging-correlation-values?utm_source=chatgpt.com
- https://www.tandfonline.com/doi/abs/10.1080/00221309809595548
"""
function average_cc(cc_data::Vector{Float64}; w::Vector{Float64}=ones(Float64, length(cc_data)))
    # Apply fisher transformation
    f_data = atanh.(cc_data)

    # Weighted mean
    mean_cc = mean(f_data, weights(w))

    # Apply inverse fisher transformation
    return tanh.(mean_cc)
end
