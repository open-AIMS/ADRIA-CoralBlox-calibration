"""
    end_of_history_dist_mean(final_tol, final_cover, k_area, initial_tol)::Matrix{Float64}

Cover-weighted mean (across locations) heat tolerance at the end of the historical run, as
`(n_groups, n_sizes)`. Forecasts start from this state, not from the initial `dist_mean`.

`final_tol` (`coral_dhw_final` from `ADRIA.run_model`) and `final_cover` (last timestep of
`raw`) are both `(n_groups, n_sizes, n_locs)`. Group/size combinations with no cover anywhere keep their `initial_tol`.
"""
function end_of_history_dist_mean(
    final_tol::AbstractArray{<:Real,3},
    final_cover::AbstractArray{<:Real,3},
    k_area::AbstractVector{<:Real},
    initial_tol::AbstractMatrix{<:Real}
)::Matrix{Float64}
    weights = final_cover .* reshape(k_area, 1, 1, :)
    w_sum = dropdims(sum(weights; dims=3); dims=3)
    mean_tol = dropdims(sum(final_tol .* weights; dims=3); dims=3) ./ w_sum
    return ifelse.(w_sum .> 0.0, mean_tol, Float64.(initial_tol))
end
