"""
    dist_mean_order_penalty(
        group_means::AbstractVector{<:Real}; margin::Real=0.0, scale::Real=DIST_MEAN_ORDER_SCALE
    )::Float64

Dimensionless penalty for initial mean heat tolerances that break the expected order:
tabular Acropora lowest, corymbose Acropora and corymbose non-Acropora in the middle,
small and large massives highest (Hughes et al. 2018, Extended Data Fig. 5). The order within
the middle pair and within the massives is free. Returns 0.0 when the order holds.

`group_means` must be in `ADRIA.functional_group_names()` order. `margin` (DHW-weeks) is the
minimum required gap between tiers, and `scale` (DHW-weeks) normalises the violation.
"""
function dist_mean_order_penalty(
    group_means::AbstractVector{<:Real};
    margin::Real=0.0,
    scale::Real=DIST_MEAN_ORDER_SCALE
)::Float64
    length(group_means) == N_TAXA ||
        throw(ArgumentError("Expected $N_TAXA group means, got $(length(group_means))"))
    ta, ca, cna, sm, lm = group_means
    violation = max(0.0, ta - min(ca, cna) + margin) +
                max(0.0, max(ca, cna) - min(sm, lm) + margin)
    return Float64(violation / scale)
end
