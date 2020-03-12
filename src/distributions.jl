export dnorm, dgumbel, dgamma, dlognorm

"""
    dnorm(x::Real, μ::Real, σ::Real)

Normal distribution density function
"""
dnorm(x::Real, μ::Real, σ::Real) = exp(-0.5 * ((x - μ)/σ)^2)/(σ*sqrt2π)

"""
    dgumbel(x::Real, μ::Real, σ::Real)

Gumbel density, parameterized by mean (μ) and SD (σ) of x.
"""
function dgumbel(x::Real, μ::Real, σ::Real)
    β = (√6 * σ) / π
    α = μ - (β*γ)
    z = (x - α) / β
    exp(-(z + exp(-z)))
end

"""
    dgamma(x::Real, μ::Real, σ::Real)

Gamma density, parameterized by mean (μ) and sigma(σ) of x.
"""
function dgamma(x::Real, μ::Real, σ::Real)
    α = μ^2 / σ^2
    β = μ / σ^2
    (x^(α - 1) * exp(-β * x))  * (β^α / Γ(α))
end

"""
    dlognorm(x::Real, μ::Real, σ::Real)

Lognormal density, parameterized by mean (μ) and SD (σ) of x, NOT of log(x).
"""
function dlognorm(x::Real, μ::Real, σ::Real)
    lμ = log(μ^2 / sqrt(1 + σ^2 / μ^2))
    lσ = √(log(1 + σ^2 / μ^2))
    dnorm(log(x), lμ, lσ) / x
end
