export em_step!, em_run!, mixfit, densfit

"""
    em_step!(est::MixModel, x::Vector{<:Real})

Run a single EM step.

See also: [`em_run!`](@em_run!), [`mixfit`](@mixfit)
"""
function em_step!(est::MixModel, x::Vector{<:Real})
    p = Matrix{Float32}(undef, length(est.μ), length(x))
    for i = 1:length(est.μ)
        p[i,:] = est.α[i] * dnorm.(x, est.μ[i], est.σ[i]) ./
            sum(est.α[j] * dnorm.(x, est.μ[j], est.σ[j]) for j=1:length(est.μ))
    end
    for i = 1:length(est.μ)
        est.μ[i] = sum(p[i,:] .* x) / sum(p[i,:])
        est.σ[i] = √(sum(p[i,:] .* (x .- est.μ[i,:]).^2) / sum(p[i,:]))
        est.α[i] = mean(p[i,:])
    end
end

"""
    em_run!(est::MixModel, x::Vector{<:Real}, rtol::AbstractFloat = 0.00001, maxiter::Int = 0)

Run EM steps until the percent increase in log-likelihood is below `rtol` or the number of iterations is greater than `maxiter`.
If `maxiter` is set to zero, EM steps will continue until the increase is below `rtol` regardless of the number of iterations.

See also: [`em_run!`](@em_run!), [`mixfit`](@mixfit)
"""
function em_run!(est::MixModel, x::Vector{<:Real}, rtol::AbstractFloat = 0.00001, maxiter::Int = 0)
    ll::Float32 = LL(x, est)
    oll::Float32 = -ll
    i::Int = 1
    while true
        em_step!(est, x)
        oll = ll
        ll = LL(x, est)
        if abs((oll - ll) / oll) < rtol
            break
        end
        if i == maxiter
            @warn "Hit maximum number of iterations"
        end
        i += 1
    end
end

"""
    mixfit(x::Vector{<:Real},
            m::Int;
            rtol::AbstractFloat = 0.00001,
            α::Vector{<:Real} = fill(1/m, m),
            μ::Vector{<:Real} = quantile!(x, (1:m)/m),
            σ::Vector{<:Real} = fill(std(x) / √(m), m),
            maxswap::Int = 5 * m^2,
            maxiter_inner::Int = 0,
            maxiter::Int = 0,
            silent::Bool = false,
            kernel::Function = dnorm)

Get the maximum likelihood estimate of an `m`-component mixture model using random-swap EM.
Random-swap EM avoids local optimums by randomly replacing components and
using the result with the maximum likelihood [^1].  The component distributions
are given by `kernel`.  The maximum number of swaps is given by `maxswap`.
If results vary across runs, then `maxswap` is too low - the default is ``5m^2``.
maxiter_inner is the maximum number of iterations for the estimates that go through
the swapping process, and maxiter is the maximum for the final estimate.  Similarly,
`rtol` is for the final estimate, while the inner estimates use 0.1.  To supress output,
simply set `silent` to true.  Starting values can be provided via the α, μ, and σ
arguments, but this shouldn't be necassary due to the use of random swapping.

[^1]: Zhao, Q., Hautamäki, V., Kärkkäinen, I., & Fränti, P. (2012). Random swap EM algorithm for Gaussian mixture models. Pattern Recognition Letters, 33(16), 2120-2126.

See also: [`densfit`](@densfit), [`em_run!`](@em_run!)
"""
function mixfit(x::Vector{<:Real},
                m::Int;
                rtol::AbstractFloat = 0.00001,
                t::Int = 0,
                α::Union{Vector{<:Real}, Nothing} = nothing,
                μ::Union{Vector{<:Real}, Nothing} = nothing,
                σ::Union{Vector{<:Real}, Nothing} = nothing,
                maxswap::Int = 0,
                maxiter_inner::Int = 0,
                maxiter::Int = 0,
                silent::Bool = false,
                kernel::Function = dnorm)
    if maxswap == 0
        maxswap = 5 * m^2
    end

    est = MixModel(fill(1/m, m), quantile!(x, (1:m)/m), fill(std(x) / √(m), m), kernel)

    if μ != nothing
        est.μ = μ
    end
    if σ != nothing
        est.σ = σ
    end
    if α != nothing
        est.α = α
    end

    em_run!(est, x, 0.1, maxiter_inner)
    for i ∈ 1:maxswap
        est_s = copy(est)
        choice = rand(1:m)
        est_s.μ[choice] = rand(x)
        est_s.α = fill(1/m, m)
        est_s.σ[choice] = std(x) / √(m)
        em_run!(est_s, x, 0.1, maxiter_inner)
        if LL(x, est_s) > LL(x, est)
            est = est_s
        end
    end
    em_run!(est, x, rtol, maxiter)

    if silent == false
        describe(est, data = x)
    end
    return est
end

"""
    densfit(x::Vector{<:Real};
            wait::Int = 3,
            rtol_em::AbstractFloat = 0.00001,
            criterion::Function = AIC,
            silent::Bool = false,
            maxiter::Int = 0,
            maxiter_inner::Int = 0,
            maxswap::Int = 0,
            kernel::Function = dnorm)

Estimate the density of `x` by a mixture model.  Successive mixture model estimates
are done via random swap EM with increasing number of clusters until `criterion`
decreases for 3 iterations.  The relative tolerance for EM convergence is given
by `rtol_em`.  To disable output, set `silent` to true.

See also: [`densfit`](@densfit), [`em_run!`](@em_run!)
"""

function densfit(x::Vector{<:Real};
                wait::Int = 3,
                rtol_em::AbstractFloat = 0.00001,
                criterion::Function = AIC,
                silent::Bool = false,
                maxiter::Int = 0,
                maxiter_inner::Int = 0,
                maxswap::Int = 0,
                kernel::Function = dnorm)
    new::MixModel = mixfit(x, 1, silent = true)
    current::MixModel = new
    m::Int = 2
    nbad::Int = 0

    while nbad < wait
        new = mixfit(x, m,
                    silent = true,
                    maxiter = maxiter,
                    maxiter_inner = maxiter_inner,
                    maxswap = maxswap,
                    kernel = kernel)

        if (criterion(x, new) >= criterion(x, current))
            nbad += 1
        else
            nbad = 0
            current = new
        end

        m += 1
    end

    if silent == false
        describe(current, data = x)
    end

    return current
end
