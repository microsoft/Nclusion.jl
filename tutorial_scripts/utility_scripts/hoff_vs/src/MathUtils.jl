module MathUtils

using Distributions
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics

export  StickBreaking, jitter!, calc_r, calc_p, calc_dispersion, t_test,norm_weights,normToProb,maxk!,norm_weights3,norm_weights3!, normToProb3!,log_norm_pdf3,vectorfloor!,sigmoidNorm!,normalizelogweights

function StickBreaking(β_prime)
    one_minus_β_prime = 1. .- β_prime
    cumprod_one_minus_β_prime = cumprod(one_minus_β_prime, dims = 1)
    β_prime_copy = deepcopy(β_prime)
    push!(β_prime_copy, 1.0)
    pushfirst!(cumprod_one_minus_β_prime, 1.0)
    β = β_prime_copy .* cumprod_one_minus_β_prime
    return β
end

function jitter!(a::Array, factor=1.0)
    @assert eltype(a) <: AbstractFloat
    a .+= rand(size(a)...) .* factor
end

function calc_r(μ,σ²)
    r = μ.^2 ./ (σ² .- μ)
    return r
end

function calc_p(μ,σ²)
    r = calc_r(μ,σ²)
    p = r ./ (r .+ μ)
    return p
end
function calc_dispersion(μ,a_)
    # r = calc_r(μ,σ²)
    σ² = μ .+ a_ .* μ .^ 2 
    return σ²
end

function t_test(x; conf_level=0.95)
    alpha = (1 - conf_level)
    tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
    SE = std(x)./sqrt(length(x))

    lo, hi = mean(x) .+ [-1, 1] .* tstar * SE
    # "($lo, $hi)"
    return lo, hi
end

function maxk!(ix, a, k; initialized=false)
    partialsortperm!(ix, a, 1:k, rev=true, initialized=initialized)
    @views collect(zip(ix[1:k], a[ix[1:k]]))
end

function norm_weights(p)
    psum = StatsFuns.logsumexp(p)
    w = exp.(p .- psum)
    return w
end

function normToProb(p)
    psum = sum(p)
    w = p ./ psum
    return w
end

function normalizelogweights(p)
    pmax =  maximum(p)
    w = exp.(p .- pmax)
    return w ./sum(w)
end

function norm_weights3(p;float_type=nothing)
    K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    psum = convert(float_type,StatsFuns.logsumexp(p))
    w = Vector{float_type}(undef,K)
    for k in 1:K
        w[k] = exp(p[k] - psum)
    end
    
    return w
end
function norm_weights3!(p;float_type=nothing)
    K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    psum = convert(float_type,StatsFuns.logsumexp(p))
    # w = Vector{float_type}(undef,K)
    for k in 1:K
        p[k] = exp(p[k] - psum)
    end
    
    return p
end

function norm_weights3!(K,p;float_type=nothing)
    # K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    psum = @views convert(float_type,StatsFuns.logsumexp(p[1:K]))
    # w = Vector{float_type}(undef,K)
    for k in 1:K
        p[k] = exp(p[k] - psum)
    end
    
    return p
end
function normToProb3!(p;float_type=nothing)
    K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    psum = convert(float_type,sum(p))
    for k in 1:K
        p[k] =  p[k] / psum
    end
    return p
end
function sigmoidNorm!(p;float_type=nothing)
    K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    for k in 1:K
        p[k] =  StatsFuns.logistic(p[k])
    end
    return p
end
function sigmoidNorm!(K,p;float_type=nothing)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    for k in 1:K
        p[k] =  StatsFuns.logistic(p[k])
    end
    return p
end
function normToProb3!(K,p;float_type=nothing)
    # K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    psum = @views convert(float_type,sum(p[1:K]))
    for k in 1:K
        p[k] =  p[k] / psum
    end
    return p
end

function log_norm_pdf3(x,μk,a0k,b0k)
    float_type = typeof(x)
    logpi = convert(float_type,log(2π))
    return 1/2 * log(a0k/b0k) -  1/2 * ((x-μk)^2 * (a0k/b0k)) - 1/2 * logpi
end
function vectorfloor!(vec_vals)
    for i in eachindex(vec_vals)
        vec_vals[i] = floor(vec_vals[i])
    end
    return vec_vals
end


end