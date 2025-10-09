@doc raw"""
    t_test(x; conf_level=0.95)
This function calculates the upper and lower confidence interval of a population of parameters using a t-test
"""
function t_test(x; conf_level=0.95)
    alpha = (1 - conf_level)
    tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
    SE = std(x)./sqrt(length(x))

    lo, hi = mean(x) .+ [-1, 1] .* tstar * SE

    return lo, hi
end

@doc raw"""
    norm_weights(p)
This function normalizes a vector of values on the log scale.
```math
\pi_i=\exp(x_i− \text{logsumexp}(x)) 

```
where ``\text{logsumexp}(x)=b+\log \sum_{j=1}^n \exp(x_j-b)`` and ``\pi_i \in [0,1]`` and ``b = \max([x_1,...,x_n])``.
"""
function norm_weights(p::Vector{U})  where {U <: AbstractFloat}
    psum = StatsFuns.logsumexp(p)
    w = exp.(p .- psum)
    return w
end

@doc raw"""
    normToProb(p)
This function normalizes a vector of values
```math
w_i = \frac{x_i}{\sum_{j=1}^n x_j}
```
"""
function normToProb(p::Vector{U}) where {U <: AbstractFloat}
    psum = sum(p)
    w = p ./ psum
    return w
end


@doc raw"""
    norm_weights3(p;float_type=nothing)
This function normalizes a vector of values on the log scale by precallocating an its output.
```math
\pi_i=\exp(x_i− \text{logsumexp}(x)) 

```
where ``\text{logsumexp}(x)=b+\log \sum_{j=1}^n \exp(x_j-b)`` and ``\pi_i \in [0,1]`` and ``b = \max([x_1,...,x_n])``.
"""
function norm_weights3(p::Vector{U};float_type=nothing)  where {U <: AbstractFloat}
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

@doc raw"""
    norm_weights3!(p;float_type=nothing)
This function normalizes a vector of values on the log scale by precallocating an its output and performs operations in place.
```math
\pi_i=\exp(x_i− \text{logsumexp}(x)) 

```
where ``\text{logsumexp}(x)=b+\log \sum_{j=1}^n \exp(x_j-b)`` and ``\pi_i \in [0,1]`` and ``b = \max([x_1,...,x_n])``.
"""
function norm_weights3!(p::Vector{U};float_type=nothing)  where {U <: AbstractFloat}
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

@doc raw"""
    norm_weights3!(K,p;float_type=nothing)
This function normalizes a vector of values on the log scale by precallocating an its output and performs operations in place. Normalization only occurs up until the Kth element in the vector
```math
\pi_i=\exp(x_i− \text{logsumexp}(x)) 

```
where ``\text{logsumexp}(x)=b+\log \sum_{j=1}^n \exp(x_j-b)`` and ``\pi_i \in [0,1]`` and ``i \in \{1,..,K\}`` and ``b = \max([x_1,...,x_n])``
"""
function norm_weights3!(K::Int,p::Vector{U};float_type=nothing) where {U <: AbstractFloat}
    if isnothing(float_type)
        float_type =eltype(p)
    end
    psum = @views convert(float_type,StatsFuns.logsumexp(p[1:K]))
    for k in 1:K
        p[k] = exp(p[k] - psum)
    end
    
    return p
end

@doc raw"""
    normToProb3!(p;float_type=nothing)
This function normalizes a vector of values by precallocating an its output and performs operations in place.
```math
w_i = \frac{x_i}{\sum_{j=1}^n x_j}
```
"""
function normToProb3!(p::Vector{U};float_type=nothing) where {U <: AbstractFloat}
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

@doc raw"""
    sigmoidNorm!(p;float_type=nothing)
This function normalizes a vector of values using a logistic function by precallocating an its output and performs operations in place.
```math
f(x_i) = \frac{1}{1 + e^{-x_i}}
```
"""
function sigmoidNorm!(p::Vector{U};float_type=nothing)  where {U <: AbstractFloat}
    K = length(p)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    for k in 1:K
        p[k] =  StatsFuns.logistic(p[k])
    end
    return p
end

@doc raw"""
    sigmoidNorm!(K,p;float_type=nothing)
This function normalizes a vector of values using a logistic function by precallocating an its output and performs operations in place. Normalization only occurs up until the Kth element in the vector
```math
f(x_i) = \frac{1}{1 + e^{-x_i}} \text{ for } i \in \{1,..,K\}
```
"""
function sigmoidNorm!(K::Int,p::Vector{U};float_type=nothing)  where {U <: AbstractFloat}
    if isnothing(float_type)
        float_type =eltype(p)
    end
    for k in 1:K
        p[k] =  StatsFuns.logistic(p[k])
    end
    return p
end

@doc raw"""
    normToProb3!(K,p;float_type=nothing)
This function normalizes a vector of values by precallocating an its output and performs operations in place. Normalization only occurs up until the Kth element in the vector
```math
w_i = \frac{x_i}{\sum_{j=1}^K x_j} \text{ for } i \in \{1,..,K\}
```
"""
function normToProb3!(K::Int,p::Vector{U};float_type=nothing)  where {U <: AbstractFloat}
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


##########################################################
######## SIMPLE DISTRIBUTION NORMALIZER FUNCTIONS ########
##########################################################
@doc raw"""
    ln_Gamma_distribution_normalizer(a::AbstractFloat,b::AbstractFloat)
This function computes the log normalizer of a Gamma distribution
```math
\ln(Z) = a*\log(b) - \log\Gamma(a)
```
where ``\Gamma(a)`` is the Gamma function
"""
function ln_Gamma_distribution_normalizer(a::AbstractFloat,b::AbstractFloat)
    return a*log(b) - loggamma(a)
end

@doc raw"""
    ln_Beta_distribution_normalizer(a::AbstractFloat,b::AbstractFloat)
This function computes the log normalizer of a Beta distribution

```math
\ln(Z) = \log \text{B}(a,b)
```
where ``\text{B}(a,b)`` is the Beta function
"""
function ln_Beta_distribution_normalizer(a::AbstractFloat,b::AbstractFloat)
    return logbeta(a,b)
end

@doc raw"""
    ln_Dirichlet_distribution_nomralizer(a::Vector{U}) where {U <: AbstractFloat}
This function computes the log normalizer of a Dirichlet distribution
```math
\ln(Z) = \sum_{i=1}^K \log\Gamma(a_i) - \log\Gamma(\sum_{i=1}^K a_i)
```
where ``\Gamma(a)`` is the Gamma function
"""
function ln_Dirichlet_distribution_nomralizer(a::Vector{U}) where {U <: AbstractFloat}
    return sum(loggamma.(a)) - loggamma(sum(a))
end

############################################
############################################
############################################