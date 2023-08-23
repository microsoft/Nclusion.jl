"""
    t_test(x; conf_level=0.95)
This function calculates the upper and lower confidence interval of a population of parameters using a t-test
"""
function t_test(x; conf_level=0.95)
    alpha = (1 - conf_level)
    tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
    SE = std(x)./sqrt(length(x))

    lo, hi = mean(x) .+ [-1, 1] .* tstar * SE
    # "($lo, $hi)"
    return lo, hi
end

"""
    norm_weights(p)
This function normalizes a vector of values on the log scale.
```math
π_i=exp(x_i− logsumexp(x)) where logsumexp(x)=b+log∑_{j=1}^n exp(x_j−b)
π_i in [0,1]
```
"""
function norm_weights(p)
    psum = StatsFuns.logsumexp(p)
    w = exp.(p .- psum)
    return w
end

"""
    normToProb(p)
This function normalizes a vector of values
```math
w_i = frac{x_i}{∑_{j=1}^n x_j}
```
"""
function normToProb(p)
    psum = sum(p)
    w = p ./ psum
    return w
end


"""
    norm_weights3(p;float_type=nothing)
This function normalizes a vector of values on the log scale by precallocating an its output.
```math
π_i=exp(x_i− logsumexp(x)) where logsumexp(x)=b+log∑_{j=1}^n exp(x_j−b)
π_i in [0,1]
```
"""
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

"""
    norm_weights3!(p;float_type=nothing)
This function normalizes a vector of values on the log scale by precallocating an its output and performs operations in place.
```math
π_i=exp(x_i− logsumexp(x)) where logsumexp(x)=b+log∑_{j=1}^n exp(x_j−b)
π_i in [0,1]
```
"""
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

"""
    norm_weights3!(K,p;float_type=nothing)
This function normalizes a vector of values on the log scale by precallocating an its output and performs operations in place. Normalization only occurs up until the Kth element in the vector
```math
π_i=exp(x_i− logsumexp(x)) text{ where } logsumexp(x)=b+log∑_{j=1}^K exp(x_j−b)
text{ for } π_i in [0,1] text{and} i in {1,..,K}
```
"""
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

"""
    normToProb3!(p;float_type=nothing)
This function normalizes a vector of values by precallocating an its output and performs operations in place.
```math
w_i = frac{x_i}{∑_{j=1}^n x_j}
```
"""
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

"""
    sigmoidNorm!(p;float_type=nothing)
This function normalizes a vector of values using a logistic function by precallocating an its output and performs operations in place.
```math
f(x_i) = frac{1}{1 + e^{-x_i}}
```
"""
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

"""
    sigmoidNorm!(K,p;float_type=nothing)
This function normalizes a vector of values using a logistic function by precallocating an its output and performs operations in place. Normalization only occurs up until the Kth element in the vector
```math
f(x_i) = frac{1}{1 + e^{-x_i}} text{ for } i in {1,..,K}
```
"""
function sigmoidNorm!(K,p;float_type=nothing)
    if isnothing(float_type)
        float_type =eltype(p)
    end
    for k in 1:K
        p[k] =  StatsFuns.logistic(p[k])
    end
    return p
end

"""
    normToProb3!(K,p;float_type=nothing)
This function normalizes a vector of values by precallocating an its output and performs operations in place. Normalization only occurs up until the Kth element in the vector
```math
w_i = frac{x_i}{∑_{j=1}^K x_j} text{ for } i in {1,..,K} 
```
"""
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


