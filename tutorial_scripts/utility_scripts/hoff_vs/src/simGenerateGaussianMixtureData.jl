
function fake_mvGausssian_data_for_testing(G,C_t,Ktrue;μ =nothing,τ = nothing, mix_prob =nothing,same_prob_t = true,dynamic = false)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* 10.0 .^(power[k]) .* ones(G)#rand(G)
        end
    end
    if isnothing(τ)
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            x[t][c] = rand(MultivariateNormal(μ[z[t][c]],diagm((τ[z[t][c]]) .^(-1)) ))
        end
    end
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob,μ,τ 
end

function fake_mvGausssian_data_largeG_indepG_for_testing(G,C_t,Ktrue;μ =nothing,τ = nothing, mix_prob =nothing,same_prob_t = true,dynamic = false)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* 100.0 .^(power[k]) .* rand(G)
        end
    end
    if isnothing(τ)
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            x[t][c] = rand(arraydist(Normal.(μ[z[t][c]],(τ[z[t][c]]) .^(-1) )))
        end
    end
        # τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob,μ,τ 
end



function fake_mvGausssian_data_indepG_for_testing(G,C_t,Ktrue;μ =nothing,τ = nothing, mix_prob =nothing,same_prob_t = true,dynamic = false,μ_magnitude = 10.0,largeG = true)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* μ_magnitude .^(power[k]) .* rand(G)
        end
    elseif typeof(μ) <: VecOrMat && eltype(μ) <: Number
        μ = [μ[k] .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(τ)
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    elseif typeof(τ) <: VecOrMat && eltype(τ) <: Number
        τ = [τ[k] .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            if largeG 
                x[t][c] = rand(arraydist(Normal.(μ[z[t][c]],(τ[z[t][c]]) .^(-1) )))
            else
                x[t][c] = rand(MultivariateNormal(μ[z[t][c]],diagm((τ[z[t][c]]) .^(-1)) ))
            end
            
        end
    end
        # τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob,μ,τ 
end
