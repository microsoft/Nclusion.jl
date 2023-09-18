@model function timeseries_indep_dp_pmm1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # basically the same as dp_pmm_sb2 but now we can change Hyperparameters a_θ,b_θ ,a,b
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    # v =  rand(arraydist(filldist.(StickBreakingProcess.(crm), 12 - 1)))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))
    # z = rand.(filldist.(Categorical.(π_), nobs_t))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # for t in 1:T
    #     x[t] ~ arraydist(Poisson.(λ_z[t]))
    # end 

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
            else
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
            end
        end
        # k[t] = k_t
    end
end

@model function timeseries_indep_dp_pmm2(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # Faster because of less loop and type definitions
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)
    z = Vector{Vector}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
        λ_z = λ[z[t]]
        x[t] ~ arraydist(Poisson.(λ_z))
    end 
end

@model function timeseries_indep_dp_pmm2FactorialDesign(x, T, K,factorLevels;a_θ = 1,b_θ = 1,a=6,b = 1) # Faster because of less loop and type definitions
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    if factorLevels["θt_constant"]
        θ_t ~ arraydist(Dirac.(ones(T)))
    else
        θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    end

    crm = DirichletProcess.(θ_t)



    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T) 
    π_ = Vector{Vector}(undef, T)
    if factorLevels["likelihood_form"] == "no filldist"
        z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    else
        z = Vector{Vector}(undef, T)#
    end
   

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        if factorLevels["likelihood_form"] == "MixtureModel likelihood"
            x[t] .~ MixtureModel(Poisson.(λ), π_[t])
        else
            if factorLevels["likelihood_form"] == "no filldist"
                z[t] .~ Categorical(π_[t])
            else
                z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])#
            end
            
            # λ_z = λ[z[t]]
            if factorLevels["likelihood_form"] == "no arraydist"
                λ_z = λ[z[t]]
                x[t] .~ Poisson.(λ_z)
            else
                if factorLevels["likelihood_form"] == "direct index of λ with z "
                    x[t] ~ arraydist(Poisson.(λ[z[t]]))
                else
                    λ_z = λ[z[t]]
                    x[t] ~ arraydist(Poisson.(λ_z))
                end
                
            end
            # x[t] ~ arraydist(Poisson.(λ_z))
        end

        
    end 
end
@model function timeseries_indep_dp_pmm3(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # Work in Progress; Tryna make it faster
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    # crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)
    z = Vector{Vector}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))

    for t in 1:T
        v[t] ~ filldist(Beta(1,θ_t[t]), K - 1)#filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
        λ_z = λ[z[t]]
        x[t] ~ arraydist(Poisson.(λ_z))
    end 
end
@model function timeseries_indep_dp_pmm_MixtureModelLiklihood1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        x[t] .~ MixtureModel(Poisson.(λ), π_[t])
    end 
end

@model function timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood1(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    
    #OUTPUT IS A MATRIX
    λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

    H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        x[t] .~ MixtureModel(H, π_[t])
    end 
end
@model function timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood2(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    
    #OUTPUT IS A MATRIX
    # λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)
    # H = Product.(map(j -> Poisson.(λ)[:,j],1:K))

    λ = Vector(map(k -> tzeros(Float64, G),1:K))
    λ .~ map(k-> filldist(Gamma(ak[k],1/bk[k]),G)  ,1:K)
    H = Product.(map(j -> Poisson.(λ[j]),1:K))


    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        x[t] .~ MixtureModel(H, π_[t])
    end 
end
@model function timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood22(x, T, K,G,a,b_ak,b_bk;a_θ = 1,b_θ = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))

    # filldist.(Gamma.(a, 1.0 ./b), G)
    
    a=a
    b = TArray(Float64,K)#Vector{Float64}(undef,K)
    # λ =  TArray(TArray,K)
    # H = Vector{Distribution}(undef,K)
    println("DUM BDU")
    # println(λ)
    for k in 1:K
        b[k] ~  Gamma(b_ak[k], b_bk[k])
        # λ[k] ~  filldist(Gamma(a, 1.0 ./ b[k]), G) #Vector{Vector}(undef,K)
    end
    println("1")
    λ ~ arraydist(filldist.(Gamma.(1,b),G))
    println(λ)
    H =  Product.(map(k -> Poisson.(λ[:,k]),1:K))
    println("First Class")
    # b=b
    
    #OUTPUT IS A MATRIX
    # println("here1")

    # H = Product.(map(j -> Poisson.(λ[j]),1:K))

    #filldist(Gamma(a,1/b), K)
    # println("here2")S
    
    # println("here3")
    #  Vector(map(k -> filldist(Gamma(a[k],1/b[k]),G),1:K)) #Vector(map(k -> Vector{Float64}(undef,G),1:K))  #  Vector(map(k -> TArray(Float64,G),1:K))#
    # cluster_gamma_prior = Gamma.(a, 1 ./ b)
    # # λ .~ vec(map(k -> filldist(Gamma(a[k],1/b[k]),G),1:K))
    # # λ .~ map(k-> filldist(Gamma(a[k],1/b[k]),G)  ,1:K)
    # # # println(typeof(λ))
    # for k in 1:K
    #     λ[k] ~ filldist(cluster_gamma_prior[k],G)
    # end
    
    
    
    # tmpPrior =  Vector(map(dist-> filldist(dist,G), Gamma.(a, 1 ./ b)))
    # println(H)
    # # println(size(tmpPrior))
    # # println(typeof(tmpPrior))
    # λ .~ Vector(map(dist-> filldist(dist,G), Gamma.(a, 1 ./ b)))
    # println("here4")
    # H = Product.(map(j -> Poisson.(λ[j]),1:K))
    # println(typeof(H))
    # println("Here5")
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] = stickbreak(v[t])
        # println("Here6")
        x[t] .~ MixtureModel(H, π_[t])
    end 
end


@model function timeseries_hdp_pmm_MultivariateMixtureModelLiklihood1(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1, γ = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    γ = γ
    vβ ~  filldist(Beta(1,γ), K - 1)
    β = stickbreak(vβ)
    #OUTPUT IS A MATRIX
    λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

    H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        c =  θ_t[t] .* β
        d =  θ_t[t] .*( 1 .-  cumsum(β))
        c = c[1:end-1]
        d = d[1:end-1]
        # v[t] = tzeros(Float64, K-1)
        v[t] ~ arraydist(Beta.(c,d))
        π_[t] = stickbreak(v[t])
        x[t] .~ MixtureModel(H, π_[t])
    end 
end
@model function timeseries_hdp_pmm_MultivariateMixtureModelLiklihood2(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1, γ = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    γ = γ
    vβ ~  filldist(Beta(1,γ), K - 1)
    β = stickbreak(vβ)
    #OUTPUT IS A MATRIX
    λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

    H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        π_[t] ~  Dirichlet(θ_t[t] .* β)
        x[t] .~ MixtureModel(H, π_[t])
    end 
end
@model function timeseries_hdp_pmm_MultivariateMixtureModelLiklihood3(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1, γ = 1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    # crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    γ = γ
    β ~  Dirichlet(γ .* ones(K) ./ K)
    #OUTPUT IS A MATRIX
    λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

    H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        π_[t] ~  Dirichlet(θ_t[t] .* β)
        x[t] .~ MixtureModel(H, π_[t])
    end 
end

@model function timeseries_hdp_gmm_MultivariateMixtureModelLiklihood1(x, T, K,G;a_θ = 1,b_θ = 1,a=1,b = 1, γ = 1, λ0 = 1, μ0=1) # s
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)



    a=a
    b=b
    γ = γ
    vβ ~  filldist(Beta(1,γ), K - 1)
    β = stickbreak(vβ)
    
    τk ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)
    μk ~ arraydist(Normal.(μ0,(λ0 .* τk) .^(-1.)))
    H = map(j -> MultivariateNormal(μk[:,j],diagm(τk[:,j] .^(-1.) ) ),1:K)

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)

    for t in 1:T
        π_[t] ~  Dirichlet(θ_t[t] .* β)
        x[t] .~ MixtureModel(H, π_[t])
    end 
end







@model function timeseries_indep_dp_gmm1(x, T, K;a_θ = 1,b_θ = 1) 
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)

    s² ~ filldist(InverseGamma(2, 3), K)
    m ~ arraydist(Normal.(0, sqrt.(s²)))
    H =  Normal.(m, sqrt.(s²))

    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ H[z[t][c]]#MixtureModel(Poisson,λ, π_[t])
            else
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ H[z[t][c]]#MixtureModel(Poisson,λ, π_[t])
            end
        end
        # k[t] = k_t
    end
end





# One time Point

@model function _indep_dp_gmm1(x, K;a_θ = 1,b_θ = 1,a=6,b = 1) 
    nobs = length(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ Gamma(a_θ,1/b_θ) 
    crm = DirichletProcess(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=a
    b=b
    s² ~ filldist(InverseGamma(2, 3), K)
    sqrt_s = sqrt.(s²)
    m ~ arraydist(Normal.(0, sqrt.(s²)))
    H =  Normal.(m, sqrt.(s²))
    
    v ~ filldist(StickBreakingProcess(crm), K - 1)
    π_ =  stickbreak(v)
    z = tzeros(Int, nobs)
    for c in 1:nobs
        z[c] ~ Categorical(π_)
        x[c] ~  H[z[c]]#Normal(m[z[c]], sqrt.(s²)[z[c]])#MixtureModel(Poisson,λ, π_[t])
    end
    return x,z
end
@model function _indep_dp_pmm1(x, K;a_θ = 1,b_θ = 1,a=6,b = 1) 
    nobs = length(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ Gamma(a_θ,1/b_θ) 
    crm = DirichletProcess(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    v ~ filldist(StickBreakingProcess(crm), K - 1)
    π_ =  stickbreak(v)
    z = tzeros(Int, nobs)
    for c in 1:nobs
        z[c] ~ Categorical(π_)
        x[c] ~ Poisson(λ[z[c]])#MixtureModel(Poisson,λ, π_[t])
    end
    return x,z
end


@model function timeseries_indep_dp_pmm1_fixedλ(x, T, K,true_λ;a_θ = 1,b_θ = 1,a=6,b = 1) # basically the same as dp_pmm_sb2 but now we can change Hyperparameters a_θ,b_θ ,a,b
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=a
    b=b
    λ  ~ arraydist(Dirac.(true_λ))
    v = Vector{Vector}(undef, T)
    π_ = Vector{Vector}(undef, T)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
            else
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
            end
        end
        # k[t] = k_t
    end
end

@model function infiniteTimeSeriesPMM(x,T)
    # Hyper-parameters, i.e. concentration parameter and parameters of H.

    # α = 1.0
    # μ0 = 0.0
    # σ0 = 1.0
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)

    a=6
    b=1
    
    


    rpm = DirichletProcess.(θ_t)
    
    # Define the base distribution, i.e. expected value of the Dirichlet process.
    H = Gamma(a,1/b)
    
    # Latent assignment.
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
        
    # Locations of the infinitely many clusters.
    λ = tzeros(Float64,0)
    nk = Vector{Vector{Int}}(undef, T)
    for t in 1:T
        for c in 1:nobs_t[t]
            # Number of clusters.
            K = maximum(vcat(z[1:t]...))
            nk[t] = Vector{Int}(map(k -> sum(z[t] .== k), 1:K))

            # Draw the latent assignment.
            if t==1
                n__ = nk[t]#z[t][c] ~ ChineseRestaurantProcess(rpm[t], nk[t])
            else
                n__ =  nk[t] .+ nk[t-1]#z[t][c] ~ ChineseRestaurantProcess(rpm[t], nk[t] .+ nk[t-1])
            end
            z[t][c] ~ ChineseRestaurantProcess(rpm[t], n__)
            
        
            # Create a new cluster?
            if z[t][c] > K
                if t != 1
                    for i in 1:t-1
                        push!(nk[i],zero(eltype(nk[i])))
                    end
                end
                push!(λ, 0.0)

                # Draw location of new cluster.
                # println("Here at $K")
                λ[z[t][c]] ~ H
            end
                    
            # Draw observation.
            x[t][c] ~ Poisson(λ[z[t][c]])
        end
    end
end

# m3 = infiniteTimeSeriesPMM(x,T)
# iterations = 1000
# @time smc_chain = begin
#     burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
#     n_samples = 1000
#     iterations = burn + n_samples
 
#     chain = sample(m3, SMC(), iterations);
# end

@model dp_pmm_sb2(x, T, K) = begin
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=6
    b=1
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, K-1)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                z[t][c] ~ Categorical(eta)
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                z[t][c] ~ Categorical(eta)
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
            end
        end
        # k[t] = k_t
    end

    
    
end
# burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
# n_samples = 500
# iterations = burn + n_samples
# K = 25
# T = 4
# m2=dp_pmm_sb2(x, T, K)
# # pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))
# tchain = sample(m2, PG(100),iterations);
# k = map(
#     t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :k), :].value))),
#     1:iterations
# );
# ids = findall(map(name -> occursin("λ", string(name)), names(tchain)));
# p = plot(tchain[:, ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter)


# DP PMM model under stick-breaking construction
@model dp_pmm_sb1(x, T, K) = begin
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=6
    b=1
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, K-1)
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            end
        end
    end

    
    
end

# x = [[data_dict[(i,t)][1] for i in 1:6] for t in 1:4]
# true_clusters = [[truth_dict[(i,t)] for i in 1:6] for t in 1:4]
# coordinates = [[(i,t) for i in 1:6] for t in 1:4]
# K = 5
# T = 4
# m=dp_pmm_sb1(x, T, K)

# # Set random seed for reproducibility
# Random.seed!(0);

# # Compile time approx. 32s.
# # Run time approx. 70s.

# @time hmc_chain = begin
#     burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
#     n_samples = 500
#     iterations = burn + n_samples
#     n_components = 10
#     stepsize = 0.01
#     nleapfrog = floor(Int, 1 / stepsize)
 
#     chain = sample(dp_pmm_sb1(x, T, K), 
#            HMC(stepsize, nleapfrog),
#            iterations)
# end
# ;


#Hierarchical 
#Work in Progress

@model function timeseries_Hdpsb_dp_pmm1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1)
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    # crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))
    γ ~ Gamma(1,1)
    # vβ ~ filldist(Beta(1,γ), K - 1)
    β ~ Dirichlet(γ .* ones(Float64, K) )#stickbreak(vβ)
    # minus_β_cumsum = 1 .- cumsum(β)
    # vπ = Vector{Vector}(undef, T)
    π_ = Vector{Vector{Float64}}(undef, T)
    z = Vector{Vector{Int64}}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    # rand.(Beta.(alpha[1] .* beta,alpha[1] .* minus_beta_cumsum)[1:K-1])
    # println(iszero.(minus_β_cumsum))
    for t in 1:T
        # Beta.(θ_t[t] .* β, θ_t[t] .* minus_β_cumsum)
        # vπ[t] ~ arraydist(Beta.(θ_t[t] .* β[1:K-1], θ_t[t] .* minus_β_cumsum[1:K-1]))#filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] ~ Dirichlet(θ_t[t] .* β) # = stickbreak(vπ[t])
        z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
        λ_z = λ[z[t]]
        x[t] ~ arraydist(Poisson.(λ_z))
    end 
end

@model function timeseries_Hdpsb_dp_pmm2(x, T, K, G;a_θ = 1,b_θ = 1,a=6,b = 1)
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    # crm = DirichletProcess.(θ_t)

    # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
    # println(size(v))
    # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
    # println(size(π_))
    # z  .~ filldist.(Categorical.(π_), nobs_t)
    # println(size(z))


    a=a
    b=b
    a_g = a .* ones(Float64,G)
    b_g = b .* ones(Float64,G)
    g = Gamma.(a_g,1 ./ b_g)
    λ ~ filldist(Product(g), K)
    # println(size(λ))
    # λ_z = Vector(map(t -> λ[z[t]], 1:T))
    # println(size(λ_z))
    γ ~ Gamma(1,1)
    # vβ ~ filldist(Beta(1,γ), K - 1)
    β ~ Dirichlet(γ .* ones(Float64, K) )#stickbreak(vβ)
    # minus_β_cumsum = 1 .- cumsum(β)
    # vπ = Vector{Vector}(undef, T)
    π_ = Vector{Vector{Float64}}(undef, T)
    z = Vector{Vector{Int64}}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    # rand.(Beta.(alpha[1] .* beta,alpha[1] .* minus_beta_cumsum)[1:K-1])
    # println(iszero.(minus_β_cumsum))
    for t in 1:T
        # Beta.(θ_t[t] .* β, θ_t[t] .* minus_β_cumsum)
        # vπ[t] ~ arraydist(Beta.(θ_t[t] .* β[1:K-1], θ_t[t] .* minus_β_cumsum[1:K-1]))#filldist(StickBreakingProcess(crm[t]), K - 1)
        π_[t] ~ Dirichlet(θ_t[t] .* β) # = stickbreak(vπ[t])
        z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
        λ_z = λ[z[t]]
        #Uniform.(rand(10), 1))
        x[t] .~ arraydist.(map(el -> Poisson.(el), λ_z))
    end 
end

# NOT COMPLETE
@model function timeseries_ddp_pmm(x,T,K,W)
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_0 ~ Gamma(a_θ,1/b_θ)
    # θ_t ~ filldist(θ_0, T)
    θ_t = Vector{Float64}(undef, T)

    a=6
    b=1
    H_0 = Gamma(a,1/b)
    λ = Vector{Vector}(undef, T)
    β = Vector{Vector}(undef, T)
    η = Vector{Vector}(undef, T) # or Vector{Float64}(undef, T)
    # λ ~ filldist(H_0, K)


    v = Vector{Vector}(undef, K-1)
    # crm = DirichletProcess.(θ_t)
    crm = Vector{DirichletProcess{Float64}}()
    
    τ = Vector{Float64}(undef, T)
    ϕ = Vector{Float64}(undef, T)
    W[1] = 0.0


   
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    
    
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            end
        end
    end
end
    
# module turingInferenceModels
#     using Random
#     using Distributions
#     using Turing
#     using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     using Test
#     import Debugger
#     using CSV,DataFrames
#     using LinearAlgebra

#     export timeseries_indep_dp_pmm1, timeseries_ddp_pmm, dp_pmm_sb1, dp_pmm_sb2,infiniteTimeSeriesPMM,timeseries_indep_dp_pmm1_fixedλ,timeseries_indep_dp_pmm2,timeseries_indep_dp_pmm_MixtureModelLiklihood1,timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood1
#     @model function timeseries_indep_dp_pmm1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # basically the same as dp_pmm_sb2 but now we can change Hyperparameters a_θ,b_θ ,a,b
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)
#         z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         # v =  rand(arraydist(filldist.(StickBreakingProcess.(crm), 12 - 1)))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))
#         # z = rand.(filldist.(Categorical.(π_), nobs_t))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # for t in 1:T
#         #     x[t] ~ arraydist(Poisson.(λ_z[t]))
#         # end 

#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             # k_t = Vector{Int}(undef, nobs_t[t])
#             for c in 1:nobs_t[t]
#                 if t ==1 
#                     π_[t] = stickbreak(v[t])
#                     z[t][c] ~ Categorical(π_[t])
#                     x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
#                 else
#                     π_[t] = stickbreak(v[t])
#                     z[t][c] ~ Categorical(π_[t])
#                     x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
#                 end
#             end
#             # k[t] = k_t
#         end
#     end
    
#     @model function timeseries_indep_dp_pmm2(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # Faster because of less loop and type definitions
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)
#         z = Vector{Vector}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))

#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
#             λ_z = λ[z[t]]
#             x[t] ~ arraydist(Poisson.(λ_z))
#         end 
#     end
    
#     @model function timeseries_indep_dp_pmm2FactorialDesign(x, T, K,factorLevels;a_θ = 1,b_θ = 1,a=6,b = 1) # Faster because of less loop and type definitions
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         if factorLevels["θt_constant"]
#             θ_t ~ arraydist(Dirac.(ones(T)))
#         else
#             θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         end
    
#         crm = DirichletProcess.(θ_t)



#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T) 
#         π_ = Vector{Vector}(undef, T)
#         if factorLevels["likelihood_form"] == "no filldist"
#             z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         else
#             z = Vector{Vector}(undef, T)#
#         end
       

#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             if factorLevels["likelihood_form"] == "MixtureModel likelihood"
#                 x[t] .~ MixtureModel(Poisson.(λ), π_[t])
#             else
#                 if factorLevels["likelihood_form"] == "no filldist"
#                     z[t] .~ Categorical(π_[t])
#                 else
#                     z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])#
#                 end
                
#                 # λ_z = λ[z[t]]
#                 if factorLevels["likelihood_form"] == "no arraydist"
#                     λ_z = λ[z[t]]
#                     x[t] .~ Poisson.(λ_z)
#                 else
#                     if factorLevels["likelihood_form"] == "direct index of λ with z "
#                         x[t] ~ arraydist(Poisson.(λ[z[t]]))
#                     else
#                         λ_z = λ[z[t]]
#                         x[t] ~ arraydist(Poisson.(λ_z))
#                     end
                    
#                 end
#                 # x[t] ~ arraydist(Poisson.(λ_z))
#             end

            
#         end 
#     end
#     @model function timeseries_indep_dp_pmm3(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # Work in Progress; Tryna make it faster
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         # crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)
#         z = Vector{Vector}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))

#         for t in 1:T
#             v[t] ~ filldist(Beta(1,θ_t[t]), K - 1)#filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
#             λ_z = λ[z[t]]
#             x[t] ~ arraydist(Poisson.(λ_z))
#         end 
#     end
#     @model function timeseries_indep_dp_pmm_MixtureModelLiklihood1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             x[t] .~ MixtureModel(Poisson.(λ), π_[t])
#         end 
#     end

#     @model function timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood1(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
        
#         #OUTPUT IS A MATRIX
#         λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

#         H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end
#     @model function timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood2(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
        
#         #OUTPUT IS A MATRIX
#         # λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)
#         # H = Product.(map(j -> Poisson.(λ)[:,j],1:K))

#         λ = Vector(map(k -> tzeros(Float64, G),1:K))
#         λ .~ map(k-> filldist(Gamma(ak[k],1/bk[k]),G)  ,1:K)
#         H = Product.(map(j -> Poisson.(λ[j]),1:K))


#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end
#     @model function timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood22(x, T, K,G,a,b_ak,b_bk;a_θ = 1,b_θ = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)
    
#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))
    
#         # filldist.(Gamma.(a, 1.0 ./b), G)
        
#         a=a
#         b = TArray(Float64,K)#Vector{Float64}(undef,K)
#         # λ =  TArray(TArray,K)
#         # H = Vector{Distribution}(undef,K)
#         println("DUM BDU")
#         # println(λ)
#         for k in 1:K
#             b[k] ~  Gamma(b_ak[k], b_bk[k])
#             # λ[k] ~  filldist(Gamma(a, 1.0 ./ b[k]), G) #Vector{Vector}(undef,K)
#         end
#         println("1")
#         λ ~ arraydist(filldist.(Gamma.(1,b),G))
#         println(λ)
#         H =  Product.(map(k -> Poisson.(λ[:,k]),1:K))
#         println("First Class")
#         # b=b
        
#         #OUTPUT IS A MATRIX
#         # println("here1")
    
#         # H = Product.(map(j -> Poisson.(λ[j]),1:K))
    
#         #filldist(Gamma(a,1/b), K)
#         # println("here2")S
        
#         # println("here3")
#         #  Vector(map(k -> filldist(Gamma(a[k],1/b[k]),G),1:K)) #Vector(map(k -> Vector{Float64}(undef,G),1:K))  #  Vector(map(k -> TArray(Float64,G),1:K))#
#         # cluster_gamma_prior = Gamma.(a, 1 ./ b)
#         # # λ .~ vec(map(k -> filldist(Gamma(a[k],1/b[k]),G),1:K))
#         # # λ .~ map(k-> filldist(Gamma(a[k],1/b[k]),G)  ,1:K)
#         # # # println(typeof(λ))
#         # for k in 1:K
#         #     λ[k] ~ filldist(cluster_gamma_prior[k],G)
#         # end
        
        
        
#         # tmpPrior =  Vector(map(dist-> filldist(dist,G), Gamma.(a, 1 ./ b)))
#         # println(H)
#         # # println(size(tmpPrior))
#         # # println(typeof(tmpPrior))
#         # λ .~ Vector(map(dist-> filldist(dist,G), Gamma.(a, 1 ./ b)))
#         # println("here4")
#         # H = Product.(map(j -> Poisson.(λ[j]),1:K))
#         # println(typeof(H))
#         # println("Here5")
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))
    
#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)
    
#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] = stickbreak(v[t])
#             # println("Here6")
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end
#     @model function timeseries_hdp_pmm_MultivariateMixtureModelLiklihood1(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1, γ = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         γ = γ
#         vβ ~  filldist(Beta(1,γ), K - 1)
#         β = stickbreak(vβ)
#         #OUTPUT IS A MATRIX
#         λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

#         H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             c =  θ_t[t] .* β
#             d =  θ_t[t] .*( 1 .-  cumsum(β))
#             c = c[1:end-1]
#             d = d[1:end-1]
#             # v[t] = tzeros(Float64, K-1)
#             v[t] ~ arraydist(Beta.(c,d))
#             π_[t] = stickbreak(v[t])
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end

#     @model function timeseries_hdp_pmm_MultivariateMixtureModelLiklihood2(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1, γ = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         γ = γ
#         vβ ~  filldist(Beta(1,γ), K - 1)
#         β = stickbreak(vβ)
#         #OUTPUT IS A MATRIX
#         λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

#         H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             π_[t] ~  Dirichlet(θ_t[t] .* β)
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end

#     @model function timeseries_hdp_gmm_MultivariateMixtureModelLiklihood1(x, T, K,G;a_θ = 1,b_θ = 1,a=1,b = 1, γ = 1, λ0 = 1, μ0=1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)



#         a=a
#         b=b
#         γ = γ
#         vβ ~  filldist(Beta(1,γ), K - 1)
#         β = stickbreak(vβ)
        
#         τk ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)
#         μk ~ arraydist(Normal.(μ0,(λ0 .* τk) .^(-1.)))
#         H = map(j -> MultivariateNormal(μk[:,j],diagm(τk[:,j] .^(-1.) ) ),1:K)

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             π_[t] ~  Dirichlet(θ_t[t] .* β)
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end
    

#     @model function timeseries_hdp_pmm_MultivariateMixtureModelLiklihood3(x, T, K,G;a_θ = 1,b_θ = 1,a=6,b = 1, γ = 1) # s
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         # crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         γ = γ
#         β ~  Dirichlet(γ .* ones(K) ./ K)
#         #OUTPUT IS A MATRIX
#         λ ~  filldist(filldist(Gamma(a,1/b), G),K)#filldist(Gamma(a,1/b), K)

#         H = Product.(map(j -> Poisson.(λ[:,j]),1:K)) # vs Product.(map(j -> Poisson.(λ)[:,j],1:K))

#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)

#         for t in 1:T
#             π_[t] ~  Dirichlet(θ_t[t] .* β)
#             x[t] .~ MixtureModel(H, π_[t])
#         end 
#     end

    

#     export timeseries_indep_dp_gmm1
#     @model function timeseries_indep_dp_gmm1(x, T, K;a_θ = 1,b_θ = 1) 
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)

#         s² ~ filldist(InverseGamma(2, 3), K)
#         m ~ arraydist(Normal.(0, sqrt.(s²)))
#         H =  Normal.(m, sqrt.(s²))

#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)
#         z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             # k_t = Vector{Int}(undef, nobs_t[t])
#             for c in 1:nobs_t[t]
#                 if t ==1 
#                     π_[t] = stickbreak(v[t])
#                     z[t][c] ~ Categorical(π_[t])
#                     x[t][c] ~ H[z[t][c]]#MixtureModel(Poisson,λ, π_[t])
#                 else
#                     π_[t] = stickbreak(v[t])
#                     z[t][c] ~ Categorical(π_[t])
#                     x[t][c] ~ H[z[t][c]]#MixtureModel(Poisson,λ, π_[t])
#                 end
#             end
#             # k[t] = k_t
#         end
#     end
    
#     #Hierarchical 
#     #Work in Progress
#     export timeseries_Hdpsb_dp_pmm1,timeseries_Hdpsb_dp_pmm2
#     @model function timeseries_Hdpsb_dp_pmm1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1)
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         # crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))
#         γ ~ Gamma(1,1)
#         # vβ ~ filldist(Beta(1,γ), K - 1)
#         β ~ Dirichlet(γ .* ones(Float64, K) )#stickbreak(vβ)
#         # minus_β_cumsum = 1 .- cumsum(β)
#         # vπ = Vector{Vector}(undef, T)
#         π_ = Vector{Vector{Float64}}(undef, T)
#         z = Vector{Vector{Int64}}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         # rand.(Beta.(alpha[1] .* beta,alpha[1] .* minus_beta_cumsum)[1:K-1])
#         # println(iszero.(minus_β_cumsum))
#         for t in 1:T
#             # Beta.(θ_t[t] .* β, θ_t[t] .* minus_β_cumsum)
#             # vπ[t] ~ arraydist(Beta.(θ_t[t] .* β[1:K-1], θ_t[t] .* minus_β_cumsum[1:K-1]))#filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] ~ Dirichlet(θ_t[t] .* β) # = stickbreak(vπ[t])
#             z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
#             λ_z = λ[z[t]]
#             x[t] ~ arraydist(Poisson.(λ_z))
#         end 
#     end
    
#     @model function timeseries_Hdpsb_dp_pmm2(x, T, K, G;a_θ = 1,b_θ = 1,a=6,b = 1)
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         # crm = DirichletProcess.(θ_t)

#         # v ~  arraydist(filldist.(Beta.(1,θ_t), K - 1))
#         # println(size(v))
#         # π_ =  Vector(map(t -> stickbreak(v[:,t]), 1:T))#stickbreak.(v)
#         # println(size(π_))
#         # z  .~ filldist.(Categorical.(π_), nobs_t)
#         # println(size(z))


#         a=a
#         b=b
#         a_g = a .* ones(Float64,G)
#         b_g = b .* ones(Float64,G)
#         g = Gamma.(a_g,1 ./ b_g)
#         λ ~ filldist(Product(g), K)
#         # println(size(λ))
#         # λ_z = Vector(map(t -> λ[z[t]], 1:T))
#         # println(size(λ_z))
#         γ ~ Gamma(1,1)
#         # vβ ~ filldist(Beta(1,γ), K - 1)
#         β ~ Dirichlet(γ .* ones(Float64, K) )#stickbreak(vβ)
#         # minus_β_cumsum = 1 .- cumsum(β)
#         # vπ = Vector{Vector}(undef, T)
#         π_ = Vector{Vector{Float64}}(undef, T)
#         z = Vector{Vector{Int64}}(undef, T)#Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         # rand.(Beta.(alpha[1] .* beta,alpha[1] .* minus_beta_cumsum)[1:K-1])
#         # println(iszero.(minus_β_cumsum))
#         for t in 1:T
#             # Beta.(θ_t[t] .* β, θ_t[t] .* minus_β_cumsum)
#             # vπ[t] ~ arraydist(Beta.(θ_t[t] .* β[1:K-1], θ_t[t] .* minus_β_cumsum[1:K-1]))#filldist(StickBreakingProcess(crm[t]), K - 1)
#             π_[t] ~ Dirichlet(θ_t[t] .* β) # = stickbreak(vπ[t])
#             z[t] ~ filldist(Categorical(π_[t]), nobs_t[t])
#             λ_z = λ[z[t]]
#             #Uniform.(rand(10), 1))
#             x[t] .~ arraydist.(map(el -> Poisson.(el), λ_z))
#         end 
#     end



#     # One time Point
#     export _indep_dp_gmm1,_indep_dp_pmm1
#     @model function _indep_dp_gmm1(x, K;a_θ = 1,b_θ = 1,a=6,b = 1) 
#         nobs = length(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ Gamma(a_θ,1/b_θ) 
#         crm = DirichletProcess(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
#         a=a
#         b=b
#         s² ~ filldist(InverseGamma(2, 3), K)
#         sqrt_s = sqrt.(s²)
#         m ~ arraydist(Normal.(0, sqrt.(s²)))
#         H =  Normal.(m, sqrt.(s²))
        
#         v ~ filldist(StickBreakingProcess(crm), K - 1)
#         π_ =  stickbreak(v)
#         z = tzeros(Int, nobs)
#         for c in 1:nobs
#             z[c] ~ Categorical(π_)
#             x[c] ~  H[z[c]]#Normal(m[z[c]], sqrt.(s²)[z[c]])#MixtureModel(Poisson,λ, π_[t])
#         end
#         return x,z
#     end
#     @model function _indep_dp_pmm1(x, K;a_θ = 1,b_θ = 1,a=6,b = 1) 
#         nobs = length(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ Gamma(a_θ,1/b_θ) 
#         crm = DirichletProcess(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
#         a=a
#         b=b
#         λ ~ filldist(Gamma(a,1/b), K)
#         v ~ filldist(StickBreakingProcess(crm), K - 1)
#         π_ =  stickbreak(v)
#         z = tzeros(Int, nobs)
#         for c in 1:nobs
#             z[c] ~ Categorical(π_)
#             x[c] ~ Poisson(λ[z[c]])#MixtureModel(Poisson,λ, π_[t])
#         end
#         return x,z
#     end

#     export timeseries_indep_dp_pmm1_fixedλ
#     @model function timeseries_indep_dp_pmm1_fixedλ(x, T, K,true_λ;a_θ = 1,b_θ = 1,a=6,b = 1) # basically the same as dp_pmm_sb2 but now we can change Hyperparameters a_θ,b_θ ,a,b
#         nobs_t = length.(x)
#         a_θ = a_θ
#         b_θ = b_θ
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
#         a=a
#         b=b
#         λ  ~ arraydist(Dirac.(true_λ))
#         v = Vector{Vector}(undef, T)
#         π_ = Vector{Vector}(undef, T)
#         z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             # k_t = Vector{Int}(undef, nobs_t[t])
#             for c in 1:nobs_t[t]
#                 if t ==1 
#                     π_[t] = stickbreak(v[t])
#                     z[t][c] ~ Categorical(π_[t])
#                     x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
#                 else
#                     π_[t] = stickbreak(v[t])
#                     z[t][c] ~ Categorical(π_[t])
#                     x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
#                 end
#             end
#             # k[t] = k_t
#         end
#     end

#     @model function infiniteTimeSeriesPMM(x,T)
#         # Hyper-parameters, i.e. concentration parameter and parameters of H.
    
#         # α = 1.0
#         # μ0 = 0.0
#         # σ0 = 1.0
#         nobs_t = length.(x)
#         a_θ = 1
#         b_θ = 1
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    
#         a=6
#         b=1
        
        
    
    
#         rpm = DirichletProcess.(θ_t)
        
#         # Define the base distribution, i.e. expected value of the Dirichlet process.
#         H = Gamma(a,1/b)
        
#         # Latent assignment.
#         z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
            
#         # Locations of the infinitely many clusters.
#         λ = tzeros(Float64,0)
#         nk = Vector{Vector{Int}}(undef, T)
#         for t in 1:T
#             for c in 1:nobs_t[t]
#                 # Number of clusters.
#                 K = maximum(vcat(z[1:t]...))
#                 nk[t] = Vector{Int}(map(k -> sum(z[t] .== k), 1:K))
    
#                 # Draw the latent assignment.
#                 if t==1
#                     n__ = nk[t]#z[t][c] ~ ChineseRestaurantProcess(rpm[t], nk[t])
#                 else
#                     n__ =  nk[t] .+ nk[t-1]#z[t][c] ~ ChineseRestaurantProcess(rpm[t], nk[t] .+ nk[t-1])
#                 end
#                 z[t][c] ~ ChineseRestaurantProcess(rpm[t], n__)
                
            
#                 # Create a new cluster?
#                 if z[t][c] > K
#                     if t != 1
#                         for i in 1:t-1
#                             push!(nk[i],zero(eltype(nk[i])))
#                         end
#                     end
#                     push!(λ, 0.0)
    
#                     # Draw location of new cluster.
#                     # println("Here at $K")
#                     λ[z[t][c]] ~ H
#                 end
                        
#                 # Draw observation.
#                 x[t][c] ~ Poisson(λ[z[t][c]])
#             end
#         end
#     end
    
#     # m3 = infiniteTimeSeriesPMM(x,T)
#     # iterations = 1000
#     # @time smc_chain = begin
#     #     burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
#     #     n_samples = 1000
#     #     iterations = burn + n_samples
     
#     #     chain = sample(m3, SMC(), iterations);
#     # end

#     @model dp_pmm_sb2(x, T, K) = begin
#         nobs_t = length.(x)
#         a_θ = 1
#         b_θ = 1
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
#         a=6
#         b=1
#         λ ~ filldist(Gamma(a,1/b), K)
#         v = Vector{Vector}(undef, K-1)
#         z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             # k_t = Vector{Int}(undef, nobs_t[t])
#             for c in 1:nobs_t[t]
#                 if t ==1 
#                     eta = stickbreak(v[t])
#                     z[t][c] ~ Categorical(eta)
#                     x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
#                 else
#                     eta = stickbreak(v[t])
#                     z[t][c] ~ Categorical(eta)
#                     x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
#                 end
#             end
#             # k[t] = k_t
#         end

        
        
#     end
#     # burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
#     # n_samples = 500
#     # iterations = burn + n_samples
#     # K = 25
#     # T = 4
#     # m2=dp_pmm_sb2(x, T, K)
#     # # pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))
#     # tchain = sample(m2, PG(100),iterations);
#     # k = map(
#     #     t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :k), :].value))),
#     #     1:iterations
#     # );
#     # ids = findall(map(name -> occursin("λ", string(name)), names(tchain)));
#     # p = plot(tchain[:, ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter)


#     # DP PMM model under stick-breaking construction
#     @model dp_pmm_sb1(x, T, K) = begin
#         nobs_t = length.(x)
#         a_θ = 1
#         b_θ = 1
#         θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
#         crm = DirichletProcess.(θ_t)
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
#         a=6
#         b=1
#         λ ~ filldist(Gamma(a,1/b), K)
#         v = Vector{Vector}(undef, K-1)
#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             for c in 1:nobs_t[t]
#                 if t ==1 
#                     eta = stickbreak(v[t])
#                     x[t][c] ~ MixtureModel(Poisson,λ, eta)
#                 else
#                     eta = stickbreak(v[t])
#                     x[t][c] ~ MixtureModel(Poisson,λ, eta)
#                 end
#             end
#         end

        
        
#     end
    
#     # x = [[data_dict[(i,t)][1] for i in 1:6] for t in 1:4]
#     # true_clusters = [[truth_dict[(i,t)] for i in 1:6] for t in 1:4]
#     # coordinates = [[(i,t) for i in 1:6] for t in 1:4]
#     # K = 5
#     # T = 4
#     # m=dp_pmm_sb1(x, T, K)
    
#     # # Set random seed for reproducibility
#     # Random.seed!(0);
    
#     # # Compile time approx. 32s.
#     # # Run time approx. 70s.
    
#     # @time hmc_chain = begin
#     #     burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
#     #     n_samples = 500
#     #     iterations = burn + n_samples
#     #     n_components = 10
#     #     stepsize = 0.01
#     #     nleapfrog = floor(Int, 1 / stepsize)
     
#     #     chain = sample(dp_pmm_sb1(x, T, K), 
#     #            HMC(stepsize, nleapfrog),
#     #            iterations)
#     # end
#     # ;


#     # NOT COMPLETE
#     @model function timeseries_ddp_pmm(x,T,K,W)
#         nobs_t = length.(x)
#         a_θ = 1
#         b_θ = 1
#         θ_0 ~ Gamma(a_θ,1/b_θ)
#         # θ_t ~ filldist(θ_0, T)
#         θ_t = Vector{Float64}(undef, T)
    
#         a=6
#         b=1
#         H_0 = Gamma(a,1/b)
#         λ = Vector{Vector}(undef, T)
#         β = Vector{Vector}(undef, T)
#         η = Vector{Vector}(undef, T) # or Vector{Float64}(undef, T)
#         # λ ~ filldist(H_0, K)
    
    
#         v = Vector{Vector}(undef, K-1)
#         # crm = DirichletProcess.(θ_t)
#         crm = Vector{DirichletProcess{Float64}}()
        
#         τ = Vector{Float64}(undef, T)
#         ϕ = Vector{Float64}(undef, T)
#         W[1] = 0.0
    
    
       
#         # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
        
        
#         for t in 1:T
#             v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
#             for c in 1:nobs_t[t]
#                 if t ==1 
#                     eta = stickbreak(v[t])
#                     x[t][c] ~ MixtureModel(Poisson,λ, eta)
#                 else
#                     eta = stickbreak(v[t])
#                     x[t][c] ~ MixtureModel(Poisson,λ, eta)
#                 end
#             end
#         end
#     end
    
    
# end



