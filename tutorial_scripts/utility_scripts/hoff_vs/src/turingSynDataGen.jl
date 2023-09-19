
@model function gen1Phen_SimpleCase1_NB_TimeInvar(x,T,k,C_t,μ,σ²,α,mixing_prob)
    if x === missing
        x = Vector{Vector}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = tzeros(Int, cells)
        end
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end



    if μ >= σ²
        throw(DomainError(μ, "μ must be less than σ²"))
    end
    if k !== maximum([length(μ), length(σ²)])
        error("must ha μ and σ² for each of the k components")
    end


    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²)

    H = NegativeBinomial.(r,p)
    z = Vector{Vector}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z

    
end  


@model function gen1Phen_SimpleCase1_Poisson_TimeInvar(x,T,k,C_t,μ,α,mixing_prob)
    if x === missing
        x = Vector{Vector}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = tzeros(Int, cells)
        end
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end

    if k !== length(μ)
        error("must have μ for each of the k components")
    end


    # r = calc_r(μ,σ²)
    # p = calc_p(μ,σ²)

    H = Poisson.(μ)
    z = Vector{Vector}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z

    
end  

@model function gen1Phen_SimpleCase1_PoissonGamma_TimeInvar(x,T,k,C_t,α,mixing_prob)
    if x === missing
        x = Vector{Vector}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = tzeros(Int, cells)
        end
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end


    a =  6

    b = 1
    λ ~ filldist(Gamma(a,1/b), k)

    # r = calc_r(μ,σ²)
    # p = calc_p(μ,σ²)

    H = Poisson.(λ)
    z = Vector{Vector}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z

    
end 

@model function gen1Phen_SimpleCase1_Normal_TimeInvar(x,T,k,C_t,α,mixing_prob)
    if x === missing
        x = Vector{Vector}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = tzeros(Float64, cells)
        end
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end



    # if μ >= σ²
    #     throw(DomainError(μ, "μ must be less than σ²"))
    # end
    # if k !== maximum([length(μ), length(σ²)])
    #     error("must ha μ and σ² for each of the k components")
    # end


    s² ~ filldist(InverseGamma(2, 3), k)
    m ~ arraydist(Normal.(0, sqrt.(s²)))

    H = Normal.(m,s²)
    z = Vector{Vector}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z

    
end  


@model function genGPhen_SimpleCase1_NB_TimeInvar(x,T,G,k,C_t,μ,σ²,α,mixing_prob)
    if x === missing
        x = Vector{Vector{Vector{Int64}}}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = map(c -> tzeros(Int, G), 1:cells)
        end
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end



    if all(μ .>= σ²)
        throw(DomainError(μ, "μ must be less than σ²"))
    end
    if k !== maximum([length(μ), length(σ²)])
        error("must ha μ and σ² for each of the k components")
    end


    r = calc_r.(μ,σ²)
    p = calc_p.(μ,σ²)



    H =  vec(map((n,m) -> filldist(NegativeBinomial(n,m), G), r,p))#NegativeBinomial.(r,p)
    z = Vector{Vector{Int}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z

    
end  


@model function NBMixtureModel(x,k,cells,μ,σ²,α,mixing_prob)
    if x === missing
        x = tzeros(Int, cells)
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end



    if μ >= σ²
        throw(DomainError(μ, "μ must be less than σ²"))
    end
    if k !== maximum([length(μ), length(σ²)])
        error("must ha μ and σ² for each of the k components")
    end


    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²)

    H = NegativeBinomial.(r,p)
    z = Vector{Int}(undef,cells)
    for c in 1:cells
        z[c] ~ Categorical(_prob)
        x[c] ~ H[z[c]]
    end 
    return x,z

    
end  

@model function NormalMixtureModel(x,k,cells,μ,σ²,α,mixing_prob)
    if x === missing
        x = tzeros(Int, cells)
    end
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end



    if μ >= σ²
        throw(DomainError(μ, "μ must be less than σ²"))
    end
    if k !== maximum([length(μ), length(σ²)])
        error("must ha μ and σ² for each of the k components")
    end


    # r = calc_r(μ,σ²)
    # p = calc_p(μ,σ²)

    H = Normal.(μ,σ²)
    z = Vector{Int}(undef,cells)
    for c in 1:cells
        z[c] ~ Categorical(_prob)
        x[c] ~ H[z[c]]
    end 
    return x,z

    
end  


@model function gen4gene_Corr_PGMvLN(x,G,C,Ω;a=0.5,b=0.6)
    if x === missing
        x = Vector{Vector}(undef,C)
    end
    Σ = exp.(Ω) .-  ones(Int,G) * transpose(ones(Int,G))
    λ = Vector{Vector}(undef,C)
    u = Vector{Vector}(undef,C)
    E = MvLogNormal(ones(Int,G),Σ)
    for c in 1:C
        u[c] ~ E
        λ[c] ~ filldist(Gamma(a,b),G)
        x[c] ~ arraydist(Poisson.(u[c] .* λ[c]))
    end

    
end

# module turingSytheticDataGeneration
#     # include("MathUtils.jl")
#     # using .MathUtils


#     # using Random
#     # using Distributions
#     # using Turing
#     # using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     # using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     # using Test
#     # import Debugger
#     # using CSV,DataFrames

#     export  gen1Phen_SimpleCase1_NB_TimeInvar , gen4gene_Corr_PGMvLN,gen1Phen_SimpleCase1_Normal_TimeInvar,NBMixtureModel,NormalMixtureModel,genGPhen_SimpleCase1_NB_TimeInvar,gen1Phen_SimpleCase1_Poisson_TimeInvar,gen1Phen_SimpleCase1_PoissonGamma_TimeInvar

#     @model function gen1Phen_SimpleCase1_NB_TimeInvar(x,T,k,C_t,μ,σ²,α,mixing_prob)
#         if x === missing
#             x = Vector{Vector}(undef,T)
#             for t in 1:T
#                 cells = C_t[t]
#                 x[t] = tzeros(Int, cells)
#             end
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end



#         if μ >= σ²
#             throw(DomainError(μ, "μ must be less than σ²"))
#         end
#         if k !== maximum([length(μ), length(σ²)])
#             error("must ha μ and σ² for each of the k components")
#         end


#         r = calc_r(μ,σ²)
#         p = calc_p(μ,σ²)

#         H = NegativeBinomial.(r,p)
#         z = Vector{Vector}(undef,T)
#         for t in 1:T
#             cells = C_t[t]
#             z[t] = tzeros(Int, cells)
#             for c in 1:cells
#                 z[t][c] ~ Categorical(_prob)
#                 x[t][c] ~ H[z[t][c]]
#             end 
            
#         end
#         return x,z

        
#     end  

    
#     @model function gen1Phen_SimpleCase1_Poisson_TimeInvar(x,T,k,C_t,μ,α,mixing_prob)
#         if x === missing
#             x = Vector{Vector}(undef,T)
#             for t in 1:T
#                 cells = C_t[t]
#                 x[t] = tzeros(Int, cells)
#             end
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end

#         if k !== length(μ)
#             error("must have μ for each of the k components")
#         end


#         # r = calc_r(μ,σ²)
#         # p = calc_p(μ,σ²)

#         H = Poisson.(μ)
#         z = Vector{Vector}(undef,T)
#         for t in 1:T
#             cells = C_t[t]
#             z[t] = tzeros(Int, cells)
#             for c in 1:cells
#                 z[t][c] ~ Categorical(_prob)
#                 x[t][c] ~ H[z[t][c]]
#             end 
            
#         end
#         return x,z

        
#     end  

#     @model function gen1Phen_SimpleCase1_PoissonGamma_TimeInvar(x,T,k,C_t,α,mixing_prob)
#         if x === missing
#             x = Vector{Vector}(undef,T)
#             for t in 1:T
#                 cells = C_t[t]
#                 x[t] = tzeros(Int, cells)
#             end
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end

    
#         a =  6

#         b = 1
#         λ ~ filldist(Gamma(a,1/b), k)

#         # r = calc_r(μ,σ²)
#         # p = calc_p(μ,σ²)

#         H = Poisson.(λ)
#         z = Vector{Vector}(undef,T)
#         for t in 1:T
#             cells = C_t[t]
#             z[t] = tzeros(Int, cells)
#             for c in 1:cells
#                 z[t][c] ~ Categorical(_prob)
#                 x[t][c] ~ H[z[t][c]]
#             end 
            
#         end
#         return x,z

        
#     end 

#     @model function gen1Phen_SimpleCase1_Normal_TimeInvar(x,T,k,C_t,α,mixing_prob)
#         if x === missing
#             x = Vector{Vector}(undef,T)
#             for t in 1:T
#                 cells = C_t[t]
#                 x[t] = tzeros(Float64, cells)
#             end
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end



#         # if μ >= σ²
#         #     throw(DomainError(μ, "μ must be less than σ²"))
#         # end
#         # if k !== maximum([length(μ), length(σ²)])
#         #     error("must ha μ and σ² for each of the k components")
#         # end


#         s² ~ filldist(InverseGamma(2, 3), k)
#         m ~ arraydist(Normal.(0, sqrt.(s²)))

#         H = Normal.(m,s²)
#         z = Vector{Vector}(undef,T)
#         for t in 1:T
#             cells = C_t[t]
#             z[t] = tzeros(Int, cells)
#             for c in 1:cells
#                 z[t][c] ~ Categorical(_prob)
#                 x[t][c] ~ H[z[t][c]]
#             end 
            
#         end
#         return x,z

        
#     end  
    

#     @model function genGPhen_SimpleCase1_NB_TimeInvar(x,T,G,k,C_t,μ,σ²,α,mixing_prob)
#         if x === missing
#             x = Vector{Vector{Vector{Int64}}}(undef,T)
#             for t in 1:T
#                 cells = C_t[t]
#                 x[t] = map(c -> tzeros(Int, G), 1:cells)
#             end
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end



#         if all(μ .>= σ²)
#             throw(DomainError(μ, "μ must be less than σ²"))
#         end
#         if k !== maximum([length(μ), length(σ²)])
#             error("must ha μ and σ² for each of the k components")
#         end


#         r = calc_r.(μ,σ²)
#         p = calc_p.(μ,σ²)



#         H =  vec(map((n,m) -> filldist(NegativeBinomial(n,m), G), r,p))#NegativeBinomial.(r,p)
#         z = Vector{Vector{Int}}(undef,T)
#         for t in 1:T
#             cells = C_t[t]
#             z[t] = tzeros(Int, cells)
#             for c in 1:cells
#                 z[t][c] ~ Categorical(_prob)
#                 x[t][c] ~ H[z[t][c]]
#             end 
            
#         end
#         return x,z

        
#     end  


#     @model function NBMixtureModel(x,k,cells,μ,σ²,α,mixing_prob)
#         if x === missing
#             x = tzeros(Int, cells)
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end
    
    
    
#         if μ >= σ²
#             throw(DomainError(μ, "μ must be less than σ²"))
#         end
#         if k !== maximum([length(μ), length(σ²)])
#             error("must ha μ and σ² for each of the k components")
#         end
    
    
#         r = calc_r(μ,σ²)
#         p = calc_p(μ,σ²)
    
#         H = NegativeBinomial.(r,p)
#         z = Vector{Int}(undef,cells)
#         for c in 1:cells
#             z[c] ~ Categorical(_prob)
#             x[c] ~ H[z[c]]
#         end 
#         return x,z
    
        
#     end  
    
#     @model function NormalMixtureModel(x,k,cells,μ,σ²,α,mixing_prob)
#         if x === missing
#             x = tzeros(Int, cells)
#         end
#         if mixing_prob === missing
#             _prob ~ Dirichlet(α .* ones(k) ./ k)
#         else
#             _prob ~ arraydist(Dirac.(mixing_prob))
#         end
    
    
    
#         if μ >= σ²
#             throw(DomainError(μ, "μ must be less than σ²"))
#         end
#         if k !== maximum([length(μ), length(σ²)])
#             error("must ha μ and σ² for each of the k components")
#         end
    
    
#         # r = calc_r(μ,σ²)
#         # p = calc_p(μ,σ²)
    
#         H = Normal.(μ,σ²)
#         z = Vector{Int}(undef,cells)
#         for c in 1:cells
#             z[c] ~ Categorical(_prob)
#             x[c] ~ H[z[c]]
#         end 
#         return x,z
    
        
#     end  

    
#     @model function gen4gene_Corr_PGMvLN(x,G,C,Ω;a=0.5,b=0.6)
#         if x === missing
#             x = Vector{Vector}(undef,C)
#         end
#         Σ = exp.(Ω) .-  ones(Int,G) * transpose(ones(Int,G))
#         λ = Vector{Vector}(undef,C)
#         u = Vector{Vector}(undef,C)
#         E = MvLogNormal(ones(Int,G),Σ)
#         for c in 1:C
#             u[c] ~ E
#             λ[c] ~ filldist(Gamma(a,b),G)
#             x[c] ~ arraydist(Poisson.(u[c] .* λ[c]))
#         end

        
#     end
    
# end
