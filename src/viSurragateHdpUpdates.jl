"""
"""
function genterate_Delta_mk(rhok_hat_vec, omegak_hat_vec)
    K = length(rhok_hat_vec)
    e_βk = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    Δ_mk =  Matrix{Float64}(undef,K,K+1)
    for k in 1:K+1
        for m in 1:K
            if m < k
                Δ_mk[m,k] = -1.0 /(1 - rhok_hat_vec[m]) * e_βk[k]
            elseif m == k
                Δ_mk[m,k] = 1.0 / rhok_hat_vec[m] * e_βk[k]
            elseif m > 0
                Δ_mk[m,k] = 0.0
            end
        end
    end
    return Δ_mk
end

"""
"""
function SurragateLowerBound_util(vec_args,T,γ,α_0,T_k)
    rho_hat = vec_args[1, :]
    omega_hat = vec_args[2, :]
    calc_SurragateLowerBound(rho_hat,omega_hat,T,γ,α_0,T_k)
end

"""
"""
function SurragateLowerBound_unconstrained_util(vec_args,T,γ,α_0,T_k)
    c = vec_args[1, :] # rho_hat
    d = vec_args[2, :] # omega_hat
    calc_SurragateLowerBound_unconstrained(c,d,T,γ,α_0,T_k)
end

"""
"""
function g_constrained!(GD, x, T,γ,α_0,T_k)
    rho_hat =  x[1,:]
    omega_hat = x[2,:]

    K = length(x[1,:])
    k_vec = collect(1:K)
    rho_omega_hat =  rho_hat.*omega_hat 
    minusrho_omega_hat =  (1 .- rho_hat) .* omega_hat 
    minusrho_hat = (1 .- rho_hat)
    Δ_mk = genterate_Delta_mk(rho_hat, omega_hat)
    Δ_prod_T =  Δ_mk[:,1:K] .* T_k
    a_Δ_prod_T = α_0 .* permutedims(sum(Δ_prod_T,dims=1))
    #take Negative to find max  

    GD[1,:] = 1.0 .* ( omega_hat .* (T .+ 1 .-  rho_omega_hat) .* polygamma.(1,rho_omega_hat) .- omega_hat .*( T .*(K .+ 1.0 .- k_vec)  .+ γ .- minusrho_omega_hat) .* polygamma.(1,minusrho_omega_hat) .+ a_Δ_prod_T) #  rho
    #-2.0 .* (1.0 .- x[1,:]) .- 400.0 .* (x[2,:] .- x[1,:].^2) .* x[1,:]
    GD[2,:] = 1.0 .* ((T .+ 1.0 .- rho_omega_hat) .* (rho_hat .* polygamma.(1,rho_omega_hat) .- polygamma.(1,omega_hat)) .+ ( T .* (K .+ 1.0 .- k_vec) .+ γ .-  minusrho_omega_hat) .* ( minusrho_hat .* polygamma.(1, minusrho_omega_hat)  .- polygamma.(1,omega_hat) )) #  omega
     #(x[2,:] .- x[1,:].^2)
end

"""
"""
function g_unconstrained!(GD, x, T,γ,α_0,T_k)
    rho_hat =  x[1,:] # c
    omega_hat = x[2,:] # d

    K = length(x[1,:])
    k_vec = collect(1:K)
    rho_omega_hat =  rho_hat.*omega_hat 
    minusrho_omega_hat =  (1 .- rho_hat) .* omega_hat 
    minusrho_hat = (1 .- rho_hat)
    Δ_mk = genterate_Delta_mk(rho_hat,omega_hat)
    Δ_prod_T =  Δ_mk[:,1:K] .* T_k[1:K]
    a_Δ_prod_T = α_0 .* permutedims(sum(Δ_prod_T,dims=1))
    #take Negative to find max  

    GD[1,:] = 1.0 .*(rho_hat .* minusrho_hat .* (omega_hat .* (T .+ 1 .-  rho_omega_hat) .* polygamma.(1,rho_omega_hat) .- omega_hat .*( T .*(K .+ 1.0 .- k_vec)  .+ γ .- minusrho_omega_hat) .* polygamma.(1,minusrho_omega_hat) .+ a_Δ_prod_T)) # c and rho
    #-2.0 .* (1.0 .- x[1,:]) .- 400.0 .* (x[2,:] .- x[1,:].^2) .* x[1,:] 
    GD[2,:] = 1.0 .* ( omega_hat .* ((T .+ 1.0 .- rho_omega_hat) .* (rho_hat .* polygamma.(1,rho_omega_hat) .- polygamma.(1,omega_hat)) .+ ( T .* (K .+ 1.0 .- k_vec) .+ γ .-  minusrho_omega_hat) .* ( minusrho_hat .* polygamma.(1, minusrho_omega_hat)  .- polygamma.(1,omega_hat) ))) # d and omega
     #(x[2,:] .- x[1,:].^2)
end

"""
"""
g_constrained_closure!(T,γ,α_0,T_k) = (GD, x)  ->  g_constrained!(GD, x, T,γ,α_0,T_k)

"""
"""
g_unconstrained_closure!(T,γ,α_0,T_k) = (GD, x)  ->  g_unconstrained!(GD, x, T,γ,α_0,T_k)


"""
"""
SurragateLowerBound_closure(T,γ,α_0,T_k) = vec_args ->  SurragateLowerBound_util(vec_args,T,γ,α_0,T_k)

"""
"""
SurragateLowerBound_unconstrained_closure(T,γ,α_0,T_k) = vec_args ->  SurragateLowerBound_unconstrained_util(vec_args,T,γ,α_0,T_k)

