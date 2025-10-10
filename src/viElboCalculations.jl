"""
    Lp_data(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculate the expected log likelihood of the data given the variational parameters.
"""
function Lp_data(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        @simd for j in 1:J
            @fastmath @inbounds  LB += ( -0.5 * clusters[k].Nk[1]* dataparams.logpi - 0.5 *clusters[k].Nk[1]*E_ln_sigma_sq(j,clusters[k]) )- 0.5*E_one_over_sigma_sq(j,clusters[k]) *( clusters[k].x_hat_sq[j] -2*clusters[k].m_mu[j]* clusters[k].y[j]*clusters[k].x_hat[j]-2*clusters[k].m_nu[j]*clusters[k].x_hat[j]+2*clusters[k].Nk[1]*clusters[k].m_mu[j]* clusters[k].y[j]*clusters[k].m_nu[j] + clusters[k].Nk[1]* clusters[k].y[j]*(clusters[k].m_mu[j]^2 + clusters[k].s_sq_mu[j]) + clusters[k].Nk[1]* (clusters[k].m_nu[j]^2 + clusters[k].s_sq_nu[j]) ) #  
        end
    end
    return LB
end


"""
    Lp_pi_tau_phi(cells::Vector{CellFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters pi, tau, and phi.
"""
function Lp_pi_tau_phi(cells,conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0    
    for n in 1:N    
        t = cells[n].t #dataparams.LinearAddress[n][2]
        i = cells[n].i #dataparams.LinearAddress[n][1]
        for k in 1:K
            for tt in 1:t
                LB += cells[n].r[k] * cells[n].c[tt] * E_ln_pi(k,conditions[sum(T[1:i-1])+tt])
            end
        end
    end
    return LB
end

"""
    Lp_tau_omega(cells::Vector{CellFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}  
Calculates the log likelihood of parameters tau and omega.
"""
function Lp_tau_omega(cells,conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for n in 1:N
        t = cells[n].t #dataparams.LinearAddress[n][2]
        i = cells[n].i #dataparams.LinearAddress[n][1]
        for tt in 1:t
            minus_ln_omega_lb = 0.0
            minus_ln_omega_lb += recursive_cumsum_E_ln_minusomega(i,t,tt+1,T,conditions)
            LB += cells[n].c[tt] * (E_ln_omega(conditions[sum(T[1:i-1])+tt]) + minus_ln_omega_lb)
        end
    end
    return LB
end


"""
    Lp_pi_chi(clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;use_log=true) where {U <: AbstractFloat, W <: Int64}   
Calculates the log likelihood of parameters pi and chi.
"""
function Lp_pi_chi(clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;use_log=true) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    T_all = sum(T)
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for it in 1:T_all
        condition_i = conditions[it].i
        condition_t = conditions[it].t
        for k in 1:Kplus
            LB +=  modelparams.alpha0[condition_i][condition_t] * expectation_SBk(k,clusters,modelparams;use_log= use_log) * E_ln_pi(k,conditions[it])
        end
    end
    return LB
end


"""
    Lp_surrogate(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of the surrogate bound.
"""
function Lp_surrogate( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    alpha0s = recursive_flatten(modelparams.alpha0)
    Kplus = K + 1
    LB = 0.0
    for alpha0 in alpha0s
        LB += K * log(alpha0)
    end
    for k in 1:K
        LB += E_ln_chi(clusters[k]) + (K + 1 - k)*E_ln_minus_chi(clusters[k])
    end
    return LB
end

"""
    Lp_chi(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters chi.
"""
function Lp_chi(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        LB += -ln_Beta_distribution_normalizer(1.0,modelparams.gamma0[1]) + (1.0 - 1)*E_ln_chi(clusters[k]) + (modelparams.gamma0[1]- 1)*E_ln_minus_chi(clusters[k])
    end
    return LB
end


"""
    Lp_omega(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters omega.
""" 
function Lp_omega(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    T_all = sum(T)
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for it in 1:T_all
        condition_t = conditions[it].t
        if condition_t ==1
            continue
        end
        LB += -ln_Beta_distribution_normalizer(modelparams.phi1[1],modelparams.phi2[1]) + (modelparams.phi1[1] - 1)*E_ln_omega(conditions[it]) + (modelparams.phi2[1]- 1)*E_ln_minusomega(conditions[it])
    end
    return LB
end


"""
    Lp_sigma_sq(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;sigma_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters sigma squared.
"""
function Lp_sigma_sq(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;sigma_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    if sigma_update_mode == "Clusterwise"
        error("sigma_update_mode not implemented for 'Clusterwise'")
    elseif sigma_update_mode == "Genewise"
        for j in 1:J
            LB += ln_Gamma_distribution_normalizer(modelparams.xi1[1],modelparams.xi2[1]) - (modelparams.xi1[1] - 1)*E_ln_sigma_sq(j,clusters[1]) - modelparams.xi2[1]*E_one_over_sigma_sq(j,clusters[1])
        end
    elseif sigma_update_mode == "Global"
         error("sigma_update_mode not implemented for 'Global'")
    elseif sigma_update_mode == "Local"
        for k in 1:K
            for j in 1:J
                LB += ln_Gamma_distribution_normalizer(modelparams.xi1[1],modelparams.xi2[1]) - (modelparams.xi1[1] - 1)*E_ln_sigma_sq(j,clusters[k]) - modelparams.xi2[1]*E_one_over_sigma_sq(j,clusters[k])
            end
        end
    else
        error("Invalid sigma_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return LB
end

"""
    Lp_mu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}  
Calculates the log likelihood of parameters mu.
"""
function Lp_mu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        for j in 1:J
            LB += -0.5 * dataparams.logpi - 0.5 * E_ln_lambda(j,clusters[k]) - 0.5 * E_ln_sigma_sq(j,clusters[k]) - 0.5 * E_one_over_sigma_sq(j,clusters[k]) * E_mu_sq(j,clusters[k])
        end
    end
    return LB
end

"""
    Lp_nu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}  
Calculates the log likelihood of parameters nu.
"""
function Lp_nu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for j in 1:J
        LB += -0.5 * dataparams.logpi - 0.5 * log(modelparams.sigma_sq_nu[j]) - 0.5 * 1/(modelparams.sigma_sq_nu[j]) * E_nu_diff_sq(j,clusters[1],modelparams)
    end
    return LB
end


"""
    Lp_lambda(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;lambda_update_mode="Local")  where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters lambda.
"""
function Lp_lambda(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;lambda_update_mode="Local")  where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    if lambda_update_mode == "Clusterwise"
        for k in 1:K
            LB += ln_Gamma_distribution_normalizer(modelparams.kappa1[1],modelparams.kappa2[1]) - (modelparams.kappa1[1] - 1)*E_ln_lambda(1,clusters[k]) - modelparams.kappa2[1]*E_one_over_lambda(1,clusters[k])
        end
    elseif lambda_update_mode == "Genewise"
        for j in 1:J
            LB += ln_Gamma_distribution_normalizer(modelparams.kappa1[1],modelparams.kappa2[1]) - (modelparams.kappa1[1] - 1)*E_ln_lambda(j,clusters[1]) - modelparams.kappa2[1]*E_one_over_lambda(j,clusters[1])
        end
    elseif lambda_update_mode == "Global"
        LB += ln_Gamma_distribution_normalizer(modelparams.kappa1[1],modelparams.kappa2[1]) - (modelparams.kappa1[1] - 1)*E_ln_lambda(1,clusters[1]) - modelparams.kappa2[1]*E_one_over_lambda(1,clusters[1])
    elseif lambda_update_mode == "Local"
        for k in 1:K
            for j in 1:J
                LB += ln_Gamma_distribution_normalizer(modelparams.kappa1[1],modelparams.kappa2[1]) - (modelparams.kappa1[1] - 1)*E_ln_lambda(j,clusters[k]) - modelparams.kappa2[1]*E_one_over_lambda(j,clusters[k])
            end
        end
    else
        error("Invalid lambda_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return LB
end


"""
    Lp_rho( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters rho.
"""
function Lp_rho( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        for j in 1:J
            LB += xlogy(clusters[k].y[j], exp(E_ln_eta(j,clusters[k]))) + xlogy(1 - clusters[k].y[j], exp(E_ln_minus_eta(j,clusters[k])))
        end
    end
    return LB
end


"""
    Lp_eta( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;eta_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of parameters eta.
"""
function Lp_eta( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;eta_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    if eta_update_mode == "Clusterwise"
        for k in 1:K
            LB += -ln_Beta_distribution_normalizer(modelparams.varphi1[1],modelparams.varphi2[1]) + (modelparams.varphi1[1] - 1)*E_ln_eta(1,clusters[k]) + (modelparams.varphi2[1]- 1)*E_ln_minus_eta(1,clusters[k])
        end
    elseif  eta_update_mode == "Genewise"
        for j in 1:J
            LB += -ln_Beta_distribution_normalizer(modelparams.varphi1[1],modelparams.varphi2[1]) + (modelparams.varphi1[1] - 1)*E_ln_eta(j,clusters[1]) + (modelparams.varphi2[1]- 1)*E_ln_minus_eta(j,clusters[1])
        end
    elseif  eta_update_mode == "Global"
        LB += -ln_Beta_distribution_normalizer(modelparams.varphi1[1],modelparams.varphi2[1]) + (modelparams.varphi1[1] - 1)*E_ln_eta(1,clusters[1]) + (modelparams.varphi2[1]- 1)*E_ln_minus_eta(1,clusters[1])
    elseif eta_update_mode == "Local"
        for k in 1:K
            for j in 1:J
                LB += -ln_Beta_distribution_normalizer(modelparams.varphi1[1],modelparams.varphi2[1]) + (modelparams.varphi1[1] - 1)*E_ln_eta(j,clusters[k]) + (modelparams.varphi2[1]- 1)*E_ln_minus_eta(j,clusters[k])
            end
        end
    else
        error("Invalid eta_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return LB
end

"""
    Lq_r(cells::Vector{CellFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of the variational parameters r.
"""
function Lq_r(cells,dataparams::DataFeature,modelparams::ModelParameterFeature) # Flipping the sign for etropies so i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for n in 1:N
        LB  += -entropy(cells[n].r)
    end
    return LB
end
# ::Vector{CellFeature{U,W,J}}
"""
    Lq_c(cells::Vector{CellFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
Calculates the log likelihood of the variational parameters c.
"""
function Lq_c(cells,dataparams::DataFeature,modelparams::ModelParameterFeature) # Flipping the sign for etropies so i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    # LB = Threads.Atomic{float_type}(0.0)
    for n in 1:N
        # Threads.@threads for n in 1:N
        # Threads.atomic_add!(LB,-entropy(cells[n].c))
        LB+=-entropy(cells[n].c)
    end
    # return LB[]
    return LB
end


"""
    Lq_d(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}# Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters d.
"""
function Lq_d(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}# Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    T_all = sum(T)
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for it in 1:T_all
        sum_loggamma = 0.0
        loggamma_sum = 0.0
        E_Lq_d_sum = 0.0
        for k in 1:Kplus
            sum_loggamma += loggamma(conditions[it].d[k])
            loggamma_sum += conditions[it].d[k]
            E_Lq_d_sum += (conditions[it].d[k] - 1) * E_ln_pi(k,conditions[it])
        end
        LB += (sum_loggamma - loggamma(loggamma_sum)) + E_Lq_d_sum
    end
    return LB
end

"""
    Lq_w1w2(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters w1 and w2.
"""
function Lq_w1w2(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    T_all = sum(T)
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for it in 1:T_all
        condition_t = conditions[it].t
        if condition_t ==1
            continue
        end
        LB += -ln_Beta_distribution_normalizer(conditions[it].w1[1],conditions[it].w2[1]) + (conditions[it].w1[1] - 1)*E_ln_omega(conditions[it]) + (conditions[it].w2[1] - 1)*E_ln_minusomega(conditions[it])
    end
    return LB
end

"""
    Lq_g1g2(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters g1 and g2.
"""
function Lq_g1g2(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        alpha_part = clusters[k].g1[1]*clusters[k].g2[1]
        beta_part = (1 - clusters[k].g1[1])*clusters[k].g2[1]
        LB += -ln_Beta_distribution_normalizer(alpha_part,beta_part) + (beta_part - 1)*E_ln_chi(clusters[k]) + (beta_part - 1)*E_ln_minus_chi(clusters[k])
    end
    return LB
end

"""
    Lq_ab(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; sigma_update_mode="Local") where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters a and b.
"""
function Lq_ab(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; sigma_update_mode="Local") where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    if sigma_update_mode == "Clusterwise"
        error("sigma_update_mode not implemented for 'Clusterwise'")
    elseif sigma_update_mode == "Genewise"
        for j in 1:J
            LB += ln_Gamma_distribution_normalizer(clusters[1].a[j],clusters[1].b[j]) - (clusters[1].a[j] - 1)*E_ln_sigma_sq(j,clusters[1]) - clusters[1].b[j]*E_one_over_sigma_sq(j,clusters[1])
        end
    elseif sigma_update_mode == "Global"
         error("sigma_update_mode not implemented for 'Global'")
    elseif sigma_update_mode == "Local"
        for k in 1:K
            for j in 1:J
                LB += ln_Gamma_distribution_normalizer(clusters[k].a[j],clusters[k].b[j]) - (clusters[k].a[j] - 1)*E_ln_sigma_sq(j,clusters[k]) - clusters[k].b[j]*E_one_over_sigma_sq(j,clusters[k])
            end
        end
    else
        error("Invalid sigma_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return LB
end

"""
    Lq_ms_sq_mu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}# Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters m and s squared for mu.
"""
function Lq_ms_sq_mu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}# Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        for j in 1:J
            LB += -0.5 * dataparams.logpi - 0.5 * log(clusters[k].s_sq_mu[j]) - 0.5 * 1/clusters[k].s_sq_mu[j] * E_prior_sq_diff_mu(j,clusters[k])
        end
    end
    return LB
end

"""
    Lq_ms_sq_nu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}# Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters m and s squared for nu.
"""
function Lq_ms_sq_nu(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}# Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for j in 1:J
        LB += -0.5 * dataparams.logpi - 0.5 * log(clusters[1].s_sq_nu[j]) - 0.5 * 1/clusters[1].s_sq_nu[j] * E_prior_sq_diff_nu(j,clusters[1])
    end
    return LB
end

"""
    Lq_uv(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; lambda_update_mode="Local")  where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters u and v.
"""
function Lq_uv(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;  lambda_update_mode="Local")  where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    if lambda_update_mode == "Clusterwise"
        for k in 1:K
            LB += ln_Gamma_distribution_normalizer(clusters[k].u[1],clusters[k].v[1]) - (clusters[k].u[1] - 1)*E_ln_lambda(1,clusters[k]) - clusters[k].v[1]*E_one_over_lambda(1,clusters[k])
        end
    elseif lambda_update_mode == "Genewise"
        for j in 1:J
            LB += ln_Gamma_distribution_normalizer(clusters[1].u[j],clusters[1].v[j]) - (clusters[1].u[j] - 1)*E_ln_lambda(j,clusters[1]) - clusters[1].v[j]*E_one_over_lambda(j,clusters[1])
        end
    elseif lambda_update_mode == "Global"
        LB += ln_Gamma_distribution_normalizer(clusters[1].u[1],clusters[1].v[1]) - (clusters[1].u[1] - 1)*E_ln_lambda(1,clusters[1]) - clusters[1].v[1]*E_one_over_lambda(1,clusters[1])
    elseif lambda_update_mode == "Local"
        for k in 1:K
            for j in 1:J
                LB += ln_Gamma_distribution_normalizer(clusters[k].u[j],clusters[k].v[j]) - (clusters[k].u[j] - 1)*E_ln_lambda(j,clusters[k]) - clusters[k].v[j]*E_one_over_lambda(j,clusters[k])
            end
        end
    else
        error("Invalid lambda_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return LB
end



"""
    Lq_y(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters y.
"""
function Lq_y(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    for k in 1:K
        for j in 1:J
            LB += xlogy(clusters[k].y[j], clusters[k].y[j]) + xlogy(1 - clusters[k].y[j], 1 - clusters[k].y[j])
        end
    end
    return LB
end

"""
    Lq_h1h2(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;eta_update_mode="Local") where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
Calculates the log likelihood of the variational parameters h1 and h2.
"""
function Lq_h1h2(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;eta_update_mode="Local") where {U <: AbstractFloat, W <: Int64} # Keep signs the same so that i can subtract from elbo
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    LB = 0.0
    if eta_update_mode == "Clusterwise"
        for k in 1:K
            LB += -ln_Beta_distribution_normalizer(clusters[k].h1[1],clusters[k].h2[1]) + (clusters[k].h1[1] - 1)*E_ln_eta(1,clusters[k]) + (clusters[k].h2[1]- 1)*E_ln_minus_eta(1,clusters[k])
        end
    elseif  eta_update_mode == "Genewise"
        for j in 1:J
            LB += -ln_Beta_distribution_normalizer(clusters[1].h1[j],clusters[1].h2[j]) + (clusters[1].h1[j] - 1)*E_ln_eta(j,clusters[1]) + (clusters[1].h2[j]- 1)*E_ln_minus_eta(j,clusters[1])
        end
    elseif  eta_update_mode == "Global"
        LB += -ln_Beta_distribution_normalizer(clusters[1].h1[1],clusters[1].h2[1]) + (clusters[1].h1[1] - 1)*E_ln_eta(1,clusters[1]) + (clusters[1].h2[1]- 1)*E_ln_minus_eta(1,clusters[1])
    elseif eta_update_mode == "Local"
        for k in 1:K
            for j in 1:J
                LB += -ln_Beta_distribution_normalizer(clusters[k].h1[1],clusters[k].h2[1]) + (clusters[k].h1[1] - 1)*E_ln_eta(j,clusters[k]) + (clusters[k].h2[1]- 1)*E_ln_minus_eta(j,clusters[k])
            end
        end
    else
        error("Invalid eta_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end

    # for k in 1:K
    #     LB += -ln_Beta_distribution_normalizer(clusters[k].h1[1],clusters[k].h2[1]) + (clusters[k].h1[1] - 1)*E_ln_eta(clusters[k]) + (clusters[k].h2[1]- 1)*E_ln_minus_eta(clusters[k])
    # end
    return LB
end


"""
    ELBO(cells::Vector{CellFeature{U,W}},clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; use_log::Bool=true,update_clusterwise::Bool = false,eta_update_mode="Local", lambda_update_mode="Local", sigma_update_mode = "Local") where {U <: AbstractFloat, W <: Int64}
Calculates the Evidence Lower Bound (ELBO) for the entire model.            
"""
function ELBO(cells,clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; use_log::Bool=true,update_clusterwise::Bool = false,eta_update_mode="Local", lambda_update_mode="Local", sigma_update_mode = "Local") where {U <: AbstractFloat, W <: Int64}
    LB = 0.0
    LB += Lp_data(clusters,dataparams,modelparams) 
    LB += Lp_pi_tau_phi(cells,conditions,dataparams,modelparams) 
    LB += Lp_tau_omega(cells,conditions,dataparams,modelparams) 
    LB += Lp_pi_chi(clusters,conditions,dataparams,modelparams;use_log = use_log) 
    LB += Lp_surrogate(clusters,dataparams,modelparams) 
    LB += Lp_chi(clusters,dataparams,modelparams) 
    LB += Lp_omega(conditions,dataparams,modelparams) 
    LB += Lp_sigma_sq(clusters,dataparams,modelparams; sigma_update_mode = sigma_update_mode ) 
    LB += Lp_mu(clusters,dataparams,modelparams) 
    LB += Lp_nu(clusters,dataparams,modelparams)
    LB += Lp_lambda(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode) 
    LB += Lp_rho(clusters,dataparams,modelparams) 
    LB += Lp_eta(clusters,dataparams,modelparams;eta_update_mode=eta_update_mode) 
    LB -= Lq_r(cells,dataparams,modelparams) 
    LB -= Lq_c(cells,dataparams,modelparams) 
    LB -= Lq_d(conditions,dataparams,modelparams) 
    LB -= Lq_w1w2(conditions,dataparams,modelparams) 
    LB -= Lq_g1g2(clusters,dataparams,modelparams) 
    LB -= Lq_ab(clusters,dataparams,modelparams; sigma_update_mode = sigma_update_mode ) 
    LB -= Lq_ms_sq_mu(clusters,dataparams,modelparams) 
    LB -= Lq_ms_sq_nu(clusters,dataparams,modelparams)
    LB -= Lq_uv(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode) 
    LB -= Lq_y(clusters,dataparams,modelparams) 
    LB -= Lq_h1h2(clusters,dataparams,modelparams;eta_update_mode=eta_update_mode)
    return LB
end

"""
    calc_SurragateLowerBound(rho_hat,omega_hat,T,γ,α0,Tk)
 Calculates the Surrogate Bound
"""
function calc_SurragateLowerBound(rho_hat,omega_hat,T,γ,α0,Tk)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K]
    k_vec = collect(1:K)
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α0 .* e_βk .* Tk
    lb_lg = sum(lb_lg_k)
    return lb_lg
end

"""
    calc_SurragateLowerBound_unconstrained(c,d,T,γ,α0,Tk)
Calculates the unconstrained Surrogate Bound 
"""
function calc_SurragateLowerBound_unconstrained(c,d,T,γ,α0,Tk)
    rho_hat = sigmoid.(c)
    omega_hat = exp.(d)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α0 .* e_βk .* Tk[1:K]
    lb_lg = sum(lb_lg_k)
    return lb_lg
end


"""
    c_Ga(a0, b0)
Calculates the log of the Gamma function
"""
function c_Ga(a0, b0)
    a0 .* log.(b0) .- loggamma.(a0)
end

"""
    c_Beta(a0, b0)
Calculates the log of the Beta function
"""
function c_Beta(a0, b0)
    - logbeta.(a0,b0) 
end

"""
    calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
Calculates the current iterations elbo. Wrapper function that calls the component elbo calculation functions
"""
function calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
    K = modelparams.K
    for k in 1:K
        elbolog.per_k_elbo[k,iter] = 0.0
    end
    elbo_val = 0.0
    dataelbo,elbolog = calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    # dataelbo = 0.0
    zentropy = calc_Hz_fast3(cellpop,clusters,dataparams)
    lg_elbo,elbolog =calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    # lg_elbo = 0.0
    w_elbo = calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    sentropy = calc_HsElbo(conditionparams, dataparams, modelparams)
    elbo_val +=  dataelbo + zentropy + lg_elbo + w_elbo + sentropy
    for k in 1:K
        elbolog.per_k_elbo[k,iter] += zentropy + w_elbo + sentropy
    end
    return elbo_val,elbolog
end


"""
    calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
Calculates the current iterations data elbo.
"""
function calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    one_half_const = 1/2
    data_elbo = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        clusters[k].cache .= 0.0
        yjk_entropy_perK = 0.0
        perK_data_elbo = 0.0
        for j in 1:G
            perK_data_elbo += -1*one_half_const * clusters[k].Nk[1] * log(2π)
            perK_data_elbo += -1*one_half_const * clusters[k].Nk[1]  * log(clusters[k].σ_sq_k_hat[j])
            perK_data_elbo += -1*one_half_const * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat_sq[j] 
            perK_data_elbo +=  1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat[j] * clusters[k].κk_hat[j]
            perK_data_elbo +=  -1*one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].var_muk[j]
            perK_data_elbo += -1* one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 - one_half_const * clusters[k].yjk_hat[j] * log(geneparams[j].λ_sq[1]) 
            perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].var_muk[j] 
            perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 
            perK_data_elbo += clusters[k].yjk_hat[j] * log(modelparams.ηk[1]) + (1 - clusters[k].yjk_hat[j]) * log((1-modelparams.ηk[1]))
            perK_data_elbo += one_half_const * clusters[k].yjk_hat[j] * log(clusters[k].v_sq_k_hat[j])
            perK_data_elbo +=  one_half_const * clusters[k].yjk_hat[j]
        end
    

        yjk_entropy_perK += entropy(clusters[k].yjk_hat)
        yjk_entropy_perK = -yjk_entropy_perK
        perK_ebloval =   perK_data_elbo + yjk_entropy_perK
        elbolog.per_k_elbo[k,iter] += perK_ebloval
        data_elbo += perK_ebloval
    end
    return data_elbo,elbolog
end

"""
    calc_Hz_fast3(cellpop,clusters,dataparams)
Fast calculation of the current iterations cell cluster assignment entropy.
"""
function calc_Hz_fast3(cellpop,clusters,dataparams)
    float_type = dataparams.BitType
    N = dataparams.N
    z_entropy = 0.0
    @fastmath @inbounds @simd for i in 1:N
        z_entropy += entropy(cellpop[i].rtik)
    end
    return  -z_entropy
end

"""
    calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
Calculates the current iterations surrogate priors' elbo.
"""
function calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0


    lb_lg = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        perK_lg_ebloval = 0.0
        c_B = beta(clusters[k].gk_hat[1] * clusters[k].hk_hat[1] , (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1])
        e_βk = expectation_βk(k,clusters,modelparams)
        e_logUk = expectation_logUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
        e_log1minusUk = expectation_log1minusUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
        perK_lg_ebloval += -1.0 * c_B +  (T + 1. - clusters[k].gk_hat[1] * clusters[k].hk_hat[1]) *e_logUk  +  (T * ( K + 1. - k) + γ0 - (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1]) * e_log1minusUk  +  e_βk * α0 * Tk[k]
        elbolog.per_k_elbo[k,iter] += perK_lg_ebloval
        lb_lg += perK_lg_ebloval
    end


    return lb_lg,elbolog
end

"""
    calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
Calculates the current iterations dynamic priors' elbo.
"""
function calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    T = dataparams.T
    wAlloc_elbo = 0.0



    @fastmath @inbounds for t in 2:T
        b_cttprime = 0.0
        @fastmath @inbounds for t_prime_b in t:T
            @fastmath @inbounds @simd for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                b_cttprime += conditionparams[t_prime_b].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        c_Beta_p = -logbeta(1, modelparams.ϕ0)
        c_Beta_q = -(-logbeta(1, conditionparams[t-1].st_hat[1]))
        c_Beta_pq = c_Beta_p + c_Beta_q
        ϕ0_st_hat_vec = modelparams.ϕ0- conditionparams[t-1].st_hat[1]
        e_log_tilde_wt = expectation_log_tilde_wtt(1, conditionparams[t-1].st_hat[1])
        e_log_minus_tilde_wt = expectation_log_minus_tilde_wtt(1, conditionparams[t-1].st_hat[1])


        wAlloc_elbo_t = c_Beta_pq + (1) * e_log_tilde_wt + (ϕ0_st_hat_vec + b_cttprime) * e_log_minus_tilde_wt
        wAlloc_elbo += wAlloc_elbo_t
    end
    return wAlloc_elbo
end

"""
    calc_HsElbo(conditionparams, dataparams, modelparams)
Calculation of the current iterations conditions assignment entropy.
"""
function calc_HsElbo(conditionparams, dataparams, modelparams)
    float_type = dataparams.BitType
    T = dataparams.T


    s_entropy = 0.0
    @fastmath @inbounds @simd for t in 1:T
        s_entropy += entropy(conditionparams[t].c_tt_prime)
    end


    s_entropy = -s_entropy

    return s_entropy
end


