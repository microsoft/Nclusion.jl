
############################################
###### SIMPLE EXPECTATION FUNCTIONS ########
############################################
function e_ln_sigma_sq(a::AbstractFloat,b::AbstractFloat)
    return log(b) - digamma(a)
end

function e_one_over_sigma_sq(a::AbstractFloat,b::AbstractFloat)
    return a/b
end

function e_ll_sq_diff_mu(x_sq::AbstractFloat,x::AbstractFloat,m_nu::AbstractFloat,m_mu::AbstractFloat,y::AbstractFloat,s_sq_mu::AbstractFloat,s_sq_nu::AbstractFloat)
    return x_sq - 2*y*x*m_mu - 2*x*m_nu + 2*y*m_mu*m_nu + y*m_mu^2 + y*s_sq_mu + s_sq_nu + m_nu^2
    # return x_sq - 2*x*m_mu - 2*x*m_nu + 2*m_mu*m_nu + m_mu^2 + s_sq_mu + s_sq_nu + m_nu^2
end

function e_ln_pi(d::AbstractFloat,dsum::AbstractFloat)
    return digamma(d) - digamma(dsum)
end

function e_ln_omega(w1::AbstractFloat,w2::AbstractFloat)
    return digamma(w1) - digamma(w1 + w2)
end

function e_ln_minusomega(w1::AbstractFloat,w2::AbstractFloat)
    return digamma(w2) - digamma(w1 + w2)
end

function e_omega(w1::AbstractFloat,w2::AbstractFloat)
    return w1 /(w1 + w2)
end

function e_minusomega(w1::AbstractFloat,w2::AbstractFloat)
    return w2 /(w1 + w2)
end

function e_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return alpha_part/(alpha_part + beta_part)
end

function e_minus_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return beta_part/(alpha_part + beta_part)
end

function e_ln_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return digamma(alpha_part) - digamma(alpha_part + beta_part)
end

function e_ln_minus_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return digamma(beta_part) - digamma(alpha_part + beta_part)
end

function e_ln_lambda(u::AbstractFloat,v::AbstractFloat)
    return log(v) -  digamma(u)
end

function e_one_over_lambda(u::AbstractFloat,v::AbstractFloat)
    return u/v
end

function e_mu_sq(m_mu::AbstractFloat,s_sq_mu::AbstractFloat)
    return m_mu^2 + s_sq_mu
end

function e_mu(m_mu::AbstractFloat)
    return m_mu
end

function e_prior_sq_diff_mu(m_mu::AbstractFloat,s_sq_mu::AbstractFloat)
    return m_mu^2 -2*m_mu*m_mu + s_sq_mu + m_mu^2
end

function e_nu_sq(m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
    return m_nu^2 + s_sq_nu
end

function e_nu_diff_sq(m_nu::AbstractFloat,s_sq_nu::AbstractFloat,nu0::AbstractFloat)
    return m_nu^2 + s_sq_nu -2*m_nu*nu0 + nu0^2
end

function e_nu(m_nu::AbstractFloat)
    return m_nu
end

function e_prior_sq_diff_nu(m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
    return m_nu^2 -2*m_nu*m_nu + s_sq_nu + m_nu^2
end

function e_ln_eta(h1::AbstractFloat,h2::AbstractFloat)
    return digamma(h1) - digamma(h1 + h2)
end

function e_ln_minus_eta(h1::AbstractFloat,h2::AbstractFloat)
    return digamma(h2) - digamma(h1 + h2)
end

function e_eta(h1::AbstractFloat,h2::AbstractFloat)
    return h1/(h1 + h2)
end

function e_minus_eta(h1::AbstractFloat,h2::AbstractFloat)
    return h2/(h1 + h2)
end

"""
"""
function recursive_minus_e_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,g1::Vector{U},g2::Vector{U}) where {U <: AbstractFloat}# formerly recursive_minus_e_uk_cumprod
    if iszero(k)
        return cummulative_prod
    else
        cummulative_prod *= e_minus_chi(g1[k],g2[k])
        k -= 1
        recursive_minus_e_chi_cumprod(k,cummulative_prod,g1,g2)
    end
end

"""
"""
function log_of_recursive_minus_e_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,g1::Vector{U},g2::Vector{U})  where {U <: AbstractFloat} #  Logging for stability?
    if iszero(k)
        return cummulative_prod
    else
        cummulative_prod += log(e_minus_chi(g1[k],g2[k]))
        k -= 1
        log_of_recursive_minus_e_chi_cumprod(k,cummulative_prod,g1,g2)
    end
end

"""
"""
function expectation_sbk(k::Int,K::Int,g1::Vector{U},g2::Vector{U};use_log=false)  where {U <: AbstractFloat}# formerly expectation_βk
    Kplus = K + 1
    if k == Kplus
        e_chi_k = 1.0
    else
        e_chi_k = e_chi(g1[k],g2[k])
    end
    if isone(k)
        cumprod_e_minus_chi_k = 1.0
        if use_log
            cumprod_e_minus_chi_k = log(cumprod_e_minus_chi_k)          
        end
    else
        if use_log
            cumprod_e_minus_chi_k = log_of_recursive_minus_e_chi_cumprod(k-1,0.0,g1,g2)
        else
            cumprod_e_minus_chi_k = recursive_minus_e_chi_cumprod(k-1,1.0,g1,g2)
        end
    end
    # println("($e_uk,$cumprod_minus_e_uk)")
    if use_log
        e_sbk = e_chi_k * exp(cumprod_e_minus_chi_k)
    else
        e_sbk = e_chi_k * cumprod_e_minus_chi_k
    end
    return e_sbk
end

############################################
############################################
############################################

TATION FUNCTIONS ########
#################################################

function E_ln_sigma_sq(j,clusterfeature::ClusterFeature)
    return e_ln_sigma_sq(clusterfeature.a[j],clusterfeature.b[j])
end

function E_one_over_sigma_sq(j,clusterfeature::ClusterFeature)
    return e_one_over_sigma_sq(clusterfeature.a[j],clusterfeature.b[j])
end


function E_ll_sq_diff_mu(j,cellfeature::CellFeature,clusterfeature::ClusterFeature)
    return e_ll_sq_diff_mu(cellfeature.xsq[j],cellfeature.x[j],clusterfeature.m_mu[j],clusterfeature.y[j],clusterfeature.s_sq_mu[j],clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j])
end

function E_ln_pi(k,conditionfeature::ConditionFeature)
    return e_ln_pi(conditionfeature.d[k],conditionfeature.d_sum[1])
end

function E_ln_omega(conditionfeature::ConditionFeature)
    return e_ln_omega(conditionfeature.w1[1],conditionfeature.w2[1])
end

function E_ln_minusomega(conditionfeature::ConditionFeature)
    return e_ln_minusomega(conditionfeature.w1[1],conditionfeature.w2[1])
end

function E_omega(conditionfeature::ConditionFeature)
    return e_omega(conditionfeature.w1[1],conditionfeature.w2[1])
end

function E_minusomega(conditionfeature::ConditionFeature)
    return e_minusomega(conditionfeature.w1[1],conditionfeature.w2[1])
end

function E_chi(clusterfeature::ClusterFeature)
    return e_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

function E_minus_chi(clusterfeature::ClusterFeature)
    return e_minus_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

function E_ln_chi(clusterfeature::ClusterFeature)
    return e_ln_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

function E_ln_minus_chi(clusterfeature::ClusterFeature)
    return e_ln_minus_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

function E_ln_lambda(j,clusterfeature::ClusterFeature)
    return e_ln_lambda(clusterfeature.u[j],clusterfeature.v[j])
end

function E_one_over_lambda(j,clusterfeature::ClusterFeature)
    return e_one_over_lambda(clusterfeature.u[j],clusterfeature.v[j])
end

function E_mu_sq(j,clusterfeature::ClusterFeature)
    return e_mu_sq(clusterfeature.m_mu[j],clusterfeature.s_sq_mu[j])
end

function E_mu(j,clusterfeature::ClusterFeature)
    return e_mu(clusterfeature.m_mu[j])
end

function E_prior_sq_diff_mu(j,clusterfeature::ClusterFeature)
    return e_prior_sq_diff_mu(clusterfeature.m_mu[j],clusterfeature.s_sq_mu[j])
end

function E_nu_sq(j,clusterfeature::ClusterFeature)
    return e_nu_sq(clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j])
end

function E_nu(j,clusterfeature::ClusterFeature)
    return e_nu(clusterfeature.m_nu[j])
end

function E_prior_sq_diff_nu(j,clusterfeature::ClusterFeature)
    return e_prior_sq_diff_nu(clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j])
end

function E_nu_diff_sq(j,clusterfeature::ClusterFeature,modelparams::ModelParameterFeature)
    return e_nu_diff_sq(clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j],modelparams.nu0[j])
end

function E_ln_eta(j,clusterfeature::ClusterFeature)
    return e_ln_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

function E_ln_minus_eta(j,clusterfeature::ClusterFeature)
    return e_ln_minus_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

function E_eta(j,clusterfeature::ClusterFeature)
    return e_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

function E_minus_eta(j,clusterfeature::ClusterFeature)
    return e_minus_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

#copy_cells = deepcopy(cells); copy_precomputed_genefeatures_cells = deepcopy(cells); for n in 1:dataparams.N copy_cells[n].cache .= 0.0; copy_precomputed_genefeatures_cells[n].cache .= 0.0; adjust_E_ln_pi!(copy_cells[n],conditions,dataparams);adjust_E_ln_pi!(copy_precomputed_genefeatures_cells[n],conditions,dataparams); for k in 1:modelparams.K E_log_normal_l_j!(copy_cells[n],clusters[k], dataparams); E_log_normal_l_j!(preupdated_genefeatures,copy_precomputed_genefeatures_cells[n],clusters[k], dataparams); end; end; copy_cells_z_argmax = [argmax(copy_cells[n].cache) for n in 1:dataparams.N]; copy_precomputed_genefeatures_cells_z_argmax = [argmax(copy_precomputed_genefeatures_cells[n].cache) for n in 1:dataparams.N];println(Clustering.randindex(cell_cluster_labels,cell_cluster_labels)[1]); println(countmap(copy_cells_z_argmax)); println(Clustering.randindex(cell_cluster_labels,copy_cells_z_argmax)[1]); println(countmap(copy_precomputed_genefeatures_cells_z_argmax)); println(Clustering.randindex(cell_cluster_labels,copy_precomputed_genefeatures_cells_z_argmax)[1]); println(all([el in collect(1:dataparams.N)[cell_cluster_labels .!= copy_cells_z_argmax]  for el in collect(1:dataparams.N)[cell_cluster_labels .!= copy_precomputed_genefeatures_cells_z_argmax]]));

"""
"""
function recursive_minus_E_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,clusters::Vector{ClusterFeature{U,W}}) where {U <: AbstractFloat, W <: Int64} # formerly recursive_minus_e_uk_cumprod
    if iszero(k)
        return cummulative_prod
    else
        cummulative_prod *= E_minus_chi(clusters[k])
        k -= 1
        recursive_minus_E_chi_cumprod(k,cummulative_prod,clusters)
    end
end

"""
"""
function log_of_recursive_minus_E_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,clusters::Vector{ClusterFeature{U,W}}) where {U <: AbstractFloat, W <: Int64} #  Logging for stability?
    if iszero(k)
        return cummulative_prod
    else
        cummulative_prod += log(E_minus_chi(clusters[k]))
        k -= 1
        log_of_recursive_minus_E_chi_cumprod(k,cummulative_prod,clusters)
    end
end

"""
"""
function expectation_SBk(k::Int,clusters::Vector{ClusterFeature{U,W}},modelparams::ModelParameterFeature;use_log=false) where {U <: AbstractFloat, W <: Int64} # formerly expectation_βk
    K = modelparams.K
    Kplus = K + 1
    if k == Kplus
        e_chi_k = 1.0
    else
        e_chi_k = E_chi(clusters[k])
    end
    if isone(k)
        cumprod_e_minus_chi_k = 1.0
        if use_log
            cumprod_e_minus_chi_k = log(cumprod_e_minus_chi_k)          
        end
    else
        if use_log
            cumprod_e_minus_chi_k = log_of_recursive_minus_E_chi_cumprod(k-1,0.0,clusters)
        else
            cumprod_e_minus_chi_k = recursive_minus_E_chi_cumprod(k-1,1.0,clusters)
        end
    end
    # println("($e_uk,$cumprod_minus_e_uk)")
    if use_log
        e_SBk = e_chi_k * exp(cumprod_e_minus_chi_k)
    else
        e_SBk = e_chi_k * cumprod_e_minus_chi_k
    end
    return e_SBk
end


"""
"""
function βk_expected_value(rho_hat_vec, omega_hat_vec)
    K = length(rho_hat_vec)
    e_uk_vec = uk_expected_value(rho_hat_vec, omega_hat_vec)
    minus_e_uk_vec = 1. .- e_uk_vec
    cumprod_minus_e_uk_vec= cumprod(minus_e_uk_vec)
    app_e_uk_vec= deepcopy(e_uk_vec)
    app_cumprod_minus_e_uk_vec = deepcopy(cumprod_minus_e_uk_vec)
    append!(app_e_uk_vec,1.0) 
    insert!(app_cumprod_minus_e_uk_vec,1,1)
    e_βk_vec = Vector{Float64}(undef,K+1)
    for k in 1:K+1
        e_βk_vec[k] = app_e_uk_vec[k] * app_cumprod_minus_e_uk_vec[k]
    end
    return e_βk_vec
end

"""
"""
function uk_expected_value(rho_hat_vec, omega_hat_vec)
    minus_rho_hat_vec = 1 .- rho_hat_vec
    rho_omega_hat_vec = rho_hat_vec .* omega_hat_vec
    minus_rho_omega_hat_vec = minus_rho_hat_vec .*omega_hat_vec
    e_uk_vec = rho_omega_hat_vec ./ (rho_omega_hat_vec .+ minus_rho_omega_hat_vec)
    return e_uk_vec
end

"""
"""
function logUk_expected_value(rho_hat,omega_hat)
    return digamma.(rho_hat .* omega_hat) .- digamma.(omega_hat)
end

"""
"""
function log1minusUk_expected_value(rho_hat,omega_hat)
    return digamma.((1.0 .- rho_hat) .* omega_hat) .- digamma.(omega_hat)
end
#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################

"""
"""
function adjust_e_log_π_tk3!(t,conditionparams)
    Kplus = length(conditionparams[1].d_hat_t)
    # conditions_pis_sums = 0.0
    # e_log_π_t_cache = Vector{Float64}(undef,Kplus)
    for k in 1:Kplus
        conditions_pis_sums = 0.0
        for tt in 1:t
            conditions_pis_sums += conditionparams[t].c_tt_prime[tt] * (digamma(conditionparams[tt].d_hat_t[k]) - digamma(conditionparams[tt].d_hat_t_sum[1]))
        end
        conditionparams[t].e_log_π_t_cache[k] = conditions_pis_sums
    end
    return conditionparams
end

"""
"""
function log_π_expected_value_fast3!(conditionparams,dataparams,modelparams)
    T = dataparams.T
    Kplus = modelparams.K + 1
    for t in 1:T
        for k in 1:Kplus
            conditionparams[t].e_log_π_t_cache[k] = expectation_log_π_tk(conditionparams[t].d_hat_t[k], conditionparams[t].d_hat_t_sum[1] )
        end
    end
    return conditionparams
end

"""
"""
function expectation_βk(k,clusters,modelparams)
    K = modelparams.K
    Kplus = K + 1
    if k == Kplus
        e_uk = 1.0
    else
        e_uk = expectation_uk(clusters[k].gk_hat[1], clusters[k].hk_hat[1])
    end
    if isone(k)
        cumprod_minus_e_uk = 1.0
    else
        cumprod_minus_e_uk = recursive_minus_e_uk_cumprod(k-1,1.0,clusters)
    end
    # println("($e_uk,$cumprod_minus_e_uk)")
    e_βk = e_uk * cumprod_minus_e_uk
    return e_βk
end

"""
"""
function recursive_minus_e_uk_cumprod(k,cummulative_prod,clusters)
    if iszero(k)
        return cummulative_prod
    else
        cummulative_prod *= 1-expectation_uk(clusters[k].gk_hat[1], clusters[k].hk_hat[1])
        k -= 1
        recursive_minus_e_uk_cumprod(k,cummulative_prod,clusters)
    end
end

"""
"""
function expectation_log_normal_l_j!(cluster,cell, dataparams)
    G = dataparams.G
    for j in 1:G
        cell.cache[j]  =  -0.5 * (log(cluster.σ_sq_k_hat[j]) + 1/cluster.σ_sq_k_hat[j] * ( cell.xsq[j] - 2 * cell.x[j] * cluster.κk_hat[j] + cluster.yjk_hat[j] * (cluster.mk_hat[j] ^2 + cluster.v_sq_k_hat[j])  ))
        # cell.cache[j]  =  -0.5 * (log(cluster.σ_sq_k_hat[j]) + 1/cluster.σ_sq_k_hat[j] * ( cell.xsq[j] - 2 * cell.x[j] * (cluster.yjk_hat[j] * cluster.mk_hat[j]) + cluster.yjk_hat[j] * (cluster.mk_hat[j] ^2 + cluster.v_sq_k_hat[j])  ))
    end
    return cluster
end

"""
"""
function expectation_log_π_tk(θ_hat_tk,θ_hat_t_sum)
    e_log_π_tk =  digamma(θ_hat_tk) - digamma(θ_hat_t_sum)
    return e_log_π_tk 
end

"""
"""
function expectation_log_tilde_wtt(awt_hat,bwt_hat)
    return  digamma(awt_hat) - digamma(awt_hat + bwt_hat)
end

"""
"""
function expectation_log_minus_tilde_wtt(awt_hat,bwt_hat)
    return  digamma(bwt_hat) - digamma(awt_hat + bwt_hat)
end

"""
"""
function expectation_uk(rho_hat, omega_hat)
    # minus_rho_hat = 1 - rho_hat
    # rho_omega_hat = rho_hat * omega_hat
    # minus_rho_omega_hat = (1 - rho_hat) * omega_hat
    e_uk = (rho_hat * omega_hat) / ((rho_hat * omega_hat)+ ((1 - rho_hat) * omega_hat))
    return e_uk
end

"""
"""
function expectation_αt(a_α,b_α)
    return a_α / b_α
end

"""
"""
function expectation_log_αt(a_α,b_α)
    return digamma(a_α) - log(b_α)
end

"""
"""
function expectation_logUk(rho_hat,omega_hat)
    return digamma(rho_hat * omega_hat) - digamma(omega_hat)
end

"""
"""
function expectation_log1minusUk(rho_hat,omega_hat)
    return digamma((1.0 - rho_hat) * omega_hat) - digamma(omega_hat)
end



