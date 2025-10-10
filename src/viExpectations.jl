
############################################
###### SIMPLE EXPECTATION FUNCTIONS ########
############################################
"""
    e_ln_sigma_sq(a::AbstractFloat,b::AbstractFloat)
Expected value of the logarithm of a variable following a Gamma distribution with shape `a` and rate `b`.
"""
function e_ln_sigma_sq(a::AbstractFloat,b::AbstractFloat)
    return log(b) - digamma(a)
end

"""
    e_one_over_sigma_sq(a::AbstractFloat,b::AbstractFloat)
Expected value of the inverse of a variable following a Gamma distribution with shape `a` and rate `b`.
"""
function e_one_over_sigma_sq(a::AbstractFloat,b::AbstractFloat)
    return a/b
end

"""
    e_ll_sq_diff_mu(x_sq::AbstractFloat,x::AbstractFloat,m_mu::AbstractFloat,y::AbstractFloat,s_sq_mu::AbstractFloat,m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
Expected value of the squared difference between a data point `x` and the sum of two variables `mu` and `nu`, where `mu` follows a Normal distribution with mean `m_mu` and variance `s_sq_mu`, and `nu` follows a Normal distribution with mean `m_nu` and variance `s_sq_nu`. The variable `y` is a scaling factor.
"""
function e_ll_sq_diff_mu(x_sq::AbstractFloat,x::AbstractFloat,m_nu::AbstractFloat,m_mu::AbstractFloat,y::AbstractFloat,s_sq_mu::AbstractFloat,s_sq_nu::AbstractFloat)
    return x_sq - 2*y*x*m_mu - 2*x*m_nu + 2*y*m_mu*m_nu + y*m_mu^2 + y*s_sq_mu + s_sq_nu + m_nu^2
    # return x_sq - 2*x*m_mu - 2*x*m_nu + 2*m_mu*m_nu + m_mu^2 + s_sq_mu + s_sq_nu + m_nu^2
end

"""
    e_ln_pi(d::AbstractFloat,dsum::AbstractFloat)
Expected value of the logarithm of a variable following a Dirichlet distribution with parameter `d` and sum of parameters `dsum`.
"""
function e_ln_pi(d::AbstractFloat,dsum::AbstractFloat)
    return digamma(d) - digamma(dsum)
end

"""
    e_ln_omega(w1::AbstractFloat,w2::AbstractFloat)
Expected value of the logarithm of a variable following a Beta distribution with parameters `w1` and `w2`.
"""
function e_ln_omega(w1::AbstractFloat,w2::AbstractFloat)
    return digamma(w1) - digamma(w1 + w2)
end

"""
    e_ln_minusomega(w1::AbstractFloat,w2::AbstractFloat)
Expected value of the logarithm of one minus a variable following a Beta distribution with parameters `w1` and `w2`.
"""
function e_ln_minusomega(w1::AbstractFloat,w2::AbstractFloat)
    return digamma(w2) - digamma(w1 + w2)
end

"""
    e_omega(w1::AbstractFloat,w2::AbstractFloat)
Expected value of a variable following a Beta distribution with parameters `w1` and `w2`.
"""
function e_omega(w1::AbstractFloat,w2::AbstractFloat)
    return w1 /(w1 + w2)
end

"""
    e_minusomega(w1::AbstractFloat,w2::AbstractFloat)
Expected value of one minus a variable following a Beta distribution with parameters `w1` and `w2`.
"""
function e_minusomega(w1::AbstractFloat,w2::AbstractFloat)
    return w2 /(w1 + w2)
end

"""
    e_chi(g1::AbstractFloat,g2::AbstractFloat)
Expected value of a variable following a Beta distribution with parameters `g1` and `g2`.
"""
function e_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return alpha_part/(alpha_part + beta_part)
end

"""
    e_minus_chi(g1::AbstractFloat,g2::AbstractFloat)
Expected value of one minus a variable following a Beta distribution with parameters `g1` and `g2`.
"""
function e_minus_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return beta_part/(alpha_part + beta_part)
end

"""
    e_ln_chi(g1::AbstractFloat,g2::AbstractFloat)
Expected value of the logarithm of a variable following a Beta distribution with parameters `g1` and `g2`.
"""
function e_ln_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return digamma(alpha_part) - digamma(alpha_part + beta_part)
end

"""
    e_ln_minus_chi(g1::AbstractFloat,g2::AbstractFloat)
Expected value of the logarithm of one minus a variable following a Beta distribution with parameters `g1` and `g2`.
"""
function e_ln_minus_chi(g1::AbstractFloat,g2::AbstractFloat)
    alpha_part = g1*g2
    beta_part = (1-g1)*g2
    return digamma(beta_part) - digamma(alpha_part + beta_part)
end

"""
    e_ln_lambda(u::AbstractFloat,v::AbstractFloat)
Expected value of the logarithm of a variable following a Gamma distribution with shape `u` and rate `v`.
"""
function e_ln_lambda(u::AbstractFloat,v::AbstractFloat)
    return log(v) -  digamma(u)
end

"""
    e_one_over_lambda(u::AbstractFloat,v::AbstractFloat)
Expected value of the inverse of a variable following a Gamma distribution with shape `u` and rate `v`.
"""
function e_one_over_lambda(u::AbstractFloat,v::AbstractFloat)
    return u/v
end

"""
    e_mu_sq(m_mu::AbstractFloat,s_sq_mu::AbstractFloat)
Expected value of the square of a variable following a Normal distribution with mean `m_mu` and variance `s_sq_mu`.
"""
function e_mu_sq(m_mu::AbstractFloat,s_sq_mu::AbstractFloat)
    return m_mu^2 + s_sq_mu
end

"""
    e_mu(m_mu::AbstractFloat)
Expected value of a variable following a Normal distribution with mean `m_mu`.
"""
function e_mu(m_mu::AbstractFloat)
    return m_mu
end

"""
    e_prior_sq_diff_mu(m_mu::AbstractFloat,s_sq_mu::AbstractFloat)
Expected value of the squared difference between a variable following a Normal distribution with mean `m_mu` and variance `s_sq_mu`, and its prior mean `m_mu`.
"""
function e_prior_sq_diff_mu(m_mu::AbstractFloat,s_sq_mu::AbstractFloat)
    return m_mu^2 -2*m_mu*m_mu + s_sq_mu + m_mu^2
end

"""
    e_nu_sq(m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
Expected value of the square of a variable following a Normal distribution with mean `m_nu` and variance `s_sq_nu`.
"""
function e_nu_sq(m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
    return m_nu^2 + s_sq_nu
end

"""
    e_nu_diff_sq(m_nu::AbstractFloat,s_sq_nu::AbstractFloat,nu0::AbstractFloat)
Expected value of the squared difference between a variable following a Normal distribution with mean `m_nu` and variance `s_sq_nu`, and its prior mean `nu0`.
"""
function e_nu_diff_sq(m_nu::AbstractFloat,s_sq_nu::AbstractFloat,nu0::AbstractFloat)
    return m_nu^2 + s_sq_nu -2*m_nu*nu0 + nu0^2
end

"""
    e_nu(m_nu::AbstractFloat)
Expected value of a variable following a Normal distribution with mean `m_nu`.
"""
function e_nu(m_nu::AbstractFloat)
    return m_nu
end

"""
    e_prior_sq_diff_nu(m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
Expected value of the squared difference between a variable following a Normal distribution with mean `m_nu` and variance `s_sq_nu`, and its prior mean `m_nu`.
"""
function e_prior_sq_diff_nu(m_nu::AbstractFloat,s_sq_nu::AbstractFloat)
    return m_nu^2 -2*m_nu*m_nu + s_sq_nu + m_nu^2
end

"""
    e_ln_eta(h1::AbstractFloat,h2::AbstractFloat)
Expected value of the logarithm of a variable following a Beta distribution with parameters `h1` and `h2`.
"""
function e_ln_eta(h1::AbstractFloat,h2::AbstractFloat)
    return digamma(h1) - digamma(h1 + h2)
end

"""
    e_ln_minus_eta(h1::AbstractFloat,h2::AbstractFloat)
Expected value of the logarithm of one minus a variable following a Beta distribution with parameters `h1` and `h2`.
"""
function e_ln_minus_eta(h1::AbstractFloat,h2::AbstractFloat)
    return digamma(h2) - digamma(h1 + h2)
end

"""
    e_eta(h1::AbstractFloat,h2::AbstractFloat)
Expected value of a variable following a Beta distribution with parameters `h1` and `h2`.
"""
function e_eta(h1::AbstractFloat,h2::AbstractFloat)
    return h1/(h1 + h2)
end

"""
    e_minus_eta(h1::AbstractFloat,h2::AbstractFloat)
Expected value of one minus a variable following a Beta distribution with parameters `h1` and `h2`.
"""
function e_minus_eta(h1::AbstractFloat,h2::AbstractFloat)
    return h2/(h1 + h2)
end

"""
    recursive_minus_e_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,g1::Vector{U},g2::Vector{U}) where {U <: AbstractFloat}# formerly recursive_minus_e_uk_cumprod
Recursively computes the cumulative product of (1 - e_chi) from index `k` down to 1, where `e_chi` is the expected value of a variable following a Beta distribution with parameters from vectors `g1` and `g2`. The recursion stops when `k` reaches 0, returning the cumulative product.
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
    log_of_recursive_minus_e_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,g1::Vector{U},g2::Vector{U})  where {U <: AbstractFloat} #  Logging for stability?
Recursively computes the logarithm of the cumulative product of (1 - e_chi) from index `k` down to 1, where `e_chi` is the expected value of a variable following a Beta distribution with parameters from vectors `g1` and `g2`. The recursion stops when `k` reaches 0, returning the cumulative sum of the logarithms.
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
    expectation_sbk(k::Int,K::Int,g1::Vector{U},g2::Vector{U};use_log=false)  where {U <: AbstractFloat}# formerly expectation_βk
Calculates the expected value of the stick-breaking weight `SB_k` for a given index `k` in a stick-breaking process with `K` components. The function uses parameters from vectors `g1` and `g2` to compute the expected values of the Beta-distributed variables involved in the stick-breaking process. If `use_log` is set to true, the function computes the logarithm of the cumulative product for numerical stability.
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

#################################################
###### CUSTOM TYPE EXPECTATION FUNCTIONS ########
#################################################

"""
    E_ln_sigma_sq(j,clusterfeature::ClusterFeature)
Calculates the expected value of the logarithm of the sigma squared parameter for a given gene index `j` in a cluster feature.
""" 
function E_ln_sigma_sq(j,clusterfeature::ClusterFeature)
    return e_ln_sigma_sq(clusterfeature.a[j],clusterfeature.b[j])
end

"""
    E_one_over_sigma_sq(j,clusterfeature::ClusterFeature)   
Calculates the expected value of the reciprocal of the sigma squared parameter for a given gene index `j` in a cluster feature.
""" 
function E_one_over_sigma_sq(j,clusterfeature::ClusterFeature)
    return e_one_over_sigma_sq(clusterfeature.a[j],clusterfeature.b[j])
end


""" 
    E_ll_sq_diff_mu(j,cellfeature::CellFeature,clusterfeature::ClusterFeature)
Calculates the expected value of the squared difference between the observed expression and the mean expression for a given gene index `j` in a cell feature and a cluster feature.
"""     
function E_ll_sq_diff_mu(j,cellfeature::CellFeature,clusterfeature::ClusterFeature)
    return e_ll_sq_diff_mu(cellfeature.xsq[j],cellfeature.x[j],clusterfeature.m_mu[j],clusterfeature.y[j],clusterfeature.s_sq_mu[j],clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j])
end


"""
    E_ln_pi(k,conditionfeature::ConditionFeature)
Calculates the expected value of the logarithm of the pi parameter for a given condition index `k` in a condition feature.
"""     
function E_ln_pi(k,conditionfeature::ConditionFeature)
    return e_ln_pi(conditionfeature.d[k],conditionfeature.d_sum[1])
end


"""    
    E_ln_omega(conditionfeature::ConditionFeature)
Calculates the expected value of the logarithm of the omega parameter for a given condition feature.
""" 
function E_ln_omega(conditionfeature::ConditionFeature)
    return e_ln_omega(conditionfeature.w1[1],conditionfeature.w2[1])
end

"""    
    E_ln_minusomega(conditionfeature::ConditionFeature)
Calculates the expected value of the logarithm of the minus omega parameter for a given condition feature.
""" 
function E_ln_minusomega(conditionfeature::ConditionFeature)
    return e_ln_minusomega(conditionfeature.w1[1],conditionfeature.w2[1])
end

"""
    E_omega(conditionfeature::ConditionFeature) 
Calculates the expected value of the omega parameter for a given condition feature.
"""    
function E_omega(conditionfeature::ConditionFeature)
    return e_omega(conditionfeature.w1[1],conditionfeature.w2[1])
end

"""
    E_minusomega(conditionfeature::ConditionFeature)
Calculates the expected value of the minus omega parameter for a given condition feature.
"""    
function E_minusomega(conditionfeature::ConditionFeature)
    return e_minusomega(conditionfeature.w1[1],conditionfeature.w2[1])
end

"""    
    E_chi(clusterfeature::ClusterFeature)
Calculates the expected value of the chi parameter for a given cluster feature.
"""    
function E_chi(clusterfeature::ClusterFeature)
    return e_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

"""    
    E_minus_chi(clusterfeature::ClusterFeature) 
Calculates the expected value of the minus chi parameter for a given cluster feature.
"""    
function E_minus_chi(clusterfeature::ClusterFeature)
    return e_minus_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

"""    
    E_ln_chi(clusterfeature::ClusterFeature)    
Calculates the expected value of the logarithm of the chi parameter for a given cluster feature.
"""    
function E_ln_chi(clusterfeature::ClusterFeature)
    return e_ln_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

"""    
    E_ln_minus_chi(clusterfeature::ClusterFeature)
Calculates the expected value of the logarithm of the minus chi parameter for a given cluster feature.
"""     
function E_ln_minus_chi(clusterfeature::ClusterFeature)
    return e_ln_minus_chi(clusterfeature.g1[1],clusterfeature.g2[1])
end

"""    
    E_ln_lambda(j,clusterfeature::ClusterFeature)
Calculates the expected value of the logarithm of the lambda parameter for a given gene index `j` in a cluster feature.
"""    
function E_ln_lambda(j,clusterfeature::ClusterFeature)
    return e_ln_lambda(clusterfeature.u[j],clusterfeature.v[j])
end

"""    
    E_one_over_lambda(j,clusterfeature::ClusterFeature)
Calculates the expected value of the one over lambda parameter for a given gene index `j` in a cluster feature.
"""    
function E_one_over_lambda(j,clusterfeature::ClusterFeature)
    return e_one_over_lambda(clusterfeature.u[j],clusterfeature.v[j])
end

"""    
    E_mu_sq(j,clusterfeature::ClusterFeature)   
Calculates the expected value of the squared mu parameter for a given gene index `j` in a cluster feature.
"""    
function E_mu_sq(j,clusterfeature::ClusterFeature)
    return e_mu_sq(clusterfeature.m_mu[j],clusterfeature.s_sq_mu[j])
end

"""    
    E_mu(j,clusterfeature::ClusterFeature)  
Calculates the expected value of the mu parameter for a given gene index `j` in a cluster feature.
"""    
function E_mu(j,clusterfeature::ClusterFeature)
    return e_mu(clusterfeature.m_mu[j])
end

"""    
    E_prior_sq_diff_mu(j,clusterfeature::ClusterFeature)    
Calculates the expected value of the prior squared difference of the mu parameter for a given gene index `j` in a cluster feature.
"""
function E_prior_sq_diff_mu(j,clusterfeature::ClusterFeature)
    return e_prior_sq_diff_mu(clusterfeature.m_mu[j],clusterfeature.s_sq_mu[j])
end

"""    
    E_nu_sq(j,clusterfeature::ClusterFeature)
Calculates the expected value of the squared nu parameter for a given gene index `j` in a cluster feature.
"""
function E_nu_sq(j,clusterfeature::ClusterFeature)
    return e_nu_sq(clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j])
end

"""    
    E_nu(j,clusterfeature::ClusterFeature)
Calculates the expected value of the nu parameter for a given gene index `j` in a cluster feature.
"""
function E_nu(j,clusterfeature::ClusterFeature)
    return e_nu(clusterfeature.m_nu[j])
end

"""    
    E_prior_sq_diff_nu(j,clusterfeature::ClusterFeature)
Calculates the expected value of the prior squared difference of the nu parameter for a given gene index `j` in a cluster feature.
"""
function E_prior_sq_diff_nu(j,clusterfeature::ClusterFeature)
    return e_prior_sq_diff_nu(clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j])
end

"""    
    E_nu_diff_sq(j,clusterfeature::ClusterFeature,modelparams::ModelParameterFeature)   
Calculates the expected value of the squared difference of the nu parameter for a given gene index `j` in a cluster feature.
"""
function E_nu_diff_sq(j,clusterfeature::ClusterFeature,modelparams::ModelParameterFeature)
    return e_nu_diff_sq(clusterfeature.m_nu[j],clusterfeature.s_sq_nu[j],modelparams.nu0[j])
end

"""    
    E_ln_eta(j,clusterfeature::ClusterFeature)  
Calculates the expected value of the logarithm of the eta parameter for a given gene index `j` in a cluster feature.
"""
function E_ln_eta(j,clusterfeature::ClusterFeature)
    return e_ln_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

"""    
    E_ln_minus_eta(j,clusterfeature::ClusterFeature)    
Calculates the expected value of the logarithm of the complement of the eta parameter for a given gene index `j` in a cluster feature.
"""
function E_ln_minus_eta(j,clusterfeature::ClusterFeature)
    return e_ln_minus_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

"""    
    E_eta(j,clusterfeature::ClusterFeature)
Calculates the expected value of the eta parameter for a given gene index `j` in a cluster feature.
"""
function E_eta(j,clusterfeature::ClusterFeature)
    return e_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

"""    
    E_minus_eta(j,clusterfeature::ClusterFeature)   
Calculates the expected value of the complement of the eta parameter for a given gene index `j` in a cluster feature.
"""
function E_minus_eta(j,clusterfeature::ClusterFeature)
    return e_minus_eta(clusterfeature.h1[j],clusterfeature.h2[j])
end

#copy_cells = deepcopy(cells); copy_precomputed_genefeatures_cells = deepcopy(cells); for n in 1:dataparams.N copy_cells[n].cache .= 0.0; copy_precomputed_genefeatures_cells[n].cache .= 0.0; adjust_E_ln_pi!(copy_cells[n],conditions,dataparams);adjust_E_ln_pi!(copy_precomputed_genefeatures_cells[n],conditions,dataparams); for k in 1:modelparams.K E_log_normal_l_j!(copy_cells[n],clusters[k], dataparams); E_log_normal_l_j!(preupdated_genefeatures,copy_precomputed_genefeatures_cells[n],clusters[k], dataparams); end; end; copy_cells_z_argmax = [argmax(copy_cells[n].cache) for n in 1:dataparams.N]; copy_precomputed_genefeatures_cells_z_argmax = [argmax(copy_precomputed_genefeatures_cells[n].cache) for n in 1:dataparams.N];println(Clustering.randindex(cell_cluster_labels,cell_cluster_labels)[1]); println(countmap(copy_cells_z_argmax)); println(Clustering.randindex(cell_cluster_labels,copy_cells_z_argmax)[1]); println(countmap(copy_precomputed_genefeatures_cells_z_argmax)); println(Clustering.randindex(cell_cluster_labels,copy_precomputed_genefeatures_cells_z_argmax)[1]); println(all([el in collect(1:dataparams.N)[cell_cluster_labels .!= copy_cells_z_argmax]  for el in collect(1:dataparams.N)[cell_cluster_labels .!= copy_precomputed_genefeatures_cells_z_argmax]]));

"""
    recursive_minus_E_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,clusters::Vector{ClusterFeature{U,W}}) where {U <: AbstractFloat, W <: Int64} # formerly recursive_minus_e_uk_cumprod
Recursively computes the cumulative product of (1 - E_chi) from index `k` down to 1, where `E_chi` is the expected value of a variable following a Beta distribution with parameters from a vector of cluster features. The recursion stops when `k` reaches 0, returning the cumulative product.
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
    log_of_recursive_minus_E_chi_cumprod(k::Int,cummulative_prod::AbstractFloat,clusters::Vector{ClusterFeature{U,W}}) where {U <: AbstractFloat, W <: Int64} #  Logging for stability?
Recursively computes the cumulative sum of the logarithm of (1 - E_chi) from index `k` down to 1, where `E_chi` is the expected value of a variable following a Beta distribution with parameters from a vector of cluster features. The recursion stops when `k` reaches 0, returning the cumulative sum.
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
    expectation_SBk(k::Int,clusters::Vector{ClusterFeature{U,W}},modelparams::ModelParameterFeature;use_log=false) where {U <: AbstractFloat, W <: Int64} # formerly expectation_βk
Calculates the expected value of the stick-breaking weight `SB_k` for a given index `k` in a stick-breaking process with `K` components, using a vector of cluster features to obtain the necessary parameters. If `use_log` is set to true, the function computes the logarithm of the cumulative product for numerical stability.
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


# sum(T[1:i-1])+tt
"""
    recursive_cumsum_E_ln_minusomega(i::Int,t::Int,tt::Int,T::Vector{Int},conditions::Vector{ConditionFeature{U,W}}) where {U <: AbstractFloat, W <: Int64}
Recursively computes the cumulative sum of the expected value of the logarithm of one minus omega for a given condition feature, iterating from `tt` to `t` for the `i`-th condition. The recursion stops when `tt` exceeds `t`, returning 0.0.
"""
function recursive_cumsum_E_ln_minusomega(i::Int,t::Int,tt::Int,T::Vector{Int},conditions::Vector{ConditionFeature{U,W}}) where {U <: AbstractFloat, W <: Int64}
    if t < tt
        return 0.0
    else
        if tt == 1 && t == 1
            return 0.0
        else
            return E_ln_minusomega(conditions[sum(T[1:i-1])+tt]) + recursive_cumsum_E_ln_minusomega(i,t,tt+1,T,conditions)
        end
    end
end

############################################
############################################
############################################

"""
    βk_expected_value(rho_hat_vec, omega_hat_vec)
Calculates the expected value of the stick-breaking weight `β_k` for a given index `k` in a stick-breaking process with `K` components, using a vector of cluster features to obtain the necessary parameters. If `use_log` is set to true, the function computes the logarithm of the cumulative product for numerical stability.
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
    uk_expected_value(rho_hat_vec, omega_hat_vec)
Calculates the expected value of the stick-breaking variable `u_k` for each component in a stick-breaking process, given vectors of parameters `rho_hat_vec` and `omega_hat_vec`.
"""
function uk_expected_value(rho_hat_vec, omega_hat_vec)
    minus_rho_hat_vec = 1 .- rho_hat_vec
    rho_omega_hat_vec = rho_hat_vec .* omega_hat_vec
    minus_rho_omega_hat_vec = minus_rho_hat_vec .*omega_hat_vec
    e_uk_vec = rho_omega_hat_vec ./ (rho_omega_hat_vec .+ minus_rho_omega_hat_vec)
    return e_uk_vec
end

"""
    logUk_expected_value(rho_hat,omega_hat)
Calculates the expected value of the logarithm of `u_k` for each component in a stick-breaking process, given vectors of parameters `rho_hat` and `omega_hat`.
"""
function logUk_expected_value(rho_hat,omega_hat)
    return digamma.(rho_hat .* omega_hat) .- digamma.(omega_hat)
end

"""
    log1minusUk_expected_value(rho_hat,omega_hat)
Calculates the expected value of the logarithm of `1 - u_k` for each component in a stick-breaking process, given vectors of parameters `rho_hat` and `omega_hat`.
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
    adjust_e_log_π_tk3!(t,conditionparams)
Adjusts the cached expected logarithm of the mixing proportions `π_tk` for a specific time point `t` across all conditions, based on the current estimates of the Dirichlet parameters. This function updates the `e_log_π_t_cache` for the condition at time `t` by summing contributions from all previous time points weighted by their respective transition probabilities.
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
    log_π_expected_value_fast3!(conditionparams,dataparams,modelparams)
Calculates the expected logarithm of the mixing proportions `π_tk` for all time points and conditions, storing the results in the `e_log_π_t_cache` of each condition parameters. This function iterates over all time points and components, computing the expected logarithm based on the current estimates of the Dirichlet parameters.
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
    expectation_βk(k,clusters,modelparams)
Calculates the expected value of the stick-breaking variable `β_k` for a specific component `k` in a stick-breaking process, given the current estimates of the parameters.
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
    recursive_minus_e_uk_cumprod(k::Int,cummulative_prod::AbstractFloat,clusters::Vector{ClusterFeature{U,W}}) where {U <: AbstractFloat, W <: Int64} # formerly recursive_minus_e_uk_cumprod
Recursively computes the cumulative product of (1 - e_uk) from index `k` down to 1, where `e_uk` is the expected value of a variable following a Beta distribution with parameters from a vector of cluster features. The recursion stops when `k` reaches 0, returning the cumulative product.
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
    expectation_log_normal_l_j!(cluster,cell, dataparams)
Calculates the expected logarithm of the likelihood of the observed data for a given cell and cluster, updating the cell's cache with the computed values for each gene. This function uses the current estimates of the cluster parameters to compute the expected log-likelihood.
"""
function expectation_log_normal_l_j!(cluster,cell, dataparams)
    G = dataparams.G
    for j in 1:G
        cell.cache[j]  =  -0.5 * (log(cluster.σ_sq_k_hat[j]) + 1/cluster.σ_sq_k_hat[j] * ( cell.xsq[j] - 2 * cell.x[j] * cluster.κk_hat[j] + cluster.yjk_hat[j] * (cluster.mk_hat[j] ^2 + cluster.v_sq_k_hat[j])  ))
        # cell.cache[j]  =  -0.5 * (log(cluster.σ_sq_k_hat[j]) + 1/cluster.σ_sq_k_hat[j] * ( cell.xsq[j] - 2 * cell.x[j] * (cluster.yjk_hat[j] * cluster.mk_hat[j]) + cluster.yjk_hat[j] * (cluster.mk_hat[j] ^2 + cluster.v_sq_k_hat[j])  ))
    end
    return cluster
end


function expectation_log_π_tk(θ_hat_tk,θ_hat_t_sum)
    e_log_π_tk =  digamma(θ_hat_tk) - digamma(θ_hat_t_sum)
    return e_log_π_tk 
end



function expectation_log_tilde_wtt(awt_hat,bwt_hat)
    return  digamma(awt_hat) - digamma(awt_hat + bwt_hat)
end


function expectation_log_minus_tilde_wtt(awt_hat,bwt_hat)
    return  digamma(bwt_hat) - digamma(awt_hat + bwt_hat)
end


function expectation_uk(rho_hat, omega_hat)
    # minus_rho_hat = 1 - rho_hat
    # rho_omega_hat = rho_hat * omega_hat
    # minus_rho_omega_hat = (1 - rho_hat) * omega_hat
    e_uk = (rho_hat * omega_hat) / ((rho_hat * omega_hat)+ ((1 - rho_hat) * omega_hat))
    return e_uk
end


function expectation_αt(a_α,b_α)
    return a_α / b_α
end


function expectation_log_αt(a_α,b_α)
    return digamma(a_α) - log(b_α)
end


function expectation_logUk(rho_hat,omega_hat)
    return digamma(rho_hat * omega_hat) - digamma(omega_hat)
end


function expectation_log1minusUk(rho_hat,omega_hat)
    return digamma((1.0 - rho_hat) * omega_hat) - digamma(omega_hat)
end



