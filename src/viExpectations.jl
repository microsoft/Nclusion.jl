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



