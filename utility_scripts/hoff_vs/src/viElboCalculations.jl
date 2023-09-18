function calc_Hz(rtik)
    z_entropy = 0.0
    T = length(rtik)
    C_t = [length(el) for el in rtik]
    K = length(rtik[1][1])
    for t in 1:T
        for i in 1:C_t[t]
            z_entropy += entropy(rtik[t][i][1:K])
        end
    end
    return  -z_entropy
end
function calc_Hv(v_tikj)
    v_entropy = 0.0
    T = length(v_tikj)
    C_t = [length(el) for el in v_tikj]
    K = length(v_tikj[1][1])
    G = length(v_tikj[1][1][1])
    for t in 1:T
        for i in 1:C_t[t]
            for k in 1:K
                for j in 1:G
                    v_entropy += entropy(v_tikj[t][i][k][j])
                end
            end
        end
    end
    return  -v_entropy
end
function calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
    K = length(b_ηkj_hat)
    G = length(b_ηkj_hat[1])
    T = length(v_tikj)
    C_t = length.(v_tikj)
    a_eta = [a_η for j in 1:G]
    b_eta = [b_η for j in 1:G]
    v_true =  [[[[v_tikj[t][i][k][j][1] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    v_kj = sum(sum(v_true))
    v_false =  [[[[v_tikj[t][i][k][j][2] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    minus_v_kj = sum(sum(v_false))
    cBeta_sum =  sum([-c_Beta(a_eta,b_eta) .+ c_Beta(a_ηkj_hat[k],b_ηkj_hat[k])  for k in 1:K])
    logη_sum =sum([(a_eta .- a_ηkj_hat[k] + v_kj[k]) .*e_log_ηkj[k] for k in 1:K]) 
    logminusη_sum = sum([(b_eta .- b_ηkj_hat[k] + minus_v_kj[k]) .*e_log_minus_ηkj[k] for k in 1:K]) 

    imp_lb  = sum(cBeta_sum + logη_sum +logminusη_sum)
    return imp_lb
end
function calc_ImportanceElbo(v_tikj,ηkj_prior)
    K = length(v_tikj[1][1])
    G = length(v_tikj[1][1][1])
    T = length(v_tikj)
    C_t = length.(v_tikj)

    imp_lb =  [[[[v_tikj[t][i][k][j][1]*log(ηkj_prior[k][j]) + v_tikj[t][i][k][j][2]*log(1 - ηkj_prior[k][j])  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    imp_lb  = sum(sum(sum(sum(sum.(imp_lb)))))
    return imp_lb
end
function calc_DataElbo7(x,rtik,v_tikj,mk_hat_vec,μ0_vec,μ0_err_vec,m_err_hat_vec,λ0k_hat_vec,λ0_vec,λ0_err_vec,λ0_err_hat_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
    data_lb_sum = 0.0
    K = length(rtik[1][1][1])
    N_signal,N_error = update_N(rtik,v_tikj);
    # Nj_error = update_errorNj(N_error);
    # Nkj_signal = update_signalNkj(N_signal);       
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
    x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
    x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
    x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
    x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
    sum_ti_rtik_v_tikj = sum(sum.(N_signal))
    sum_ti_rtik_minusv_tikj = sum(sum.(N_error))
    Nk = update_Nk(rtik)
    E_q_τ_μ_sq(l0,m0,a0,b0) = 1 ./ l0 .+ m0 .^2 .* a0 ./ b0
    E_q_τ_μ(m0,a0,b0) = m0 .* a0 ./ b0
    E_q_τ(a0,b0) =  a0 ./ b0
    for k in 1:K
        log_λ_sum = 0.5*(log.(λ0_vec) .+ log.(λ0_err_vec) .- log.(λ0k_hat_vec[k]) .- log.(λ0_err_hat_vec))
        log_c_Gamma_sum = c_Ga.(a0_vec, b0_vec) .+  c_Ga.(a0_err_vec, b0_err_vec) .- c_Ga.(a0k_hat_vec[k], b0k_hat_vec[k]) .-  c_Ga.(a0_err_hat_vec, b0_err_hat_vec) 
        e_log_τ_kj = log_τ_kj_expected_value(a0k_hat_vec[k], b0k_hat_vec[k])
        e_log_τ_kj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
        a0_sum = (a0_vec - a0k_hat_vec[k] + 0.5* sum_ti_rtik_v_tikj[k]) .* e_log_τ_kj .+ (a0_err_vec .- a0_err_hat_vec .+ sum_ti_rtik_minusv_tikj[k]) .*e_log_τ_kj_err
        log2π_sum = -0.5 .* log(2π) .* Nk[k]
        e_q_τ_μ_sq_sum = 0.5 .*(λ0k_hat_vec[k] .- λ0_vec - sum_ti_rtik_v_tikj[k]) .* E_q_τ_μ_sq(λ0k_hat_vec[k],mk_hat_vec[k],a0k_hat_vec[k], b0k_hat_vec[k])
        e_q_τ_μ_sq_err_sum = 0.5 .*(λ0_err_hat_vec .- λ0_err_vec - sum_ti_rtik_minusv_tikj[k]) .* E_q_τ_μ_sq(λ0_err_hat_vec,m_err_hat_vec,a0_err_hat_vec, b0_err_hat_vec)
        e_q_τ_μ_sum = (λ0_vec .* μ0_vec - λ0k_hat_vec[k] .* mk_hat_vec[k] +  x_hat_k[k]) .*  E_q_τ_μ(mk_hat_vec[k],a0k_hat_vec[k], b0k_hat_vec[k])
        e_q_τ_μ_err_sum = (λ0_err_vec .* μ0_err_vec - λ0_err_hat_vec .* m_err_hat_vec +  x_hat_err) .*  E_q_τ_μ(m_err_hat_vec,a0_err_hat_vec, b0_err_hat_vec)
        e_q_τ_sum = ((b0k_hat_vec[k] .- b0_vec) .+ 0.5 .*(λ0k_hat_vec[k] .* mk_hat_vec[k] .^2 .-  λ0_vec .* μ0_vec .^2 .- x_hat_sq_k[k])) .* E_q_τ(a0k_hat_vec[k], b0k_hat_vec[k])
        e_q_τ_err_sum = ((b0_err_hat_vec .-  b0_err_vec) .+ 0.5 .*(λ0_err_hat_vec .* m_err_hat_vec .^2 -  λ0_err_vec .* μ0_err_vec .^2 - x_hat_sq_err ) ) .* E_q_τ(a0_err_hat_vec, b0_err_hat_vec)
        val = sum(log_λ_sum .+ log_c_Gamma_sum .+ a0_sum .+ log2π_sum .+ e_q_τ_μ_sq_sum .+ e_q_τ_μ_sq_err_sum .+ e_q_τ_μ_sum .+ e_q_τ_μ_err_sum .+ e_q_τ_sum .+ e_q_τ_err_sum)
        data_lb_sum += val
    end
    return data_lb_sum
end

function calc_DataElbo12(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
    data_lb_sum = 0.0
    N_signal,N_error = update_N(rtik,v_tikj);
    K = length(rtik[1][1][1])
    # Nj_error = update_errorNj(N_error);
    # Nkj_signal = update_signalNkj(N_signal);       
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
    x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
    x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
    x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
    x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
    sum_ti_rtik_v_tikj = sum(sum.(N_signal))
    sum_ti_rtik_minusv_tikj = sum(sum.(N_error))
    Nk = update_Nk(rtik)
    E_q_τ_μ_sq(l0,m0,a0,b0) = 1 ./ l0 .+ m0 .^2 .* a0 ./ b0
    E_q_τ_μ(m0,a0,b0) = m0 .* a0 ./ b0
    E_q_τ(a0,b0) =  a0 ./ b0
    for k in 1:K
        log_λ_sum = 0.5*(log.(λ0_vec)  .- log.(λ0k_hat_vec[k]))
        log_c_Gamma_sum = c_Ga.(a0_vec, b0_vec) .+  c_Ga.(a0_err_vec, b0_err_vec) .- c_Ga.(a0k_hat_vec[k], b0k_hat_vec[k]) .-  c_Ga.(a0_err_hat_vec, b0_err_hat_vec) 
        e_log_τ_kj = log_τ_kj_expected_value(a0k_hat_vec[k], b0k_hat_vec[k])
        e_log_τ_kj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
        a0_sum = (a0_vec - a0k_hat_vec[k] + 0.5* sum_ti_rtik_v_tikj[k]) .* e_log_τ_kj .+ (a0_err_vec .- a0_err_hat_vec .+ sum_ti_rtik_minusv_tikj[k]) .*e_log_τ_kj_err
        log2π_sum = -0.5 .* log(2π) .* Nk[k]
        e_q_τ_μ_sq_sum = 0.5 .*(λ0k_hat_vec[k] .- λ0_vec - sum_ti_rtik_v_tikj[k]) .* E_q_τ_μ_sq(λ0k_hat_vec[k],mk_hat_vec[k],a0k_hat_vec[k], b0k_hat_vec[k])
        # e_q_τ_μ_sq_err_sum = 0.5 .*(λ0_err_hat_vec .- λ0_err_vec - sum_ti_rtik_minusv_tikj[k]) .* E_q_τ_μ_sq(λ0_err_hat_vec,m_err_hat_vec,a0_err_hat_vec, b0_err_hat_vec)
        e_q_τ_μ_sum = (λ0_vec .* μ0_vec - λ0k_hat_vec[k] .* mk_hat_vec[k] +  x_hat_k[k]) .*  E_q_τ_μ(mk_hat_vec[k],a0k_hat_vec[k], b0k_hat_vec[k])
        # e_q_τ_μ_err_sum = (λ0_err_vec .* μ0_err_vec - λ0_err_hat_vec .* m_err_hat_vec +  x_hat_err) .*  E_q_τ_μ(m_err_hat_vec,a0_err_hat_vec, b0_err_hat_vec)
        e_q_τ_sum = ((b0k_hat_vec[k] .- b0_vec) .+ 0.5 .*(λ0k_hat_vec[k] .* mk_hat_vec[k] .^2 .-  λ0_vec .* μ0_vec .^2 .- x_hat_sq_k[k])) .* E_q_τ(a0k_hat_vec[k], b0k_hat_vec[k])
        e_q_τ_err_sum = ((b0_err_hat_vec .-  b0_err_vec) .- 0.5 .*(x_hat_sq_err ) ) .* E_q_τ(a0_err_hat_vec, b0_err_hat_vec)
        val = sum(log_λ_sum .+ log_c_Gamma_sum .+ a0_sum .+ log2π_sum .+ e_q_τ_μ_sq_sum .+ e_q_τ_μ_sum .+ e_q_τ_sum .+ e_q_τ_err_sum)
        data_lb_sum += val
    end
    return data_lb_sum
end


function calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(a0k_hat_vec)
    G = length(a0k_hat_vec[1])
    # halfNklogπ = 1/2 .* Nk .*log(2π) # K by 1
    halfNklogπ = -1/2 .* Nk .*log(2π) # K by 1
    halfNkoverλ0k_hat = 1/2 .* [ Nk[k] ./λ0k_hat_vec[k] for k in 1:K] # K by G
    a0b0_c_Ga =  [c_Ga.(a0_vec, b0_vec) for k in 1:K] # K by G [zeros(G) for k in 1:K]#
    # a0b0_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#
    a0kb0k_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#

    
    
    halfLogλ0 = 1/2 .* [ log.(λ0_vec) for k in 1:K] # K by G
    halfLogλ0k = 1/2 .* [ log.(λ0k_hat_vec[k]) for k in 1:K] # K by G
    halfλ0μ0_sq = 1/2 .* λ0_vec.* μ0_vec .^ 2 #G By 1
    a0kOverb0k_hat = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K] # K by G
    bs_e_τ = [(b0k_hat_vec[k]  .- b0_vec .- halfλ0μ0_sq) .*a0kOverb0k_hat[k]  for k in 1:K] #K by G
    # digamma.(a_0kj_hat) .- log.(b_0kj_hat) 
    # e_log_τkj = log_τ_expected_value.(a0k_hat_vec, b0k_hat_vec) #K by G
    as_e_log_τkj =[(1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])) for k in 1:K] #K by G
    λ0μ0 = λ0_vec .* μ0_vec #G By 1
    λ0μ0mk_a0kOverb0k_hat = [ λ0μ0 .* mk_hat_vec[k] .* a0kOverb0k_hat[k] for k in 1:K] #K by G
    halfλ0mk_sq_a0kOverb0k_hat = [1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* a0kOverb0k_hat[k]  for k in 1:K] #K by G
    one_half_const = 1/2
    halfλ0overλ0k_hat = [1/2 .* λ0_vec ./ λ0k_hat_vec[k] for k in 1:K] #K by G
    weighted_ss_kjti = [[[[rtik[t][i][k] .*(x[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
    weighted_ss_kjt = [[[sum(weighted_ss_kjti[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    weighted_ss_kj = [[sum(weighted_ss_kjt[k][j]) for j in 1:G] for k in 1:K] #K by G
    halfa0kOverb0k_hat_weighted_ss_kj =  [ 1/2 .* a0kOverb0k_hat[k] .* weighted_ss_kj[k]  for k in 1:K] #K by G

    # data_elbo_kj = [ halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .+ a0b0_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .- one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    data_elbo_kj = [halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .- a0kb0k_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .+ one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    data_elbo_j  = sum(data_elbo_kj)
    data_elbo  = sum(data_elbo_j)
    #     
    return data_elbo
end

function calc_SurragateLowerBound(rho_hat,omega_hat,T,γ,α0,Tk)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α0 .* e_βk .* Tk
    #take Negative to find max  
    lb_lg = sum(lb_lg_k)
    return lb_lg
end
function calc_SurragateLowerBound_unconstrained(c,d,T,γ,α0,Tk)
    rho_hat = sigmoid.(c)
    omega_hat = exp.(d)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α0 .* e_βk .* Tk[1:K]
    #take Negative to find max  
    lb_lg = sum(lb_lg_k)
    return lb_lg
end
function calc_SurragateLowerBound_unconstrained(c,d,T,e_γ,Tαk)
    rho_hat = sigmoid.(c)
    omega_hat = exp.(d)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ e_γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  e_βk .* Tαk[1:K]
    #take Negative to find max  
    lb_lg = sum(lb_lg_k)
    return lb_lg
end
function c_Ga(a0, b0)
    # try
    #     a0 .* log.(b0) .- loggamma.(a0) 
    # catch
    #     println("a0: $a0")
    #     println("b0: $b0")
    # end
    a0 .* log.(b0) .- loggamma.(a0)
end
function c_Beta(a0, b0)
    - logbeta.(a0,b0) 
end

function calc_Hs(c_ttprime)
    s_entropy = 0.0
    T = length(c_ttprime)
    for t in 1:T
        s_entropy += entropy(c_ttprime[t])
    end
    return -s_entropy
end

function calc_wAllocationsLowerBound(c_ttprime, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    T = length(c_ttprime)
    wAlloc_elbo = 0.0
    c_Beta_p = c_Beta.(adot_w,bdot_w)
    c_Beta_q = -c_Beta.(awt_hat_vec,bwt_hat_vec)
    c_Beta_pq = c_Beta_p .+ c_Beta_q
    e_log_tilde_wt = log_tilde_wt_expected_value.(awt_hat_vec,bwt_hat_vec)
    e_log_minus_tilde_wt = log_minus_tilde_wt_expected_value.(awt_hat_vec,bwt_hat_vec)
    adot_w_awt_hat_vec = adot_w .- awt_hat_vec
    bdot_w_bwt_hat_vec = bdot_w .- bwt_hat_vec
    for t in 2:T
        a_cttprime = 0.0
        for t_prime_a  in t:T
            a_cttprime += c_ttprime[t_prime_a][t]
        end 

        b_cttprime = 0.0
        for t_prime_b  in t:T
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                b_cttprime += c_ttprime[t_prime_b][l]
                # sum_string *= c_string
            end
        end
        wAlloc_elbo_t = c_Beta_pq[t-1] .+ (adot_w_awt_hat_vec[t-1] .+ a_cttprime) .* e_log_tilde_wt[t-1] .+ (bdot_w_bwt_hat_vec[t-1] .+ b_cttprime) .* e_log_minus_tilde_wt[t-1]
        wAlloc_elbo += wAlloc_elbo_t
    end
    return wAlloc_elbo
end

function calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    c_Gamma_γ_p = c_Ga(a_γ,b_γ)
    c_Gamma_γ_q = -c_Ga(a_γ_hat,b_γ_hat)
    e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
    e_log_γ  = log_γ_expected_value(a_γ_hat,b_γ_hat)
    γ_elbo = c_Gamma_γ_p + c_Gamma_γ_q+(a_γ- a_γ_hat)*e_log_γ - (b_γ-b_γ_hat) * e_γ
    return γ_elbo
end


function calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    T = length(a_αt_hat_vec)
    α_elbo = 0.0
    a_αt_vec = [a_α for t in 1:T]
    b_αt_vec = [b_α for t in 1:T]
    c_Gamma_α_p = c_Ga.(a_αt_vec,b_αt_vec)
    c_Gamma_α_q = -c_Ga(a_αt_hat_vec,b_αt_hat_vec)
    e_αt = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    e_log_αt  = log_αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    for t in 1:T
        αt_elbo = c_Gamma_α_p[t] .+ c_Gamma_α_q[t] .+ (a_αt_vec[t] .- a_αt_hat_vec[t]) .*e_log_αt[t] .- (b_αt_vec[t] .-b_αt_hat_vec[t]) .* e_αt[t]
        α_elbo += αt_elbo
    end
    return α_elbo
end



#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################
function calculate_elbo(Tk,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
    elbo_val = 0.0
    elbo_val += calc_DataElbo(cellpop,clusters,geneparams,dataparams,modelparams)
    elbo_val += calc_Hz_fast3(cellpop,clusters,dataparams)
    elbo_val += calc_HyjkSurragateLowerBound_unconstrained(Tk,clusters,dataparams,modelparams)
    elbo_val += calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    # elbo_val += calc_alphaElbo_fast3(conditionparams,dataparams,modelparams)
    # elbo_val += calc_HsGammaAlphaElbo_fast3(a_γ,b_γ,a_γ_hat,b_γ_hat,conditionparams,dataparams,modelparams)
    elbo_val +=  calc_HsElbo(conditionparams, dataparams, modelparams)
    return elbo_val
end

function calculate_elbo_perK(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
    K = modelparams.K
    for k in 1:K
        elbolog.per_k_elbo[k,iter] = 0.0
    end
    elbo_val = 0.0
    dataelbo,elbolog = calc_DataElbo_perK(cellpop,clusters,geneparams,elbolog,dataparams,modelparams,iter)
    zentropy = calc_Hz_fast3(cellpop,clusters,dataparams)
    lg_elbo,elbolog =calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    w_elbo = calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    sentropy = calc_HsElbo(conditionparams, dataparams, modelparams)
    elbo_val +=  dataelbo + zentropy + lg_elbo + w_elbo + sentropy
    for k in 1:K
        elbolog.per_k_elbo[k,iter] += zentropy + w_elbo + sentropy
    end
    return elbo_val,elbolog
end
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

# data_elbo =  calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec);
function calc_DataElbo(cellpop,clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    one_half_const = 1/2
    λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
    data_elbo = 0.0
    @fastmath @inbounds @simd for k in 1:K
        clusters[k].cache .= 0.0
        # @fastmath @inbounds @simd for i in 1:N
        #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
        # end
        # try
        #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
        # catch
        #     println("K =  $k")
        #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
        #     println("λ_sq_vec =  $(λ_sq_vec)")
        #     println("ηk =  $(modelparams.ηk)")
        #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
        #     println("var_muk =  $(clusters[k].var_muk )")
        #     println("mk_hat =  $(clusters[k].mk_hat)")

        # end

        
    #    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
       data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1] ) .+ (1 .- clusters[k].yjk_hat) .* log((1-(modelparams.ηk[1]) )) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
    end
    #     
    return data_elbo
end

function calc_DataElbo_perK(cellpop,clusters,geneparams,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    one_half_const = 1/2
    λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
    data_elbo = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        clusters[k].cache .= 0.0
        yjk_entropy_perK = 0.0
        # @fastmath @inbounds @simd for i in 1:N
        #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
        # end
        # try
        #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
        # catch
        #     println("K =  $k")
        #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
        #     println("λ_sq_vec =  $(λ_sq_vec)")
        #     println("ηk =  $(modelparams.ηk)")
        #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
        #     println("var_muk =  $(clusters[k].var_muk )")
        #     println("mk_hat =  $(clusters[k].mk_hat)")

        # end

        
    #    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
    perK_data_elbo = sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1]) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk[1])) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
    yjk_entropy_perK += entropy(clusters[k].yjk_hat)
    yjk_entropy_perK = -yjk_entropy_perK
    perK_ebloval =   perK_data_elbo + yjk_entropy_perK
    elbolog.per_k_elbo[k,iter] += perK_ebloval
    data_elbo += perK_ebloval
    end
    #     
    return data_elbo,elbolog
end
function calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    one_half_const = 1/2
    # λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
    data_elbo = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        clusters[k].cache .= 0.0
        yjk_entropy_perK = 0.0
        # @fastmath @inbounds @simd for i in 1:N
        #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
        # end
        # try
        #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
        # catch
        #     println("K =  $k")
        #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
        #     println("λ_sq_vec =  $(λ_sq_vec)")
        #     println("ηk =  $(modelparams.ηk)")
        #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
        #     println("var_muk =  $(clusters[k].var_muk )")
        #     println("mk_hat =  $(clusters[k].mk_hat)")

        # end

        
    ####    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


    # perK_data_elbo = sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1]) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk[1])) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


    # perK_data_elbo = sum( .+ . .+  .+)

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
    #     
    return data_elbo,elbolog
end
function calc_DataElbo_mpu2(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    one_half_const = 1/2
    λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
    data_elbo = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        clusters[k].cache .= 0.0
        yjk_entropy_perK = 0.0
        # @fastmath @inbounds @simd for i in 1:N
        #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
        # end
        # try
        #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
        # catch
        #     println("K =  $k")
        #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
        #     println("λ_sq_vec =  $(λ_sq_vec)")
        #     println("ηk =  $(modelparams.ηk)")
        #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
        #     println("var_muk =  $(clusters[k].var_muk )")
        #     println("mk_hat =  $(clusters[k].mk_hat)")

        # end

        
    ####    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


    perK_data_elbo = sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1]) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk[1])) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


    # perK_data_elbo = sum( .+ . .+  .+)

    # perK_data_elbo = 0.0
    # for j in 1:G
    #     perK_data_elbo += -1*one_half_const * clusters[k].Nk[1] * log(2π)
    #     perK_data_elbo += -1*one_half_const * clusters[k].Nk[1]  * log(clusters[k].σ_sq_k_hat[j])
    #     perK_data_elbo += -1*one_half_const * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat_sq[j] 
    #     perK_data_elbo +=  1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat[j] * clusters[k].κk_hat[j]
    #     perK_data_elbo +=  -1*one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].var_muk[j] 
    #     perK_data_elbo += -1* one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 - one_half_const * clusters[k].yjk_hat[j] * log(geneparams[j].λ_sq[1]) 
    #     perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].var_muk[j]
    #     perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 
    #     perK_data_elbo += clusters[k].yjk_hat[j] * log(modelparams.ηk[1]) + (1 - clusters[k].yjk_hat[j]) * log((1-modelparams.ηk[1]))
    #     perK_data_elbo += one_half_const * clusters[k].yjk_hat[j] * log(clusters[k].v_sq_k_hat[j])
    #     perK_data_elbo +=  one_half_const * clusters[k].yjk_hat[j]
    # end
    

    yjk_entropy_perK += entropy(clusters[k].yjk_hat)
    yjk_entropy_perK = -yjk_entropy_perK
    perK_ebloval =   perK_data_elbo + yjk_entropy_perK
    elbolog.per_k_elbo[k,iter] += perK_ebloval
    data_elbo += perK_ebloval
    end
    #     
    return data_elbo,elbolog
end
# assgn_entropy =  calc_Hz(rtik) 
function calc_Hz_fast3(cellpop,clusters,dataparams)
    float_type = dataparams.BitType
    N = dataparams.N
    z_entropy = 0.0
    @fastmath @inbounds @simd for i in 1:N
        z_entropy += entropy(cellpop[i].rtik)
    end
    return  -z_entropy
end

# dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
# pip_entropy = calc_Hpip(pip_kj);
# pip_entropy + dHDP_surragate_elbo

function calc_HyjkSurragateLowerBound_unconstrained(Tk,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0
    # rho_hat = sigmoid.(c)
    # omega_hat = exp.(d)
    # c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    # e_logUk = logUk_expected_value(rho_hat,omega_hat)
    # e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    # K = length(rho_hat)
    # e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    # k_vec = collect(1:K)
    # # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    # lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ e_γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  e_βk .* Tαk[1:K]
    #take Negative to find max  
    # lb_lg = sum(lb_lg_k)

    lb_lg = 0.0
    yjk_entropy = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        c_B = beta(clusters[k].gk_hat[1] * clusters[k].hk_hat[1] , (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1])
        e_βk = expectation_βk(k,clusters,modelparams)
        e_logUk = expectation_logUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
        e_log1minusUk = expectation_log1minusUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
        lb_lg += -1.0 * c_B +  (T + 1. - clusters[k].gk_hat[1] * clusters[k].hk_hat[1]) *e_logUk  +  (T * ( K + 1. - k) + γ0 - (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1]) * e_log1minusUk  +  e_βk * α0 * Tk[k]
        yjk_entropy += entropy(clusters[k].yjk_hat)
    end
    yjk_entropy = -yjk_entropy 


    return lb_lg + yjk_entropy
end

function calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0
    # rho_hat = sigmoid.(c)
    # omega_hat = exp.(d)
    # c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    # e_logUk = logUk_expected_value(rho_hat,omega_hat)
    # e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    # K = length(rho_hat)
    # e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    # k_vec = collect(1:K)
    # # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    # lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ e_γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  e_βk .* Tαk[1:K]
    #take Negative to find max  
    # lb_lg = sum(lb_lg_k)

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


# wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
function calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # T = length(c_ttprime)
    wAlloc_elbo = 0.0



    @fastmath @inbounds for t in 2:T
        # a_cttprime = 0.0
        # @fastmath @inbounds @simd for t_prime_a in t:T
        #     a_cttprime += conditionparams[t_prime_a].ctt_prime[t]
        # end

        b_cttprime = 0.0
        @fastmath @inbounds for t_prime_b in t:T
            @fastmath @inbounds @simd for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                b_cttprime += conditionparams[t_prime_b].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        c_Beta_p = -logbeta(1, modelparams.ϕ0)#c_Beta(adot_w,bdot_w)
        c_Beta_q = -(-logbeta(1, conditionparams[t-1].st_hat[1]))#c_Beta(conditionparams[t-1].awt_hat[1],conditionparams[t-1].bwt_hat[1])
        c_Beta_pq = c_Beta_p + c_Beta_q
        # adot_w_awt_hat_vec = modelparams.adot_w - conditionparams[t-1].awt_hat[1]
        ϕ0_st_hat_vec = modelparams.ϕ0- conditionparams[t-1].st_hat[1]
        e_log_tilde_wt = expectation_log_tilde_wtt(1, conditionparams[t-1].st_hat[1])
        e_log_minus_tilde_wt = expectation_log_minus_tilde_wtt(1, conditionparams[t-1].st_hat[1])


        wAlloc_elbo_t = c_Beta_pq + (1) * e_log_tilde_wt + (ϕ0_st_hat_vec + b_cttprime) * e_log_minus_tilde_wt
        wAlloc_elbo += wAlloc_elbo_t
    end
    return wAlloc_elbo
end

# α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
function calc_alphaElbo_fast3(conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # N = dataparams.N
    # K = modelparams.K
    # T = length(a_αt_hat_vec)
    α_elbo = 0.0
    # a_αt_vec = [a_α for t in 1:T]
    # b_αt_vec = [b_α for t in 1:T]
    # c_Gamma_α_p = c_Ga.(a_αt_vec,b_αt_vec)
    # c_Gamma_α_q = -c_Ga(a_αt_hat_vec,b_αt_hat_vec)
    # e_αt = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    # e_log_αt  = log_αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    @simd for t in 1:T
        αt_elbo = c_Ga(modelparams.a_α,modelparams.b_α) - c_Ga(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) + (modelparams.a_α - conditionparams[t].a_αt_hat[1]) *expectation_log_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) - (modelparams.b_α -conditionparams[t].b_αt_hat[1]) * expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
        α_elbo += αt_elbo
    end
    return α_elbo
end

# α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
# s_entropy = calc_Hs(c_ttprime_vec)
# γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
# α_elbo + s_entropy + γ_elbo

function calc_HsGammaAlphaElbo_fast3(a_γ,b_γ,a_γ_hat,b_γ_hat,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # N = dataparams.N
    # K = modelparams.K
    # T = length(a_αt_hat_vec)
    α_elbo = 0.0
    s_entropy = 0.0
    @fastmath @inbounds @simd for t in 1:T
        αt_elbo = c_Ga(modelparams.a_α,modelparams.b_α) - c_Ga(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) + (modelparams.a_α - conditionparams[t].a_αt_hat[1]) *expectation_log_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) - (modelparams.b_α -conditionparams[t].b_αt_hat[1]) * expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
        α_elbo += αt_elbo
        s_entropy += entropy(conditionparams[t].c_tt_prime)
    end

        
    s_entropy = -s_entropy

    c_Gamma_γ_p = c_Ga(a_γ,b_γ)
    c_Gamma_γ_q = -c_Ga(a_γ_hat,b_γ_hat)
    e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
    e_log_γ  = log_γ_expected_value(a_γ_hat,b_γ_hat)
    γ_elbo = c_Gamma_γ_p + c_Gamma_γ_q+(a_γ- a_γ_hat)*e_log_γ - (b_γ-b_γ_hat) * e_γ


    return α_elbo + γ_elbo + s_entropy
end


function calc_HsElbo(conditionparams, dataparams, modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # N = dataparams.N
    # K = modelparams.K
    # T = length(a_αt_hat_vec)

    s_entropy = 0.0
    @fastmath @inbounds @simd for t in 1:T
        # αt_elbo = c_Ga(modelparams.a_α, modelparams.b_α) - c_Ga(conditionparams[t].a_αt_hat[1], conditionparams[t].b_αt_hat[1]) + (modelparams.a_α - conditionparams[t].a_αt_hat[1]) * expectation_log_αt(conditionparams[t].a_αt_hat[1], conditionparams[t].b_αt_hat[1]) - (modelparams.b_α - conditionparams[t].b_αt_hat[1]) * expectation_αt(conditionparams[t].a_αt_hat[1], conditionparams[t].b_αt_hat[1])
        # α_elbo += αt_elbo
        s_entropy += entropy(conditionparams[t].c_tt_prime)
    end


    s_entropy = -s_entropy

    # c_Gamma_γ_p = c_Ga(a_γ, b_γ)
    # c_Gamma_γ_q = -c_Ga(a_γ_hat, b_γ_hat)
    # e_γ = γ_expected_value(a_γ_hat, b_γ_hat)
    # e_log_γ = log_γ_expected_value(a_γ_hat, b_γ_hat)
    # γ_elbo = c_Gamma_γ_p + c_Gamma_γ_q + (a_γ - a_γ_hat) * e_log_γ - (b_γ - b_γ_hat) * e_γ


    return s_entropy
end






function tidy_calc_Hz(rtik_mat)
    T = length(unique(rtik_mat[:,1]));
    N = length(unique(rtik_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat);
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    entropies = [entropy(rtik_mat[(i-1)*K+1:i*K,end]) for i in 1:N]
    z_entropy = sum([sum(entropies[st:en])  for (st,en) in timeranges])
    return  -z_entropy
end
function tidy_calc_Hs(c_ttprime_mat)
    T = length(unique(c_ttprime_mat[:,1]));
    s_entropies = [entropy(c_ttprime_mat[(t-1)*T+1:t*T,end]) for t in 1:T]
    s_entropy = sum(s_entropies)
    return  -s_entropy
end
function tidy_calc_Hxi(ηtkj_mat)
    T = length(unique(ηtkj_mat[:,1]));
    K = length(unique(ηtkj_mat[:,2]));
    G = length(unique(ηtkj_mat[:,3]));
    xi_entropies = [entropy(ηtkj_mat[(t-1)*2+1:t*2,end]) for t in 1:T*K*G]
    xi_entropy = sum(xi_entropies)
    return  -xi_entropy
end
function tidy_calc_SurragateLowerBound_unconstrained(ρkωkckdk_hat_mat,T,e_γ_mat,Tαk_mat)
    c_hat_vec_ = ρkωkckdk_hat_mat[:,4]
    d_hat_vec = ρkωkckdk_hat_mat[:,5]
    Tαk_ = Tαk_mat[:,2]
    e_γ_ = e_γ_mat[1,2]
    rho_hat = sigmoid.(c_hat_vec_)
    omega_hat = exp.(d_hat_vec)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ e_γ_ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  e_βk .* Tαk_[1:K]
    #take Negative to find max  
    lb_lg = sum(lb_lg_k)
    return lb_lg
end
function tidy_calc_GammaElbo(aαbαawbwaγbγ_mat,aγbγ_hat_mat)
    a_γ = aαbαawbwaγbγ_mat[1,6]
    b_γ = aαbαawbwaγbγ_mat[1,7]
    a_γ_hat =aγbγ_hat_mat[1,2]
    b_γ_hat  = aγbγ_hat_mat[1,3]
    c_Gamma_γ_p = c_Ga(a_γ,b_γ)
    c_Gamma_γ_q = -c_Ga(a_γ_hat,b_γ_hat)
    e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
    e_log_γ  = log_γ_expected_value(a_γ_hat,b_γ_hat)
    γ_elbo = c_Gamma_γ_p + c_Gamma_γ_q+(a_γ- a_γ_hat)*e_log_γ - (b_γ-b_γ_hat) * e_γ
    return γ_elbo
end
function tidy_calc_wAllocationsLowerBound(c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(c_ttprime_mat[:,1]))
    adot_w_= aαbαawbwaγbγ_mat[1,4]
    bdot_w_ = aαbαawbwaγbγ_mat[1,5]
    awt_hat_vec_ = aαtbαtawtbwt_hat_mat[2:end,4]
    bwt_hat_vec_ = aαtbαtawtbwt_hat_mat[2:end,5]
    c_ttprime_ = [c_ttprime_mat[(t-1)*T+1:t*T,end] for t in 1:T]
    wAlloc_elbo = 0.0
    c_Beta_p = c_Beta.(adot_w_,bdot_w_)
    c_Beta_q = -c_Beta.(awt_hat_vec_,bwt_hat_vec_)
    c_Beta_pq = c_Beta_p .+ c_Beta_q
    e_log_tilde_wt = log_tilde_wt_expected_value.(awt_hat_vec_,bwt_hat_vec_)
    e_log_minus_tilde_wt = log_minus_tilde_wt_expected_value.(awt_hat_vec_,bwt_hat_vec_)
    adot_w_awt_hat_vec = adot_w_ .- awt_hat_vec_
    bdot_w_bwt_hat_vec = bdot_w_ .- bwt_hat_vec_
    for t in 2:T
        a_cttprime = 0.0
        for t_prime_a  in t:T
            a_cttprime += c_ttprime_[t_prime_a][t]
        end 

        b_cttprime = 0.0
        for t_prime_b  in t:T
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                b_cttprime += c_ttprime_[t_prime_b][l]
                # sum_string *= c_string
            end
        end
        wAlloc_elbo_t = c_Beta_pq[t-1] .+ (adot_w_awt_hat_vec[t-1] .+ a_cttprime) .* e_log_tilde_wt[t-1] .+ (bdot_w_bwt_hat_vec[t-1] .+ b_cttprime) .* e_log_minus_tilde_wt[t-1]
        wAlloc_elbo += wAlloc_elbo_t
    end
    return wAlloc_elbo
end
function tidy_calc_alphaElbo(aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)
    T = length(unique(c_ttprime_mat[:,1]));
    
    a_α = aαbαawbwaγbγ_mat[1,2]
    b_α = aαbαawbwaγbγ_mat[1,3]

    a_αt_hat_vec = aαtbαtawtbwt_hat_mat[:,2]
    b_αt_hat_vec = aαtbαtawtbwt_hat_mat[:,3]
    α_elbo = 0.0
    a_αt_vec = [a_α for t in 1:T]
    b_αt_vec = [b_α for t in 1:T]
    c_Gamma_α_p = c_Ga.(a_αt_vec,b_αt_vec)
    c_Gamma_α_q = -c_Ga(a_αt_hat_vec,b_αt_hat_vec)
    e_αt = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    e_log_αt  = log_αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    for t in 1:T
        αt_elbo = c_Gamma_α_p[t] .+ c_Gamma_α_q[t] .+ (a_αt_vec[t] .- a_αt_hat_vec[t]) .*e_log_αt[t] .- (b_αt_vec[t] .-b_αt_hat_vec[t]) .* e_αt[t]
        α_elbo += αt_elbo
    end
    return α_elbo
end
function tidy_calc_DataElbo(xmat,rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(rtik_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,K)
    Nk_mat = tidy_update_Nk(rtik_mat)
    λ0_vec = λ0μ0a0b0_mat[:,2]
    μ0_vec = λ0μ0a0b0_mat[:,3]
    a0_vec = λ0μ0a0b0_mat[:,4]
    b0_vec = λ0μ0a0b0_mat[:,5]

    λ0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,3] for k in states_id]
    mk_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,4] for k in states_id]
    a0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,5] for k in states_id]
    b0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,6] for k in states_id]


    Nk = Nk_mat[:,end]
    x_ =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    rtik_ =  [rtik_mat[(i-1)*K+1:(i)*K,end] for i in 1:N];

    lin_timeranges = zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)

    halfNklogπ = 1/2 .* Nk .*log(2π) # K by 1

    halfNkoverλ0k_hat = 1/2 .* [ Nk[k] ./λ0k_hat_vec[k] for k in 1:K] # K by G
    a0b0_c_Ga =  [c_Ga.(a0_vec, b0_vec) for k in 1:K] # K by G [zeros(G) for k in 1:K]#
    a0b0_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#
    
    
    halfLogλ0 = 1/2 .* [ log.(λ0_vec) for k in 1:K] # K by G
    halfLogλ0k = 1/2 .* [ log.(λ0k_hat_vec[k]) for k in 1:K] # K by G
    halfλ0μ0_sq = 1/2 .* λ0_vec.* μ0_vec .^ 2 #G By 1
    a0kOverb0k_hat = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K] # K by G
    bs_e_τ = [(b0k_hat_vec[k]  .- b0_vec .- halfλ0μ0_sq) .*a0kOverb0k_hat[k]  for k in 1:K] #K by G
    # digamma.(a_0kj_hat) .- log.(b_0kj_hat) 
    # e_log_τkj = log_τ_expected_value.(a0k_hat_vec, b0k_hat_vec) #K by G
    as_e_log_τkj =[(1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])) for k in 1:K] #K by G
    λ0μ0 = λ0_vec .* μ0_vec #G By 1
    λ0μ0mk_a0kOverb0k_hat = [ λ0μ0 .* mk_hat_vec[k] .* a0kOverb0k_hat[k] for k in 1:K] #K by G
    halfλ0mk_sq_a0kOverb0k_hat = [1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* a0kOverb0k_hat[k]  for k in 1:K] #K by G
    one_half_const = 1/2
    halfλ0overλ0k_hat = [1/2 .* λ0_vec ./ λ0k_hat_vec[k] for k in 1:K] #K by G

    weighted_ss_kjti = [[rtik_[i][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id]
    # weighted_ss_kjti = [[[[rtik_[t][i][k] .*(x_[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
    # aa = [[[[rtik[t][i][k] .*(x_to_use[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K];
    # bb  = [[rtik_[i][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id];
    # aa2 = [[[sum(aa[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    # aa3 = [[sum(aa2[k][j]) for j in 1:G] for k in 1:K]
    #@test sum(sum(sum(bb))) == sum(sum(sum(sum(aa))))
    #@test all(sum.(sum.(bb)) .== sum.(aa3))
    #@test all(all.([ uu .≈ vv for (uu,vv) in zip(sum.(bb),aa3)]))
    #@test all(all.([ uu .== vv for (uu,vv) in zip([sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in states_id],aa3)])) 
    # [sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in state_ids]
    # weighted_ss_kjt = [[[sum(weighted_ss_kjti[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    weighted_ss_kj =  sum.(weighted_ss_kjti)#[[sum(weighted_ss_kjt[k][j]) for j in 1:G] for k in 1:K] #K by G
    halfa0kOverb0k_hat_weighted_ss_kj =  [ 1/2 .* a0kOverb0k_hat[k] .* weighted_ss_kj[k]  for k in 1:K] #K by G

    data_elbo_kj = [ halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .+ a0b0_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .- one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    data_elbo_j  = sum(data_elbo_kj)
    data_elbo  = sum(data_elbo_j)
    #     
    return data_elbo
end 
#TODO CHECK THIS MATH!!!!
function tidy_calc_DataElbo_SparseVS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(rtik_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,K)
    Nk_mat = tidy_update_Nk(rtik_mat)
    λ0_vec = λ0μ0a0b0_mat[:,2]
    μ0_vec = λ0μ0a0b0_mat[:,3]
    a0_vec = λ0μ0a0b0_mat[:,4]
    b0_vec = λ0μ0a0b0_mat[:,5]

    λ0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,3] for k in states_id]
    mk_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,4] for k in states_id]
    a0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,5] for k in states_id]
    b0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,6] for k in states_id]


    Nk = Nk_mat[:,end]
    x_ =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    rtik_ =  [rtik_mat[(i-1)*K+1:(i)*K,end] for i in 1:N];
    ηtkj_mat_tr = ηtkj_mat[1:2:end,:];
    ηtkj_vec = [ηtkj_mat_tr[(t-1)*G+1:(t)*G,end] for t in 1:T*K];
    ηtkj_ = [ηtkj_vec[t:T:end] for t in 1:T]

    lin_timeranges = zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    time_to_cellid(id) = collect(1:T)[collect(0:T-1).*  N_t  .+ 1 .<= id .<= collect(1:T).*  N_t][1]
    
    halfNklogπ = 1/2 .* Nk .*log(2π) # K by 1

    halfNkoverλ0k_hat = 1/2 .* [ Nk[k] ./λ0k_hat_vec[k] for k in 1:K] # K by G
    a0b0_c_Ga =  [c_Ga.(a0_vec, b0_vec) for k in 1:K] # K by G [zeros(G) for k in 1:K]#
    a0b0_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#
    
    
    halfLogλ0 = 1/2 .* [ log.(λ0_vec) for k in 1:K] # K by G
    halfLogλ0k = 1/2 .* [ log.(λ0k_hat_vec[k]) for k in 1:K] # K by G
    halfλ0μ0_sq = 1/2 .* λ0_vec.* μ0_vec .^ 2 #G By 1
    a0kOverb0k_hat = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K] # K by G
    bs_e_τ = [(b0k_hat_vec[k]  .- b0_vec .- halfλ0μ0_sq) .*a0kOverb0k_hat[k]  for k in 1:K] #K by G
    # digamma.(a_0kj_hat) .- log.(b_0kj_hat) 
    # e_log_τkj = log_τ_expected_value.(a0k_hat_vec, b0k_hat_vec) #K by G
    as_e_log_τkj =[(1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])) for k in 1:K] #K by G
    λ0μ0 = λ0_vec .* μ0_vec #G By 1
    λ0μ0mk_a0kOverb0k_hat = [ λ0μ0 .* mk_hat_vec[k] .* a0kOverb0k_hat[k] for k in 1:K] #K by G
    halfλ0mk_sq_a0kOverb0k_hat = [1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* a0kOverb0k_hat[k]  for k in 1:K] #K by G
    one_half_const = 1/2
    halfλ0overλ0k_hat = [1/2 .* λ0_vec ./ λ0k_hat_vec[k] for k in 1:K] #K by G

    weighted_ss_kjti = [[rtik_[i][k] .* ηtkj_[time_to_cellid(i)][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id]
    # weighted_ss_kjti = [[[[rtik_[t][i][k] .*(x_[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
    # aa = [[[[rtik[t][i][k] .*(x_to_use[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K];
    # bb  = [[rtik_[i][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id];
    # aa2 = [[[sum(aa[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    # aa3 = [[sum(aa2[k][j]) for j in 1:G] for k in 1:K]
    #@test sum(sum(sum(bb))) == sum(sum(sum(sum(aa))))
    #@test all(sum.(sum.(bb)) .== sum.(aa3))
    #@test all(all.([ uu .≈ vv for (uu,vv) in zip(sum.(bb),aa3)]))
    #@test all(all.([ uu .== vv for (uu,vv) in zip([sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in states_id],aa3)])) 
    # [sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in state_ids]
    # weighted_ss_kjt = [[[sum(weighted_ss_kjti[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    weighted_ss_kj =  sum.(weighted_ss_kjti)#[[sum(weighted_ss_kjt[k][j]) for j in 1:G] for k in 1:K] #K by G
    halfa0kOverb0k_hat_weighted_ss_kj =  [ 1/2 .* a0kOverb0k_hat[k] .* weighted_ss_kj[k]  for k in 1:K] #K by G

    data_elbo_kj = [ halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .+ a0b0_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .- one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    data_elbo_j  = sum(data_elbo_kj)
    data_elbo  = sum(data_elbo_j)
    #     
    return data_elbo
end 
function tidy_calc_DataElbo_VS1(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(rtik_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,K)
    Nk_mat = tidy_update_Nk(rtik_mat)
    λ0_vec = λ0μ0a0b0_mat[:,2]
    μ0_vec = λ0μ0a0b0_mat[:,3]
    a0_vec = λ0μ0a0b0_mat[:,4]
    b0_vec = λ0μ0a0b0_mat[:,5]

    λ0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,3] for k in states_id]
    mk_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,4] for k in states_id]
    a0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,5] for k in states_id]
    b0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,6] for k in states_id]


    Nk = Nk_mat[:,end]
    x_ =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    rtik_ =  [rtik_mat[(i-1)*K+1:(i)*K,end] for i in 1:N];
    ηtkj_mat_tr = ηtkj_mat[1:2:end,:];
    ηtkj_vec = [ηtkj_mat_tr[(t-1)*G+1:(t)*G,end] for t in 1:T*K];
    ηtkj_ = [ηtkj_vec[t:T:end] for t in 1:T]

    lin_timeranges = zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    time_to_cellid(id) = collect(1:T)[collect(0:T-1).*  N_t  .+ 1 .<= id .<= collect(1:T).*  N_t][1]
    
    halfNklogπ = 1/2 .* Nk .*log(2π) # K by 1

    halfNkoverλ0k_hat = 1/2 .* [ Nk[k] ./λ0k_hat_vec[k] for k in 1:K] # K by G
    a0b0_c_Ga =  [c_Ga.(a0_vec, b0_vec) for k in 1:K] # K by G [zeros(G) for k in 1:K]#
    a0b0_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#
    
    
    halfLogλ0 = 1/2 .* [ log.(λ0_vec) for k in 1:K] # K by G
    halfLogλ0k = 1/2 .* [ log.(λ0k_hat_vec[k]) for k in 1:K] # K by G
    halfλ0μ0_sq = 1/2 .* λ0_vec.* μ0_vec .^ 2 #G By 1
    a0kOverb0k_hat = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K] # K by G
    bs_e_τ = [(b0k_hat_vec[k]  .- b0_vec .- halfλ0μ0_sq) .*a0kOverb0k_hat[k]  for k in 1:K] #K by G
    # digamma.(a_0kj_hat) .- log.(b_0kj_hat) 
    # e_log_τkj = log_τ_expected_value.(a0k_hat_vec, b0k_hat_vec) #K by G
    as_e_log_τkj =[(1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])) for k in 1:K] #K by G
    λ0μ0 = λ0_vec .* μ0_vec #G By 1
    λ0μ0mk_a0kOverb0k_hat = [ λ0μ0 .* mk_hat_vec[k] .* a0kOverb0k_hat[k] for k in 1:K] #K by G
    halfλ0mk_sq_a0kOverb0k_hat = [1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* a0kOverb0k_hat[k]  for k in 1:K] #K by G
    one_half_const = 1/2
    halfλ0overλ0k_hat = [1/2 .* λ0_vec ./ λ0k_hat_vec[k] for k in 1:K] #K by G

    weighted_ss_kjti = [[rtik_[i][k] .* ηtkj_[time_to_cellid(i)][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id]
    # weighted_ss_kjti = [[[[rtik_[t][i][k] .*(x_[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
    # aa = [[[[rtik[t][i][k] .*(x_to_use[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K];
    # bb  = [[rtik_[i][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id];
    # aa2 = [[[sum(aa[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    # aa3 = [[sum(aa2[k][j]) for j in 1:G] for k in 1:K]
    #@test sum(sum(sum(bb))) == sum(sum(sum(sum(aa))))
    #@test all(sum.(sum.(bb)) .== sum.(aa3))
    #@test all(all.([ uu .≈ vv for (uu,vv) in zip(sum.(bb),aa3)]))
    #@test all(all.([ uu .== vv for (uu,vv) in zip([sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in states_id],aa3)])) 
    # [sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in state_ids]
    # weighted_ss_kjt = [[[sum(weighted_ss_kjti[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    weighted_ss_kj =  sum.(weighted_ss_kjti)#[[sum(weighted_ss_kjt[k][j]) for j in 1:G] for k in 1:K] #K by G
    halfa0kOverb0k_hat_weighted_ss_kj =  [ 1/2 .* a0kOverb0k_hat[k] .* weighted_ss_kj[k]  for k in 1:K] #K by G

    data_elbo_kj = [ halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .+ a0b0_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .- one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    data_elbo_j  = sum(data_elbo_kj)
    data_elbo  = sum(data_elbo_j)
    #     
    return data_elbo
end 
function tidy_calc_DataElbo_VS2(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(rtik_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,K)
    Nk_mat = tidy_update_Nk(rtik_mat)
    λ0_vec = λ0μ0a0b0_mat[:,2]
    μ0_vec = λ0μ0a0b0_mat[:,3]
    a0_vec = λ0μ0a0b0_mat[:,4]
    b0_vec = λ0μ0a0b0_mat[:,5]

    λ0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,3] for k in states_id]
    mk_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,4] for k in states_id]
    a0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,5] for k in states_id]
    b0k_hat_vec = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k*G,6] for k in states_id]


    Nk = Nk_mat[:,end]
    x_ =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    rtik_ =  [rtik_mat[(i-1)*K+1:(i)*K,end] for i in 1:N];
    ηtkj_mat_tr = ηtkj_mat[1:2:end,:];
    ηtkj_vec = [ηtkj_mat_tr[(t-1)*G+1:(t)*G,end] for t in 1:T*K];
    ηtkj_ = [ηtkj_vec[t:T:end] for t in 1:T]

    lin_timeranges = zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    time_to_cellid(id) = collect(1:T)[collect(0:T-1).*  N_t  .+ 1 .<= id .<= collect(1:T).*  N_t][1]
    
    halfNklogπ = 1/2 .* Nk .*log(2π) # K by 1

    halfNkoverλ0k_hat = 1/2 .* [ Nk[k] ./λ0k_hat_vec[k] for k in 1:K] # K by G
    a0b0_c_Ga =  [c_Ga.(a0_vec, b0_vec) for k in 1:K] # K by G [zeros(G) for k in 1:K]#
    a0b0_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#
    
    
    halfLogλ0 = 1/2 .* [ log.(λ0_vec) for k in 1:K] # K by G
    halfLogλ0k = 1/2 .* [ log.(λ0k_hat_vec[k]) for k in 1:K] # K by G
    halfλ0μ0_sq = 1/2 .* λ0_vec.* μ0_vec .^ 2 #G By 1
    a0kOverb0k_hat = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K] # K by G
    bs_e_τ = [(b0k_hat_vec[k]  .- b0_vec .- halfλ0μ0_sq) .*a0kOverb0k_hat[k]  for k in 1:K] #K by G
    # digamma.(a_0kj_hat) .- log.(b_0kj_hat) 
    # e_log_τkj = log_τ_expected_value.(a0k_hat_vec, b0k_hat_vec) #K by G
    as_e_log_τkj =[(1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])) for k in 1:K] #K by G
    λ0μ0 = λ0_vec .* μ0_vec #G By 1
    λ0μ0mk_a0kOverb0k_hat = [ λ0μ0 .* mk_hat_vec[k] .* a0kOverb0k_hat[k] for k in 1:K] #K by G
    halfλ0mk_sq_a0kOverb0k_hat = [1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* a0kOverb0k_hat[k]  for k in 1:K] #K by G
    one_half_const = 1/2
    halfλ0overλ0k_hat = [1/2 .* λ0_vec ./ λ0k_hat_vec[k] for k in 1:K] #K by G

    weighted_ss_kjti = [[rtik_[i][k] .* ηtkj_[time_to_cellid(i)][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id]
    # weighted_ss_kjti = [[[[rtik_[t][i][k] .*(x_[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
    # aa = [[[[rtik[t][i][k] .*(x_to_use[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K];
    # bb  = [[rtik_[i][k] .* (x_[i] .- mk_hat_vec[k]) .^2 for i in cell_ids] for k in states_id];
    # aa2 = [[[sum(aa[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    # aa3 = [[sum(aa2[k][j]) for j in 1:G] for k in 1:K]
    #@test sum(sum(sum(bb))) == sum(sum(sum(sum(aa))))
    #@test all(sum.(sum.(bb)) .== sum.(aa3))
    #@test all(all.([ uu .≈ vv for (uu,vv) in zip(sum.(bb),aa3)]))
    #@test all(all.([ uu .== vv for (uu,vv) in zip([sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in states_id],aa3)])) 
    # [sum([sum(bb[k][st:en])  for (st,en) in lin_timeranges]) for k in state_ids]
    # weighted_ss_kjt = [[[sum(weighted_ss_kjti[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    weighted_ss_kj =  sum.(weighted_ss_kjti)#[[sum(weighted_ss_kjt[k][j]) for j in 1:G] for k in 1:K] #K by G
    halfa0kOverb0k_hat_weighted_ss_kj =  [ 1/2 .* a0kOverb0k_hat[k] .* weighted_ss_kj[k]  for k in 1:K] #K by G

    data_elbo_kj = [ halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .+ a0b0_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .- one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    data_elbo_j  = sum(data_elbo_kj)
    data_elbo  = sum(data_elbo_j)
    #     
    return data_elbo
end 
function tidy_elbo_calc(xmat,rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(xmat[:,1]))
    # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    data_elbo = tidy_calc_DataElbo(xmat,rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    # assgn_entropy =  calc_Hz(rtik)
    assgn_entropy = tidy_calc_Hz(rtik_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    e_γ_mat =  tidy_γ_expected_value(aγbγ_hat_mat)
    # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    dHDP_surragate_elbo = tidy_calc_SurragateLowerBound_unconstrained(ρkωkckdk_hat_mat,T,e_γ_mat,Tαk_mat)##NOT EXACTLY THE SAME (different optimization function)
    # s_entropy = calc_Hs(c_ttprime_vec)
    s_entropy = tidy_calc_Hs(c_ttprime_mat)
    # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    wAlloc_elbo = tidy_calc_wAllocationsLowerBound(c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    γ_elbo = tidy_calc_GammaElbo(aαbαawbwaγbγ_mat,aγbγ_hat_mat)
    # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    α_elbo = tidy_calc_alphaElbo(aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)
    elbo_val = data_elbo +  dHDP_surragate_elbo + assgn_entropy + s_entropy + wAlloc_elbo + γ_elbo + α_elbo
    return elbo_val
end
#CURRENTLY SETS DATA ELBO TO 0!!!!
function tidy_elbo_calc_SparseVS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(xmat[:,1]))
    # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    data_elbo = tidy_calc_DataElbo_SparseVS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    # assgn_entropy =  calc_Hz(rtik)
    assgn_entropy = tidy_calc_Hz(rtik_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    e_γ_mat =  tidy_γ_expected_value(aγbγ_hat_mat)
    # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    dHDP_surragate_elbo = tidy_calc_SurragateLowerBound_unconstrained(ρkωkckdk_hat_mat,T,e_γ_mat,Tαk_mat)##NOT EXACTLY THE SAME (different optimization function)
    # s_entropy = calc_Hs(c_ttprime_vec)
    s_entropy = tidy_calc_Hs(c_ttprime_mat)
    # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    wAlloc_elbo = tidy_calc_wAllocationsLowerBound(c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    γ_elbo = tidy_calc_GammaElbo(aαbαawbwaγbγ_mat,aγbγ_hat_mat)
    # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    α_elbo = tidy_calc_alphaElbo(aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)
    importance_elbo = tidy_calc_Hxi(ηtkj_mat)
    data_elbo = 0
    elbo_val = data_elbo +  dHDP_surragate_elbo + assgn_entropy + s_entropy + wAlloc_elbo + γ_elbo + α_elbo + importance_elbo
    return elbo_val
end
function tidy_elbo_calc_VS1(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(xmat[:,1]))
    # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    data_elbo = tidy_calc_DataElbo_VS1(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    # assgn_entropy =  calc_Hz(rtik)
    assgn_entropy = tidy_calc_Hz(rtik_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    e_γ_mat =  tidy_γ_expected_value(aγbγ_hat_mat)
    # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    dHDP_surragate_elbo = tidy_calc_SurragateLowerBound_unconstrained(ρkωkckdk_hat_mat,T,e_γ_mat,Tαk_mat)##NOT EXACTLY THE SAME (different optimization function)
    # s_entropy = calc_Hs(c_ttprime_vec)
    s_entropy = tidy_calc_Hs(c_ttprime_mat)
    # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    wAlloc_elbo = tidy_calc_wAllocationsLowerBound(c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    γ_elbo = tidy_calc_GammaElbo(aαbαawbwaγbγ_mat,aγbγ_hat_mat)
    # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    α_elbo = tidy_calc_alphaElbo(aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)
    importance_elbo = tidy_calc_Hxi(ηtkj_mat)
    data_elbo = 0
    elbo_val = data_elbo +  dHDP_surragate_elbo + assgn_entropy + s_entropy + wAlloc_elbo + γ_elbo + α_elbo + importance_elbo
    return elbo_val
end
function tidy_elbo_calc_VS2(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(xmat[:,1]))
    # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    data_elbo = tidy_calc_DataElbo_VS2(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    # assgn_entropy =  calc_Hz(rtik)
    assgn_entropy = tidy_calc_Hz(rtik_mat) ##NOT EXACTLY THE SAME (Numerical stability?)
    e_γ_mat =  tidy_γ_expected_value(aγbγ_hat_mat)
    # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    dHDP_surragate_elbo = tidy_calc_SurragateLowerBound_unconstrained(ρkωkckdk_hat_mat,T,e_γ_mat,Tαk_mat)##NOT EXACTLY THE SAME (different optimization function)
    # s_entropy = calc_Hs(c_ttprime_vec)
    s_entropy = tidy_calc_Hs(c_ttprime_mat)
    # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    wAlloc_elbo = tidy_calc_wAllocationsLowerBound(c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
    # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    γ_elbo = tidy_calc_GammaElbo(aαbαawbwaγbγ_mat,aγbγ_hat_mat)
    # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    α_elbo = tidy_calc_alphaElbo(aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)
    importance_elbo = tidy_calc_Hxi(ηtkj_mat)
    data_elbo = 0
    elbo_val = data_elbo +  dHDP_surragate_elbo + assgn_entropy + s_entropy + wAlloc_elbo + γ_elbo + α_elbo + importance_elbo
    return elbo_val
end

