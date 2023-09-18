# function update_Nc_t(c_ttprime)
#     T = length(c_ttprime)
#     Nc_t = Vector{Float64}(undef,T+1)
#     for t in 1:T
#         Nc_t[t] = sum(c_ttprime[t])
#     end
#     Nc_t[T+1] = 0
#     return Nc_t
# end
function update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ)
    T = length(e_log_π)
    K = length(e_log_τ)
    C_t = [length(el) for el in e_τ_μ]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = [[[e_log_π[t][k] - 1/2 * Glog + 1/2*e_log_τ[k] - 1/2*e_τ_μ[t][i][k]   for k in 1:K] for i in 1:C_t[t] ] for t in 1:T]

    for t in 1:T
        numcells = C_t[t]
        rtik[t] = Vector{Vector{Float64}}(undef,numcells)
        for i in 1:numcells
            rtik[t][i] = Vector{Float64}(undef,K)
            val_sum = StatsFuns.logsumexp(ptik_tilde[t][i])
            val = exp.(ptik_tilde[t][i] .- val_sum)
            # shifted_val = val .+ eps(1.0)
            # rtik[t][i] = shifted_val ./ sum(shifted_val)#val#
            rtik[t][i] = val#
        end
    end

    return rtik
end
function update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τ)
    C_t = [length(el) for el in e_τ_μ]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] - 1/2 * Glog + 1/2*e_log_τ[k] - 1/2*e_τ_μ[t][i][k]
            end
            ptik_tilde_ti[i] = ptik_tilde_tik
        end
        ptik_tilde[t] = ptik_tilde_ti
    end
    # ptik_tilde = [[[e_log_π[t][k] - 1/2 * Glog + 1/2*e_log_τ[k] - 1/2*e_τ_μ[t][i][k]   for k in 1:K] for i in 1:C_t[t] ] for t in 1:T]

    for t in 1:T
        numcells = C_t[t]
        rtik[t] = Vector{Vector{Float64}}(undef,numcells)
        for i in 1:numcells
            rtik[t][i] = Vector{Float64}(undef,K)
            val_sum = StatsFuns.logsumexp(ptik_tilde[t][i])
            val = exp.(ptik_tilde[t][i] .- val_sum)
            # shifted_val = val .+ eps(1.0)
            # rtik[t][i] = shifted_val ./ sum(shifted_val)#val#
            rtik[t][i] = val#
        end
    end

    return rtik
end
function update_Ntk(rtik)
    T = length(rtik)
    C_t = [length(el) for el in rtik]
    K = length(rtik[1][1])
    Ntk =  Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        Nt_vec[1:K] =  sum(rtik[t]) 
        Ntk[t] = Nt_vec
    end
    return Ntk
end
function update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,α0)
    T = length(Ntk)
    e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    θ_hat = [α0 .* e_βk_vec .+ Ntk[t] for t in 1:T]
    return θ_hat
end
function update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime)
    T = length(Ntk)
    K = length(rhok_hat_vec)
    e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    e_αt_vec = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    θ_hat = Vector{Vector{Float64}}(undef,T) 
    for t in 1:T
        updated_Ntk = Ntk .*  c_ttprime[t]
        sum_updated_Ntk = sum(updated_Ntk)
        e_αt_βk = e_αt_vec[t] .* e_βk_vec
        θ_hat[t] = e_αt_βk .+ sum_updated_Ntk
    end
    return θ_hat
end
function update_Nk(rtik)
    Nk = sum(sum.(rtik))#sum(Ntk)[1:end-1]
    return Nk
end
function update_N(rtik,v_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(v_tikj[1][1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error = Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Ntik_signal = Vector{Float64}(undef,G)
                Ntik_error = Vector{Float64}(undef,G)
                for j in 1:G
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * v_tikj[t][i][k][j][2]
                end
                Nti_signal[k] = Ntik_signal
                Nti_error[k] = Ntik_error
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_x_hat_error_vs_forloops(x,N_error)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N_error[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    push!(value, x[t][i][j] * N_error[t][i][k][j])  
                end
            end
        end
        x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_err
end
function update_x_hatk_signal_vs_forloops(x,N_signal)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N_signal[1][1])
    G = length(x[1][1])
    x_hatk_signal = Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                push!(value, x[t][i] .* N_signal[t][i][k])
            end
        end
        x_hatk_signal[k] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hatk_signal
end
function update_x_hat_sq_error_vs_forloops(x,N_error)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N_error[1][1])
    G = length(x[1][1])
    x_hat_sq_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    push!(value, x[t][i][j]^2 * N_error[t][i][k][j])  
                end
            end
        end
        x_hat_sq_err[j] = sum(value)
    end
    return x_hat_sq_err
end
function update_x_hatk_sq_signal_vs_forloops(x,N_signal)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N_signal[1][1])
    G = length(x[1][1])
    x_hatk_sq_signal = Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                push!(value, x[t][i] .^2 .* N_signal[t][i][k])
            end
        end
        x_hatk_sq_signal[k] = sum(value)
    end
    return x_hatk_sq_signal
end
function update_xbar_k(x,rtik)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    xbar_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    xbar_tk = [[ sum(xbar_tik[k][t]) for t in 1:T] for k in 1:K]
    xbar_k = [sum(xbar_tk[k]) for k in 1:K]
    return xbar_k
end
function update_sk(x,xbar_k,rtik)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])

    s_tik = [[[ rtik[t][i][k] .* (x[t][i] .- xbar_k[k]) .^2 for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    s_tk = [[ sum(s_tik[k][t]) for t in 1:T] for k in 1:K]
    s_k = [sum(s_tk[k]) for k in 1:K]
    return s_k
end
function update_x_hat_sq_k(x,rtik)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])

    x_hat_sq_tik = [[[ rtik[t][i][k] .* (x[t][i]) .^2 for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_sq_tk = [[ sum(x_hat_sq_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_sq_k = [sum(x_hat_sq_tk[k]) for k in 1:K]
    return x_hat_sq_k
end
function update_x_hat_k(x,rtik)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_k
end
function update_λ0k_hat(λ0_vec,Nk) 
    K = length(Nk)
    λ0k_hat_vec = [λ0_vec .+ Nk[k] for k in 1:K]
    return λ0k_hat_vec
end
function update_a0k_hat_usingXhat(a0_vec,Nk)
    K = length(Nk) 
    a0k_hat_vec = [ a0_vec .+ 1/2 * (Nk[k]+1) for k in 1:K]
    return a0k_hat_vec
end
function update_a0k_hat(a0_vec,Nk)
    K = length(Nk) 
    a0k_hat_vec = [ a0_vec .+ 1/2 * Nk[k] for k in 1:K]
    return a0k_hat_vec
end
function update_mk_hat(λ0_vec,μ0_vec, Nk,xbar_k)
    K = length(Nk)
    Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_vec .* μ0_vec
    denom = [λ0_vec .+ Nk[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    mkj_hat = [ (λ0_μ0 .+ Nk_xbar_k[k]) ./denom[k] for k in 1:K]
    return mkj_hat
end
function update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
    K = length(Nk)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_vec .* μ0_vec
    denom = [λ0_vec .+ Nk[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    mkj_hat = [ (λ0_μ0 .+ x_hat_k[k]) ./denom[k] for k in 1:K]
    return mkj_hat
end
function update_b0k_hat(b0_vec,λ0_vec,μ0_vec, Nk,xbar_k,sk)
    K = length(Nk)
    denom = [λ0_vec .+ Nk[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    # numer = [(λ0_vec .* Nk[k] .* (xbar_k[k] .- μ0_vec) .^2) ./denom[k]  for k in 1:K ]
    numer = [(λ0_vec .* Nk[k] .* (xbar_k[k] .- μ0_vec) .^2) for k in 1:K ]
    ssd = [numer[k] ./ denom[k] for k in 1:K]
    half_sk_ssd =  1/2 .* [sk[k] .+ ssd[k] for k in 1:K] 
    # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    return  b0k_hat_vec
end
function update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)
    K = length(Nk)
    denom = [λ0_vec .+ Nk[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    μ0_sq_vec = μ0_vec .^2
    μ0λ0_vec =  λ0_vec .* μ0_vec
    μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    numer = [(x_hat_k[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    ssd = [numer[k] ./ denom[k] for k in 1:K]
    half_sk_ssd =  1/2 .* [x_hat_sq_k[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    return  b0k_hat_vec
end
function update_Tk(θ_hat)
    e_log_π = log_π_expected_value(θ_hat)
    Tk = sum(e_log_π)
    return Tk
end
function update_rho_omega_hat(rhok_hat_vec,omegak_hat_vec,T,γ,α0,Tk;optim_max_iter=100000)
    c_hat = StatsFuns.logit.(rhok_hat_vec)
    d_hat =  log.(omegak_hat_vec)
    cd_vec_args = [c_hat , d_hat]
    cd_vec_args0  = permutedims(reduce(hcat,cd_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,γ,α0,Tk)
    gg_uncon! = g_unconstrained_closure!(T,γ,α0,Tk)
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(), Optim.Options(iterations = optim_max_iter))
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter)) #,linesearch=Optim.LineSearches.BackTracking() 
    # optimize(LB_LG_unconstrained,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))

    # @debug Optim.converged(lb_lg_results)
    new_rho_hat = sigmoid.(lb_lg_results.res.minimizer[1,:])
    new_omega_hat = exp.(lb_lg_results.res.minimizer[2,:])

    new_c_hat = StatsFuns.logit.(new_rho_hat)
    new_d_hat = log.(new_omega_hat)
    return new_rho_hat,new_omega_hat, new_c_hat,new_d_hat
end
function update_rho_omega_hat(rhok_hat_vec,omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=100000)
    c_hat = StatsFuns.logit.(rhok_hat_vec)
    d_hat =  log.(omegak_hat_vec)
    cd_vec_args = [c_hat , d_hat]
    cd_vec_args0  = permutedims(reduce(hcat,cd_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,e_γ,Tαk)
    gg_uncon! = g_unconstrained_closure!(T,e_γ,Tαk)
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(), Optim.Options(iterations = optim_max_iter))
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter)) #,linesearch=Optim.LineSearches.BackTracking() 
    # optimize(LB_LG_unconstrained,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))

    # @debug Optim.converged(lb_lg_results)
    new_rho_hat = sigmoid.(lb_lg_results.res.minimizer[1,:])
    new_omega_hat = exp.(lb_lg_results.res.minimizer[2,:])

    new_c_hat = StatsFuns.logit.(new_rho_hat)
    new_d_hat = log.(new_omega_hat)
    return new_rho_hat,new_omega_hat, new_c_hat,new_d_hat
end
function update_awt_hat(adot_w, c_ttprime)
    T = length(c_ttprime)
    awt_hat_vec = Vector{Float64}(undef,T-1)
    for t in 2:T
        awt_hat = 0.0
        awt_hat += adot_w
        # sum_string = ""
        # sum_string *= "adot_w "
        for t_prime in t:T
            # c_string = "+ c$(t_prime)$(t) "
            awt_hat += c_ttprime[t_prime][t]
            # sum_string *= c_string
        end
        awt_hat_vec[t-1]=awt_hat
        # println(sum_string)
    end
    return awt_hat_vec
end
function update_bwt_hat(bdot_w, c_ttprime)
    #Test with c_ttprime = [[1.0,0.0,0.0,0.0],[0.5,0.5,0.0,0.0],[0.25,0.25,0.5,0.0],[0.25,0.25,0.25,0.25]]
    T = length(c_ttprime)
    bwt_hat_vec = Vector{Float64}(undef,T-1)
    for t in 2:T
        bwt_hat = 0.0
        bwt_hat += bdot_w
        # sum_string = ""
        # sum_string *= "bdot_w "
        for t_prime in t:T
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                bwt_hat += c_ttprime[t_prime][l]
                # sum_string *= c_string
            end
        end
        bwt_hat_vec[t-1]=bwt_hat
        # println(sum_string)
    end
    return bwt_hat_vec
end
function update_c_ttprime(awt_hat,bwt_hat,rtik,θ_hat)
    e_w_ttprime = log_w_ttprime_expected_value(awt_hat,bwt_hat)
    e_log_π = log_π_expected_value(θ_hat)
    Ntk = update_Ntk(rtik)#sum.(rtik)
    T = length(Ntk)
    # adjusted_e_log_π
    c_ttprime = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        adjusted_e_log_π_k = [Ntk[t] .* el for el in e_log_π[1:t]]
        adjusted_e_log_π = sum.(adjusted_e_log_π_k)
        p_c_ttprime = e_w_ttprime[t] .+ adjusted_e_log_π
        # c_ttprime[t] = norm_weights(p_c_ttprime)
        c_ttprime[t] = vcat(norm_weights(p_c_ttprime),zeros(Float64,T-t))
        # c_full = push!(c_trunc,zeros(Float64,T-t))
        # c_ttprime[t] = c_full
    end
    return c_ttprime
end
function update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat)
    K = length(rhok_hat_vec)
    T = length(θ_hat)
    a_α_hat = a_α + K
    a_α_hat_vec = Vector{Float64}(undef,T)
    b_α_hat_vec = Vector{Float64}(undef,T)
    e_log_π = log_π_expected_value(θ_hat)
    e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    for t in 1:T
        a_α_hat_vec[t] = a_α_hat
        e_βk_log_π = [e_βk_vec[k] .* e_log_π[t][k]  for k in 1:K+1]
        b_α_hat_vec[t] = b_α - sum(e_βk_log_π)
    end
    return a_α_hat_vec,b_α_hat_vec
end
function update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    K = length(rhok_hat_vec)
    e_log_minus_uk = log1minusUk_expected_value(rhok_hat_vec, omegak_hat_vec)
    a_γ_hat = a_γ + K
    b_γ_hat = b_γ - sum(e_log_minus_uk)
    return a_γ_hat, b_γ_hat
end
function update_Tαk(θ_hat,a_α_hat_vec,b_α_hat_vec)
    e_log_π = log_π_expected_value(θ_hat)
    e_αt_vec = αt_expected_value.(a_α_hat_vec,b_α_hat_vec)
    e_αt_log_π =  e_αt_vec .*  e_log_π
    Tαk = sum(e_αt_log_π)
    return Tαk
end



#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################

function get_gene_PIP25_fast3(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10,float_type=nothing)#,N_k=nothing
    if isnothing(float_type)
        float_type =eltype(x[1][1])
    end
    G = length(x[1][1])
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    pip_kj  = Vector{Vector{float_type}}(undef,K)
    N_k = Vector{float_type}(undef,K)
    # unnormalized_weights  = Vector{Vector{Float64}}(undef,K)
    # # println("here")
    # for k in 1:K
    #     unnormalized_weights[k]=zeros(Float64,G)
    # end
    gene_vec = Vector{float_type}(undef,G)
    zeros_vec = zeros(float_type,G)
    m_null = zeros(float_type,G)
    a0_null = ones(float_type,G)
    for j in 1:G
        a0_null[j] = null_precision * a0_null[j]
    end
    b0_null = ones(float_type,G)
    for k in 1:K
        unnormalized_weights=copy(zeros_vec)
        n_tik = 0.0
        for t in 1:T
            for i in 1:C_t[t]
                # @views unnormalized_weights .+= norm_weights(cell_ll_scores[t][i][k] .- null_cell_ll_scores[t][i][k]) .* rtik[t][i][k]


                # cell_unnormalized_weights = log_norm_pdf2(x[t][i],mk_hat_vec[k],a0k_hat_vec[k] , b0k_hat_vec[k]) - log_norm_pdf2(x[t][i],m_null,a0_null , b0_null)
                # unnormalized_weights += norm_weights(cell_unnormalized_weights) * rtik[t][i][k]

                # cell_unnormalized_weights = 
                # unnormalized_weights .+= norm_weights(log_norm_pdf2(x[t][i],mk_hat_vec[k],a0k_hat_vec[k] , b0k_hat_vec[k]) .- log_norm_pdf2(x[t][i],m_null,a0_null , b0_null)) .* rtik[t][i][k]
                unnormalized_weights = accumulate_deviance_log_norm_pdf2!(unnormalized_weights,x[t][i],mk_hat_vec[k],a0k_hat_vec[k],b0k_hat_vec[k],m_null,a0_null,b0_null, rtik[t][i][k])
                n_tik += rtik[t][i][k]
            end
        end
        N_k[k] = n_tik
        normalized_weights = gene_vec
        for j in 1:G
            normalized_weights[j] = unnormalized_weights[j] / N_k[k]
            if isnan(normalized_weights[j])
                normalized_weights[j] = 0.0
            end
        end
        if iszero(sum(normalized_weights))
            for j in 1:G
                normalized_weights[j] = 1.0
            end
        end
        pip_kj[k] = normToProb(fix_nan_or_allzero2!(normalized_weights))
    end
    return pip_kj
end
function get_gene_PIP25_fast3!(cellpop,clusters,dataparams,modelparams)#,N_k=nothing,nulldist::NullDistributionFeatures,
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K

    # N_t = dataparams.N_t
    # lin_address = dataparams.LinearAddress


    # null_precision = modelparams.null_precision
    # T = length(x)
    # C_t = [length(el) for el in x]
    # K = length(rtik[1][1])
    # pip_kj  = Vector{Vector{float_type}}(undef,K)
    # N_k = Vector{float_type}(undef,K)
    # unnormalized_weights  = Vector{Vector{Float64}}(undef,K)
    # # println("here")
    # for k in 1:K
    #     unnormalized_weights[k]=zeros(Float64,G)
    # end

    # gene_vec = Vector{float_type}(undef,G)
    # zeros_vec = zeros(float_type,G)


    # m_null = nulldist.m_null
    # a0_null = nulldist.a0_null
    # for j in 1:G
    #     a0_null[j] = null_precision * a0_null[j]
    # end
    # b0_null = ones(float_type,G)
    for k in 1:K
        # clusters[k]._reset!(clusters[k].pip_k, clusters[k].BitType)
        clusters[k].pip_k .= 0.0
        n_tik = 0.0
        for i in 1:N
            # println("here")
            # println(typeof(cellpop[i].rtik))
            # println(typeof(clusters[k].pip_k))
            # println(typeof(clusters[k].k))
            # println(typeof(clusters[k].cache))
            # println(typeof(clusters[k]))
            # println(typeof(cellpop[i]))
            _accumulate_deviance_log_norm_pdf3!(clusters[k],cellpop[i],dataparams)
            n_tik += cellpop[i].rtik[clusters[k].k]
        end
        clusters[k].Nk[1] = n_tik
        for j in 1:G
            clusters[k].pip_k[j] = clusters[k].pip_k[j] / clusters[k].Nk[1]
            if isnan(clusters[k].pip_k[j])
                clusters[k].pip_k[j] = 0.0
            end
        end
        if iszero(sum(clusters[k].pip_k))
            for j in 1:G
                clusters[k].pip_k[j] = 1.0
            end
        end
        normToProb3!(clusters[k].pip_k)
        # clusters[k].pip_k= normToProb(clusters[k].pip_k)
    end
    return clusters
end
function get_gene_PIP25_fast3_mt!(cellpop,clusters,dataparams,modelparams)#,N_k=nothing,nulldist::NullDistributionFeatures,
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K

    # N_t = dataparams.N_t
    # lin_address = dataparams.LinearAddress


    # null_precision = modelparams.null_precision
    # T = length(x)
    # C_t = [length(el) for el in x]
    # K = length(rtik[1][1])
    # pip_kj  = Vector{Vector{float_type}}(undef,K)
    # N_k = Vector{float_type}(undef,K)
    # unnormalized_weights  = Vector{Vector{Float64}}(undef,K)
    # # println("here")
    # for k in 1:K
    #     unnormalized_weights[k]=zeros(Float64,G)
    # end

    # gene_vec = Vector{float_type}(undef,G)
    # zeros_vec = zeros(float_type,G)


    # m_null = nulldist.m_null
    # a0_null = nulldist.a0_null
    # for j in 1:G
    #     a0_null[j] = null_precision * a0_null[j]
    # end
    # b0_null = ones(float_type,G)
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].pip_k, clusters[k].BitType)
        clusters[k].pip_k .= 0.0
        n_tik = 0.0
        for i in 1:N
            # println("here")
            # println(typeof(cellpop[i].rtik))
            # println(typeof(clusters[k].pip_k))
            # println(typeof(clusters[k].k))
            # println(typeof(clusters[k].cache))
            # println(typeof(clusters[k]))
            # println(typeof(cellpop[i]))
            _accumulate_deviance_log_norm_pdf3!(clusters[k],cellpop[i],dataparams)
            n_tik += cellpop[i].rtik[clusters[k].k]
        end
        clusters[k].Nk[1] = n_tik
        for j in 1:G
            clusters[k].pip_k[j] = clusters[k].pip_k[j] / clusters[k].Nk[1]
            if isnan(clusters[k].pip_k[j])
                clusters[k].pip_k[j] = 0.0
            end
        end
        if iszero(sum(clusters[k].pip_k))
            for j in 1:G
                clusters[k].pip_k[j] = 1.0
            end
        end
        normToProb3!(clusters[k].pip_k)
        # clusters[k].pip_k= normToProb(clusters[k].pip_k)
    end
    return clusters
end
function accumulate_deviance_log_norm_pdf3!(unnormalized_weights,x,mk_hat_k,a0k_hat_k, b0k_hat_k,m_null,a0_null , b0_null,rtik;float_type=nothing)
    G = length(x)
    if isnothing(float_type)
        float_type =eltype(x)
    end
    # logpi =  log(2π)
    log_norm_pdf_ = Vector{float_type}(undef,G)
    for j in 1:G
        # log_norm_pdf_[j] = (1/2 * log(a0k_hat_k[j]/b0k_hat_k[j]) -  1/2 * ((x[j]-mk_hat_k[j])^2 * (a0k_hat_k[j]/b0k_hat_k[j])) - 1/2 * logpi) - (1/2 * log(a0_null[j]/b0_null[j]) -  1/2 * ((x[j]-m_null[j])^2 * (a0_null[j]/b0_null[j])) - 1/2 * logpi) 
        log_norm_pdf_[j] = log_norm_pdf3(x[j],mk_hat_k[j],a0k_hat_k[j],b0k_hat_k[j])- log_norm_pdf3(x[j],m_null[j],a0_null[j],b0_null[j])
    end
    normed_weight = norm_weights3(log_norm_pdf_)
    for j in 1:G
        unnormalized_weights[j] +=  normed_weight[j]* rtik
    end
    return unnormalized_weights
end
function _accumulate_deviance_log_norm_pdf3!(cluster,cell,dataparams)
    # println("here")

    G = dataparams.G
    # float_type = dataparams.BitType
    for j in 1:G
        # log_norm_pdf_[j] = (1/2 * log(a0k_hat_k[j]/b0k_hat_k[j]) -  1/2 * ((x[j]-mk_hat_k[j])^2 * (a0k_hat_k[j]/b0k_hat_k[j])) - 1/2 * logpi) - (1/2 * log(a0_null[j]/b0_null[j]) -  1/2 * ((x[j]-m_null[j])^2 * (a0_null[j]/b0_null[j])) - 1/2 * logpi) 
        cluster.cache[j] = log_norm_pdf3(cell.x[j],cluster.mk_hat[j],cluster.a0k_hat[j],cluster.b0k_hat[j])- cell.x_null[j]
    end
    
    norm_weights3!(cluster.cache)
    for j in 1:G
        cluster.pip_k[j] +=  cluster.cache[j]* cell.rtik[cluster.k]
    end
    # cluster._reset!(cluster.cache, cluster.BitType)
    return cluster
end

function update_yjk!_depracated(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]  

    for k in 1:K

        clusters[k].yjk_hat .= 0.0
        logodds = log10((modelparams.ηk[1])  / (1 - (modelparams.ηk[1])))
        logp10 = log(10)*logodds 
        clusters[k].yjk_hat .=  logp10 .+ log.(sqrt.(clusters[k].v_sq_k_hat) ./ (sqrt.(λ_sq_vec))) .+ 0.5 .* (clusters[k].mk_hat) .^2 ./clusters[k].v_sq_k_hat 

        sigmoidNorm!(clusters[k].yjk_hat) 
    end

end

function update_yjk!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    # λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]  

    for k in 1:K

        clusters[k].yjk_hat .= 0.0
        logodds = log10((modelparams.ηk[1])  / (1 - (modelparams.ηk[1])))
        logp10 = log(10)*logodds 
        for j in 1:G
            clusters[k].yjk_hat[j] =  logp10 + log(sqrt(clusters[k].v_sq_k_hat[j]) / (sqrt(geneparams[j].λ_sq[1]))) + 0.5 * (clusters[k].mk_hat[j]) ^2 /clusters[k].v_sq_k_hat[j] 
        end
        sigmoidNorm!(clusters[k].yjk_hat) 
    end

end
function update_yjk_mt!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    # λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]  

    Threads.@threads for k in 1:K

        clusters[k].yjk_hat .= 0.0
        logodds = log10((modelparams.ηk[1])  / (1 - (modelparams.ηk[1])))
        logp10 = log(10)*logodds 
        for j in 1:G
            clusters[k].yjk_hat[j] =  logp10 + log(sqrt(clusters[k].v_sq_k_hat[j]) / (sqrt(geneparams[j].λ_sq[1]))) + 0.5 * (clusters[k].mk_hat[j]) ^2 /clusters[k].v_sq_k_hat[j] 
        end
        sigmoidNorm!(clusters[k].yjk_hat) 
    end

end

function update_yjk_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    # λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]  
    if mt_mode =="full"
        Threads.@threads for k in 1:K

            clusters[k].yjk_hat .= 0.0
            logodds = log10((modelparams.ηk[1])  / (1 - (modelparams.ηk[1])))
            logp10 = log(10)*logodds 
            for j in 1:G
                clusters[k].yjk_hat[j] =  logp10 + log(sqrt(clusters[k].v_sq_k_hat[j]) / (sqrt(geneparams[j].λ_sq[1]))) + 0.5 * (clusters[k].mk_hat[j]) ^2 /clusters[k].v_sq_k_hat[j] 
            end
            sigmoidNorm!(clusters[k].yjk_hat) 
        end
    else
        for k in 1:K

            clusters[k].yjk_hat .= 0.0
            logodds = log10((modelparams.ηk[1])  / (1 - (modelparams.ηk[1])))
            logp10 = log(10)*logodds 
            for j in 1:G
                clusters[k].yjk_hat[j] =  logp10 + log(sqrt(clusters[k].v_sq_k_hat[j]) / (sqrt(geneparams[j].λ_sq[1]))) + 0.5 * (clusters[k].mk_hat[j]) ^2 /clusters[k].v_sq_k_hat[j] 
            end
            sigmoidNorm!(clusters[k].yjk_hat) 
        end
    end


end
function update_ηk!(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K


    for k in 1:K

        modelparams.ηk[1]  = 0.0
        modelparams.ηk[1]  = sum(clusters[k].yjk_hat) / sum(1 .- clusters[k].yjk_hat)
        sigmoidNorm!(modelparams.ηk) 
    end

end
# .* sqrt.(clusters[k].σ_sq_k_hat)
# clusters[k].yjk_hat .= log(modelparams.ηk / (1 - modelparams.ηk)) .+ log.(sqrt.(clusters[k].v_sq_k_hat) ./ sqrt.(λ_sq_vec)) .- 0.5 .* (clusters[k].mk_hat) .^2 ./clusters[k].v_sq_k_hat  .+ 0.5

function update_rtik!(cellpop,clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    for i in 1:N
        t = first(dataparams.LinearAddress[i])
        adjust_e_log_π_tk3!(t,conditionparams)
        cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
        for k in 1:K
            expectation_log_normal_l_j!(clusters[k],cellpop[i], dataparams)
            cell_gene_sums  = sum(cellpop[i].cache)
            cellpop[i].rtik[k] = conditionparams[t].e_log_π_t_cache[k] + cell_gene_sums
        end
        norm_weights3!(cellpop[i].rtik)
    end
    return cellpop
end

function update_rtik_mt!(cellpop,clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for i in 1:N
        t = first(dataparams.LinearAddress[i])
        adjust_e_log_π_tk3!(t,conditionparams)
        cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
        for k in 1:K
            expectation_log_normal_l_j!(clusters[k],cellpop[i], dataparams)
            cell_gene_sums  = sum(cellpop[i].cache)
            cellpop[i].rtik[k] = conditionparams[t].e_log_π_t_cache[k] + cell_gene_sums
        end
        norm_weights3!(cellpop[i].rtik)
    end
    return cellpop
end

function update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams;mt_mode=nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    if mt_mode == "full" || mt_mode == "optimal"
        # println("Here")
        Threads.@threads for i in 1:N
            t = first(dataparams.LinearAddress[i])
            adjust_e_log_π_tk3!(t,conditionparams)
            cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
            # cellpop[i].rtik .= 0.0
            for k in 1:K
                expectation_log_normal_l_j!(clusters[k],cellpop[i], dataparams)
                cell_gene_sums = 0.0
                for el in cellpop[i].cache
                    cell_gene_sums+=el
                end
                # cell_gene_sums  = sum(cellpop[i].cache)
                cellpop[i].rtik[k] = conditionparams[t].e_log_π_t_cache[k] + cell_gene_sums
            end
            norm_weights3!(cellpop[i].rtik)
        end
    else
        # float_type = dataparams.BitType
        # G = dataparams.G
        # T = dataparams.T
        # N = dataparams.N
        # K = modelparams.K
        for i in 1:N
            t = first(dataparams.LinearAddress[i])
            adjust_e_log_π_tk3!(t,conditionparams)
            cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
            # cellpop[i].rtik .= 0.0
            for k in 1:K
                expectation_log_normal_l_j!(clusters[k],cellpop[i], dataparams)
                cell_gene_sums = 0.0
                for el in cellpop[i].cache
                    cell_gene_sums+=el
                end
                # cell_gene_sums  = sum(cellpop[i].cache)
                cellpop[i].rtik[k] = conditionparams[t].e_log_π_t_cache[k] + cell_gene_sums
            end
            norm_weights3!(cellpop[i].rtik)
        end
    end

    return cellpop
end

function update_rtik_vs25_fast3_mt!(cellpop,clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    # emptycells = [Vector{Vector{Float64}}(undef,C_t[t]) for t in 1:T]
    # emptyclusters =  Vector{Float64}(undef,K)
    # ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    # θ_hat_t_sum = Vector{Float64}(undef,T)

    # for t in 1:T
    #     conditionparams[t].θ_hat_t_sum[1] = sum(conditionparams[t].θ_hat_t )
    # end
    # logpi = Glog/G
    Threads.@threads for i in 1:N
        t = first(dataparams.LinearAddress[i])
        adjust_e_log_π_tk3!(t,conditionparams)
        cellpop[i].rtik .= 0.0#_reset!(cellpop[i].rtik,cellpop[i].BitType)
        
        # rtik_atom = Threads.Atomic{Float64}(0)
        for k in 1:K
            # cellpop[i].cache .= 0.0
            expectation_log_normal_l_j3!(clusters[k],cellpop[i], dataparams)
            cell_gene_sums  = sum(cellpop[i].cache)
            cellpop[i].rtik[k] += conditionparams[t].e_log_π_t_cache[k] + cell_gene_sums
            # Threads.atomic_add!(rtik_atom, conditionparams[t].e_log_π_t_cache[k])
            # Threads.atomic_add!(rtik_atom, cell_gene_sums)
            # cellpop[i].rtik[k] = convert(float_type,rtik_atom)
        end
        norm_weights3!(cellpop[i].rtik)
    end
    return cellpop
end


# Ntk = update_Ntk(rtik)
function update_Ntk!(cellpop,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # Ntk =  Vector{Vector{Float64}}(undef,T)
    #@fastmath @inbounds
    @inbounds for t in 1:T
        _reset!(conditionparams[t].Ntk,float_type)
        st,en = dataparams.TimeRanges[t]
        # Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        # Nt_vec[1:K] = 
        @inbounds for i in st:en
            @inbounds for k in 1:K #@fastmath 
                conditionparams[t].Ntk[k] += cellpop[i].rtik[k]
            end
            # conditionparams[t].Ntk[1:K] += cellpop[i].rtik
        end
        # @views conditionparams[t].Ntk[1:K] .= sum([el.rtik for el in cellpop[st:en]])
        # 
    end
    return conditionparams
end
function update_Ntk_fast3_mt!(cellpop,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # Ntk =  Vector{Vector{Float64}}(undef,T)
    #@fastmath @inbounds
    Threads.@threads  for t in 1:T
        _reset!(conditionparams[t].Ntk,float_type)
        st,en = dataparams.TimeRanges[t]
        # Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        # Nt_vec[1:K] = 
        for i in st:en
            for k in 1:K #@fastmath 
                conditionparams[t].Ntk[k] += cellpop[i].rtik[k]
            end
            # conditionparams[t].Ntk[1:K] += cellpop[i].rtik
        end
        # @views conditionparams[t].Ntk[1:K] .= sum([el.rtik for el in cellpop[st:en]])
        # 
    end
    return conditionparams
end
function update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams;mt_mode=nothing)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    if mt_mode == "full"
        Threads.@threads  for t in 1:T
            _reset!(conditionparams[t].Ntk,float_type)
            st,en = dataparams.TimeRanges[t]
            # Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
            # Nt_vec[1:K] = 
            for i in st:en
                for k in 1:K #@fastmath 
                    conditionparams[t].Ntk[k] += cellpop[i].rtik[k]
                end
                # conditionparams[t].Ntk[1:K] += cellpop[i].rtik
            end
            # @views conditionparams[t].Ntk[1:K] .= sum([el.rtik for el in cellpop[st:en]])
            # 
        end
    else
        @inbounds for t in 1:T
            _reset!(conditionparams[t].Ntk,float_type)
            st,en = dataparams.TimeRanges[t]
            # Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
            # Nt_vec[1:K] = 
            @inbounds for i in st:en
                @inbounds for k in 1:K #@fastmath 
                    conditionparams[t].Ntk[k] += cellpop[i].rtik[k]
                end
                # conditionparams[t].Ntk[1:K] += cellpop[i].rtik
            end
            # @views conditionparams[t].Ntk[1:K] .= sum([el.rtik for el in cellpop[st:en]])
            # 
        end
    end

    return conditionparams
end
# mod_Ntk = zeros(float_type,Kplus) 
# old_Ntk = [zeros(float_type,Kplus) for t in 1:T]
# for t in 1:T
#     conditionparams[t].Ntk .= copy(update_Ntk(rtik_)[t])
# end
# print([conditionparams[t].Ntk for t in 1:T][1])
# print(old_Ntk[1])
# print(mod_Ntk)
function calculate_sparse_restart!(old_Ntk,mod_Ntk,conditionparams,dataparams,sparse_restart_smallest_n)
    T = dataparams.T
    Kplus = length(conditionparams[1].Ntk)
    for t in 1:T
        mod_Ntk .= 0.0
        old_Ntk[t] .= copy(conditionparams[t].Ntk)
        mod_Ntk .= vectorfloor!(copy(conditionparams[t].Ntk))
        non_zero_indices = broadcast(!,iszero.(mod_Ntk))
        zero_indices = iszero.(mod_Ntk)
        ordered_clusters = sortperm(conditionparams[t].Ntk[non_zero_indices],rev=true)
        conditionparams[t].Ntk[zero_indices] .= 0.0
        if length(ordered_clusters) > sparse_restart_smallest_n
            @views conditionparams[t].Ntk[non_zero_indices][ordered_clusters][end-sparse_restart_smallest_n+1:end] .= 0.0
        end
    end
    return old_Ntk,conditionparams
end
function reset_Ntk_sparse_restart!(old_Ntk,conditionparams,dataparams)
    T = dataparams.T
    for t in 1:T
        conditionparams[t].Ntk  .= copy(old_Ntk[t])
        old_Ntk[t] .=0.0
    end
    return old_Ntk,conditionparams
end

function update_c_ttprime!(conditionparams, dataparams, modelparams)
    float_type = dataparams.BitType
    # e_w_ttprime = log_w_ttprime_expected_value(awt_hat,bwt_hat)
    log_π_expected_value_fast3!(conditionparams, dataparams, modelparams)
    # Ntk = update_Ntk(rtik)#sum.(rtik)
    T = dataparams.T
    Kplus = modelparams.K + 1
    # adjusted_e_log_π
    # c_ttprime = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        _reset!(conditionparams[t].c_tt_prime, float_type)
        for tt in 1:t
            ctt_k_accum = 0.0
            for k in 1:Kplus
                ctt_k_accum += conditionparams[t].Ntk[k] * conditionparams[tt].e_log_π_t_cache[k]
            end
            ctt_w_accum = 0.0
            if !isone(t)
                if isone(tt)
                    e_log_tilde_wt = 0.0
                else
                    e_log_tilde_wt = log_tilde_wt_expected_value(1.0, conditionparams[tt-1].st_hat[1])#e_log_tilde_wt_vec[tt-1]
                end
                sum_e_log_minus_tilde_wt = 0.0
                if tt != t
                    for tprime in tt:t-1
                        sum_e_log_minus_tilde_wt += expectation_log_minus_tilde_wtt(1.0, conditionparams[tprime].st_hat[1])
                    end
                end
                ctt_w_accum += sum_e_log_minus_tilde_wt + e_log_tilde_wt
            end
            conditionparams[t].c_tt_prime[tt] = ctt_k_accum + ctt_w_accum
        end
        norm_weights3!(t, conditionparams[t].c_tt_prime)
        # for tt_prime in t+1:T
        #     conditionparams[t].c_tt_prime[tt_prime] = 0.0
        # end

        # adjusted_e_log_π_k = [Ntk[t] .* el for el in e_log_π[1:t]]
        # adjusted_e_log_π = sum.(adjusted_e_log_π_k)
        # p_c_ttprime = e_w_ttprime[t] .+ adjusted_e_log_π
        # # c_ttprime[t] = norm_weights(p_c_ttprime)
        # c_ttprime[t] = vcat(norm_weights(p_c_ttprime),zeros(Float64,T-t))
        # c_full = push!(c_trunc,zeros(Float64,T-t))
        # c_ttprime[t] = c_full
    end
    return conditionparams
end
function update_d_hat!(clusters, conditionparams, dataparams, modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    # e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    # e_αt_vec = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    # θ_hat = Vector{Vector{Float64}}(undef,T) 
    @inbounds for t in 1:T
        _reset!(conditionparams[t].d_hat_t, float_type)
        @inbounds for k in 1:Kplus
            e_βk = expectation_βk(k, clusters, modelparams)
            updated_Ntk_sum = 0.0
            @inbounds for tt in 1:T
                updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k]
            end
            conditionparams[t].d_hat_t[k] = modelparams.α0 * e_βk + updated_Ntk_sum
        end
        # updated_Ntk = Ntk .*  c_ttprime[t]
        # sum_updated_Ntk = sum(updated_Ntk)
        # e_αt_βk = e_αt_vec[t] .* e_βk_vec
        # θ_hat[t] = e_αt_βk .+ sum_updated_Ntk
    end
    return conditionparams
end
function update_d_hat_sum!(conditionparams, dataparams)
    float_type = dataparams.BitType
    T = dataparams.T
    for t in 1:T
        conditionparams[t].d_hat_t_sum[1] = sum(conditionparams[t].d_hat_t)
    end
    return conditionparams
end
function update_d_hat_mt!(clusters, conditionparams, dataparams, modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    # e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    # e_αt_vec = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    # θ_hat = Vector{Vector{Float64}}(undef,T) 
    Threads.@threads for t in 1:T
        _reset!(conditionparams[t].d_hat_t, float_type)
        @inbounds for k in 1:Kplus
            e_βk = expectation_βk(k, clusters, modelparams)
            updated_Ntk_sum = 0.0
            @inbounds for tt in 1:T
                updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k]
            end
            conditionparams[t].d_hat_t[k] = modelparams.α0 * e_βk + updated_Ntk_sum
        end
        # updated_Ntk = Ntk .*  c_ttprime[t]
        # sum_updated_Ntk = sum(updated_Ntk)
        # e_αt_βk = e_αt_vec[t] .* e_βk_vec
        # θ_hat[t] = e_αt_βk .+ sum_updated_Ntk
    end
    return conditionparams
end
function update_d_hat_sum_mt!(conditionparams, dataparams)
    float_type = dataparams.BitType
    T = dataparams.T
    Threads.@threads for t in 1:T
        conditionparams[t].d_hat_t_sum[1] = sum(conditionparams[t].d_hat_t)
    end
    return conditionparams
end

function update_d_hat_mpu!(clusters, conditionparams, dataparams, modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1

    if mt_mode == "full"
        Threads.@threads for t in 1:T
            _reset!(conditionparams[t].d_hat_t, float_type)
            @inbounds for k in 1:Kplus
                e_βk = expectation_βk(k, clusters, modelparams)
                updated_Ntk_sum = 0.0
                @inbounds for tt in 1:T
                    updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k]
                end
                conditionparams[t].d_hat_t[k] = modelparams.α0 * e_βk + updated_Ntk_sum
            end
        end
    else
        @inbounds for t in 1:T
            _reset!(conditionparams[t].d_hat_t, float_type)
            @inbounds for k in 1:Kplus
                e_βk = expectation_βk(k, clusters, modelparams)
                updated_Ntk_sum = 0.0
                @inbounds for tt in 1:T
                    updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k]
                end
                conditionparams[t].d_hat_t[k] = modelparams.α0 * e_βk + updated_Ntk_sum
            end
        end
    end
    return conditionparams
end
function update_d_hat_sum_mpu!(conditionparams, dataparams;mt_mode = nothing)
    float_type = dataparams.BitType
    T = dataparams.T
    if mt_mode == "full"
        Threads.@threads for t in 1:T
            conditionparams[t].d_hat_t_sum[1] = sum(conditionparams[t].d_hat_t)
        end
    else
        for t in 1:T
            conditionparams[t].d_hat_t_sum[1] = sum(conditionparams[t].d_hat_t)
        end
    end
    return conditionparams
end
function update_c_ttprime_fast3_mt!(conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # e_w_ttprime = log_w_ttprime_expected_value(awt_hat,bwt_hat)
    log_π_expected_value_fast3!(conditionparams,dataparams,modelparams)
    # Ntk = update_Ntk(rtik)#sum.(rtik)
    T = dataparams.T
    Kplus = modelparams.K + 1
    # adjusted_e_log_π
    # c_ttprime = Vector{Vector{Float64}}(undef,T)
    Threads.@threads for t in 1:T
        _reset!(conditionparams[t].c_tt_prime,float_type)
        for tt in 1:t
            ctt_k_accum = 0.0
            for k in 1:Kplus
                ctt_k_accum += conditionparams[t].Ntk[k] * conditionparams[tt].e_log_π_t_cache[k]
            end
            ctt_w_accum = 0.0 
            if !isone(t)
                if isone(tt)
                    e_log_tilde_wt = 0.0
                else
                    e_log_tilde_wt = log_tilde_wt_expected_value(conditionparams[tt-1].awt_hat[1],conditionparams[tt-1].bwt_hat[1])#e_log_tilde_wt_vec[tt-1]
                end
                sum_e_log_minus_tilde_wt = 0.0
                if tt != t
                    for tprime in tt:t-1
                        sum_e_log_minus_tilde_wt += expectation_log_minus_tilde_wtt(conditionparams[tprime].awt_hat[1],conditionparams[tprime].bwt_hat[1])
                    end
                end
                ctt_w_accum += sum_e_log_minus_tilde_wt + e_log_tilde_wt
            end
            conditionparams[t].c_tt_prime[tt] = ctt_k_accum + ctt_w_accum
        end
        norm_weights3!(t,conditionparams[t].c_tt_prime)
        # for tt_prime in t+1:T
        #     conditionparams[t].c_tt_prime[tt_prime] = 0.0
        # end
        
        # adjusted_e_log_π_k = [Ntk[t] .* el for el in e_log_π[1:t]]
        # adjusted_e_log_π = sum.(adjusted_e_log_π_k)
        # p_c_ttprime = e_w_ttprime[t] .+ adjusted_e_log_π
        # # c_ttprime[t] = norm_weights(p_c_ttprime)
        # c_ttprime[t] = vcat(norm_weights(p_c_ttprime),zeros(Float64,T-t))
        # c_full = push!(c_trunc,zeros(Float64,T-t))
        # c_ttprime[t] = c_full
    end
    return conditionparams
end
function update_c_ttprime_mpu!(conditionparams,dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    log_π_expected_value_fast3!(conditionparams,dataparams,modelparams)
    T = dataparams.T
    Kplus = modelparams.K + 1
    

    if mt_mode == "full"
        Threads.@threads for t in 1:T
            _reset!(conditionparams[t].c_tt_prime,float_type)
            for tt in 1:t
                ctt_k_accum = 0.0
                for k in 1:Kplus
                    ctt_k_accum += conditionparams[t].Ntk[k] * conditionparams[tt].e_log_π_t_cache[k]
                end
                ctt_w_accum = 0.0 
                if !isone(t)
                    if isone(tt)
                        e_log_tilde_wt = 0.0
                    else
                        e_log_tilde_wt = log_tilde_wt_expected_value(conditionparams[tt-1].awt_hat[1],conditionparams[tt-1].bwt_hat[1])#e_log_tilde_wt_vec[tt-1]
                    end
                    sum_e_log_minus_tilde_wt = 0.0
                    if tt != t
                        for tprime in tt:t-1
                            sum_e_log_minus_tilde_wt += expectation_log_minus_tilde_wtt(conditionparams[tprime].awt_hat[1],conditionparams[tprime].bwt_hat[1])
                        end
                    end
                    ctt_w_accum += sum_e_log_minus_tilde_wt + e_log_tilde_wt
                end
                conditionparams[t].c_tt_prime[tt] = ctt_k_accum + ctt_w_accum
            end
            norm_weights3!(t,conditionparams[t].c_tt_prime)
        end
    else
        for t in 1:T
            _reset!(conditionparams[t].c_tt_prime, float_type)
            for tt in 1:t
                ctt_k_accum = 0.0
                for k in 1:Kplus
                    ctt_k_accum += conditionparams[t].Ntk[k] * conditionparams[tt].e_log_π_t_cache[k]
                end
                ctt_w_accum = 0.0
                if !isone(t)
                    if isone(tt)
                        e_log_tilde_wt = 0.0
                    else
                        e_log_tilde_wt = log_tilde_wt_expected_value(1.0, conditionparams[tt-1].st_hat[1])#e_log_tilde_wt_vec[tt-1]
                    end
                    sum_e_log_minus_tilde_wt = 0.0
                    if tt != t
                        for tprime in tt:t-1
                            sum_e_log_minus_tilde_wt += expectation_log_minus_tilde_wtt(1.0, conditionparams[tprime].st_hat[1])
                        end
                    end
                    ctt_w_accum += sum_e_log_minus_tilde_wt + e_log_tilde_wt
                end
                conditionparams[t].c_tt_prime[tt] = ctt_k_accum + ctt_w_accum
            end
            norm_weights3!(t, conditionparams[t].c_tt_prime)
        end
    end
    return conditionparams
end

# θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
function update_θ_hat_fast3!(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    # e_αt_vec = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    # θ_hat = Vector{Vector{Float64}}(undef,T) 
    @inbounds for t in 1:T
        _reset!(conditionparams[t].θ_hat_t,float_type)
        e_αt = expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
        @inbounds for k in 1:Kplus
            e_βk = expectation_βk(k,clusters,modelparams)
            updated_Ntk_sum = 0.0
            @inbounds for tt in 1:T
                updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k] 
            end
            conditionparams[t].θ_hat_t[k] = e_αt*e_βk + updated_Ntk_sum
        end
        # updated_Ntk = Ntk .*  c_ttprime[t]
        # sum_updated_Ntk = sum(updated_Ntk)
        # e_αt_βk = e_αt_vec[t] .* e_βk_vec
        # θ_hat[t] = e_αt_βk .+ sum_updated_Ntk
    end
    return conditionparams
end
function update_θ_hat_fast3_mt!(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    # e_αt_vec = αt_expected_value.(a_αt_hat_vec,b_αt_hat_vec)
    # θ_hat = Vector{Vector{Float64}}(undef,T) 
    Threads.@threads for t in 1:T
        _reset!(conditionparams[t].θ_hat_t,float_type)
        e_αt = expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
        @inbounds for k in 1:Kplus
            e_βk = expectation_βk(k,clusters,modelparams)
            updated_Ntk_sum = 0.0
            @inbounds for tt in 1:T
                updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k] 
            end
            conditionparams[t].θ_hat_t[k] = e_αt*e_βk + updated_Ntk_sum
        end
        # updated_Ntk = Ntk .*  c_ttprime[t]
        # sum_updated_Ntk = sum(updated_Ntk)
        # e_αt_βk = e_αt_vec[t] .* e_βk_vec
        # θ_hat[t] = e_αt_βk .+ sum_updated_Ntk
    end
    return conditionparams
end
function update_θ_hat_sum_fast3!(conditionparams,dataparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    # G = dataparams.G
    T = dataparams.T
    # N = dataparams.N
    # K = modelparams.K
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    # emptycells = [Vector{Vector{Float64}}(undef,C_t[t]) for t in 1:T]
    # emptyclusters =  Vector{Float64}(undef,K)
    # ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    # θ_hat_t_sum = Vector{Float64}(undef,T)
    for t in 1:T
        conditionparams[t].θ_hat_t_sum[1] = sum(conditionparams[t].θ_hat_t )
    end
    return conditionparams
end
function update_θ_hat_sum_fast3_mt!(conditionparams,dataparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    # G = dataparams.G
    T = dataparams.T
    # N = dataparams.N
    # K = modelparams.K
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    # emptycells = [Vector{Vector{Float64}}(undef,C_t[t]) for t in 1:T]
    # emptyclusters =  Vector{Float64}(undef,K)
    # ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    # θ_hat_t_sum = Vector{Float64}(undef,T)
    Threads.@threads for t in 1:T
        conditionparams[t].θ_hat_t_sum[1] = sum(conditionparams[t].θ_hat_t )
    end
    return conditionparams
end



function update_Nk!(cellpop,clusters,dataparams,modelparams)#,N_k=nothing,nulldist::NullDistributionFeatures,
    float_type = dataparams.BitType
    N = dataparams.N
    K = modelparams.K

    for k in 1:K
        clusters[k].Nk .= 0.0
        Nk_sum = 0.0
        for i in 1:N
            Nk_sum += cellpop[i].rtik[k]
        end
        clusters[k].Nk[1] = Nk_sum
    end
    return clusters
end
function update_Nk_mt!(cellpop,clusters,dataparams,modelparams)#,N_k=nothing,nulldist::NullDistributionFeatures,
    float_type = dataparams.BitType
    N = dataparams.N
    K = modelparams.K

    Threads.@threads for k in 1:K
        clusters[k].Nk .= 0.0
        Nk_sum = 0.0
        for i in 1:N
            Nk_sum += cellpop[i].rtik[k]
        end
        clusters[k].Nk[1] = Nk_sum
    end
    return clusters
end
function update_Nk_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    N = dataparams.N
    K = modelparams.K
    if mt_mode == "full"
        Threads.@threads for k in 1:K
            clusters[k].Nk .= 0.0
            Nk_sum = 0.0
            for i in 1:N
                Nk_sum += cellpop[i].rtik[k]
            end
            clusters[k].Nk[1] = Nk_sum
        end
    else
        for k in 1:K
            clusters[k].Nk .= 0.0
            Nk_sum = 0.0
            for i in 1:N
                Nk_sum += cellpop[i].rtik[k]
            end
            clusters[k].Nk[1] = Nk_sum
        end
    end
    return clusters
end

# Nkj = update_Nkj25_fast(rtik, pip_kj);
function update_Nkj25_fast3!(cellpop,clusters,dataparams,modelparams)#,N_k=nothing,nulldist::NullDistributionFeatures,
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    # T = dataparams.T
    N = dataparams.N
    K = modelparams.K

    for k in 1:K
        # _reset!(clusters[k].Nkj, float_type)
        clusters[k].Nkj.= 0.0
        Nk_sum = 0.0
        for i in 1:N
        # for t in 1:T
        #     # t = first(dataparams.LinearAddress[i])
        #     st,en = dataparams.TimeRanges[t]
        #     #     # Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        #     #     # Nt_vec[1:K] = 
        #     @inbounds for i in st:en
            Nk_sum += cellpop[i].rtik[k]
        end
        # clusters[k].Nkj .+=   Nk_sum .* clusters[k].pip_k
        for j in 1:G
            clusters[k].Nkj[j] +=   Nk_sum * clusters[k].pip_k[j]
        end
            # end

    end
    return clusters
end

function update_Nkj25_fast3_mt!(cellpop,clusters,dataparams,modelparams)#,N_k=nothing,nulldist::NullDistributionFeatures,
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    # T = dataparams.T
    N = dataparams.N
    K = modelparams.K

    Threads.@threads for k in 1:K
        # _reset!(clusters[k].Nkj, float_type)
        clusters[k].Nkj.= 0.0
        Nk_sum = 0.0
        for i in 1:N
        # for t in 1:T
        #     # t = first(dataparams.LinearAddress[i])
        #     st,en = dataparams.TimeRanges[t]
        #     #     # Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        #     #     # Nt_vec[1:K] = 
        #     @inbounds for i in st:en
            Nk_sum += cellpop[i].rtik[k]
        end
        # clusters[k].Nkj .+=   Nk_sum .* clusters[k].pip_k
        for j in 1:G
            clusters[k].Nkj[j] +=   Nk_sum * clusters[k].pip_k[j]
        end
            # end

    end
    return clusters
end
function update_x_hat_k!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    @inbounds for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k]
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_k_mt!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k]
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode=nothing)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    if mt_mode == "full" || mt_mode == "optimal"
        Threads.@threads for k in 1:K
            clusters[k].x_hat .= 0.0
            @inbounds @fastmath for i in 1:N
                clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k]
            end
    
        end
    else
        @inbounds for k in 1:K
            clusters[k].x_hat .= 0.0
            @inbounds @fastmath for i in 1:N
                clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k]
            end
        end
    end

    return clusters
    # return x_hat_k
end
# x_hat_k = update_x_hat_k25_fast(x,rtik,pip_kj);
function update_x_hat_k25_fast3!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    @inbounds for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k] .* clusters[k].pip_k
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_k25_fast3_mt!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k] .* clusters[k].pip_k
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_sq_k!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    @inbounds for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat_sq .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k]
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_sq_k_mt!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat_sq .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k]
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode=nothing)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    if mt_mode == "full" || mt_mode == "optimal"
        Threads.@threads for k in 1:K
            clusters[k].x_hat_sq .= 0.0
            @inbounds @fastmath for i in 1:N
                clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k]
            end
    
        end
    else
        @inbounds for k in 1:K
            clusters[k].x_hat_sq .= 0.0
            @inbounds @fastmath for i in 1:N
                clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k]
            end
    
        end
    end
    return clusters
end
# x_hat_sq_k = update_x_hat_sq_k25_fast(x,rtik,pip_kj);  
function update_x_hat_sq_k25_fast3!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    @inbounds for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat_sq .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k] .* clusters[k].pip_k
        end

    end
    return clusters
    # return x_hat_k
end
function update_x_hat_sq_k25_fast3_mt!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].x_hat, float_type)
        clusters[k].x_hat_sq .= 0.0
        @inbounds @fastmath for i in 1:N
            # @inbounds for j in 1:G
            #     clusters[k].x_hat[j] +=   cellpop[i].x[j] * cellpop[i].rtik[k] * clusters[k].pip_k[j]
            # end
            clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k] .* clusters[k].pip_k
        end

    end
    return clusters
    # return x_hat_k
end
# λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
function update_λ0k_hat_fast3!(clusters,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].λ0k_hat.= 0.0
        # for j in 1:G
            # clusters[k].λ0k_hat[j] +=  modelparams.λ0_vec[j] + clusters[k].Nkj[j]
        # end
        clusters[k].λ0k_hat .+=  modelparams.λ0_vec .+ clusters[k].Nkj
    end
    # λ0k_hat_vec = [λ0_vec .+ Nk[k] for k in 1:K]
    return clusters
end
function update_λ0k_hat_fast3_mt!(clusters,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].λ0k_hat.= 0.0
        # for j in 1:G
            # clusters[k].λ0k_hat[j] +=  modelparams.λ0_vec[j] + clusters[k].Nkj[j]
        # end
        clusters[k].λ0k_hat .+=  modelparams.λ0_vec .+ clusters[k].Nkj
    end
    # λ0k_hat_vec = [λ0_vec .+ Nk[k] for k in 1:K]
    return clusters
end
# a0k_hat_vec = update_a0k_hat_usingXhat25_fast(a0_vec,Nkj) 
function update_a0k_hat_usingXhat25_fast3!(clusters,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    # a0k_hat_vec = [ a0_vec .+ 1/2 .* (Nkj[k] .+1 ) for k in 1:K]
    for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].a0k_hat.= 0.0

        clusters[k].a0k_hat .+=  modelparams.a0_vec .+ 1/2 .*  (clusters[k].Nkj .+1 )
    end

    return clusters
end
function update_a0k_hat_usingXhat25_fast3_mt!(clusters,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    # a0k_hat_vec = [ a0_vec .+ 1/2 .* (Nkj[k] .+1 ) for k in 1:K]
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].a0k_hat.= 0.0

        clusters[k].a0k_hat .+=  modelparams.a0_vec .+ 1/2 .*  (clusters[k].Nkj .+1 )
    end

    return clusters
end

  
function update_mk_hat!_depracated(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]  
    for k in 1:K
        clusters[k].mk_hat .= 0.0
        clusters[k].mk_hat .+=  (λ_sq_vec .* clusters[k].x_hat) ./  (clusters[k].Nk .* λ_sq_vec  .+ clusters[k].σ_sq_k_hat  )
    end

    return clusters
end

function update_mk_hat!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].mk_hat .= 0.0
        for j in 1:G
            clusters[k].mk_hat[j] +=  (geneparams[j].λ_sq[1] * clusters[k].x_hat[j]) /  (clusters[k].Nk[1] * geneparams[j].λ_sq[1]  + clusters[k].σ_sq_k_hat[j]  )
        end
    end

    return clusters
end
function update_mk_hat_mt!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].mk_hat .= 0.0
        for j in 1:G
            clusters[k].mk_hat[j] +=  (geneparams[j].λ_sq[1] * clusters[k].x_hat[j]) /  (clusters[k].Nk[1] * geneparams[j].λ_sq[1]  + clusters[k].σ_sq_k_hat[j]  )
        end
    end

    return clusters
end
function update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode=nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    if mt_mode == "full"
        Threads.@threads for k in 1:K
            clusters[k].mk_hat .= 0.0
            for j in 1:G
                clusters[k].mk_hat[j] +=  (geneparams[j].λ_sq[1] * clusters[k].x_hat[j]) /  (clusters[k].Nk[1] * geneparams[j].λ_sq[1]  + clusters[k].σ_sq_k_hat[j]  )
            end
        end
    else
        for k in 1:K
            clusters[k].mk_hat .= 0.0
            for j in 1:G
                clusters[k].mk_hat[j] +=  (geneparams[j].λ_sq[1] * clusters[k].x_hat[j]) /  (clusters[k].Nk[1] * geneparams[j].λ_sq[1]  + clusters[k].σ_sq_k_hat[j]  )
            end
        end
    end

    return clusters
end
function update_v_sq_k_hat!_depracated(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]  
    for k in 1:K
        clusters[k].v_sq_k_hat .= 0.0
        clusters[k].v_sq_k_hat .+=   (λ_sq_vec .* clusters[k].σ_sq_k_hat ) ./ (clusters[k].Nk .* λ_sq_vec .+ clusters[k].σ_sq_k_hat  )  
    end

    return clusters
end
function update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].v_sq_k_hat .= 0.0
        for j in 1:G
            clusters[k].v_sq_k_hat[j] +=   (geneparams[j].λ_sq[1] * clusters[k].σ_sq_k_hat[j] ) / (clusters[k].Nk[1] * geneparams[j].λ_sq[1] + clusters[k].σ_sq_k_hat[j])
        end  
    end

    return clusters
end
function update_v_sq_k_hat_mt!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].v_sq_k_hat .= 0.0
        for j in 1:G
            clusters[k].v_sq_k_hat[j] +=   (geneparams[j].λ_sq[1] * clusters[k].σ_sq_k_hat[j] ) / (clusters[k].Nk[1] * geneparams[j].λ_sq[1] + clusters[k].σ_sq_k_hat[j])
        end  
    end
    return clusters
end
function update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode=nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    if mt_mode =="full"
        Threads.@threads for k in 1:K
            clusters[k].v_sq_k_hat .= 0.0
            for j in 1:G
                clusters[k].v_sq_k_hat[j] +=   (geneparams[j].λ_sq[1] * clusters[k].σ_sq_k_hat[j] ) / (clusters[k].Nk[1] * geneparams[j].λ_sq[1] + clusters[k].σ_sq_k_hat[j])
            end  
        end
    else
        for k in 1:K
            clusters[k].v_sq_k_hat .= 0.0
            for j in 1:G
                clusters[k].v_sq_k_hat[j] +=   (geneparams[j].λ_sq[1] * clusters[k].σ_sq_k_hat[j] ) / (clusters[k].Nk[1] * geneparams[j].λ_sq[1] + clusters[k].σ_sq_k_hat[j])
            end  
        end
    end
    return clusters
end
function update_var_muk_hat!_depracated(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].var_muk .= 0.0

        clusters[k].var_muk .+=  clusters[k].yjk_hat .* (clusters[k].mk_hat .^2  + clusters[k].v_sq_k_hat)   .-  clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2
    end

    return clusters
end
function update_var_muk_hat!(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].var_muk .= 0.0
        for j in 1:G
            clusters[k].var_muk[j] +=  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j] ^2  + clusters[k].v_sq_k_hat[j])   -  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end

    end

    return clusters
end
function update_var_muk_hat_mt!(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].var_muk .= 0.0
        for j in 1:G
            clusters[k].var_muk[j] +=  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j] ^2  + clusters[k].v_sq_k_hat[j])   -  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end

    end
    return clusters
end
function update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    if mt_mode == "full"
        Threads.@threads for k in 1:K
            clusters[k].var_muk .= 0.0
            for j in 1:G
                clusters[k].var_muk[j] +=  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j] ^2  + clusters[k].v_sq_k_hat[j])   -  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
            end
        end
    else
        for k in 1:K
            clusters[k].var_muk .= 0.0
            for j in 1:G
                clusters[k].var_muk[j] +=  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j] ^2  + clusters[k].v_sq_k_hat[j])   -  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
            end
    
        end
    end
    return clusters
end
function update_κk_hat!(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].κk_hat .= 0.0

        clusters[k].κk_hat .+=  clusters[k].yjk_hat .* clusters[k].mk_hat
    end

    return clusters
end

function update_κk_hat_mt!(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].κk_hat .= 0.0

        clusters[k].κk_hat .+=  clusters[k].yjk_hat .* clusters[k].mk_hat
    end

    return clusters
end
function update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K


    if mt_mode == "full"
        Threads.@threads for k in 1:K
            clusters[k].κk_hat .= 0.0
    
            clusters[k].κk_hat .+=  clusters[k].yjk_hat .* clusters[k].mk_hat
        end
    else
        for k in 1:K
            clusters[k].κk_hat .= 0.0
    
            clusters[k].κk_hat .+=  clusters[k].yjk_hat .* clusters[k].mk_hat
        end    
    end

    return clusters
end
function update_σ_sq_k_hat!(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].σ_sq_k_hat .= 0.0

        clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 10.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ 10.0)
        # if clusters[k].Nk[1] != 0.0
        #     clusters[k].σ_sq_k_hat .= 0.0

        #     clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 1.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) ) 
        # else
        #     clusters[k].σ_sq_k_hat .= 0.0

        #     clusters[k].σ_sq_k_hat .+=   1 ./ (1e-10 .+ 1.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) ) 
        # end
    end

    return clusters
end
function update_σ_sq_k_hat_mt!(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].σ_sq_k_hat .= 0.0

        clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 10.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ 10.0)
        # if clusters[k].Nk[1] != 0.0
        #     clusters[k].σ_sq_k_hat .= 0.0

        #     clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 1.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) ) 
        # else
        #     clusters[k].σ_sq_k_hat .= 0.0

        #     clusters[k].σ_sq_k_hat .+=   1 ./ (1e-10 .+ 1.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) ) 
        # end
    end

    return clusters
end
function update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    if mt_mode == "full"
        Threads.@threads for k in 1:K
            clusters[k].σ_sq_k_hat .= 0.0
    
            clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 10.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ 10.0)
        end
    else
        for k in 1:K
            clusters[k].σ_sq_k_hat .= 0.0
    
            clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 10.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ 10.0)
        end
    end

    return clusters
end
#### 
####
function findmaxindicies(vvec;boolvec=nothing)
    K = length(vvec)
    vmax = maximum(vvec)
    if isnothing(boolvec)
        boolvec = Vector{Bool}(undef,K)
        boolvec .= false
    end

    for k in 1:K
        if vvec[k] == vmax
            boolvec[k] = true
        else
            boolvec[k] = false
        end
    end
    return boolvec
end
function allocate_to_clusters(x,argmax_r)
    
end
function update_σ_sq_k_hat!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    clusters_cells = [[] for k in 1:K]
    for i in 1:N
        argmax_r = findmaxindicies(cellpop[i].rtik)
        for k in 1:K
            if argmax_r[k]
                push!(clusters_cells[k], collect(cellpop[i].x))
            end
        end
    end
    for k in 1:K
        if isempty(clusters_cells[k])
            xbar = zeros(G)
        else
            xbar =  vec(mean(hcat(clusters_cells[k]...),dims=2))
        end
        for j in 1:G
            if clusters[k].Nk[1] != 0.0
                clusters[k].σ_sq_k_hat[j] = 0.0
                # sum_sq = sum([cellpop[i].rtik[k] * (cellpop[i].x[j] -  clusters[k].κk_hat[j]) ^2  for i in 1:N]) 
                # clusters[k].σ_sq_k_hat[j] +=   1 /clusters[k].Nk[1] * ( sum_sq  +  clusters[k].Nk[1]  * clusters[k].yjk_hat[j] * clusters[k].v_sq_k_hat[j] ) 
                clusters[k].σ_sq_k_hat[j] +=   1 /clusters[k].Nk[1] * ( (xbar[j] -   clusters[k].κk_hat[j]) ^2  +  clusters[k].Nk[1]  * clusters[k].yjk_hat[j] * clusters[k].v_sq_k_hat[j] )
            else
                clusters[k].σ_sq_k_hat[j] = 0.0
                # sum_sq = sum([cellpop[i].rtik[k] * (cellpop[i].x[j] -  clusters[k].κk_hat[j]) ^2  for i in 1:N])
                # clusters[k].σ_sq_k_hat[j] +=   1 ./1e-10 * ( sum_sq  +  clusters[k].Nk[1]  * clusters[k].yjk_hat[j] * clusters[k].v_sq_k_hat[j] ) 
                clusters[k].σ_sq_k_hat[j] +=   1 /1e-10  * ( (xbar[j] -   clusters[k].κk_hat[j]) ^2  +  clusters[k].Nk[1]  * clusters[k].yjk_hat[j] * clusters[k].v_sq_k_hat[j] )
            end
            # if clusters[k].Nk[1] != 0.0
            #     clusters[k].σ_sq_k_hat .= 0.0
            #     sum_sq = sum([cellpop[i].rtik[k] .* (cellpop[i].x .-  clusters[k].κk_hat) .^2  for i in 1:N]) 
            #     clusters[k].σ_sq_k_hat .+=   1 ./clusters[k].Nk .* ( sum_sq  .+  clusters[k].Nk  .* clusters[k].yjk_hat .* clusters[k].v_sq_k_hat ) 
            # else
            #     clusters[k].σ_sq_k_hat .= 0.0
            #     sum_sq = sum([cellpop[i].rtik[k] .* (cellpop[i].x .-  clusters[k].κk_hat) .^2  for i in 1:N])
            #     clusters[k].σ_sq_k_hat .+=   1 ./1e-10 .* ( sum_sq  .+  clusters[k].Nk  .* clusters[k].yjk_hat .* clusters[k].v_sq_k_hat ) 
            # end
        end
    end

    return clusters
end

function update_σ_sq_k_hat!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    clusters_cells = [[] for k in 1:K]
    for i in 1:N
        argmax_r = findmaxindicies(cellpop[i].rtik)
        for k in 1:K
            if argmax_r[k]
                push!(clusters_cells[k], collect(cellpop[i].x))
            end
        end
    end
    for k in 1:K
        if isempty(clusters_cells[k])
            xbar = zeros(G)
        else
            xbar =  vec(mean(hcat(clusters_cells[k]...),dims=2))
        end
        for j in 1:G
            if clusters[k].Nk[1] != 0.0
                clusters[k].σ_sq_k_hat[j] = 0.0
                clusters[k].σ_sq_k_hat[j] +=   1 /clusters[k].Nk[1] * ( (xbar[j] -   clusters[k].κk_hat[j]) ^2  +  clusters[k].Nk[1]  * clusters[k].yjk_hat[j] * clusters[k].v_sq_k_hat[j] )
            else
                clusters[k].σ_sq_k_hat[j] = 0.0
                clusters[k].σ_sq_k_hat[j] +=   1 /1e-10  * ( (xbar[j] -   clusters[k].κk_hat[j]) ^2  +  clusters[k].Nk[1]  * clusters[k].yjk_hat[j] * clusters[k].v_sq_k_hat[j] )
            end

        end
    end

    return clusters
end


function update_λ_sq_k_hat!_depracated(geneparams,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    K = modelparams.K
    for j in 1:G
        geneparams[j].cache[1] = 0.0
        yjk_sum = sum([clusters[k].yjk_hat[j] for k in 1:K])
        for k in 1:K
            geneparams[j].cache[1] += 10. + clusters[k].var_muk[j] + clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end
        geneparams[j].λ_sq[1] = geneparams[j].cache[1] ./ (10. + yjk_sum)
    end

    return geneparams
end
function update_λ_sq_k_hat!(geneparams,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    K = modelparams.K
    for j in 1:G
        geneparams[j].cache[1] = 0.0
        yjk_sum = 0.0 
        for k in 1:K
            yjk_sum += clusters[k].yjk_hat[j]
        end
        for k in 1:K
            geneparams[j].cache[1] += 10. + clusters[k].var_muk[j] + clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end
        geneparams[j].λ_sq[1] = geneparams[j].cache[1] ./ (10. + yjk_sum)
    end

    return geneparams
end
function update_λ_sq_k_hat_mt!(geneparams,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    K = modelparams.K
    Threads.@threads for j in 1:G
        geneparams[j].cache[1] = 0.0
        yjk_sum = 0.0 
        for k in 1:K
            yjk_sum += clusters[k].yjk_hat[j]
        end
        for k in 1:K
            geneparams[j].cache[1] += 10. + clusters[k].var_muk[j] + clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end
        geneparams[j].λ_sq[1] = geneparams[j].cache[1] ./ (10. + yjk_sum)
    end
    return geneparams
end
function update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    G = dataparams.G
    K = modelparams.K

    if mt_mode == "full"
        Threads.@threads for j in 1:G
            geneparams[j].cache[1] = 0.0
            yjk_sum = 0.0 
            for k in 1:K
                yjk_sum += clusters[k].yjk_hat[j]
            end
            for k in 1:K
                geneparams[j].cache[1] += 10. + clusters[k].var_muk[j] + clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
            end
            geneparams[j].λ_sq[1] = geneparams[j].cache[1] ./ (10. + yjk_sum)
        end
    else
        for j in 1:G
            geneparams[j].cache[1] = 0.0
            yjk_sum = 0.0 
            for k in 1:K
                yjk_sum += clusters[k].yjk_hat[j]
            end
            for k in 1:K
                geneparams[j].cache[1] += 10. + clusters[k].var_muk[j] + clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
            end
            geneparams[j].λ_sq[1] = geneparams[j].cache[1] ./ (10. + yjk_sum)
        end
    end
    return geneparams
end

# mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)  
function update_mk_hat_usingXhat_fast3!(clusters,dataparams,modelparams)
    # K = length(Nk)
    # # Nk_xbar_k = Nk .* xbar_k
    # λ0_μ0 =  λ0_vec .* μ0_vec
    # denom = [λ0_vec .+ Nk[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    # mkj_hat = [ (λ0_μ0 .+ x_hat_k[k]) ./denom[k] for k in 1:K]
    # return mkj_hat
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].mk_hat.= 0.0

        clusters[k].mk_hat .+=  (modelparams.λ0_vec .* modelparams.μ0_vec .+ clusters[k].x_hat) ./  (clusters[k].Nkj .+ modelparams.λ0_vec  )
    end

    return clusters
end
function update_mk_hat_usingXhat_fast3_mt!(clusters,dataparams,modelparams)
    # K = length(Nk)
    # # Nk_xbar_k = Nk .* xbar_k
    # λ0_μ0 =  λ0_vec .* μ0_vec
    # denom = [λ0_vec .+ Nk[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    # mkj_hat = [ (λ0_μ0 .+ x_hat_k[k]) ./denom[k] for k in 1:K]
    # return mkj_hat
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].mk_hat.= 0.0

        clusters[k].mk_hat .+=  (modelparams.λ0_vec .* modelparams.μ0_vec .+ clusters[k].x_hat) ./  (clusters[k].Nkj .+ modelparams.λ0_vec  )
    end

    return clusters
end
# b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k) 
function update_b0k_hat_usingXhat_fast3!(clusters,dataparams,modelparams)
    # K = length(Nk)
    # denom = [λ0_vec .+ Nk[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    # μ0_sq_vec = μ0_vec .^2
    # μ0λ0_vec =  λ0_vec .* μ0_vec
    # μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    # numer = [(x_hat_k[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    # ssd = [numer[k] ./ denom[k] for k in 1:K]
    # half_sk_ssd =  1/2 .* [x_hat_sq_k[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    # b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].b0k_hat.= 0.0

        clusters[k].b0k_hat .+=  modelparams.b0_vec .+ 1/2 .* (clusters[k].x_hat_sq .+  (modelparams.μ0_vec .^2 .* modelparams.λ0_vec ) .- ( (clusters[k].x_hat .-  modelparams.λ0_vec .* modelparams.μ0_vec ) .^ 2 ./  (clusters[k].Nkj .+ modelparams.λ0_vec  )))
    end

    return clusters
end
function update_b0k_hat_usingXhat_fast3_mt!(clusters,dataparams,modelparams)
    # K = length(Nk)
    # denom = [λ0_vec .+ Nk[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    # μ0_sq_vec = μ0_vec .^2
    # μ0λ0_vec =  λ0_vec .* μ0_vec
    # μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    # numer = [(x_hat_k[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    # ssd = [numer[k] ./ denom[k] for k in 1:K]
    # half_sk_ssd =  1/2 .* [x_hat_sq_k[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    # b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        # clusters[k]._reset!(clusters[k].λ0k_hat, float_type)
        clusters[k].b0k_hat.= 0.0

        clusters[k].b0k_hat .+=  modelparams.b0_vec .+ 1/2 .* (clusters[k].x_hat_sq .+  (modelparams.μ0_vec .^2 .* modelparams.λ0_vec ) .- ( (clusters[k].x_hat .-  modelparams.λ0_vec .* modelparams.μ0_vec ) .^ 2 ./  (clusters[k].Nkj .+ modelparams.λ0_vec  )))
    end

    return clusters
end

# a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec) 
function update_αt_fast3!(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # K = length(rhok_hat_vec)
    # T = length(θ_hat)
    # a_α_hat = a_α + K
    # a_α_hat_vec = Vector{Float64}(undef,T)
    # b_α_hat_vec = Vector{Float64}(undef,T)
    # e_log_π = log_π_expected_value(θ_hat)
    # e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    @inbounds for t in 1:T
        # a_α_hat_vec[t] = a_α_hat
        # e_βk_log_π = [e_βk_vec[k] .* e_log_π[t][k]  for k in 1:K+1]
        # b_α_hat_vec[t] = b_α - sum(e_βk_log_π)
        conditionparams[t].a_αt_hat[1] = modelparams.a_α + K 
        e_βk_log_π = 0.0
        @inbounds @fastmath for k in 1:Kplus
            e_βk = expectation_βk(k,clusters,modelparams)
            e_log_π = expectation_log_π_tk(conditionparams[t].θ_hat_t[k],conditionparams[t].θ_hat_t_sum[1])
            e_βk_log_π += e_βk * e_log_π
        end
        conditionparams[t].b_αt_hat[1] = modelparams.b_α - e_βk_log_π
    end
    return conditionparams
end

function update_αt_fast3_mt!(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # K = length(rhok_hat_vec)
    # T = length(θ_hat)
    # a_α_hat = a_α + K
    # a_α_hat_vec = Vector{Float64}(undef,T)
    # b_α_hat_vec = Vector{Float64}(undef,T)
    # e_log_π = log_π_expected_value(θ_hat)
    # e_βk_vec = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    Threads.@threads for t in 1:T
        # a_α_hat_vec[t] = a_α_hat
        # e_βk_log_π = [e_βk_vec[k] .* e_log_π[t][k]  for k in 1:K+1]
        # b_α_hat_vec[t] = b_α - sum(e_βk_log_π)
        conditionparams[t].a_αt_hat[1] = modelparams.a_α + K 
        e_βk_log_π = 0.0
        @inbounds @fastmath for k in 1:Kplus
            e_βk = expectation_βk(k,clusters,modelparams)
            e_log_π = expectation_log_π_tk(conditionparams[t].θ_hat_t[k],conditionparams[t].θ_hat_t_sum[1])
            e_βk_log_π += e_βk * e_log_π
        end
        conditionparams[t].b_αt_hat[1] = modelparams.b_α - e_βk_log_π
    end
    return conditionparams
end
# awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
# bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec) 
function update_abwt_hat_fast3!(conditionparams,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # T = length(c_ttprime)
    # awt_hat_vec = Vector{Float64}(undef,T-1)
    # bwt_hat_vec = Vector{Float64}(undef,T-1)
    for t in 2:T
        awt_hat = 0.0
        awt_hat += modelparams.adot_w
        bwt_hat = 0.0
        bwt_hat += modelparams.bdot_w
        # sum_string = ""
        # sum_string *= "adot_w "
        for t_prime in t:T
            # c_string = "+ c$(t_prime)$(t) "
            awt_hat += conditionparams[t_prime].c_tt_prime[t]
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                bwt_hat += conditionparams[t_prime].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        conditionparams[t-1].awt_hat[1]=awt_hat
        conditionparams[t-1].bwt_hat[1]=bwt_hat
        # println(sum_string)
    end
    conditionparams[T].awt_hat[1]=0.0
    conditionparams[T].bwt_hat[1]=0.0
    return conditionparams
end
function update_abwt_hat_fast3_mt!(conditionparams,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # T = length(c_ttprime)
    # awt_hat_vec = Vector{Float64}(undef,T-1)
    # bwt_hat_vec = Vector{Float64}(undef,T-1)
    Threads.@threads for t in 2:T
        awt_hat = 0.0
        awt_hat += modelparams.adot_w
        bwt_hat = 0.0
        bwt_hat += modelparams.bdot_w
        # sum_string = ""
        # sum_string *= "adot_w "
        for t_prime in t:T
            # c_string = "+ c$(t_prime)$(t) "
            awt_hat += conditionparams[t_prime].c_tt_prime[t]
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                bwt_hat += conditionparams[t_prime].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        conditionparams[t-1].awt_hat[1]=awt_hat
        conditionparams[t-1].bwt_hat[1]=bwt_hat
        # println(sum_string)
    end
    conditionparams[T].awt_hat[1]=0.0
    conditionparams[T].bwt_hat[1]=0.0
    return conditionparams
end
function update_st_hat_fast3!(conditionparams,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # T = length(c_ttprime)
    # awt_hat_vec = Vector{Float64}(undef,T-1)
    # bwt_hat_vec = Vector{Float64}(undef,T-1)
    for t in 2:T
        st_hat = 0.0
        st_hat += modelparams.ϕ0
        # sum_string = ""
        # sum_string *= "adot_w "
        for t_prime in t:T
            # c_string = "+ c$(t_prime)$(t) "
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                st_hat += conditionparams[t_prime].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        conditionparams[t-1].st_hat[1]=st_hat
        # println(sum_string)
    end
    conditionparams[T].st_hat[1]=0.0
    return conditionparams
end
function update_st_hat_fast3_mt!(conditionparams,dataparams,modelparams) 
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # T = length(c_ttprime)
    # awt_hat_vec = Vector{Float64}(undef,T-1)
    # bwt_hat_vec = Vector{Float64}(undef,T-1)
    Threads.@threads for t in 2:T
        st_hat = 0.0
        st_hat += modelparams.ϕ0
        # sum_string = ""
        # sum_string *= "adot_w "
        for t_prime in t:T
            # c_string = "+ c$(t_prime)$(t) "
            for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                st_hat += conditionparams[t_prime].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        conditionparams[t-1].st_hat[1]=st_hat
        # println(sum_string)
    end
    conditionparams[T].st_hat[1]=0.0
    return conditionparams
end


function update_st_hat_mpu!(conditionparams,dataparams,modelparams;mt_mode = nothing) 
    float_type = dataparams.BitType
    T = dataparams.T

    if mt_mode == "full"
        Threads.@threads for t in 2:T
            st_hat = 0.0
            st_hat += modelparams.ϕ0
            for t_prime in t:T
                for l in 1:t-1
                    st_hat += conditionparams[t_prime].c_tt_prime[l]
                end
            end
            conditionparams[t-1].st_hat[1]=st_hat
        end
        conditionparams[T].st_hat[1]=0.0
    else
        for t in 2:T
            st_hat = 0.0
            st_hat += modelparams.ϕ0
            for t_prime in t:T
                for l in 1:t-1
                    st_hat += conditionparams[t_prime].c_tt_prime[l]
                end
            end
            conditionparams[t-1].st_hat[1]=st_hat
        end
        conditionparams[T].st_hat[1]=0.0
    end
    return conditionparams
end


# a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec) 
function update_γ_fast3(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    K = modelparams.K
    # K = length(rhok_hat_vec)
    e_log_minus_uk_sum = 0.0 #expectation_log1minusUk(rhok_hat_vec, omegak_hat_vec)
    a_γ_hat = modelparams.a_γ + K
    for k in 1:K
        e_log_minus_uk_sum += expectation_log1minusUk(clusters[k].rhok_hat[1],clusters[k].omegak_hat[1])
    end
    b_γ_hat = modelparams.b_γ - e_log_minus_uk_sum
    return a_γ_hat, b_γ_hat
end
function update_γ_fast3_mt(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    K = modelparams.K
    # K = length(rhok_hat_vec)
    e_log_minus_uk_sum = 0.0 #expectation_log1minusUk(rhok_hat_vec, omegak_hat_vec)
    a_γ_hat = modelparams.a_γ + K
    Threads.@threads for k in 1:K
        e_log_minus_uk_sum += expectation_log1minusUk(clusters[k].rhok_hat[1],clusters[k].omegak_hat[1])
    end
    b_γ_hat = modelparams.b_γ - e_log_minus_uk_sum
    return a_γ_hat, b_γ_hat
end

function update_Tk(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # e_log_π = log_π_expected_value(θ_hat)
    # e_αt_vec = αt_expected_value.(a_α_hat_vec,b_α_hat_vec)
    # e_αt_log_π =  e_αt_vec .*  e_log_π
    # Tαk = sum(e_αt_log_π)
    Tk = zeros(Float64,Kplus)
    @inbounds for k in 1:Kplus
        @inbounds @fastmath for t in 1:T
            e_log_π = expectation_log_π_tk(conditionparams[t].d_hat_t[k],conditionparams[t].d_hat_t_sum[1])
            Tk[k] +=  e_log_π
        end
    end
    return Tk
end
function update_Tk_mt(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # e_log_π = log_π_expected_value(θ_hat)
    # e_αt_vec = αt_expected_value.(a_α_hat_vec,b_α_hat_vec)
    # e_αt_log_π =  e_αt_vec .*  e_log_π
    # Tαk = sum(e_αt_log_π)
    Tk = zeros(Float64,Kplus)
    Threads.@threads  for k in 1:Kplus
        @inbounds @fastmath for t in 1:T
            e_log_π = expectation_log_π_tk(conditionparams[t].d_hat_t[k],conditionparams[t].d_hat_t_sum[1])
            Tk[k] +=  e_log_π
        end
    end
    return Tk
end
function update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams;mt_mode = nothing)
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    Tk .= 0.0
    if mt_mode == "full"
        Threads.@threads  for k in 1:Kplus
            @inbounds @fastmath for t in 1:T
                e_log_π = expectation_log_π_tk(conditionparams[t].d_hat_t[k],conditionparams[t].d_hat_t_sum[1])
                Tk[k] +=  e_log_π
            end
        end
    else
        @inbounds for k in 1:Kplus
            @inbounds @fastmath for t in 1:T
                e_log_π = expectation_log_π_tk(conditionparams[t].d_hat_t[k],conditionparams[t].d_hat_t_sum[1])
                Tk[k] +=  e_log_π
            end
        end
    end
end
# e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
# Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)  Threads.@threads
function update_Tαk_fast3(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # e_log_π = log_π_expected_value(θ_hat)
    # e_αt_vec = αt_expected_value.(a_α_hat_vec,b_α_hat_vec)
    # e_αt_log_π =  e_αt_vec .*  e_log_π
    # Tαk = sum(e_αt_log_π)
    Tαk = zeros(Float64,Kplus)
    @inbounds for k in 1:Kplus
        @inbounds @fastmath for t in 1:T
            e_αt = expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
            # e_βk = expectation_βk(k,clusters,modelparams)
            e_log_π = expectation_log_π_tk(conditionparams[t].θ_hat_t[k],conditionparams[t].θ_hat_t_sum[1])
            Tαk[k] += e_αt * e_log_π
        end
    end
    return Tαk
end
function update_Tαk_fast3_mt(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    # e_log_π = log_π_expected_value(θ_hat)
    # e_αt_vec = αt_expected_value.(a_α_hat_vec,b_α_hat_vec)
    # e_αt_log_π =  e_αt_vec .*  e_log_π
    # Tαk = sum(e_αt_log_π)
    Tαk = zeros(Float64,Kplus)
    Threads.@threads for k in 1:Kplus
        @inbounds @fastmath for t in 1:T
            e_αt = expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
            # e_βk = expectation_βk(k,clusters,modelparams)
            e_log_π = expectation_log_π_tk(conditionparams[t].θ_hat_t[k],conditionparams[t].θ_hat_t_sum[1])
            Tαk[k] += e_αt * e_log_π
        end
    end
    return Tαk
end
# rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
function update_rho_omega_hat_fast3!(clusters,dataparams,modelparams,e_γ,Tαk;optim_max_iter=100000)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    K = modelparams.K
    # rhok_hat_vec = [clusters[k].rhok_hat[1] for k in 1:K ]
    # omegak_hat_vec = [clusters[k].omegak_hat[1] for k in 1:K]
    c_hat = [clusters[k].c_hat[1] for k in 1:K ]
    d_hat =  [clusters[k].d_hat[1] for k in 1:K ]
    cd_vec_args = [c_hat , d_hat]
    cd_vec_args0  = permutedims(reduce(hcat,cd_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,e_γ,Tαk)
    gg_uncon! = g_unconstrained_closure!(T,e_γ,Tαk)
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(), Optim.Options(iterations = optim_max_iter))
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter)) #,linesearch=Optim.LineSearches.BackTracking() 
    # optimize(LB_LG_unconstrained,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))

    # @debug Optim.converged(lb_lg_results)
    # new_rho_hat = sigmoid.(lb_lg_results.res.minimizer[1,:])
    # new_omega_hat = exp.(lb_lg_results.res.minimizer[2,:])

    # new_c_hat = StatsFuns.logit.(new_rho_hat)
    # new_d_hat = log.(new_omega_hat)
    for k in 1:K
        clusters[k].rhok_hat[1] = sigmoid(lb_lg_results.res.minimizer[1,k])#new_rho_hat[k]
        clusters[k].omegak_hat[1] = exp(lb_lg_results.res.minimizer[2,k])#new_omega_hat[k]
        clusters[k].c_hat[1] = StatsFuns.logit(clusters[k].rhok_hat[1]) #new_c_hat[k]
        clusters[k].d_hat[1] = log(clusters[k].omegak_hat[1])#new_d_hat[k]
    end
    return clusters
end
function update_gh_hat!(clusters,dataparams,modelparams,Tk;optim_max_iter=100000)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0
    # rhok_hat_vec = [clusters[k].rhok_hat[1] for k in 1:K ]
    # omegak_hat_vec = [clusters[k].omegak_hat[1] for k in 1:K]
    a_hat = [clusters[k].ak_hat[1] for k in 1:K ]
    b_hat =  [clusters[k].bk_hat[1] for k in 1:K ]
    ab_vec_args = [a_hat , b_hat]
    ab_vec_args0  = permutedims(reduce(hcat,ab_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,γ0,α0,Tk)
    gg_uncon! = g_unconstrained_closure!(T,γ0,α0,Tk)
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(), Optim.Options(iterations = optim_max_iter))
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,ab_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter)) #,linesearch=Optim.LineSearches.BackTracking() 
    # optimize(LB_LG_unconstrained,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))

    # @debug Optim.converged(lb_lg_results)
    # new_rho_hat = sigmoid.(lb_lg_results.res.minimizer[1,:])
    # new_omega_hat = exp.(lb_lg_results.res.minimizer[2,:])

    # new_c_hat = StatsFuns.logit.(new_rho_hat)
    # new_d_hat = log.(new_omega_hat)
    for k in 1:K
        clusters[k].gk_hat[1] = sigmoid(lb_lg_results.res.minimizer[1,k])#new_rho_hat[k]
        clusters[k].hk_hat[1] = exp(lb_lg_results.res.minimizer[2,k])#new_omega_hat[k]
        clusters[k].ak_hat[1] = StatsFuns.logit(clusters[k].gk_hat[1]) #new_c_hat[k]
        clusters[k].bk_hat[1] = log(clusters[k].hk_hat[1])#new_d_hat[k]
    end
    return clusters
end
function update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=100000)
    float_type = dataparams.BitType
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0
    a_hat = [clusters[k].ak_hat[1] for k in 1:K ]
    b_hat =  [clusters[k].bk_hat[1] for k in 1:K ]
    ab_vec_args = [a_hat , b_hat]
    ab_vec_args0  = permutedims(reduce(hcat,ab_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,γ0,α0,Tk)
    gg_uncon! = g_unconstrained_closure!(T,γ0,α0,Tk)
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(), Optim.Options(iterations = optim_max_iter))
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,ab_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter)) #,linesearch=Optim.LineSearches.BackTracking() 
    # optimize(LB_LG_unconstrained,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    for k in 1:K
        clusters[k].gk_hat[1] = sigmoid(lb_lg_results.res.minimizer[1,k])#new_rho_hat[k]
        clusters[k].hk_hat[1] = exp(lb_lg_results.res.minimizer[2,k])#new_omega_hat[k]
        clusters[k].ak_hat[1] = StatsFuns.logit(clusters[k].gk_hat[1]) #new_c_hat[k]
        clusters[k].bk_hat[1] = log(clusters[k].hk_hat[1])#new_d_hat[k]
    end
    return clusters
end

#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################
function tidy_update_rtik(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat; rtik_mat_init = nothing )
    T = length(unique(e_log_π_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
    Glog = G*log(2π)
    Kplus = K + 1
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))
    N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
    if isnothing(rtik_mat_init)
        timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))
        N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
        nrows = N*K
        ncols = 4
        rtik_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        # println("The start")
        # println("K: $K")
        # println("N_t: $N_t")
        # println("timepoints: $timepoints")
        time = innermelt(timepoints, K .* N_t)
        # println("The mid")
        cells = innermelt(cell_ids,K)
        # println("The End")
        states = outermelt(state_ids,N)
        rtik_mat[:,1] = time
        rtik_mat[:,2] = cells
        rtik_mat[:,3] = states
    else
        rtik_mat = rtik_mat_init
    end
    # println("We outside")
    e_τ_μ_counter = 0
    time_indx_starts = collect(0:T-1).*  N_t .+ 1
    new_t = 1
    t = 0
    for i in cell_ids
        ptik_tilde = Vector{Float64}(undef,K)
        if i <= time_indx_starts[new_t]
            t += 1
            new_t += 1
            if new_t > T
                new_t = T
            end
            if t > T
                t = T
            end
        end
        time_index = (t-1)*T + 1:(t)*T
        c_tl = c_ttprime_mat[time_index,3]
        states_index(n) = (n-1)*Kplus + 1:(n)*Kplus
        # e_log_π_tk = e_log_π[states_index, 3] states_index(tt)
        adjusted_e_log_π_tk = sum([c_tl[tt] .* e_log_π_mat[states_index(tt), 3] for tt in 1:t]) # (tt-1)*Kplus + 1:(tt)*Kplus
        for k in state_ids
            e_τ_μ_counter +=1
            τ_kj_index = (k-1)*G + 1:(k)*G
            τ_μ_index = (e_τ_μ_counter-1)*G + 1:(e_τ_μ_counter)*G
            e_log_τ = sum(e_log_τ_kj_mat[τ_kj_index,end])
            e_log_τ_μ = sum(e_τ_μ_tikj_mat[τ_μ_index,end])
            ptik_tilde[k] = adjusted_e_log_π_tk[k] .+ 0.5 .* (e_log_τ .- Glog .- e_log_τ_μ) # 
        end
        rtik = norm_weights(ptik_tilde)
        cell_index = (i-1)*K + 1:(i)*K
        rtik_mat[cell_index,end] = rtik
    end

    return rtik_mat
end
function tidy_update_rtik_SparseVS(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat,ηtkj_mat; rtik_mat_init = nothing )
    T = length(unique(e_log_π_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
 
    Glog = G*log(2π)
    logpi= Glog/G
    Kplus = K + 1
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))

    N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
    
    if isnothing(rtik_mat_init)
        timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))
        N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
        nrows = N*K
        ncols = 4
        rtik_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* N_t)
        cells = innermelt(cell_ids,K)
        states = outermelt(state_ids,N)
        rtik_mat[:,1] = time
        rtik_mat[:,2] = cells
        rtik_mat[:,3] = states
    else
        rtik_mat = rtik_mat_init
    end

    e_τ_μ_counter = 0
    time_indx_starts = collect(0:T-1).*  N_t .+ 1
    new_t = 1
    t = 0

    for i in cell_ids
        ptik_tilde = Vector{Float64}(undef,K)
        if i <= time_indx_starts[new_t]
            t += 1
            new_t += 1
            if new_t > T
                new_t = T
            end
            if t > T
                t = T
            end
        end
        time_index = (t-1)*T + 1:(t)*T
        @views c_tl = c_ttprime_mat[time_index,3]
        states_index(n) = (n-1)*Kplus + 1:(n)*Kplus
        # e_log_π_tk = e_log_π[states_index, 3] states_index(tt)
        adjusted_e_log_π_tk = sum([c_tl[tt] .* e_log_π_mat[states_index(tt), 3] for tt in 1:t]) # (tt-1)*Kplus + 1:(tt)*Kplus
        η_time_index = (t-1)*K*G*2 + 1:(t)*K*G*2
        @views ηkj_mat = ηtkj_mat[η_time_index,:]
        for k in state_ids
            e_τ_μ_counter +=1
            e_τ_μ_err_index = (i-1)*G + 1:(i)*G
            τ_kj_index = (k-1)*G + 1:(k)*G
            η_state_index = (k-1)*2*G+1:(k)*2*G
            τ_μ_index = (e_τ_μ_counter-1)*G + 1:(e_τ_μ_counter)*G
            @views e_log_τ_kj = e_log_τ_kj_mat[τ_kj_index,end]
            @views e_log_τ_μ_tikj = e_τ_μ_tikj_mat[τ_μ_index,end]
            @views ηj_mat = ηkj_mat[η_state_index,:]
            @views η_true = ηj_mat[1:2:end,5]
            ptik_tilde[k] = adjusted_e_log_π_tk[k] .+ sum(0.5 .* η_true .* (e_log_τ_kj .- logpi .- e_log_τ_μ_tikj))# 
        end
        rtik = norm_weights(ptik_tilde)
        cell_index = (i-1)*K + 1:(i)*K
        @views rtik_mat[cell_index,end] = rtik
    end

    return rtik_mat
end
function tidy_update_rtik_VS1(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat,e_log_τj_err_mat,e_τ_μ_tij_err_mat,ηtkj_mat; rtik_mat_init = nothing )
    T = length(unique(e_log_π_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
 
    Glog = G*log(2π)
    logpi= Glog/G
    Kplus = K + 1
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))

    N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
    
    if isnothing(rtik_mat_init)
        timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))
        N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
        nrows = N*K
        ncols = 4
        rtik_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* N_t)
        cells = innermelt(cell_ids,K)
        states = outermelt(state_ids,N)
        rtik_mat[:,1] = time
        rtik_mat[:,2] = cells
        rtik_mat[:,3] = states
    else
        rtik_mat = rtik_mat_init
    end

    e_τ_μ_counter = 0
    time_indx_starts = collect(0:T-1).*  N_t .+ 1
    new_t = 1
    t = 0

    for i in cell_ids
        ptik_tilde = Vector{Float64}(undef,K)
        if i <= time_indx_starts[new_t]
            t += 1
            new_t += 1
            if new_t > T
                new_t = T
            end
            if t > T
                t = T
            end
        end
        time_index = (t-1)*T + 1:(t)*T
        @views c_tl = c_ttprime_mat[time_index,3]
        states_index(n) = (n-1)*Kplus + 1:(n)*Kplus
        # e_log_π_tk = e_log_π[states_index, 3] states_index(tt)
        adjusted_e_log_π_tk = sum([c_tl[tt] .* e_log_π_mat[states_index(tt), 3] for tt in 1:t]) # (tt-1)*Kplus + 1:(tt)*Kplus
        η_time_index = (t-1)*K*G*2 + 1:(t)*K*G*2
        @views ηkj_mat = ηtkj_mat[η_time_index,:]
        for k in state_ids
            e_τ_μ_counter +=1
            e_τ_μ_err_index = (i-1)*G + 1:(i)*G
            τ_kj_index = (k-1)*G + 1:(k)*G
            η_state_index = (k-1)*2*G+1:(k)*2*G
            τ_μ_index = (e_τ_μ_counter-1)*G + 1:(e_τ_μ_counter)*G
            @views e_log_τ = e_log_τ_kj_mat[τ_kj_index,end]
            @views e_log_τ_μ = e_τ_μ_tikj_mat[τ_μ_index,end]
            @views e_log_τj_err = e_log_τj_err_mat[:,end]
            @views e_τ_μ_tij_err = e_τ_μ_tij_err_mat[e_τ_μ_err_index,end] 
            @views ηj_mat = ηkj_mat[η_state_index,:]
            @views η_true = ηj_mat[1:2:end,5]
            @views η_false = ηj_mat[2:2:end,5]
            ptik_tilde[k] = adjusted_e_log_π_tk[k] .+ sum(0.5 .* η_true .* (e_log_τ .- logpi .- e_log_τ_μ) + 0.5 .* η_false .* (e_log_τj_err .- logpi .- e_τ_μ_tij_err) )# 
        end
        rtik = norm_weights(ptik_tilde)
        cell_index = (i-1)*K + 1:(i)*K
        @views rtik_mat[cell_index,end] = rtik
    end

    return rtik_mat
end
function tidy_update_rtik_VS2(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat,ηtkj_mat; rtik_mat_init = nothing )
    # e_log_τj_err_mat,e_τ_μ_tij_err_mat,
    T = length(unique(e_log_π_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
 
    Glog = G*log(2π)
    logpi= Glog/G
    Kplus = K + 1
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))

    N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
    
    if isnothing(rtik_mat_init)
        timepoint_freq = countmap(Int.(e_τ_μ_tikj_mat[1:G:end,1]))
        N_t = Int.([timepoint_freq[key] ./K for key in sort(collect(keys(timepoint_freq)))])
        nrows = N*K
        ncols = 4
        rtik_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* N_t)
        cells = innermelt(cell_ids,K)
        states = outermelt(state_ids,N)
        rtik_mat[:,1] = time
        rtik_mat[:,2] = cells
        rtik_mat[:,3] = states
    else
        rtik_mat = rtik_mat_init
    end

    e_τ_μ_counter = 0
    time_indx_starts = collect(0:T-1).*  N_t .+ 1
    new_t = 1
    t = 0

    for i in cell_ids
        ptik_tilde = Vector{Float64}(undef,K)
        if i <= time_indx_starts[new_t]
            t += 1
            new_t += 1
            if new_t > T
                new_t = T
            end
            if t > T
                t = T
            end
        end
        time_index = (t-1)*T + 1:(t)*T
        @views c_tl = c_ttprime_mat[time_index,3]
        states_index(n) = (n-1)*Kplus + 1:(n)*Kplus
        # e_log_π_tk = e_log_π[states_index, 3] states_index(tt)
        adjusted_e_log_π_tk = sum([c_tl[tt] .* e_log_π_mat[states_index(tt), 3] for tt in 1:t]) # (tt-1)*Kplus + 1:(tt)*Kplus
        η_time_index = (t-1)*K*G*2 + 1:(t)*K*G*2
        @views ηkj_mat = ηtkj_mat[η_time_index,:]
        for k in state_ids
            e_τ_μ_counter +=1
            
            τ_kj_index = (k-1)*G + 1:(k)*G
            η_state_index = (k-1)*2*G+1:(k)*2*G
            τ_μ_index = (e_τ_μ_counter-1)*G + 1:(e_τ_μ_counter)*G
            @views e_log_τ = e_log_τ_kj_mat[τ_kj_index,end]
            @views e_log_τ_μ = e_τ_μ_tikj_mat[τ_μ_index,end] 
            @views ηj_mat = ηkj_mat[η_state_index,:]
            @views η_true = ηj_mat[1:2:end,5]
            # e_τ_μ_err_index = (i-1)*G + 1:(i)*G
            # @views e_log_τj_err = e_log_τj_err_mat[:,end]
            # @views e_τ_μ_tij_err = e_τ_μ_tij_err_mat[e_τ_μ_err_index,end]
            # @views η_false = ηj_mat[2:2:end,5]
            ptik_tilde[k] = adjusted_e_log_π_tk[k] .+ sum(0.5 .* η_true .* (e_log_τ .- logpi .- e_log_τ_μ) )# + 0.5 .* η_false .* (e_log_τj_err .- logpi .- e_τ_μ_tij_err)  
        end
        rtik = norm_weights(ptik_tilde)
        cell_index = (i-1)*K + 1:(i)*K
        @views rtik_mat[cell_index,end] = rtik
    end

    return rtik_mat
end
function tidy_update_θ_hat(Ntk_mat,ρkωkckdk_hat_mat,aαtbαtawtbwt_hat_mat,c_ttprime_hat_mat)
    T = length(unique(Ntk_mat[:,1]))
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    Kplus = K+1
    e_βk_mat = tidy_βk_expected_value(ρkωkckdk_hat_mat);
    e_αt_mat =  tidy_αt_expected_value(aαtbαtawtbwt_hat_mat);
    state_ids = collect(1:Kplus)
    timepoints = collect(1:T)
    states = outermelt(state_ids,T)
    time = innermelt(timepoints,Kplus)
    nrows = length(time)
    ncols = 3
    θ_hat_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    θ_hat_mat[:,1] = time
    θ_hat_mat[:,2] = states
    θ_hat = Vector{Vector{Float64}}(undef,T) 
    ntk_indx(t) = (t-1)*K + 1:(t)*K
    Ntk = [push!(deepcopy(Ntk_mat[ntk_indx(t),end]),0.0) for t in 1:T]
    for t in 1:T
        c_indx = (t-1)*T+1:(t)*T # = 
        updated_Ntk = Ntk .*  c_ttprime_hat_mat[c_indx,3]
        sum_updated_Ntk = sum(updated_Ntk)
        e_αt_βk = e_αt_mat[t,2] .* e_βk_mat[:,2]
        θ_hat[t] = e_αt_βk .+ sum_updated_Ntk
    end
    θ_hat_mat[:,3] = recursive_flatten(θ_hat)
    return θ_hat_mat
end
function tidy_update_c_ttprime(aαtbαtawtbwt_hat_mat,rtik_mat,θ_hat_mat; e_log_π_mat_init=nothing)
    T = length(unique(rtik_mat[:,1]))
    Kplus = length(unique(θ_hat_mat[:,2]))
    K = Kplus -1
    e_log_w_ttprime_mat  = tidy_log_w_ttprime_expected_value(aαtbαtawtbwt_hat_mat);
    e_log_π_mat =  tidy_log_π_expected_value(θ_hat_mat; e_log_π_mat_init=e_log_π_mat_init);
    Ntk_mat = tidy_update_Ntk(rtik_mat);
    timepoints = collect(1:T)
    timepoints_prime = collect(1:T)
    time = innermelt(timepoints,T)
    timeprime = outermelt(timepoints_prime,T)
    nrows = length(time)
    ncols = 3
    
    c_ttprime = Vector{Vector{Float64}}(undef,T)
    ntk_indx(t) = (t-1)*K + 1:(t)*K
    Ntk = [push!(deepcopy(Ntk_mat[ntk_indx(t),end]),0.0) for t in 1:T]
    e_log_π_indx(t) = (t-1)*Kplus + 1:(t)*Kplus
    e_log_π = [deepcopy(e_log_π_mat[e_log_π_indx(t),end]) for t in 1:T]
    ts = 0
    for t in 1:T
        w_indx = ts+1:(ts+t)
        ts = ts+t
        adjusted_e_log_π_k = [Ntk[t] .* el for el in e_log_π[1:t]]
        adjusted_e_log_π = sum.(adjusted_e_log_π_k)
        p_c_ttprime = e_log_w_ttprime_mat[w_indx,end] .+ adjusted_e_log_π
        # c_ttprime[t] = norm_weights(p_c_ttprime)
        c_ttprime[t] = vcat(norm_weights(p_c_ttprime),zeros(Float64,T-t))
        # c_full = push!(c_trunc,zeros(Float64,T-t))
        # c_ttprime[t] = c_full
    end
    c_ttprime_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    c_ttprime_mat[:,1] = time
    c_ttprime_mat[:,2] = timeprime
    c_ttprime_mat[:,3] = recursive_flatten(c_ttprime)
    return c_ttprime_mat
end
function tidy_update_Nk(rtik_mat)
    K = length(unique(rtik_mat[:,3]))
    T = length(unique(rtik_mat[:,1]));
    timepoints = collect(1:T);
    state_ids = collect(1:K);
    nrows = K
    ncols = 2
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    Nk_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    Nk_mat[:,1] = state_ids
    Nk_mat[:,2] = [sum([sum(rtik_mat[k:K:end,end][st:en]) for (st,en) in timeranges]) for k in 1:K] #[sum(rtik_mat[k:K:end,end]) for k in 1:K]
    return Nk_mat
end
function tidy_update_Nkj(rtik_mat,ηtkj_mat)
    K = length(unique(rtik_mat[:,3]))
    T = length(unique(rtik_mat[:,1]));
    G = length(unique(ηtkj_mat[:,3]));
    timepoints = collect(1:T);
    state_ids = collect(1:K);
    gene_ids = collect(1:G);
    nrows = K*G
    ncols = 3
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    Nkj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    Nkj_vec = Vector{Vector{Float64}}(undef,K);

    Nkj_mat[:,1] = innermelt(state_ids,G)
    Nkj_mat[:,2] = outermelt(gene_ids,K)
    ηtkj_mat_tr = ηtkj_mat[1:2:end,:];
    ηtkj_vec_tr = [ηtkj_mat_tr[(i-1)*G+1:(i)*G,end] for i in 1:T*K];
    for k in state_ids
        rtik_k =  rtik_mat[k:K:end,end];
        ηtkj_vec_tr_k = [repeat([el],val) for (el,val) in zip(ηtkj_vec_tr[k:K:end],N_t)];
        ηtkj_vec_tr_k = vcat(ηtkj_vec_tr_k...)
        Nkj_ = [r .*h for (r,h) in zip(rtik_k,ηtkj_vec_tr_k)]#rtik .* x_sq_vec
        Nkj_vec[k] = sum([sum(Nkj_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end

    Nkj_mat[:,3] = recursive_flatten(Nkj_vec)
    return Nkj_mat
end
function tidy_update_Nej(rtik_mat,ηtkj_mat)
    K = length(unique(rtik_mat[:,3]))
    T = length(unique(rtik_mat[:,1]));
    G = length(unique(ηtkj_mat[:,3]));
    timepoints = collect(1:T);
    state_ids = collect(1:1);
    gene_ids = collect(1:G);
    nrows = 1*G
    ncols = 3
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    Nej_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    Nej_vec = Vector{Vector{Float64}}(undef,K);

    Nej_mat[:,1] = innermelt(state_ids,G)
    Nej_mat[:,2] = outermelt(gene_ids,1)
    ηtkj_mat_fs = ηtkj_mat[2:2:end,:];
    ηtkj_vec_fs = [ηtkj_mat_fs[(i-1)*G+1:(i)*G,end] for i in 1:T*K];

    for k in 1:K
        rtik_k =  rtik_mat[k:K:end,end];
        ηtkj_vec_fs_k = [repeat([el],val) for (el,val) in zip(ηtkj_vec_fs[k:K:end],N_t)];
        ηtkj_vec_fs_k = vcat(ηtkj_vec_fs_k...)
        Nej_ = [r .*h for (r,h) in zip(rtik_k,ηtkj_vec_fs_k)]#rtik .* x_sq_vec
        Nej_vec[k] = sum([sum(Nej_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end

    Nej_mat[:,3] = sum(Nej_vec)
    return Nej_mat
end
function tidy_update_Ntk(rtik_mat)
    T = length(unique(rtik_mat[:,1]))
    N = length(unique(rtik_mat[:,2]))
    K = length(unique(rtik_mat[:,3]))
    # Ntk =  Vector{Vector{Float64}}(undef,T)
    nrows = T*K
    ncols = 3
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    time = innermelt(timepoints,K)
    states = outermelt(state_ids,T)
    Ntk_mat =  Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    Ntk_mat[:,1]= time
    Ntk_mat[:,2]= states
    # time_cell_starts_indx = collect(0:T-1).*  N_t .*K .+ 1
    # time_cell_ends_indx =  collect(1:T).*  N_t*K
    # for t in 1:T
    #     time_indx = (t-1)*K + 1:(t)*K
    #     ts = time_cell_starts_indx[t]
    #     te = time_cell_ends_indx[t]
    #     Ntk_vec = [sum(rtik_mat[ts+k:K:te,end] ) for k in 0:K-1] # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
    #     Ntk_mat[time_indx,end] = Ntk_vec
    # end
    rtik_vec = [rtik_mat[(i-1)*K+1:(i)*K,end] for i in 1:N]
    timeranges = tidy_get_timeranges(N_t) #zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    t_indx = 0
    for (st,en) in timeranges
        t_indx +=1
        time_indx = (t_indx-1)*K + 1:(t_indx)*K
        Ntk_vec = sum(rtik_vec[st:en] )# K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        Ntk_mat[time_indx,end] = Ntk_vec
    end
    return Ntk_mat
end
function tidy_update_x_hat_k(xmat,rtik_mat)
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
    nrows = K*G; 
    ncols = 3;
    x_hat_k_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_k_vec = Vector{Vector{Float64}}(undef,K);
    x_vec =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    # all([all(a .== b) for (a,b) in zip(x_vec,vcat(x_to_use...)) ])
    # all([all(rtik_mat[k:K:end,end] .== [el[k] for el in vcat(rtik...)]) for k in 1:K])
    
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in states_id
        rtik =  rtik_mat[k:K:end,end];
        x_hat_ = rtik .* x_vec
        x_hat_k_vec[k] = sum([sum(x_hat_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_k_mat[:,1] = states
    x_hat_k_mat[:,2] = genes
    x_hat_k_mat[:,end] = recursive_flatten(x_hat_k_vec)
    return x_hat_k_mat
    # all(sum(x_hat_[1:600]) .== x_hat_tk[1][1])
    # all(sum(x_hat_[601:1200]) .== x_hat_tk[1][2])
    # all(sum(x_hat_[1201:1800]) .== x_hat_tk[1][3])
    # all(sum(x_hat_[2401:3000]) .== x_hat_tk[1][5])
end
function tidy_update_x_hat_sq_k(xmat,rtik_mat)
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
    nrows = K*G; 
    ncols = 3;
    x_hat_sq_k_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_sq_k_vec = Vector{Vector{Float64}}(undef,K);
    x_sq_vec =  [(xmat[(i-1)*G+1:(i)*G,end]) .^ 2 for i in 1:N];

    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in states_id
        rtik =  rtik_mat[k:K:end,end];
        x_hat_sq_ = rtik .* x_sq_vec
        x_hat_sq_k_vec[k] = sum([sum(x_hat_sq_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_sq_k_mat[:,1] = states
    x_hat_sq_k_mat[:,2] = genes
    x_hat_sq_k_mat[:,end] = recursive_flatten(x_hat_sq_k_vec)
    return x_hat_sq_k_mat

end
function tidy_update_x_hat_k(xmat,rtik_mat,ηtkj_mat)
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
    nrows = K*G; 
    ncols = 3;
    x_hat_k_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_k_vec = Vector{Vector{Float64}}(undef,K);
    x_vec =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    ηtkj_mat_tr = ηtkj_mat[1:2:end,:];
    ηtkj_vec_tr = [ηtkj_mat_tr[(i-1)*G+1:(i)*G,end] for i in 1:T*K];
    # all([all(a .== b) for (a,b) in zip(x_vec,vcat(x_to_use...)) ])
    # all([all(rtik_mat[k:K:end,end] .== [el[k] for el in vcat(rtik...)]) for k in 1:K])
    
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in states_id
        rtik_k =  rtik_mat[k:K:end,end];
        ηtkj_vec_tr_k = [repeat([el],val) for (el,val) in zip(ηtkj_vec_tr[k:K:end],N_t)]
        ηtkj_vec_tr_k = vcat(ηtkj_vec_tr_k...)
        x_hat_ = [r .*h .*d for (r,h,d) in zip(rtik_k,ηtkj_vec_tr_k,x_vec)]
        x_hat_k_vec[k] = sum([sum(x_hat_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_k_mat[:,1] = states
    x_hat_k_mat[:,2] = genes
    x_hat_k_mat[:,end] = recursive_flatten(x_hat_k_vec)
    return x_hat_k_mat
    # all(sum(x_hat_[1:600]) .== x_hat_tk[1][1])
    # all(sum(x_hat_[601:1200]) .== x_hat_tk[1][2])
    # all(sum(x_hat_[1201:1800]) .== x_hat_tk[1][3])
    # all(sum(x_hat_[2401:3000]) .== x_hat_tk[1][5])
end
function tidy_update_x_hat_e(xmat,rtik_mat,ηtkj_mat)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(rtik_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:1);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,1)
    nrows = 1*G; 
    ncols = 3;
    x_hat_e_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_k_vec = Vector{Vector{Float64}}(undef,K);
    x_vec =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    ηtkj_mat_fs = ηtkj_mat[2:2:end,:];
    ηtkj_vec_fs = [ηtkj_mat_fs[(i-1)*G+1:(i)*G,end] for i in 1:T*K];
    # all([all(a .== b) for (a,b) in zip(x_vec,vcat(x_to_use...)) ])
    # all([all(rtik_mat[k:K:end,end] .== [el[k] for el in vcat(rtik...)]) for k in 1:K])
    
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in 1:K
        rtik_k =  rtik_mat[k:K:end,end];
        ηtkj_vec_fs_k = [repeat([el],val) for (el,val) in zip(ηtkj_vec_fs[k:K:end],N_t)]
        ηtkj_vec_fs_k = vcat(ηtkj_vec_fs_k...)
        x_hat_ = [r .*h .*d for (r,h,d) in zip(rtik_k,ηtkj_vec_fs_k,x_vec)]
        x_hat_k_vec[k] = sum([sum(x_hat_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_e_mat[:,1] = states
    x_hat_e_mat[:,2] = genes
    x_hat_e_mat[:,end] = sum(x_hat_k_vec)
    return x_hat_e_mat
    # all(sum(x_hat_[1:600]) .== x_hat_tk[1][1])
    # all(sum(x_hat_[601:1200]) .== x_hat_tk[1][2])
    # all(sum(x_hat_[1201:1800]) .== x_hat_tk[1][3])
    # all(sum(x_hat_[2401:3000]) .== x_hat_tk[1][5])
end
function tidy_update_x_hat_sq_k(xmat,rtik_mat,ηtkj_mat)
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
    nrows = K*G; 
    ncols = 3;
    x_hat_sq_k_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_sq_k_vec = Vector{Vector{Float64}}(undef,K);
    x_sq_vec =  [(xmat[(i-1)*G+1:(i)*G,end]) .^ 2 for i in 1:N];
    ηtkj_mat_tr = ηtkj_mat[1:2:end,:];
    ηtkj_vec_tr = [ηtkj_mat_tr[(i-1)*G+1:(i)*G,end] for i in 1:T*K];

    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in states_id
        # rtik =  rtik_mat[k:K:end,end];
        rtik_k =  rtik_mat[k:K:end,end];
        ηtkj_vec_tr_k = [repeat([el],val) for (el,val) in zip(ηtkj_vec_tr[k:K:end],N_t)]
        ηtkj_vec_tr_k = vcat(ηtkj_vec_tr_k...)
        x_hat_sq_ = [r .*h .*d for (r,h,d) in zip(rtik_k,ηtkj_vec_tr_k,x_sq_vec)]#rtik .* x_sq_vec
        x_hat_sq_k_vec[k] = sum([sum(x_hat_sq_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_sq_k_mat[:,1] = states
    x_hat_sq_k_mat[:,2] = genes
    x_hat_sq_k_mat[:,end] = recursive_flatten(x_hat_sq_k_vec)
    return x_hat_sq_k_mat

end
function tidy_update_x_hat_sq_e(xmat,rtik_mat,ηtkj_mat)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(rtik_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:1);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,1)
    nrows = 1*G; 
    ncols = 3;
    x_hat_sq_e_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_sq_k_vec = Vector{Vector{Float64}}(undef,K);
    x_sq_vec =  [(xmat[(i-1)*G+1:(i)*G,end]) .^ 2 for i in 1:N];
    ηtkj_mat_fs = ηtkj_mat[2:2:end,:];
    ηtkj_vec_fs = [ηtkj_mat_fs[(i-1)*G+1:(i)*G,end] for i in 1:T*K];

    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in 1:K
        # rtik =  rtik_mat[k:K:end,end];
        rtik_k =  rtik_mat[k:K:end,end];
        ηtkj_vec_fs_k = [repeat([el],val) for (el,val) in zip(ηtkj_vec_fs[k:K:end],N_t)]
        ηtkj_vec_fs_k = vcat(ηtkj_vec_fs_k...)
        x_hat_sq_ = [r .*h .*d for (r,h,d) in zip(rtik_k,ηtkj_vec_fs_k,x_sq_vec)]#rtik .* x_sq_vec
        x_hat_sq_k_vec[k] = sum([sum(x_hat_sq_[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_sq_e_mat[:,1] = states
    x_hat_sq_e_mat[:,2] = genes
    x_hat_sq_e_mat[:,end] = sum(x_hat_sq_k_vec)
    return x_hat_sq_e_mat

end
function tidy_update_λ0k_hat!(rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    λ0_vec = λ0μ0a0b0_mat[:,2]
    Nk_mat = tidy_update_Nk(rtik_mat);
    Nk_vec = innermelt(Nk_mat[:,end],G)
    λ0_vec_vec = outermelt(λ0_vec,K)
    @views λ0kmka0kb0k_hat_mat[:,3] = λ0_vec_vec .+ Nk_vec
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_λ0k_hat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    λ0_vec = λ0μ0a0b0_mat[:,2]
    Nkj_mat = tidy_update_Nkj(rtik_mat,ηtkj_mat);
    λ0_vec_vec = outermelt(λ0_vec,K)
    @views λ0kmka0kb0k_hat_mat[:,3] = λ0_vec_vec .+ Nkj_mat[:,end]
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_λ0_err_hat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    G = length(unique(λ0ma0b0_err_hat_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    λ0_err_vec = λ0μ0a0b0_err_mat[:,2]
    Nej_mat = tidy_update_Nej(rtik_mat,ηtkj_mat);
    λ0_vec_vec = outermelt(λ0_err_vec,1)
    @views λ0ma0b0_err_hat_mat[:,3] = λ0_vec_vec .+ Nej_mat[:,end]
    return λ0ma0b0_err_hat_mat
end
function tidy_update_a0k_hat_usingXhat!(rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    Nk_mat = tidy_update_Nk(rtik_mat);
    Nk_vec = innermelt(Nk_mat[:,end],G)
    a0_vec = λ0μ0a0b0_mat[:,4]
    a0_vec_vec = outermelt(a0_vec,K)
    @views λ0kmka0kb0k_hat_mat[:,5] = a0_vec_vec .+  0.5 .* (Nk_vec .+ 1)
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_a0k_hat_usingXhat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    Nkj_mat = tidy_update_Nkj(rtik_mat,ηtkj_mat);
    a0_vec = λ0μ0a0b0_mat[:,4]
    a0_vec_vec = outermelt(a0_vec,K)
    @views λ0kmka0kb0k_hat_mat[:,5] = a0_vec_vec .+  0.5 .* (Nkj_mat[:,end] .+ 1)
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_a0_err_hat_usingXhat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    G = length(unique(λ0ma0b0_err_hat_mat[:,2]));
    K = length(unique(rtik_mat[:,3]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    Nej_mat = tidy_update_Nej(rtik_mat,ηtkj_mat);
    a0_err_vec = λ0μ0a0b0_err_mat[:,4]
    a0_vec_vec = outermelt(a0_err_vec,1)
    @views λ0ma0b0_err_hat_mat[:,5] = a0_vec_vec .+  0.5 .* (Nej_mat[:,end] .+ 1)
    return λ0ma0b0_err_hat_mat
end
function tidy_update_mk_hat_usingXhat!(x_hat_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(x_hat_k_mat[:,1]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    μ0_vec = λ0μ0a0b0_mat[:,3]
    λ0_vec = λ0μ0a0b0_mat[:,2]
    λ0_μ0 =  λ0_vec .* μ0_vec
    λ0k_hat_vec = λ0kmka0kb0k_hat_mat[:,3]
    # all(outermelt(λ0_μ0,K)[1:10] .== outermelt(λ0_μ0,K)[end-9:end])
    λ0_μ0_vec = outermelt(λ0_μ0,K)
    @views λ0kmka0kb0k_hat_mat[:,4] = ( λ0_μ0_vec .+ x_hat_k_mat[:,end] ) ./ λ0k_hat_vec
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_mk_hat_usingXhat!(x_hat_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(x_hat_k_mat[:,1]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    μ0_vec = λ0μ0a0b0_mat[:,3]
    λ0_vec = λ0μ0a0b0_mat[:,2]
    λ0_μ0 =  λ0_vec .* μ0_vec
    λ0k_hat_vec = λ0kmka0kb0k_hat_mat[:,3]
    # all(outermelt(λ0_μ0,K)[1:10] .== outermelt(λ0_μ0,K)[end-9:end])
    λ0_μ0_vec = outermelt(λ0_μ0,K)
    @views λ0kmka0kb0k_hat_mat[:,4] = ( λ0_μ0_vec .+ x_hat_k_mat[:,end] ) ./ λ0k_hat_vec
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_m_err_hat_usingXhat!(x_hat_e_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    G = length(unique(λ0ma0b0_err_hat_mat[:,2]));
    K = length(unique(x_hat_e_mat[:,1]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    μ0_err_vec = λ0μ0a0b0_err_mat[:,3]
    λ0_err_vec = λ0μ0a0b0_err_mat[:,2]
    λ0_μ0_err =  λ0_err_vec .* μ0_err_vec
    λ0_err_hat_vec = λ0ma0b0_err_hat_mat[:,3]
    # all(outermelt(λ0_μ0,K)[1:10] .== outermelt(λ0_μ0,K)[end-9:end])
    λ0_μ0_err_vec = outermelt(λ0_μ0_err,1)
    @views λ0ma0b0_err_hat_mat[:,4] = ( λ0_μ0_err_vec .+ x_hat_e_mat[:,end] ) ./ λ0_err_hat_vec
    return λ0ma0b0_err_hat_mat
end

function tidy_update_b0k_hat_usingXhat!(x_hat_sq_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(x_hat_sq_k_mat[:,1]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    μ0_vec = λ0μ0a0b0_mat[:,3]
    λ0_vec = λ0μ0a0b0_mat[:,2]
    b0_vec = λ0μ0a0b0_mat[:,5]
    μ0_sq_vec = μ0_vec .^2
    λ0_μ0_sq = μ0_sq_vec .*λ0_vec
    λ0k_hat_vec = λ0kmka0kb0k_hat_mat[:,3]
    mk_hat_vec = λ0kmka0kb0k_hat_mat[:,4]
    # all(outermelt(λ0_μ0,K)[1:10] .== outermelt(λ0_μ0,K)[end-9:end])
    λ0_μ0_sq_vec = outermelt(λ0_μ0_sq,K)
    b0_vec_vec = outermelt(b0_vec,K)
    @views λ0kmka0kb0k_hat_mat[:,6] = b0_vec_vec + 0.5 * (x_hat_sq_k_mat[:,end] .+ λ0_μ0_sq_vec .- λ0k_hat_vec .* (mk_hat_vec) .^2)
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_b0k_hat_usingXhat!(x_hat_sq_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    K = length(unique(x_hat_sq_k_mat[:,1]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    μ0_vec = λ0μ0a0b0_mat[:,3]
    λ0_vec = λ0μ0a0b0_mat[:,2]
    b0_vec = λ0μ0a0b0_mat[:,5]
    μ0_sq_vec = μ0_vec .^2
    λ0_μ0_sq = μ0_sq_vec .*λ0_vec
    λ0k_hat_vec = λ0kmka0kb0k_hat_mat[:,3]
    mk_hat_vec = λ0kmka0kb0k_hat_mat[:,4]
    # all(outermelt(λ0_μ0,K)[1:10] .== outermelt(λ0_μ0,K)[end-9:end])
    λ0_μ0_sq_vec = outermelt(λ0_μ0_sq,K)
    b0_vec_vec = outermelt(b0_vec,K)
    @views λ0kmka0kb0k_hat_mat[:,6] = b0_vec_vec + 0.5 * (x_hat_sq_k_mat[:,end] .+ λ0_μ0_sq_vec .- λ0k_hat_vec .* (mk_hat_vec) .^2)
    return λ0kmka0kb0k_hat_mat
end
function tidy_update_b0_err_hat_usingXhat!(x_hat_sq_e_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    G = length(unique(λ0ma0b0_err_hat_mat[:,2]));
    K = length(unique(x_hat_sq_e_mat[:,1]));
    states_id = collect(1:K)
    gene_id = collect(1:G)
    μ0_err_vec = λ0μ0a0b0_err_mat[:,3]
    λ0_err_vec = λ0μ0a0b0_err_mat[:,2]
    b0_err_vec = λ0μ0a0b0_err_mat[:,5]
    μ0_sq_err_vec = μ0_err_vec .^2
    λ0_μ0_sq_err = μ0_sq_err_vec .*λ0_err_vec
    λ0_err_hat_vec = λ0ma0b0_err_hat_mat[:,3]
    m_err_hat_vec = λ0ma0b0_err_hat_mat[:,4]
    # all(outermelt(λ0_μ0,K)[1:10] .== outermelt(λ0_μ0,K)[end-9:end])
    λ0_μ0_sq_err_vec = outermelt(λ0_μ0_sq_err,1)
    b0_vec_vec = outermelt(b0_err_vec,1)
    @views λ0ma0b0_err_hat_mat[:,6] = b0_vec_vec + 0.5 * (x_hat_sq_e_mat[:,end] .+ λ0_μ0_sq_err_vec .- λ0_err_hat_vec .* (m_err_hat_vec) .^2)
    return λ0ma0b0_err_hat_mat
end
function tidy_update_ηtkj_SparseVS(rtik_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,ηtkj_prior_mat)
    T = length(unique(e_τ_μ_tikj_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
    Glog = G*log(2π)
    logpi= Glog/G
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    gene_ids = collect(1:G)
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)


    rtik_cell = innermelt(rtik_mat[:,end],G)
    e_log_τ_kj_cell = outermelt(e_log_τ_kj_mat[:,end],N)
    e_τ_μ_tikj_cell = e_τ_μ_tikj_mat[:,end]
    cell_signal_sq = rtik_cell .* (0.5 .* (e_log_τ_kj_cell .- logpi .-  e_τ_μ_tikj_cell))
    cell_signal_sq_lin = [cell_signal_sq[(i-1)*K*G+1:(i)*K*G,end] for i in 1:N]
    prior_odds = [log.(ηtkj_prior_mat[1:2:end,end])[(t-1)*K*G+1:(t)*K*G] - log.(ηtkj_prior_mat[2:2:end,end])[(t-1)*K*G+1:(t)*K*G] for t in 1:T]
    cell_signal_imp = [sum(cell_signal_sq_lin[st:en]) for (st,en) in timeranges] .+ prior_odds
    # cell_signal_imp = recursive_flatten(cell_signal_imp)
    # e_log_τj_err_cell = outermelt(e_log_τj_err_mat[:,end],N*K)
    # e_τ_μ_tij_err_cell = recursive_flatten(outermelt.([e_τ_μ_tij_err_mat[(i-1)*G+1:i*G,end] for i in 1:N],K))
    # cell_err_sq = rtik_cell .* (0.5 .* (e_log_τj_err_cell .- logpi .-  e_τ_μ_tij_err_cell))
    # cell_err_sq_lin = [cell_err_sq[(i-1)*K*G+1:(i)*K*G,end] for i in 1:N ]
    # cell_signal_notimp = [sum(cell_err_sq_lin[st:en]) for (st,en) in timeranges] .+ [log.(ηtkj_prior_mat[2:2:end,end])[(t-1)*K*G+1:(t)*K*G] for t in 1:T]
    # cell_signal_notimp = recursive_flatten(cell_signal_notimp)







    log_ηtkj_tilde = recursive_flatten(cell_signal_imp) #[[imp,notimp] for (imp,notimp) in zip(cell_signal_imp,cell_signal_notimp)]




    ## OLD 
    # ηtkj_imp = logistic.(log_ηtkj_tilde)#[norm_weights(el) for el in log_ηtkj_tilde]
    # ηtkj_notimp = 1 .- ηtkj_imp
    # ηtkj_ = [[imp,notimp] for (imp,notimp) in zip(ηtkj_imp,ηtkj_notimp)]
    # ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    # log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],log_ηtkj_tilde  ])

    ## NEW (Multinomial)
    log_ηtkj_tilde_vec = [log_ηtkj_tilde[(i-1)*G+1:i*G,end] for i in 1:T*K]
    ηtkj_vec_ = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        tjk_mat = Matrix{Float64}(undef,K,G)
        for j in 1:G
            tjk_mat[:,j] = norm_weights([el[j] for  el in log_ηtkj_tilde_vec[(t-1)*K+1:t*K]])
            
        end
        ηtkj_vec_[t] = [vec(row) for row in eachrow(tjk_mat)]
    end
    ηtkj_imp = recursive_flatten(ηtkj_vec_)
    ηtkj_notimp = 1 .- ηtkj_imp
    ηtkj_ = [[imp,notimp] for (imp,notimp) in zip(ηtkj_imp,ηtkj_notimp)]
    ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],log_ηtkj_tilde  ])

    #NEW Alternative
    # log_ηtkj_tilde_vec
    # log_ηtkj_tilde_vec = [log_ηtkj_tilde[(i-1)*G+1:i*G,end] for i in 1:T*K]
    # ηtkj_vec_ = norm_weights.(log_ηtkj_tilde_vec)
    # ηtkj_imp = recursive_flatten(ηtkj_vec_)
    # ηtkj_notimp = 1 .- ηtkj_imp
    # ηtkj_ = [[imp,notimp] for (imp,notimp) in zip(ηtkj_imp,ηtkj_notimp)]
    # ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    # log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],log_ηtkj_tilde  ])

    return ηtkj_mat,log_ηtkj_tilde_mat
end
function tidy_update_ηtkj_VS(rtik_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,e_log_τj_err_mat,e_τ_μ_tij_err_mat,ηtkj_prior_mat)
    T = length(unique(e_τ_μ_tikj_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
    Glog = G*log(2π)
    logpi= Glog/G
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    gene_ids = collect(1:G)
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)


    rtik_cell = innermelt(rtik_mat[:,end],G)
    e_log_τ_kj_cell = outermelt(e_log_τ_kj_mat[:,end],N)
    e_τ_μ_tikj_cell = e_τ_μ_tikj_mat[:,end]
    e_log_τj_err_cell =  outermelt(e_log_τj_err_mat[:,end],K*N)
    e_τ_μ_tij_err_vec = [e_τ_μ_tij_err_mat[(i-1)*G+1:i*G,5] for i in 1:N]
    e_τ_μ_tij_err_cell = recursive_flatten([repeat([el],K) for el in e_τ_μ_tij_err_vec])

    #OLD
    cell_signal_sq = rtik_cell .* (0.5 .* (e_log_τ_kj_cell  .-  e_τ_μ_tikj_cell .- e_log_τj_err_cell .+  e_τ_μ_tij_err_cell))
    cell_signal_sq_vec = [cell_signal_sq[(i-1)*G+1:i*G] for i in 1:N*K]
    ηtkj_prior_mat_tr = ηtkj_prior_mat[1:2:end,:]
    ηtkj_prior_mat_fs = ηtkj_prior_mat[2:2:end,:]
    ηtkj_prior_logged_trfs_time_state = log.(ηtkj_prior_mat_tr[:,end]) .- log.(ηtkj_prior_mat_fs[:,end])
    ηtkj_prior_logged_trfs_vec = [ηtkj_prior_logged_trfs_time_state[(i-1)*G+1:i*G] for i in 1:K*T]

    ηtkj_prior_timeranges = tidy_get_timeranges(ones(Int,T) .* K) 
    ηtkj_prior_timeranges_tuple_vec = [(st,en) for (st,en) in ηtkj_prior_timeranges]
    t_indx = 0
    eta_sum = nothing
    ηtkj_vec = Vector{Vector{Vector}}(undef,T)
    for (st,en) in timeranges
        t_indx +=1
        ts = first(ηtkj_prior_timeranges_tuple_vec[t_indx])
        te = last(ηtkj_prior_timeranges_tuple_vec[t_indx])
        eta_sum = [sum(cell_signal_sq_vec[k:K:end][st:en] )  .+ ηtkj_prior_logged_trfs_vec[ts:te][k:K:end][1] for k in state_ids]
        ηtkj_vec[t_indx] = eta_sum 
    end
    log_ηtkj_tilde = recursive_flatten(ηtkj_vec)




    ηtkj_imp = logistic.(log_ηtkj_tilde)
    ηtkj_notimp = 1 .- ηtkj_imp
    ηtkj_ = [[imp,notimp] for (imp,notimp) in zip(ηtkj_imp,ηtkj_notimp)]
    ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],log_ηtkj_tilde  ])


    ########New
    # cell_signal_sq = rtik_cell .* (0.5 .* (e_log_τ_kj_cell .- logpi .-  e_τ_μ_tikj_cell))
    # cell_noise_sq = rtik_cell .* (0.5 .* (e_log_τj_err_cell .-logpi .-  e_τ_μ_tij_err_cell))
    # cell_signal_sq_vec = [cell_signal_sq[(i-1)*G+1:i*G] for i in 1:N*K]
    # cell_noise_sq_vec = [cell_noise_sq[(i-1)*G+1:i*G] for i in 1:N*K]
    # ηtkj_prior_mat_tr = ηtkj_prior_mat[1:2:end,:]
    # ηtkj_prior_mat_fs = ηtkj_prior_mat[2:2:end,:]
    # ηtkj_prior_logged_tr_time_state = log.(ηtkj_prior_mat_tr[:,end]) 
    # ηtkj_prior_logged_fs_time_state = log.(ηtkj_prior_mat_fs[:,end])
    # ηtkj_prior_logged_tr_vec = [ηtkj_prior_logged_tr_time_state[(i-1)*G+1:i*G] for i in 1:K*T]
    # ηtkj_prior_logged_fs_vec = [ηtkj_prior_logged_fs_time_state[(i-1)*G+1:i*G] for i in 1:K*T]

    # ηtkj_prior_timeranges = tidy_get_timeranges(ones(Int,T) .* K) 
    # ηtkj_prior_timeranges_tuple_vec = [(st,en) for (st,en) in ηtkj_prior_timeranges]
    # t_indx = 0

    # eta_signal_sum = nothing
    # eta_noise_sum= nothing

    # ηtkj_vec_imp = Vector{Vector{Vector}}(undef,T)
    # ηtkj_vec_notimp = Vector{Vector{Vector}}(undef,T)
    # for (st,en) in timeranges
    #     t_indx +=1
    #     ts = first(ηtkj_prior_timeranges_tuple_vec[t_indx])
    #     te = last(ηtkj_prior_timeranges_tuple_vec[t_indx])

    #     eta_signal_sum = [sum(cell_signal_sq_vec[k:K:end][st:en] )  .+ ηtkj_prior_logged_tr_vec[ts:te][k:K:end][1] for k in state_ids]
    #     eta_noise_sum = [sum(cell_noise_sq_vec[k:K:end][st:en] )  .+ ηtkj_prior_logged_fs_vec[ts:te][k:K:end][1] for k in state_ids]
    #     ηtkj_vec_imp[t_indx] = eta_signal_sum
    #     ηtkj_vec_notimp[t_indx] = eta_noise_sum
    # end
    # log_ηtkj_tilde_imp = recursive_flatten(ηtkj_vec_imp)
    # log_ηtkj_tilde_notimp = recursive_flatten(ηtkj_vec_notimp)
    # # ηtkj_imp = logistic.(log_ηtkj_tilde)
    # # ηtkj_notimp = 1 .- ηtkj_imp
    # ηtkj_ = [norm_weights([imp,notimp]) for (imp,notimp) in zip(log_ηtkj_tilde_imp,log_ηtkj_tilde_notimp)]
    # ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    # log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],log_ηtkj_tilde  ])
    

    # NEW Alternative
    # log_ηtkj_tilde_vec = [log_ηtkj_tilde[(i-1)*G+1:i*G,end] for i in 1:T*K]
    # ηtkj_vec_ = Vector{Vector{Vector{Float64}}}(undef,T)
    # for t in 1:T
    #     tjk_mat = Matrix{Float64}(undef,K,G)
    #     for j in 1:G
    #         tjk_mat[:,j] = norm_weights([el[j] for  el in log_ηtkj_tilde_vec[(t-1)*K+1:t*K]])
            
    #     end
    #     ηtkj_vec_[t] = [vec(row) for row in eachrow(tjk_mat)]
    # end
    # ηtkj_imp = recursive_flatten(ηtkj_vec_)
    # ηtkj_notimp = 1 .- ηtkj_imp
    # ηtkj_ = [[imp,notimp] for (imp,notimp) in zip(ηtkj_imp,ηtkj_notimp)]
    # ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    # log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],log_ηtkj_tilde  ])
    return ηtkj_mat,log_ηtkj_tilde_mat
end
function tidy_update_αt!(aαbαawbwaγbγ_mat,ρkωkckdk_hat_mat,θ_hat_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(θ_hat_mat[:,1]))
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    Kplus = K + 1
    a_α = aαbαawbwaγbγ_mat[1,2]
    b_α = aαbαawbwaγbγ_mat[1,3]
    a_α_hat = a_α + K
    
    
    e_log_π_mat =tidy_log_π_expected_value(θ_hat_mat)
    e_βk_mat = tidy_βk_expected_value(ρkωkckdk_hat_mat)
    e_βk_vec_ = outermelt(e_βk_mat[:,end],T)
    e_log_π_vec_ = e_log_π_mat[:,end]
    e_βk_log_π_ = e_βk_vec_ .* e_log_π_vec_

    # a_α_hat_vec = outermelt(a_α_hat,T)
    # b_α_hat_vec =  b_α .- [sum(e_βk_log_π_[(t-1)*Kplus+1:t*Kplus]) for t in 1:T]
    @views aαtbαtawtbwt_hat_mat[:,2] = outermelt(a_α_hat,T)
    @views aαtbαtawtbwt_hat_mat[:,3] =  b_α .- [sum(e_βk_log_π_[(t-1)*Kplus+1:t*Kplus]) for t in 1:T]
    return aαtbαtawtbwt_hat_mat

end
function tidy_update_wtildet!(aαbαawbwaγbγ_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(c_ttprime_mat[:,1]))
    adot_w = aαbαawbwaγbγ_mat[1,4]
    bdot_w = aαbαawbwaγbγ_mat[1,5]
    c_ttprime_ = [c_ttprime_mat[(t-1)*T+1:t*T,end] for t in 1:T]
    @views aαtbαtawtbwt_hat_mat[2:end,4] = update_awt_hat(adot_w, c_ttprime_)
    @views aαtbαtawtbwt_hat_mat[2:end,5] = update_bwt_hat(bdot_w, c_ttprime_)
    return aαtbαtawtbwt_hat_mat
end
function tidy_update_γ(ρkωkckdk_hat_mat,aαbαawbwaγbγ_mat)
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    e_log_minus_uk_mat = tidy_log1minusUk_expected_value(ρkωkckdk_hat_mat);
    a_γ = aαbαawbwaγbγ_mat[1,6]
    b_γ = aαbαawbwaγbγ_mat[1,7]
    a_γ_hat = a_γ + K
    b_γ_hat = b_γ - sum(e_log_minus_uk_mat[:,end])
    nrows = 1
    ncols = 3
    aγbγ_hat_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    aγbγ_hat_mat[1,1] = 1
    aγbγ_hat_mat[1,2] = a_γ_hat
    aγbγ_hat_mat[1,3] = b_γ_hat
    return aγbγ_hat_mat
end
function tidy_update_Tαk(θ_hat_mat,aαtbαtawtbwt_hat_mat)
    T = length(unique(aαtbαtawtbwt_hat_mat[:,1]))
    Kplus = length(unique(θ_hat_mat[:,2]))
    timepoints = collect(1:T)
    kplus_id = collect(1:Kplus)
    e_log_π_mat = tidy_log_π_expected_value(θ_hat_mat)
    e_αt_mat = tidy_αt_expected_value(aαtbαtawtbwt_hat_mat)
    e_αt_vec_ = innermelt(e_αt_mat[:,2],Kplus)
    e_αt_log_π_ = e_αt_vec_ .* e_log_π_mat[:,end]
    nrows = Kplus
    ncols = 2
    Tαk_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    Tαk_mat[:,1] = kplus_id
    Tαk_mat[:,2] = sum([e_αt_log_π_[(t-1)*Kplus+1:t*Kplus] for t in timepoints])

    return Tαk_mat
    
end
function tidy_update_rho_omega_hat!(T,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat;optim_max_iter=100000) ###SLIGHTLY DIFFERENT
    rhok_hat_vec_ = ρkωkckdk_hat_mat[:,2]
    omegak_hat_vec_ = ρkωkckdk_hat_mat[:,3]
    Tαk_ = Tαk_mat[:,2]
    e_γ_ = e_γ_mat[1,2]
    c_hat = StatsFuns.logit.(rhok_hat_vec_)
    d_hat =  log.(omegak_hat_vec_)
    cd_vec_args = [c_hat , d_hat]
    cd_vec_args0  = permutedims(reduce(hcat,cd_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,e_γ_,Tαk_)
    gg_uncon! = g_unconstrained_closure!(T,e_γ_,Tαk_)
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(), Optim.Options(iterations = optim_max_iter))
    # lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter)) #,linesearch=Optim.LineSearches.BackTracking() 
    # optimize(LB_LG_unconstrained,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))

    # @debug Optim.converged(lb_lg_results)
    new_rho_hat = sigmoid.(lb_lg_results.res.minimizer[1,:])
    new_omega_hat = exp.(lb_lg_results.res.minimizer[2,:])

    new_c_hat = StatsFuns.logit.(new_rho_hat)
    new_d_hat = log.(new_omega_hat)
    @views ρkωkckdk_hat_mat[:,2] = new_rho_hat
    @views ρkωkckdk_hat_mat[:,3] = new_omega_hat
    @views ρkωkckdk_hat_mat[:,4] = new_c_hat
    @views ρkωkckdk_hat_mat[:,5] = new_d_hat 
    return ρkωkckdk_hat_mat
end
#####################
#####################
function tidy_local_update(xmat,θ_hat_mat,λ0kmka0kb0k_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    Glog = G*log(2π)
    e_log_π_mat =  tidy_log_π_expected_value(θ_hat_mat; e_log_π_mat_init=e_log_π_mat_init)
    e_log_τ_kj_mat = tidy_log_τ_kj_expected_value(λ0kmka0kb0k_hat_mat; e_log_τ_kj_mat_init = e_log_τ_kj_mat_init )
    e_τ_μ_tikj_mat = tidy_τ_μ_expected_value(xmat,λ0kmka0kb0k_hat_mat);
    rtik_mat = tidy_update_rtik(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat; rtik_mat_init = rtik_mat_init);
    Ntk_mat = tidy_update_Ntk(rtik_mat);
    c_ttprime_mat = tidy_update_c_ttprime(aαtbαtawtbwt_hat_mat,rtik_mat,θ_hat_mat);
    θ_hat_mat =tidy_update_θ_hat(Ntk_mat,ρkωkckdk_hat_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)

    return rtik_mat, θ_hat_mat, c_ttprime_mat
end
function tidy_global_update(xmat,rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = 10000)
    
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    x_hat_k_mat =tidy_update_x_hat_k(xmat,rtik_mat);
    x_hat_sq_k_mat = tidy_update_x_hat_sq_k(xmat,rtik_mat);
    tidy_update_λ0k_hat!(rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_a0k_hat_usingXhat!(rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_mk_hat_usingXhat!(x_hat_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_b0k_hat_usingXhat!(x_hat_sq_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_αt!(aαbαawbwaγbγ_mat,ρkωkckdk_hat_mat,θ_hat_mat,aαtbαtawtbwt_hat_mat)
    tidy_update_wtildet!(aαbαawbwaγbγ_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat)
    aγbγ_hat_mat = tidy_update_γ(ρkωkckdk_hat_mat,aαbαawbwaγbγ_mat)
    e_γ_mat = tidy_γ_expected_value(aγbγ_hat_mat)
    Tαk_mat = tidy_update_Tαk(θ_hat_mat,aαtbαtawtbwt_hat_mat)
    tidy_update_rho_omega_hat!(T,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat;optim_max_iter=optim_max_iter)

    return x_hat_k_mat, x_hat_sq_k_mat,λ0kmka0kb0k_hat_mat, aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat
end
function tidy_local_update_SparseVS(xmat,rtik_mat,ηtkj_mat,θ_hat_mat,λ0kmka0kb0k_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat,ηtkj_prior_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    Glog = G*log(2π)
    e_log_π_mat =  tidy_log_π_expected_value(θ_hat_mat; e_log_π_mat_init=e_log_π_mat_init)
    e_log_τ_kj_mat = tidy_log_τ_kj_expected_value(λ0kmka0kb0k_hat_mat; e_log_τ_kj_mat_init = e_log_τ_kj_mat_init )
    e_τ_μ_tikj_mat = tidy_τ_μ_expected_value(xmat,λ0kmka0kb0k_hat_mat);
    ηtkj_mat,log_ηtkj_tilde_mat = tidy_update_ηtkj_SparseVS(rtik_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,ηtkj_prior_mat)
    rtik_mat = tidy_update_rtik_SparseVS(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat,ηtkj_mat; rtik_mat_init = rtik_mat_init );
    Ntk_mat = tidy_update_Ntk(rtik_mat);
    c_ttprime_mat = tidy_update_c_ttprime(aαtbαtawtbwt_hat_mat,rtik_mat,θ_hat_mat);
    θ_hat_mat =tidy_update_θ_hat(Ntk_mat,ρkωkckdk_hat_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)

    return rtik_mat, ηtkj_mat, θ_hat_mat, c_ttprime_mat,log_ηtkj_tilde_mat
end
function tidy_global_update_SparseVS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = 10000)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    x_hat_k_mat =tidy_update_x_hat_k(xmat,rtik_mat,ηtkj_mat);
    x_hat_sq_k_mat = tidy_update_x_hat_sq_k(xmat,rtik_mat,ηtkj_mat);
    tidy_update_λ0k_hat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_a0k_hat_usingXhat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_mk_hat_usingXhat!(x_hat_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_b0k_hat_usingXhat!(x_hat_sq_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_αt!(aαbαawbwaγbγ_mat,ρkωkckdk_hat_mat,θ_hat_mat,aαtbαtawtbwt_hat_mat)
    tidy_update_wtildet!(aαbαawbwaγbγ_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat)
    aγbγ_hat_mat = tidy_update_γ(ρkωkckdk_hat_mat,aαbαawbwaγbγ_mat)
    e_γ_mat = tidy_γ_expected_value(aγbγ_hat_mat)
    Tαk_mat = tidy_update_Tαk(θ_hat_mat,aαtbαtawtbwt_hat_mat)
    tidy_update_rho_omega_hat!(T,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat;optim_max_iter=optim_max_iter)

    return x_hat_k_mat, x_hat_sq_k_mat,λ0kmka0kb0k_hat_mat, aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat
end
function tidy_local_update_VS1(xmat,rtik_mat,ηtkj_mat,θ_hat_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat,ηtkj_prior_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing, e_log_τj_err_mat_init = nothing,e_τ_μ_tij_err_mat_init = nothing)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    Glog = G*log(2π)
    e_log_π_mat =  tidy_log_π_expected_value(θ_hat_mat; e_log_π_mat_init=e_log_π_mat_init)
    e_log_τ_kj_mat = tidy_log_τ_kj_expected_value(λ0kmka0kb0k_hat_mat; e_log_τ_kj_mat_init = e_log_τ_kj_mat_init )
    e_τ_μ_tikj_mat = tidy_τ_μ_expected_value(xmat,λ0kmka0kb0k_hat_mat);
    e_log_τj_err_mat = tidy_log_τj_err_expected_value(λ0ma0b0_err_hat_mat; e_log_τj_err_mat_init = e_log_τj_err_mat_init )
    e_τ_μ_tij_err_mat = tidy_τ_μ_err_expected_value(xmat,λ0ma0b0_err_hat_mat; e_τ_μ_tij_err_mat_init = e_τ_μ_tij_err_mat_init)
    ηtkj_mat,log_ηtkj_tilde_mat = tidy_update_ηtkj_VS(rtik_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,e_log_τj_err_mat,e_τ_μ_tij_err_mat,ηtkj_prior_mat)
    rtik_mat = tidy_update_rtik_VS1(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat,e_log_τj_err_mat,e_τ_μ_tij_err_mat,ηtkj_mat; rtik_mat_init = nothing );
    Ntk_mat = tidy_update_Ntk(rtik_mat);
    c_ttprime_mat = tidy_update_c_ttprime(aαtbαtawtbwt_hat_mat,rtik_mat,θ_hat_mat);
    θ_hat_mat =tidy_update_θ_hat(Ntk_mat,ρkωkckdk_hat_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)

    return rtik_mat, ηtkj_mat, θ_hat_mat, c_ttprime_mat,log_ηtkj_tilde_mat
end
function tidy_local_update_VS2(xmat,rtik_mat,ηtkj_mat,θ_hat_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat,ηtkj_prior_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing, e_log_τj_err_mat_init = nothing,e_τ_μ_tij_err_mat_init = nothing)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    Glog = G*log(2π)
    e_log_π_mat =  tidy_log_π_expected_value(θ_hat_mat; e_log_π_mat_init=e_log_π_mat_init)
    e_log_τ_kj_mat = tidy_log_τ_kj_expected_value(λ0kmka0kb0k_hat_mat; e_log_τ_kj_mat_init = e_log_τ_kj_mat_init )
    e_τ_μ_tikj_mat = tidy_τ_μ_expected_value(xmat,λ0kmka0kb0k_hat_mat);
    e_log_τj_err_mat = tidy_log_τj_err_expected_value(λ0ma0b0_err_hat_mat; e_log_τj_err_mat_init = e_log_τj_err_mat_init )
    e_τ_μ_tij_err_mat = tidy_τ_μ_err_expected_value(xmat,λ0ma0b0_err_hat_mat; e_τ_μ_tij_err_mat_init = e_τ_μ_tij_err_mat_init)
    ηtkj_mat,log_ηtkj_tilde_mat = tidy_update_ηtkj_VS(rtik_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,e_log_τj_err_mat,e_τ_μ_tij_err_mat,ηtkj_prior_mat)
    rtik_mat = tidy_update_rtik_VS2(e_log_π_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,c_ttprime_mat,ηtkj_mat; rtik_mat_init = nothing );
    Ntk_mat = tidy_update_Ntk(rtik_mat);
    c_ttprime_mat = tidy_update_c_ttprime(aαtbαtawtbwt_hat_mat,rtik_mat,θ_hat_mat);
    θ_hat_mat =tidy_update_θ_hat(Ntk_mat,ρkωkckdk_hat_mat,aαtbαtawtbwt_hat_mat,c_ttprime_mat)

    return rtik_mat, ηtkj_mat, θ_hat_mat, c_ttprime_mat,log_ηtkj_tilde_mat
end
function tidy_global_update_VS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = 10000)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    x_hat_k_mat =tidy_update_x_hat_k(xmat,rtik_mat,ηtkj_mat);
    x_hat_sq_k_mat = tidy_update_x_hat_sq_k(xmat,rtik_mat,ηtkj_mat);

    Nej_mat = tidy_update_Nej(rtik_mat,ηtkj_mat);
    x_hat_e_mat = tidy_update_x_hat_e(xmat,rtik_mat,ηtkj_mat)
    x_hat_sq_e_mat = tidy_update_x_hat_sq_e(xmat,rtik_mat,ηtkj_mat)


    tidy_update_λ0k_hat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_a0k_hat_usingXhat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_mk_hat_usingXhat!(x_hat_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)
    tidy_update_b0k_hat_usingXhat!(x_hat_sq_k_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat)

    tidy_update_λ0_err_hat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    tidy_update_a0_err_hat_usingXhat!(rtik_mat,ηtkj_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    tidy_update_m_err_hat_usingXhat!(x_hat_e_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)
    tidy_update_b0_err_hat_usingXhat!(x_hat_sq_e_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat)

    tidy_update_αt!(aαbαawbwaγbγ_mat,ρkωkckdk_hat_mat,θ_hat_mat,aαtbαtawtbwt_hat_mat)
    tidy_update_wtildet!(aαbαawbwaγbγ_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat)
    aγbγ_hat_mat = tidy_update_γ(ρkωkckdk_hat_mat,aαbαawbwaγbγ_mat)
    e_γ_mat = tidy_γ_expected_value(aγbγ_hat_mat)
    Tαk_mat = tidy_update_Tαk(θ_hat_mat,aαtbαtawtbwt_hat_mat)
    tidy_update_rho_omega_hat!(T,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat;optim_max_iter=optim_max_iter)

    return x_hat_k_mat, x_hat_sq_k_mat,x_hat_e_mat,x_hat_sq_e_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat
end


#####################
#####################
#####################################################
############# DEPRACATED FUNCTIONS ##################
#####################################################
function depracated_tidy_update_N(rtik_mat,ηtkj_mat)
    T = length(unique(rtik_mat[:,1]))
    N = length(unique(rtik_mat[:,2]))
    K = length(unique(rtik_mat[:,3]))
    G = length(unique(ηtkj_mat[:,3]))
    
    nrows = 2*G*K*N
    ncols = 6
    N_mat = Matrix{Union{Float64,Int64}}(undef,nrows,ncols)
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    gene_ids = collect(1:G)
    cell_ids = collect(1:N)
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    choice_ids = collect(1:2)
    time = innermelt(timepoints, 2*G*K*N_t)
    cells = innermelt(cell_ids,2*G*K)
    states = outermelt(innermelt(state_ids,2*G),N)
    genes = outermelt(outermelt(innermelt(gene_ids,2),K),N)
    importance_choice = outermelt(choice_ids,G*K*N)
    rtik_cell = recursive_flatten(innermelt.([rtik_mat[(i-1)*K+1:i*K,end] for i in 1:N],2*G))#recursive_flatten([rtik_mat[(i-1)*K+1:i*K,end] for i in 1:N for j in 1:2*G])#
    ηtkj_cell = recursive_flatten(outermelt.([ηtkj_mat[(t-1)*2*G*K+1:t*2*G*K,5] for t in 1:T],N_t))#recursive_flatten([[ηtkj_mat[(t-1)*2*G*K+1:t*2*G*K,5] for i in 1:N_t[t]] for t in 1:T])#
    
    @views N_mat[:,1] = time
    @views N_mat[:,2] = cells
    @views N_mat[:,3] = states
    @views N_mat[:,4] = genes
    @views N_mat[:,5] = importance_choice
    @views N_mat[:,6] = rtik_cell .* ηtkj_cell

    return N_mat
    
end
function depracated_tidy_update_N2(rtik_mat,ηtkj_mat)
    T = length(unique(rtik_mat[:,1]))
    N = length(unique(rtik_mat[:,2]))
    K = length(unique(rtik_mat[:,3]))
    G = length(unique(ηtkj_mat[:,3]))
    
    nrows = G*K*N
    ncols = 6
    N_mat = Matrix{Union{Float64,Int64}}(undef,nrows,ncols)
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    gene_ids = collect(1:G)
    cell_ids = collect(1:N)
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    # choice_ids = collect(1:2)
    time = innermelt(timepoints, G*K*N_t)
    cells = innermelt(cell_ids,G*K)
    states = outermelt(innermelt(state_ids,G),N)
    genes = outermelt(gene_ids,K*N)
    # importance_choice = outermelt(choice_ids,G*K*N)
    rtik_cell = [[rtik_mat[(t-1)*K*N_t[t]+1:t*K*N_t[t],end][(i-1)*K+1:i*K] for i in 1:N_t[t]] for t in 1:T]#innermelt(rtik_mat[:,end],G) #recursive_flatten([rtik_mat[(i-1)*K+1:i*K,end] for i in 1:N for j in 1:G])
    ηtkj_imp = ηtkj_mat[1:2:end,:]
    ηtkj_notimp = ηtkj_mat[2:2:end,:]
    recursive_flatten(outermelt.([ηtkj_mat[(t-1)*2*G*K+1:t*2*G*K,5] for t in 1:T],N_t))#recursive_flatten([[ηtkj_mat[(t-1)*2*G*K+1:t*2*G*K,5] for i in 1:N_t[t]] for t in 1:T])#
    ηtkj_imp_cell = [[ηtkj_imp[(t-1)*K*G+1:t*G*K,5][(k-1)*G+1:k*G] for k in 1:K] for t in 1:T]
    ηtkj_notimp_cell = [[ηtkj_notimp[(t-1)*K*G+1:t*G*K,5][(k-1)*G+1:k*G] for k in 1:K] for t in 1:T]
    @views N_mat[:,1] = time
    @views N_mat[:,2] = cells
    @views N_mat[:,3] = states
    @views N_mat[:,4] = genes
    @views N_mat[:,5] = recursive_flatten([rtik_cell[t][i][k] .* ηtkj_imp_cell[t][k]  for t in 1:T for i in 1:N_t[t] for k in 1:K])
    @views N_mat[:,6] = recursive_flatten([rtik_cell[t][i][k] .* ηtkj_notimp_cell[t][k]  for t in 1:T for i in 1:N_t[t] for k in 1:K])

    return N_mat
    
end
function depracated_tidy_update_errorNj(N_mat)
    G = length(unique(N_mat[:,4]))
    nrows = G
    ncols = 2
    Nj_err_mat = Matrix{Union{Float64,Int64}}(undef,nrows,ncols)
    gene_ids = collect(1:G)
    genes = gene_ids
    @views N_err_mat = N_mat[2:2:end,:]
    @views Nj_err = [sum(N_err_mat[j:G:end,end]) for j in 1:G]
    @views Nj_err_mat[:,1] = genes
    @views Nj_err_mat[:,2] = Nj_err

    return Nj_err_mat
end
function depracated_tidy_update_signalNkj(N_mat)
    K = length(unique(N_mat[:,3]))
    G = length(unique(N_mat[:,4]))
    nrows = K*G
    ncols = 3
    Nkj_signal_mat = Matrix{Union{Float64,Int64}}(undef,nrows,ncols)
    gene_ids = collect(1:G)
    state_ids = collect(1:K)
    genes =outermelt(gene_ids,K)
    states =innermelt(state_ids,G)
    @views N_signal_mat = N_mat[1:2:end,:]
    @views Nkj_signal = [sum([N_signal_mat[(i-1)*K*G+1:K*G*i,:][(k-1)*G+1:k*G,end] for i in 1:N]) for k in 1:K]
    @views Nkj_signal_mat[:,1] = states
    @views Nkj_signal_mat[:,2] = genes
    @views Nkj_signal_mat[:,3] = recursive_flatten(Nkj_signal)

    return Nkj_signal_mat
end
function depracated_tidy_update_x_hat_error(xmat,N_mat)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(N_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,1)
    nrows = G; 
    ncols = 3;
    N_error_mat = N_mat[2:2:end,:]
    x_hat_err_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_err_vec = Vector{Vector{Float64}}(undef,K);
    x_vec =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    x_vec =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N for k in 1:K];
    # all([all(a .== b) for (a,b) in zip(x_vec,vcat(x_to_use...)) ])
    # all([all(rtik_mat[k:K:end,end] .== [el[k] for el in vcat(rtik...)]) for k in 1:K])
    # [sum([N_signal_mat[(i-1)*K*G+1:K*G*i,:][(k-1)*G+1:k*G,end] for i in 1:N]) for k in 1:K]
    #     [[N_error_mat[(i-1)*K*G+1:K*G*i,end][(k-1)*G+1:k*G] for i in 1:N] for k in 1:K]

    N_error_ = [N_error_mat[(i-1)*G+1:i*G,6] for i in 1:N*K]
    # Nj_error_k = [N_error_[k:K:end,end] for k in 1:K];

    # timeranges = zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    
    # for k in states_id
    #     Nj_error_ =  Nj_error_k[k];
    #     x_hat_err_ = broadcast( .*, Nj_error_, x_vec)
    #     x_hat_err_vec[k] = sum([x_hat_err_[st:en]  for (st,en) in timeranges])#sum(x_hat_err)# 
    # end
    x_hat_err_vec = broadcast( .*, N_error_, x_vec)
    x_hat_err_vec = sum(x_hat_err_vec)
    x_hat_err_mat[:,1] .= 1
    x_hat_err_mat[:,2] = genes
    x_hat_err_mat[:,end] = recursive_flatten(x_hat_err_vec)
    return x_hat_err_mat
end
function depracated_tidy_update_x_hatk_signal(xmat,N_mat)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = length(unique(N_mat[:,3]));
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,K)
    nrows = K*G; 
    ncols = 3;
    @views N_signal_mat = N_mat[1:2:end,:]
    x_hatk_signal_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hatk_signal_vec = Vector{Vector{Float64}}(undef,K);
    x_vec =  [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
    # all([all(a .== b) for (a,b) in zip(x_vec,vcat(x_to_use...)) ])
    # all([all(rtik_mat[k:K:end,end] .== [el[k] for el in vcat(rtik...)]) for k in 1:K])
    N_signal_ = [N_signal_mat[(i-1)*G+1:i*G,6] for i in 1:N*K]
    Nj_signal_k = [N_signal_[k:K:end,end] for k in 1:K];

    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    for k in states_id
        Nj_signal_ =  Nj_signal_k[k];
        x_hat_err = broadcast( .*, Nj_signal_, x_vec)
        x_hatk_signal_vec[k] = sum(x_hat_err)#sum(x_hat_) 
    end
    x_hatk_signal_mat[:,1] = states
    x_hatk_signal_mat[:,2] = genes
    x_hatk_signal_mat[:,end] = recursive_flatten(x_hatk_signal_vec)
    return x_hatk_signal_mat
    # all(sum(x_hat_[1:600]) .== x_hat_tk[1][1])
    # all(sum(x_hat_[601:1200]) .== x_hat_tk[1][2])
    # all(sum(x_hat_[1201:1800]) .== x_hat_tk[1][3])
    # all(sum(x_hat_[2401:3000]) .== x_hat_tk[1][5])
end
function depracated_tidy_update_x_hat_sq_error(xmat,N_mat)
    T = length(unique(xmat[:,1]));
    N = length(unique(xmat[:,2]));
    G = length(unique(xmat[:,3]));
    K = 1;
    timepoints = collect(1:T);
    states_id = collect(1:K);
    cell_ids = collect(1:N);
    gene_ids = collect(1:G);
    N_t = tidy_get_Nt_from_xmat(xmat);
    states = innermelt(states_id,G)
    genes = outermelt(gene_ids,K)
    nrows = G; 
    ncols = 3;
    @views N_error_mat = N_mat[2:2:end,:]
    N_t = tidy_get_Nt_from_rtikmat(xmat)
    x_hat_sq_err_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols);
    x_hat_sq_err_vec = Vector{Vector{Float64}}(undef,K);
    x_sq_vec =  [(xmat[(i-1)*G+1:(i)*G,end]) .^2 for i in 1:N];
    # all([all(a .== b) for (a,b) in zip(x_vec,vcat(x_to_use...)) ])
    # all([all(rtik_mat[k:K:end,end] .== [el[k] for el in vcat(rtik...)]) for k in 1:K])
    # [sum([N_signal_mat[(i-1)*K*G+1:K*G*i,:][(k-1)*G+1:k*G,end] for i in 1:N]) for k in 1:K]
    N_error_ = [N_error_mat[(i-1)*G+1:i*G,6] for i in 1:N*K]
    Nj_error_k = [N_error_[k:K:end,end] for k in 1:K];

    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    
    for k in states_id
        Nj_error_ =  Nj_error_k[k];
        x_hat_sq_err = broadcast( .*, Nj_error_, x_sq_vec)
        x_hat_sq_err_vec[k] = sum(x_hat_sq_err)#sum([sum(x_hat_err[st:en])  for (st,en) in timeranges])#sum(x_hat_) 
    end
    x_hat_sq_err_vec = sum(x_hat_sq_err_vec)
    x_hat_sq_err_mat[:,1] = states
    x_hat_sq_err_mat[:,2] = genes
    x_hat_sq_err_mat[:,end] = recursive_flatten(x_hat_sq_err_vec)
    return x_hat_sq_err_mat
end
function deprecated_tidy_update_ηtkj(rtik_mat,e_log_τ_kj_mat,e_τ_μ_tikj_mat,e_log_τj_err_mat,e_τ_μ_tij_err_mat,ηtkj_prior_mat)
    T = length(unique(e_τ_μ_tikj_mat[:,1]))
    K = length(unique(e_log_τ_kj_mat[:,1]))
    G = length(unique(e_log_τ_kj_mat[:,2]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
    Glog = G*log(2π)
    logpi= Glog/G
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
    gene_ids = collect(1:G)
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)


    rtik_cell = innermelt(rtik_mat[:,end],G)
    e_log_τ_kj_cell = outermelt(e_log_τ_kj_mat[:,end],N)
    e_τ_μ_tikj_cell = e_τ_μ_tikj_mat[:,end]
    cell_signal_sq = rtik_cell .* (0.5 .* (e_log_τ_kj_cell .- logpi .-  e_τ_μ_tikj_cell))
    cell_signal_sq_lin = [cell_signal_sq[(i-1)*K*G+1:(i)*K*G,end] for i in 1:N]
    cell_signal_imp = [sum(cell_signal_sq_lin[st:en]) for (st,en) in timeranges] .+ [log.(ηtkj_prior_mat[1:2:end,end])[(t-1)*K*G+1:(t)*K*G] - log.(ηtkj_prior_mat[2:2:end,end])[(t-1)*K*G+1:(t)*K*G] for t in 1:T]
    # cell_signal_imp = recursive_flatten(cell_signal_imp)
    # e_log_τj_err_cell = outermelt(e_log_τj_err_mat[:,end],N*K)
    # e_τ_μ_tij_err_cell = recursive_flatten(outermelt.([e_τ_μ_tij_err_mat[(i-1)*G+1:i*G,end] for i in 1:N],K))
    # cell_err_sq = rtik_cell .* (0.5 .* (e_log_τj_err_cell .- logpi .-  e_τ_μ_tij_err_cell))
    # cell_err_sq_lin = [cell_err_sq[(i-1)*K*G+1:(i)*K*G,end] for i in 1:N ]
    # cell_signal_notimp = [sum(cell_err_sq_lin[st:en]) for (st,en) in timeranges] .+ [log.(ηtkj_prior_mat[2:2:end,end])[(t-1)*K*G+1:(t)*K*G] for t in 1:T]
    # cell_signal_notimp = recursive_flatten(cell_signal_notimp)
    log_ηtkj_tilde = recursive_flatten(cell_signal_imp) #[[imp,notimp] for (imp,notimp) in zip(cell_signal_imp,cell_signal_notimp)]
    ηtkj_imp = logistic.(log_ηtkj_tilde)#[norm_weights(el) for el in log_ηtkj_tilde]
    ηtkj_notimp = 1 .- ηtkj_imp
    ηtkj_ = [[imp,notimp] for (imp,notimp) in zip(ηtkj_imp,ηtkj_notimp)]
    ηtkj_mat = reduce(hcat,[ηtkj_prior_mat[:,1:4],recursive_flatten(ηtkj_) ])
    log_ηtkj_tilde_mat = reduce(hcat,[ηtkj_prior_mat[1:2:end,1:4],recursive_flatten(ηtkj_) ])
    return ηtkj_mat,log_ηtkj_tilde_mat
end
#####################################################
