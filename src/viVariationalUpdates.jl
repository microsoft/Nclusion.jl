"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
function update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams;mt_mode = nothing)
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

"""
"""
function update_d_hat_sum_mpu!(conditionparams,dataparams;mt_mode = nothing)
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

"""
"""
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

"""
"""
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

"""
"""
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
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,ab_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter))
    for k in 1:K
        clusters[k].gk_hat[1] = sigmoid(lb_lg_results.res.minimizer[1,k])#new_rho_hat[k]
        clusters[k].hk_hat[1] = exp(lb_lg_results.res.minimizer[2,k])#new_omega_hat[k]
        clusters[k].ak_hat[1] = StatsFuns.logit(clusters[k].gk_hat[1]) #new_c_hat[k]
        clusters[k].bk_hat[1] = log(clusters[k].hk_hat[1])#new_d_hat[k]
    end
    return clusters
end

"""
"""
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

"""
"""
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

"""
"""
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

"""
"""
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
