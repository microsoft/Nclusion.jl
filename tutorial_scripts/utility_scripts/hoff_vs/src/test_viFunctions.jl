
# variational_inference_dynamicHDP_vs,variational_inference_dynamicHDP_vs1,variational_inference_dynamicHDP_vs2,variational_inference_dynamicHDP_vs3,variational_inference_dynamicHDP_vs3,variational_inference_dynamicHDP_vs4,variational_inference_dynamicHDP_vs5,variational_inference_dynamicHDP_vs6,variational_inference_dynamicHDP_vs6_testFixedClusters,variational_inference_dynamicHDP_vs7,variational_inference_dynamicHDP_vs7_testFixedClusters
#####################################
function update_η_tikj(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,2) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(η_prior[t][i][k][j]) 
                    log_η_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - η_prior[t][i][k][j])
                    log_η_tikj = norm_weights(log_η_tikj_tilde)
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
function update_η_tikj_3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,2) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(η_prior[t][i][k][j]) 
                    log_η_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - η_prior[t][i][k][j])
                    log_η_tikj = norm_weights(log_η_tikj_tilde)
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
#####################################
function update_η_tikj1(Glog,rtik,e_log_τkj,e_τ_μ_tikj,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi = Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,2) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(η_prior[t][i][k][j]) 
                    log_η_tikj_tilde[2] = log(1 - η_prior[t][i][k][j])#rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * Glog - 0.5 * e_τ_μj_err[t][i][j]) + 
                    log_η_tikj = norm_weights(log_η_tikj_tilde)
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
function update_η_tikj1_2(Glog,rtik,e_log_τkj,e_τ_μ_tikj,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi = Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,2) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(η_prior[t][i][k][j]) 
                    log_η_tikj_tilde[2] = log(1 - η_prior[t][i][k][j])#rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * Glog - 0.5 * e_τ_μj_err[t][i][j]) + 
                    log_η_tikj = norm_weights(log_η_tikj_tilde)
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
#####################################

function update_η_tikj2(Glog,rtik,e_log_τ,e_τ_μ,e_log_τj_err,e_τ_μ_tij_err,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_k = Vector{Vector{Float64}}(undef,K)
    η_k_tilde = Vector{Vector{Float64}}(undef,K)
    for  k in 1:K
        η_tilde = Vector{Float64}(undef,2)
        η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
        η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
        η_k_tilde[k] = η_tilde
        η_k[k] = norm_weights(η_tilde);
    end

    return η_k,η_k_tilde
end
function update_η_tikj2_2(Glog,rtik,e_log_τ,e_τ_μ,e_log_τj_err,e_τ_μ_tij_err,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_k = Vector{Vector{Float64}}(undef,K)
    η_k_tilde = Vector{Vector{Float64}}(undef,K)
    for  k in 1:K
        η_tilde = Vector{Float64}(undef,2)
        η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
        η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
        η_k_tilde[k] = η_tilde
        η_k[k] = norm_weights(η_tilde);
    end

    return η_k,η_k_tilde
end
#####################################
function update_η_tikj3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    # η_k = Vector{Vector{Float64}}(undef,K)
    # η_k_tilde = Vector{Vector{Float64}}(undef,K)
    # for  k in 1:K
    #     η_tilde = Vector{Float64}(undef,2)
    #     η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
    #     η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
    #     η_k_tilde[k] = η_tilde
    #     η_k[k] = norm_weights(η_tilde);
    # end


    ###############################
    # [[[sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]
    # [[[sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]

    # [[[sum([rtik[t][i][k] for i in  1:C_t[t]]) for k in 1:K] for t in 1:T] for j in 1:G]

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,2)
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) 
                log_η_tkj_tilde[2] = sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j])
                log_η_tkj = norm_weights(log_η_tkj_tilde)
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end
function update_η_tikj3_2(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    # η_k = Vector{Vector{Float64}}(undef,K)
    # η_k_tilde = Vector{Vector{Float64}}(undef,K)
    # for  k in 1:K
    #     η_tilde = Vector{Float64}(undef,2)
    #     η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
    #     η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
    #     η_k_tilde[k] = η_tilde
    #     η_k[k] = norm_weights(η_tilde);
    # end


    ###############################
    # [[[sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]
    # [[[sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]

    # [[[sum([rtik[t][i][k] for i in  1:C_t[t]]) for k in 1:K] for t in 1:T] for j in 1:G]

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,2)
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) 
                log_η_tkj_tilde[2] = sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j])
                log_η_tkj = norm_weights(log_η_tkj_tilde)
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end

function update_η_tikj3_3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,2)
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j][1]) 
                log_η_tkj_tilde[2] = sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(ηprior[t][k][j][2])
                log_η_tkj = norm_weights(log_η_tkj_tilde)
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end

function update_η_tkj_sigmoid6(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    # η_k = Vector{Vector{Float64}}(undef,K)
    # η_k_tilde = Vector{Vector{Float64}}(undef,K)
    # for  k in 1:K
    #     η_tilde = Vector{Float64}(undef,2)
    #     η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
    #     η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
    #     η_k_tilde[k] = η_tilde
    #     η_k[k] = norm_weights(η_tilde);
    # end


    ###############################
    # [[[sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]
    # [[[sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]

    # [[[sum([rtik[t][i][k] for i in  1:C_t[t]]) for k in 1:K] for t in 1:T] for j in 1:G]

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,1)
                # log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j] /(1 - ηprior[t][k][j])) 
                # log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j] /(1 - ηprior[t][k][j])) 
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * e_τ_μ_tikj[t][i][k][j] - 0.5 * e_log_τj_err[j] + 0.5 * e_τ_μ_tij_err[t][i][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j] /(1 - ηprior[t][k][j])) 
                # log_η_tkj_tilde[2] = sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j])
                log_η_tkj[1] = logistic.(log_η_tkj_tilde)
                log_η_tkj[2] = 1 .- log_η_tkj[1]
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end


#####################################
function update_η_tikj4(Glog,rtik,e_log_τkj,e_τ_μ_tikj,spike_logpdf,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    # η_k = Vector{Vector{Float64}}(undef,K)
    # η_k_tilde = Vector{Vector{Float64}}(undef,K)
    # for  k in 1:K
    #     η_tilde = Vector{Float64}(undef,2)
    #     η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
    #     η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
    #     η_k_tilde[k] = η_tilde
    #     η_k[k] = norm_weights(η_tilde);
    # end

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,2)
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) 
                log_η_tkj_tilde[2] = sum([ rtik[t][i][k] .* spike_logpdf[t][i][j] for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j])
                # if t ==1 && k == 1
                #     println(log_η_tkj_tilde)
                # end
                log_η_tkj = norm_weights(log_η_tkj_tilde)
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end
function update_η_tikj4_2(Glog,rtik,e_log_τkj,e_τ_μ_tikj,spike_logpdf,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    # η_k = Vector{Vector{Float64}}(undef,K)
    # η_k_tilde = Vector{Vector{Float64}}(undef,K)
    # for  k in 1:K
    #     η_tilde = Vector{Float64}(undef,2)
    #     η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
    #     η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
    #     η_k_tilde[k] = η_tilde
    #     η_k[k] = norm_weights(η_tilde);
    # end

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,2)
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) 
                log_η_tkj_tilde[2] = sum([ rtik[t][i][k] .* spike_logpdf[t][i][j] for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j])
                # if t ==1 && k == 1
                #     println(log_η_tkj_tilde)
                # end
                log_η_tkj = norm_weights(log_η_tkj_tilde)
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end
#####################################
function update_η_tikj5(Glog,rtik,e_log_τkj,e_τ_μ_tikj,spike_logpdf,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,2) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(η_prior[t][i][k][j]) 
                    log_η_tikj_tilde[2] = rtik[t][i][k] *  (spike_logpdf[t][i][j]) + log(1 - η_prior[t][i][k][j])
                    log_η_tikj = norm_weights(log_η_tikj_tilde)
                    # if t ==1 && k == 1 && i == cells_
                    #     println(log_η_tikj_tilde)
                    # end
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
function update_η_tikj5(Glog,rtik,e_log_τkj,e_τ_μ_tikj,spike_logpdf,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,2) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(η_prior[t][i][k][j]) 
                    log_η_tikj_tilde[2] = rtik[t][i][k] *  (spike_logpdf[t][i][j]) + log(1 - η_prior[t][i][k][j])
                    log_η_tikj = norm_weights(log_η_tikj_tilde)
                    # if t ==1 && k == 1 && i == cells_
                    #     println(log_η_tikj_tilde)
                    # end
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
#####################################

function update_η_tikj7sigmoid(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    η_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    η_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        η_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        η_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_η_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_η_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_η_tik = Vector{Vector{Float64}}(undef,G)
                log_η_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_η_tikj = Vector{Float64}(undef,2)
                    log_η_tikj_tilde = Vector{Float64}(undef,1) 
                    log_η_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) - rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j])  + log(η_prior[t][i][k][j]) - log(1 - η_prior[t][i][k][j])
                    log_η_tikj = [logistic(log_η_tikj_tilde[1]), 1- logistic(log_η_tikj_tilde[1]) ] ##norm_weights(log_η_tikj_tilde)
                    log_η_tik_tilde[j] = log_η_tikj_tilde
                    log_η_tik[j] = log_η_tikj
                    # println(" not broke")
                end
                log_η_ti[k] = log_η_tik
                log_η_ti_tilde[k] = log_η_tik_tilde
            end
            η_t[i] = log_η_ti
            η_t_tilde[i] = log_η_ti_tilde
        end
        η_tikj[t] = η_t
        η_tikj_tilde[t] = η_t_tilde
    end
    return η_tikj,η_tikj_tilde
end
######################################

# function τ_μ_error_expected_value(x,λ0_err_vec,m_err_vec,a0_err_vec, b0_err_vec)
#     T = length(x)
#     C_t = [length(el) for el in x]
#     G = length(a0_err_vec)
#     e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
#     e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
#     for t in 1:T
#         cells = C_t[t]
#         e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
#         e_τ_μ_13 =  Vector{Float64}(undef,cells)
#         for i in 1:cells
#             e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i] .- m_err_vec) .^2 .+ 1 ./λ0_err_vec
#             e_τ_μ_23 =  sum(e_τ_μ_kjti3)
#             e_τ_μ_kjt3[i] = e_τ_μ_kjti3
#             e_τ_μ_13[i] = e_τ_μ_23
#         end
#         e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
#         e_τ_μ_true3[t] =e_τ_μ_13
#     end

#     return e_τ_μ_kj_true3,e_τ_μ_true3
# end
# function τ_μ_error_expected_value(x,λ0_err_vec,m_err_vec,a0_err_vec, b0_err_vec)
#     T = length(x)
#     C_t = [length(el) for el in x]
#     G = length(a0_err_vec)
#     e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
#     e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
#     for t in 1:T
#         cells = C_t[t]
#         e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
#         e_τ_μ_13 =  Vector{Float64}(undef,cells)
#         for i in 1:cells
#             e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i] .- m_err_vec) .^2 .+ 1 ./λ0_err_vec
#             e_τ_μ_23 =  sum(e_τ_μ_kjti3)
#             e_τ_μ_kjt3[i] = e_τ_μ_kjti3
#             e_τ_μ_13[i] = e_τ_μ_23
#         end
#         e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
#         e_τ_μ_true3[t] =e_τ_μ_13
#     end

#     return e_τ_μ_kj_true3,e_τ_μ_true3
# end

function initialize_η_prior(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[[pct_important for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return η_prior
end
function initialize_η_prior2(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [pct_important for k in 1:K] 
    return η_prior
end
function initialize_η_prior3(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[pct_important for j in 1:G] for k in 1:K] for t in 1:T] 
    return η_prior
end
function initialize_η_prior4(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[pct_important for j in 1:G] for k in 1:K] for t in 1:T] 
    return η_prior
end
function initialize_η_prior5(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[[pct_important for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return η_prior
end

# function log_τ_k_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
#     e_log_τ_kj_vec = log_τ_kj_expected_value(a0_err_hat_vec, b0_err_hat_vec)
#     e_log_τ_k = sum(e_log_τ_kj_vec)
#     return e_log_τ_k
# end
# function log_τ_k_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
#     e_log_τ_kj_vec = log_τ_kj_expected_value(a0_err_hat_vec, b0_err_hat_vec)
#     e_log_τ_k = sum(e_log_τ_kj_vec)
#     return e_log_τ_k
# end


# function log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
#     e_log_τ_kj_vec = digamma.(a0_err_hat_vec) .- log.(b0_err_hat_vec)
    
#     return e_log_τ_kj_vec
# end
# function log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
#     e_log_τ_kj_vec = digamma.(a0_err_hat_vec) .- log.(b0_err_hat_vec)
    
#     return e_log_τ_kj_vec
# end
function update_N(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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

function update_N_forloops(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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

function update_N_forloops2(rtik,η_k)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    N_signal = Vector{Vector{Vector{Float64}}}(undef,T)
    N_error = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Float64}}(undef,cells_)
        Nt_error =  Vector{Vector{Float64}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Float64}(undef,K)
            Nti_error = Vector{Float64}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] * η_k[k][1]
                Nti_error[k] = rtik[t][i][k] * η_k[k][2]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_N_forloops3(rtik,η_tkj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tkj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal =Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_tkj[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_tkj[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_N_forloops3(rtik,η_jkt)
    T = length(rtik)
    K = length(rtik[1][1])
    G = length(η_jkt[1][1])
    C_t = [length(el) for el in rtik]
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_jkt[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_jkt[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_N_forloops_sigmoid(rtik,η_tkj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tkj[1][1])
    N = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt =Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti[k] = rtik[t][i][k] .* [η_tkj[t][k][j][1] for j in 1:G]
            end
            Nt[i] = Nti
        end
        N[t] = Nt
    end
    return N
end
function update_N_forloops4(rtik,η_jkt)
    T = length(rtik)
    K = length(rtik[1][1])
    G = length(η_jkt[1][1])
    C_t = [length(el) for el in rtik]
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_jkt[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_jkt[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end

function update_N_forloops5(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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
function update_N_broadcast(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)

    #main idea broadcast(*,vcat(rtik...)[1] ,vcat(η_tikj...)[1] )
end 
function update_errorNj(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_errorNj(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_errorNj(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_errorNj(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_errorNj_forloops(N_error) # SLOWER and MORE MEMORY when compared to @benchmark  update_errorNj(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    Nj_error = Vector{Float64}(undef, G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                  push!(value,N_error[t][i][k][j])  
                end
            end
        end
        Nj_error[j] = sum(value)
    end
    return Nj_error
end
function update_errorNj_forloops(N_error) # SLOWER and MORE MEMORY when compared to @benchmark  update_errorNj(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    Nj_error = Vector{Float64}(undef, G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                  push!(value,N_error[t][i][k][j])  
                end
            end
        end
        Nj_error[j] = sum(value)
    end
    return Nj_error
end
function update_errorNj_forloops(N_error) # SLOWER and MORE MEMORY when compared to @benchmark  update_errorNj(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    Nj_error = Vector{Float64}(undef, G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                  push!(value,N_error[t][i][k][j])  
                end
            end
        end
        Nj_error[j] = sum(value)
    end
    return Nj_error
end
function update_signalNkj(N_signal)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_signal = reduce(vcat,N_signal)

    Nkj_signal = sum(perCell_linerize_N_signal)
    return Nkj_signal
end
function update_Nkj(N)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_signal = reduce(vcat,N)

    Nkj = sum(perCell_linerize_N_signal)
    return Nkj
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
function update_x_hatk_signal_vs_forloops22(x,N_signal)
    # with η_k
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
function update_x_hatk_signal_vs_forloops3(x,N_signal)
    # with η_k
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
function update_x_hatk_sigmoid(x,N)
    # with η_k
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N[1][1])
    G = length(x[1][1])
    x_hatk = Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                push!(value, x[t][i] .* N[t][i][k])
            end
        end
        x_hatk[k] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hatk
end
function update_x_hatk_signal_vs_forloops4(x,N_signal)
    # with η_k
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
function update_x_hatk_signal_vs_forloops5(x,N_signal)
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
function update_x_hat_error_vs_forloops2(x,rtik,η_tikj) #IDK WHY THIS IS NOT THE SAME AS THE OTHER TWO VERSIONS OTHER THAN THE FACT THAT I DONT PREALLOCATE????
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
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

function update_x_hat_error_vs_forloops(x,rtik,η_tikj)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = zeros(sum(C_t)*K)
        counter = 0
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    counter+=1 
                    value[Int(counter)] =  x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2]
                    # push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
                end
            end
        end
        # println(length(value))
        x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_err
end
function update_x_hat_error_vs_forloops22(x,rtik,η_k)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = zeros(sum(C_t)*K)
        counter = 0
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    counter+=1 
                    value[Int(counter)] =  x[t][i][j] * rtik[t][i][k] * η_k[k][2]
                    # push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
                end
            end
        end
        # println(length(value))
        x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_err
end

function update_x_hat_error_vs(x,N_error)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    # x_hat_err = Vector{Float64}(undef,G)
    # perCell_rtik = reduce(vcat,rtik)
    # perCell_η_tikj = reduce(vcat,η_tikj)
    # perCell_rη = [ r .* η for (r,η) in zip(perCell_rtik,perCell_η_tikj)]
    # perCell_rηx = [[broadcast(*, vv,el) for el in bb] for  (vv,bb) in zip(vcat(x...),perCell_rη)]
    perCell_x = reduce(vcat,x)
    perCell_N_error = reduce(vcat,N_error)
    perCell_xNerror = [[broadcast(*, j,el) for el in k] for  (j,k) in zip(perCell_x,perCell_N_error)]
    x_hat_err = sum(vcat(perCell_xNerror...))
    return x_hat_err
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
function update_x_hatk_sq_signal_vs_forloops2(x,N_signal)
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
function update_x_hatk_sq_signal_vs_forloops3(x,N_signal)
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

function update_x_hatk_sq_sigmoid(x,N)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N[1][1])
    G = length(x[1][1])
    x_hatk_sq= Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                push!(value, x[t][i] .^2 .* N[t][i][k])
            end
        end
        x_hatk_sq[k] = sum(value)
    end
    return x_hatk_sq
end

function update_x_hatk_sq_signal_vs_forloops4(x,N_signal)
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
function update_x_hatk_sq_signal_vs_forloops5(x,N_signal)
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
function update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error .+1)
    return a0_err_hat_vec
end
function update_λ0_err_hat(λ0_err_vec,Nj_error) 
    λ0_err_hat_vec = λ0_err_vec .+ Nj_error 
    return λ0_err_hat_vec
end
function update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_err_vec .* μ0_err_vec
    denom = λ0_err_vec .+  Nj_error 
    m_err_hat_vec = (λ0_μ0 .+ x_hat_err) ./ denom #[ (λ0_μ0 .+ x_hat_k[k]) ./denom[k] for k in 1:K]
    return m_err_hat_vec
end

function update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)
    denom = λ0_err_vec .+  Nj_error 
    μ0_sq_vec = μ0_err_vec .^2
    μ0λ0_vec =  λ0_err_vec .* μ0_err_vec
    μ0_sq_λ0_vec = λ0_err_vec .* μ0_sq_vec
    numer = (x_hat_err.- μ0λ0_vec) .^2 
    ssd = numer ./ denom
    half_sk_ssd =  1/2 .* (x_hat_sq_err .+ μ0_sq_λ0_vec .- ssd)
    b0_err_hat_vec = b0_err_vec .+ half_sk_ssd
    return  b0_err_hat_vec
end


# function update_λ0k_signal_hat(λ0_vec,Nkj_signal) 
#     K = length(Nkj_signal)
#     λ0k_hat_vec = [λ0_vec .+ Nkj_signal[k] for k in 1:K]
#     return λ0k_hat_vec
# end
function update_λ0k_sigmoid(λ0_vec,Nkj) 
    K = length(Nkj)
    λ0k_hat_vec = [λ0_vec .+ Nkj[k] for k in 1:K]
    return λ0k_hat_vec
end


# function update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
#     K = length(Nkj_signal) 
#     a0k_hat_vec = [ a0_vec .+ 1/2 * (Nkj_signal[k] .+ 1) for k in 1:K]
#     return a0k_hat_vec
# end

function update_a0k_sigmoid(a0_vec,Nkj)
    K = length(Nkj) 
    a0k_hat_vec = [ a0_vec .+ 1/2 * (Nkj[k] .+ 1) for k in 1:K]
    return a0k_hat_vec
end

# function update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
#     K = length(Nkj_signal)
#     # Nk_xbar_k = Nk .* xbar_k
#     λ0_μ0 =  λ0_vec .* μ0_vec
#     denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
#     mkj_hat = [ (λ0_μ0 .+ x_hatk_signal[k]) ./denom[k] for k in 1:K]
#     return mkj_hat
# end
function update_mk_sigmoid(λ0_vec,μ0_vec, Nkj,x_hatk)
    K = length(Nkj)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_vec .* μ0_vec
    denom = [λ0_vec .+ Nkj[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    mkj_hat = [ (λ0_μ0 .+ x_hatk[k]) ./denom[k] for k in 1:K]
    return mkj_hat
end

# function update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
#     K = length(Nkj_signal)
#     denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
#     μ0_sq_vec = μ0_vec .^2
#     μ0λ0_vec =  λ0_vec .* μ0_vec
#     μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
#     numer = [(x_hatk_signal[k] .- μ0λ0_vec) .^2 for k in 1:K ]
#     ssd = [numer[k] ./ denom[k] for k in 1:K]
#     half_sk_ssd =  1/2 .* [x_hatk_sq_signal[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
#     # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
#     b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
#     # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
#     # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
#     return  b0k_hat_vec
# end

function update_b0k_sigmoid(b0_vec,λ0_vec,μ0_vec, Nkj,x_hatk,x_hatk_sq)
    K = length(Nkj)
    denom = [λ0_vec .+ Nkj[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    μ0_sq_vec = μ0_vec .^2
    μ0λ0_vec =  λ0_vec .* μ0_vec
    μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    numer = [(x_hatk[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    ssd = [numer[k] ./ denom[k] for k in 1:K]
    half_sk_ssd =  1/2 .* [x_hatk_sq[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    return  b0k_hat_vec
end

function update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][i][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][i][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  logpi .-e_τ_μ_tij_err[t][i]))
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
function update_rtik_vs1(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][i][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][i][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  Glog .- e_τ_μ_tikj[t][i][k]))
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

function update_rtik_vs2(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_k,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = η_k[k][1]
                η_false = η_k[k][2]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  Glog .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  Glog .-e_τ_μ_tij_err[t][i]))
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

function update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tkj[t][k][j][1] for j in 1:G]
                η_false = [η_tkj[t][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  logpi .-e_τ_μ_tij_err[t][i]))
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
function update_rtik_vs3_sigmoid(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tkj[t][k][j][1] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .- e_τ_μ_tikj[t][i][k] .- e_log_τj_err .+ e_τ_μ_tij_err[t][i]) .+ (e_log_τj_err .-  logpi .-e_τ_μ_tij_err[t][i]))
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


function update_rtik_vs4(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,spike_logpdf,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* spike_logpdf[t][i])
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
function update_rtik_vs5(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, spike_logpdf,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][i][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][i][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (spike_logpdf[t][i]))
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

function update_rtik_vs7(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                # η_true = [η_tikj[t][i][k][j][1] for j in 1:G]
                # η_false = [η_tikj[t][i][k][j][2] for j in 1:G]
                # ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  logpi .-e_τ_μ_tij_err[t][i])) #
                log_like_gene_vec = Vector{Float64}(undef,G)
                for j in 1:G
                    η_true = η_tikj[t][i][k][j][1]
                    η_false = η_tikj[t][i][k][j][2]
                    log_like_gene = 0.5 * η_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j]) + 0.5 * η_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    # log_like_gene = 0.5 * η_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) 
                    # # + 0.5 * η_false * ()
                    log_like_gene_vec[j] =  log_like_gene
                end
                sum_log_like_gene = sum(log_like_gene_vec)
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] + sum_log_like_gene
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

function variational_inference_dynamicHDP_vs(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,η_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    
    


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = Dict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["xbar_k_beforeNorm"]= []
        debug_val["xbar_k_afterNorm"]= []
        debug_val["sk_beforeNorm"]= []
        debug_val["sk_afterNorm"]= []
        debug_val["Tk"]= []
        debug_val["data_elbo"]= []
        debug_val["assgn_entropy"]= []
        debug_val["HDP_surragate_elbo"]= []
        
    end
    #init debug dict initial values
    if debugme
        push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
        push!(debug_val["mk_hat_vec"],mk_hat_vec)
        push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
        push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat)
        push!(debug_val["rtik"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["xbar_k_beforeNorm"],[])
        push!(debug_val["xbar_k_afterNorm"],[])
        push!(debug_val["sk_beforeNorm"],[])
        push!(debug_val["sk_afterNorm"],[])
        push!(debug_val["Tk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x_to_use,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tikj,c_ttprime_vec) #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            η_tikj,_  = update_η_tikj(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_prior);

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
            if debugme
                push!(debug_val["θ_hat"],θ_hat)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
            end
        end

        Nk = update_Nk(rtik)
        if debugme
            push!(debug_val["Nk"],Nk)
        end

        x_hat_k = update_x_hat_k(x,rtik)
        if debugme
            push!(debug_val["xbar_k_beforeNorm"],xbar_k)
        end
        
        # xbar_k = 1 ./ Nk .* xbar_k
        if debugme
            push!(debug_val["xbar_k_afterNorm"],xbar_k)
        end
        
        x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        if debugme
            push!(debug_val["sk_beforeNorm"],sk)
        end
        
        # sk = 1 ./ Nk .* sk
        if debugme
            push!(debug_val["sk_afterNorm"],sk)
        end

        N_signal,N_error = update_N(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x_to_use,N_signal)


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x_to_use,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)
        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["Tk"],Tk)
            push!(debug_val["data_elbo"],data_elbo)
            push!(debug_val["assgn_entropy"],assgn_entropy)
            push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,η_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    
    


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = Dict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["xbar_k_beforeNorm"]= []
        debug_val["xbar_k_afterNorm"]= []
        debug_val["sk_beforeNorm"]= []
        debug_val["sk_afterNorm"]= []
        debug_val["Tk"]= []
        debug_val["data_elbo"]= []
        debug_val["assgn_entropy"]= []
        debug_val["HDP_surragate_elbo"]= []
        
    end
    #init debug dict initial values
    if debugme
        push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
        push!(debug_val["mk_hat_vec"],mk_hat_vec)
        push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
        push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat)
        push!(debug_val["rtik"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["xbar_k_beforeNorm"],[])
        push!(debug_val["xbar_k_afterNorm"],[])
        push!(debug_val["sk_beforeNorm"],[])
        push!(debug_val["sk_afterNorm"],[])
        push!(debug_val["Tk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            η_tikj,_  = update_η_tikj(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior);

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
            if debugme
                push!(debug_val["θ_hat"],θ_hat)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
            end
        end

        Nk = update_Nk(rtik)
        if debugme
            push!(debug_val["Nk"],Nk)
        end

        x_hat_k = update_x_hat_k(x,rtik)
        if debugme
            push!(debug_val["xbar_k_beforeNorm"],xbar_k)
        end
        
        # xbar_k = 1 ./ Nk .* xbar_k
        if debugme
            push!(debug_val["xbar_k_afterNorm"],xbar_k)
        end
        
        x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        if debugme
            push!(debug_val["sk_beforeNorm"],sk)
        end
        
        # sk = 1 ./ Nk .* sk
        if debugme
            push!(debug_val["sk_afterNorm"],sk)
        end

        N_signal,N_error = update_N(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)
        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["Tk"],Tk)
            push!(debug_val["data_elbo"],data_elbo)
            push!(debug_val["assgn_entropy"],assgn_entropy)
            push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

function variational_inference_dynamicHDP_vs1(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,η_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            η_tikj,_  = update_η_tikj(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior);

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
            if debugme
                push!(debug_val["θ_hat"],θ_hat)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
            end
        end

        Nk = update_Nk(rtik)
        if debugme
            push!(debug_val["Nk"],Nk)
        end

        x_hat_k = update_x_hat_k(x,rtik)
        if debugme
            push!(debug_val["xbar_k_beforeNorm"],xbar_k)
        end
        
        # xbar_k = 1 ./ Nk .* xbar_k
        if debugme
            push!(debug_val["xbar_k_afterNorm"],xbar_k)
        end
        
        x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        if debugme
            push!(debug_val["sk_beforeNorm"],sk)
        end
        
        # sk = 1 ./ Nk .* sk
        if debugme
            push!(debug_val["sk_afterNorm"],sk)
        end

        N_signal,N_error = update_N(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs2(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_k_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_k_vec_init) && rand_init
        η_k_vec_init = [rand(Dirichlet(ones(2) ./2)) for k in 1:K]
    elseif isnothing(η_k_vec_init) && !rand_init
        η_k_vec_init = [ones(2) ./2  for k in 1:K]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_k = η_k_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs2(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_k,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            η_k,_  = update_η_tikj2(Glog,rtik,e_log_τ,e_τ_μ,e_log_τj_err,e_τ_μ_tij_err,ηprior);
            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        
        N_signal,N_error = update_N_forloops2(rtik,η_k);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops22(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error) update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops2(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_k,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end


function variational_inference_dynamicHDP_vs3(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            η_tkj,_  = update_η_tikj3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);


            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        

        N_signal,N_error = update_N_forloops3(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops3(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error) update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops3(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

function variational_inference_dynamicHDP_vs4(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,σ²_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    spike_logpdf = [[[logpdf(Normal(0,sqrt(σ²_err)),x[t][i][j]) for j in 1:G] for i in 1:C_t[t]] for t in 1:T]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            # e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            # e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs4(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,spike_logpdf,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            η_tkj,_  = update_η_tikj4(Glog,rtik,e_log_τkj,e_τ_μ_tikj,spike_logpdf,ηprior);

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        
        

        N_signal,N_error = update_N_forloops4(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops4(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error) update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops4(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs5(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,σ²_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    spike_logpdf = [[[logpdf(Normal(0,sqrt(σ²_err)),x[t][i][j]) for j in 1:G] for i in 1:C_t[t]] for t in 1:T]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            # e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            # e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs5(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,spike_logpdf,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            η_tikj,_  = update_η_tikj5(Glog,rtik,e_log_τkj,e_τ_μ_tikj,spike_logpdf,ηprior);

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        
        

        N_signal,N_error = update_N_forloops5(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops5(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error) update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops5(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

### THE SAME A Version 3 but we update error Distributions

function variational_inference_dynamicHDP_vs6(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tkj,_  = update_η_tikj3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);
            rtik = update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        
        

        N_signal,N_error = update_N_forloops3(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops3(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error) #update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops3(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

function variational_inference_dynamicHDP_vs6_testFixedClusters(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true,rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    
    

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tkj,_  = update_η_tikj3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);
            for t in 1:T
                η_tkj[t][1] = [[0.0,1.0] for j in 1:G]
                η_tkj[t][2] = [ j<=3 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
                η_tkj[t][3] = [ G-j<3 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
            end
            rtik = update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        
        

        N_signal,N_error = update_N_forloops3(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops3(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error) #update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops3(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

function variational_inference_dynamicHDP_vs6sigmoid(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            # e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            # e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tkj,_  = update_η_tkj_sigmoid6(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);
            rtik = update_rtik_vs3_sigmoid(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end
        # println(e_τ_μ_tikj[1][1])
        # println(e_log_τkj)
        
        
        

        N_signal,N_error = update_N_forloops3(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)

        

    


        # x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops3(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error) #update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops3(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_sigmoid(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_sigmoid(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_sigmoid(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_sigmoid(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end


###### THE SAME AS version VS but we update the Distributions

function variational_inference_dynamicHDP_vs7(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,η_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    


    # mk_hat_vec = mk_hat_vec_init 
    # λ0k_hat_vec = λ0k_hat_vec_init
    # a0k_hat_vec = a0k_hat_vec_init
    # b0k_hat_vec = b0k_hat_vec_init
    # rhok_hat_vec = rhok_hat_vec_init
    # omegak_hat_vec = omegak_hat_vec_init
    # a_γ_hat = a_γ_hat_init 
    # b_γ_hat = b_γ_hat_init

     
    
    # awt_hat_vec = awt_hat_vec_init 
    # bwt_hat_vec = bwt_hat_vec_init
    # a_αt_hat_vec = a_αt_hat_vec_init 
    # b_αt_hat_vec = b_αt_hat_vec_init
    # θ_hat_vec = θ_hat_vec_init
    # c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    #init debug dict initial values

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tikj,_  = update_η_tikj(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior);

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs7(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 

            
        end


        
        
        
        # sk = 1 ./ Nk .* sk

        

        N_signal,N_error = update_N(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

function variational_inference_dynamicHDP_vs7sigmoid(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,η_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    


    # mk_hat_vec = mk_hat_vec_init 
    # λ0k_hat_vec = λ0k_hat_vec_init
    # a0k_hat_vec = a0k_hat_vec_init
    # b0k_hat_vec = b0k_hat_vec_init
    # rhok_hat_vec = rhok_hat_vec_init
    # omegak_hat_vec = omegak_hat_vec_init
    # a_γ_hat = a_γ_hat_init 
    # b_γ_hat = b_γ_hat_init

     
    
    # awt_hat_vec = awt_hat_vec_init 
    # bwt_hat_vec = bwt_hat_vec_init
    # a_αt_hat_vec = a_αt_hat_vec_init 
    # b_αt_hat_vec = b_αt_hat_vec_init
    # θ_hat_vec = θ_hat_vec_init
    # c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    #init debug dict initial values

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tikj,_  = update_η_tikj7sigmoid(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior);

            rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 

            
        end


        
        
        
        # sk = 1 ./ Nk .* sk

        

        N_signal,N_error = update_N(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end


function variational_inference_dynamicHDP_vs7_testFixedClusters(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,η_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tikj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tikj_vec_init) && rand_init
        η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(η_tikj_vec_init) && !rand_init
        η_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tikj = η_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    # for t in 1:T
    #     for i in 1:C_t[t]
    #         η_tikj[t][i][1] = [ 4<=j<=6 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
    #         η_tikj[t][i][2] = [ j<=3 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
    #         η_tikj[t][i][3] = [ G-j<3 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
    #     end
    # end
    


    # mk_hat_vec = mk_hat_vec_init 
    # λ0k_hat_vec = λ0k_hat_vec_init
    # a0k_hat_vec = a0k_hat_vec_init
    # b0k_hat_vec = b0k_hat_vec_init
    # rhok_hat_vec = rhok_hat_vec_init
    # omegak_hat_vec = omegak_hat_vec_init
    # a_γ_hat = a_γ_hat_init 
    # b_γ_hat = b_γ_hat_init

     
    
    # awt_hat_vec = awt_hat_vec_init 
    # bwt_hat_vec = bwt_hat_vec_init
    # a_αt_hat_vec = a_αt_hat_vec_init 
    # b_αt_hat_vec = b_αt_hat_vec_init
    # θ_hat_vec = θ_hat_vec_init
    # c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    #init debug dict initial values

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tikj,_  = update_η_tikj(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_prior);
            for t in 1:T
                for i in 1:C_t[t]
                    η_tikj[t][i][1] = [ j<=3 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
                    η_tikj[t][i][2] = [ 4<=j<=6 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
                    η_tikj[t][i][3] = [ G-j<3 ? [1.0,0.0] : [0.0,1.0] for j in 1:G]
                end
                
            end

            rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 

            
        end


        
        
        
        # sk = 1 ./ Nk .* sk

        

        N_signal,N_error = update_N(rtik,η_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end


function variational_inference_dynamicHDP_vs3_sigmoid(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,rtik_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);
    mk_hat_vec_init = init_mk_hat!(mk_hat_vec_init,x,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_vec_init = init_λ0k_hat_vec!(λ0k_hat_vec_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_vec_init = init_a0k_hat_vec!(a0k_hat_vec_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_vec_init = init_b0k_hat_vec!(b0k_hat_vec_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_vec_init,omegak_hat_vec_init = init_ρωk_hat_vec!(rhok_hat_vec_init,omegak_hat_vec_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_vec_init,dk_hat_vec_init = rhok_hat_vec_init,omegak_hat_vec_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    awt_hat_vec_init = init_awt_hat_vec!(awt_hat_vec_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_vec_init =init_bwt_hat_vec!(bwt_hat_vec_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_vec_init = init_a_αt_hat_vec!(a_αt_hat_vec_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_vec_init =  init_b_αt_hat_vec!(b_αt_hat_vec_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_vec_init = init_c_ttprime_hat_vec!(c_ttprime_vec_init,T;rand_init = rand_init);
    θ_hat_vec_init = init_θ_hat_vec!(θ_hat_vec_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_vec_init, omegak_hat_init= omegak_hat_vec_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)


    η_tkj_vec_init = init_η_tkj_vec!(η_tkj_vec_init,G,K,T;rand_init = rand_init)
    

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            # e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            # e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs3_sigmoid(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            η_tkj,_  = update_η_tkj_sigmoid(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);


            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        

        N = update_N_forloops_sigmoid(rtik,η_tkj);
        # Nj_error = update_errorNj(N_error)
        Nkj = update_Nkj(N)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk =  update_x_hatk_sigmoid(x,N) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error) update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq = update_x_hatk_sq_sigmoid(x,N)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_sigmoid(λ0_vec,Nkj)
        a0k_hat_vec = update_a0k_sigmoid(a0_vec,Nkj)
        mk_hat_vec= update_mk_sigmoid(λ0_vec,μ0_vec, Nkj,x_hatk)
        b0k_hat_vec = update_b0k_sigmoid(b0_vec,λ0_vec,μ0_vec, Nkj,x_hatk,x_hatk_sq)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,initDict,debug_val
end

function initialize_η_prior(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[[pct_important for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return η_prior
end


function update_N(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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

function update_N_forloops(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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
function update_N_broadcast(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)

    #main idea broadcast(*,vcat(rtik...)[1] ,vcat(η_tikj...)[1] )
end 

function update_signalNkj(N_signal)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_signal = reduce(vcat,N_signal)

    Nkj_signal = sum(perCell_linerize_N_signal)
    return Nkj_signal
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
function update_x_hat_error_vs_forloops2(x,rtik,η_tikj) #IDK WHY THIS IS NOT THE SAME AS THE OTHER TWO VERSIONS OTHER THAN THE FACT THAT I DONT PREALLOCATE????
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
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

function update_x_hat_error_vs_forloops(x,rtik,η_tikj)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = zeros(sum(C_t)*K)
        counter = 0
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    counter+=1 
                    value[Int(counter)] =  x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2]
                    # push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
                end
            end
        end
        # println(length(value))
        x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_err
end

function update_x_hat_error_vs(x,N_error)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    # x_hat_err = Vector{Float64}(undef,G)
    # perCell_rtik = reduce(vcat,rtik)
    # perCell_η_tikj = reduce(vcat,η_tikj)
    # perCell_rη = [ r .* η for (r,η) in zip(perCell_rtik,perCell_η_tikj)]
    # perCell_rηx = [[broadcast(*, vv,el) for el in bb] for  (vv,bb) in zip(vcat(x...),perCell_rη)]
    perCell_x = reduce(vcat,x)
    perCell_N_error = reduce(vcat,N_error)
    perCell_xNerror = [[broadcast(*, j,el) for el in k] for  (j,k) in zip(perCell_x,perCell_N_error)]
    x_hat_err = sum(vcat(perCell_xNerror...))
    return x_hat_err
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
function update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error .+1)
    return a0_err_hat_vec
end
function update_λ0_err_hat(λ0_err_vec,Nj_error) 
    λ0_err_hat_vec = λ0_err_vec .+ Nj_error 
    return λ0_err_hat_vec
end
function update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_err_vec .* μ0_err_vec
    denom = λ0_err_vec .+  Nj_error 
    m_err_hat_vec = (λ0_μ0 .+ x_hat_err) ./ denom #[ (λ0_μ0 .+ x_hat_k[k]) ./denom[k] for k in 1:K]
    return m_err_hat_vec
end

function update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)
    denom = λ0_err_vec .+  Nj_error 
    μ0_sq_vec = μ0_err_vec .^2
    μ0λ0_vec =  λ0_err_vec .* μ0_err_vec
    μ0_sq_λ0_vec = λ0_err_vec .* μ0_sq_vec
    numer = (x_hat_err.- μ0λ0_vec) .^2 
    ssd = numer ./ denom
    half_sk_ssd =  1/2 .* x_hat_sq_err .+ μ0_sq_λ0_vec .- ssd
    b0_err_hat_vec = b0_err_vec .+ half_sk_ssd
    return  b0_err_hat_vec
end

function update_λ0k_signal_hat(λ0_vec,Nkj_signal) 
    K = length(Nkj_signal)
    λ0k_hat_vec = [λ0_vec .+ Nkj_signal[k] for k in 1:K]
    return λ0k_hat_vec
end

function update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
    K = length(Nkj_signal) 
    a0k_hat_vec = [ a0_vec .+ 1/2 * (Nkj_signal[k] .+ 1) for k in 1:K]
    return a0k_hat_vec
end

function update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
    K = length(Nkj_signal)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_vec .* μ0_vec
    denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    mkj_hat = [ (λ0_μ0 .+ x_hatk_signal[k]) ./denom[k] for k in 1:K]
    return mkj_hat
end


function update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
    K = length(Nkj_signal)
    denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    μ0_sq_vec = μ0_vec .^2
    μ0λ0_vec =  λ0_vec .* μ0_vec
    μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    numer = [(x_hatk_signal[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    ssd = [numer[k] ./ denom[k] for k in 1:K]
    half_sk_ssd =  1/2 .* [x_hatk_sq_signal[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    return  b0k_hat_vec
end

function update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][i][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][i][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  Glog .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  Glog .-e_τ_μ_tij_err[t][i]))
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

##############################################


function initialize_η_prior(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[[pct_important for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return η_prior
end
function initialize_η_prior2(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [pct_important for k in 1:K] 
    return η_prior
end
function initialize_η_prior3(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[pct_important for j in 1:G] for k in 1:K] for t in 1:T] 
    return η_prior
end
function initialize_η_prior4(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[pct_important for j in 1:G] for k in 1:K] for t in 1:T] 
    return η_prior
end
function initialize_η_prior5(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[[[pct_important for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return η_prior
end

function update_N(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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

function update_N_forloops(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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

function update_N_forloops2(rtik,η_k)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    N_signal = Vector{Vector{Vector{Float64}}}(undef,T)
    N_error = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Float64}}(undef,cells_)
        Nt_error =  Vector{Vector{Float64}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Float64}(undef,K)
            Nti_error = Vector{Float64}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] * η_k[k][1]
                Nti_error[k] = rtik[t][i][k] * η_k[k][2]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_N_forloops3(rtik,η_tkj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tkj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal =Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_tkj[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_tkj[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_N_forloops3(rtik,η_jkt)
    T = length(rtik)
    K = length(rtik[1][1])
    G = length(η_jkt[1][1])
    C_t = [length(el) for el in rtik]
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_jkt[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_jkt[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_N_forloops4(rtik,η_jkt)
    T = length(rtik)
    K = length(rtik[1][1])
    G = length(η_jkt[1][1])
    C_t = [length(el) for el in rtik]
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_jkt[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_jkt[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end

function update_N_forloops5(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
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
                    Ntik_signal[j] = rtik[t][i][k] * η_tikj[t][i][k][j][1]
                    Ntik_error[j] = rtik[t][i][k] * η_tikj[t][i][k][j][2]
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
function update_N_broadcast(rtik,η_tikj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tikj[1][1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)

    #main idea broadcast(*,vcat(rtik...)[1] ,vcat(η_tikj...)[1] )
end 

function update_signalNkj(N_signal)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_signal = reduce(vcat,N_signal)

    Nkj_signal = sum(perCell_linerize_N_signal)
    return Nkj_signal
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
function update_x_hatk_signal_vs_forloops22(x,N_signal)
    # with η_k
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
function update_x_hatk_signal_vs_forloops3(x,N_signal)
    # with η_k
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
function update_x_hatk_signal_vs_forloops4(x,N_signal)
    # with η_k
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
function update_x_hatk_signal_vs_forloops5(x,N_signal)
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
function update_x_hat_error_vs_forloops2(x,rtik,η_tikj) #IDK WHY THIS IS NOT THE SAME AS THE OTHER TWO VERSIONS OTHER THAN THE FACT THAT I DONT PREALLOCATE????
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = []
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
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

function update_x_hat_error_vs_forloops(x,rtik,η_tikj)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = zeros(sum(C_t)*K)
        counter = 0
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    counter+=1 
                    value[Int(counter)] =  x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2]
                    # push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
                end
            end
        end
        # println(length(value))
        x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_err
end
function update_x_hat_error_vs_forloops22(x,rtik,η_k)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    for j in 1:G
        value = zeros(sum(C_t)*K)
        counter = 0
        for t in 1:T
            for i in 1:C_t[t]
                for k in 1:K
                    counter+=1 
                    value[Int(counter)] =  x[t][i][j] * rtik[t][i][k] * η_k[k][2]
                    # push!(value, x[t][i][j] * rtik[t][i][k] * η_tikj[t][i][k][j][2])  
                end
            end
        end
        # println(length(value))
        x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_err
end

function update_x_hat_error_vs(x,N_error)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    # x_hat_err = Vector{Float64}(undef,G)
    # perCell_rtik = reduce(vcat,rtik)
    # perCell_η_tikj = reduce(vcat,η_tikj)
    # perCell_rη = [ r .* η for (r,η) in zip(perCell_rtik,perCell_η_tikj)]
    # perCell_rηx = [[broadcast(*, vv,el) for el in bb] for  (vv,bb) in zip(vcat(x...),perCell_rη)]
    perCell_x = reduce(vcat,x)
    perCell_N_error = reduce(vcat,N_error)
    perCell_xNerror = [[broadcast(*, j,el) for el in k] for  (j,k) in zip(perCell_x,perCell_N_error)]
    x_hat_err = sum(vcat(perCell_xNerror...))
    return x_hat_err
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
function update_x_hatk_sq_signal_vs_forloops2(x,N_signal)
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
function update_x_hatk_sq_signal_vs_forloops3(x,N_signal)
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
function update_x_hatk_sq_signal_vs_forloops4(x,N_signal)
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
function update_x_hatk_sq_signal_vs_forloops5(x,N_signal)
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
function update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error .+1)
    return a0_err_hat_vec
end
function update_λ0_err_hat(λ0_err_vec,Nj_error) 
    λ0_err_hat_vec = λ0_err_vec .+ Nj_error 
    return λ0_err_hat_vec
end
function update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_err_vec .* μ0_err_vec
    denom = λ0_err_vec .+  Nj_error 
    m_err_hat_vec = (λ0_μ0 .+ x_hat_err) ./ denom #[ (λ0_μ0 .+ x_hat_k[k]) ./denom[k] for k in 1:K]
    return m_err_hat_vec
end

function update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)
    denom = λ0_err_vec .+  Nj_error 
    μ0_sq_vec = μ0_err_vec .^2
    μ0λ0_vec =  λ0_err_vec .* μ0_err_vec
    μ0_sq_λ0_vec = λ0_err_vec .* μ0_sq_vec
    numer = (x_hat_err.- μ0λ0_vec) .^2 
    ssd = numer ./ denom
    half_sk_ssd =  1/2 .* (x_hat_sq_err .+ μ0_sq_λ0_vec .- ssd)
    b0_err_hat_vec = b0_err_vec .+ half_sk_ssd
    return  b0_err_hat_vec
end
function update_λ0k_signal_hat(λ0_vec,Nkj_signal) 
    K = length(Nkj_signal)
    λ0k_hat_vec = [λ0_vec .+ Nkj_signal[k] for k in 1:K]
    return λ0k_hat_vec
end
function update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
    K = length(Nkj_signal) 
    a0k_hat_vec = [ a0_vec .+ 1/2 * (Nkj_signal[k] .+ 1) for k in 1:K]
    return a0k_hat_vec
end

function update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
    K = length(Nkj_signal)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_vec .* μ0_vec
    denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    mkj_hat = [ (λ0_μ0 .+ x_hatk_signal[k]) ./denom[k] for k in 1:K]
    return mkj_hat
end


function update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
    K = length(Nkj_signal)
    denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    μ0_sq_vec = μ0_vec .^2
    μ0λ0_vec =  λ0_vec .* μ0_vec
    μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    numer = [(x_hatk_signal[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    ssd = [numer[k] ./ denom[k] for k in 1:K]
    half_sk_ssd =  1/2 .* [x_hatk_sq_signal[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    return  b0k_hat_vec
end

function update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  logpi .-e_τ_μ_tij_err[t][i]))
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
function update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tikj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tikj[t][k][j][1] for j in 1:G]
                η_false = [η_tikj[t][k][j][2] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]))
                # ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]) .+ 0.5 .* η_false .* (e_log_τj_err .-  logpi .-e_τ_μ_tij_err[t][i]))
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

function variational_inference_dynamicHDP_vs3(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            rtik = update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            η_tkj,_  = update_η_tikj3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);


            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        

        N_signal,N_error = update_N_forloops3(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        # x_hat_err = update_x_hat_error_vs_forloops(x_to_use,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops3(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x_to_use,N_error) update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops3(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end

function variational_inference_dynamicHDP_vs6(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    
    if isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_vec_init) && rand_init
        λ0_err_hat_vec_init = λ0_err_vec
    end

    if isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_vec_init) && rand_init
        m_err_hat_vec_init =μ0_err_vec
    end

    if isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = rand(Uniform(0,1),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_vec_init) && rand_init
        a0_err_hat_vec_init = a0_err_vec
    end
    if isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = rand(Uniform(0,1),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_vec_init) && rand_init
        b0_err_hat_vec_init = b0_err_vec
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err = nothing,nothing,nothing,nothing
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μ_tij_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            η_tkj,_  = update_η_tikj3(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,ηprior);
            rtik = update_rtik_vs3(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end

        
        
        

        N_signal,N_error = update_N_forloops3(rtik,η_tkj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)


        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal =  update_x_hatk_signal_vs_forloops3(x,N_signal) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error) #update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops3(x,N_signal)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end



function update_N_forloops3(rtik,η_tkj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(η_tkj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    N_error = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal =Vector{Vector{Vector{Float64}}}(undef,cells_)
        Nt_error =  Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            Nti_error = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Nti_signal[k] = rtik[t][i][k] .* [η_tkj[t][k][j][1] for j in 1:G]
                Nti_error[k] = rtik[t][i][k] .* [η_tkj[t][k][j][2] for j in 1:G]
            end
            Nt_signal[i] = Nti_signal
            Nt_error[i] = Nti_error 
        end
        N_signal[t] = Nt_signal
        N_error[t] = Nt_error
    end
    return N_signal,N_error
end
function update_signalNkj(N_signal)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_signal = reduce(vcat,N_signal)

    Nkj_signal = sum(perCell_linerize_N_signal)
    return Nkj_signal
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
function update_x_hat_error_vs_forloops2(x,N_error)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(N_error[1][1])
    G = length(x[1][1])
    x_hat_err = Vector{Float64}(undef,G)
    value = []
    for t in 1:T
        for i in 1:C_t[t]
            for k in 1:K
                push!(value, x[t][i] .* N_error[t][i][k])  
            end
        end
        # x_hat_err[j] = sum(value)
    end
    # x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    # x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return sum(value)
end
function update_x_hatk_signal_vs_forloops3(x,N_signal)
    # with η_k
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
function variational_inference_dynamicHDP_vs6logistic(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,ηprior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, η_tkj_vec_init = nothing,rtik_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


    if isnothing(mk_hat_vec_init) && rand_init
        mk_hat_vec_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_vec_init) && !rand_init
        mk_hat_vec_init = [μ0_vec for k in 1:K]
    end 
    if isnothing(λ0k_hat_vec_init) && rand_init
        λ0k_hat_vec_init = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    elseif isnothing(λ0k_hat_vec_init) && !rand_init
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]
    end
    if isnothing(a0k_hat_vec_init) && rand_init
        a0k_hat_vec_init = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    elseif isnothing(a0k_hat_vec_init) && !rand_init
        a0k_hat_vec_init = [a0_vec for k in 1:K] #
    end
    if isnothing(b0k_hat_vec_init) && rand_init
        b0k_hat_vec_init =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    elseif isnothing(b0k_hat_vec_init) && !rand_init
        b0k_hat_vec_init =  [b0_vec for k in 1:K] #
    end 
    if isnothing(rhok_hat_vec_init) || isnothing(omegak_hat_vec_init)
        if rand_init
            rhok_hat_vec_init = rand(Uniform(0,1), (K,));
            omegak_hat_vec_init = rand(Uniform(0,2), (K,));
        else
            rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K)
        end
    end
    if isnothing(a_γ_hat_init) && rand_init
        a_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(a_γ_hat_init) && !rand_init
        a_γ_hat_init = a_γ
    end
    if isnothing(b_γ_hat_init) && rand_init
        b_γ_hat_init = rand(Uniform(0,10))
    elseif isnothing(b_γ_hat_init) && !rand_init
        b_γ_hat_init = b_γ
    end


    # DYNAMIC PARAMETERS
    if isnothing(η_tkj_vec_init) && rand_init
        η_tkj_vec_init = [[[rand(Dirichlet(ones(2) ./2)) for j in 1:G] for k in 1:K] for t in 1:T]
    elseif isnothing(η_tkj_vec_init) && !rand_init
        η_tkj_vec_init =[ [[ones(2) ./2 for j in 1:G] for k in 1:K] for t in 1:T]
    end

    if isnothing(awt_hat_vec_init) && rand_init
        awt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(awt_hat_vec_init) && !rand_init
        awt_hat_vec_init = [adot_w for t in 1:T]
    end
    if isnothing(bwt_hat_vec_init) && rand_init
        bwt_hat_vec_init = [rand(Uniform(0,1)) for t in 1:T]
    elseif isnothing(bwt_hat_vec_init) && !rand_init
        bwt_hat_vec_init = [bdot_w for t in 1:T]
    end
    if isnothing(a_αt_hat_vec_init) && rand_init
        a_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(a_αt_hat_vec_init) && !rand_init
        a_αt_hat_vec_init = [a_α for t in 1:T]
    end
    if isnothing(b_αt_hat_vec_init) && rand_init
        b_αt_hat_vec_init = [rand(Uniform(0,10)) for t in 1:T]
    elseif isnothing(b_αt_hat_vec_init) && !rand_init
        b_αt_hat_vec_init = [b_α for t in 1:T]
    end
    if isnothing(c_ttprime_vec_init) && rand_init
        c_ttprime_vec_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_vec_init) && !rand_init
        c_ttprime_vec_init = [ones(T) ./T  for t in 1:T]
    end
    if isnothing(θ_hat_vec_init)
        if uniform_theta_init
            θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                θ_hat_vec_init = [rand(K+1) for t in 1:T]
            else
                θ_hat_vec_init = init_θ_hat_tk(T,rhok_hat_vec_init, omegak_hat_vec_init);
            end
        end
    end


    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:C_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T]
    end

    rtik = rtik_init

    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     
    η_tkj = η_tkj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init


    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
 
    #init debug dict initial values
 
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K


            η_tkj,_  = update_η_tkj_logistic6(Glog,rtik,e_log_τkj,e_τ_μ_tikj,ηprior);
            rtik = update_rtik_vs_logistic6(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,η_tkj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            
        end
        # println(e_τ_μ_tikj[1][1])
        # println(e_log_τkj)
        
        
        

        Nkj = update_Nkj_logistic6(rtik,η_tkj);
       

        

    


        # x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk =  update_x_hatk_logistic6(x,rtik,η_tkj) #update_x_hatk_signal_vs_forloops(x,N_signal) 


        # x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error) #update_x_hat_error_vs_forloops22(x,rtik,η_k)
        x_hatk_sq = update_x_hatk_sq_logistic6(x,rtik,η_tkj)# update_x_hatk_sq_signal_vs_forloops(x,N_signal)


        # a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        # b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_hat_logistic6(λ0_vec,Nkj)
        a0k_hat_vec = update_a0k_hat_usingXhat_logistic6(a0_vec,Nkj)
        mk_hat_vec= update_mk_hat_usingXhat_logistic6(λ0_vec,μ0_vec, Nkj,x_hatk)
        b0k_hat_vec = update_b0k_hat_usingXhat_logistic6(b0_vec,λ0_vec,μ0_vec, Nkj,x_hatk,x_hatk_sq)

        # update_λ0k_hat(λ0_vec,Nk)
        # update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        # update_a0k_hat_usingXhat(a0_vec,Nk)
        # update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

        # data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        # assgn_entropy =  calc_Hz(rtik) 
        # dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        # s_entropy = calc_Hs(c_ttprime_vec)
        # wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        # γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        # α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
      


        # iter = Int64(iter)
        # # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        # elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        # elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        # if iter > 2
        #     delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
        #     if delta_elbo <= elbo_ep || iter>=num_iter
        #         converged_bool = true
        #     end
        # end
        # iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,η_tkj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,initDict,debug_val
end
function update_λ0k_hat_logistic6(λ0_vec,Nkj_signal) 
    K = length(Nkj_signal)
    λ0k_hat_vec = [λ0_vec .+ Nkj_signal[k] for k in 1:K]
    return λ0k_hat_vec
end
function update_a0k_hat_usingXhat_logistic6(a0_vec,Nkj_signal)
    K = length(Nkj_signal) 
    a0k_hat_vec = [ a0_vec .+ 1/2 * (Nkj_signal[k] .+ 1) for k in 1:K]
    return a0k_hat_vec
end
function update_mk_hat_usingXhat_logistic6(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
    K = length(Nkj_signal)
    # Nk_xbar_k = Nk .* xbar_k
    λ0_μ0 =  λ0_vec .* μ0_vec
    denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]#update_λ0k_hat(λ0_vec,Nk)
    mkj_hat = [ (λ0_μ0 .+ x_hatk_signal[k]) ./denom[k] for k in 1:K]
    return mkj_hat
end
function update_Nkj_logistic6(rtik,η_tkj)
    T = length(rtik)
    C_t = [length(el) for el in rtik]
    K = length(rtik[1][1])
    G = length(η_tkj[1][1])
    N_rη_tkj  = [[[rtik[t][i][k] .* [η_tkj[t][k][j][1] for j in 1:G] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    N_rη_kj = [[sum(N_rη_tkj[k][t]) for t in 1:T] for k in 1:K]
    Nkj = [sum(N_rη_kj[k]) for k in 1:K]
    #sum(Ntk)[1:end-1]
    return Nkj
end
function update_b0k_hat_usingXhat_logistic6(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
    K = length(Nkj_signal)
    denom = [λ0_vec .+ Nkj_signal[k] for k in 1:K]# update_λ0k_hat(λ0_vec,Nk)
    μ0_sq_vec = μ0_vec .^2
    μ0λ0_vec =  λ0_vec .* μ0_vec
    μ0_sq_λ0_vec = λ0_vec .* μ0_sq_vec
    numer = [(x_hatk_signal[k] .- μ0λ0_vec) .^2 for k in 1:K ]
    ssd = [numer[k] ./ denom[k] for k in 1:K]
    half_sk_ssd =  1/2 .* [x_hatk_sq_signal[k] .+ μ0_sq_λ0_vec .- ssd[k] for k in 1:K] 
    # half_sk_ssd =  1/2 .* [Nk[k] .* sk[k] .+ ssd[k] for k in 1:K]
    b0k_hat_vec = [b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0_vec .+ half_sk_ssd[k] for k in 1:K]
    # b0k_hat_vec = [1.0 ./b0k_hat_vec[k] for k in 1:K]
    return  b0k_hat_vec
end
function update_rtik_vs_logistic6(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,η_tkj,c_ttprime)
    T = length(e_log_π)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    G = length(e_τ_μ_tikj[1][1][1])
    logpi = Glog/G
    for t in 1:T
        cells_ = C_t[t]
        ptik_tilde_ti = Vector{Vector{Float64}}(undef,cells_)
        adjusted_e_log_π_tk = sum([c_ttprime[t][tt] .* e_log_π[tt] for tt in 1:t])#sum(c_ttprime[t] .* e_log_π[1:t]) #sum([c_ttprime[t] .* el for el in e_log_π[1:t]])
        # println(adjusted_e_log_π_tk)
        for i in 1:cells_
            ptik_tilde_tik = Vector{Float64}(undef,K)
            for k in 1:K
                η_true = [η_tkj[t][k][j][1] for j in 1:G]
                ptik_tilde_tik[k] =  adjusted_e_log_π_tk[k] .+ sum( 0.5 .* η_true .* (e_log_τkj[k] .-  logpi .- e_τ_μ_tikj[t][i][k]))
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
function update_x_hatk_sq_logistic6(x,rtik,η_tkj)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])

    x_hat_sq_tikj = [[[ rtik[t][i][k] .* [η_tkj[t][k][j][1] for j in 1:G] .* (x[t][i]) .^2 for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_sq_tkj = [[ sum(x_hat_sq_tikj[k][t]) for t in 1:T] for k in 1:K]
    x_hat_sq_kj = [sum(x_hat_sq_tkj[k]) for k in 1:K]
    return x_hat_sq_kj
end
function update_x_hatk_logistic6(x,rtik,η_tkj)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    G = length(x[1][1])
    x_hat_tikj = [[[ rtik[t][i][k].* [η_tkj[t][k][j][1] for j in 1:G] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_tkj = [[ sum(x_hat_tikj[k][t]) for t in 1:T] for k in 1:K]
    x_hat_kj = [sum(x_hat_tkj[k]) for k in 1:K]
    return x_hat_kj
end
function update_η_tkj_logistic6(Glog,rtik,e_log_τkj,e_τ_μ_tikj,ηprior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    # η_k = Vector{Vector{Float64}}(undef,K)
    # η_k_tilde = Vector{Vector{Float64}}(undef,K)
    # for  k in 1:K
    #     η_tilde = Vector{Float64}(undef,2)
    #     η_tilde[1] = log(ηprior[k]) + sum([rtik[t][i][k] * 0.5 * (e_log_τ[k] .- Glog .- e_τ_μ[t][i][k]) for t in 1:T for i in 1:C_t[t]])
    #     η_tilde[2] = log(1 - ηprior[k])  + sum([rtik[t][i][k] * 0.5 * (sum(e_log_τj_err) .- Glog .- sum(e_τ_μ_tij_err[t][i])) for t in 1:T for i in 1:C_t[t]])
    #     η_k_tilde[k] = η_tilde
    #     η_k[k] = norm_weights(η_tilde);
    # end


    ###############################
    # [[[sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]
    # [[[sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j]) for k in 3:3] for t in 1:T] for j in 1:G]

    # [[[sum([rtik[t][i][k] for i in  1:C_t[t]]) for k in 1:K] for t in 1:T] for j in 1:G]

    η_tkj = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    η_tkj_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        η_t = Vector{Vector{Vector{Float64}}}(undef,K)
        η_t_tilde =  Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            log_η_tk = Vector{Vector{Float64}}(undef,G)
            log_η_tk_tilde = Vector{Vector{Float64}}(undef,G)
            for j in 1:G
                log_η_tkj = Vector{Float64}(undef,2)
                log_η_tkj_tilde =Vector{Float64}(undef,1)
                log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) + log(ηprior[t][k][j] /(1 - ηprior[t][k][j])) 
                # log_η_tkj_tilde[1] = sum([rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) for i in  1:C_t[t]]) - sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(ηprior[t][k][j] /(1 - ηprior[t][k][j])) 
                
                # log_η_tkj_tilde[2] = sum([ rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μ_tij_err[t][i][j]) for i in 1:C_t[t]]) + log(1 - ηprior[t][k][j])
                log_η_tkj[1] = logistic(log_η_tkj_tilde[1])
                log_η_tkj[2] = 1 - log_η_tkj[1]
                
                log_η_tk[j] = log_η_tkj
                log_η_tk_tilde[j]= log_η_tkj_tilde
            end
            η_t[k] = log_η_tk
            η_t_tilde[k] = log_η_tk_tilde
        end
        η_tkj[t] = η_t
        η_tkj_tilde[t] = η_t_tilde
    end
    return η_tkj,η_tkj_tilde
end







