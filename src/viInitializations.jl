"""
"""
function init_params_states(K)
    rho_hat_vec = 0.25 .* ones(Float64,K)
    omega_hat_vec = 2 .* ones(Float64,K)
    return rho_hat_vec, omega_hat_vec
end

######################################################

"""
"""
function init_mk_hat!(mk_hat_init,x,K,G;rand_init = false)
    μ0_vec = ones(G)
    if isnothing(mk_hat_init) && rand_init
        mk_hat_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_init) && !rand_init
        mk_hat_init = [μ0_vec for k in 1:K]
    end
    return mk_hat_init
end

"""
"""
function init_λ_sq_vec!(λ_sq_init,G;rand_init = false, lo=0,hi=1)
    λ_sq_vec = ones(G)
    if isnothing(λ_sq_init) && rand_init
        λ_sq_init = rand(Uniform(lo,hi),length(λ_sq_vec))
    elseif isnothing(λ_sq_init) && !rand_init
        λ_sq_init = λ_sq_vec
    end
    return λ_sq_init
end

"""
"""
function init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = false, lo=0,hi=1)
    σ_sq_k_vec = ones(G)
    if isnothing(σ_sq_k_init) && rand_init
        σ_sq_k_init = [rand(Uniform(lo,hi),length(σ_sq_k_vec)) for k in 1:K]
    elseif isnothing(σ_sq_k_init) && !rand_init
        σ_sq_k_init = [σ_sq_k_vec for k in 1:K] #
    end
    return σ_sq_k_init
end

"""
"""
function init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = false, lo=0,hi=1)
    v_sq_k_vec = ones(G)
    if isnothing(v_sq_k_hat_init) && rand_init
        v_sq_k_hat_init =  [rand(Uniform(lo,hi),length(v_sq_k_vec)) for k in 1:K]
    elseif isnothing(v_sq_k_hat_init) && !rand_init
        v_sq_k_hat_init =  [v_sq_k_vec for k in 1:K] #
    end 
    return v_sq_k_hat_init
end

"""
"""
function init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = false, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2)
    if isnothing(gk_hat_init) || isnothing(hk_hat_init)
        if rand_init
            gk_hat_init = rand(Uniform(g_lo,g_hi), (K,));
            hk_hat_init = rand(Uniform(h_lo,h_hi), (K,));
        else
            gk_hat_init, hk_hat_init = init_params_states(K)
        end
    end
    return gk_hat_init,hk_hat_init
end

"""
"""
function init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = false)
    if isnothing(c_ttprime_init) && rand_init
        c_ttprime_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_init) && !rand_init
        c_ttprime_init = [ones(T) ./T  for t in 1:T]
    end
    
    return c_ttprime_init
end

"""
"""
function init_d_hat_vec!(d_hat_init,K,T;rand_init = false,uniform_theta_init=false, gk_hat_init = nothing, hk_hat_init= nothing)
    if isnothing(d_hat_init)
        if uniform_theta_init
            d_hat_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                d_hat_init = [rand(K+1) for t in 1:T]
            else
                d_hat_init = init_d_hat_tk(T,gk_hat_init, hk_hat_init);
            end
        end
    end

    return d_hat_init
end

"""
"""
function init_d_hat_tk(T,g_hat_vec, h_hat_vec)
    d_hat_vec = [βk_expected_value(g_hat_vec, h_hat_vec) for t in 1:T]
    return d_hat_vec
end

"""
"""
function init_yjk_vec!(yjk_init,G,K;rand_init = false)
    if isnothing(yjk_init) && rand_init
        yjk_init = [[rand(Beta(1.,1.)) for j in 1:G] for k in 1:K] 
    elseif isnothing(yjk_init) && !rand_init
        yjk_init =[[0.5 for j in 1:G] for k in 1:K] 
    end

    return yjk_init
end

"""
"""
function init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    if isnothing(st_hat_init) && rand_init
        st_hat_init = [rand(Uniform(lo,hi)) for t in 1:T]
    elseif isnothing(st_hat_init) && !rand_init
        st_hat_init = [ϕ0 for t in 1:T]
    end
    return st_hat_init
end

"""
"""
function init_rtik_vec!(rtik_init,K,T,N_t;rand_init = false)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:N_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:N_t[t]] for t in 1:T]
    end
    

    return rtik_init
end

#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################

"""
        get_unique_time_id()
    This is an example of Docstring. This function receives two 
    numbers x and y and returns the sum of the squares.
    ```math

    ```
"""
function initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    N_t = dataparams.N_t
    K = modelparams.K
    for j in 1:G
        geneparams[j].λ_sq[1] = λ_sq_init[j]
    end
    for k in 1:K
        clusters[k].mk_hat .= copy(mk_hat_init[k])
        clusters[k].v_sq_k_hat .= copy(v_sq_k_hat_init[k])
        clusters[k].σ_sq_k_hat .= copy(σ_sq_k_init[k])
        clusters[k].var_muk .= copy(yjk_init[k] .* (mk_hat_init[k].^2 .+ v_sq_k_hat_init[k]) .- yjk_init[k] .* mk_hat_init[k].^2)
        clusters[k].κk_hat .= copy(yjk_init[k] .* mk_hat_init[k])
        clusters[k].yjk_hat .= copy(yjk_init[k])
        clusters[k].gk_hat[1] = gk_hat_init[k] 
        clusters[k].hk_hat[1] = hk_hat_init[k]
        clusters[k].ak_hat[1] = StatsFuns.logit(gk_hat_init[k])
        clusters[k].bk_hat[1] = log(hk_hat_init[k])
    end
    n = 0
    for t in 1:T
        for i in 1:N_t[t]
        n += 1
        cellpop[n].rtik .= copy(rtik_init[t][i])
        end
        conditionparams[t].c_tt_prime .= copy(c_ttprime_init[t])
        conditionparams[t].d_hat_t .= copy(d_hat_init[t])
        conditionparams[t].d_hat_t_sum[1] = sum(d_hat_init[t])
        if t ==T
            conditionparams[t].st_hat[1] = 0.0
        else
            conditionparams[t].st_hat[1] = st_hat_init[t]
        end
    end

    return cellpop,clusters,conditionparams,dataparams,modelparams,geneparams
end
