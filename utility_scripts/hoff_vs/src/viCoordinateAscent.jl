################################################################
################################################################
################# FAST FUNCTIONS (w/ Variable Selection) #######
################################################################
################################################################
################################################################
################################################################


function variational_inference_dynamicHDP_vshoff(inputs_dict;mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, ηk,α0,γ0,ϕ0, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,G;rand_init = rand_init);
    v_sq_k_hat_init = init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = rand_init, lo=0,hi=1);
    λ_sq_init = init_λ_sq_vec!(λ_sq_init,G;rand_init = rand_init, lo=0,hi=1) ;
    σ_sq_k_init = init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = rand_init, lo=0,hi=1);
    gk_hat_init,hk_hat_init = init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = rand_init, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2);
    # DYNAMIC PARAMETERS
    st_hat_init = init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    d_hat_init = init_d_hat_vec!(d_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, gk_hat_init = gk_hat_init, hk_hat_init= hk_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    yjk_init = init_yjk_vec!(yjk_init,G,K;rand_init = rand_init)

    mk_hat_vec = mk_hat_init 
    v_sq_k_hat_vec = v_sq_k_hat_init
    λ_sq_hat_vec = λ_sq_init
    σ_sq_k_vec = σ_sq_k_init
    gk_hat_vec = gk_hat_init
    hk_hat_vec = hk_hat_init
    d_hat_vec = d_hat_init
    st_hat_vec = st_hat_init
    c_ttprime_vec = c_ttprime_init
    rtik = rtik_init
    yjk=  yjk_init
    



    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i]) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,ηk,α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];
    geneparams = [GeneFeatures(j) for j in 1:G];
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init);

    Tk = nothing;
    chain_dict = nothing;

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    converged_bool = false
    is_converged = false

    update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)

    while !converged_bool #for iter in 1:num_iter
        update_rtik!(cellpop,clusters,conditionparams,dataparams,modelparams)
        # if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
        #     println("r's nan here @ iter $iter")
        # end
        update_Ntk!(cellpop,conditionparams,dataparams,modelparams)
        update_Nk!(cellpop,clusters,dataparams,modelparams)
        update_mk_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
        #     println("m's nan here @ iter $iter")
        # end
        update_c_ttprime!(conditionparams,dataparams,modelparams)
        update_yjk!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
        #     println("y's nan here @ iter $iter")
        # end
        Tk = update_Tk(clusters,conditionparams,dataparams,modelparams)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        # if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
        #     println("zero summation here @ iter $iter")
        # end
        # if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
        #     println("zero Nk here @ iter $iter")
        # end

        iter = Int64(iter)
        elbo_iter =  calculate_elbo(Tk,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
        elbo_[iter] = elbo_iter

        update_x_hat_k!(cellpop,clusters,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat's nan here @ iter $iter v1")
        # end
        update_x_hat_sq_k!(cellpop,clusters,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat_sq's nan here @ iter $iter v1")
        # end
        # if iter ==431
        #     println("*********************************")
        #     println("yjk_hat")
        #     println("*********************************")
        #     println.([clusters[k].yjk_hat for k in 1:K])
        #     println("*********************************")
        #     println("mk_hat")
        #     println("*********************************")
        #     println.([clusters[k].mk_hat for k in 1:K])
        #     println("*********************************")
        #     println("vk_sq_hat")
        #     println("*********************************")
        #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # end
        update_var_muk_hat!(clusters, dataparams,modelparams)
        # if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("var_muk's nan here @ iter $iter v1")
        # end
        update_κk_hat!(clusters, dataparams,modelparams)
        # if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("κk_hat's nan here @ iter $iter v1")
        # end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat!(clusters,conditionparams,dataparams,modelparams)
        # if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
        #     println("d's nan here @ iter $iter")
        # end
        update_d_hat_sum!(conditionparams,dataparams)

        

        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v1")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end
        update_σ_sq_k_hat!(clusters,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
        #     println("σ's nan here @ iter $iter")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end

        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v2")
        # end
        update_λ_sq_k_hat!(geneparams,clusters,dataparams,modelparams)
        # if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
        #     println("λ's nan here @ iter $iter")
        # end
        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v3")
        # end


        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams);
    Tk_ = Tk;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value;#initDict,
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value];#initDict,
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function variational_inference_dynamicHDP_vshoff_perK(inputs_dict;mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false,update_η_bool= false)
    x, K, ηk,α0,γ0,ϕ0,elbolog, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,G;rand_init = rand_init);
    v_sq_k_hat_init = init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = rand_init, lo=0,hi=1);
    λ_sq_init = init_λ_sq_vec!(λ_sq_init,G;rand_init = rand_init, lo=0,hi=1) ;
    σ_sq_k_init = init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = rand_init, lo=0,hi=1);
    gk_hat_init,hk_hat_init = init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = rand_init, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2);
    # DYNAMIC PARAMETERS
    st_hat_init = init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    d_hat_init = init_d_hat_vec!(d_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, gk_hat_init = gk_hat_init, hk_hat_init= hk_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    yjk_init = init_yjk_vec!(yjk_init,G,K;rand_init = rand_init)

    mk_hat_vec = mk_hat_init 
    v_sq_k_hat_vec = v_sq_k_hat_init
    λ_sq_hat_vec = λ_sq_init
    σ_sq_k_vec = σ_sq_k_init
    gk_hat_vec = gk_hat_init
    hk_hat_vec = hk_hat_init
    d_hat_vec = d_hat_init
    st_hat_vec = st_hat_init
    c_ttprime_vec = c_ttprime_init
    rtik = rtik_init
    yjk=  yjk_init
    



    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i]) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,ηk,α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];
    geneparams = [GeneFeatures(j) for j in 1:G];
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init);

    Tk = nothing;
    chain_dict = nothing;

    # elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    converged_bool = false
    is_converged = false

    update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
    ηk_trend_vec = []
    while !converged_bool #for iter in 1:num_iter
        update_rtik!(cellpop,clusters,conditionparams,dataparams,modelparams)
        # if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
        #     println("r's nan here @ iter $iter")
        # end
        update_Ntk!(cellpop,conditionparams,dataparams,modelparams)
        update_Nk!(cellpop,clusters,dataparams,modelparams)
        update_mk_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
        #     println("m's nan here @ iter $iter")
        # end
        update_c_ttprime!(conditionparams,dataparams,modelparams)
        update_yjk!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
        #     println("y's nan here @ iter $iter")
        # end
        Tk = update_Tk(clusters,conditionparams,dataparams,modelparams)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        # if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
        #     println("zero summation here @ iter $iter")
        # end
        # if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
        #     println("zero Nk here @ iter $iter")
        # end

        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_perK(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter

        update_x_hat_k!(cellpop,clusters,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat's nan here @ iter $iter v1")
        # end
        update_x_hat_sq_k!(cellpop,clusters,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat_sq's nan here @ iter $iter v1")
        # end
        # if iter ==431
        #     println("*********************************")
        #     println("yjk_hat")
        #     println("*********************************")
        #     println.([clusters[k].yjk_hat for k in 1:K])
        #     println("*********************************")
        #     println("mk_hat")
        #     println("*********************************")
        #     println.([clusters[k].mk_hat for k in 1:K])
        #     println("*********************************")
        #     println("vk_sq_hat")
        #     println("*********************************")
        #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # end
        update_var_muk_hat!(clusters, dataparams,modelparams)
        # if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("var_muk's nan here @ iter $iter v1")
        # end
        update_κk_hat!(clusters, dataparams,modelparams)
        # if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("κk_hat's nan here @ iter $iter v1")
        # end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat!(clusters,conditionparams,dataparams,modelparams)
        # if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
        #     println("d's nan here @ iter $iter")
        # end
        update_d_hat_sum!(conditionparams,dataparams)
        update_st_hat_fast3!(conditionparams,dataparams,modelparams) 
        

        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v1")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end
        update_σ_sq_k_hat!(clusters,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
        #     println("σ's nan here @ iter $iter")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end

        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v2")
        # end
        update_λ_sq_k_hat!(geneparams,clusters,dataparams,modelparams)
        # if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
        #     println("λ's nan here @ iter $iter")
        # end
        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v3")
        # end
        if update_η_bool
            push!(ηk_trend_vec,modelparams.ηk[1])
            update_ηk!(clusters,dataparams,modelparams)
        end
        

        if iter > 2
            delta_elbo = abs(elbolog.elbo_[iter] - elbolog.elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbolog.elbo_) .|| isnan.(elbolog.elbo_)) 
    elbo_ = elbolog.elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog);
    Tk_,ηk_trend_vec_ = Tk,ηk_trend_vec;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ;#initDict,
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ];#initDict,
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function variational_inference_dynamicHDP_vshoff_mpu(inputs_dict;mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false,update_η_bool= false,mt_mode = nothing)
    x, K, ηk,α0,γ0,ϕ0,elbolog, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,G;rand_init = rand_init);
    v_sq_k_hat_init = init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = rand_init, lo=0,hi=1);
    λ_sq_init = init_λ_sq_vec!(λ_sq_init,G;rand_init = rand_init, lo=0,hi=1) ;
    σ_sq_k_init = init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = rand_init, lo=0,hi=1);
    gk_hat_init,hk_hat_init = init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = rand_init, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2);
    # DYNAMIC PARAMETERS
    st_hat_init = init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    d_hat_init = init_d_hat_vec!(d_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, gk_hat_init = gk_hat_init, hk_hat_init= hk_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    yjk_init = init_yjk_vec!(yjk_init,G,K;rand_init = rand_init)

    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i]) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,ηk,α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];
    geneparams = [GeneFeatures(j) for j in 1:G];
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init);

    Tk = Vector{Float64}(undef,K+1);
    chain_dict = nothing;

    # elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    converged_bool = false
    is_converged = false

    update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    ηk_trend_vec = []
    while !converged_bool #for iter in 1:num_iter
        update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
        #     println("r's nan here @ iter $iter")
        # end
        update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_Nk_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        ### Cause of INF: Solution is to move it up! (BEGIN) #######
        update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat's nan here @ iter $iter v1")
        # end
        update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat_sq's nan here @ iter $iter v1")
        # end
        # if iter ==431
        #     println("*********************************")
        #     println("yjk_hat")
        #     println("*********************************")
        #     println.([clusters[k].yjk_hat for k in 1:K])
        #     println("*********************************")
        #     println("mk_hat")
        #     println("*********************************")
        #     println.([clusters[k].mk_hat for k in 1:K])
        #     println("*********************************")
        #     println("vk_sq_hat")
        #     println("*********************************")
        #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # end
        ### Cause of INF: Solution is to move it up! (END) #######


        update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
        #     println("m's nan here @ iter $iter")
        # end
        update_c_ttprime_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_yjk_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
        #     println("y's nan here @ iter $iter")
        # end
        update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        # if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
        #     println("zero summation here @ iter $iter")
        # end
        # if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
        #     println("zero Nk here @ iter $iter")
        # end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter

        ### Cause of INF (BEGIN) #######
        # update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        # #     println("x_hat's nan here @ iter $iter v1")
        # # end
        # update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        # #     println("x_hat_sq's nan here @ iter $iter v1")
        # # end
        # # if iter ==431
        # #     println("*********************************")
        # #     println("yjk_hat")
        # #     println("*********************************")
        # #     println.([clusters[k].yjk_hat for k in 1:K])
        # #     println("*********************************")
        # #     println("mk_hat")
        # #     println("*********************************")
        # #     println.([clusters[k].mk_hat for k in 1:K])
        # #     println("*********************************")
        # #     println("vk_sq_hat")
        # #     println("*********************************")
        # #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # # end

        ### Cause of INF (END) #######

        ####COMMENTED RECENTLY
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)

        # if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("var_muk's nan here @ iter $iter v1")
        # end
        
        # if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("κk_hat's nan here @ iter $iter v1")
        # end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
        #     println("d's nan here @ iter $iter")
        # end
        update_d_hat_sum_mpu!(conditionparams,dataparams;mt_mode = mt_mode)
        update_st_hat_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode) 
        

        # update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)


        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v1")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
        #     println("σ's nan here @ iter $iter")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end

        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v2")
        # end
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
        #     println("λ's nan here @ iter $iter")
        # end
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v3")
        # end
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        if update_η_bool
            push!(ηk_trend_vec,modelparams.ηk[1])
            update_ηk!(clusters,dataparams,modelparams)
        end
        

        if iter > 2
            delta_elbo = abs(elbolog.elbo_[iter] - elbolog.elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbolog.elbo_) .|| isnan.(elbolog.elbo_)) 
    elbo_ = elbolog.elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog);
    Tk_,ηk_trend_vec_ = Tk,ηk_trend_vec;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ;#initDict,
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ];#initDict,
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function convert_float_type(val)
    val_type = typeof(val)
    

end
function variational_inference_dynamicHDP_vshoff_mpu_test(inputs_dict;mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false,update_η_bool= false,mt_mode = nothing)
    x, K, ηk,α0,γ0,ϕ0,elbolog, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    # ηk = ηk_L[19]

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,G;rand_init = rand_init);
    v_sq_k_hat_init = init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = rand_init, lo=0,hi=1);
    λ_sq_init = init_λ_sq_vec!(λ_sq_init,G;rand_init = rand_init, lo=0,hi=1) ;
    σ_sq_k_init = init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = rand_init, lo=0,hi=1);
    gk_hat_init,hk_hat_init = init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = rand_init, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2);
    # DYNAMIC PARAMETERS
    st_hat_init = init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    d_hat_init = init_d_hat_vec!(d_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, gk_hat_init = gk_hat_init, hk_hat_init= hk_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    yjk_init = init_yjk_vec!(yjk_init,G,K;rand_init = rand_init)

    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[Int(t)][Int(i)]) for t in 1:T for i in 1:C_t[t]];#cellpop = [CellFeatures(t,i,K,x[t][i]) for t in 1:T for i in 1:C_t[t]]; #
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,ηk,α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];
    geneparams = [GeneFeatures(j) for j in 1:G];
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init);

    Tk = Vector{Float64}(undef,K+1);
    chain_dict = nothing;

    # elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    converged_bool = false
    is_converged = false

    update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    ηk_trend_vec = []
    any(isnan.(recursive_flatten([el.rtik for el in cellpop])))
    any(isnan.(recursive_flatten([el.yjk_hat for el in clusters])))
    any(isnan.(recursive_flatten([el.yjk_hat for el in clusters])))
    any(isnan.(recursive_flatten([el.yjk_hat for el in clusters])))
    for iter in 12:20
        recursive_flatten([el.rtik for el in cellpop])
        update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
            println("r's nan here @ iter $iter")
        end
        recursive_flatten([el.rtik for el in cellpop])

        update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        println("\n\n********************************************")
        println("Before")
        println(recursive_flatten([el.Nk for el in clusters]))
        println("********************************************")
        update_Nk_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        println("\n\n********************************************")
        println("After")
        println(recursive_flatten([el.Nk for el in clusters]))
        println("********************************************")

        update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("x_hat's nan here @ iter $iter v1")
        end
        println("\n\n********************************************")
        println("mean of x_hat")
        println(mean(recursive_flatten([el.x_hat for el in clusters])))
        println("********************************************")
        println("\n\n********************************************")
        println("max of x_hat")
        println(maximum(recursive_flatten([el.x_hat for el in clusters])))
        println("********************************************")
        println("\n\n********************************************")
        println("mean of x_hat (Per K)")
        println([mean(el.x_hat) for el in clusters])
        println("********************************************")
        println("\n\n********************************************")
        println("max of x_hat (Per K)")
        println([maximum(el.x_hat) for el in clusters])
        println("********************************************")
        # println("\n\n********************************************")
        # println("min of x_hat (Per K)")
        # println([minimum(el.x_hat) for el in clusters])
        # println("********************************************")

        update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("x_hat_sq's nan here @ iter $iter v1")
        end




        recursive_flatten([el.mk_hat for el in clusters])
        update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        
        if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
            println("m's nan here @ iter $iter")
        end
        Threads.@threads for k in 1:K
            for j in 1:G
                if typeof(clusters[k].mk_hat[j]) <: Float64 && clusters[k].mk_hat[j] > 1e154
                    # println("Here #1")
                    clusters[k].mk_hat[j] = BigFloat(clusters[k].mk_hat[j])
                end
                if typeof(clusters[k].mk_hat[j]) <: BigFloat && clusters[k].mk_hat[j] < 1e155
                    println("Here #2")
                    clusters[k].mk_hat[j] = Float64(clusters[k].mk_hat[j])
                end
            end
        end
        recursive_flatten([el.mk_hat for el in clusters])
        println("\n\n********************************************")
        println("mean of mk_hat")
        println(mean(recursive_flatten([el.mk_hat for el in clusters])))
        println("********************************************")
        println("\n\n********************************************")
        println("max of mk_hat")
        println(maximum(recursive_flatten([el.mk_hat for el in clusters])))
        println("********************************************")

        println("\n\n********************************************")
        println("mean of mk_hat (Per K)")
        println([mean(el.mk_hat) for el in clusters])
        println("********************************************")
        println("\n\n********************************************")
        println("max of mk_hat (Per K)")
        println([maximum(el.mk_hat) for el in clusters])
        println("********************************************")
        println("\n\n********************************************")
        println("min of mk_hat (Per K)")
        println([minimum(el.mk_hat) for el in clusters])
        println("********************************************")


        update_c_ttprime_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode)


        recursive_flatten([el.yjk_hat for el in clusters])
        update_yjk_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
            println("y's nan here @ iter $iter")
        end
        recursive_flatten([el.yjk_hat for el in clusters])
        println("\n\n********************************************")
        println("mean of yjk_hat")
        println(mean(recursive_flatten([el.yjk_hat for el in clusters])))
        println("********************************************")


        update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
            println("zero summation here @ iter $iter")
        end
        if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
            println("zero Nk here @ iter $iter")
        end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter

        # update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat's nan here @ iter $iter v1")
        # end
        # println("\n\n********************************************")
        # println("mean of x_hat")
        # println(mean(recursive_flatten([el.x_hat for el in clusters])))
        # println("********************************************")
        # println("\n\n********************************************")
        # println("max of x_hat")
        # println(maximum(recursive_flatten([el.x_hat for el in clusters])))
        # println("********************************************")
        # println("\n\n********************************************")
        # println("mean of x_hat (Per K)")
        # println([mean(el.x_hat) for el in clusters])
        # println("********************************************")
        # println("\n\n********************************************")
        # println("max of x_hat (Per K)")
        # println([maximum(el.x_hat) for el in clusters])
        # println("********************************************")
        # # println("\n\n********************************************")
        # # println("min of x_hat (Per K)")
        # # println([minimum(el.x_hat) for el in clusters])
        # # println("********************************************")

        # update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat_sq's nan here @ iter $iter v1")
        # end




        if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("yjk_hat's nan here @ iter $iter v1")
        end
        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v_sq_k_hat's nan here @ iter $iter v1")
        end
        if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("σ_sq_k_hat's nan here @ iter $iter v1")
        end
        if any(iszero.([clusters[k].σ_sq_k_hat for k in 1:K]))
            println("zero σ_sq_k_hat here @ iter $iter")
        end
        if any(any.([clusters[k].σ_sq_k_hat .< 0 for k in 1:K]))
            println("v_sq_k_hat < 0  here @ iter $iter")
        end
        if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("mk_hat's nan here @ iter $iter v1")
        end
        if any(iszero.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("zero v_sq_k_hat here @ iter $iter")
        end
        if any(any.([clusters[k].v_sq_k_hat .< 0 for k in 1:K]))
            println("v_sq_k_hat < 0  here @ iter $iter")
        end
        # [any(iszero.(geneparams[j].λ_sq)) for j in 1:G]
        if any([any(iszero.(geneparams[j].λ_sq)) for j in 1:G])
            println("zero λ_sq here @ iter $iter")
        end
        if any(any.([geneparams[j].λ_sq .< 0 for j in 1:G]))
            println("λ_sq < 0  here @ iter $iter")
        end
        # if iter ==10
        #     println("*********************************")
        #     println("yjk_hat")
        #     println("*********************************")
        #     println.([clusters[k].yjk_hat for k in 1:K])
        #     println("*********************************")
        #     println("mk_hat")
        #     println("*********************************")
        #     println.([clusters[k].mk_hat for k in 1:K])
        #     println("*********************************")
        #     println("vk_sq_hat")
        #     println("*********************************")
        #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # end

        ####COMMENTED RECENTLY
        recursive_flatten([el.var_muk for el in clusters])
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        recursive_flatten([el.var_muk for el in clusters])



        recursive_flatten([el.κk_hat for el in clusters])
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        recursive_flatten([el.κk_hat for el in clusters])


        if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("var_muk's nan here @ iter $iter v1")
        end
        
        if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("κk_hat's nan here @ iter $iter v1")
        end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
            println("d's nan here @ iter $iter")
        end
        update_d_hat_sum_mpu!(conditionparams,dataparams;mt_mode = mt_mode)
        update_st_hat_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode) 
        

        # update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)


        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v's nan here @ iter $iter v1")
        end

        recursive_flatten([el.σ_sq_k_hat for el in clusters])
        update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
            println("σ's nan here @ iter $iter")
        end
        recursive_flatten([el.σ_sq_k_hat for el in clusters])
        println("\n\n********************************************")
        println("mean of σ_sq_k_hat")
        println(mean(recursive_flatten([el.σ_sq_k_hat for el in clusters])))
        println("********************************************")
        println("\n\n********************************************")
        println("max of σ_sq_k_hat")
        println(maximum(recursive_flatten([el.σ_sq_k_hat for el in clusters])))
        println("********************************************")
        println("\n\n********************************************")
        println("mean of σ_sq_k_hat (Per K)")
        println([mean(el.σ_sq_k_hat) for el in clusters])
        println("********************************************")
        println("\n\n********************************************")
        println("max of σ_sq_k_hat (Per K)")
        println([maximum(el.σ_sq_k_hat) for el in clusters])
        println("********************************************")
        println("\n\n********************************************")
        println("min of σ_sq_k_hat (Per K)")
        println([minimum(el.σ_sq_k_hat) for el in clusters])
        println("********************************************")

        recursive_flatten([el.v_sq_k_hat for el in clusters])
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v's nan here @ iter $iter v2")
        end
        recursive_flatten([el.v_sq_k_hat for el in clusters])


        
        recursive_flatten([el.var_muk for el in clusters])
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        recursive_flatten([el.var_muk for el in clusters])



        recursive_flatten([el.κk_hat for el in clusters])
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        recursive_flatten([el.κk_hat for el in clusters])


        recursive_flatten([el.λ_sq for el in geneparams])
        update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
            println("λ's nan here @ iter $iter")
        end
        recursive_flatten([el.λ_sq for el in geneparams])
        println("\n\n********************************************")
        println("mean of λ_sq")
        println(mean(recursive_flatten([el.λ_sq for el in geneparams])))
        println("********************************************")
        println("\n\n********************************************")
        println("max of λ_sq")
        println(maximum(recursive_flatten([el.λ_sq for el in geneparams])))
        println("********************************************")
        println("\n\n********************************************")
        println("min of λ_sq")
        println(minimum(recursive_flatten([el.λ_sq for el in geneparams])))
        println("********************************************")

        recursive_flatten([el.v_sq_k_hat for el in clusters])
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v's nan here @ iter $iter v3")
        end
        recursive_flatten([el.v_sq_k_hat for el in clusters])

        
        recursive_flatten([el.var_muk for el in clusters])
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        recursive_flatten([el.var_muk for el in clusters])



        recursive_flatten([el.κk_hat for el in clusters])
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        recursive_flatten([el.κk_hat for el in clusters])
        
        
        if update_η_bool
            push!(ηk_trend_vec,modelparams.ηk[1])
            update_ηk!(clusters,dataparams,modelparams)
        end

        recursive_flatten([el.rtik for el in cellpop])
        recursive_flatten([el.Nk for el in clusters])
        recursive_flatten([el.yjk_hat for el in clusters])
        recursive_flatten([el.v_sq_k_hat for el in clusters])
        recursive_flatten([el.mk_hat for el in clusters])
        recursive_flatten([el.var_muk for el in clusters])
        recursive_flatten([el.κk_hat for el in clusters])
        recursive_flatten([el.σ_sq_k_hat for el in clusters])
        recursive_flatten([el.λ_sq for el in geneparams])

        println(any(isnan.(recursive_flatten([el.yjk_hat for el in clusters]))))
    end
    while !converged_bool #for iter in 1:num_iter
        update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
        #     println("r's nan here @ iter $iter")
        # end
        update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_Nk_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
        #     println("m's nan here @ iter $iter")
        # end
        update_c_ttprime_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_yjk_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
        #     println("y's nan here @ iter $iter")
        # end
        update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        # if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
        #     println("zero summation here @ iter $iter")
        # end
        # if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
        #     println("zero Nk here @ iter $iter")
        # end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter

        update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat's nan here @ iter $iter v1")
        # end
        update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat_sq's nan here @ iter $iter v1")
        # end
        # if iter ==431
        #     println("*********************************")
        #     println("yjk_hat")
        #     println("*********************************")
        #     println.([clusters[k].yjk_hat for k in 1:K])
        #     println("*********************************")
        #     println("mk_hat")
        #     println("*********************************")
        #     println.([clusters[k].mk_hat for k in 1:K])
        #     println("*********************************")
        #     println("vk_sq_hat")
        #     println("*********************************")
        #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # end

        ####COMMENTED RECENTLY
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)

        # if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("var_muk's nan here @ iter $iter v1")
        # end
        
        # if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("κk_hat's nan here @ iter $iter v1")
        # end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
        #     println("d's nan here @ iter $iter")
        # end
        update_d_hat_sum_mpu!(conditionparams,dataparams;mt_mode = mt_mode)
        update_st_hat_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode) 
        

        # update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)


        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v1")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
        #     println("σ's nan here @ iter $iter")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end

        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v2")
        # end
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
        #     println("λ's nan here @ iter $iter")
        # end
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v3")
        # end
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        if update_η_bool
            push!(ηk_trend_vec,modelparams.ηk[1])
            update_ηk!(clusters,dataparams,modelparams)
        end
        

        if iter > 2
            delta_elbo = abs(elbolog.elbo_[iter] - elbolog.elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
        #any(isnan.(recursive_flatten([el.yjk_hat for el in clusters])))
    end

    
    
    nonemptychain_indx = broadcast(!,ismissing.(elbolog.elbo_) .|| isnan.(elbolog.elbo_)) 
    elbo_ = elbolog.elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog);
    Tk_,ηk_trend_vec_ = Tk,ηk_trend_vec;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ;#initDict,
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ];#initDict,
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function variational_inference_dynamicHDP_vshoff_lowmem_mpu(inputs_dict;ep = 0.001,elbo_ep = 10^(-6),record_chain = false,update_η_bool= false,mt_mode = nothing)

    float_type, cellpop, clusters,dataparams,modelparams,conditionparams, geneparams,ηk,Tk,elbolog, num_iter, num_local_iter = (; inputs_dict...)
    chain_dict = nothing;

    # elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    converged_bool = false
    is_converged = false

    update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    ηk_trend_vec = []
    while !converged_bool #for iter in 1:num_iter
        update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
        #     println("r's nan here @ iter $iter")
        # end
        update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_Nk_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
        #     println("m's nan here @ iter $iter")
        # end
        update_c_ttprime_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_yjk_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
        #     println("y's nan here @ iter $iter")
        # end
        update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        # if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
        #     println("zero summation here @ iter $iter")
        # end
        # if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
        #     println("zero Nk here @ iter $iter")
        # end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter

        update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat's nan here @ iter $iter v1")
        # end
        update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("x_hat_sq's nan here @ iter $iter v1")
        # end
        # if iter ==431
        #     println("*********************************")
        #     println("yjk_hat")
        #     println("*********************************")
        #     println.([clusters[k].yjk_hat for k in 1:K])
        #     println("*********************************")
        #     println("mk_hat")
        #     println("*********************************")
        #     println.([clusters[k].mk_hat for k in 1:K])
        #     println("*********************************")
        #     println("vk_sq_hat")
        #     println("*********************************")
        #     println.([clusters[k].v_sq_k_hat for k in 1:K])
        # end

        ####COMMENTED RECENTLY
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)

        # if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("var_muk's nan here @ iter $iter v1")
        # end
        
        # if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("κk_hat's nan here @ iter $iter v1")
        # end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
        #     println("d's nan here @ iter $iter")
        # end
        update_d_hat_sum_mpu!(conditionparams,dataparams;mt_mode = mt_mode)
        update_st_hat_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode) 
        

        # update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)


        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v1")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end
        # update_var_muk_hat!(clusters, dataparams,modelparams)
        # update_κk_hat!(clusters, dataparams,modelparams)
        update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
        #     println("σ's nan here @ iter $iter")
        # end
        # if iter == 29
        #     println.([clusters[k].σ_sq_k_hat for k in 1:K])
        #     println("*****************************************")
        # end

        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v2")
        # end
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
        #     println("λ's nan here @ iter $iter")
        # end
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        # if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
        #     println("v's nan here @ iter $iter v3")
        # end
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        if update_η_bool
            push!(ηk_trend_vec,modelparams.ηk[1])
            update_ηk!(clusters,dataparams,modelparams)
        end
        

        if iter > 2
            delta_elbo = abs(elbolog.elbo_[iter] - elbolog.elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbolog.elbo_) .|| isnan.(elbolog.elbo_)) 
    elbo_ = elbolog.elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog);
    Tk_,ηk_trend_vec_ = Tk,ηk_trend_vec;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ;#initDict,
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ];#initDict,
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function variational_inference_dynamicHDP_vshoff2(inputs_dict;mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, ηk,α0,γ0,ϕ0, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,G;rand_init = rand_init);
    v_sq_k_hat_init = init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = rand_init, lo=0,hi=1);
    λ_sq_init = init_λ_sq_vec!(λ_sq_init,G;rand_init = rand_init, lo=0,hi=1) ;
    σ_sq_k_init = init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = rand_init, lo=0,hi=1);
    gk_hat_init,hk_hat_init = init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = rand_init, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2);
    # DYNAMIC PARAMETERS
    st_hat_init = init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    d_hat_init = init_d_hat_vec!(d_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, gk_hat_init = gk_hat_init, hk_hat_init= hk_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    yjk_init = init_yjk_vec!(yjk_init,G,K;rand_init = rand_init)

    mk_hat_vec = mk_hat_init 
    v_sq_k_hat_vec = v_sq_k_hat_init
    λ_sq_hat_vec = λ_sq_init
    σ_sq_k_vec = σ_sq_k_init
    gk_hat_vec = gk_hat_init
    hk_hat_vec = hk_hat_init
    d_hat_vec = d_hat_init
    st_hat_vec = st_hat_init
    c_ttprime_vec = c_ttprime_init
    rtik = rtik_init
    yjk=  yjk_init
    



    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i]) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,ηk,α0,γ0,ϕ0,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];
    geneparams = [GeneFeatures(j) for j in 1:G];
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init);

    Tk = nothing;
    chain_dict = nothing;

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    converged_bool = false
    is_converged = false

    update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)

    while !converged_bool #for iter in 1:num_iter
        update_rtik!(cellpop,clusters,conditionparams,dataparams,modelparams)
        if any([any(isnan.(cellpop[i].rtik)) for i in 1:sum(C_t)])
            println("r's nan here @ iter $iter")
        end
        update_Ntk!(cellpop,conditionparams,dataparams,modelparams)
        update_Nk!(cellpop,clusters,dataparams,modelparams)
        update_mk_hat!(clusters,geneparams,dataparams,modelparams)
        if any([any(isnan.(clusters[k].mk_hat)) for k in 1:K])
            println("m's nan here @ iter $iter")
        end
        update_c_ttprime!(conditionparams,dataparams,modelparams)
        update_yjk!(clusters,geneparams,dataparams,modelparams)
        if any([any(isnan.(clusters[k].yjk_hat)) for k in 1:K])
            println("y's nan here @ iter $iter")
        end
        Tk = update_Tk(clusters,conditionparams,dataparams,modelparams)
        # [clusters[k].mk_hat for k in 1:K]
        # [clusters[k].yjk_hat for k in 1:K]
        
        if any(iszero.(sum([clusters[k].yjk_hat for k in 1:K])))
            println("zero summation here @ iter $iter")
        end
        if any(iszero.([clusters[k].Nk[1] for k in 1:K]))
            println("zero Nk here @ iter $iter")
        end

        iter = Int64(iter)
        elbo_iter =  calculate_elbo(Tk,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
        elbo_[iter] = elbo_iter

        update_x_hat_k!(cellpop,clusters,dataparams,modelparams)
        if any([any(isnan.(clusters[k].x_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("x_hat's nan here @ iter $iter v1")
        end
        update_x_hat_sq_k!(cellpop,clusters,dataparams,modelparams)
        if any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("x_hat_sq's nan here @ iter $iter v1")
        end
        if iter ==431
            println("*********************************")
            println("yjk_hat")
            println("*********************************")
            println.([clusters[k].yjk_hat for k in 1:K])
            println("*********************************")
            println("mk_hat")
            println("*********************************")
            println.([clusters[k].mk_hat for k in 1:K])
            println("*********************************")
            println("vk_sq_hat")
            println("*********************************")
            println.([clusters[k].v_sq_k_hat for k in 1:K])
        end
        update_var_muk_hat!(clusters, dataparams,modelparams)
        if any([any(isnan.(clusters[k].var_muk)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("var_muk's nan here @ iter $iter v1")
        end
        update_κk_hat!(clusters, dataparams,modelparams)
        if any([any(isnan.(clusters[k].κk_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("κk_hat's nan here @ iter $iter v1")
        end
        # [clusters[k].x_hat for k in 1:K]
        # [clusters[k].x_hat_sq for k in 1:K]
        # [clusters[k].var_muk for k in 1:K]
        # [clusters[k].κk_hat for k in 1:K]

        update_gh_hat!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat!(clusters,conditionparams,dataparams,modelparams)
        # if any([any(isnan.(conditionparams[t].d_hat_t)) for t in 1:T])
        #     println("d's nan here @ iter $iter")
        # end
        update_d_hat_sum!(conditionparams,dataparams)

        

        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v's nan here @ iter $iter v1")
        end
        if iter == 29
            println.([clusters[k].σ_sq_k_hat for k in 1:K])
            println("*****************************************")
        end
        update_σ_sq_k_hat!(clusters,dataparams,modelparams)
        if any([any(isnan.(clusters[k].σ_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].σ_sq_k_hat for k in 1:K]))
            println("σ's nan here @ iter $iter")
        end
        if iter == 29
            println.([clusters[k].σ_sq_k_hat for k in 1:K])
            println("*****************************************")
        end

        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v's nan here @ iter $iter v2")
        end
        update_λ_sq_k_hat!(geneparams,clusters,dataparams,modelparams)
        if any([any(isnan.(geneparams[j].λ_sq)) for j in 1:G])#any(isnan.([geneparams[j].λ_sq for j in 1:G]))
            println("λ's nan here @ iter $iter")
        end
        update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)
        if any([any(isnan.(clusters[k].v_sq_k_hat)) for k in 1:K])#any(isnan.([clusters[k].v_sq_k_hat for k in 1:K]))
            println("v's nan here @ iter $iter v3")
        end


        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams);
    Tk_ = Tk;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value;#initDict,
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value];#initDict,
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function variational_inference_dynamicHDP_vs25_fast3_mtall(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,pip_kj_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, null_precision, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    # DYNAMIC PARAMETERS
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= omegak_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    pip_kj_init = init_pip_kj_vec!(pip_kj_init,G,K;rand_init = rand_init)

    mk_hat_vec = mk_hat_init 
    λ0k_hat_vec = λ0k_hat_init
    a0k_hat_vec = a0k_hat_init
    b0k_hat_vec = b0k_hat_init
    rhok_hat_vec = rhok_hat_init
    omegak_hat_vec = omegak_hat_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init
    awt_hat_vec = awt_hat_init 
    bwt_hat_vec = bwt_hat_init
    a_αt_hat_vec = a_αt_hat_init 
    b_αt_hat_vec = b_αt_hat_init
    θ_hat_vec = θ_hat_init
    c_ttprime_vec = c_ttprime_init
    rtik = rtik_init
    pip_kj =  pip_kj_init
    
    e_γ,Tαk = nothing,nothing


    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i],null_precision) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,a0,b0,λ0,μ0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,null_precision,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];

    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,pip_kj_init)
    if prior_cluster_membership == true
        update_Ntk_fast3!(cellpop,conditionparams,dataparams,modelparams)
        update_c_ttprime_fast3!(conditionparams,dataparams,modelparams)
        update_θ_hat_fast3!(clusters,conditionparams,dataparams,modelparams)
        update_θ_hat_sum_fast3!(conditionparams,dataparams) 
        global_cluster_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
    end

    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain,prior_cluster_membership ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain,prior_cluster_membership ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    chain_dict = nothing

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    # Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            local_updates_fast3_mtall!(cellpop,clusters,conditionparams,dataparams,modelparams)
        end
        # rpip = update_N_rpip25_fast(rtik,pip_kj)
        global_cluster_updates_fast3_mtall!(cellpop,clusters,conditionparams,dataparams,modelparams)
        global_condition_updates_fast3_mtall!(clusters,conditionparams,dataparams,modelparams)
        a_γ_hat,b_γ_hat = update_γ_fast3(clusters,dataparams,modelparams)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk_fast3(clusters,conditionparams,dataparams,modelparams)
        update_rho_omega_hat_fast3!(clusters,dataparams,modelparams,e_γ,Tαk;optim_max_iter=1000)
        iter = Int64(iter)
        elbo_iter =  calculate_elbo(a_γ,b_γ,a_γ_hat,b_γ_hat,e_γ,Tαk,cellpop,clusters,conditionparams,dataparams,modelparams)
        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,conditionparams,dataparams,modelparams);
    a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = a_γ_hat,b_γ_hat,e_γ,Tαk;
    output_str_list2 = @name a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end
function variational_inference_dynamicHDP_vs25_fast3_mtbest(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,pip_kj_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, null_precision, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    # DYNAMIC PARAMETERS
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= omegak_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    pip_kj_init = init_pip_kj_vec!(pip_kj_init,G,K;rand_init = rand_init)

    mk_hat_vec = mk_hat_init 
    λ0k_hat_vec = λ0k_hat_init
    a0k_hat_vec = a0k_hat_init
    b0k_hat_vec = b0k_hat_init
    rhok_hat_vec = rhok_hat_init
    omegak_hat_vec = omegak_hat_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init
    awt_hat_vec = awt_hat_init 
    bwt_hat_vec = bwt_hat_init
    a_αt_hat_vec = a_αt_hat_init 
    b_αt_hat_vec = b_αt_hat_init
    θ_hat_vec = θ_hat_init
    c_ttprime_vec = c_ttprime_init
    rtik = rtik_init
    pip_kj =  pip_kj_init
    
    e_γ,Tαk = nothing,nothing


    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i],null_precision) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,a0,b0,λ0,μ0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,null_precision,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];

    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,pip_kj_init)

    if prior_cluster_membership == true
        update_Ntk_fast3!(cellpop,conditionparams,dataparams,modelparams)
        update_c_ttprime_fast3!(conditionparams,dataparams,modelparams)
        update_θ_hat_fast3!(clusters,conditionparams,dataparams,modelparams)
        update_θ_hat_sum_fast3!(conditionparams,dataparams) 
        global_cluster_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
    end

    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain,prior_cluster_membership ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain,prior_cluster_membership ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    chain_dict = nothing

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    # Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            local_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
            # local_updates_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)
        end
        # rpip = update_N_rpip25_fast(rtik,pip_kj)
        global_cluster_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
        # global_cluster_updates_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)
        global_condition_updates_fast3_mtbest!(clusters,conditionparams,dataparams,modelparams)
        a_γ_hat,b_γ_hat = update_γ_fast3(clusters,dataparams,modelparams)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk_fast3(clusters,conditionparams,dataparams,modelparams)
        update_rho_omega_hat_fast3!(clusters,dataparams,modelparams,e_γ,Tαk;optim_max_iter=1000)
        iter = Int64(iter)
        elbo_iter =  calculate_elbo(a_γ,b_γ,a_γ_hat,b_γ_hat,e_γ,Tαk,cellpop,clusters,conditionparams,dataparams,modelparams)
        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,conditionparams,dataparams,modelparams);
    a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = a_γ_hat,b_γ_hat,e_γ,Tαk;
    output_str_list2 = @name a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end

function variational_inference_dynamicHDP_vs25_fast3_mtbest_sparserestart(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,pip_kj_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-2),record_chain = false,sparse_restart=false,sparse_restart_smallest_n=5,sparse_restart_itr_budget=2,local_iter_ep = 10^(-4))
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, null_precision, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)

    if typeof(K) <: AbstractFloat
        K = Int(round(K))
    end
    if !isnothing(rtik_init)
        prior_cluster_membership = true
    else
        prior_cluster_membership = false

    end
    mk_hat_init = init_mk_hat!(mk_hat_init,x,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    # DYNAMIC PARAMETERS
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= omegak_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
    pip_kj_init = init_pip_kj_vec!(pip_kj_init,G,K;rand_init = rand_init)

    mk_hat_vec = mk_hat_init 
    λ0k_hat_vec = λ0k_hat_init
    a0k_hat_vec = a0k_hat_init
    b0k_hat_vec = b0k_hat_init
    rhok_hat_vec = rhok_hat_init
    omegak_hat_vec = omegak_hat_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init
    awt_hat_vec = awt_hat_init 
    bwt_hat_vec = bwt_hat_init
    a_αt_hat_vec = a_αt_hat_init 
    b_αt_hat_vec = b_αt_hat_init
    θ_hat_vec = θ_hat_init
    c_ttprime_vec = c_ttprime_init
    rtik = rtik_init
    pip_kj =  pip_kj_init
    
    e_γ,Tαk = nothing,nothing


    float_type=eltype(x[1][1])
    cellpop = [CellFeatures(t,i,K,x[t][i],null_precision) for t in 1:T for i in 1:C_t[t]];
    clusters = [ClusterFeatures(k,G;float_type=float_type) for k in 1:K];
    dataparams = DataFeatures(x);
    modelparams = ModelParameterFeatures(x,K,a0,b0,λ0,μ0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,null_precision,num_iter,num_local_iter,uniform_theta_init,rand_init);
    conditionparams = [ConditionFeatures(t,K,T;float_type=float_type) for t in 1:T];
    Kplus = K + 1
    if sparse_restart

    end
    mod_Ntk = zeros(float_type,Kplus) 
    old_Ntk = [zeros(float_type,Kplus) for t in 1:T]
    curr_Ntk = [zeros(float_type,Kplus) for t in 1:T]
    prev_Ntk = [zeros(float_type,Kplus) for t in 1:T]
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,pip_kj_init)

    if prior_cluster_membership == true
        update_Ntk_fast3!(cellpop,conditionparams,dataparams,modelparams)
        update_c_ttprime_fast3!(conditionparams,dataparams,modelparams)
        update_θ_hat_fast3!(clusters,conditionparams,dataparams,modelparams)
        update_θ_hat_sum_fast3!(conditionparams,dataparams) 
        global_cluster_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
    end

    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain,prior_cluster_membership ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain,prior_cluster_membership ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    chain_dict = nothing

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    num_accepted = 0.0
    num_rejected = 0.0
    # Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            local_updates_fast3_mtbest_sparserestart!(cellpop,clusters,conditionparams,dataparams,modelparams,iter,loc_iter,old_Ntk,mod_Ntk,sparse_restart,sparse_restart_smallest_n,sparse_restart_itr_budget)
            # local_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
            # local_updates_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)
            curr_Ntk .= [conditionparams[t].Ntk for t in 1:T]
            if loc_iter >= 2

                if all([all(abs.(prev_ .- curr_) .<= local_iter_ep ) for (prev_,curr_) in zip(prev_Ntk,curr_Ntk)])
                    # println(loc_iter)
                    break
                    
                else
                    prev_Ntk .= copy(curr_Ntk)
                    # println(loc_iter)
                end
            else
                prev_Ntk .= copy(curr_Ntk)
            end
        end

        # rpip = update_N_rpip25_fast(rtik,pip_kj)
        global_cluster_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
        # global_cluster_updates_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)
        global_condition_updates_fast3_mtbest!(clusters,conditionparams,dataparams,modelparams)
        a_γ_hat,b_γ_hat = update_γ_fast3(clusters,dataparams,modelparams)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk_fast3(clusters,conditionparams,dataparams,modelparams)
        update_rho_omega_hat_fast3!(clusters,dataparams,modelparams,e_γ,Tαk;optim_max_iter=1000)
        iter = Int64(iter)
        elbo_iter =  calculate_elbo(a_γ,b_γ,a_γ_hat,b_γ_hat,e_γ,Tαk,cellpop,clusters,conditionparams,dataparams,modelparams)
        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if sparse_restart && isone(iter % sparse_restart_itr_budget) && !isone(iter)
            if (elbo_[iter] < elbo_[iter-sparse_restart_itr_budget]) && (elbo_[iter] < elbo_[iter-1]) 
                println("Sparse Restart REJECTED")
                num_rejected += 1
                reset_Ntk_sparse_restart!(old_Ntk,conditionparams,dataparams)

            else
                println("Sparse Restart ACCEPTED")
                num_accepted += 1
            end
        end
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1
    restarts_acceptance_rate = num_accepted / (num_accepted + num_rejected)
    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    

    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];


    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,conditionparams,dataparams,modelparams);
    a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = a_γ_hat,b_γ_hat,e_γ,Tαk;
    output_str_list2 = @name a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value,restarts_acceptance_rate;
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value,restarts_acceptance_rate];
    addToDict!(outputs_dict,output_key_list2,output_var_list2);

    return outputs_dict
end


function local_updates_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # rtik = update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec);
    update_rtik!(cellpop,clusters,conditionparams,dataparams,modelparams)

    # Ntk = update_Ntk(rtik)
    update_Ntk!(cellpop,conditionparams,dataparams,modelparams)

    # c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    update_c_ttprime!(conditionparams,dataparams,modelparams)

    # θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    update_d_hat!(clusters,conditionparams,dataparams,modelparams)
    update_d_hat_sum!(conditionparams,dataparams)
    return cellpop,clusters,conditionparams,dataparams,modelparams
 
end
function local_updates_fast3_mtall!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # rtik = update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec);
    update_rtik_vs25_fast3_mt!(cellpop,clusters,conditionparams,dataparams,modelparams)

    # Ntk = update_Ntk(rtik)
    update_Ntk_fast3_mt!(cellpop,conditionparams,dataparams,modelparams)

    # c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    update_c_ttprime_fast3_mt!(conditionparams,dataparams,modelparams)

    # θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
    update_θ_hat_fast3_mt!(clusters,conditionparams,dataparams,modelparams)
    update_θ_hat_sum_fast3_mt!(conditionparams,dataparams)
    return cellpop,clusters,conditionparams,dataparams,modelparams

end
function local_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # rtik = update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec);
    update_rtik_vs25_fast3_mt!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # update_rtik_vs25_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)

    # Ntk = update_Ntk(rtik)
    update_Ntk_fast3!(cellpop,conditionparams,dataparams,modelparams)

    # c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    update_c_ttprime_fast3!(conditionparams,dataparams,modelparams)

    # θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    update_θ_hat_fast3!(clusters,conditionparams,dataparams,modelparams)
    update_θ_hat_sum_fast3!(conditionparams,dataparams) 
    return cellpop,clusters,conditionparams,dataparams,modelparams

end
function local_updates_fast3_mtbest_sparserestart!(cellpop,clusters,conditionparams,dataparams,modelparams,iter,loc_iter,old_Ntk,mod_Ntk,sparse_restart,sparse_restart_smallest_n,sparse_restart_itr_budget)
    # rtik = update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec);
    update_rtik_vs25_fast3_mt!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # update_rtik_vs25_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)

    # Ntk = update_Ntk(rtik)
    update_Ntk_fast3!(cellpop,conditionparams,dataparams,modelparams)
    if sparse_restart && iszero(iter % sparse_restart_itr_budget) && isone(loc_iter)
        println("Doing a sparse restart")
        calculate_sparse_restart!(old_Ntk,mod_Ntk,conditionparams,dataparams,sparse_restart_smallest_n)
    end
    # c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    update_c_ttprime_fast3!(conditionparams,dataparams,modelparams)

    # θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    update_θ_hat_fast3!(clusters,conditionparams,dataparams,modelparams)
    update_θ_hat_sum_fast3!(conditionparams,dataparams) 
    return cellpop,clusters,conditionparams,dataparams,modelparams,old_Ntk,mod_Ntk

end 
function global_cluster_updates_expectation!(cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
    update_mk_hat!(clusters,geneparams,dataparams,modelparams)

    update_v_sq_k_hat!(clusters,geneparams,dataparams,modelparams)

    update_yjk!(clusters,geneparams,dataparams,modelparams)

    return cellpop,clusters,conditionparams,dataparams,geneparams,modelparams
end
function global_cluster_updates_statistics!(cellpop,clusters,conditionparams,dataparams,modelparams)
    update_Nk!(cellpop,clusters,dataparams,modelparams)

    update_x_hat_k!(cellpop,clusters,dataparams,modelparams)

    update_x_hat_sq_k!(cellpop,clusters,dataparams,modelparams)

    update_var_muk_hat!(clusters, dataparams,modelparams)

    update_κk_hat!(clusters, dataparams,modelparams)

    return cellpop,clusters,conditionparams,dataparams,modelparams
end
function global_cluster_updates_maximization!(cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
    
    update_σ_sq_k_hat!(clusters,dataparams,modelparams)
    # update_σ_sq_k_hat!(cellpop,clusters,dataparams,modelparams)
    
    update_λ_sq_k_hat!(geneparams,clusters,dataparams,modelparams)

    return cellpop,clusters,conditionparams,dataparams,geneparams,modelparams
end
function global_cluster_updates_fast3!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # Nkj = update_Nkj25_fast(rtik, pip_kj);
    update_Nkj25_fast3!(cellpop,clusters,dataparams,modelparams)

    # x_hat_k = update_x_hat_k25_fast(x,rtik,pip_kj);
    update_x_hat_k25_fast3!(cellpop,clusters,dataparams,modelparams)

    # x_hat_sq_k = update_x_hat_sq_k25_fast(x,rtik,pip_kj);
    update_x_hat_sq_k25_fast3!(cellpop,clusters,dataparams,modelparams)

    # λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
    update_λ0k_hat_fast3!(clusters,dataparams,modelparams)

    # mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    update_mk_hat_usingXhat_fast3!(clusters,dataparams,modelparams)

    # a0k_hat_vec = update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    update_a0k_hat_usingXhat25_fast3!(clusters,dataparams,modelparams)

    # b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    update_b0k_hat_usingXhat_fast3!(clusters,dataparams,modelparams)

    # pip_kj =  get_gene_PIP25_fast2(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=null_precision);
    get_gene_PIP25_fast3!(cellpop,clusters,dataparams,modelparams)
    return cellpop,clusters,conditionparams,dataparams,modelparams
end
function global_cluster_updates_fast3_mtall!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # Nkj = update_Nkj25_fast(rtik, pip_kj);
    update_Nkj25_fast3_mt!(cellpop,clusters,dataparams,modelparams)

    # x_hat_k = update_x_hat_k25_fast(x,rtik,pip_kj);
    update_x_hat_k25_fast3_mt!(cellpop,clusters,dataparams,modelparams)

    # x_hat_sq_k = update_x_hat_sq_k25_fast(x,rtik,pip_kj);
    update_x_hat_sq_k25_fast3_mt!(cellpop,clusters,dataparams,modelparams)

    # λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
    update_λ0k_hat_fast3_mt!(clusters,dataparams,modelparams)

    # mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    update_mk_hat_usingXhat_fast3_mt!(clusters,dataparams,modelparams)

    # a0k_hat_vec = update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    update_a0k_hat_usingXhat25_fast3_mt!(clusters,dataparams,modelparams)

    # b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    update_b0k_hat_usingXhat_fast3_mt!(clusters,dataparams,modelparams)

    # pip_kj =  get_gene_PIP25_fast2(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=null_precision);
    get_gene_PIP25_fast3_mt!(cellpop,clusters,dataparams,modelparams)
    return cellpop,clusters,conditionparams,dataparams,modelparams
end
function global_cluster_updates_fast3_mtbest!(cellpop,clusters,conditionparams,dataparams,modelparams)
    # Nkj = update_Nkj25_fast(rtik, pip_kj);
    update_Nkj25_fast3!(cellpop,clusters,dataparams,modelparams)

    # x_hat_k = update_x_hat_k25_fast(x,rtik,pip_kj);
    update_x_hat_k25_fast3_mt!(cellpop,clusters,dataparams,modelparams)

    # x_hat_sq_k = update_x_hat_sq_k25_fast(x,rtik,pip_kj);
    update_x_hat_sq_k25_fast3_mt!(cellpop,clusters,dataparams,modelparams)

    # λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
    update_λ0k_hat_fast3!(clusters,dataparams,modelparams)

    # mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    update_mk_hat_usingXhat_fast3!(clusters,dataparams,modelparams)

    # a0k_hat_vec = update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    update_a0k_hat_usingXhat25_fast3!(clusters,dataparams,modelparams)

    # b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    update_b0k_hat_usingXhat_fast3!(clusters,dataparams,modelparams)

    # pip_kj =  get_gene_PIP25_fast2(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=null_precision);
    get_gene_PIP25_fast3_mt!(cellpop,clusters,dataparams,modelparams)
    return cellpop,clusters,conditionparams,dataparams,modelparams
end

function global_condition_updates_fast3!(clusters,conditionparams,dataparams,modelparams)
    # a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    update_αt_fast3!(clusters,conditionparams,dataparams,modelparams)

    # awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
    # bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
    update_abwt_hat_fast3!(conditionparams,dataparams,modelparams)
    return clusters,conditionparams,dataparams,modelparams

end
function global_condition_updates_fast3_mtall!(clusters,conditionparams,dataparams,modelparams)
    # a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    update_αt_fast3_mt!(clusters,conditionparams,dataparams,modelparams)

    # awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
    # bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
    update_abwt_hat_fast3_mt!(conditionparams,dataparams,modelparams)
    return clusters,conditionparams,dataparams,modelparams

end
function global_condition_updates_fast3_mtbest!(clusters,conditionparams,dataparams,modelparams)
    # a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    update_αt_fast3!(clusters,conditionparams,dataparams,modelparams)

    # awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
    # bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
    update_abwt_hat_fast3!(conditionparams,dataparams,modelparams)
    return clusters,conditionparams,dataparams,modelparams
end

function extract_cluster_paramter(paramname,clusters,modelparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    K = modelparams.K

    param =[ isone(length(getfield(clusters[k],paramname))) ? getfield(clusters[k],paramname)[1] : getfield(clusters[k],paramname) for k in 1:K] 
    return param
end
function extract_condition_paramter(paramname,conditionparams,dataparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    T = dataparams.T

    param = [ isone(length(getfield(conditionparams[t],paramname))) ? getfield(conditionparams[t],paramname)[1] : getfield(conditionparams[t],paramname) for t in 1:T]
    return param
end

function extract_gene_paramter(paramname,geneparams,dataparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    G = dataparams.G

    param = [ isone(length(getfield(geneparams[j],paramname))) ? getfield(geneparams[j],paramname)[1] : getfield(geneparams[j],paramname) for j in 1:G]
    return param
end
function extract_elbo_vals_perK(paramname,elbolog)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    # G = dataparams.G

    vals = getfield(elbolog,paramname)
    return vals
end


function extract_rtik_paramter(cellpop,dataparams)

    paramname = :rtik
    T = dataparams.T
    N_t = dataparams.N_t
    rtik = Vector{Vector{Vector{cellpop[1].BitType}}}(undef,T)
    n = 0
    for t in 1:T
        rtik[t] = Vector{Vector{cellpop[1].BitType}}(undef,N_t[t])
        for i in 1:N_t[t]
        n += 1
        rtik[t][i] = cellpop[n].rtik  
        end
    end

    return rtik
end

function extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
    cluster_params_of_interest = [:yjk_hat, :mk_hat, :v_sq_k_hat, :σ_sq_k_hat, :var_muk,:κk_hat, :Nk, :gk_hat, :hk_hat, :ak_hat, :bk_hat, :x_hat,:x_hat_sq]
    condition_params_of_interest = [:d_hat_t, :c_tt_prime, :st_hat]
    gene_params_of_interest = [:λ_sq]

    outputs_dict[:rtik_] = extract_rtik_paramter(cellpop, dataparams)
    for fn in cluster_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_cluster_paramter(fn, clusters, modelparams)
    end
    for fn in condition_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_condition_paramter(fn, conditionparams, dataparams)
    end
    for fn in gene_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_gene_paramter(fn, geneparams, dataparams)
    end
end
function extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog)
    cluster_params_of_interest = [:yjk_hat, :mk_hat, :v_sq_k_hat, :σ_sq_k_hat, :var_muk,:κk_hat, :Nk, :gk_hat, :hk_hat, :ak_hat, :bk_hat, :x_hat,:x_hat_sq]
    condition_params_of_interest = [:d_hat_t, :c_tt_prime, :st_hat]
    gene_params_of_interest = [:λ_sq]
    perK_elbos = [:per_k_elbo]
    model_params = [:ηk]

    outputs_dict[:rtik_] = extract_rtik_paramter(cellpop, dataparams)
    for fn in cluster_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_cluster_paramter(fn, clusters, modelparams)
    end
    for fn in condition_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_condition_paramter(fn, conditionparams, dataparams)
    end
    for fn in gene_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_gene_paramter(fn, geneparams, dataparams)
    end
    for fn in perK_elbos
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_elbo_vals_perK(fn, elbolog)
    end
    for fn in model_params
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = getfield(modelparams,fn)
    end
end
#####################################################
#####################################################
################# ORIGINAL FUNCTIONS ################
#####################################################
#####################################################
#####################
#####################
function variational_inference_HDP(x, G,K,γ,α0,λ0,μ0,a0,b0,num_iter, num_local_iter;mk_hat_vec_init=nothing,λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing,b0k_hat_vec_init=nothing,θ_hat_vec_init = nothing,rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    # if rand_init
    #     λ0k_hat_vec = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]#[λ0_vec for k in 1:K] # 
    #     mk_hat_vec = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]#[μ0_vec for k in 1:K]
    #     a0k_hat_vec = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]#[a0_vec for k in 1:K] #
    #     b0k_hat_vec =  [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]#[b0_vec for k in 1:K] #
    # else
    #     λ0k_hat_vec = [λ0_vec for k in 1:K]
    #     mk_hat_vec = [μ0_vec for k in 1:K]
    #     a0k_hat_vec = [a0_vec for k in 1:K]
    #     b0k_hat_vec = [b0_vec for k in 1:K]
    # end
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

    mk_hat_vec = mk_hat_vec_init;
    λ0k_hat_vec = λ0k_hat_vec_init;
    a0k_hat_vec = a0k_hat_vec_init; 
    b0k_hat_vec = b0k_hat_vec_init;
    rhok_hat_vec = rhok_hat_vec_init;
    omegak_hat_vec = omegak_hat_vec_init;
    θ_hat_vec = θ_hat_vec_init;

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,θ_hat_vec_init;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,θ_hat_vec_init];

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
    is_converged = false
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ)
            Ntk = update_Ntk(rtik)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,α0) 
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
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        a0k_hat_vec = update_a0k_hat_usingXhat(a0_vec,Nk)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        Tk = update_Tk(θ_hat_vec)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,γ,α0,Tk;optim_max_iter=1000)
        data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,γ,α0,Tk)
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
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter =  data_elbo + assgn_entropy  + HDP_surragate_elbo
        # println(iter)
        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            # println("Here #1")
            if delta_elbo <= elbo_ep || iter>=num_iter
                # println("Here #2")
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    return elbo_, rtik,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,initDict,is_converged,debug_val
end

function variational_inference_dynamicHDP(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    # if rand_init
    #     λ0k_hat_vec = [rand(Uniform(0,1),length(λ0_vec)) for k in 1:K]
    #     mk_hat_vec = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    #     a0k_hat_vec = [rand(Uniform(0,1),length(a0_vec)) for k in 1:K]
    #     b0k_hat_vec = [rand(Uniform(0,1),length(b0_vec)) for k in 1:K]
    #     awt_hat_vec = [rand(Uniform(0,1),length(adot_w)) for t in 1:T]
    #     bwt_hat_vec = [rand(Uniform(0,1),length(bdot_w)) for t in 1:T]
    #     a_αt_hat_vec = [rand(Uniform(0,10),length(a_α)) for t in 1:T]
    #     b_αt_hat_vec = [rand(Uniform(0,10),length(b_α)) for t in 1:T]
    #     a_γ_hat = rand(Uniform(0,10))
    #     b_γ_hat = rand(Uniform(0,10))
    # else
    #     λ0k_hat_vec = [λ0_vec for k in 1:K]
    #     mk_hat_vec = [μ0_vec for k in 1:K]
    #     a0k_hat_vec = [a0_vec for k in 1:K]
    #     b0k_hat_vec = [b0_vec for k in 1:K]
    #     awt_hat_vec = [adot_w for t in 1:T]
    #     bwt_hat_vec = [bdot_w for t in 1:T]
    #     a_αt_hat_vec = [a_α for t in 1:T]
    #     b_αt_hat_vec = [b_α for t in 1:T]
    #     a_γ_hat,b_γ_hat = a_γ,b_γ
    # end

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
    
    


    mk_hat_vec = mk_hat_vec_init 
    λ0k_hat_vec = λ0k_hat_vec_init
    a0k_hat_vec = a0k_hat_vec_init
    b0k_hat_vec = b0k_hat_vec_init
    rhok_hat_vec = rhok_hat_vec_init
    omegak_hat_vec = omegak_hat_vec_init
    a_γ_hat = a_γ_hat_init 
    b_γ_hat = b_γ_hat_init

     

    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init ];

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
    is_converged = false
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
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
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)

        a0k_hat_vec = update_a0k_hat_usingXhat(a0_vec,Nk)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)

        e_γ = γ_expected_value(a_γ_hat,b_γ_hat) # γ_expected_value(a_γ,b_γ)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        
        a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        
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
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo
        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
    end
    
    return elbo_, rtik,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,initDict,is_converged,debug_val
end

#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################
function tidy_variational_inference_dynamicHDP(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing, rhok_hat_init=nothing, omegak_hat_init=nothing, rtik_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),optim_max_iter = 100)

    #ORDER HERE MATTERS######
    # K,a0,a0_err,a_α, a_γ, adot_w, b0, b0_err, b_α, b_γ, bdot_w, num_iter, num_local_iter, x_input, λ0, λ0_err, μ0, μ0_err = (; sort(inputs_dict)...)
    x_input,K,a0,b0,μ0,λ0,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,num_iter,num_local_iter = (; inputs_dict...)

    T = length(x_input)
    G = length(x_input[1][1])
    N_t = [length(el) for el in x_input]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    # λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    mk_hat_init = init_mk_hat!(mk_hat_init,x_input,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= rhok_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,N_t;rand_init = rand_init)


    # η_tkj_init = init_η_tkj_vec!(η_tkj_init,G,K,T;rand_init = rand_init)
    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    # a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    # b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)
    

    xmat = tidify_data(x_input);

    λ0μ0a0b0_mat = tidify_state_hyperprior_params(λ0_vec, μ0_vec, a0_vec, b0_vec);
    aαbαawbwaγbγ_mat = tidify_scalar_hyperprior_params(a_α, b_α,adot_w,bdot_w,a_γ,b_γ);

    
    θ_hat_mat = tidify_θ(θ_hat_init);
    λ0kmka0kb0k_hat_mat = tidify_state_params(λ0k_hat_init,mk_hat_init,a0k_hat_init, b0k_hat_init);
    c_ttprime_mat=tidify_cttprime(c_ttprime_init);
    rtik_mat = tidify_rtik(rtik_init);
    aαtbαtawtbwt_hat_mat = tidify_time_params(a_αt_hat_init, b_αt_hat_init,awt_hat_init,bwt_hat_init);
    ρkωkckdk_hat_mat = tidify_clustering_params(rhok_hat_init,omegak_hat_init,ck_hat_init,dk_hat_init );
    

    
    s2e = tidy_get_next_chain_rows_indices;


    θ_hat_chain_mat = tidy_make_chain(num_iter,θ_hat_mat)
    λ0kmka0kb0k_hat_chain_mat = tidy_make_chain(num_iter,λ0kmka0kb0k_hat_mat)
    c_ttprime_chain_mat = tidy_make_chain(num_iter,c_ttprime_mat)
    
    aαtbαtawtbwt_hat_chain_mat = tidy_make_chain(num_iter,aαtbαtawtbwt_hat_mat)
    ρkωkckdk_hat_chain_mat = tidy_make_chain(num_iter,ρkωkckdk_hat_mat)
    rtik_chain_mat = tidy_make_chain(num_iter,rtik_mat)
    θ_hat_chain_mat[s2e(1,θ_hat_mat),2:end] = θ_hat_mat
    λ0kmka0kb0k_hat_chain_mat[s2e(1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
    c_ttprime_chain_mat[s2e(1,c_ttprime_mat),2:end] = c_ttprime_mat
    rtik_chain_mat[s2e(1,rtik_mat),2:end] = rtik_mat
    aαtbαtawtbwt_hat_chain_mat[s2e(1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_chain_mat[s2e(1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat

    x_hat_k_mat = nothing
    x_hat_sq_k_mat = nothing
    aγbγ_hat_mat = nothing
    e_γ_mat = nothing 
    Tαk_mat = nothing
    x_hat_k_chain_mat = nothing
    x_hat_sq_k_chain_mat = nothing
    aγbγ_hat_chain_mat = nothing
    e_γ_chain_mat = nothing
    Tαk_chain_mat= nothing
    

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    chn_iter = 1
    converged_bool = false
    is_converged = false
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            rtik_mat, θ_hat_mat, c_ttprime_mat = tidy_local_update(xmat,θ_hat_mat,λ0kmka0kb0k_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing)
        end
        x_hat_k_mat, x_hat_sq_k_mat,λ0kmka0kb0k_hat_mat, aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat = tidy_global_update(xmat,rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = optim_max_iter)

        if chn_iter == 1
            x_hat_k_chain_mat = tidy_make_chain(num_iter,x_hat_k_mat;start = 1)
            x_hat_sq_k_chain_mat = tidy_make_chain(num_iter,x_hat_sq_k_mat;start = 1)
            aγbγ_hat_chain_mat = tidy_make_chain(num_iter,aγbγ_hat_mat;start = 1)
            e_γ_chain_mat = tidy_make_chain(num_iter,e_γ_mat;start = 1)
            Tαk_chain_mat = tidy_make_chain(num_iter,Tαk_mat;start = 1)
        end

        elbo_iter = tidy_elbo_calc(xmat,rtik_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
        iter = Int64(iter)


        ## SAVE RESULTS IN CHAIN
        θ_hat_chain_mat[s2e(chn_iter+1,θ_hat_mat),2:end] = θ_hat_mat
        λ0kmka0kb0k_hat_chain_mat[s2e(chn_iter+1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
        c_ttprime_chain_mat[s2e(chn_iter+1,c_ttprime_mat),2:end] = c_ttprime_mat
        rtik_chain_mat[s2e(chn_iter+1,rtik_mat),2:end] = rtik_mat
        aαtbαtawtbwt_hat_chain_mat[s2e(chn_iter+1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
        ρkωkckdk_hat_chain_mat[s2e(chn_iter+1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
    
        x_hat_k_chain_mat[s2e(chn_iter,x_hat_k_mat),2:end] = x_hat_k_mat
        x_hat_sq_k_chain_mat[s2e(chn_iter,x_hat_sq_k_mat),2:end]  = x_hat_sq_k_mat
        aγbγ_hat_chain_mat[s2e(chn_iter,aγbγ_hat_mat),2:end]  = aγbγ_hat_mat
        e_γ_chain_mat[s2e(chn_iter,e_γ_mat),2:end]  = e_γ_mat
        Tαk_chain_mat[s2e(chn_iter,Tαk_mat),2:end]  = Tαk_mat

        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
        chn_iter += 1
    end

    nonemptychain_indx = broadcast(!,ismissing.(elbo_))
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_)

    θ_hat_chain_mat_ = tidy_truncate_chain(θ_hat_chain_mat,θ_hat_mat,truncation_value)
    λ0kmka0kb0k_hat_chain_mat_ =  tidy_truncate_chain(λ0kmka0kb0k_hat_chain_mat,λ0kmka0kb0k_hat_mat,truncation_value)
    c_ttprime_chain_mat_=  tidy_truncate_chain(c_ttprime_chain_mat,c_ttprime_mat,truncation_value)
    rtik_chain_mat_=  tidy_truncate_chain(rtik_chain_mat,rtik_mat,truncation_value)
    aαtbαtawtbwt_hat_chain_mat_=  tidy_truncate_chain(aαtbαtawtbwt_hat_chain_mat,aαtbαtawtbwt_hat_mat,truncation_value)
    ρkωkckdk_hat_chain_mat_=  tidy_truncate_chain(ρkωkckdk_hat_chain_mat,ρkωkckdk_hat_mat,truncation_value)
    # println("Failing here")
    x_hat_k_chain_mat_=  tidy_truncate_chain(x_hat_k_chain_mat,x_hat_k_mat,truncation_value)
    # println("now here")
    x_hat_sq_k_chain_mat_=  tidy_truncate_chain(x_hat_sq_k_chain_mat,x_hat_sq_k_mat,truncation_value)
    aγbγ_hat_chain_mat_=  tidy_truncate_chain(aγbγ_hat_chain_mat,aγbγ_hat_mat,truncation_value)
    e_γ_chain_mat_=  tidy_truncate_chain(e_γ_chain_mat,e_γ_mat,truncation_value)
    Tαk_chain_mat_=  tidy_truncate_chain(Tαk_chain_mat,Tαk_mat,truncation_value)
    
    
    θ_hat_mat_ = θ_hat_mat
    λ0kmka0kb0k_hat_mat_ = λ0kmka0kb0k_hat_mat
    c_ttprime_mat_ = c_ttprime_mat
    rtik_mat_ = rtik_mat
    aαtbαtawtbwt_hat_mat_ = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_mat_ = ρkωkckdk_hat_mat
    x_hat_k_mat_ = x_hat_k_mat
    x_hat_sq_k_mat_ =x_hat_sq_k_mat
    aγbγ_hat_mat_ = aγbγ_hat_mat
    e_γ_mat_ = e_γ_mat
    Tαk_mat_ = Tαk_mat

    λ0μ0a0b0_mat_ = λ0μ0a0b0_mat
    aαbαawbwaγbγ_mat_ = aαbαawbwaγbγ_mat

    output_str_list = @name elbo_,rtik_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0kmka0kb0k_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_, is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_,rtik_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0kmka0kb0k_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_,is_converged,truncation_value];


    outputs_dict = Dict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);
    return outputs_dict
    
end
function tidy_variational_inference_dynamicHDP_SparseVS(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing, rhok_hat_init=nothing, omegak_hat_init=nothing, ηtkj_init = nothing,rtik_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),optim_max_iter = 100)
    # ,λ0_err_hat_init=nothing,m_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing
    #ORDER HERE MATTERS######
    # K,a0,a0_err,a_α, a_γ, adot_w, b0, b0_err, b_α, b_γ, bdot_w, num_iter, num_local_iter, x_input,ηtkj_prior, λ0, λ0_err, μ0, μ0_err = (; sort(inputs_dict)...)
    x_input, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, ηtkj_prior, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x_input)
    G = length(x_input[1][1])
    N_t = [length(el) for el in x_input]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    # λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    mk_hat_init = init_mk_hat!(mk_hat_init,x_input,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= rhok_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,N_t;rand_init = rand_init)



    ηtkj_init = init_η_tkj_vec!(ηtkj_init,G,K,T;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    # a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    # b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)
    

    xmat = tidify_data(x_input);

    λ0μ0a0b0_mat = tidify_state_hyperprior_params(λ0_vec, μ0_vec, a0_vec, b0_vec);
    aαbαawbwaγbγ_mat = tidify_scalar_hyperprior_params(a_α, b_α,adot_w,bdot_w,a_γ,b_γ);

    
    θ_hat_mat = tidify_θ(θ_hat_init);
    λ0kmka0kb0k_hat_mat = tidify_state_params(λ0k_hat_init,mk_hat_init,a0k_hat_init, b0k_hat_init);
    c_ttprime_mat=tidify_cttprime(c_ttprime_init);
    rtik_mat = tidify_rtik(rtik_init);
    aαtbαtawtbwt_hat_mat = tidify_time_params(a_αt_hat_init, b_αt_hat_init,awt_hat_init,bwt_hat_init);
    ρkωkckdk_hat_mat = tidify_clustering_params(rhok_hat_init,omegak_hat_init,ck_hat_init,dk_hat_init );
    ηtkj_mat = tidify_ηtkj(ηtkj_init)
    ηtkj_prior_mat = tidify_ηtkj(ηtkj_prior)

    
    s2e = tidy_get_next_chain_rows_indices;


    θ_hat_chain_mat = tidy_make_chain(num_iter,θ_hat_mat)
    λ0kmka0kb0k_hat_chain_mat = tidy_make_chain(num_iter,λ0kmka0kb0k_hat_mat)
    c_ttprime_chain_mat = tidy_make_chain(num_iter,c_ttprime_mat)
    
    aαtbαtawtbwt_hat_chain_mat = tidy_make_chain(num_iter,aαtbαtawtbwt_hat_mat)
    ρkωkckdk_hat_chain_mat = tidy_make_chain(num_iter,ρkωkckdk_hat_mat)
    rtik_chain_mat = tidy_make_chain(num_iter,rtik_mat)
    ηtkj_chain_mat = tidy_make_chain(num_iter,ηtkj_mat)

    θ_hat_chain_mat[s2e(1,θ_hat_mat),2:end] = θ_hat_mat
    λ0kmka0kb0k_hat_chain_mat[s2e(1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
    c_ttprime_chain_mat[s2e(1,c_ttprime_mat),2:end] = c_ttprime_mat
    rtik_chain_mat[s2e(1,rtik_mat),2:end] = rtik_mat
    aαtbαtawtbwt_hat_chain_mat[s2e(1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_chain_mat[s2e(1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
    ηtkj_chain_mat[s2e(1,ηtkj_mat),2:end] = ηtkj_mat

    x_hat_k_mat = nothing
    x_hat_sq_k_mat = nothing
    aγbγ_hat_mat = nothing
    e_γ_mat = nothing 
    Tαk_mat = nothing
    x_hat_k_chain_mat = nothing
    x_hat_sq_k_chain_mat = nothing
    aγbγ_hat_chain_mat = nothing
    e_γ_chain_mat = nothing
    Tαk_chain_mat= nothing
    

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    chn_iter = 1
    converged_bool = false
    is_converged = false
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            rtik_mat, ηtkj_mat, θ_hat_mat, c_ttprime_mat,log_ηtkj_tilde_mat = tidy_local_update_SparseVS(xmat,rtik_mat,ηtkj_mat,θ_hat_mat,λ0kmka0kb0k_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat,ηtkj_prior_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing)
        end
        x_hat_k_mat, x_hat_sq_k_mat,λ0kmka0kb0k_hat_mat, aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat = tidy_global_update_SparseVS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = optim_max_iter)

        if chn_iter == 1
            x_hat_k_chain_mat = tidy_make_chain(num_iter,x_hat_k_mat;start = 1)
            x_hat_sq_k_chain_mat = tidy_make_chain(num_iter,x_hat_sq_k_mat;start = 1)
            aγbγ_hat_chain_mat = tidy_make_chain(num_iter,aγbγ_hat_mat;start = 1)
            e_γ_chain_mat = tidy_make_chain(num_iter,e_γ_mat;start = 1)
            Tαk_chain_mat = tidy_make_chain(num_iter,Tαk_mat;start = 1)
        end

        elbo_iter = tidy_elbo_calc_SparseVS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
        iter = Int64(iter)


        ## SAVE RESULTS IN CHAIN
        θ_hat_chain_mat[s2e(chn_iter+1,θ_hat_mat),2:end] = θ_hat_mat
        λ0kmka0kb0k_hat_chain_mat[s2e(chn_iter+1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
        c_ttprime_chain_mat[s2e(chn_iter+1,c_ttprime_mat),2:end] = c_ttprime_mat
        rtik_chain_mat[s2e(chn_iter+1,rtik_mat),2:end] = rtik_mat
        aαtbαtawtbwt_hat_chain_mat[s2e(chn_iter+1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
        ρkωkckdk_hat_chain_mat[s2e(chn_iter+1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
        ηtkj_chain_mat[s2e(chn_iter+1,ηtkj_mat),2:end] = ηtkj_mat
    
        x_hat_k_chain_mat[s2e(chn_iter,x_hat_k_mat),2:end] = x_hat_k_mat
        x_hat_sq_k_chain_mat[s2e(chn_iter,x_hat_sq_k_mat),2:end]  = x_hat_sq_k_mat
        aγbγ_hat_chain_mat[s2e(chn_iter,aγbγ_hat_mat),2:end]  = aγbγ_hat_mat
        e_γ_chain_mat[s2e(chn_iter,e_γ_mat),2:end]  = e_γ_mat
        Tαk_chain_mat[s2e(chn_iter,Tαk_mat),2:end]  = Tαk_mat

        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
        chn_iter += 1
    end

    nonemptychain_indx = broadcast(!,ismissing.(elbo_))
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_)

    θ_hat_chain_mat_ = tidy_truncate_chain(θ_hat_chain_mat,θ_hat_mat,truncation_value)
    λ0kmka0kb0k_hat_chain_mat_ =  tidy_truncate_chain(λ0kmka0kb0k_hat_chain_mat,λ0kmka0kb0k_hat_mat,truncation_value)
    c_ttprime_chain_mat_=  tidy_truncate_chain(c_ttprime_chain_mat,c_ttprime_mat,truncation_value)
    rtik_chain_mat_=  tidy_truncate_chain(rtik_chain_mat,rtik_mat,truncation_value)
    ηtkj_chain_mat_=  tidy_truncate_chain(ηtkj_chain_mat,ηtkj_mat,truncation_value)
    aαtbαtawtbwt_hat_chain_mat_=  tidy_truncate_chain(aαtbαtawtbwt_hat_chain_mat,aαtbαtawtbwt_hat_mat,truncation_value)
    ρkωkckdk_hat_chain_mat_=  tidy_truncate_chain(ρkωkckdk_hat_chain_mat,ρkωkckdk_hat_mat,truncation_value)
    # println("Failing here")
    x_hat_k_chain_mat_=  tidy_truncate_chain(x_hat_k_chain_mat,x_hat_k_mat,truncation_value)
    # println("now here")
    x_hat_sq_k_chain_mat_=  tidy_truncate_chain(x_hat_sq_k_chain_mat,x_hat_sq_k_mat,truncation_value)
    aγbγ_hat_chain_mat_=  tidy_truncate_chain(aγbγ_hat_chain_mat,aγbγ_hat_mat,truncation_value)
    e_γ_chain_mat_=  tidy_truncate_chain(e_γ_chain_mat,e_γ_mat,truncation_value)
    Tαk_chain_mat_=  tidy_truncate_chain(Tαk_chain_mat,Tαk_mat,truncation_value)
    
    
    θ_hat_mat_ = θ_hat_mat
    λ0kmka0kb0k_hat_mat_ = λ0kmka0kb0k_hat_mat
    c_ttprime_mat_ = c_ttprime_mat
    rtik_mat_ = rtik_mat
    ηtkj_mat_ = ηtkj_mat
    aαtbαtawtbwt_hat_mat_ = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_mat_ = ρkωkckdk_hat_mat
    x_hat_k_mat_ = x_hat_k_mat
    x_hat_sq_k_mat_ =x_hat_sq_k_mat
    aγbγ_hat_mat_ = aγbγ_hat_mat
    e_γ_mat_ = e_γ_mat
    Tαk_mat_ = Tαk_mat

    λ0μ0a0b0_mat_ = λ0μ0a0b0_mat
    aαbαawbwaγbγ_mat_ = aαbαawbwaγbγ_mat

    output_str_list = @name elbo_,rtik_mat_,ηtkj_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0kmka0kb0k_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,ηtkj_chain_mat,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_, ηtkj_chain_mat_, is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_,rtik_mat_,ηtkj_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0kmka0kb0k_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,ηtkj_chain_mat,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_,ηtkj_chain_mat_,is_converged,truncation_value];


    outputs_dict = Dict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);
    return outputs_dict
    
end

function tidy_variational_inference_dynamicHDP_VS1(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,m_err_hat_init=nothing, λ0_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing, rhok_hat_init=nothing, omegak_hat_init=nothing, ηtkj_init = nothing,rtik_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),optim_max_iter = 100)
    # ,λ0_err_hat_init=nothing,m_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing
    #ORDER HERE MATTERS######
    # K,a0,a0_err,a_α, a_γ, adot_w, b0, b0_err, b_α, b_γ, bdot_w, num_iter, num_local_iter, x_input,ηtkj_prior, λ0, λ0_err, μ0, μ0_err = (; sort(inputs_dict)...)
    x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x_input)
    G = length(x_input[1][1])
    N_t = [length(el) for el in x_input]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    mk_hat_init = init_mk_hat!(mk_hat_init,x_input,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= rhok_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,N_t;rand_init = rand_init)



    ηtkj_init = init_η_tkj_vec!(ηtkj_init,G,K,T;rand_init = rand_init)


    m_err_hat_init = init_m_err_hat!(m_err_hat_init,x_input,μ0_err_vec;rand_init = rand_init)
    λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)
    

    xmat = tidify_data(x_input);

    λ0μ0a0b0_mat = tidify_state_hyperprior_params(λ0_vec, μ0_vec, a0_vec, b0_vec);
    λ0μ0a0b0_err_mat = tidify_state_hyperprior_params(λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec);
    aαbαawbwaγbγ_mat = tidify_scalar_hyperprior_params(a_α, b_α,adot_w,bdot_w,a_γ,b_γ);

    
    θ_hat_mat = tidify_θ(θ_hat_init);
    λ0kmka0kb0k_hat_mat = tidify_state_params(λ0k_hat_init,mk_hat_init,a0k_hat_init, b0k_hat_init);
    λ0ma0b0_err_hat_mat = tidify_state_params(λ0_err_hat_init,m_err_hat_init,a0_err_hat_init, b0_err_hat_init);
    c_ttprime_mat=tidify_cttprime(c_ttprime_init);
    rtik_mat = tidify_rtik(rtik_init);
    aαtbαtawtbwt_hat_mat = tidify_time_params(a_αt_hat_init, b_αt_hat_init,awt_hat_init,bwt_hat_init);
    ρkωkckdk_hat_mat = tidify_clustering_params(rhok_hat_init,omegak_hat_init,ck_hat_init,dk_hat_init );
    ηtkj_mat = tidify_ηtkj(ηtkj_init)
    ηtkj_prior_mat = tidify_ηtkj(ηtkj_prior)

    
    s2e = tidy_get_next_chain_rows_indices;


    θ_hat_chain_mat = tidy_make_chain(num_iter,θ_hat_mat)
    λ0kmka0kb0k_hat_chain_mat = tidy_make_chain(num_iter,λ0kmka0kb0k_hat_mat)
    λ0ma0b0_err_hat_chain_mat = tidy_make_chain(num_iter,λ0ma0b0_err_hat_mat)
    c_ttprime_chain_mat = tidy_make_chain(num_iter,c_ttprime_mat)
    
    aαtbαtawtbwt_hat_chain_mat = tidy_make_chain(num_iter,aαtbαtawtbwt_hat_mat)
    ρkωkckdk_hat_chain_mat = tidy_make_chain(num_iter,ρkωkckdk_hat_mat)
    rtik_chain_mat = tidy_make_chain(num_iter,rtik_mat)
    ηtkj_chain_mat = tidy_make_chain(num_iter,ηtkj_mat)

    θ_hat_chain_mat[s2e(1,θ_hat_mat),2:end] = θ_hat_mat
    λ0kmka0kb0k_hat_chain_mat[s2e(1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
    λ0ma0b0_err_hat_chain_mat[s2e(1,λ0ma0b0_err_hat_mat),2:end] = λ0ma0b0_err_hat_mat
    c_ttprime_chain_mat[s2e(1,c_ttprime_mat),2:end] = c_ttprime_mat
    rtik_chain_mat[s2e(1,rtik_mat),2:end] = rtik_mat
    aαtbαtawtbwt_hat_chain_mat[s2e(1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_chain_mat[s2e(1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
    ηtkj_chain_mat[s2e(1,ηtkj_mat),2:end] = ηtkj_mat

    x_hat_k_mat = nothing
    x_hat_sq_k_mat = nothing
    x_hat_e_mat = nothing
    x_hat_sq_e_mat = nothing
    aγbγ_hat_mat = nothing
    e_γ_mat = nothing 
    Tαk_mat = nothing
    x_hat_k_chain_mat = nothing
    x_hat_sq_k_chain_mat = nothing
    x_hat_e_chain_mat = nothing
    x_hat_sq_e_chain_mat = nothing
    aγbγ_hat_chain_mat = nothing
    e_γ_chain_mat = nothing
    Tαk_chain_mat= nothing
    

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    chn_iter = 1
    converged_bool = false
    is_converged = false
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            rtik_mat, ηtkj_mat, θ_hat_mat, c_ttprime_mat,log_ηtkj_tilde_mat = tidy_local_update_VS1(xmat,rtik_mat,ηtkj_mat,θ_hat_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat,ηtkj_prior_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing, e_log_τj_err_mat_init = nothing,e_τ_μ_tij_err_mat_init = nothing)
        end
        x_hat_k_mat, x_hat_sq_k_mat,x_hat_e_mat,x_hat_sq_e_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat = tidy_global_update_VS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = optim_max_iter)
        if chn_iter == 1
            x_hat_k_chain_mat = tidy_make_chain(num_iter,x_hat_k_mat;start = 1)
            x_hat_sq_k_chain_mat = tidy_make_chain(num_iter,x_hat_sq_k_mat;start = 1)
            x_hat_e_chain_mat = tidy_make_chain(num_iter,x_hat_e_mat;start = 1)
            x_hat_sq_e_chain_mat = tidy_make_chain(num_iter,x_hat_sq_e_mat;start = 1)
            aγbγ_hat_chain_mat = tidy_make_chain(num_iter,aγbγ_hat_mat;start = 1)
            e_γ_chain_mat = tidy_make_chain(num_iter,e_γ_mat;start = 1)
            Tαk_chain_mat = tidy_make_chain(num_iter,Tαk_mat;start = 1)
        end

        elbo_iter = tidy_elbo_calc_VS1(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
        iter = Int64(iter)


        ## SAVE RESULTS IN CHAIN
        θ_hat_chain_mat[s2e(chn_iter+1,θ_hat_mat),2:end] = θ_hat_mat
        λ0kmka0kb0k_hat_chain_mat[s2e(chn_iter+1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
        λ0ma0b0_err_hat_chain_mat[s2e(chn_iter+1,λ0ma0b0_err_hat_mat),2:end] = λ0ma0b0_err_hat_mat
        c_ttprime_chain_mat[s2e(chn_iter+1,c_ttprime_mat),2:end] = c_ttprime_mat
        rtik_chain_mat[s2e(chn_iter+1,rtik_mat),2:end] = rtik_mat
        aαtbαtawtbwt_hat_chain_mat[s2e(chn_iter+1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
        ρkωkckdk_hat_chain_mat[s2e(chn_iter+1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
        ηtkj_chain_mat[s2e(chn_iter+1,ηtkj_mat),2:end] = ηtkj_mat
    
        x_hat_k_chain_mat[s2e(chn_iter,x_hat_k_mat),2:end] = x_hat_k_mat
        x_hat_sq_k_chain_mat[s2e(chn_iter,x_hat_sq_k_mat),2:end]  = x_hat_sq_k_mat
        x_hat_e_chain_mat[s2e(chn_iter,x_hat_e_mat),2:end] = x_hat_e_mat
        x_hat_sq_e_chain_mat[s2e(chn_iter,x_hat_sq_e_mat),2:end]  = x_hat_sq_e_mat
        aγbγ_hat_chain_mat[s2e(chn_iter,aγbγ_hat_mat),2:end]  = aγbγ_hat_mat
        e_γ_chain_mat[s2e(chn_iter,e_γ_mat),2:end]  = e_γ_mat
        Tαk_chain_mat[s2e(chn_iter,Tαk_mat),2:end]  = Tαk_mat

        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
        chn_iter += 1
    end

    nonemptychain_indx = broadcast(!,ismissing.(elbo_))
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_)

    θ_hat_chain_mat_ = tidy_truncate_chain(θ_hat_chain_mat,θ_hat_mat,truncation_value)
    λ0kmka0kb0k_hat_chain_mat_ =  tidy_truncate_chain(λ0kmka0kb0k_hat_chain_mat,λ0kmka0kb0k_hat_mat,truncation_value)
    λ0ma0b0_err_hat_chain_mat_ =  tidy_truncate_chain(λ0ma0b0_err_hat_chain_mat,λ0ma0b0_err_hat_mat,truncation_value)
    c_ttprime_chain_mat_=  tidy_truncate_chain(c_ttprime_chain_mat,c_ttprime_mat,truncation_value)
    rtik_chain_mat_=  tidy_truncate_chain(rtik_chain_mat,rtik_mat,truncation_value)
    ηtkj_chain_mat_=  tidy_truncate_chain(ηtkj_chain_mat,ηtkj_mat,truncation_value)
    aαtbαtawtbwt_hat_chain_mat_=  tidy_truncate_chain(aαtbαtawtbwt_hat_chain_mat,aαtbαtawtbwt_hat_mat,truncation_value)
    ρkωkckdk_hat_chain_mat_=  tidy_truncate_chain(ρkωkckdk_hat_chain_mat,ρkωkckdk_hat_mat,truncation_value)
    # println("Failing here")
    x_hat_k_chain_mat_=  tidy_truncate_chain(x_hat_k_chain_mat,x_hat_k_mat,truncation_value)
    # println("now here")
    x_hat_sq_k_chain_mat_=  tidy_truncate_chain(x_hat_sq_k_chain_mat,x_hat_sq_k_mat,truncation_value)
    x_hat_sq_e_chain_mat_=  tidy_truncate_chain(x_hat_sq_e_chain_mat,x_hat_sq_e_mat,truncation_value)
    x_hat_e_chain_mat_=  tidy_truncate_chain(x_hat_e_chain_mat,x_hat_e_mat,truncation_value)

    
    aγbγ_hat_chain_mat_=  tidy_truncate_chain(aγbγ_hat_chain_mat,aγbγ_hat_mat,truncation_value)
    e_γ_chain_mat_=  tidy_truncate_chain(e_γ_chain_mat,e_γ_mat,truncation_value)
    Tαk_chain_mat_=  tidy_truncate_chain(Tαk_chain_mat,Tαk_mat,truncation_value)
    
    
    θ_hat_mat_ = θ_hat_mat
    λ0kmka0kb0k_hat_mat_ = λ0kmka0kb0k_hat_mat
    λ0ma0b0_err_hat_mat_ = λ0ma0b0_err_hat_mat
    c_ttprime_mat_ = c_ttprime_mat
    rtik_mat_ = rtik_mat
    ηtkj_mat_ = ηtkj_mat
    aαtbαtawtbwt_hat_mat_ = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_mat_ = ρkωkckdk_hat_mat
    x_hat_k_mat_ = x_hat_k_mat
    x_hat_sq_k_mat_ =x_hat_sq_k_mat
    x_hat_e_mat_ = x_hat_e_mat
    x_hat_sq_e_mat_ =x_hat_sq_e_mat
    aγbγ_hat_mat_ = aγbγ_hat_mat
    e_γ_mat_ = e_γ_mat
    Tαk_mat_ = Tαk_mat

    λ0μ0a0b0_mat_ = λ0μ0a0b0_mat
    λ0μ0a0b0_err_mat_=λ0μ0a0b0_err_mat
    aαbαawbwaγbγ_mat_ = aαbαawbwaγbγ_mat

    output_str_list = @name elbo_,rtik_mat_,ηtkj_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0μ0a0b0_err_mat_,λ0kmka0kb0k_hat_mat_,λ0ma0b0_err_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,x_hat_e_mat_,x_hat_sq_e_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,λ0ma0b0_err_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,ηtkj_chain_mat,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,x_hat_e_chain_mat_,x_hat_sq_e_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_, ηtkj_chain_mat_, is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_,rtik_mat_,ηtkj_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0μ0a0b0_err_mat_,λ0kmka0kb0k_hat_mat_,λ0ma0b0_err_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,x_hat_e_mat_,x_hat_sq_e_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,λ0ma0b0_err_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,ηtkj_chain_mat,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,x_hat_e_chain_mat_,x_hat_sq_e_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_, ηtkj_chain_mat_, is_converged,truncation_value];


    outputs_dict = Dict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);
    return outputs_dict
    
end


function tidy_variational_inference_dynamicHDP_VS2(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,m_err_hat_init=nothing, λ0_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing, rhok_hat_init=nothing, omegak_hat_init=nothing, ηtkj_init = nothing,rtik_init = nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),optim_max_iter = 100)
    # ,λ0_err_hat_init=nothing,m_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing
    #ORDER HERE MATTERS######
    # K,a0,a0_err,a_α, a_γ, adot_w, b0, b0_err, b_α, b_γ, bdot_w, num_iter, num_local_iter, x_input,ηtkj_prior, λ0, λ0_err, μ0, μ0_err = (; sort(inputs_dict)...)
    x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x_input)
    G = length(x_input[1][1])
    N_t = [length(el) for el in x_input]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    mk_hat_init = init_mk_hat!(mk_hat_init,x_input,K,μ0_vec;rand_init = rand_init);
    λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1) ;
    a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
    b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
    rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    ck_hat_init,dk_hat_init = rhok_hat_init,omegak_hat_init;
    a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
    b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
    awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
    bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
    a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
    b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
    c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
    θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= rhok_hat_init)
    rtik_init = init_rtik_vec!(rtik_init,K,T,N_t;rand_init = rand_init)



    ηtkj_init = init_η_tkj_vec!(ηtkj_init,G,K,T;rand_init = rand_init)


    m_err_hat_init = init_m_err_hat!(m_err_hat_init,x_input,μ0_err_vec;rand_init = rand_init)
    λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)
    

    xmat = tidify_data(x_input);

    λ0μ0a0b0_mat = tidify_state_hyperprior_params(λ0_vec, μ0_vec, a0_vec, b0_vec);
    λ0μ0a0b0_err_mat = tidify_state_hyperprior_params(λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec);
    aαbαawbwaγbγ_mat = tidify_scalar_hyperprior_params(a_α, b_α,adot_w,bdot_w,a_γ,b_γ);

    
    θ_hat_mat = tidify_θ(θ_hat_init);
    λ0kmka0kb0k_hat_mat = tidify_state_params(λ0k_hat_init,mk_hat_init,a0k_hat_init, b0k_hat_init);
    λ0ma0b0_err_hat_mat = tidify_state_params(λ0_err_hat_init,m_err_hat_init,a0_err_hat_init, b0_err_hat_init);
    c_ttprime_mat=tidify_cttprime(c_ttprime_init);
    rtik_mat = tidify_rtik(rtik_init);
    aαtbαtawtbwt_hat_mat = tidify_time_params(a_αt_hat_init, b_αt_hat_init,awt_hat_init,bwt_hat_init);
    ρkωkckdk_hat_mat = tidify_clustering_params(rhok_hat_init,omegak_hat_init,ck_hat_init,dk_hat_init );
    ηtkj_mat = tidify_ηtkj(ηtkj_init)
    ηtkj_prior_mat = tidify_ηtkj(ηtkj_prior)

    
    s2e = tidy_get_next_chain_rows_indices;


    θ_hat_chain_mat = tidy_make_chain(num_iter,θ_hat_mat)
    λ0kmka0kb0k_hat_chain_mat = tidy_make_chain(num_iter,λ0kmka0kb0k_hat_mat)
    λ0ma0b0_err_hat_chain_mat = tidy_make_chain(num_iter,λ0ma0b0_err_hat_mat)
    c_ttprime_chain_mat = tidy_make_chain(num_iter,c_ttprime_mat)
    
    aαtbαtawtbwt_hat_chain_mat = tidy_make_chain(num_iter,aαtbαtawtbwt_hat_mat)
    ρkωkckdk_hat_chain_mat = tidy_make_chain(num_iter,ρkωkckdk_hat_mat)
    rtik_chain_mat = tidy_make_chain(num_iter,rtik_mat)
    ηtkj_chain_mat = tidy_make_chain(num_iter,ηtkj_mat)

    θ_hat_chain_mat[s2e(1,θ_hat_mat),2:end] = θ_hat_mat
    λ0kmka0kb0k_hat_chain_mat[s2e(1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
    λ0ma0b0_err_hat_chain_mat[s2e(1,λ0ma0b0_err_hat_mat),2:end] = λ0ma0b0_err_hat_mat
    c_ttprime_chain_mat[s2e(1,c_ttprime_mat),2:end] = c_ttprime_mat
    rtik_chain_mat[s2e(1,rtik_mat),2:end] = rtik_mat
    aαtbαtawtbwt_hat_chain_mat[s2e(1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_chain_mat[s2e(1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
    ηtkj_chain_mat[s2e(1,ηtkj_mat),2:end] = ηtkj_mat

    x_hat_k_mat = nothing
    x_hat_sq_k_mat = nothing
    x_hat_e_mat = nothing
    x_hat_sq_e_mat = nothing
    aγbγ_hat_mat = nothing
    e_γ_mat = nothing 
    Tαk_mat = nothing
    x_hat_k_chain_mat = nothing
    x_hat_sq_k_chain_mat = nothing
    x_hat_e_chain_mat = nothing
    x_hat_sq_e_chain_mat = nothing
    aγbγ_hat_chain_mat = nothing
    e_γ_chain_mat = nothing
    Tαk_chain_mat= nothing
    

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    iter = 1
    chn_iter = 1
    converged_bool = false
    is_converged = false
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            rtik_mat, ηtkj_mat, θ_hat_mat, c_ttprime_mat,log_ηtkj_tilde_mat = tidy_local_update_VS2(xmat,rtik_mat,ηtkj_mat,θ_hat_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,c_ttprime_mat,aαtbαtawtbwt_hat_mat,ρkωkckdk_hat_mat,ηtkj_prior_mat; e_log_π_mat_init= nothing,e_log_τ_kj_mat_init = nothing,rtik_mat_init = nothing, e_log_τj_err_mat_init = nothing,e_τ_μ_tij_err_mat_init = nothing)
        end
        x_hat_k_mat, x_hat_sq_k_mat,x_hat_e_mat,x_hat_sq_e_mat,λ0kmka0kb0k_hat_mat,λ0ma0b0_err_hat_mat,aγbγ_hat_mat,e_γ_mat,Tαk_mat,ρkωkckdk_hat_mat = tidy_global_update_VS(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,λ0μ0a0b0_err_mat,λ0ma0b0_err_hat_mat,ρkωkckdk_hat_mat,θ_hat_mat,c_ttprime_mat,aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat;optim_max_iter = optim_max_iter)
        if chn_iter == 1
            x_hat_k_chain_mat = tidy_make_chain(num_iter,x_hat_k_mat;start = 1)
            x_hat_sq_k_chain_mat = tidy_make_chain(num_iter,x_hat_sq_k_mat;start = 1)
            x_hat_e_chain_mat = tidy_make_chain(num_iter,x_hat_e_mat;start = 1)
            x_hat_sq_e_chain_mat = tidy_make_chain(num_iter,x_hat_sq_e_mat;start = 1)
            aγbγ_hat_chain_mat = tidy_make_chain(num_iter,aγbγ_hat_mat;start = 1)
            e_γ_chain_mat = tidy_make_chain(num_iter,e_γ_mat;start = 1)
            Tαk_chain_mat = tidy_make_chain(num_iter,Tαk_mat;start = 1)
        end

        elbo_iter = tidy_elbo_calc_VS2(xmat,rtik_mat,ηtkj_mat,λ0μ0a0b0_mat,λ0kmka0kb0k_hat_mat,aγbγ_hat_mat,ρkωkckdk_hat_mat,Tαk_mat,c_ttprime_mat, aαbαawbwaγbγ_mat,aαtbαtawtbwt_hat_mat)
        iter = Int64(iter)


        ## SAVE RESULTS IN CHAIN
        θ_hat_chain_mat[s2e(chn_iter+1,θ_hat_mat),2:end] = θ_hat_mat
        λ0kmka0kb0k_hat_chain_mat[s2e(chn_iter+1,λ0kmka0kb0k_hat_mat),2:end] = λ0kmka0kb0k_hat_mat
        λ0ma0b0_err_hat_chain_mat[s2e(chn_iter+1,λ0ma0b0_err_hat_mat),2:end] = λ0ma0b0_err_hat_mat
        c_ttprime_chain_mat[s2e(chn_iter+1,c_ttprime_mat),2:end] = c_ttprime_mat
        rtik_chain_mat[s2e(chn_iter+1,rtik_mat),2:end] = rtik_mat
        aαtbαtawtbwt_hat_chain_mat[s2e(chn_iter+1,aαtbαtawtbwt_hat_mat),2:end] = aαtbαtawtbwt_hat_mat
        ρkωkckdk_hat_chain_mat[s2e(chn_iter+1,ρkωkckdk_hat_mat),2:end] = ρkωkckdk_hat_mat
        ηtkj_chain_mat[s2e(chn_iter+1,ηtkj_mat),2:end] = ηtkj_mat
    
        x_hat_k_chain_mat[s2e(chn_iter,x_hat_k_mat),2:end] = x_hat_k_mat
        x_hat_sq_k_chain_mat[s2e(chn_iter,x_hat_sq_k_mat),2:end]  = x_hat_sq_k_mat
        x_hat_e_chain_mat[s2e(chn_iter,x_hat_e_mat),2:end] = x_hat_e_mat
        x_hat_sq_e_chain_mat[s2e(chn_iter,x_hat_sq_e_mat),2:end]  = x_hat_sq_e_mat
        aγbγ_hat_chain_mat[s2e(chn_iter,aγbγ_hat_mat),2:end]  = aγbγ_hat_mat
        e_γ_chain_mat[s2e(chn_iter,e_γ_mat),2:end]  = e_γ_mat
        Tαk_chain_mat[s2e(chn_iter,Tαk_mat),2:end]  = Tαk_mat

        elbo_[iter] = elbo_iter# push!(elbo_, elbo_iter)
        if iter > 2
            delta_elbo = abs(elbo_[iter] - elbo_[iter-1])
            if delta_elbo <= elbo_ep || iter>=num_iter
                converged_bool = true
                if iter>=num_iter && delta_elbo > elbo_ep
                    is_converged = false

                else
                    is_converged = true
                end
            end
        end
        iter += 1
        chn_iter += 1
    end

    nonemptychain_indx = broadcast(!,ismissing.(elbo_))
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_)

    θ_hat_chain_mat_ = tidy_truncate_chain(θ_hat_chain_mat,θ_hat_mat,truncation_value)
    λ0kmka0kb0k_hat_chain_mat_ =  tidy_truncate_chain(λ0kmka0kb0k_hat_chain_mat,λ0kmka0kb0k_hat_mat,truncation_value)
    λ0ma0b0_err_hat_chain_mat_ =  tidy_truncate_chain(λ0ma0b0_err_hat_chain_mat,λ0ma0b0_err_hat_mat,truncation_value)
    c_ttprime_chain_mat_=  tidy_truncate_chain(c_ttprime_chain_mat,c_ttprime_mat,truncation_value)
    rtik_chain_mat_=  tidy_truncate_chain(rtik_chain_mat,rtik_mat,truncation_value)
    ηtkj_chain_mat_=  tidy_truncate_chain(ηtkj_chain_mat,ηtkj_mat,truncation_value)
    aαtbαtawtbwt_hat_chain_mat_=  tidy_truncate_chain(aαtbαtawtbwt_hat_chain_mat,aαtbαtawtbwt_hat_mat,truncation_value)
    ρkωkckdk_hat_chain_mat_=  tidy_truncate_chain(ρkωkckdk_hat_chain_mat,ρkωkckdk_hat_mat,truncation_value)
    # println("Failing here")
    x_hat_k_chain_mat_=  tidy_truncate_chain(x_hat_k_chain_mat,x_hat_k_mat,truncation_value)
    # println("now here")
    x_hat_sq_k_chain_mat_=  tidy_truncate_chain(x_hat_sq_k_chain_mat,x_hat_sq_k_mat,truncation_value)
    x_hat_sq_e_chain_mat_=  tidy_truncate_chain(x_hat_sq_e_chain_mat,x_hat_sq_e_mat,truncation_value)
    x_hat_e_chain_mat_=  tidy_truncate_chain(x_hat_e_chain_mat,x_hat_e_mat,truncation_value)

    
    aγbγ_hat_chain_mat_=  tidy_truncate_chain(aγbγ_hat_chain_mat,aγbγ_hat_mat,truncation_value)
    e_γ_chain_mat_=  tidy_truncate_chain(e_γ_chain_mat,e_γ_mat,truncation_value)
    Tαk_chain_mat_=  tidy_truncate_chain(Tαk_chain_mat,Tαk_mat,truncation_value)
    
    
    θ_hat_mat_ = θ_hat_mat
    λ0kmka0kb0k_hat_mat_ = λ0kmka0kb0k_hat_mat
    λ0ma0b0_err_hat_mat_ = λ0ma0b0_err_hat_mat
    c_ttprime_mat_ = c_ttprime_mat
    rtik_mat_ = rtik_mat
    ηtkj_mat_ = ηtkj_mat
    aαtbαtawtbwt_hat_mat_ = aαtbαtawtbwt_hat_mat
    ρkωkckdk_hat_mat_ = ρkωkckdk_hat_mat
    x_hat_k_mat_ = x_hat_k_mat
    x_hat_sq_k_mat_ =x_hat_sq_k_mat
    x_hat_e_mat_ = x_hat_e_mat
    x_hat_sq_e_mat_ =x_hat_sq_e_mat
    aγbγ_hat_mat_ = aγbγ_hat_mat
    e_γ_mat_ = e_γ_mat
    Tαk_mat_ = Tαk_mat

    λ0μ0a0b0_mat_ = λ0μ0a0b0_mat
    λ0μ0a0b0_err_mat_=λ0μ0a0b0_err_mat
    aαbαawbwaγbγ_mat_ = aαbαawbwaγbγ_mat

    output_str_list = @name elbo_,rtik_mat_,ηtkj_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0μ0a0b0_err_mat_,λ0kmka0kb0k_hat_mat_,λ0ma0b0_err_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,x_hat_e_mat_,x_hat_sq_e_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,λ0ma0b0_err_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,ηtkj_chain_mat,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,x_hat_e_chain_mat_,x_hat_sq_e_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_, ηtkj_chain_mat_, is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_,rtik_mat_,ηtkj_mat_,θ_hat_mat_,λ0μ0a0b0_mat_,λ0μ0a0b0_err_mat_,λ0kmka0kb0k_hat_mat_,λ0ma0b0_err_hat_mat_,c_ttprime_mat_,aαbαawbwaγbγ_mat_,aαtbαtawtbwt_hat_mat_,ρkωkckdk_hat_mat_,x_hat_k_mat_,x_hat_sq_k_mat_,x_hat_e_mat_,x_hat_sq_e_mat_,aγbγ_hat_mat_,e_γ_mat_,Tαk_mat_, θ_hat_chain_mat_, λ0kmka0kb0k_hat_chain_mat_,λ0ma0b0_err_hat_chain_mat_,c_ttprime_chain_mat_,rtik_chain_mat_,aαtbαtawtbwt_hat_chain_mat_,ρkωkckdk_hat_chain_mat_,ηtkj_chain_mat,x_hat_k_chain_mat_,x_hat_sq_k_chain_mat_,x_hat_e_chain_mat_,x_hat_sq_e_chain_mat_,aγbγ_hat_chain_mat_,e_γ_chain_mat_,Tαk_chain_mat_, ηtkj_chain_mat_, is_converged,truncation_value];


    outputs_dict = Dict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);
    return outputs_dict
    
end


function depracated_variational_inference(x, G,K,γ,α0,λ0,μ0,a0,b0,num_iter, num_local_iter; rand_init = false,ep = 0.001)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0k_hat_vec = [λ0_vec for k in 1:K]
    mk_hat_vec = [μ0_vec for k in 1:K]
    a0k_hat_vec = [a0_vec for k in 1:K]
    b0k_hat_vec = [b0_vec for k in 1:K]
    rhok_hat_vec, omegak_hat_vec = init_params_states(K)
    θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)#[ones(K+1) ./K  for t in 1:T]#
    
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Float64}(undef,num_iter)
    for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat) # T by K
            e_log_τ = sum.(depracated_log_τ_expected_value.(a0k_hat_vec, b0k_hat_vec)) # K by 1
            e_τ_μ_tikj = depracated_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G
            e_τ_μ = [[[ sum(e_τ_μ_tikj[t][i][k]) for k in 1:K] for i in 1:C_t[t]] for t in 1:T] # T by C_t by K
            rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ)
            Ntk = update_Ntk(rtik)
            θ_hat = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,α0) 
        end

        Nk = update_Nk(rtik)
        xbar_k = update_xbar_k(x,rtik)
        # xbar_k = 1 ./ Nk .* xbar_k
        sk = update_sk_GIndepGaussian(x,xbar_k,rtik)
        # sk = 1 ./ Nk .* sk
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        mk_hat_vec= update_mk_hat(λ0_vec,μ0_vec, Nk,xbar_k)

        a0k_hat_vec =  update_a0k_hat(a0_vec,Nk)
        b0k_hat_vec = update_b0k_hat(b0_vec,λ0_vec,μ0_vec, Nk,xbar_k,sk)

        Tk = update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk;optim_max_iter=1000)
        data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,G,γ,α0,Tk)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_[iter] = data_elbo + assgn_entropy  + HDP_surragate_elbo
        
    end
    return elbo_, rtik,θ_hat, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec
end

function depracated_variational_inference_notusingXhat(x, G,K,γ,α0,λ0,μ0,a0,b0,num_iter, num_local_iter; rand_init = false,ep = 0.001,debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0k_hat_vec = [λ0_vec for k in 1:K]
    mk_hat_vec = [μ0_vec for k in 1:K]
    a0k_hat_vec = [a0_vec for k in 1:K]
    b0k_hat_vec = [b0_vec for k in 1:K]
    rhok_hat_vec, omegak_hat_vec = init_params_states(K)
    θ_hat = [ones(K+1) ./(K+1)  for t in 1:T]#init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
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
    elbo_ = Vector{Float64}(undef,num_iter)
    for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ)
            Ntk = update_Ntk(rtik)
            θ_hat = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,α0) 
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
        
        xbar_k = update_xbar_k(x,rtik)
        if debugme
            push!(debug_val["xbar_k_beforeNorm"],xbar_k)
        end
        
        # xbar_k = 1 ./ Nk .* xbar_k
        if debugme
            push!(debug_val["xbar_k_afterNorm"],xbar_k)
        end
        
        sk = update_sk(x,xbar_k,rtik)
        if debugme
            push!(debug_val["sk_beforeNorm"],sk)
        end
        
        # sk = 1 ./ Nk .* sk
        if debugme
            push!(debug_val["sk_afterNorm"],sk)
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        mk_hat_vec= update_mk_hat(λ0_vec,μ0_vec, Nk,xbar_k)

        a0k_hat_vec =  update_a0k_hat(a0_vec,Nk)
        b0k_hat_vec = update_b0k_hat(b0_vec,λ0_vec,μ0_vec, Nk,xbar_k,sk)

        Tk = update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk;optim_max_iter=1000)
        data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,γ,α0,Tk)
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
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_[iter] = data_elbo + assgn_entropy  + HDP_surragate_elbo
        
    end
    
    return elbo_, rtik,θ_hat, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,debug_val
end
