"""
        variational_inference_dynamicHDP_vshoff_mpu()
    This is the main variational inference function. It perfroms coordinate ascent to infer the paramters of the NCLUSION Model

"""
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
    iter = 1
    converged_bool = false
    is_converged = false
    
    update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    ηk_trend_vec = []
    while !converged_bool
        update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_Nk_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        update_c_ttprime_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_yjk_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);
        update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams;mt_mode = mt_mode)
        update_d_hat_sum_mpu!(conditionparams,dataparams;mt_mode = mt_mode)
        update_st_hat_mpu!(conditionparams,dataparams,modelparams;mt_mode = mt_mode) 
        update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams;mt_mode = mt_mode)
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams;mt_mode = mt_mode)
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
        update_κk_hat_mpu!(clusters, dataparams,modelparams;mt_mode = mt_mode)
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
    outputs_dict = OrderedDict{Symbol,Any}();
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog);
    Tk_,ηk_trend_vec_ = Tk,ηk_trend_vec;
    output_str_list2 = @name Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ;
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,chain_dict,is_converged,truncation_value,ηk_trend_vec_ ];
    addToDict!(outputs_dict,output_key_list2,output_var_list2);
    return outputs_dict
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

"""
        extract_cluster_paramter(paramname,clusters,modelparams)
    This function extracts the cluster specific parameters from the ClusterFeatures object used in inference
"""
function extract_cluster_paramter(paramname,clusters,modelparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    K = modelparams.K

    param =[ isone(length(getfield(clusters[k],paramname))) ? getfield(clusters[k],paramname)[1] : getfield(clusters[k],paramname) for k in 1:K] 
    return param
end

"""
        extract_condition_paramter(paramname,conditionparams,dataparams)
    This function extracts the condition specific parameters from the ConditionFeatures object used in inference
"""
function extract_condition_paramter(paramname,conditionparams,dataparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    T = dataparams.T

    param = [ isone(length(getfield(conditionparams[t],paramname))) ? getfield(conditionparams[t],paramname)[1] : getfield(conditionparams[t],paramname) for t in 1:T]
    return param
end

"""
        extract_gene_paramter(paramname,geneparams,dataparams)
    This function extracts the gene specific parameters from the GeneFeatures object used in inference
"""
function extract_gene_paramter(paramname,geneparams,dataparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    G = dataparams.G

    param = [ isone(length(getfield(geneparams[j],paramname))) ? getfield(geneparams[j],paramname)[1] : getfield(geneparams[j],paramname) for j in 1:G]
    return param
end

"""
        extract_elbo_vals_perK(paramname,elbolog)
    This function extracts the cluster specific elbo values from the ElboFeatures object used in inference
"""
function extract_elbo_vals_perK(paramname,elbolog)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    # G = dataparams.G

    vals = getfield(elbolog,paramname)
    return vals
end

"""
        extract_rtik_paramter(cellpop,dataparams)
    This function extracts the cell-level cluster probability vector for each cell in the CellFeatures object
"""
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

"""
        extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
    This function extracts all parameters from custom objects and adds them to the previously instantiated output dictionary.
"""
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

"""
        extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog)
    This function extracts all parameters from custom objects and adds them to the previously instantiated output dictionary.
"""
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
