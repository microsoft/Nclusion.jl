
"""
        cavi()
    This is the main variational inference function. It perfroms coordinate ascent to infer the paramters of the NCLUSION Model

"""
function cavi(inputs;elbo_ep = 10^(-6),update_η_bool= false)

    cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,Tk,elbolog  = (; inputs...);
    iter = 1
    converged_bool = false
    is_converged = false
    num_iter=modelparams.num_iter
    # update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    ηk_trend_vec = []
    while !converged_bool
        # E-STEP
        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams) #1 
        update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams) #2
        update_yjk_mpu!(clusters,geneparams,dataparams,modelparams)#3

        update_c_ttprime_mpu!(conditionparams,dataparams,modelparams) # DONT NEED TIME
        update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams) # DONT NEED TIME

        update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams)#4
        update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams) #DETERMINISTIC
        update_Nk_mpu!(cellpop,clusters,dataparams,modelparams) #DETERMINISTIC
        update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams) #DETERMINISTIC
        update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams) #DETERMINISTIC


        # Calculate ELBO
        iter = Int64(iter)
        elbo_iter,elbolog =  calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
        elbolog.elbo_[iter] = elbo_iter
        
        # M-STEP
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams) # TO REMOVE #DETERMINISTIC
        update_κk_hat_mpu!(clusters, dataparams,modelparams) # TO REMOVE #DETERMINISTIC
        update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams)#5

        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams) # TO REMOVE #DETERMINISTIC
        update_κk_hat_mpu!(clusters, dataparams,modelparams) # TO REMOVE #DETERMINISTIC

        update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams)#6

        update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
        update_var_muk_hat_mpu!(clusters, dataparams,modelparams) # TO REMOVE #DETERMINISTIC
        update_κk_hat_mpu!(clusters, dataparams,modelparams) # TO REMOVE #DETERMINISTIC

        update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=10000);#7
        update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams) #8 
        update_d_hat_sum_mpu!(conditionparams,dataparams) #DETERMINISTIC

        update_st_hat_mpu!(conditionparams,dataparams,modelparams) # DONT NEED TIME
        
        
        
        
        
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
    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];
    outputs_dict = OrderedDict{Symbol,Any}();
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog);
    Tk_,ηk_trend_vec_ = Tk,ηk_trend_vec;
    output_str_list2 = @name Tk_,is_converged,truncation_value,ηk_trend_vec_ ;
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [Tk_,is_converged,truncation_value,ηk_trend_vec_ ];
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
