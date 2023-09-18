function f0(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    # bdot_w = 1
    
    input_str_list = @name x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter_tune, num_local_iter_tune];


    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end

function f12(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f14(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end

function f15(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f16(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    # bdot_w = 1

    input_str_list = @name x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter_tune, num_local_iter_tune];


    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end


function f17(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f18(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    # bdot_w = 1

    input_str_list = @name x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter_tune, num_local_iter_tune];


    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f19(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f20(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f21(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f22(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f23(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f24(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    a0_err = a0#1*10^(0.);
    # b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    # λ0_err = 1*10^(0.);
    # bdot_w = 1

    ηkj_prior = initialize_η_kj_prior(x_opt,KMax; pct_important=pct_important);
    # x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior
    input_str_list = @name x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior, num_iter_tune, num_local_iter_tune];

    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f25(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w,null_precision)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    # bdot_w = 1

    input_str_list = @name x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune];


    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f25_fast3series(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w,null_precision)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    # bdot_w = 1

    input_str_list = @name x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune];


    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_, pip_kj_, a0k_hat_vec_, b0k_hat_vec_, mk_hat_vec_, λ0k_hat_vec_, Nk_, rhok_hat_vec_, omegak_hat_vec_, θ_hat_vec_, c_ttprime_vec_, a_αt_hat_vec_, b_αt_hat_vec_, awt_hat_vec_, bwt_hat_vec_, a_γ_hat_, b_γ_hat_, e_γ_, Tαk_, chain_dict, initDict, is_converged, truncation_value  = (; outputs_dict...);


    num_posterior_samples = 1000;
    try 
        vi_make_z_post_s(rtik_, S=num_posterior_samples);
    catch e
        return 0
    end
    z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
    ari_RandIndices_summary_invar = calc_time_invariant_ARI_summarization(time_invariant_ari(z_opt,z_post_s))
    return ari_RandIndices_summary_invar[1,2]
end
function f25_fast3seriesCH(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w,null_precision)
    G = length(x_opt[1][1])
    μ0= mean([el[1]  for el in mean.(x_opt,dims=1)],dims = 1)[1];
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    # bdot_w = 1

    input_str_list = @name x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [x_opt, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune];


    inputs_dict = OrderedDict()
    addToDict!(inputs_dict,input_key_list,input_var_list);
    outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);
    
    elbo_, rtik_, pip_kj_, a0k_hat_vec_, b0k_hat_vec_, mk_hat_vec_, λ0k_hat_vec_, Nk_, rhok_hat_vec_, omegak_hat_vec_, θ_hat_vec_, c_ttprime_vec_, a_αt_hat_vec_, b_αt_hat_vec_, awt_hat_vec_, bwt_hat_vec_, a_γ_hat_, b_γ_hat_, e_γ_, Tαk_, chain_dict, initDict, is_converged, truncation_value  = (; outputs_dict...);

    xmat = hcat(vcat(x_opt...)...)
    z_arg_max = [argmax.(el) for el in rtik_]
    # cvi_criterion_value_b_vec = Vector{Float64}(undef, num_posterior_samples)
    # cvi  = CH()
    cvi_b = CH()
    z_infer_vec = vcat(z_arg_max...)
    cvi_criterion_value_b_vec = get_cvi!(cvi_b, xmat, z_infer_vec)
    return cvi_criterion_value_b_vec
end
function f25_fast3seriesCH_kfolds(vi_func,num_iter_tune,x_opt,z_opt,KMax,a0,b0,λ0,b_γ,b_α,bdot_w,null_precision;min_pct_subset=0.75,num_fold=4,seed=12345)
    num_local_iter_tune = 1;
    # a0 = 1
    a_γ = 1
    a_α = 1
    adot_w =1
    

    println("Subsampling data...")
    tune_fold_splits = Vector{Vector{Bool}}(undef,num_fold)
    test_fold_splits = Vector{Vector{Bool}}(undef,num_fold)
    ch_vec = Vector{Float64}(undef,num_fold)
    for fold in 1:num_fold
        tune_subsampled_cells_bool, test_subsampled_cells_bool= subsample_tuning_set(vcat(z_opt...);min_pct_subset= min_pct_subset);
        tune_fold_splits[fold] = tune_subsampled_cells_bool
        test_fold_splits[fold] = test_subsampled_cells_bool
        _,data_dict_tune,_ = generate_subsample_data_processed_1condition_anndata_data(x_opt,z_opt,1111,tune_subsampled_cells_bool;cell_timepoints_index=nothing,used_x="x_shifted",is_tuning_set = true,seed = seed);
        x_fold,_,z_tune,_,_,π_tune = (; data_dict_tune...);
        μ0= mean([el[1]  for el in mean.(x_fold,dims=1)],dims = 1)[1];
        input_str_list = @name x_fold, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune;
        input_key_list = Symbol.(naming_vec(input_str_list));
        input_var_list = [x_fold, KMax, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w,null_precision, num_iter_tune, num_local_iter_tune];


        inputs_dict = OrderedDict()
        addToDict!(inputs_dict,input_key_list,input_var_list);
        xmat = hcat(vcat(x_fold...)...)
        println("Tuning data on fold $fold on $(sum(tune_subsampled_cells_bool)) cells...")

        outputs_dict = vi_func(inputs_dict;uniform_theta_init = true, rand_init = false);

        elbo_, rtik_, pip_kj_, a0k_hat_vec_, b0k_hat_vec_, mk_hat_vec_, λ0k_hat_vec_, Nk_, rhok_hat_vec_, omegak_hat_vec_, θ_hat_vec_, c_ttprime_vec_, a_αt_hat_vec_, b_αt_hat_vec_, awt_hat_vec_, bwt_hat_vec_, a_γ_hat_, b_γ_hat_, e_γ_, Tαk_, chain_dict, initDict, is_converged, truncation_value  = (; outputs_dict...);

        xmat = hcat(vcat(x_opt...)...)
        z_arg_max = [argmax.(el) for el in rtik_]
        # cvi_criterion_value_b_vec = Vector{Float64}(undef, num_posterior_samples)
        # cvi  = CH()
        cvi_b = CH()
        z_infer_vec = vcat(z_arg_max...)
        cvi_criterion_value_b_vec = get_cvi!(cvi_b, xmat, z_infer_vec)
        # ch_summary_invar = calc_time_invariant_CVI_summarization([cvi_criterion_value_b_vec];conf_level=0.95)
        ch_vec[fold]  = cvi_criterion_value_b_vec
    end
    return mean(ch_vec)
end


function select_and_run_optimization_procedure(vi_function_used,num_iter_tune,x_tune,z_tune,model_indx,avail_models,avail_opt_models;num_opt_samples = 1000)
    ho = nothing
    G= length(x_tune[1][1])
    if vi_function_used == "original"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K = collect(5:25),
            a0 = collect(1:1),
            b0 = exp.(LinRange(-6,0,1000)),
            λ0  = exp.(LinRange(-8,0,1000)),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1)#exp.(LinRange(-9,9,100))
        print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5), ", λ0: ", round(λ0; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,λ0,b_γ,b_α,bdot_w)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,λ0,b_γ,b_α,bdot_w): $best_params")
        println("ARI @ Best K: $max_f")
    elseif  vi_function_used == "vs12"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = exp.(LinRange(-8,0,1000)),
            λ0_err  = exp.(LinRange(-8,0,1000)),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif  vi_function_used == "vs14"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = exp.(LinRange(-8,0,1000)),
            λ0_err  = exp.(LinRange(-8,0,1000)),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif  vi_function_used == "vs15"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")

    elseif vi_function_used == "vs17"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif vi_function_used == "vs19"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif vi_function_used == "vs20"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif vi_function_used == "vs21"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif vi_function_used == "vs22"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif vi_function_used == "vs23"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif vi_function_used == "vs24"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K =  collect(5:25),
            a0 = collect(1:1),#exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-8,0,1000)),
            b0_err = exp.(LinRange(-8,0,1000)),
            λ0  = collect(1:1),
            λ0_err  = collect(1:1),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            pct_important = LinRange(0,1,1000)
            print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5),", b0_err: ", round(b0_err; digits = 5), ", λ0: ", round(λ0; digits = 5), ", λ0_err: ", round(λ0_err; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", %imp: ", round(pct_important; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,b0_err,λ0,λ0_err,b_γ,b_α,bdot_w,pct_important): $best_params")
        println("ARI @ Best K: $max_f")
    elseif  vi_function_used == "vs16"  
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K = collect(5:25),
            a0 = exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-6,6,1000)),
            λ0  = exp.(LinRange(-8,0,1000)),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1)#exp.(LinRange(-9,9,100))
        print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5), ", λ0: ", round(λ0; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,λ0,b_γ,b_α,bdot_w)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,λ0,b_γ,b_α,bdot_w): $best_params")
        println("ARI @ Best K: $max_f")
    elseif   vi_function_used == "vs18"  
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K = collect(5:25),
            a0 = exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-6,6,1000)),
            λ0  = exp.(LinRange(-8,0,1000)),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1)#exp.(LinRange(-9,9,100))
        print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5), ", λ0: ", round(λ0; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,λ0,b_γ,b_α,bdot_w)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,λ0,b_γ,b_α,bdot_w): $best_params")
        println("ARI @ Best K: $max_f")
    
    elseif   vi_function_used == "vs25"   ||  vi_function_used == "vs25_fast" ||  vi_function_used == "vs25_fast2" || vi_function_used == "vs25_fast3" ||  vi_function_used == "vs25_fast3_mtall" ||  vi_function_used == "vs25_fast3_mtbest"
        ho = @hyperopt for i=num_opt_samples,
            sampler = RandomSampler(), # This is default if none provided
            K = collect(5:25),
            a0 = exp.(LinRange(-6,6,1000)),
            b0 = exp.(LinRange(-6,6,1000)),
            λ0  = exp.(LinRange(-8,0,1000)),
            b_γ  = exp.(LinRange(-6,6,1000)),
            b_α  = exp.(LinRange(-6,6,1000)),
            bdot_w  = collect(1:1),#exp.(LinRange(-9,9,100))
            null_precision  = [0.1,1.0,10.0,100.0]
        print("\n Itr: ", i, ", K: ", K, ", a0: ", round(a0; digits = 5), ", b0: ", round(b0; digits = 5), ", λ0: ", round(λ0; digits = 5), ", b_γ: ", round(b_γ; digits = 5), ", b_α: ", round(b_α; digits = 5), ", bdot_w: ", round(bdot_w; digits = 5),", null_precision: ", round(null_precision; digits = 5), " \t")
        vi_func = avail_models[model_indx]
        f = avail_opt_models[model_indx]
        @show f(vi_func,num_iter_tune,x_tune,z_tune,K,a0,b0,λ0,b_γ,b_α,bdot_w,null_precision)
        end
        best_params, max_f = ho.maximizer , ho.maximum
        println("Best (K,a0,b0,λ0,b_γ,b_α,bdot_w,null_precision): $best_params")
        println("ARI @ Best K: $max_f")

    
    
    end

    
    vi_func = avail_models[model_indx]
   return ho,vi_func
end

function setup_avail_models()
    avail_models= [variational_inference_dynamicHDP_dev,variational_inference_dynamicHDP_vs12,variational_inference_dynamicHDP_vs14,variational_inference_dynamicHDP_vs15,variational_inference_dynamicHDP_vs16,variational_inference_dynamicHDP_vs17,variational_inference_dynamicHDP_vs18,variational_inference_dynamicHDP_vs19,variational_inference_dynamicHDP_vs20,variational_inference_dynamicHDP_vs21,variational_inference_dynamicHDP_vs22,variational_inference_dynamicHDP_vs23,variational_inference_dynamicHDP_vs24,variational_inference_dynamicHDP_vs25,variational_inference_dynamicHDP_vs25_fast,variational_inference_dynamicHDP_vs25_fast2,variational_inference_dynamicHDP_vs25_fast3,variational_inference_dynamicHDP_vs25_fast3_mtall,variational_inference_dynamicHDP_vs25_fast3_mtbest,variational_inference_dynamicHDP_vs25_fast3_mtbest,variational_inference_dynamicHDP_vs25_fast3_mtbest] 
    return avail_models
end
function setup_avail_opt_functions()
    avail_opt_models= [f0,f12,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f25,f25,f25_fast3series,f25_fast3series,f25_fast3series,f25_fast3seriesCH,f25_fast3seriesCH_kfolds]
    return avail_opt_models
end
function setup_avail_model_names()
    avail_models_names = ["variational_inference_dynamicHDP_dev","variational_inference_dynamicHDP_vs12","variational_inference_dynamicHDP_vs14","variational_inference_dynamicHDP_vs15","variational_inference_dynamicHDP_vs16","variational_inference_dynamicHDP_vs17","variational_inference_dynamicHDP_vs18","variational_inference_dynamicHDP_vs19","variational_inference_dynamicHDP_vs20","variational_inference_dynamicHDP_vs21","variational_inference_dynamicHDP_vs22","variational_inference_dynamicHDP_vs23","variational_inference_dynamicHDP_vs24","variational_inference_dynamicHDP_vs25","variational_inference_dynamicHDP_vs25_fast","variational_inference_dynamicHDP_vs25_fast2","variational_inference_dynamicHDP_vs25_fast3","variational_inference_dynamicHDP_vs25_fast3_mtall","variational_inference_dynamicHDP_vs25_fast3_mtbest","variational_inference_dynamicHDP_vs25_fast3_mtbest","variational_inference_dynamicHDP_vs25_fast3_mtbest"]
    return avail_models_names
end
function setup_avail_opt_functions_names()
    avail_opt_models_names =  ["f0","f12","f14","f15","f16","f17","f18","f19","f20","f21","f22","f23","f24","f25","f25","f25","f25_fast3series","f25_fast3series","f25_fast3series","f25_fast3seriesCH","f25_fast3seriesCH_kfolds"]
    return avail_opt_models_names
end
function setup_model_indexing()
    model_indx_dict = OrderedDict("original"=>1,"vs12"=>2,"vs14"=>3,"vs15"=>4,"vs16"=>5,"vs17"=>6,"vs18"=>7,"vs19"=>8,"vs20"=>9,"vs21"=>10,"vs22"=>11,"vs23"=>12,"vs24"=>13,"vs25"=>14,"vs25_fast"=>15,"vs25_fast2"=>16,"vs25_fast3"=>17,"vs25_fast3_mtall"=>18,"vs25_fast3_mtbest"=>19,"vs25_fast3_mtbest_CHmp"=>20,"vs25_fast3_mtbest_CHmp_kfold"=>21)#,"f7"=>3,"f7_2"=>4,"f12_2"=>5
    




    return model_indx_dict
end
function setup_dataset_dict()
    dataset_dict = OrderedDict(1=>[[10.0,0.0,0.0,0.0,0.0], [0.0,10.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0,10.0] ],
    2=>[[10.0,0.0,10.0,0.0,0.0], [0.0,10.0,0.10,0.0,0.0], [0.0,0.0,0.0,0.0,10.0] ] ,
    3=>[[10.0,0.0,10.0,0.0,0.0], [0.0,10.0,0.10,0.0,0.0], [0.0,0.0,0.0,10.0,10.0] ],
    4=>[[10.0,0.0,10.0,0.0,0.0], [0.0,10.0,10.0,0.0,0.0], [0.0,0.0,0.0,20.0,10.0] ],
    5=>[[100.0,0.0,0.0,0.0,0.0], [0.0,100.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0,100.0] ],
    6=>[[10.0,0.0,10.0,0.0,0.0], [0.0,10.0,0.0,10.0,0.0], [0.0,0.0,0.0,0.0,10.0] ] ,
    7=>[[10.0,0.0,10.0,0.0,0.0], [0.0,10.0,0.0,10.0,0.0], [0.0,0.0,0.0,10.0,10.0] ],
    8=> [[10.0,0.0,10.0,0.0,0.0], [0.0,10.0,0.0,10.0,0.0], [0.0,0.0,0.0,20.0,10.0] ]
    )
    return dataset_dict
end
function setup_experiment_tag(experiment_filename)
    return "EXPERIMENT_$experiment_filename"
end

vi_closure_fn(x,G,λ0,μ0,a0,b0,num_iter, num_local_iter) = (γ,α0,kmax) -> variational_inference(x, G,kmax,γ,α0,λ0,μ0,a0,b0,num_iter, num_local_iter)
# (GD, x)  ->  g_constrained(GD, x, G,γ,α_0,T_k)
# variational_inference(x, G,kmax,γ,α0,λ0,μ0,a0,b0,num_iter, num_local_iter)


function tuningSqUtilizationRatioLoss(x,Ktrue, G,λ0,μ0,a0,b0,num_iter, num_local_iter)
    vi_closure = vi_closure_fn(x,G,λ0,μ0,a0,b0,num_iter, num_local_iter)
    f  
    fff3(a,b,kmax) = nonsignleton_cluster_cost( K ,update_Nk(variational_inference(x, G,kmax,a,b,1,1,1,1, 50,num_local_iter)[2]))
    # Main macro. The first argument to the for loop is always interpreted as the number of iterations (except for hyperband optimizer)
    ho = @hyperopt for i=100,
        sampler = RandomSampler(), # This is default if none provided
        γ = [0.5, 1, 10, 100],
        α0 = [0.5, 1, 10, 100],
        kmax = [3, 5, 9, 20, 100, 200],
        # λ0= [ 1, 10],
        # μ0= [ 1, 10],
        # a0= [ 1, 10],
        # b0= [ 1, 10],
        # print(sti(i) "\t", a, "\t", b, "\t", c,  "   \t")
        # println("$i \t $γ_")
        # println( "$i \t $γ_ \t $α0_ \t $λ0_ \t $μ0_ \t $a0_ \t $b0_\t")
        x=x
        @show fff3(γ,α0,kmax)
    # @show ff(data,γ_,α0_,λ0_,μ0_,a0_,b0_, G=G,Kmax=K, num_iter=num_iter, num_local_iter=num_local_iter)
    end
end

nonsignleton_cluster_cost( ktrue ,util_vec, thresh=2) = (1-sum(util_vec .> thresh)/ktrue)^2
l1norm(actual,theorectical) = sum(x -> abs(x[1] - x[2]), zip(sort(actual), theorectical))
function theoretical_utilization(Kmax,Ktrue; π_ = nothing)
    if !isnothing(π_)
        π_ = sort(π_)
    else
        π_ = 1 /Ktrue .* ones(Ktrue)
    end
    theoretical_contribution = [ k <= K ? π_[k]  : 0.0 for k in 1:Kmax]

    return theoretical_contribution
end


function norm_weights(p)
    psum = StatsFuns.logsumexp(p)
    w = exp.(p .- psum)
    return w
end

function normToProb(p)
    psum = sum(p)
    w = p ./ psum
    return w
end