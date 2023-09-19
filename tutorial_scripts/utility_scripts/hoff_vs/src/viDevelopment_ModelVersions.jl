
# seed = 1989
# Random.seed!(seed)
# unique_time_id = get_unique_time_id()
# KTrue = 3
# K = 4#KTrue
# unique_clusters = ["Cluster $k" for k in 1:KTrue]
# T = 1
# G = 5
# CperK =  100*ones(Int,T) 
# C_t = KTrue .* CperK
# μ_magnitude = 10.0
# de_magnitude = 5
# percent_important = 0.3
# num_important_features = round(Int,G*percent_important)
# use_dHDP = true
# corr_feat = false
# one_corr_mat = true
# μk = [zeros(G), generate_DEG_means(G,collect(1:num_important_features),1,de_magnitude), generate_DEG_means(G,collect(G-num_important_features+1:G),-de_magnitude,-1)] ;
# Σk= generate_rand_covariance_matrix(G, KTrue;corr_feat = corr_feat,one_corr_mat = one_corr_mat);

# x,z,π_ = fake_mvGausssian_corrdata_for_testing(C_t,KTrue,μk ,Σk; mix_prob =nothing,same_prob_t = true,dynamic = false);
# x_to_use = raghavan2021_centerAndScale(x) ;

# x_input = x_to_use

# a0 = 32;
# b0 = 36; 
# μ0= mean([el[1]  for el in mean.(x_input,dims=1)],dims = 1)[1];
# λ0 = 1*10^(-0);

# a0_err = 32;
# b0_err = 36;#1*10^(Float64(b_exp*err_prec));#
# μ0_err =  vec(median(permutedims(reduce(hcat,reduce(vcat,x_input))),dims=1))#zeros(G)#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1))#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
# λ0_err = 1*10^(0.);


# a_γ = 0.028479847173306833;
# b_γ = 4.510944280834115;

# a_α = 0.0250047695905532;
# b_α = 20.045365862745044;
# a_η = 1.0
# b_η = a_η/pct_important -a_η 

# adot_w = 0.3146934916796939; 
# bdot_w = 4.205681866641926;
# num_iter = 1000;
# num_local_iter = 1;
# ηtkj_prior = init_ηtkj_prior(x_input,K; pct_important=0.3);

# # x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err, λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior;
# x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err, λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior;


# input_str_list = @name x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err, λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior,num_iter,num_local_iter;
# input_key_list = Symbol.(naming_vec(input_str_list));
# input_var_list = [x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err, λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior,num_iter,num_local_iter];

# x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err, λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηtkj_prior;

## x_input, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter, num_local_iter
# input_str_list = @name x_input, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter, num_local_iter;
# input_key_list = Symbol.(naming_vec(input_str_list));
# input_var_list = [x_input, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter, num_local_iter];


# inputs_dict = OrderedDict()
# addToDict!(inputs_dict,input_key_list,input_var_list);

# outputs_dict = variational_inference_dynamicHDP_dev(inputs_dict;mk_hat_init=mk_hat_init, λ0k_hat_init=λ0k_hat_init,a0k_hat_init=a0k_hat_init, b0k_hat_init=b0k_hat_init,awt_hat_init=awt_hat_init,  bwt_hat_init=bwt_hat_init,a_αt_hat_init=a_αt_hat_init, b_αt_hat_init=b_αt_hat_init,a_γ_hat_init=a_γ_hat_init, b_γ_hat_init=b_γ_hat_init,θ_hat_init=θ_hat_init,c_ttprime_init = c_ttprime_init, rhok_hat_init=rhok_hat_init, omegak_hat_init=omegak_hat_init,uniform_theta_init = uniform_theta_init, rand_init = rand_init,ep = 0.001,elbo_ep = 10^(-6),record_chain = false);

# elbo_, rtik_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value = (; outputs_dict...);
# KMax, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain = (; initDict...);

# true_z = z;
# num_posterior_samples=1000;
# z_post_s = vi_make_z_post_s(rtik_, S=num_posterior_samples);
# vmeasure_vov = getVmeasure(true_z,z_post_s;beta=1.0);
# ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov  = getRandIndices(true_z,z_post_s);
# nmi_vov = getNMI(true_z,z_post_s);
# avg_counts_mat = get_average_posterior_cluster_frequency2(z_post_s,T,true_z,KMax,KTrue,C_t);
# mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_;Ns= 1000);
# gnu_plot_All_Time_Points_PCA(x_to_use,z_post_s[end],G,unique(vcat(z_post_s[end]...));unique_clusters_labels = unique(vcat(z_post_s[end]...)))
# gnu_elbo_plot(elbo_;fig_size=(1100,900))
function generate_varselect_inputs_for_testing(x_to_use,K)
    T = length(x_to_use)
    C_t = length.(x_to_use)
    G = length(x_to_use[1][1])
    KMax = K
    b_exp = -2
    err_prec = 1.0 +10^(-100) #-1.01225  -2.95 -2.00050
    a0 = 306.04242#0.08574#1*10^(0.);
    b0 = 136.85458#0.09668#1*10^(Float64(b_exp)); 
    μ0= mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1)) #mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];
    λ0 = 1.34218#1.01818#1*10^(0.);
    a_γ = 1
    a_α = 1
    adot_w =1
    bdot_w = 1
    b_γ=0.03567#0.71869
    b_α= 143.59071#0.00297
    a0_err = a0#1*10^(0.);
    b0_err = b0#1*10^(Float64(b_exp*err_prec));#1*10^(Float64(-4));#
    μ0_err = zeros(G)#[mode(permutedims(reduce(hcat,reduce(vcat,x_to_use)))[:,i]) for i in 1:G]#vec(median(permutedims(reduce(hcat,reduce(vcat,x_to_use))),dims=1)#zeros(G)#mean([el[1]  for el in mean.(x_to_use,dims=1)],dims = 1)[1];#zeros(G)#
    λ0_err = 1*10^(0.);
    # x_to_use = x_mask;
    KMax = K
    ηkj_prior = initialize_η_kj_prior(x_to_use,KMax; pct_important=0.3);
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0);
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

    ηkj_prior = initialize_η_kj_prior(x_to_use,KMax; pct_important=3/G);
    λ0k_hat_vec_init = [λ0_vec for k in 1:KMax]; # 
    mk_hat_vec_init = [μ0_vec for k in 1:KMax];
    a0k_hat_vec_init = [a0_vec for k in 1:KMax]; #
    b0k_hat_vec_init =  [b0_vec for k in 1:KMax]; #
    rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(KMax);
    a_γ_hat_init = a_γ;
    b_γ_hat_init = b_γ;
    awt_hat_vec_init = [adot_w for t in 1:T];
    bwt_hat_vec_init =  [bdot_w for t in 1:T];
    a_αt_hat_vec_init = [a_α for t in 1:T];
    b_αt_hat_vec_init = [b_α for t in 1:T];
    c_ttprime_vec_init =  [ones(T) ./T  for t in 1:T];
    θ_hat_vec_init = [ones(KMax+1) ./(KMax+1)  for t in 1:T];
    v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:KMax] for i in 1:C_t[t]] for t in 1:T];
    λ0_err_hat_vec_init = λ0_err_vec;
    m_err_hat_vec_init =μ0_err_vec;
    a0_err_hat_vec_init = a0_err_vec;
    b0_err_hat_vec_init = b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    num_iter, num_local_iter = 30, 1;
    # η_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    rtik_init = [[ones(KMax) ./KMax for i in 1:C_t[t]] for t in 1:T];
    return λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,μ0_err,λ0_err,a0_err,b0_err,ηkj_prior,num_local_iter,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,v_tikj_vec_init,λ0_err_hat_vec_init,m_err_hat_vec_init,a0_err_hat_vec_init,b0_err_hat_vec_init,rtik_init
end
function make_chain(num_iter,value_)
    value_chain_ = Vector{Union{Nothing,Missing,typeof(value_)}}(undef,num_iter)
    return value_chain_
end
function truncate_chain(full_chain_dict,truncation_value)
    iteration_min = 1
    iteration_ids = collect(iteration_min:truncation_value)
    trunc_chain_dict = OrderedDict()
    for key in keys(full_chain_dict)
        trunc_chain_dict[key] = full_chain_dict[key][iteration_ids]
    end
    return trunc_chain_dict
end
function initialize_η_kj_prior(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    η_prior = [[pct_important for j in 1:G] for k in 1:K]
    return η_prior
end
function log_ηjk_expected_value(a_ηkj,b_ηkj)
    K = length(a_ηkj)
    e_log_ηjk = Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        e_log_ηjk[k]= digamma.(a_ηkj[k]) .- digamma.(a_ηkj[k] .+ b_ηkj[k])
    end
    return e_log_ηjk
end
function log1minusηjk_expected_value(a_ηkj,b_ηkj)
    K = length(a_ηkj)
    e_log_minusηjk = Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        e_log_minusηjk[k]= digamma.(b_ηkj[k]) .- digamma.(a_ηkj[k] .+ b_ηkj[k])
    end
    return e_log_minusηjk
end
function update_a_ηkj(v_tikj,a_η)
    T = length(v_tikj)
    C_t = [length(el) for el in v_tikj]
    K = length(v_tikj[1][1])
    G = length(v_tikj[1][1][1])
    a_ηkj = Vector{Vector{Float64}}(undef,K)
    v_true =  [[[[v_tikj[t][i][k][j][1] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    v_kj = sum(sum(v_true))
    for k in 1:K
        a_ηkj[k] = v_kj[k] .+ a_η  
    end
    return a_ηkj
end
function update_b_ηkj(v_tikj,b_η)
    T = length(v_tikj)
    C_t = [length(el) for el in v_tikj]
    K = length(v_tikj[1][1])
    G = length(v_tikj[1][1][1])
    b_ηkj = Vector{Vector{Float64}}(undef,K)
    v_false =  [[[[v_tikj[t][i][k][j][2] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    minusv_kj = sum(sum(v_false))
    for k in 1:K
        b_ηkj[k] = minusv_kj[k] .+ b_η 
    end
    return b_ηkj
end
function variational_inference_dynamicHDP_dev(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


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
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        Nk_chain = make_chain(num_iter+1,Nk)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)

        
        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec

        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:Nk_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        Nk = update_Nk(rtik)
        x_hat_k = update_x_hat_k(x,rtik)
        x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        if record_chain
            chain_dict[:Nk_chain][iter + 1] = Nk
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
        a0k_hat_vec = update_a0k_hat_usingXhat(a0_vec,Nk)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec
        end

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        
        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
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
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end


function update_rtik_vs7(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,v_tikj,c_ttprime)
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
                    v_true = v_tikj[t][i][k][j][1]
                    v_false = v_tikj[t][i][k][j][2]
                    log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j]) + 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j]) #+ 2*((e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j]) - ( e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j]))
                    #0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j]) + 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
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
function update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#+ log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#+ log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    

                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function variational_inference_dynamicHDP_vs7(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,μ0_err,λ0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,λ0_err_hat_vec_init=nothing,m_err_hat_vec_init=nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

     
    v_tikj = v_tikj_vec_init
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
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        debug_val["λ0_err_hat_vec"] = []
        debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_μj_err"]= []
        debug_val["e_τ_μ_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_μj_err"],[])
        push!(debug_val["e_τ_μ_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

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

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs7(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,v_tikj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)

            

            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_μj_err"],e_τ_μj_err)
                push!(debug_val["e_τ_μ_err"],e_τ_μ_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        

        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end





        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_μj_err,_  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj7(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_μj_err,ηkj_prior);
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs7(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,m_err_hat_init=nothing, λ0_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    m_err_hat_vec= m_err_hat_init 
    λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,m_err_hat_init,λ0_err_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, v_tikj_init,num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,m_err_hat_init,λ0_err_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, v_tikj_init,num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μj_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_μ_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec_temp  = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec_temp  = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik_temp = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        # Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        # a_αt_hat_vec = Vector{Float64}()
        # b_αt_hat_vec = Vector{Float64}()
        # awt_hat_vec = Vector{Float64}()
        # bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        # v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec_temp )
        rtik_chain = make_chain(num_iter+1,rtik_temp )
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_μj_err_chain=make_chain(num_iter+1,e_τ_μj_err)
        e_τ_μ_err_chain=make_chain(num_iter+1,e_τ_μ_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec_temp )
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_μj_err_chain,e_τ_μ_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,λ0_err_chain,m_err_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_μj_err_chain,e_τ_μ_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,λ0_err_chain,m_err_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat_vec
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_μj_err_chain][1] = nothing
        chain_dict[:e_τ_μ_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
    Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            # println(rhok_hat_vec)
            # println(omegak_hat_vec)
            # println(Ntk)
            # println(a_αt_hat_vec)
            # println(b_αt_hat_vec)
            # println(c_ttprime_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                println("1")
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat_vec
                println("2")
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                println("3")
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                println("4")
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                println("5")
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                println("6")
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                println("7")
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                println("8")
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                println("9")
                chain_dict[:e_τ_μj_err_chain][iter + 1][loc_iter] = e_τ_μj_err
                println("10")
                chain_dict[:e_τ_μ_err_chain][iter + 1][loc_iter] = e_τ_μ_err
                println("12")
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        
        # λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        # mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
        # a0k_hat_vec = update_a0k_hat_usingXhat(a0_vec,Nk)
        # b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)



        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_μj_err,_  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj7(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_μj_err,ηkj_prior);
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end


        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo7(x,rtik,v_tikj,mk_hat_vec,μ0_vec,μ0_err_vec,m_err_hat_vec,λ0k_hat_vec,λ0_vec,λ0_err_vec,λ0_err_hat_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        
        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,λ0_err_hat_vec_, m_err_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,λ0_err_hat_vec_, m_err_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_, rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,λ0_err_hat_vec_, m_err_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

function variational_inference_dynamicHDP_vs7_2(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,m_err_hat_init=nothing, λ0_err_hat_init=nothing,a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,a_η,b_η,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    m_err_hat_vec= m_err_hat_init 
    λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,m_err_hat_init,λ0_err_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, v_tikj_init,num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,m_err_hat_init,λ0_err_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, v_tikj_init,num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μj_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_μ_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec_temp  = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec_temp  = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik_temp = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        # Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        # a_αt_hat_vec = Vector{Float64}()
        # b_αt_hat_vec = Vector{Float64}()
        # awt_hat_vec = Vector{Float64}()
        # bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        # v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec_temp )
        rtik_chain = make_chain(num_iter+1,rtik_temp )
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_μj_err_chain=make_chain(num_iter+1,e_τ_μj_err)
        e_τ_μ_err_chain=make_chain(num_iter+1,e_τ_μ_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec_temp )
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_μj_err_chain,e_τ_μ_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,λ0_err_chain,m_err_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_μj_err_chain,e_τ_μ_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,λ0_err_chain,m_err_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat_vec
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_μj_err_chain][1] = nothing
        chain_dict[:e_τ_μ_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
    Glog = G*log(2π)

    e_log_ηkj = [log.(el) for el in ηkj_prior]
    e_log_minus_ηkj = [log.(1 .- el) for el in ηkj_prior]
    a_ηkj_hat = 1
    b_ηkj_hat = 1
    while !converged_bool #for iter in 1:num_iter
        # println(iter)
        for loc_iter in 1:num_local_iter
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_μj_err,e_τ_μ_err  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            # println(rhok_hat_vec)
            # println(omegak_hat_vec)
            # println(Ntk)
            # println(a_αt_hat_vec)
            # println(b_αt_hat_vec)
            # println(c_ttprime_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                println("1")
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat_vec
                println("2")
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                println("3")
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                println("4")
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                println("5")
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                println("6")
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                println("7")
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                println("8")
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                println("9")
                chain_dict[:e_τ_μj_err_chain][iter + 1][loc_iter] = e_τ_μj_err
                println("10")
                chain_dict[:e_τ_μ_err_chain][iter + 1][loc_iter] = e_τ_μ_err
                println("12")
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        
        # λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
        # mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
        # a0k_hat_vec = update_a0k_hat_usingXhat(a0_vec,Nk)
        # b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)



        a0_err_hat_vec = update_a0_err_hat_usingXhat(a0_err_vec,Nj_error)
        λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat(b0_err_vec,λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_μj_err,_  = τ_μ_error_expected_value(x,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj7(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj);
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end
        
        a_ηkj_hat = update_a_ηkj(v_tikj,a_η)
        b_ηkj_hat = update_b_ηkj(v_tikj,b_η)
        e_log_ηkj = log_ηjk_expected_value(a_ηkj_hat,b_ηkj_hat)
        e_log_minus_ηkj = log1minusηjk_expected_value(a_ηkj_hat,b_ηkj_hat)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo7(x,rtik,v_tikj,mk_hat_vec,μ0_vec,μ0_err_vec,m_err_hat_vec,λ0k_hat_vec,λ0_vec,λ0_err_vec,λ0_err_hat_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        imp_elbo = calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        
        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,λ0_err_hat_vec_, m_err_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,a_ηkj_hat_,b_ηkj_hat_ = elbo_, rtik,v_tikj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk,a_ηkj_hat,b_ηkj_hat 

    output_str_list = @name elbo_,rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,λ0_err_hat_vec_, m_err_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_, rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,a_ηkj_hat_,b_ηkj_hat_ ,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,λ0_err_hat_vec_, m_err_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,a_ηkj_hat_,b_ηkj_hat_ ,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

# variational_inference_dynamicHDP_vs7
## x_input, K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter
# input_str_list = @name x_input, K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter;
# input_key_list = Symbol.(naming_vec(input_str_list));
# inputs_dict = OrderedDict()
# addToDict!(inputs_dict,input_key_list,input_var_list);
# input_var_list = [x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter];
# out7 = variational_inference_dynamicHDP_vs7(inputs_dict;mk_hat_init=mk_hat_vec_init, λ0k_hat_init=λ0k_hat_vec_init,a0k_hat_init=a0k_hat_vec_init, b0k_hat_init=b0k_hat_vec_init,m_err_hat_init=m_err_hat_vec_init, λ0_err_hat_init=λ0_err_hat_vec_init,a0_err_hat_init=a0_err_hat_vec_init, b0_err_hat_init=b0_err_hat_vec_init,awt_hat_init=awt_hat_vec_init, bwt_hat_init=bwt_hat_vec_init,a_αt_hat_init=a_αt_hat_vec_init, b_αt_hat_init=b_αt_hat_vec_init,a_γ_hat_init=a_γ_hat_init, b_γ_hat_init=b_γ_hat_init,θ_hat_init=θ_hat_vec_init,c_ttprime_init = c_ttprime_vec_init,rtik_init=rtik_init, v_tikj_init = v_tikj_vec_init,rhok_hat_init=rhok_hat_vec_init, omegak_hat_init=omegak_hat_vec_init,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false);


# variational_inference_dynamicHDP_vs7_2
## x_input, K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,a_η,b_η,ηkj_prior,num_iter,num_local_iter
# input_str_list = @name x_input, K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,a_η,b_η,ηkj_prior,num_iter,num_local_iter;
# input_key_list = Symbol.(naming_vec(input_str_list));
# input_var_list = [x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,a_η,b_η,ηkj_prior,num_iter,num_local_iter];
# inputs_dict = OrderedDict()
# addToDict!(inputs_dict,input_key_list,input_var_list);
# out7_2 = variational_inference_dynamicHDP_vs7_2(inputs_dict;mk_hat_init=mk_hat_vec_init, λ0k_hat_init=λ0k_hat_vec_init,a0k_hat_init=a0k_hat_vec_init, b0k_hat_init=b0k_hat_vec_init,m_err_hat_init=m_err_hat_vec_init, λ0_err_hat_init=λ0_err_hat_vec_init,a0_err_hat_init=a0_err_hat_vec_init, b0_err_hat_init=b0_err_hat_vec_init,awt_hat_init=awt_hat_vec_init, bwt_hat_init=bwt_hat_vec_init,a_αt_hat_init=a_αt_hat_vec_init, b_αt_hat_init=b_αt_hat_vec_init,a_γ_hat_init=a_γ_hat_init, b_γ_hat_init=b_γ_hat_init,θ_hat_init=θ_hat_vec_init,c_ttprime_init = c_ttprime_vec_init,rtik_init=rtik_init, v_tikj_init = v_tikj_vec_init,rhok_hat_init=rhok_hat_vec_init, omegak_hat_init=omegak_hat_vec_init,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false);


function update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,v_tikj,c_ttprime)
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
                    v_true = v_tikj[t][i][k][j][1]
                    v_false = v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) + 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function update_v_tikj12(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj12(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj12(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops12(x,N_error)
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
function update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value12(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end

function variational_inference_dynamicHDP_vs12(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

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

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj12(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops12(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs12(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj12(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops12(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo12(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end


function variational_inference_dynamicHDP_vs12_2(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,a_η,b_η,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    e_log_ηkj = [log.(el) for el in ηkj_prior]
    e_log_minus_ηkj = [log.(1 .- el) for el in ηkj_prior]
    a_ηkj_hat = 1
    b_ηkj_hat = 1
    Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj12(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops12(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,e_log_ηkj,e_log_minus_ηkj);
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end
        a_ηkj_hat = update_a_ηkj(v_tikj,a_η)
        b_ηkj_hat = update_b_ηkj(v_tikj,b_η)
        e_log_ηkj = log_ηjk_expected_value(a_ηkj_hat,b_ηkj_hat)
        e_log_minus_ηkj = log1minusηjk_expected_value(a_ηkj_hat,b_ηkj_hat)

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo12(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        imp_elbo = calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        
        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_,b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,a_ηkj_hat_,b_ηkj_hat_ = elbo_, rtik,v_tikj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk,a_ηkj_hat,b_ηkj_hat

    output_str_list = @name elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,a_ηkj_hat_,b_ηkj_hat_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_,b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,a_ηkj_hat_,b_ηkj_hat_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

function update_rtik_vs14(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μ_tij_err,v_tikj,c_ttprime)
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
                    # v_true = v_tikj[t][i][k][j][1]
                    # v_false = v_tikj[t][i][k][j][2]

                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] )

                    # println("last")
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
function update_v_tikj14(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj14(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj14(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops14(x,N_error)
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
function update_a0_err_hat_usingXhat14(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat14(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value14(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo14(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
    data_lb_sum = 0.0
    N_signal,N_error = update_N(rtik,v_tikj);
    K = length(rtik[1][1][1])
    # Nj_error = update_errorNj(N_error);
    # Nkj_signal = update_signalNkj(N_signal);       
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
    # x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
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

function variational_inference_dynamicHDP_vs14(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

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

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value14(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs14(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj14(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops14(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat14(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat14(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value14(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs14(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value14(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N(rtik,v_tikj);
        Nj_error = update_errorNj14(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops14(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat14(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat14(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value14(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo14(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end
## x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter
# input_str_list = @name x_input, K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter;
# input_key_list = Symbol.(naming_vec(input_str_list));
# input_var_list = [x_input,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter];


# inputs_dict = OrderedDict()
# addToDict!(inputs_dict,input_key_list,input_var_list);
# out12 = variational_inference_dynamicHDP_vs12(inputs_dict;mk_hat_init=mk_hat_vec_init, λ0k_hat_init=λ0k_hat_vec_init,a0k_hat_init=a0k_hat_vec_init, b0k_hat_init=b0k_hat_vec_init,a0_err_hat_init=a0_err_hat_vec_init, b0_err_hat_init=b0_err_hat_vec_init,awt_hat_init=awt_hat_vec_init, bwt_hat_init=bwt_hat_vec_init,a_αt_hat_init=a_αt_hat_vec_init, b_αt_hat_init=b_αt_hat_vec_init,a_γ_hat_init=a_γ_hat_init, b_γ_hat_init=b_γ_hat_init,θ_hat_init=θ_hat_vec_init,c_ttprime_init = c_ttprime_vec_init,rtik_init=rtik_init, v_tikj_init = v_tikj_vec_init,rhok_hat_init=rhok_hat_vec_init, omegak_hat_init=omegak_hat_vec_init,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false);



# elbo_original_,rtik_original_,c_ttprime_original_,θ_hat_original_, mk_hat_vec_original_, λ0k_hat_vec_original_,a0k_hat_vec_original_,b0k_hat_vec_original_,rhok_hat_vec_original_, omegak_hat_vec_original_,a_αt_hat_vec_original_,b_αt_hat_vec_original_,awt_hat_vec_original_,bwt_hat_vec_original_,a_γ_hat_original_,b_γ_hat_original_,initDict_,debug_val = variational_inference_dynamicHDP(x_input, G,KMax,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter, num_local_iter;mk_hat_vec_init=mk_hat_init, λ0k_hat_vec_init=λ0k_hat_init,a0k_hat_vec_init=a0k_hat_init, b0k_hat_vec_init=b0k_hat_init,awt_hat_vec_init=awt_hat_init,  bwt_hat_vec_init=bwt_hat_init,a_αt_hat_vec_init=a_αt_hat_init, b_αt_hat_vec_init=b_αt_hat_init,a_γ_hat_init=a_γ_hat_init, b_γ_hat_init=b_γ_hat_init,θ_hat_vec_init=θ_hat_init,c_ttprime_vec_init = c_ttprime_init, rhok_hat_vec_init=rhok_hat_init, omegak_hat_vec_init=omegak_hat_init,uniform_theta_init = uniform_theta_init, rand_init = rand_init);


# true_z = z;
# num_posterior_samples=1000;
# gnu_elbo_plot(out7_2[:elbo_];fig_size=(1100,900))
# z_post_s_original = vi_make_z_post_s(out7_2[:rtik_], S=num_posterior_samples);
# vmeasure_vov_original = getVmeasure(true_z,z_post_s_original;beta=1.0);
# ari_vov_original,ri_vov_original,mirkinindx_vov_original,hubertindx_vov_original  = getRandIndices(true_z,z_post_s_original);
# nmi_vov_original = getNMI(true_z,z_post_s_original);
# avg_counts_mat_original = get_average_posterior_cluster_frequency2(z_post_s_original,T,true_z,KMax,KTrue,C_t);
# mean_τ_post_original,mean_μ_post_original = calc_normalgamma_μ_τ_post_mean(mk_hat_vec_original_,λ0k_hat_vec_original_,a0k_hat_vec_original_,b0k_hat_vec_original_;Ns= 1000);
# gnu_plot_All_Time_Points_PCA(x_to_use,z_post_s_original[end],G,unique(vcat(z_post_s_original[end]...));unique_clusters_labels = unique(vcat(z_post_s_original[end]...)))
# gnu_elbo_plot(elbo_original_[broadcast(!,ismissing.(elbo_original_))];fig_size=(1100,900))

################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs15(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = pip_kj[k][j] * 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips15(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [a0_err_hat_vec ./ b0_err_hat_vec]#calc_gamma_τ_post_mean(;Ns= 1000);
    # mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec;Ns= 1000);
    # mean_τ_err_post = calc_gamma_τ_post_mean([a0_err_hat_vec],[b0_err_hat_vec];Ns= 1000);
    mean_μ_err_post = [zeros(G)];


    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N15(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
                    Ntik_error[j] = rtik[t][i][k] * v_tikj[t][i][k][j][2] * (pip_kj[k][j])
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

function update_v_tikj15(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj15(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj15(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops15(x,N_error)
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
function update_a0_err_hat_usingXhat15(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat15(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value15(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo15(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N15(rtik,v_tikj,pip_kj);
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
function calc_Hpip(pip_kj)
    pip_entropy = 0.0
    K = length(pip_kj)
    for k in 1:K
        pip_entropy += entropy(pip_kj[k])
    end
    return -pip_entropy
end
function variational_inference_dynamicHDP_vs15(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value15(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs15(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N15(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj15(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops15(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat15(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat14(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value15(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj15(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips15(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs15(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value15(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs15(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N15(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj15(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops15(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat15(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat15(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value15(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj15(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips15(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo15(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_Nkj16(N_rpip)
    Nkj = sum(sum.(N_rpip))#sum(Ntk)[1:end-1]
    return Nkj
end
function update_N_rpip16(rtik,pip_kj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(pip_kj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Ntik_signal =  rtik[t][i][k] .* pip_kj[k]
                # for j in 1:G
                #     Ntik_signal[j] =
                # end
                Nti_signal[k] = Ntik_signal
            end
            Nt_signal[i] = Nti_signal
        end
        N_signal[t] = Nt_signal
    end
    return N_signal
end
function calc_DataElbo16(x,rpip,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
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
    weighted_ss_kjti = [[[[rpip[t][i][k][j] .*(x[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
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
function update_rtik_vs16(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = pip_kj[k][j] * 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function update_a0k_hat_usingXhat16(a0_vec,Nkj)
    K = length(Nkj) 
    a0k_hat_vec = [ a0_vec .+ 1/2 .* (Nkj[k] .+1 ) for k in 1:K]
    return a0k_hat_vec
end

function variational_inference_dynamicHDP_vs16(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


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
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        Nk_chain = make_chain(num_iter+1,Nk)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)

        
        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec

        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:Nk_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik_vs16(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end
        rpip = update_N_rpip16(rtik,pip_kj)
        Nkj = update_Nkj16(rpip)
        x_hat_k = update_x_hat_k(x,rpip)
        x_hat_sq_k = update_x_hat_sq_k(x,rpip)
        if record_chain
            chain_dict[:Nk_chain][iter + 1] = Nk
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
        a0k_hat_vec = update_a0k_hat_usingXhat16(a0_vec,Nkj)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec
        end

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        pip_kj = get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)

        data_elbo = calc_DataElbo16(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        pip_entropy = calc_Hpip(pip_kj);

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + pip_entropy
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
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end


################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs17(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips17(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [a0_err_hat_vec ./ b0_err_hat_vec]#calc_gamma_τ_post_mean(;Ns= 1000);
    # mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec;Ns= 1000);
    # mean_τ_err_post = calc_gamma_τ_post_mean([a0_err_hat_vec],[b0_err_hat_vec];Ns= 1000);
    mean_μ_err_post = [zeros(G)];


    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N17(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
                    Ntik_error[j] = rtik[t][i][k] * v_tikj[t][i][k][j][2] * (pip_kj[k][j])
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

function update_v_tikj17(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj17(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj17(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops17(x,N_error)
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
function update_a0_err_hat_usingXhat17(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat17(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value17(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo17(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N17(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs17(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value17(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs17(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N17(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj17(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops17(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat17(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat17(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value17(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj17(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips17(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs17(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value17(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs17(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N17(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj17(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops17(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat17(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat17(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value17(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj17(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips17(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo17(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_Nkj18(N_rpip)
    Nkj = sum(sum.(N_rpip))#sum(Ntk)[1:end-1]
    return Nkj
end
function update_N_rpip18(rtik,pip_kj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(pip_kj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Ntik_signal =  rtik[t][i][k] .* pip_kj[k]
                # for j in 1:G
                #     Ntik_signal[j] =
                # end
                Nti_signal[k] = Ntik_signal
            end
            Nt_signal[i] = Nti_signal
        end
        N_signal[t] = Nt_signal
    end
    return N_signal
end
function calc_DataElbo18(x,rpip,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
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
    weighted_ss_kjti = [[[[rpip[t][i][k][j] .*(x[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
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
function update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene =  0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function update_a0k_hat_usingXhat18(a0_vec,Nkj)
    K = length(Nkj) 
    a0k_hat_vec = [ a0_vec .+ 1/2 .* (Nkj[k] .+1 ) for k in 1:K]
    return a0k_hat_vec
end
function calc_Hpip(pip_kj)
    pip_entropy = 0.0
    K = length(pip_kj)
    for k in 1:K
        pip_entropy += entropy(pip_kj[k])
    end
    return -pip_entropy
end
function variational_inference_dynamicHDP_vs18(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


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
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        Nk_chain = make_chain(num_iter+1,Nk)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)

        
        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec

        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:Nk_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            Glog = G*log(2π)
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end
        rpip = update_N_rpip18(rtik,pip_kj)
        Nkj = update_Nkj18(rpip)
        x_hat_k = update_x_hat_k(x,rpip)
        x_hat_sq_k = update_x_hat_sq_k(x,rpip)
        if record_chain
            chain_dict[:Nk_chain][iter + 1] = Nk
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
        a0k_hat_vec = update_a0k_hat_usingXhat18(a0_vec,Nkj)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec
        end

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        pip_kj = get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)

        data_elbo = calc_DataElbo18(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        pip_entropy = calc_Hpip(pip_kj);

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + pip_entropy
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
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end



################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs19(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips19(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [a0_err_hat_vec ./ b0_err_hat_vec]#calc_gamma_τ_post_mean(;Ns= 1000);
    # mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec;Ns= 1000);
    # mean_τ_err_post = calc_gamma_τ_post_mean([a0_err_hat_vec],[b0_err_hat_vec];Ns= 1000);
    mean_μ_err_post = [zeros(G)];


    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N19(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
                    Ntik_error[j] = rtik[t][i][k] * v_tikj[t][i][k][j][2] * (1 - pip_kj[k][j])
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

function update_v_tikj19(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj19(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj19(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops19(x,N_error)
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
function update_a0_err_hat_usingXhat19(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat19(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value19(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo19(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N19(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs19(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value19(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs19(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N19(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj19(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops19(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat19(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat19(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value19(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj19(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips19(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs19(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value19(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs19(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N19(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj20(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops19(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat19(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat19(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value19(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj19(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips19(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo19(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end





################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs20(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips20(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [a0_err_hat_vec ./ b0_err_hat_vec]#calc_gamma_τ_post_mean(;Ns= 1000);
    # mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec;Ns= 1000);
    # mean_τ_err_post = calc_gamma_τ_post_mean([a0_err_hat_vec],[b0_err_hat_vec];Ns= 1000);
    mean_μ_err_post = [zeros(G)];


    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N20(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
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

function update_v_tikj20(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj20(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj20(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops20(x,N_error)
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
function update_a0_err_hat_usingXhat20(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat20(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value20(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo20(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N20(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs20(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value20(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs20(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N20(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj20(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops20(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat20(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat20(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value20(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj20(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj =  get_gene_pips20(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
        
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs20(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value20(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs20(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N20(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj20(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops20(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat20(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat20(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value20(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj20(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips20(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo20(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end


################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs21(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips21(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj;null_precision=10)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [null_precision .* ones(G) for k in 1:K]
    mean_μ_err_post = [zeros(G) for k in 1:K]

    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N21(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
                    Ntik_error[j] = rtik[t][i][k] * v_tikj[t][i][k][j][2] * (pip_kj[k][j])
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

function update_v_tikj21(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj21(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj21(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops21(x,N_error)
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
function update_a0_err_hat_usingXhat21(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat21(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value21(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo21(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N21(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs21(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj =  [ones(G) ./ G  for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value21(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs21(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N21(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj21(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops21(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat21(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat21(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value21(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj21(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips21(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj);
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs21(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value21(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs21(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N21(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj21(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops21(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat21(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat21(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value21(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj21(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips21(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj);
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo21(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end



################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs22(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips22(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [a0_err_hat_vec ./ b0_err_hat_vec]#calc_gamma_τ_post_mean(;Ns= 1000);
    # mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec;Ns= 1000);
    # mean_τ_err_post = calc_gamma_τ_post_mean([a0_err_hat_vec],[b0_err_hat_vec];Ns= 1000);
    mean_μ_err_post = [zeros(G)];


    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,[mean_τ_err_post[1] for k in 1:K],[mean_μ_err_post[1] for k in 1:K]);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N22(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
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

function update_v_tikj22(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj22(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj22(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops22(x,N_error)
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
function update_a0_err_hat_usingXhat22(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat22(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value22(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo22(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N22(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs22(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value22(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs22(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N22(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj22(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops22(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat22(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat22(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value22(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj22(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj =  get_gene_pips22(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
        
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs22(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value22(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs22(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N22(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj22(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops22(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat22(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat22(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value22(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj22(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips22(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo22(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end



################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs23(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips23(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj;null_precision=10)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [null_precision .* ones(G) for k in 1:K]
    mean_μ_err_post = [zeros(G) for k in 1:K]

    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ) for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] .+ log(v_tikj[t][i][k][j][1] ))  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N23(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
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

function update_v_tikj23(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj23(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj23(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops23(x,N_error)
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
function update_a0_err_hat_usingXhat23(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat23(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value23(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo23(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N23(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs23(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value23(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs23(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N23(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj23(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops23(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat23(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat23(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value23(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj23(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj =  get_gene_pips23(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
        
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs23(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value23(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs23(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N23(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj23(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops23(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat23(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat23(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value23(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj23(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips23(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo23(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end




################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs24(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    v_true = 1#v_tikj[t][i][k][j][1]
                    v_false = 0#v_tikj[t][i][k][j][2]
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene = 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function get_gene_pips24(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj;null_precision=10)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    mean_τ_err_post = [null_precision .* ones(G) for k in 1:K]
    mean_μ_err_post = [zeros(G) for k in 1:K]

    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] )  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end

function update_N24(rtik,v_tikj,pip_kj)
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
                    Ntik_signal[j] = rtik[t][i][k] * v_tikj[t][i][k][j][1] * pip_kj[k][j]
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

function update_v_tikj24(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_v_tikj24(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,e_log_ηkj,e_log_minus_ηkj)
    T = length(e_τ_μ_tikj)
    K = length(e_log_τkj)
    C_t = [length(el) for el in e_τ_μ_tikj]
    G = length(e_τ_μ_tikj[1][1][1])
    v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    v_tikj_tilde = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T)
    logpi= Glog/G
    for t in 1:T
        
        cells_ = C_t[t]
        v_t = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        v_t_tilde = Vector{Vector{Vector{Vector{Float64}}}}(undef,cells_)
        for i in 1:cells_
            log_v_ti = Vector{Vector{Vector{Float64}}}(undef,K)
            log_v_ti_tilde = Vector{Vector{Vector{Float64}}}(undef,K)
            for k in 1:K
                log_v_tik = Vector{Vector{Float64}}(undef,G)
                log_v_tik_tilde = Vector{Vector{Float64}}(undef,G)
                for j in 1:G
                    # println(k)
                    log_v_tikj = Vector{Float64}(undef,2)
                    log_v_tikj_tilde = Vector{Float64}(undef,2) 
                    log_v_tikj_tilde[1] = rtik[t][i][k] * (0.5 * e_log_τkj[k][j] - 0.5 * logpi - 0.5 * e_τ_μ_tikj[t][i][k][j]) + e_log_ηkj[k][j]#log(ηkj_prior[k][j]) 
                    log_v_tikj_tilde[2] = rtik[t][i][k] * (0.5 * e_log_τj_err[j] - 0.5 * logpi - 0.5 * e_τ_μj_err[t][i][j]) + e_log_minus_ηkj[k][j]#log(1 - ηkj_prior[k][j])
                    log_v_tikj = norm_weights(log_v_tikj_tilde)
                    log_v_tik_tilde[j] = log_v_tikj_tilde
                    log_v_tik[j] = log_v_tikj
                    # println(" not broke")
                end
                log_v_ti[k] = log_v_tik
                log_v_ti_tilde[k] = log_v_tik_tilde
            end
            v_t[i] = log_v_ti
            v_t_tilde[i] = log_v_ti_tilde
        end
        v_tikj[t] = v_t
        v_tikj_tilde[t] = v_t_tilde
    end
    return v_tikj,v_tikj_tilde
end
function update_errorNj24(N_error) # FASTER and LESS MEMORY when compared to @benchmark  update_errorNj_forloops(N_error)
    # T = length(N_error)
    # K = length(N_error[1][1])
    # C_t = [length(el) for el in N_error]
    # G = length(N_error[1][1][1])
    perCell_linerize_N_error = reduce(vcat,N_error)
    perCell_perState_linerize_N_error = reduce(vcat,perCell_linerize_N_error)
    Nj_error = sum(perCell_perState_linerize_N_error)
    return Nj_error
end
function update_x_hat_sq_error_vs_forloops24(x,N_error)
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
function update_a0_err_hat_usingXhat24(a0_err_vec,Nj_error)
    a0_err_hat_vec = a0_err_vec .+ 1/2 .* (Nj_error ) #(Nj_error .+1)
    return a0_err_hat_vec
end
function update_b0_err_hat_usingXhat24(b0_err_vec,x_hat_sq_err)
    b0_err_hat_vec = b0_err_vec .+  1/2 .* (x_hat_sq_err)
    return  b0_err_hat_vec
end
function τ_μ_error_expected_value24(x,a0_err_vec, b0_err_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(a0_err_vec)
    e_τ_μ_kj_true3 = Vector{Vector{Vector}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector}(undef,cells)
        e_τ_μ_13 =  Vector{Float64}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i]) .^2
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function calc_DataElbo24(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec,pip_kj)
    data_lb_sum = 0.0
    N_signal,N_error = update_N24(rtik,v_tikj,pip_kj);
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

function variational_inference_dynamicHDP_vs24(x, G,K,λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,num_iter,a0_err,b0_err,ηkj_prior, num_local_iter;mk_hat_vec_init=nothing, λ0k_hat_vec_init=nothing,a0k_hat_vec_init=nothing, b0k_hat_vec_init=nothing,awt_hat_vec_init=nothing, bwt_hat_vec_init=nothing,a_αt_hat_vec_init=nothing, b_αt_hat_vec_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_vec_init=nothing,c_ttprime_vec_init = nothing, rhok_hat_vec_init=nothing, omegak_hat_vec_init=nothing, v_tikj_vec_init = nothing,rtik_init = nothing,a0_err_hat_vec_init=nothing, b0_err_hat_vec_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),debugme = false)
    T = length(x)
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    _, _, a0_err_vec, b0_err_vec = init_params_genes(G,1,0,a0_err,b0_err);

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
    if isnothing(v_tikj_vec_init) && rand_init
        v_tikj_vec_init = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    elseif isnothing(v_tikj_vec_init) && !rand_init
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
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

    
    # if isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = rand(Uniform(0,1),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    # elseif isnothing(λ0_err_hat_vec_init) && rand_init
    #     λ0_err_hat_vec_init = λ0_err_vec
    # end

    # if isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    # elseif isnothing(m_err_hat_vec_init) && rand_init
    #     m_err_hat_vec_init =μ0_err_vec
    # end

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

     
    v_tikj = v_tikj_vec_init
    awt_hat_vec = awt_hat_vec_init 
    bwt_hat_vec = bwt_hat_vec_init
    a_αt_hat_vec = a_αt_hat_vec_init 
    b_αt_hat_vec = b_αt_hat_vec_init
    θ_hat_vec = θ_hat_vec_init
    c_ttprime_vec = c_ttprime_vec_init
    # λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init = λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec;
    # λ0_err_hat_vec, m_err_hat_vec, a0_err_hat_vec, b0_err_hat_vec =  λ0_err_hat_vec_init, m_err_hat_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init;
    a0_err_hat_vec, b0_err_hat_vec =  a0_err_hat_vec_init, b0_err_hat_vec_init;

    


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

    arg_str_list_initparams = @name mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init,  a0_err_hat_vec_init, b0_err_hat_vec_init ;
    key_list_initparams = naming_vec(arg_str_list_initparams);
    var_list_initparams = [mk_hat_vec_init,λ0k_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init,omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init, bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,θ_hat_vec_init,c_ttprime_vec_init, a0_err_hat_vec_init, b0_err_hat_vec_init ];

    initDict = Dict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    debug_val = nothing
    #init debug dict
    if debugme
        num_local_iter = 1
        debug_val = OrderedDict()
        debug_val["λ0k_hat_vec"] = []
        debug_val["mk_hat_vec"]= []
        debug_val["a0k_hat_vec"]= []
        debug_val["b0k_hat_vec"]= []
        # debug_val["λ0_err_hat_vec"] = []
        # debug_val["m_err_hat_vec"]= []
        debug_val["a0_err_hat_vec"]= []
        debug_val["b0_err_hat_vec"]= []
        debug_val["rhok_hat_vec"]= []
        debug_val["omegak_hat_vec" ]= []
        debug_val["θ_hat"]= []
        debug_val["rtik"]= []
        debug_val["v_tikj"]= []
        debug_val["c_ttprime_vec"]= []
        debug_val["e_log_π"]= []
        debug_val["e_log_τ"]= []
        debug_val["e_log_τkj"]= []
        debug_val["e_τ_μ_tikj"]= []
        debug_val["e_τ_μ"]= []
        debug_val["e_log_τj_err"]= []
        debug_val["e_τ_0j_err"]= []
        debug_val["e_τ_0_err"]= [] 
        debug_val["Ntk"]= []
        debug_val["Nk"]= []
        debug_val["N_signal"]= []
        debug_val["N_error"]= []
        debug_val["Nj_error"]= []
        debug_val["Nkj_signal"]= [] 
        debug_val["x_hat_err"]=[]
        debug_val["x_hatk_signal"]=[]
        debug_val["x_hat_sq_err"]=[]
        debug_val["x_hatk_sq_signal"]=[]
        debug_val["e_γ" ]=[]
        debug_val["a_αt_hat_vec" ]=[]
        debug_val["b_αt_hat_vec" ]=[]
        debug_val["awt_hat_vec" ]=[]
        debug_val["bwt_hat_vec" ]=[]
        debug_val["a_γ_hat" ]=[]
        debug_val["b_γ_hat" ]=[]
        debug_val["Tαk"]= []
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
        # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec )
        # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
        push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
        push!(debug_val["b0_err_hat_vec"], b0_err_hat_vec)
        push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
        push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
        push!(debug_val["θ_hat"],θ_hat_vec)
        push!(debug_val["rtik"],[])
        push!(debug_val["v_tikj"],[])
        push!(debug_val["c_ttprime_vec"],[])
        push!(debug_val["e_log_π"],[])
        push!(debug_val["e_log_τ"],[])
        push!(debug_val["e_τ_μ_tikj"],[])
        push!(debug_val["e_τ_μ"],[])
        push!(debug_val["e_log_τkj"],[])
        push!(debug_val["e_log_τj_err"],[])
        push!(debug_val["e_τ_0j_err"],[])
        push!(debug_val["e_τ_0_err"],[]) 
        push!(debug_val["Ntk"],[])
        push!(debug_val["Nk"],[])
        push!(debug_val["N_signal"],[])
        push!(debug_val["N_error"],[])
        push!(debug_val["Nj_error"],[])
        push!(debug_val["Nkj_signal"],[])
        push!(debug_val["x_hat_err"],[])
        push!(debug_val["x_hatk_signal"],[])
        push!(debug_val["x_hat_sq_err"],[])
        push!(debug_val["x_hatk_sq_signal"],[])
        push!(debug_val["e_γ" ],[])
        push!(debug_val["a_αt_hat_vec" ],[])
        push!(debug_val["b_αt_hat_vec" ],[])
        push!(debug_val["awt_hat_vec" ],[])
        push!(debug_val["bwt_hat_vec" ],[])
        push!(debug_val["a_γ_hat" ],[])
        push!(debug_val["b_γ_hat" ],[])
        push!(debug_val["Tαk"],[])
        push!(debug_val["data_elbo"],[])
        push!(debug_val["assgn_entropy"],[])
        push!(debug_val["HDP_surragate_elbo"],[])
    end

    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]

    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G for k in 1:K]
    for iter in 1:num_iter
        # println("Global Iteration: $iter")
        for loc_iter in 1:num_local_iter
            # println("Local Iteration: $loc_iter")
            
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value24(x, a0_err_hat_vec, b0_err_hat_vec);

            # v_tikj,_  = update_v_tikj7(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,ηkj_prior);
            # if debugme
            #     push!(debug_val["v_tikj"],v_tikj)
            # end

            # rtik = update_rtik_vs(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_μj_err,η_tikj,c_ttprime_vec); #update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
            rtik = update_rtik_vs24(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)




            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if debugme
                push!(debug_val["θ_hat"],θ_hat_vec)
                push!(debug_val["rtik"],rtik)
                push!(debug_val["e_log_π"],e_log_π)
                push!(debug_val["e_log_τ"],e_log_τ)
                push!(debug_val["e_τ_μ_tikj"],e_τ_μ_tikj)
                push!(debug_val["e_τ_μ"],e_τ_μ)
                push!(debug_val["Ntk"],Ntk)
                # push!(debug_val["v_tikj"],v_tikj)
                push!(debug_val["c_ttprime_vec"],c_ttprime_vec)
                push!(debug_val["e_log_τj_err"],e_log_τj_err)
                push!(debug_val["e_τ_0j_err"],e_τ_0j_err)
                push!(debug_val["e_τ_0_err"],e_τ_0_err) 
            end

        end


        
        
        
        # sk = 1 ./ Nk .* sk

        
        
        

        N_signal,N_error = update_N24(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj24(N_error)
        Nkj_signal = update_signalNkj(N_signal)
        if debugme
            # push!(debug_val["Nk"],Nk)
            push!(debug_val["N_signal"],N_signal)
            push!(debug_val["N_error"],N_error)
            push!(debug_val["Nj_error"],Nj_error)
            push!(debug_val["Nkj_signal"],Nkj_signal)
        end

        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops24(x,N_error)
        x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if debugme
            push!(debug_val["x_hat_err"],x_hat_err)
            push!(debug_val["x_hatk_signal"],x_hatk_signal)
            push!(debug_val["x_hat_sq_err"],x_hat_sq_err)
            push!(debug_val["x_hatk_sq_signal"],x_hatk_sq_signal)
        end




        
        
        a0_err_hat_vec = update_a0_err_hat_usingXhat24(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat24(b0_err_vec,x_hat_sq_err)


        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)


        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value24(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj24(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj =  get_gene_pips24(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj)
        if debugme
            push!(debug_val["v_tikj"],v_tikj)
        end
        
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
        
        Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        if debugme
            push!(debug_val["λ0k_hat_vec"],λ0k_hat_vec)
            push!(debug_val["mk_hat_vec"],mk_hat_vec)
            push!(debug_val["a0k_hat_vec"],a0k_hat_vec)
            push!(debug_val["b0k_hat_vec"],b0k_hat_vec)
            # push!(debug_val["m_err_hat_vec"],m_err_hat_vec)
            # push!(debug_val["λ0_err_hat_vec"],λ0_err_hat_vec)
            push!(debug_val["a0_err_hat_vec"],a0_err_hat_vec)
            push!(debug_val["b0_err_hat_vec"],b0_err_hat_vec)
            push!(debug_val["rhok_hat_vec"],rhok_hat_vec)
            push!(debug_val["omegak_hat_vec" ],omegak_hat_vec)
            push!(debug_val["e_γ" ],e_γ)
            push!(debug_val["a_αt_hat_vec" ],a_αt_hat_vec)
            push!(debug_val["b_αt_hat_vec" ],b_αt_hat_vec)
            push!(debug_val["awt_hat_vec" ],awt_hat_vec)
            push!(debug_val["bwt_hat_vec" ],bwt_hat_vec)
            push!(debug_val["a_γ_hat" ],a_γ_hat)
            push!(debug_val["b_γ_hat" ],b_γ_hat)
            # Tαk,e_γ,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


            push!(debug_val["Tαk"],Tαk)
            # push!(debug_val["data_elbo"],data_elbo)
            # push!(debug_val["assgn_entropy"],assgn_entropy)
            # push!(debug_val["HDP_surragate_elbo"],HDP_surragate_elbo)
        end

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
    
    return elbo_, rtik,c_ttprime_vec,v_tikj, pip_kj,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat, a0_err_hat_vec, b0_err_hat_vec,initDict,debug_val
end
function variational_inference_dynamicHDP_vs24(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing, a0_err_hat_init=nothing, b0_err_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing, v_tikj_init = nothing,rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x,K,a0,b0,μ0,λ0,a0_err,b0_err,μ0_err,λ0_err,a_γ,b_γ,a_α,b_α,adot_w,bdot_w,ηkj_prior,num_iter,num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);

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

    v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)


    # m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
    # λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
    a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
    b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)

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
    # m_err_hat_vec= m_err_hat_init 
    # λ0_err_hat_vec = λ0_err_hat_init
    a0_err_hat_vec = a0_err_hat_init
    b0_err_hat_vec = b0_err_hat_init
    v_tikj = v_tikj_init

    
    rtik = rtik_init
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,a0_err_hat_init,b0_err_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init,v_tikj_init, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_log_τj_err = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_0j_err = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        e_τ_0_err = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        x_hat_err = Vector{Float64}()
        x_hat_sq_err = Vector{Float64}()

        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        N_signal = Vector{Vector{Vector{Vector{Float64}}}}()
        N_error = Vector{Vector{Vector{Vector{Float64}}}}()
        Nj_error = Vector{Float64}()
        Nkj_signal = Vector{Vector{Float64}}()
        # a0_err_hat_vec = Vector{Float64}()
        # λ0_err_hat_vec = Vector{Float64}()
        # m_err_hat_vec = Vector{Float64}()
        # b0_err_hat_vec = Vector{Float64}()
        v_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,T);
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        # λ0_err_chain = make_chain(num_iter+1,λ0_err_hat_vec)
        # m_err_chain = make_chain(num_iter+1,m_err_hat_vec)
        a0_err_chain = make_chain(num_iter+1,a0_err_hat_vec)
        b0_err_chain = make_chain(num_iter+1,b0_err_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        e_log_τj_err_chain=make_chain(num_iter+1,e_log_τj_err)
        e_τ_0j_err_chain=make_chain(num_iter+1,e_τ_0j_err)
        e_τ_0_err_chain=make_chain(num_iter+1,e_τ_0_err)
        v_tikj_chain = make_chain(num_iter+1,v_tikj)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        N_signal_chain = make_chain(num_iter+1,N_signal)
        N_error_chain = make_chain(num_iter+1,N_error)
        Nj_error_chain = make_chain(num_iter+1,Nj_error)
        Nkj_signal_chain = make_chain(num_iter+1,Nkj_signal)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        x_hat_err_chain = make_chain(num_iter+1,x_hat_err)
        x_hat_sq_err_chain = make_chain(num_iter+1,x_hat_sq_err)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)


        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,e_log_τj_err_chain,e_τ_0j_err_chain,e_τ_0_err_chain,v_tikj_chain,Ntk_chain,c_ttprime_chain,a0_err_chain,b0_err_chain,N_signal_chain,N_error_chain,Nj_error_chain,Nkj_signal_chain,x_hat_err_chain,x_hat_sq_err_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:λ0_err_chain][1] = λ0_err_hat_vec
        # chain_dict[:m_err_chain][1] = m_err_hat_vec
        chain_dict[:a0_err_chain][1] = a0_err_hat_vec
        chain_dict[:b0_err_chain][1] = b0_err_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec
        chain_dict[:v_tikj_chain][1] = nothing
        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:e_log_τj_err_chain][1] = nothing
        chain_dict[:e_τ_0j_err_chain][1] = nothing
        chain_dict[:e_τ_0_err_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:N_signal_chain][1] = nothing
        chain_dict[:N_error_chain][1] = nothing
        chain_dict[:Nj_error_chain][1] = nothing
        chain_dict[:Nkj_signal_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:x_hat_err_chain][1] = nothing
        chain_dict[:x_hat_sq_err_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
    end
    # θ_hat = init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)
    # rand_permutation= append!(shuffle(1:K),K+1)
    # θ_hat = [θ_hat[t][rand_permutation] for t in 1:T]
    # rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    elbo_ = Vector{Union{Missing,Float64}}(undef,num_iter)
    # elbo_ = convert(Array{Union{Missing,Float64}},elbo_)
    # elbo_ .= missing
    iter = 1
    converged_bool = false
    is_converged = false
    Glog = G*log(2π)
    pip_kj = [ones(G) ./ G  for k in 1:K]
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            e_log_π = log_π_expected_value(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);

            e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value24(x, a0_err_hat_vec, b0_err_hat_vec);
            rtik = update_rtik_vs24(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec);

            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:e_log_τj_err_chain][iter + 1][loc_iter] = e_log_τj_err
                chain_dict[:e_τ_0j_err_chain][iter + 1][loc_iter] = e_τ_0j_err
                chain_dict[:e_τ_0_err_chain][iter + 1][loc_iter] = e_τ_0_err
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end

        # Nk = update_Nk(rtik)
        N_signal,N_error = update_N24(rtik,v_tikj, pip_kj);
        Nj_error = update_errorNj24(N_error)
        Nkj_signal = update_signalNkj(N_signal)        
        # x_hat_k = update_x_hat_k(x,rtik)
        # x_hat_sq_k = update_x_hat_sq_k(x,rtik)
        x_hat_err = nothing#update_x_hat_error_vs_forloops(x,N_error)
        x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
        x_hat_sq_err = update_x_hat_sq_error_vs_forloops24(x,N_error)
        x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
        if record_chain
            chain_dict[:N_signal_chain][iter + 1] = N_signal
            chain_dict[:N_error_chain][iter + 1] = N_error
            chain_dict[:Nj_error_chain][iter + 1] = Nj_error
            chain_dict[:Nkj_signal_chain][iter + 1] = Nkj_signal
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
            chain_dict[:x_hat_err_chain][iter + 1] = x_hat_err
            chain_dict[:x_hat_sq_err_chain][iter + 1] = x_hat_sq_err
        end
        



        a0_err_hat_vec = update_a0_err_hat_usingXhat24(a0_err_vec,Nj_error)
        # λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
        # m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
        b0_err_hat_vec = update_b0_err_hat_usingXhat24(b0_err_vec,x_hat_sq_err)

        λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
        a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
        mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
        b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec

            # chain_dict[:λ0_err_chain][iter + 1] = λ0_err_hat_vec
            # chain_dict[:m_err_chain][iter + 1] = m_err_hat_vec
            chain_dict[:a0_err_chain][iter + 1] = a0_err_hat_vec
            chain_dict[:b0_err_chain][iter + 1] = b0_err_hat_vec
        end







        n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
        n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K

        n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
        n_e_τ_0j_err,_  = τ_μ_error_expected_value24(x, a0_err_hat_vec, b0_err_hat_vec);
        v_tikj,_  = update_v_tikj24(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
        pip_kj = get_gene_pips24(x,mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik,v_tikj)
        if record_chain
            chain_dict[:v_tikj_chain][iter + 1] = v_tikj
        end

        

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat


        data_elbo = calc_DataElbo24(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec, pip_kj)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
        v_entropy = calc_Hv(v_tikj)
        pip_entropy = calc_Hpip(pip_kj)

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + imp_elbo + v_entropy + pip_entropy
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
        # if iter == num_iter
        #     converged_bool = true
        #     is_converged = true
        # end
    end
    
    nonemptychain_indx = broadcast(!,ismissing.(elbo_)) 
    elbo_ = elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) + 1

    if record_chain
        chain_dict = truncate_chain(chain_dict,truncation_value)
    end
    
    elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_, a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,v_tikj, pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec, a0_err_hat_vec, b0_err_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,v_tikj_, pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,a0_err_hat_vec_, b0_err_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end


################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs25(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene =  0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function update_Nkj25(N_rpip)
    Nkj = sum(sum.(N_rpip))#sum(Ntk)[1:end-1]
    return Nkj
end
function update_x_hat_sq_k25(x,rtik)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])

    x_hat_sq_tik = [[[ rtik[t][i][k] .* (x[t][i]) .^2 for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_sq_tk = [[ sum(x_hat_sq_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_sq_k = [sum(x_hat_sq_tk[k]) for k in 1:K]
    return x_hat_sq_k
end
function update_x_hat_k25(x,rtik)
    T = length(rtik)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    x_hat_tik = [[[ rtik[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_k
end
function τ_μ_expected_value25(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(a0k_hat_vec)
    G = length(a0k_hat_vec[1])
    e_τ_μ_kj_true3 = Vector{Vector{Vector{Vector}}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector{Vector}}(undef,cells)
        e_τ_μ_13 =  Vector{Vector{Float64}}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = Vector{Vector{Float64}}(undef,K)
            e_τ_μ_23 =  Vector{Float64}(undef,K)
            for k in 1:K
                e_τ_μ_kjti_cell3  =  a0k_hat_vec[k] ./  b0k_hat_vec[k] .*  (x[t][i] .- mk_hat_vec[k]) .^2 .+ 1 ./λ0k_hat_vec[k]
                e_τ_μ_3_cell3  = sum(e_τ_μ_kjti_cell3)
                
                e_τ_μ_kjti3[k] = e_τ_μ_kjti_cell3
                e_τ_μ_23[k] = e_τ_μ_3_cell3
            end
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function log_τ_kj_expected_value25(a0k_hat_vec, b0k_hat_vec)
    K = length(a0k_hat_vec)
    e_log_τ_kj_vec = []
    for k in 1:K
        e_log_τ_kj = digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])
        # e_log_τ_k = sum(e_log_τ_kj)
        push!(e_log_τ_kj_vec,e_log_τ_kj)
    end
    return e_log_τ_kj_vec
end
function log_τ_k_expected_value25(a0k_hat_vec, b0k_hat_vec)
    K = length(a0k_hat_vec)
    e_log_τ_kj_vec = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec)
    e_log_τ_k_vec = []
    for k in 1:K
        # e_log_τ_kj = digamma.(a0k_hat_vec_test[k]) .- log.(b0k_hat_vec[k])
        e_log_τ_k = sum(e_log_τ_kj_vec[k])
        push!(e_log_τ_k_vec,e_log_τ_k)
    end
    return e_log_τ_k_vec
end
function log_π_expected_value25(θ_hat)
    T = length(θ_hat)
    K = Int64(length(θ_hat[1])-1)
    e_log_π = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        # digamma_sum = digamma(sum(θ_hat[t][1:K]))
        # e_log_π[t] = digamma.(θ_hat[t][1:K]) .- digamma_sum
        digamma_sum = digamma(sum(θ_hat[t]))
        e_log_π[t] = digamma.(θ_hat[t]) .- digamma_sum
    end
    
    # param_val_dict["log_π_expected_value"] = e_log_π
    return e_log_π
end
function update_N_rpip25(rtik,pip_kj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(pip_kj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Ntik_signal =  rtik[t][i][k] .* pip_kj[k]
                # for j in 1:G
                #     Ntik_signal[j] =
                # end
                Nti_signal[k] = Ntik_signal
            end
            Nt_signal[i] = Nti_signal
        end
        N_signal[t] = Nt_signal
    end
    return N_signal
end
function calc_DataElbo25(x,rpip,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
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
    weighted_ss_kjti = [[[[rpip[t][i][k][j] .*(x[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
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
function update_a0k_hat_usingXhat25(a0_vec,Nkj)
    K = length(Nkj) 
    a0k_hat_vec = [ a0_vec .+ 1/2 .* (Nkj[k] .+1 ) for k in 1:K]
    return a0k_hat_vec
end
function get_gene_PIP25(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
    G = length(x[1][1])
    T = length(x)
    C_t = length.(x)
    K = length(rtik[1][1])
    # z_post_s = vi_make_z_post_s(rtik, S=S);

    mean_τ_post = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K ]
    mean_μ_post = mk_hat_vec#calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,;Ns= 1000);
    # mean_τ_err_post = [a0_err_hat_vec ./ b0_err_hat_vec]#calc_gamma_τ_post_mean(;Ns= 1000);
    # # mean_τ_post,mean_μ_post = calc_normalgamma_μ_τ_post_mean(mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec;Ns= 1000);
    # # mean_τ_err_post = calc_gamma_τ_post_mean([a0_err_hat_vec],[b0_err_hat_vec];Ns= 1000);
    # mean_μ_err_post = [zeros(G)];
    mean_τ_err_post = [null_precision .* ones(G) for k in 1:K]
    mean_μ_err_post = [zeros(G) for k in 1:K]
    # println(K)
    # # println(length(b0k_hat_vec_))
    # # println(length(a0k_hat_vec_[1]))
    # println(G)
    # println(T)
    # println("*****************")
    try
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
        calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);
    catch e
        println(mean_τ_post)
        println("##########")
        println(mean_μ_post)
        println("##########")
        println(mean_τ_err_post)
        println("##########")
        println(mean_μ_err_post)
        println("##########")
        println(a0k_hat_vec)
        println("##########")
        println(b0k_hat_vec)
        println("##########")
        println(mk_hat_vec)
        println("##########")
        println(a0_err_hat_vec)
        println("##########")
        println(b0_err_hat_vec)
    end

    cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post);
    null_cell_ll_scores = calc_cell_normal_μ_τ_ll_scores(x,mean_τ_err_post,mean_μ_err_post);


    expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] )  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # ration2_ = [[[[ration_[t][i][j][k] .+ log( rtik[t][i][k]) for k in 1:K] for j in 1:G]  for i in 1:C_t[t]] for t in 1:T];
    # ration2_weight = [[[norm_weights(ration2_[t][i][j]) .* rtik[t][i] for j in 1:G]  for i in 1:C_t[t]] for t in 1:T];
    
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    
    N_k = sum(sum.(rtik))
    
    gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k

    # sum(sum.(ration2_weight))
    # [el ./ N_k for el in sum(sum.(ration2_weight))]
    # normToProb.([el ./ N_k for el in sum(sum.(ration2_weight))])

    
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]
    nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
    nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
    function fix_nan_or_allzero!(v)
        K = length(v)
        G = length(v[1])
        for k in 1:K
            if any(isnan.(v[k])) || any(iszero.(v[k])) 
                if all(isnan.(v[k])) || all(iszero.(v[k]))
                    v[k] .= ones(Float64,G)
                else
                    v[k][isnan.(v[k])] .= 0.0
                end
            end
        end
        return v
    end
    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])
    pip_kj = normToProb.( fix_nan_or_allzero!(deepcopy(gene_significance_weights)))
    return pip_kj
end
function variational_inference_dynamicHDP_vs25(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,pip_kj_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, null_precision, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


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
    pip_kj_init = init_pip_kj_vec!(pip_kj_init,G,K;rand_init = false)

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
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        Nk_chain = make_chain(num_iter+1,Nk)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)

        
        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec

        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:Nk_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
            e_log_π = log_π_expected_value25(θ_hat_vec) # T by K
            # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            e_log_τkj = log_τ_kj_expected_value25(a0k_hat_vec, b0k_hat_vec);
            e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value25(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik_vs25(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end
        rpip = update_N_rpip25(rtik,pip_kj)
        Nkj = update_Nkj25(rpip)
        x_hat_k = update_x_hat_k25(x,rpip)
        x_hat_sq_k = update_x_hat_sq_k25(x,rpip)
        if record_chain
            chain_dict[:Nk_chain][iter + 1] = Nk
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
        a0k_hat_vec = update_a0k_hat_usingXhat25(a0_vec,Nkj)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec
        end

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        pip_kj = get_gene_PIP25(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=null_precision)

        data_elbo = calc_DataElbo25(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        pip_entropy = calc_Hpip(pip_kj);

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + pip_entropy
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
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

################################################################################
################################################################################
################################################################################
################################################################################################################################################################
################################################################################
################################################################################
function update_rtik_vs25_fast(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime)
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
                    # println("first")
                    # println("v_true: $v_true")
                    # println("v_false: $v_false")
                    # println("e_log_τkj[k][j]: $(e_log_τkj[k][j])")
                    # println("logpi: $(logpi)")
                    # println("e_τ_μ_tikj[t][i][k][j]: $(e_τ_μ_tikj[t][i][k][j])")
                    # println("e_log_τj_err[j]: $(e_log_τj_err[j])")
                    # println("e_τ_μ_tij_err[t][i][j]: $(e_τ_μ_tij_err[t][i][j])")

                    # log_like_gene = log(v_true)  + 0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j])  # v_false * 0.5 *(e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # log_like_gene = 0.5 * v_true * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] - e_log_τj_err[j] +  logpi + e_τ_μ_tij_err[t][i][j]) + 0.5 * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    log_like_gene =  0.5 * (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])
                    
                    # log_like_gene =  0.5 * w_kj[k][j]* (e_log_τkj[k][j] -  logpi - e_τ_μ_tikj[t][i][k][j] ) #+ 0.5 * v_false * (e_log_τj_err[j] -  logpi - e_τ_μ_tij_err[t][i][j])

                    # println("last")
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
function update_Nkj25_fast(N_rpip)
    Nkj = sum(sum.(N_rpip))#sum(Ntk)[1:end-1]
    return Nkj
end
function update_x_hat_sq_k25_fast(x,rpip)
    T = length(rpip)
    C_t = [length(el) for el in x]
    K = length(rpip[1][1])

    x_hat_sq_tik = [[[ rpip[t][i][k] .* (x[t][i]) .^2 for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_sq_tk = [[ sum(x_hat_sq_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_sq_k = [sum(x_hat_sq_tk[k]) for k in 1:K]
    return x_hat_sq_k
end
function update_x_hat_k25_fast(x,rpip)
    T = length(rpip)
    C_t = [length(el) for el in x]
    K = length(rpip[1][1])
    x_hat_tik = [[[ rpip[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_k = [sum(x_hat_tk[k]) for k in 1:K]
    return x_hat_k
end

function update_x_hat_sq_k25_fast(x,rtik,pip_kj;precomputed_x_sq=false)
    @views T = length(x)
    @views C_t = [length(el) for el in x]
    @views G = length(x[1][1])
    @views  K = length(pip_kj)
    # x_hat_tik = [[[ rpip[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_sq_k = Vector{Vector{Float64}}(undef, K)
     
    for k in 1:K
        x_hat_sq_kj = zeros(Float64,G)
        for t in 1:T
            @views  num_cells = C_t[t]
            for i in 1:num_cells
                if precomputed_x_sq
                    @views x_hat_sq_kj .+= x[t][i] .* rtik[t][i][k] .* pip_kj[k]
                else
                    @views x_hat_sq_kj .+= (x[t][i]) .^2 .* rtik[t][i][k] .* pip_kj[k]
                end
            end
        end
        x_hat_sq_k[k] = x_hat_sq_kj
    end
    return x_hat_sq_k
end
function update_x_hat_k25_fast(x,rtik,pip_kj)
    @views T = length(x)
    @views C_t = [length(el) for el in x]
    @views G = length(x[1][1])
    @views  K = length(pip_kj)
    # x_hat_tik = [[[ rpip[t][i][k] .* x[t][i] for i in 1:C_t[t]] for t in 1:T] for k in 1:K]
    # x_hat_tk = [[ sum(x_hat_tik[k][t]) for t in 1:T] for k in 1:K]
    x_hat_k = Vector{Vector{Float64}}(undef, K)
    for k in 1:K
        x_hat_kj = zeros(Float64,G)
        for t in 1:T
            @views  num_cells = C_t[t]
            # x_hat_kji = zeros(Float64,G)
            for i in 1:num_cells
                @views x_hat_kj .+= x[t][i] .* rtik[t][i][k] .* pip_kj[k]
            end
            # x_hat_kj .+= x_hat_kji
        end
        x_hat_k[k] = x_hat_kj
    end
    return x_hat_k
end
function τ_μ_expected_value25_fast(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(a0k_hat_vec)
    G = length(a0k_hat_vec[1])
    e_τ_μ_kj_true3 = Vector{Vector{Vector{Vector}}}(undef,T)
    e_τ_μ_true3 = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt3 =  Vector{Vector{Vector}}(undef,cells)
        e_τ_μ_13 =  Vector{Vector{Float64}}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti3 = Vector{Vector{Float64}}(undef,K)
            e_τ_μ_23 =  Vector{Float64}(undef,K)
            for k in 1:K
                e_τ_μ_kjti_cell3  =  a0k_hat_vec[k] ./  b0k_hat_vec[k] .*  (x[t][i] .- mk_hat_vec[k]) .^2 .+ 1 ./λ0k_hat_vec[k]
                e_τ_μ_3_cell3  = sum(e_τ_μ_kjti_cell3)
                
                e_τ_μ_kjti3[k] = e_τ_μ_kjti_cell3
                e_τ_μ_23[k] = e_τ_μ_3_cell3
            end
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end
function log_τ_kj_expected_value25_fast(a0k_hat_vec, b0k_hat_vec)
    K = length(a0k_hat_vec)
    e_log_τ_kj_vec = []
    for k in 1:K
        e_log_τ_kj = digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])
        # e_log_τ_k = sum(e_log_τ_kj)
        push!(e_log_τ_kj_vec,e_log_τ_kj)
    end
    return e_log_τ_kj_vec
end
function log_τ_k_expected_value25_fast(a0k_hat_vec, b0k_hat_vec)
    K = length(a0k_hat_vec)
    e_log_τ_kj_vec = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec)
    e_log_τ_k_vec = []
    for k in 1:K
        # e_log_τ_kj = digamma.(a0k_hat_vec_test[k]) .- log.(b0k_hat_vec[k])
        e_log_τ_k = sum(e_log_τ_kj_vec[k])
        push!(e_log_τ_k_vec,e_log_τ_k)
    end
    return e_log_τ_k_vec
end
function log_π_expected_value25_fast(θ_hat)
    T = length(θ_hat)
    K = Int64(length(θ_hat[1])-1)
    e_log_π = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        # digamma_sum = digamma(sum(θ_hat[t][1:K]))
        # e_log_π[t] = digamma.(θ_hat[t][1:K]) .- digamma_sum
        digamma_sum = digamma(sum(θ_hat[t]))
        e_log_π[t] = digamma.(θ_hat[t]) .- digamma_sum
    end
    
    # param_val_dict["log_π_expected_value"] = e_log_π
    return e_log_π
end
function expectation_log_π(θ_hat_t,θ_hat_t_sum)
    @views digamma_sum = digamma(θ_hat_t_sum)
    @views e_log_π_t =  digamma.(θ_hat_t) .- digamma_sum
    return e_log_π_t 
end

function expectation_log_τ_kj(a0k_hat_k, b0k_hat_k)
    @views e_log_τ_k = digamma.(a0k_hat_k) .- log.(b0k_hat_k)
    return e_log_τ_k
end
function expectation_τ_μ_j(x,λ0k_hat_k,mk_hat_k,a0k_hat_k, b0k_hat_k)
    @views e_τ_μ_jti  =  a0k_hat_k ./  b0k_hat_k .*  (x .- mk_hat_k) .^2 .+ 1 ./λ0k_hat_k
    return e_τ_μ_jti
end
function update_rtik_vs25_fast(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec)
    @views T = length(x)
    @views C_t = length.(x)
    @views G = length(x[1][1])
    @views K = length(a0k_hat_vec)
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    # ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    @views θ_hat_t_sum = sum.(θ_hat_vec)
    logpi = Glog/G
    for t in 1:T
        num_cells = C_t[t]
        @views adjusted_e_log_π_tk = sum([c_ttprime_vec[t][tt] .* expectation_log_π(θ_hat_vec[tt], θ_hat_t_sum[tt]) for tt in 1:t])
        @views rtik[t] = Vector{Vector{Float64}}(undef,num_cells)
        # @views ptik_tilde[t] = Vector{Vector{Float64}}(undef,num_cells)
        for i in 1:num_cells
            @views rtik[t][i] = Vector{Float64}(undef,K)
            @views ptik_tilde = Vector{Float64}(undef,K)
            for k in 1:K
                @views ptik_tilde[k] = adjusted_e_log_π_tk[k] .+ sum(0.5 .* (expectation_log_τ_kj(a0k_hat_vec[k], b0k_hat_vec[k]) .- logpi .- expectation_τ_μ_j(x[t][i],λ0k_hat_vec[k],mk_hat_vec[k],a0k_hat_vec[k], b0k_hat_vec[k])))
            end
            @views rtik[t][i] = norm_weights(ptik_tilde) 
        end
    end
    return rtik
end
function update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])
    K = length(a0k_hat_vec)
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    emptycells = [Vector{Vector{Float64}}(undef,C_t[t]) for t in 1:T]
    emptyclusters =  Vector{Float64}(undef,K)
    # ptik_tilde = Vector{Vector{Vector{Float64}}}(undef,T)
    θ_hat_t_sum = Vector{Float64}(undef,T)
    for t in 1:T
        θ_hat_t_sum[t] = sum(θ_hat_vec[t])
    end
    logpi = Glog/G
    for t in 1:T
        num_cells = C_t[t]
        conditions_pis_sums = 0.0
        adjusted_e_log_π_tk = adjust_e_log_π_tk(t,c_ttprime_vec,θ_hat_vec,θ_hat_t_sum)
        #adjusted_e_log_π_tk =sum([c_ttprime_vec[t][tt] .* expectation_log_π(θ_hat_vec[tt], θ_hat_t_sum[tt]) for tt in 1:t])
        rtik[t] = emptycells[t]
        for i in 1:num_cells
            rtik[t][i] = emptyclusters
            ptik_tilde = emptyclusters
            # @views cell_gene_sums = zeros(Float64,K)
            for k in 1:K
                cell_gene_sums  = 0.0
                # e_log_τ_kj = expectation_log_τ_kj2(a0k_hat_vec[k], b0k_hat_vec[k]) 
                e_log_normal_l_j = expectation_log_normal_l_j2(x[t][i],λ0k_hat_vec[k],mk_hat_vec[k],a0k_hat_vec[k], b0k_hat_vec[k],logpi)
                cell_gene_sums =+ @views sum(e_log_normal_l_j)#sum_vec_el_in_loop(e_log_normal_l_j)
                ptik_tilde[k] = adjusted_e_log_π_tk[k] + cell_gene_sums#adjust_e_log_π_tk(t,k,c_ttprime_vec,θ_hat_vec,θ_hat_t_sum) + cell_gene_sums#
            end
            rtik[t][i] = norm_weights(ptik_tilde) 
        end
    end
    return rtik
end

function norm_weights2(p)
    K = length(p)
    psum = StatsFuns.logsumexp(p)
    w = Vector{Float64}(undef,K)
    for k in 1:K
        w[k] = exp(p[k] - psum)
    end
    
    return w
end
# Memory estimate: 83.16 MiB, allocs estimate: 800023.
function expectation_log_π2(θ_hat_t,θ_hat_t_sum)
    K = length(θ_hat_t)
    e_log_π_t = Vector{Float64}(undef,K)
    digamma_sum = digamma(θ_hat_t_sum)
    for k in 1:K
        @views e_log_π_t[k] =  digamma(θ_hat_t[k]) - digamma_sum
    end
    return e_log_π_t 
end
function adjust_e_log_π_tk(t,c_ttprime_vec,θ_hat_vec,θ_hat_t_sum)
    Kplus = length(θ_hat_vec[t])
    conditions_pis_sums = 0.0
    adjusted_e_log_π_tk = Vector{Float64}(undef,Kplus)
    for k in 1:Kplus
        conditions_pis_sums = 0.0
        for tt in 1:t
            conditions_pis_sums += c_ttprime_vec[t][tt] * (digamma(θ_hat_vec[tt][k]) - digamma(θ_hat_t_sum[tt]))
        end
        adjusted_e_log_π_tk[k] = conditions_pis_sums
    end
    return adjusted_e_log_π_tk
end
function adjust_e_log_π_tk(t,k,c_ttprime_vec,θ_hat_vec,θ_hat_t_sum)
    conditions_pis_sums = 0.0
    for tt in 1:t
        conditions_pis_sums += c_ttprime_vec[t][tt] * (digamma(θ_hat_vec[tt][k]) - digamma(θ_hat_t_sum[tt]))
    end
    return conditions_pis_sums
end
function expectation_log_τ_kj2(a0k_hat_k, b0k_hat_k)
    G = length(a0k_hat_k[1])
    e_log_τ_k  = Vector{Float64}(undef,G)
    for j in 1:G
        @views e_log_τ_k[j]  =  digamma(a0k_hat_k[j]) - log(b0k_hat_k[j]) 
    end
    return e_log_τ_k
end
function expectation_τ_μ_j2(x,λ0k_hat_k,mk_hat_k,a0k_hat_k, b0k_hat_k)
    G = length(a0k_hat_k[1])
    e_τ_μ_jti  = Vector{Float64}(undef,G)
    for j in 1:G
        @views e_τ_μ_jti[j]  =  a0k_hat_k[j]  /  b0k_hat_k[j]  *  (x[j]  - mk_hat_k[j] ) ^2 + 1 /λ0k_hat_k[j] 
    end
    return e_τ_μ_jti
end
function expectation_log_normal_l_j2(x,λ0k_hat_k,mk_hat_k,a0k_hat_k, b0k_hat_k,logpi)
    G = length(a0k_hat_k)
    # println(G)
    log_normal_l_j = Vector{Float64}(undef,G)
    for j in 1:G
        @views log_normal_l_j[j]  =  0.5 * (digamma(a0k_hat_k[j]) - log(b0k_hat_k[j]) - ( a0k_hat_k[j]  /  b0k_hat_k[j]  *  (x[j]  - mk_hat_k[j] ) ^2 + 1 /λ0k_hat_k[j]  + logpi))
    end
    return log_normal_l_j
end
function sum_vec_el_in_loop(val_vec)
    sum_val = zero(eltype(val_vec))
    num_values = length(val_vec)
    for i in 1:num_values
        sum_val += val_vec[i]
    end
    return sum_val
end
function update_N_rpip25_fast(rtik,pip_kj)
    @views T = length(rtik)
    @views K = length(rtik[1][1])
    @views C_t = length.(rtik)
    @views G = length(pip_kj[1])
    N_rpip = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        num_cells=C_t[t]
        @views N_rpip[t] = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            @views N_rpip[t][i] = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                @views N_rpip[t][i][k] = Vector{Float64}(undef,G)
            end
        end
    end
    # println("Breaking")
    # @views N_rpip = [[Vector{Vector{Float64}}(undef,K)  for i in 1:C_t[t] ]  t in 1:T]
    N_rpip = fill_N_rpip25_fast!(N_rpip, rtik, pip_kj, T,K,C_t)
    # println("weird")
    return N_rpip
end
function fill_N_rpip25_fast!(N_rpip, rtik, pip_kj, T,K,C_t)
    # @views T = length(rtik)
    # @views K = length(rtik[1][1])
    # @views C_t = length.(rtik)
    @views G = length(pip_kj[1])

    for t in 1:T
        @views num_cells=C_t[t]
        # @views N_rpip[t] = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            # @views N_rpip[t][i] = Vector{Vector{Float64}}(undef,K)
            # N_rpip[t][i] = Vector{Float64}(undef,G)
            for k in 1:K
                @views N_rpip[t][i][k] .= rtik[t][i][k] .* pip_kj[k]
                #
            end
        end
    end
    #     rtik_flat = innermelt(recursive_flatten(deepcopy(rtik)),G)
    #     pip_kj_flat = outermelt(recursive_flatten(deepcopy(pip_kj)),sum(C_t)) .* innermelt(recursive_flatten(deepcopy(rtik)),G)
    #     num_el = length(rtik_flat)
    #    for pip in pip_kj 
    #     for r in rtik
    #         for rti in r

    #         end
    #     end
    return N_rpip
end
function update_Nkj25_fast(rtik, pip_kj)
    @views T = length(rtik)
    @views K = length(rtik[1][1])
    @views C_t = length.(rtik)
    @views G = length(pip_kj[1])
    Nkj = Vector{Vector{Float64}}(undef,K)
    for k in 1:K
        @views rik_pip = zeros(Float64,G)#Vector{Vector{Float64}}(undef,num_cells)
        for t in 1:T
            num_cells=C_t[t]
            
            for i in 1:num_cells
                # @views N_rpip[t][i] = Vector{Vector{Float64}}(undef,K)
                @views rik_pip  .+= rtik[t][i][k] .* pip_kj[k]
            end
        end
        @views Nkj[k] = rik_pip#sum(rik_pip)
    end
    # Nkj = sum(sum.( update_N_rpip25_fast(rtik, pip_kj)))#sum(Ntk)[1:end-1]
    return Nkj
end

function update_N_rpip25_fast(rtik,pip_kj)
    T = length(rtik)
    K = length(rtik[1][1])
    C_t = [length(el) for el in rtik]
    G = length(pip_kj[1][1])
    N_signal = Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        cells_=C_t[t]
        Nt_signal = Vector{Vector{Vector{Float64}}}(undef,cells_)
        for i in 1:cells_
            Nti_signal = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                Ntik_signal =  rtik[t][i][k] .* pip_kj[k]
                # for j in 1:G
                #     Ntik_signal[j] =
                # end
                Nti_signal[k] = Ntik_signal
            end
            Nt_signal[i] = Nti_signal
        end
        N_signal[t] = Nt_signal
    end
    return N_signal
end
function calc_DataElbo25_fast(x,rpip,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
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
    weighted_ss_kjti = [[[[rpip[t][i][k][j] .*(x[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
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
function calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(a0k_hat_vec)
    G = length(a0k_hat_vec[1])
    # halfNklogπ = 1/2 .* Nk .*log(2π) # K by 1
    # halfNklogπ = -1/2 .* Nk .*log(2π) # K by 1
    # halfNkoverλ0k_hat = 1/2 .* [ Nk[k] ./λ0k_hat_vec[k] for k in 1:K] # K by G
    # a0b0_c_Ga =  [c_Ga.(a0_vec, b0_vec) for k in 1:K] # K by G [zeros(G) for k in 1:K]#
    # a0kb0k_c_Ga=  c_Ga.(a0k_hat_vec, b0k_hat_vec) # K by G [zeros(G) for k in 1:K]#

    
    
    # halfLogλ0 = 1/2 .* [ log.(λ0_vec) for k in 1:K] # K by G
    # halfLogλ0k = 1/2 .* [ log.(λ0k_hat_vec[k]) for k in 1:K] # K by G
    # halfλ0μ0_sq = 1/2 .* λ0_vec.* μ0_vec .^ 2 #G By 1
    # a0kOverb0k_hat = [a0k_hat_vec[k] ./ b0k_hat_vec[k] for k in 1:K] # K by G
    # bs_e_τ = [(b0k_hat_vec[k]  .- b0_vec .- ( 1/2 .* λ0_vec.* μ0_vec .^ 2 )) .*(a0k_hat_vec[k] ./ b0k_hat_vec[k] ) for k in 1:K] #K by G


    # as_e_log_τkj =[(1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])) for k in 1:K] #K by G
    # λ0μ0 = λ0_vec .* μ0_vec #G By 1
    # λ0μ0mk_a0kOverb0k_hat = [ λ0_vec .* μ0_vec .* mk_hat_vec[k] .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) for k in 1:K] #K by G
    # halfλ0mk_sq_a0kOverb0k_hat = [1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) for k in 1:K] #K by G
    one_half_const = 1/2
    # halfλ0overλ0k_hat = [1/2 .* λ0_vec ./ λ0k_hat_vec[k] for k in 1:K] #K by G
    data_elbo = 0.0
    for k in 1:K
        weighted_ss_kj = zeros(Float64,G)
        for t in 1:T
            for i in 1:C_t[t]
                @views weighted_ss_kj .+= rtik[t][i][k] .* pip_kj[k] .*(x[t][i] .-   mk_hat_vec[k]) .^ 2
            end
        end
        # weighted_ss[k] = weighted_ss_kj
        # -1/2 .* Nk[k] .*log(2π) .- 1/2 .* (Nk[k] ./λ0k_hat_vec[k]) .+ c_Ga.(a0_vec, b0_vec) .- c_Ga.(a0k_hat_vec, b0k_hat_vec)[k] .+ 1/2 .* log.(λ0_vec)  .- 1/2 .* log.(λ0k_hat_vec[k]) .+ ((b0k_hat_vec[k]  .- b0_vec .- ( 1/2 .* λ0_vec.* μ0_vec .^ 2 )) .*(a0k_hat_vec[k] ./ b0k_hat_vec[k] ))  .+ ((1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k]))) .+ λ0_vec .* μ0_vec .* mk_hat_vec[k] .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) .- 1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) .+ one_half_const .+ 1/2 .*( λ0_vec ./ λ0k_hat_vec[k]) .- 1/2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ).* weighted_ss_kj
        @views data_elbo += sum(-1/2 .* Nkj[k] .*log(2π) .- 1/2 .* (Nkj[k] ./λ0k_hat_vec[k]) .+ c_Ga.(a0_vec, b0_vec) .- c_Ga.(a0k_hat_vec, b0k_hat_vec)[k] .+ 1/2 .* log.(λ0_vec)  .- 1/2 .* log.(λ0k_hat_vec[k]) .+ ((b0k_hat_vec[k]  .- b0_vec .- ( 1/2 .* λ0_vec.* μ0_vec .^ 2 )) .*(a0k_hat_vec[k] ./ b0k_hat_vec[k] ))  .+ ((1/2 .* Nkj[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k]))) .+ λ0_vec .* μ0_vec .* mk_hat_vec[k] .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) .- 1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) .+ one_half_const .+ 1/2 .*( λ0_vec ./ λ0k_hat_vec[k]) .- 1/2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ).* weighted_ss_kj)
    end
    # weighted_ss_kjti = [[[[rtik[t][i][k] .* pip_kj[k][j] .*(x[t][i][j] .-   mk_hat_vec[k][j]) .^ 2 for i in 1:C_t[t] ] for t in 1:T] for j in 1:G] for k in 1:K]
    # weighted_ss_kjt = [[[sum(weighted_ss_kjti[k][j][t]) for t in 1:T] for j in 1:G] for k in 1:K]
    # weighted_ss_kj = [[sum(weighted_ss_kjt[k][j]) for j in 1:G] for k in 1:K] #K by G
    # halfa0kOverb0k_hat_weighted_ss_kj =  [ 1/2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ).* weighted_ss_kj[k]  for k in 1:K] #K by G

    # # data_elbo_kj = [ halfNklogπ[k] .- halfNkoverλ0k_hat[k] .+ a0b0_c_Ga[k] .+ a0b0_c_Ga[k] .+ halfLogλ0[k]  .- halfLogλ0k[k]  .+ bs_e_τ[k]  .+ as_e_log_τkj[k] .+ λ0μ0mk_a0kOverb0k_hat[k] .- halfλ0mk_sq_a0kOverb0k_hat[k] .- one_half_const .+ halfλ0overλ0k_hat[k] .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    # data_elbo_kj = [-1/2 .* Nk[k] .*log(2π) .- 1/2 .* (Nk[k] ./λ0k_hat_vec[k]) .+ c_Ga.(a0_vec, b0_vec) .- c_Ga.(a0k_hat_vec, b0k_hat_vec)[k] .+ 1/2 .* log.(λ0_vec)  .- 1/2 .* log.(λ0k_hat_vec[k]) .+ ((b0k_hat_vec[k]  .- b0_vec .- ( 1/2 .* λ0_vec.* μ0_vec .^ 2 )) .*(a0k_hat_vec[k] ./ b0k_hat_vec[k] ))  .+ ((1/2 .* Nk[k] .+  a0_vec .- a0k_hat_vec[k]) .* (digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k]))) .+ λ0_vec .* μ0_vec .* mk_hat_vec[k] .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) .- 1/2 .* λ0_vec .* (mk_hat_vec[k]) .^2 .* (a0k_hat_vec[k] ./ b0k_hat_vec[k] ) .+ one_half_const .+ 1/2 .*( λ0_vec ./ λ0k_hat_vec[k]) .- halfa0kOverb0k_hat_weighted_ss_kj[k] for k in 1:K]
    # data_elbo_j  = sum(data_elbo_kj)
    # data_elbo  = sum(data_elbo_j)
    #     
    return data_elbo
end
function update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    K = length(Nkj) 
    a0k_hat_vec = [ a0_vec .+ 1/2 .* (Nkj[k] .+1 ) for k in 1:K]
    return a0k_hat_vec
end
function log_norm_pdf(x,μk,a0k,b0k) 
    return 1/2 .* log.(a0k./b0k) .-  1/2 .* ((x.-μk).^2 .* (a0k./b0k)) .- 1/2 .* log.(2π)
end
function expectation_log_norm_pdf(x,μk,λ0k,a0k,b0k) 
    return 1/2 .* expectation_log_τ_kj(a0k, b0k) .-  1/2 .* expectation_τ_μ_j(x,λ0k,μk,a0k,b0k) .- 1/2 .* log.(2π)
end
function log_norm_pdf(x,μk,a0k,b0k) 
    return 1/2 .* log.(a0k./b0k) .-  1/2 .* ((x.-μk).^2 .* (a0k./b0k)) .- 1/2 .* log.(2π)
end
function expectation_log_norm_pdf(x,μk,λ0k,a0k,b0k) 
    return 1/2 .* expectation_log_τ_kj(a0k, b0k) .-  1/2 .* expectation_τ_μ_j(x,λ0k,μk,a0k,b0k) .- 1/2 .* log.(2π)
end

function get_gene_PIP25_fast(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
    @views G = length(x[1][1])
    @views T = length(x)
    @views C_t = length.(x)
    @views K = length(rtik[1][1])

    N_k = sum(sum.(rtik))

    # my_ll_scores = [[[log_norm_pdf(x[t][i],mk_hat_vec[k],a0k_hat_vec[k] , b0k_hat_vec[k]) for k in 1:K ] for i in 1:C_t[t]] for t in 1:T] 
    # my_expectation_ll_scores = [[[expectation_log_norm_pdf(x[t][i],mk_hat_vec[k],λ0k_hat_vec[k],a0k_hat_vec[k] , b0k_hat_vec[k]) for k in 1:K ] for i in 1:C_t[t]] for t in 1:T] 
    # all([all([all([all(cell_ll_scores[t][i][k] .≈ my_ll_scores[t][i][k]) for k in 1:K]) for i in 1:C_t[t]]) for t in 1:T])
    # expectation_log_τ_kj(a0k_hat_k, b0k_hat_k)
    # expectation_τ_μ_j(x,λ0k_hat_k,mk_hat_k,a0k_hat_k, b0k_hat_k)

    
    # expected_val_cell_cell_ll_score = [[[[cell_ll_scores[t][i][k][j] for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    # ration_ = [[[[expected_val_cell_cell_ll_score[t][i][j][k] .- (null_cell_ll_scores[t][i][k][j] )  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    # ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # ration1_weight = [[[norm_weights(ration1_[t][i][k]) .* rtik[t][i][k] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];


    # ration2_ = [[[[ration_[t][i][j][k] .+ log( rtik[t][i][k]) for k in 1:K] for j in 1:G]  for i in 1:C_t[t]] for t in 1:T];
    # ration2_weight = [[[norm_weights(ration2_[t][i][j]) .* rtik[t][i] for j in 1:G]  for i in 1:C_t[t]] for t in 1:T];
    
    # 
    # cell_ll_scores
    # gg= [[[[cell_ll_scores[t][i][k][j] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    # normalized_weights = Vector{Vector}(undef,K)
    pip_kj  = Vector{Vector}(undef,K)
    # println("here")
    for k in 1:K
        unnormalized_weights = zeros(Float64,G)
        for t in 1:T
            for i in 1:C_t[t]
                # @views unnormalized_weights .+= norm_weights(cell_ll_scores[t][i][k] .- null_cell_ll_scores[t][i][k]) .* rtik[t][i][k]
                @views unnormalized_weights .+= norm_weights(log_norm_pdf(x[t][i],mk_hat_vec[k],a0k_hat_vec[k] , b0k_hat_vec[k]) .- log_norm_pdf(x[t][i],0.0,null_precision , 1.0)) .* rtik[t][i][k]
            end
        end
        # normalized_weights[k] = unnormalized_weights ./ N_k[k]
        # if any(isnan.(normalized_weights[k])) || any(iszero.(normalized_weights[k])) 
        #     if all(isnan.(normalized_weights[k])) || all(iszero.(normalized_weights[k]))
        #         normalized_weights[k] .= ones(Float64,G)
        #     else
        #         normalized_weights[k][isnan.(normalized_weights[k])] .= 0.0
        #     end
        # end
        @views pip_kj[k] = normToProb(fix_nan_or_allzero2!(unnormalized_weights ./ N_k[k]))
    end
    # println("now here")
    
    
    # gene_significance_weights = sum(sum.(ration1_weight)) ./ N_k

    # sum(sum.(ration2_weight))
    # [el ./ N_k for el in sum(sum.(ration2_weight))]
    # normToProb.([el ./ N_k for el in sum(sum.(ration2_weight))])

    
    # gene_significance_weights_mat = hcat(gene_significance_weights...)

    # avg_gene_significance_weights_mat = mean(gene_significance_weights_mat, dims=3)
    # w_kj = [avg_gene_significance_weights_mat[:,col,1] for col in 1:size(avg_gene_significance_weights_mat)[2]]


    # w_kj = normToProb.([ nan_or_zero_to_1.(el) for el in gene_significance_weights])

    # pip_kj = normToProb.( fix_nan_or_allzero!(normalized_weights))
    return pip_kj
end
nan_or_zero_to_1(v) = nan_or_zero(v) ? 1.0 : v 
nan_or_zero(v) = iszero(v) || isnan(v)  ? true : false 
function fix_nan_or_allzero2!(v)
    G = length(v)
    if any(isnan.(v)) || any(iszero.(v)) 
        if all(isnan.(v)) || all(iszero.(v))
            v .= ones(Float64,G)
        else
            v[isnan.(v)] .= 0.0
        end
    end
    return v
end

function fix_nan_or_allzero!(v)
    K = length(v)
    G = length(v[1])
    for k in 1:K
        if any(isnan.(v[k])) || any(iszero.(v[k])) 
            if all(isnan.(v[k])) || all(iszero.(v[k]))
                v[k] .= ones(Float64,G)
            else
                v[k][isnan.(v[k])] .= 0.0
            end
        end
    end
    return v
end

function get_gene_PIP25_fast2_depracated(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)#,N_k=nothing
    G = length(x[1][1])
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    pip_kj  = Vector{Vector{Float64}}(undef,K)
    N_k = Vector{Float64}(undef,K)
    # unnormalized_weights  = Vector{Vector{Float64}}(undef,K)
    # # println("here")
    # for k in 1:K
    #     unnormalized_weights[k]=zeros(Float64,G)
    # end
    gene_vec = Vector{Float64}(undef,G)
    zeros_vec = zeros(Float64,G)
    m_null = zeros(Float64,G)
    a0_null = ones(Float64,G)
    for j in 1:G
        a0_null[j] = null_precision * a0_null[j]
    end
    b0_null = ones(Float64,G)
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
function get_gene_PIP25_fast2(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)#,N_k=nothing
    G = length(x[1][1])
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(rtik[1][1])
    pip_kj  = Vector{Vector{Float64}}(undef,K)
    N_k = Vector{Float64}(undef,K)
    # unnormalized_weights  = Vector{Vector{Float64}}(undef,K)
    # # println("here")
    # for k in 1:K
    #     unnormalized_weights[k]=zeros(Float64,G)
    # end
    gene_vec = Vector{Float64}(undef,G)
    zeros_vec = zeros(Float64,G)
    m_null = zeros(Float64,G)
    a0_null = ones(Float64,G)
    for j in 1:G
        a0_null[j] = null_precision * a0_null[j]
    end
    b0_null = ones(Float64,G)
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
        pip_kj[k] = normToProb(normalized_weights)
    end
    return pip_kj
end
function accumulate_deviance_log_norm_pdf2!(unnormalized_weights,x,mk_hat_k,a0k_hat_k, b0k_hat_k,m_null,a0_null , b0_null,rtik)
    G = length(x)
    # logpi =  log(2π)
    log_norm_pdf_ = Vector{Float64}(undef,G)
    for j in 1:G
        # log_norm_pdf_[j] = (1/2 * log(a0k_hat_k[j]/b0k_hat_k[j]) -  1/2 * ((x[j]-mk_hat_k[j])^2 * (a0k_hat_k[j]/b0k_hat_k[j])) - 1/2 * logpi) - (1/2 * log(a0_null[j]/b0_null[j]) -  1/2 * ((x[j]-m_null[j])^2 * (a0_null[j]/b0_null[j])) - 1/2 * logpi) 
        log_norm_pdf_[j] = log_norm_pdf2(x[j],mk_hat_k[j],a0k_hat_k[j],b0k_hat_k[j])- log_norm_pdf2(x[j],m_null[j],a0_null[j],b0_null[j])
    end
    normed_weight = norm_weights2(log_norm_pdf_)
    for j in 1:G
        unnormalized_weights[j] +=  normed_weight[j]* rtik
    end
    return unnormalized_weights
end
function log_norm_pdf2_depracated(x,μk,a0k,b0k)
    G = length(x)
    log_norm_pdf_ = Vector{Float64}(undef,G)
    for j in 1:G
        log_norm_pdf_[j] = 1/2 * log(a0k[j]/b0k[j]) -  1/2 * ((x[j]-μk[j])^2 * (a0k[j]/b0k[j])) - 1/2 * log(2π)
    end
    return log_norm_pdf_
end
function log_norm_pdf2(x,μk,a0k,b0k)
    return 1/2 * log(a0k/b0k) -  1/2 * ((x-μk)^2 * (a0k/b0k)) - 1/2 * log(2π)
end
function fix_nan_or_allzero22!(v)
    G = length(v)
    if any(isnan.(v)) || any(iszero.(v)) 
        if all(isnan.(v)) || all(iszero.(v))
            v .= ones(Float64,G)
        else
            v[isnan.(v)] .= 0.0
        end
    end
    return v
end


function variational_inference_dynamicHDP_vs25_fast(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,pip_kj_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, null_precision, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


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
    pip_kj_init = init_pip_kj_vec!(pip_kj_init,G,K;rand_init = false)

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
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        Nk_chain = make_chain(num_iter+1,Nk)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)

        
        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec

        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:Nk_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
    Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            # e_log_π = log_π_expected_value25_fast(θ_hat_vec) # T by K
            # # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            # e_log_τkj = log_τ_kj_expected_value25_fast(a0k_hat_vec, b0k_hat_vec);
            # e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value25_fast(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik_vs25_fast(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end
        # rpip = update_N_rpip25_fast(rtik,pip_kj)
        Nkj = update_Nkj25_fast(rtik, pip_kj);
        x_hat_k = update_x_hat_k25_fast(x,rtik,pip_kj);
        x_hat_sq_k = update_x_hat_sq_k25_fast(x,rtik,pip_kj);
        if record_chain
            chain_dict[:Nk_chain][iter + 1] = Nk
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
        a0k_hat_vec = update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec
        end

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        pip_kj =  get_gene_PIP25_fast(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=null_precision);

        data_elbo =  calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec);
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        pip_entropy = calc_Hpip(pip_kj);

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + pip_entropy
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
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end

function variational_inference_dynamicHDP_vs25_fast2(inputs_dict;mk_hat_init=nothing, λ0k_hat_init=nothing,a0k_hat_init=nothing, b0k_hat_init=nothing,awt_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,θ_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,pip_kj_init=nothing, rhok_hat_init=nothing, omegak_hat_init=nothing,uniform_theta_init = true, rand_init = false,ep = 0.001,elbo_ep = 10^(-6),record_chain = false)
    x, K, a0, b0, μ0, λ0, a_γ, b_γ, a_α, b_α, adot_w, bdot_w, null_precision, num_iter, num_local_iter = (; inputs_dict...)
    T = length(x)
    G = length(x[1][1])
    C_t = [length(el) for el in x]
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)


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
    pip_kj_init = init_pip_kj_vec!(pip_kj_init,G,K;rand_init = false)

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
    chain_dict = nothing
    arg_str_list_initparams = @name K , mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ;
    key_list_initparams = Symbol.(naming_vec(arg_str_list_initparams));
    var_list_initparams = [K, mk_hat_init,λ0k_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init,omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init, bwt_hat_init,a_αt_hat_init,b_αt_hat_init,θ_hat_init,c_ttprime_init,rtik_init, pip_kj_init,null_precision, num_iter, num_local_iter,uniform_theta_init, rand_init,ep,elbo_ep,record_chain ];
    
    initDict = OrderedDict()
    addToDict!(initDict,key_list_initparams,var_list_initparams);
    
    #init debug dict
    if record_chain
        e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        e_log_τ = Vector{Vector{Float64}}(undef,num_local_iter)
        e_τ_μ_tikj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef,num_local_iter)
        e_τ_μ = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Ntk = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        θ_hat_vec = Vector{Vector{Float64}}(undef,num_local_iter)
        c_ttprime_vec = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        rtik = Vector{Vector{Vector{Vector{Float64}}}}(undef,num_local_iter)
        Nk = Vector{Float64}()
        x_hat_k = Vector{Vector{Float64}}()
        x_hat_sq_k = Vector{Vector{Float64}}()
        a_αt_hat_vec = Vector{Float64}()
        b_αt_hat_vec = Vector{Float64}()
        awt_hat_vec = Vector{Float64}()
        bwt_hat_vec = Vector{Float64}()
        a_γ_hat,b_γ_hat = 1.0,1.0
        e_γ = 1.0
        Tαk = Vector{Float64}()
        data_elbo = 1.
        assgn_entropy = 1.
        dHDP_surragate_elbo =1. 
        s_entropy = 1
        wAlloc_elbo = 1.
        γ_elbo =1. 
        α_elbo = 1.
        λ0k_chain = make_chain(num_iter+1,λ0k_hat_vec)
        mk_chain = make_chain(num_iter+1,mk_hat_vec)
        a0k_chain = make_chain(num_iter+1,a0k_hat_vec)
        b0k_chain = make_chain(num_iter+1,b0k_hat_vec)
        rhok_chain = make_chain(num_iter+1,rhok_hat_vec)
        omegak_chain = make_chain(num_iter+1,omegak_hat_vec)
        θ_hat_chain = make_chain(num_iter+1,θ_hat_vec)
        rtik_chain = make_chain(num_iter+1,rtik)
        e_log_π_chain = make_chain(num_iter+1,e_log_π)
        e_log_τ_chain = make_chain(num_iter+1,e_log_τ)
        e_τ_μ_tikj_chain = make_chain(num_iter+1,e_τ_μ_tikj)
        e_τ_μ_chain = make_chain(num_iter+1,e_τ_μ)
        Ntk_chain = make_chain(num_iter+1,Ntk)
        c_ttprime_chain = make_chain(num_iter+1,c_ttprime_vec)
        Nk_chain = make_chain(num_iter+1,Nk)
        x_hat_k_chain = make_chain(num_iter+1,x_hat_k)
        x_hat_sq_k_chain = make_chain(num_iter+1,x_hat_sq_k)
        a_αt_hat_chain = make_chain(num_iter+1,a_αt_hat_vec)
        b_αt_hat_chain = make_chain(num_iter+1,b_αt_hat_vec)
        awt_hat_chain = make_chain(num_iter+1,awt_hat_vec)
        bwt_hat_chain = make_chain(num_iter+1,bwt_hat_vec)
        a_γ_hat_chain = make_chain(num_iter+1,a_γ_hat)
        b_γ_hat_chain = make_chain(num_iter+1,b_γ_hat)
        e_γ_chain = make_chain(num_iter+1,e_γ)
        Tαk_chain = make_chain(num_iter+1,Tαk)
        data_elbo_chain = make_chain(num_iter+1,data_elbo)
        assgn_entropy_chain = make_chain(num_iter+1,assgn_entropy)
        dHDP_surragate_elbo_chain = make_chain(num_iter+1,dHDP_surragate_elbo)
        s_entropy_chain = make_chain(num_iter+1,s_entropy)
        wAlloc_elbo_chain = make_chain(num_iter+1,wAlloc_elbo)
        γ_elbo_chain = make_chain(num_iter+1,γ_elbo)
        α_elbo_chain = make_chain(num_iter+1,α_elbo)

        
        arg_str_list_chain = @name λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        key_list_chain= Symbol.(naming_vec(arg_str_list_chain));
        var_list_chain = [λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain];

        chain_dict = OrderedDict()
        addToDict!(chain_dict,key_list_chain,var_list_chain);
        # e_log_π = Vector{Vector{Vector{Float64}}}(undef,num_local_iter)
        
    end
    #init debug dict initial values
    if record_chain
        # λ0k_chain, mk_chain,a0k_chain,b0k_chain,rhok_chain,omegak_chain,θ_hat_chain,rtik_chain,e_log_π_chain,e_log_τ_chain,e_τ_μ_tikj_chain,e_τ_μ_chain,Ntk_chain,c_ttprime_chain,Nk_chain,x_hat_k_chain,x_hat_sq_k_chain,a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain;
        chain_dict[:λ0k_chain][1] = λ0k_hat_vec
        chain_dict[:mk_chain][1] = mk_hat_vec
        chain_dict[:a0k_chain][1] = a0k_hat_vec
        chain_dict[:b0k_chain][1] = b0k_hat_vec
        chain_dict[:rhok_chain][1] = rhok_hat_vec
        chain_dict[:omegak_chain][1] = omegak_hat_vec
        chain_dict[:θ_hat_chain][1] = θ_hat
        chain_dict[:rtik_chain][1] = rtik
        chain_dict[:c_ttprime_chain][1] = c_ttprime_vec
        chain_dict[:a_αt_hat_chain][1] = a_αt_hat_vec
        chain_dict[:b_αt_hat_chain][1] = b_αt_hat_vec
        chain_dict[:a_γ_hat_chain][1] = a_γ_hat 
        chain_dict[:b_γ_hat_chain][1] = b_γ_hat
        chain_dict[:awt_hat_chain][1] = awt_hat_vec
        chain_dict[:bwt_hat_chain][1] = bwt_hat_vec

        chain_dict[:e_log_π_chain][1] = nothing
        chain_dict[:e_log_τ_chain][1] = nothing
        chain_dict[:e_τ_μ_tikj_chain][1] = nothing
        chain_dict[:e_τ_μ_chain][1] = nothing
        chain_dict[:Ntk_chain][1] = nothing
        chain_dict[:Nk_chain][1] = nothing
        chain_dict[:x_hat_k_chain][1] = nothing
        chain_dict[:x_hat_sq_k_chain][1] = nothing
        chain_dict[:e_γ_chain][1] = nothing
        chain_dict[:Tαk_chain][1] = nothing
        chain_dict[:data_elbo_chain][1] = nothing
        chain_dict[:assgn_entropy_chain][1] = nothing
        chain_dict[:dHDP_surragate_elbo_chain][1] = nothing
        chain_dict[:s_entropy_chain][1] = nothing
        chain_dict[:wAlloc_elbo_chain][1] = nothing
        chain_dict[:γ_elbo_chain][1] = nothing
        chain_dict[:α_elbo_chain][1] = nothing
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
    Glog = G*log(2π)
    while !converged_bool #for iter in 1:num_iter
        for loc_iter in 1:num_local_iter
            # e_log_π = log_π_expected_value25_fast(θ_hat_vec) # T by K
            # # e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
            # e_log_τkj = log_τ_kj_expected_value25_fast(a0k_hat_vec, b0k_hat_vec);
            # e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value25_fast(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
            rtik = update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec);
            Ntk = update_Ntk(rtik)
            c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
            θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
            if record_chain
                chain_dict[:θ_hat_chain][iter + 1][loc_iter] = θ_hat
                chain_dict[:rtik_chain][iter + 1][loc_iter] = rtik
                chain_dict[:c_ttprime_chain][iter + 1][loc_iter] = c_ttprime_vec
                chain_dict[:e_log_π_chain][iter + 1][loc_iter] = e_log_π
                chain_dict[:e_log_τ_chain][iter + 1][loc_iter] = e_log_τ
                chain_dict[:e_τ_μ_tikj_chain][iter + 1][loc_iter] = e_τ_μ_tikj
                chain_dict[:e_τ_μ_chain][iter + 1][loc_iter] = e_τ_μ
                chain_dict[:Ntk_chain][iter + 1][loc_iter] = Ntk
            end
        end
        # rpip = update_N_rpip25_fast(rtik,pip_kj)
        Nkj = update_Nkj25_fast(rtik, pip_kj);
        x_hat_k = update_x_hat_k25_fast(x,rtik,pip_kj);
        x_hat_sq_k = update_x_hat_sq_k25_fast(x,rtik,pip_kj);
        if record_chain
            chain_dict[:Nk_chain][iter + 1] = Nk
            chain_dict[:x_hat_k_chain][iter + 1] = x_hat_k
            chain_dict[:x_hat_sq_k_chain][iter + 1] = x_hat_sq_k
        end
        
        λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nkj)
        mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
        a0k_hat_vec = update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
        b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
        if record_chain
            chain_dict[:λ0k_chain][iter + 1] = λ0k_hat_vec
            chain_dict[:mk_chain][iter + 1] = mk_hat_vec
            chain_dict[:a0k_chain][iter + 1] = a0k_hat_vec
            chain_dict[:b0k_chain][iter + 1] = b0k_hat_vec
        end

        a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
        awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
        bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
        a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
        e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
        Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
        rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
        if record_chain
            chain_dict[:a_αt_hat_chain][iter + 1] = a_αt_hat_vec
            chain_dict[:b_αt_hat_chain][iter + 1] = b_αt_hat_vec
            chain_dict[:awt_hat_chain][iter + 1] = awt_hat_vec
            chain_dict[:bwt_hat_chain][iter + 1] = bwt_hat_vec
            chain_dict[:a_γ_hat_chain][iter + 1] = a_γ_hat
            chain_dict[:b_γ_hat_chain][iter + 1] = b_γ_hat
            chain_dict[:e_γ_chain][iter + 1] = e_γ
            chain_dict[:Tαk_chain][iter + 1] = Tαk
            chain_dict[:rhok_chain][iter + 1] = rhok_hat_vec
            chain_dict[:omegak_chain][iter + 1] = omegak_hat_vec
        end
        # a_αt_hat_chain,b_αt_hat_chain,awt_hat_chain,bwt_hat_chain,a_γ_hat_chain,b_γ_hat_chain,e_γ_chain,Tαk_chain,data_elbo_chain,assgn_entropy_chain,dHDP_surragate_elbo_chain,s_entropy_chain,wAlloc_elbo_chain,γ_elbo_chain,α_elbo_chain
        # a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat
        pip_kj =  get_gene_PIP25_fast2(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=null_precision);

        data_elbo =  calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec);
        assgn_entropy =  calc_Hz(rtik) 
        dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
        s_entropy = calc_Hs(c_ttprime_vec)
        wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
        γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
        α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
        pip_entropy = calc_Hpip(pip_kj);

        if record_chain
            chain_dict[:data_elbo_chain][iter + 1] = data_elbo
            chain_dict[:assgn_entropy_chain][iter + 1] = assgn_entropy
            chain_dict[:dHDP_surragate_elbo_chain][iter + 1] = dHDP_surragate_elbo
            chain_dict[:s_entropy_chain][iter + 1] = s_entropy
            chain_dict[:wAlloc_elbo_chain][iter + 1] = wAlloc_elbo
            chain_dict[:γ_elbo_chain][iter + 1] = γ_elbo
            chain_dict[:α_elbo_chain][iter + 1] = α_elbo
        end
        iter = Int64(iter)
        # HDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(rhok_hat_vec, omegak_hat_vec,G,γ,α0,Tk)
        elbo_iter = data_elbo + assgn_entropy  + dHDP_surragate_elbo + s_entropy +  wAlloc_elbo + γ_elbo + α_elbo + pip_entropy
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
    
    elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_ = elbo_, rtik,pip_kj,c_ttprime_vec,θ_hat_vec, mk_hat_vec,λ0k_hat_vec,a0k_hat_vec,b0k_hat_vec,rhok_hat_vec, omegak_hat_vec,a_αt_hat_vec,b_αt_hat_vec,awt_hat_vec,bwt_hat_vec,a_γ_hat,b_γ_hat,e_γ,Tαk

    output_str_list = @name elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value;
    output_key_list = Symbol.(naming_vec(output_str_list));
    output_var_list = [elbo_, rtik_,pip_kj_,c_ttprime_vec_,θ_hat_vec_, mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_,rhok_hat_vec_, omegak_hat_vec_,a_αt_hat_vec_,b_αt_hat_vec_,awt_hat_vec_,bwt_hat_vec_,a_γ_hat_,b_γ_hat_,e_γ_,Tαk_,chain_dict,initDict,is_converged,truncation_value];

    

    outputs_dict = OrderedDict{Symbol,Any}();#Dict{Symbol,Union{Matrix{Union{Missing, Float64, Int64}},Matrix{Union{Float64, Int64}}}}()
    addToDict!(outputs_dict,output_key_list,output_var_list);

    return outputs_dict
end


function update_cluster_parameters!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K

    @inbounds for k in 1:K
        # clusters[k].pip_k .= 0.0
        clusters[k].Nkj.= 0.0
        clusters[k].x_hat .= 0.0
        clusters[k].x_hat_sq .= 0.0
        clusters[k].λ0k_hat.= 0.0
        clusters[k].a0k_hat.= 0.0
        clusters[k].mk_hat.= 0.0
        clusters[k].b0k_hat.= 0.0
        n_tik = 0.0
        @inbounds @fastmath for i in 1:N
            # _accumulate_deviance_log_norm_pdf3!(clusters[k],cellpop[i],dataparams)
            n_tik += cellpop[i].rtik[clusters[k].k]
            clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k] .* clusters[k].pip_k
            clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k] .* clusters[k].pip_k
        end
        clusters[k].Nk[1] = n_tik
        
        @inbounds @fastmath for j in 1:G
            # clusters[k].pip_k[j] = clusters[k].pip_k[j] / clusters[k].Nk[1]
            clusters[k].Nkj[j] +=   clusters[k].Nk[1] * clusters[k].pip_k[j]
            # if isnan(clusters[k].pip_k[j])
            #     clusters[k].pip_k[j] = 0.0
            # end
        end

        clusters[k].λ0k_hat .+=  modelparams.λ0_vec .+ clusters[k].Nkj
        clusters[k].a0k_hat .+=  modelparams.a0_vec .+ 1/2 .*  (clusters[k].Nkj .+1 )
        clusters[k].mk_hat .+=  (modelparams.λ0_vec .* modelparams.μ0_vec .+ clusters[k].x_hat) ./  clusters[k].λ0k_hat#(clusters[k].Nkj .+ modelparams.λ0_vec  )
        clusters[k].b0k_hat .+=  modelparams.b0_vec .+ 1/2 .* (clusters[k].x_hat_sq .+  (modelparams.μ0_vec .^2 .* modelparams.λ0_vec ) .- ( (clusters[k].x_hat .-  modelparams.λ0_vec .* modelparams.μ0_vec ) .^ 2 ./  clusters[k].λ0k_hat)) #(clusters[k].Nkj .+ modelparams.λ0_vec  )
        # if iszero(sum(clusters[k].pip_k))
        #     for j in 1:G
        #         clusters[k].pip_k[j] = 1.0
        #     end
        # end
        # normToProb3!(clusters[k].pip_k)

        clusters[k].pip_k .= 0.0

          @inbounds @fastmath for i in 1:N
              _accumulate_deviance_log_norm_pdf3!(clusters[k],cellpop[i],dataparams)
          end

          
          @inbounds @fastmath for j in 1:G
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

    end
    return clusters
end



