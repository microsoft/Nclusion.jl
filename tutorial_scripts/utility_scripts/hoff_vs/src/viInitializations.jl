# function init_θ_hat_tk(T,K,rho_hat_vec, omega_hat_vec;rand_init = false)
#     θ_hat_vec = nothing
#     if rand_init
#         θ_hat_vec = [rand(K+1) for t in 1:T]
#     else
#         θ_hat_vec = [βk_expected_value(rho_hat_vec, omega_hat_vec) for t in 1:T]
#     end
#     return θ_hat_vec
# end

function init_params_genes(G,λ0,μ0,a0,b0)
    if typeof(λ0) <: VecOrMat && eltype(λ0) <: Number
        λ0_vec = λ0
    else
        λ0_vec = λ0 .* ones(Float64,G)
    end
    if typeof(μ0) <: VecOrMat && eltype(μ0) <: Number
        μ0_vec = μ0
    else
        μ0_vec = μ0 .* ones(Float64,G)
    end
    if typeof(a0) <: VecOrMat && eltype(a0) <: Number
        a0_vec = a0
    else
        a0_vec = a0 .* ones(Float64,G)
    end
    if typeof(b0) <: VecOrMat && eltype(b0) <: Number
        b0_vec = b0
    else
        b0_vec = b0 .* ones(Float64,G)
    end
    
    
    return λ0_vec, μ0_vec, a0_vec, b0_vec
end
function init_params_states(K)
    rho_hat_vec = 0.25 .* ones(Float64,K)
    omega_hat_vec = 2 .* ones(Float64,K)
    return rho_hat_vec, omega_hat_vec
end
######################################################
function init_ηtkj_prior(x,K; pct_important=0.5)
    T = length(x)
    C_t = [length(el) for el in x]
    G = length(x[1][1])

    ηtkj_prior = [[[[pct_important, 1-pct_important] for j in 1:G] for k in 1:K] for t in 1:T] 
    return ηtkj_prior
end
######################################################
function init_mk_hat!(mk_hat_init,x,K,G;rand_init = false)
    μ0_vec = ones(G)
    if isnothing(mk_hat_init) && rand_init
        mk_hat_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_init) && !rand_init
        mk_hat_init = [μ0_vec for k in 1:K]
    end
    return mk_hat_init
end
function init_λ_sq_vec!(λ_sq_init,G;rand_init = false, lo=0,hi=1)
    λ_sq_vec = ones(G)
    if isnothing(λ_sq_init) && rand_init
        λ_sq_init = rand(Uniform(lo,hi),length(λ_sq_vec))
    elseif isnothing(λ_sq_init) && !rand_init
        λ_sq_init = λ_sq_vec
    end
    return λ_sq_init
end
function init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = false, lo=0,hi=1)
    σ_sq_k_vec = ones(G)
    if isnothing(σ_sq_k_init) && rand_init
        σ_sq_k_init = [rand(Uniform(lo,hi),length(σ_sq_k_vec)) for k in 1:K]
    elseif isnothing(σ_sq_k_init) && !rand_init
        σ_sq_k_init = [σ_sq_k_vec for k in 1:K] #
    end
    return σ_sq_k_init
end
function init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = false, lo=0,hi=1)
    v_sq_k_vec = ones(G)
    if isnothing(v_sq_k_hat_init) && rand_init
        v_sq_k_hat_init =  [rand(Uniform(lo,hi),length(v_sq_k_vec)) for k in 1:K]
    elseif isnothing(v_sq_k_hat_init) && !rand_init
        v_sq_k_hat_init =  [v_sq_k_vec for k in 1:K] #
    end 
    return v_sq_k_hat_init
end
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
function init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = false)
    if isnothing(c_ttprime_init) && rand_init
        c_ttprime_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_init) && !rand_init
        c_ttprime_init = [ones(T) ./T  for t in 1:T]
    end
    
    return c_ttprime_init
end
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
function init_d_hat_tk(T,g_hat_vec, h_hat_vec)
    d_hat_vec = [βk_expected_value(g_hat_vec, h_hat_vec) for t in 1:T]
    return d_hat_vec
end
function init_yjk_vec!(yjk_init,G,K;rand_init = false)
    if isnothing(yjk_init) && rand_init
        yjk_init = [[rand(Beta(1.,1.)) for j in 1:G] for k in 1:K] 
    elseif isnothing(yjk_init) && !rand_init
        yjk_init =[[0.5 for j in 1:G] for k in 1:K] 
    end

    return yjk_init
end
function init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    if isnothing(st_hat_init) && rand_init
        st_hat_init = [rand(Uniform(lo,hi)) for t in 1:T]
    elseif isnothing(st_hat_init) && !rand_init
        st_hat_init = [ϕ0 for t in 1:T]
    end
    return st_hat_init
end
function init_rtik_vec!(rtik_init,K,T,N_t;rand_init = false)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:N_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:N_t[t]] for t in 1:T]
    end
    

    return rtik_init
end
function init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = false)
    if isnothing(m_err_hat_init) && rand_init
        m_err_hat_init = rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_err_vec))#[μ0_vec for k in 1:K]
    elseif isnothing(m_err_hat_init) && !rand_init
        m_err_hat_init =μ0_err_vec
    end
    return m_err_hat_init
end
function init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = false, lo=0,hi=1)
    if isnothing(λ0_err_hat_init) && rand_init
        λ0_err_hat_init = rand(Uniform(lo,hi),length(λ0_err_vec)) #[λ0_vec for k in 1:K]; # 
    elseif isnothing(λ0_err_hat_init) && !rand_init
        λ0_err_hat_init = λ0_err_vec
    end
    return λ0_err_hat_init
end
function init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = false, lo=0,hi=1)
    if isnothing(a0_err_hat_init) && rand_init
        a0_err_hat_init = rand(Uniform(lo,hi),length(a0_err_vec))#[a0_vec for k in 1:K]; #
    elseif isnothing(a0_err_hat_init) && !rand_init
        a0_err_hat_init = a0_err_vec
    end
    return a0_err_hat_init
end
function init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = false, lo=0,hi=1)
    if isnothing(b0_err_hat_init) && rand_init
        b0_err_hat_init = rand(Uniform(lo,hi),length(b0_err_vec)) #[b0_vec for k in 1:K]; #
    elseif isnothing(b0_err_hat_init) && !rand_init
        b0_err_hat_init = b0_err_vec
    end
    return b0_err_hat_init
end

#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################
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

function low_memory_initialization(inputs_dict;mk_hat_init=nothing,v_sq_k_hat_init=nothing, λ_sq_init=nothing, σ_sq_k_init=nothing,st_hat_init=nothing, bwt_hat_init=nothing,a_αt_hat_init=nothing, b_αt_hat_init=nothing,a_γ_hat_init=nothing, b_γ_hat_init=nothing,d_hat_init=nothing,c_ttprime_init = nothing,rtik_init=nothing,yjk_init=nothing, gk_hat_init=nothing, hk_hat_init=nothing,uniform_theta_init = true, rand_init = false)
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

    input_str_list = @name float_type, cellpop, clusters,dataparams,modelparams,conditionparams,geneparams,ηk,Tk, elbolog, num_iter, num_local_iter;
    input_key_list = Symbol.(naming_vec(input_str_list));
    input_var_list = [float_type, cellpop, clusters,dataparams,modelparams,conditionparams, geneparams,ηk,Tk,elbolog, num_iter, num_local_iter];
    inputs_dict2 = OrderedDict()
    addToDict!(inputs_dict2,input_key_list,input_var_list);

    return inputs_dict2
end

#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################
function tidy_init_e_τ_μ_tikj_mat(xmat,K)
    T = length(unique(xmat[:,1]))
    N = length(unique(xmat[:,2]))
    G = length(unique(xmat[:,3]))
    timepoint_freq = countmap(Int.(xmat[1:G:end,1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    nrows = N*K*G
    ncols = 5
    gene_ids = collect(1:G)
    cell_ids = collect(1:N)
    timepoints = collect(1:T)
    state_ids= collect(1:K)
    e_τ_μ_tikj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    time = innermelt(timepoints, K .* G .* N_t)
    cells = innermelt(cell_ids,K*G)
    states = outermelt(innermelt(state_ids,G),N)
    genes = outermelt(gene_ids,N*K)
    e_τ_μ_tikj_mat[:,1] = time
    e_τ_μ_tikj_mat[:,2] = cells
    e_τ_μ_tikj_mat[:,3] = states
    e_τ_μ_tikj_mat[:,4] = genes
    return e_τ_μ_tikj_mat
end
function tidy_init_e_log_τ_kj_mat(λ0kmka0kb0k_hat_mat)
    K = length(unique(λ0kmka0kb0k_hat_mat[:,1]))
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]))
    state_ids = collect(1:K)
    gene_ids = collect(1:G)
    nrows = K*G
    ncols = 3
    states = innermelt(state_ids,G)
    genes = outermelt(gene_ids,K)
    e_log_τ_kj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_τ_kj_mat[:,1] = states
    e_log_τ_kj_mat[:,2] = genes
    return e_log_τ_kj_mat
end
function tidy_init_e_log_π_mat(θ_hat_mat)
    T = length(unique(θ_hat_mat[:,1]))
    Kplus =length(unique(θ_hat_mat[:,2]))
    timepoints = collect(1:T)
    state_ids = collect(1:Kplus)
    time = innermelt(timepoints,Kplus)
    states = outermelt(state_ids,T)
    nrows = Kplus*T
    ncols = 3 
    # params_vec = Vector{Vector{Union{Float64,Int}}}(undef,num_params)
    # for i in 1:num_params
    #     params_vec[i] = vcat(args[i]...)
    # end
    e_log_π_mat_init = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_π_mat_init[:,1] = time
    e_log_π_mat_init[:,2] = states
    return e_log_π_mat_init
end
function tidy_init_rtik_mat(e_τ_μ_tikj_mat)
    T = length(unique(e_τ_μ_tikj_mat[:,1]))
    N = length(unique(e_τ_μ_tikj_mat[:,2]))
    K = length(unique(e_τ_μ_tikj_mat[:,3]))
    G = length(unique(e_τ_μ_tikj_mat[:,4]))
    timepoints = collect(1:T)
    state_ids = collect(1:K)
    cell_ids = collect(1:N)
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
    # all([all(tidify_rtik(rtik)[:,col] .== rtik_mat[:,col]) for col in 1:size(rtik_mat)[2]])
    return rtik_mat
end
function tidy_init_e_log_τj_err_mat(λ0ma0b0_err_hat_mat)
    K = 1
    G = length(unique(λ0ma0b0_err_hat_mat[:,2]))
    state_ids = collect(1:K)
    gene_ids = collect(1:G)
    nrows = K*G
    ncols = 3
    states = innermelt(state_ids,G)
    genes = outermelt(gene_ids,K)
    e_log_τj_err_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_τj_err_mat[:,1] = states
    e_log_τj_err_mat[:,2] = genes
    return e_log_τj_err_mat
end
function tidy_init_e_τ_μ_tij_err_mat(xmat)
    T = length(unique(xmat[:,1]))
    N = length(unique(xmat[:,2]))
    G = length(unique(xmat[:,3]))
    K = 1
    timepoint_freq = countmap(Int.(xmat[1:G:end,1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    nrows = N*K*G
    ncols = 5
    gene_ids = collect(1:G)
    cell_ids = collect(1:N)
    timepoints = collect(1:T)
    state_ids= collect(1:K)
    e_τ_μ_tij_err_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    time = innermelt(timepoints, K .* G .* N_t)
    cells = innermelt(cell_ids,K*G)
    states = outermelt(innermelt(state_ids,G),N)
    genes = outermelt(gene_ids,N*K)
    e_τ_μ_tij_err_mat[:,1] = time
    e_τ_μ_tij_err_mat[:,2] = cells
    e_τ_μ_tij_err_mat[:,3] = states
    e_τ_μ_tij_err_mat[:,4] = genes
    return e_τ_μ_tij_err_mat
end



