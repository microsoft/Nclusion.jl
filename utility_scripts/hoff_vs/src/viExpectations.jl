# function log_π_expected_value(θ_hat;allK=false)
#     T = length(θ_hat)
#     K = Int64(length(θ_hat[1])-1)
#     e_log_π = Vector{Vector{Float64}}(undef,T)
#     for t in 1:T
#         if allK
#             digamma_sum = digamma(sum(θ_hat[t]))
#             e_log_π[t] = digamma.(θ_hat[t]) .- digamma_sum
#         else
#             digamma_sum = digamma(sum(θ_hat[t][1:K]))
#             e_log_π[t] = digamma.(θ_hat[t][1:K]) .- digamma_sum
#         end
        
#     end
    
#     # param_val_dict["log_π_expected_value"] = e_log_π
#     return e_log_π
# end
function τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
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
function log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec)
    K = length(a0k_hat_vec)
    e_log_τ_kj_vec = []
    for k in 1:K
        e_log_τ_kj = digamma.(a0k_hat_vec[k]) .- log.(b0k_hat_vec[k])
        # e_log_τ_k = sum(e_log_τ_kj)
        push!(e_log_τ_kj_vec,e_log_τ_kj)
    end
    return e_log_τ_kj_vec
end
function log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec)
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
function log_π_expected_value(θ_hat)
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
function βk_expected_value(rho_hat_vec, omega_hat_vec)
    K = length(rho_hat_vec)
    e_uk_vec = uk_expected_value(rho_hat_vec, omega_hat_vec)
    minus_e_uk_vec = 1. .- e_uk_vec
    cumprod_minus_e_uk_vec= cumprod(minus_e_uk_vec)
    app_e_uk_vec= deepcopy(e_uk_vec)
    app_cumprod_minus_e_uk_vec = deepcopy(cumprod_minus_e_uk_vec)
    append!(app_e_uk_vec,1.0) 
    insert!(app_cumprod_minus_e_uk_vec,1,1)
    e_βk_vec = Vector{Float64}(undef,K+1)
    for k in 1:K+1
        e_βk_vec[k] = app_e_uk_vec[k] * app_cumprod_minus_e_uk_vec[k]
    end
    return e_βk_vec
end
function uk_expected_value(rho_hat_vec, omega_hat_vec)
    minus_rho_hat_vec = 1 .- rho_hat_vec
    rho_omega_hat_vec = rho_hat_vec .* omega_hat_vec
    minus_rho_omega_hat_vec = minus_rho_hat_vec .*omega_hat_vec
    e_uk_vec = rho_omega_hat_vec ./ (rho_omega_hat_vec .+ minus_rho_omega_hat_vec)
    return e_uk_vec
end
function logUk_expected_value(rho_hat,omega_hat)
    return digamma.(rho_hat .* omega_hat) .- digamma.(omega_hat)
end
function log1minusUk_expected_value(rho_hat,omega_hat)
    return digamma.((1.0 .- rho_hat) .* omega_hat) .- digamma.(omega_hat)
end
function γ_expected_value(a_γ,b_γ)
    return a_γ ./ b_γ
end
function log_γ_expected_value(a_γ,b_γ)
    return digamma.(a_γ) .- log.(b_γ)
end
function αt_expected_value(a_α,b_α)
    return a_α ./ b_α
end
function log_αt_expected_value(a_α,b_α)
    return digamma.(a_α) .- log.(b_α)
end
function log_w_ttprime_expected_value(awt_hat,bwt_hat)
    e_log_tilde_wt_vec = log_tilde_wt_expected_value.(awt_hat,bwt_hat)
    e_log_minus_tilde_wt_vec = log_minus_tilde_wt_expected_value.(awt_hat,bwt_hat)
    T = length(e_log_tilde_wt_vec) + 1
    
    e_log_w_ttprime_vec = Vector{Vector{Float64}}(undef,T)
    e_log_w_ttprime_vec[1] = [0.0]
    for t in 2:T
        e_log_w_ttprime = Vector{Float64}(undef,t)
        for t_prime in 1:t
            if isone(t_prime)
                e_log_tilde_wt = 0.0
            else
                e_log_tilde_wt = e_log_tilde_wt_vec[t_prime-1]
            end
            if t_prime == t
                sum_e_log_minus_tilde_wt = 0.0
            else
                sum_e_log_minus_tilde_wt = sum(e_log_minus_tilde_wt_vec[t_prime:t-1])
            end
            e_log_w_ttprime[t_prime] = e_log_tilde_wt + sum_e_log_minus_tilde_wt 
        end
        e_log_w_ttprime_vec[t] = e_log_w_ttprime
    end
    return e_log_w_ttprime_vec
end
function log_tilde_wt_expected_value(awt_hat,bwt_hat)
    e_log_tilde_wt = digamma.(awt_hat) .- digamma.(awt_hat .+ bwt_hat)
    return e_log_tilde_wt
end
function log_minus_tilde_wt_expected_value(awt_hat,bwt_hat)
    e_log_minus_tilde_wt = digamma.(bwt_hat) .- digamma.(awt_hat .+ bwt_hat)
    return e_log_minus_tilde_wt
end

function log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
    e_log_τ_kj_vec = digamma.(a0_err_hat_vec) .- log.(b0_err_hat_vec)
    
    return e_log_τ_kj_vec
end
function log_τ_k_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
    e_log_τ_kj_vec = log_τ_kj_expected_value(a0_err_hat_vec, b0_err_hat_vec)
    e_log_τ_k = sum(e_log_τ_kj_vec)
    return e_log_τ_k
end
function τ_μ_error_expected_value(x,λ0_err_vec,m_err_vec,a0_err_vec, b0_err_vec)
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
            e_τ_μ_kjti3 = a0_err_vec ./  b0_err_vec .*  (x[t][i] .- m_err_vec) .^2 .+ 1 ./λ0_err_vec
            e_τ_μ_23 =  sum(e_τ_μ_kjti3)
            e_τ_μ_kjt3[i] = e_τ_μ_kjti3
            e_τ_μ_13[i] = e_τ_μ_23
        end
        e_τ_μ_kj_true3[t] = e_τ_μ_kjt3 
        e_τ_μ_true3[t] =e_τ_μ_13
    end

    return e_τ_μ_kj_true3,e_τ_μ_true3
end



#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################



function adjust_e_log_π_tk3!(t,conditionparams)
    Kplus = length(conditionparams[1].d_hat_t)
    # conditions_pis_sums = 0.0
    # e_log_π_t_cache = Vector{Float64}(undef,Kplus)
    for k in 1:Kplus
        conditions_pis_sums = 0.0
        for tt in 1:t
            conditions_pis_sums += conditionparams[t].c_tt_prime[tt] * (digamma(conditionparams[tt].d_hat_t[k]) - digamma(conditionparams[tt].d_hat_t_sum[1]))
        end
        conditionparams[t].e_log_π_t_cache[k] = conditions_pis_sums
    end
    return conditionparams
end
function log_π_expected_value_fast3!(conditionparams,dataparams,modelparams)
    T = dataparams.T
    Kplus = modelparams.K + 1
    for t in 1:T
        for k in 1:Kplus
            conditionparams[t].e_log_π_t_cache[k] = expectation_log_π_tk(conditionparams[t].d_hat_t[k], conditionparams[t].d_hat_t_sum[1] )
        end
    end
    return conditionparams
end
function expectation_βk(k,clusters,modelparams)
    K = modelparams.K
    Kplus = K + 1
    if k == Kplus
        e_uk = 1.0
    else
        e_uk = expectation_uk(clusters[k].gk_hat[1], clusters[k].hk_hat[1])
    end
    if isone(k)
        cumprod_minus_e_uk = 1.0
    else
        cumprod_minus_e_uk = recursive_minus_e_uk_cumprod(k-1,1.0,clusters)
    end
    # println("($e_uk,$cumprod_minus_e_uk)")
    e_βk = e_uk * cumprod_minus_e_uk
    return e_βk
end
function recursive_minus_e_uk_cumprod(k,cummulative_prod,clusters)
    if iszero(k)
        return cummulative_prod
    else
        cummulative_prod *= 1-expectation_uk(clusters[k].gk_hat[1], clusters[k].hk_hat[1])
        k -= 1
        recursive_minus_e_uk_cumprod(k,cummulative_prod,clusters)
    end
end
function expectation_log_normal_l_j!(cluster,cell, dataparams)
    G = dataparams.G
    for j in 1:G
        cell.cache[j]  =  -0.5 * (log(cluster.σ_sq_k_hat[j]) + 1/cluster.σ_sq_k_hat[j] * ( cell.xsq[j] - 2 * cell.x[j] * cluster.κk_hat[j] + cluster.yjk_hat[j] * (cluster.mk_hat[j] ^2 + cluster.v_sq_k_hat[j])  ))
    end
    return cluster
end
function expectation_log_π_tk(θ_hat_tk,θ_hat_t_sum)
    e_log_π_tk =  digamma(θ_hat_tk) - digamma(θ_hat_t_sum)
    return e_log_π_tk 
end
function expectation_log_tilde_wtt(awt_hat,bwt_hat)
    return  digamma(awt_hat) - digamma(awt_hat + bwt_hat)
end
function expectation_log_minus_tilde_wtt(awt_hat,bwt_hat)
    return  digamma(bwt_hat) - digamma(awt_hat + bwt_hat)
end
function expectation_uk(rho_hat, omega_hat)
    # minus_rho_hat = 1 - rho_hat
    # rho_omega_hat = rho_hat * omega_hat
    # minus_rho_omega_hat = (1 - rho_hat) * omega_hat
    e_uk = (rho_hat * omega_hat) / ((rho_hat * omega_hat)+ ((1 - rho_hat) * omega_hat))
    return e_uk
end
function expectation_αt(a_α,b_α)
    return a_α / b_α
end
function expectation_log_αt(a_α,b_α)
    return digamma(a_α) - log(b_α)
end
function expectation_logUk(rho_hat,omega_hat)
    return digamma(rho_hat * omega_hat) - digamma(omega_hat)
end
function expectation_log1minusUk(rho_hat,omega_hat)
    return digamma((1.0 - rho_hat) * omega_hat) - digamma(omega_hat)
end

#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################

function tidy_τ_μ_expected_value(xmat,λ0kmka0kb0k_hat_mat; e_τ_μ_tikj_mat_init = nothing)
    T = length(unique(xmat[:,1]))
    N = length(unique(xmat[:,2]))
    G = length(unique(xmat[:,3]))
    K = length(unique(λ0kmka0kb0k_hat_mat[:,1])) 
    gene_ids = collect(1:G)
    cell_ids = collect(1:N)
    timepoints = collect(1:T)
    state_ids= collect(1:K)
    if isnothing(e_τ_μ_tikj_mat_init)
        timepoint_freq = countmap(Int.(xmat[1:G:end,1]))
        N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
        nrows = N*K*G
        ncols = 5
        e_τ_μ_tikj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* G .* N_t)
        cells = innermelt(cell_ids,K*G)
        states = outermelt(innermelt(state_ids,G),N)
        genes = outermelt(gene_ids,N*K)
        e_τ_μ_tikj_mat[:,1] = time
        e_τ_μ_tikj_mat[:,2] = cells
        e_τ_μ_tikj_mat[:,3] = states
        e_τ_μ_tikj_mat[:,4] = genes
    else
        e_τ_μ_tikj_mat = e_τ_μ_tikj_mat_init
    end

    
    counter = 0
    for i in cell_ids
        for k in state_ids
            counter += 1 
            param_indx = (k-1)*G + 1:(k)*G
            data_indx = (i-1)*G + 1:(i)*G
            indx = (counter-1)*G + 1:(counter)*G
            @views λ0k =  λ0kmka0kb0k_hat_mat[param_indx , 3]
            @views mk =   λ0kmka0kb0k_hat_mat[param_indx , 4]
            @views a0k =  λ0kmka0kb0k_hat_mat[param_indx , 5]
            @views b0k =  λ0kmka0kb0k_hat_mat[param_indx , 6]
            @views x = xmat[data_indx , 4]
            @views e_τ_μ_tikj_mat[indx,end] = a0k ./  b0k .*  (x .- mk) .^2 .+ 1 ./λ0k
        end
    end
    return e_τ_μ_tikj_mat
end
function tidy_log_τ_kj_expected_value(λ0kmka0kb0k_hat_mat; e_log_τ_kj_mat_init = nothing )
    K = length(unique(λ0kmka0kb0k_hat_mat[:,1]))
    G = length(unique(λ0kmka0kb0k_hat_mat[:,2]))
    state_ids = collect(1:K)
    gene_ids = collect(1:G)
    nrows = K*G
    ncols = 3
    if isnothing(e_log_τ_kj_mat_init)
        states = innermelt(state_ids,G)
        genes = outermelt(gene_ids,K)
        e_log_τ_kj_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        e_log_τ_kj_mat[:,1] = states
        e_log_τ_kj_mat[:,2] = genes
    else
        e_log_τ_kj_mat = e_log_τ_kj_mat_init
    end

    for k in state_ids
        param_indx =  (k-1)*G + 1:(k)*G
        @views a0k =  λ0kmka0kb0k_hat_mat[param_indx , 5]
        @views b0k =  λ0kmka0kb0k_hat_mat[param_indx , 6]
        @views e_log_τ_kj_mat[param_indx,end] = digamma.(a0k) .- log.(b0k)
    end
    return e_log_τ_kj_mat
end
function tidy_log_τj_err_expected_value(λ0ma0b0_err_hat_mat; e_log_τj_err_mat_init = nothing )
    K = 1
    G = length(unique(λ0ma0b0_err_hat_mat[:,1]))
    state_ids = collect(1:K)
    gene_ids = collect(1:G)
    nrows = K*G
    ncols = 3
    if isnothing(e_log_τj_err_mat_init)
        states = innermelt(state_ids,G)
        genes = outermelt(gene_ids,K)
        e_log_τj_err_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        e_log_τj_err_mat[:,1] = states
        e_log_τj_err_mat[:,2] = genes
    else
        e_log_τj_err_mat = e_log_τj_err_mat_init
    end

    for k in state_ids
        param_indx =  (k-1)*G + 1:(k)*G
        @views a0k =  λ0ma0b0_err_hat_mat[param_indx , 5]
        @views b0k =  λ0ma0b0_err_hat_mat[param_indx , 6]
        @views e_log_τj_err_mat[param_indx,end] = digamma.(a0k) .- log.(b0k)
    end
    return e_log_τj_err_mat
end
function tidy_τ_μ_err_expected_value(xmat,λ0ma0b0_err_hat_mat; e_τ_μ_tij_err_mat_init = nothing)
    T = length(unique(xmat[:,1]))
    N = length(unique(xmat[:,2]))
    G = length(unique(xmat[:,3]))
    K = 1 
    gene_ids = collect(1:G)
    cell_ids = collect(1:N)
    timepoints = collect(1:T)
    state_ids= collect(1:K)
    if isnothing(e_τ_μ_tij_err_mat_init)
        timepoint_freq = countmap(Int.(xmat[1:G:end,1]))
        N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
        nrows = N*K*G
        ncols = 5
        e_τ_μ_tij_err_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        time = innermelt(timepoints, K .* G .* N_t)
        cells = innermelt(cell_ids,K*G)
        states = outermelt(innermelt(state_ids,G),N)
        genes = outermelt(gene_ids,N*K)
        e_τ_μ_tij_err_mat[:,1] = time
        e_τ_μ_tij_err_mat[:,2] = cells
        e_τ_μ_tij_err_mat[:,3] = states
        e_τ_μ_tij_err_mat[:,4] = genes
    else
        e_τ_μ_tij_err_mat = e_τ_μ_tij_err_mat_init
    end

    
    counter = 0
    for i in cell_ids
        for k in state_ids
            counter += 1 
            param_indx = (k-1)*G + 1:(k)*G
            data_indx = (i-1)*G + 1:(i)*G
            indx = (counter-1)*G + 1:(counter)*G
            @views λ0k =  λ0ma0b0_err_hat_mat[param_indx , 3]
            @views mk =   λ0ma0b0_err_hat_mat[param_indx , 4]
            @views a0k =  λ0ma0b0_err_hat_mat[param_indx , 5]
            @views b0k =  λ0ma0b0_err_hat_mat[param_indx , 6]
            @views x = xmat[data_indx , 4]
            @views e_τ_μ_tij_err_mat[indx,end] = a0k ./  b0k .*  (x .- mk) .^2 .+ 1 ./λ0k
        end
    end
    return e_τ_μ_tij_err_mat
end
function tidy_log_π_expected_value(θ_hat_mat; e_log_π_mat_init=nothing)
    T = length(unique(θ_hat_mat[:,1]))
    Kplus =length(unique(θ_hat_mat[:,2]))
    timepoints = collect(1:T)
    state_ids = collect(1:Kplus)
    if isnothing(e_log_π_mat_init)
        time = innermelt(timepoints,Kplus)
        states = outermelt(state_ids,T)
        nrows = Kplus*T
        ncols = 3 
        e_log_π_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
        e_log_π_mat[:,1] = time
        e_log_π_mat[:,2] = states
    else
        e_log_π_mat = e_log_π_mat_init
    end

    for t in timepoints
        time_indx = (t-1)*Kplus + 1:(t)*Kplus
        digamma_sum = digamma(sum(θ_hat_mat[time_indx,3]))
        e_log_π_mat[time_indx,end] = digamma.(θ_hat_mat[time_indx,3]) .- digamma_sum
    end
    return e_log_π_mat
end
function tidy_βk_expected_value(ρkωkckdk_hat_mat)
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    Kplus = K+1
    nrows = Kplus
    ncols = 2
    state_ids = collect(1:Kplus)
    

    e_uk_mat = tidy_uk_expected_value(ρkωkckdk_hat_mat)
    e_uk_vec = e_uk_mat[:,2]
    minus_e_uk_vec = 1. .- e_uk_vec
    cumprod_minus_e_uk_vec= cumprod(minus_e_uk_vec)
    app_e_uk_vec= deepcopy(e_uk_vec)
    app_cumprod_minus_e_uk_vec = deepcopy(cumprod_minus_e_uk_vec)
    append!(app_e_uk_vec,1.0) 
    insert!(app_cumprod_minus_e_uk_vec,1,1)
    e_βk_vec = Vector{Float64}(undef,Kplus)
    for k in 1:Kplus
        e_βk_vec[k] = app_e_uk_vec[k] * app_cumprod_minus_e_uk_vec[k]
    end
    e_βk_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    @views e_βk_mat[:,1] = state_ids
    @views e_βk_mat[:,2] = e_βk_vec
    return e_βk_mat
end
function tidy_uk_expected_value(ρkωkckdk_hat_mat)
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    state_ids = collect(1:K)
    @views rho_hat_vec = ρkωkckdk_hat_mat[:,2];
    @views omega_hat_vec = ρkωkckdk_hat_mat[:,3];
    minus_rho_hat_vec = 1 .- rho_hat_vec
    rho_omega_hat_vec = rho_hat_vec .* omega_hat_vec
    minus_rho_omega_hat_vec = minus_rho_hat_vec .*omega_hat_vec
    e_uk_vec = rho_omega_hat_vec ./ (rho_omega_hat_vec .+ minus_rho_omega_hat_vec)
    nrows = K
    ncols = 2
    e_uk_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    @views e_uk_mat[:,1] = state_ids
    @views e_uk_mat[:,2] = e_uk_vec

    return e_uk_mat
end
function tidy_logUk_expected_value(ρkωkckdk_hat_mat)
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    state_ids = collect(1:K)
    @views rho_hat_vec = ρkωkckdk_hat_mat[:,2];
    @views omega_hat_vec = ρkωkckdk_hat_mat[:,3];
    nrows = K
    ncols = 2
    e_log_uk_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    @views e_log_uk_mat[:,1] = state_ids
    @views e_log_uk_mat[:,2] = digamma.(rho_hat_vec .* omega_hat_vec) .- digamma.(omega_hat_vec)
    return e_log_uk_mat
end
function tidy_log1minusUk_expected_value(ρkωkckdk_hat_mat)
    K = length(unique(ρkωkckdk_hat_mat[:,1]))
    state_ids = collect(1:K)
    @views rho_hat_vec = ρkωkckdk_hat_mat[:,2];
    @views omega_hat_vec = ρkωkckdk_hat_mat[:,3];
    nrows = K
    ncols = 2
    e_log_log1minusuk_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    @views e_log_log1minusuk_mat[:,1] = state_ids
    @views e_log_log1minusuk_mat[:,2] = digamma.((1.0 .- rho_hat_vec) .* omega_hat_vec) .- digamma.(omega_hat_vec)
    return e_log_log1minusuk_mat
end
function tidy_γ_expected_value(gamma_mat)
    a_γ = gamma_mat[1,2]
    b_γ = gamma_mat[1,3]
    e_γ_mat = Matrix{Union{Float64,Int}}(undef,1,2)
    e_γ_mat[1,1] = 1
    e_γ_mat[1,2] = a_γ ./ b_γ
    return e_γ_mat
end
function tidy_log_γ_expected_value(gamma_mat)
    a_γ = gamma_mat[1,2]
    b_γ = gamma_mat[1,3]
    e_log_γ_mat = Matrix{Union{Float64,Int}}(undef,1,2)
    e_log_γ_mat[1,1] = 1
    e_log_γ_mat[1,2] = digamma.(a_γ) .- log.(b_γ)
    return e_log_γ_mat
end
function tidy_αt_expected_value(aαtbαtawtbwt_hat_mat)
    T =  length(unique(aαtbαtawtbwt_hat_mat[:,1]))
    a_αt_hat_vec = aαtbαtawtbwt_hat_mat[:,2]
    b_αt_hat_vec = aαtbαtawtbwt_hat_mat[:,3]
    timepoints = collect(1:T)
    nrows = T
    ncols = 2
    e_αt_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_αt_mat[:,1] = timepoints
    e_αt_mat[:,2] = a_αt_hat_vec ./  b_αt_hat_vec
    return e_αt_mat
end
function tidy_log_αt_expected_value(aαtbαtawtbwt_hat_mat)
    T =  length(unique(aαtbαtawtbwt_hat_mat[:,1]))
    a_αt_hat_vec = aαtbαtawtbwt_hat_mat[:,2]
    b_αt_hat_vec = aαtbαtawtbwt_hat_mat[:,3]
    timepoints = collect(1:T)
    nrows = T
    ncols = 2
    e_log_αt_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_αt_mat[:,1] = timepoints
    e_log_αt_mat[:,2] = digamma.(a_αt_hat_vec) .- log.(b_αt_hat_vec)  
    return e_log_αt_mat
end
function tidy_log_tilde_wt_expected_value(aαtbαtawtbwt_hat_mat)
    T =  length(unique(aαtbαtawtbwt_hat_mat[:,1]))
    awt_hat_vec = aαtbαtawtbwt_hat_mat[:,4]
    bwt_hat_vec = aαtbαtawtbwt_hat_mat[:,5]
    timepoints = collect(1:T)
    nrows = T
    ncols = 2
    e_log_tilde_wt_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_tilde_wt_mat[:,1] = timepoints
    e_log_tilde_wt_mat[:,2] = digamma.(awt_hat_vec) .- digamma.(awt_hat_vec .+ bwt_hat_vec) 
    return e_log_tilde_wt_mat 
end
function tidy_log_minus_tilde_wt_expected_value(aαtbαtawtbwt_hat_mat)
    T =  length(unique(aαtbαtawtbwt_hat_mat[:,1]))
    awt_hat_vec = aαtbαtawtbwt_hat_mat[:,4]
    bwt_hat_vec = aαtbαtawtbwt_hat_mat[:,5]
    timepoints = collect(1:T)
    nrows = T
    ncols = 2
    e_log_minus_tilde_wt_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_minus_tilde_wt_mat[:,1] = timepoints
    e_log_minus_tilde_wt_mat[:,2] = digamma.(bwt_hat_vec) .- digamma.(awt_hat_vec .+ bwt_hat_vec) 
    return e_log_minus_tilde_wt_mat 
end
function tidy_log_w_ttprime_expected_value(aαtbαtawtbwt_hat_mat)
    
    T = length(unique(aαtbαtawtbwt_hat_mat[:,1])) + 1
    timepoints = collect(1:T)
    progressing_timepoints = [collect(1:t) for t in 1:T]
    time = innermelt(timepoints,[length(el) for el in progressing_timepoints])
    progressing_time = recursive_flatten(progressing_timepoints)
    nrows = length(progressing_time)
    ncols = 3

    e_log_tilde_wt_mat = tidy_log_tilde_wt_expected_value(aαtbαtawtbwt_hat_mat)
    e_log_minus_tilde_wt_mat = tidy_log_minus_tilde_wt_expected_value(aαtbαtawtbwt_hat_mat)
    e_log_tilde_wt_vec = e_log_tilde_wt_mat[:,end]
    e_log_minus_tilde_wt_vec =e_log_minus_tilde_wt_mat[:,end]
    e_log_w_ttprime_vec = Vector{Vector{Float64}}(undef,T)
    e_log_w_ttprime_vec[1] = [0.0]
    for t in 2:T
        e_log_w_ttprime = Vector{Float64}(undef,t)
        for t_prime in 1:t
            if isone(t_prime)
                e_log_tilde_wt = 0.0
            else
                e_log_tilde_wt = e_log_tilde_wt_vec[t_prime-1]
            end
            if t_prime == t
                sum_e_log_minus_tilde_wt = 0.0
            else
                sum_e_log_minus_tilde_wt = sum(e_log_minus_tilde_wt_vec[t_prime:t-1])
            end
            e_log_w_ttprime[t_prime] = e_log_tilde_wt + sum_e_log_minus_tilde_wt 
        end
        e_log_w_ttprime_vec[t] = e_log_w_ttprime
    end
    e_log_w_ttprime_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    e_log_w_ttprime_mat[:,1] = time
    e_log_w_ttprime_mat[:,2] = progressing_time 
    e_log_w_ttprime_mat[:,3] = recursive_flatten(e_log_w_ttprime_vec) 
    return e_log_w_ttprime_mat
end

#####################################################
############# DEPRACATED FUNCTIONS ##################
#####################################################
function depracated_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    T = length(x)
    C_t = [length(el) for el in x]
    K = length(a0k_hat_vec)
    G = length(a0k_hat_vec[1])
    e_τ_μ_tikj = [[[[(a0k_hat_vec[k][j]/b0k_hat_vec[k][j])*(x[t][i][j] - mk_hat_vec[k][j])^2 + (λ0k_hat_vec[k][j])^(-1.0)  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return e_τ_μ_tikj
end
function depracated_log_τ_expected_value(a_0kj_hat, b_0kj_hat)
    e_log_τkj = digamma.(a_0kj_hat) .- log.(b_0kj_hat) 
    return e_log_τkj
end
#####################################################


