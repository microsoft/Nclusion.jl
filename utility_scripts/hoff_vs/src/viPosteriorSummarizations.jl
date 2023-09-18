function sample_z_post(rtik)
    T = length(rtik)
    z_post = [rand(arraydist(Distributions.Categorical.(rtik[t]))) for t in 1:T]
    return z_post
end

function vi_make_z_post_s(rtik; S=100)
    z_post_s = [sample_z_post(rtik) for s in 1:S]
    return z_post_s
end
function sample_st_post(c_ttprime)
    # T = length(rtik)
    st_post = rand(arraydist(Distributions.Categorical.(c_ttprime)))
    return st_post
end

function vi_make_st_post_s(c_ttprime; S=100)
    st_post_s = [sample_st_post(c_ttprime) for s in 1:S]
    return st_post_s
end

function sample_π_post(θ_hat)
    # T = length(rtik)
    π_post = permutedims(rand(arraydist(Dirichlet.(θ_hat))))
    return π_post
end

function vi_make_π_post_s(θ_hat; S=100)
    st_post_s = [sample_π_post(θ_hat) for s in 1:S]
    return st_post_s
end
function calc_π_post_mean(π_post_s)
    S = length(π_post_s)
    KMax_and1 = size(π_post_s[1])[2]
    T = size(π_post_s[1])[1]
    π_post_s_mat = Array{Float64}(undef,T,KMax_and1,S);
    for s in 1:S
        π_post_s_mat[:,:,s] = π_post_s[s]
    end
    π_post_mean = mean(π_post_s_mat,dims=3)[:,:,1]
    return π_post_mean
end


function sample_normalgamma_μ_τ_post(mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_;Ns= 1000)
    KMax =  length(b0k_hat_vec_)
    τ_post_s = [rand(arraydist(Gamma.(a0k_hat_vec_[k], 1 ./ b0k_hat_vec_[k])),Ns) for k in 1: KMax];
    λτ_inv_post_s = [(λ0k_hat_vec_[k] .*  τ_post_s[k]) .^ (-1) for k in 1: KMax];
    μ_post_s = [rand(arraydist(Normal.(mk_hat_vec_[k], sqrt.(λτ_inv_post_s[k])))) for k in 1:KMax];
    return τ_post_s,μ_post_s
end
function calc_normalgamma_μ_τ_post_mean(mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_;Ns= 1000)
    τ_post_s,μ_post_s = sample_normalgamma_μ_τ_post(mk_hat_vec_,λ0k_hat_vec_,a0k_hat_vec_,b0k_hat_vec_;Ns= Ns);
    mean_τ_post = vec.(mean.(τ_post_s,dims=2));
    mean_μ_post = vec.(mean.(μ_post_s,dims=2));
    return mean_τ_post,mean_μ_post
end
function calc_normalgamma_μ_τ_post_mean(τ_post_s,μ_post_s)
    mean_τ_post = vec.(mean.(τ_post_s,dims=2));
    mean_μ_post = vec.(mean.(μ_post_s,dims=2));
    return mean_τ_post,mean_μ_post
end

function sample_gamma_τ_post(a0k_hat_vec_,b0k_hat_vec_;Ns= 1000)
    KMax =  length(b0k_hat_vec_)
    τ_post_s = [rand(arraydist(Gamma.(a0k_hat_vec_[k], 1 ./ b0k_hat_vec_[k])),Ns) for k in 1: KMax];
    
    return τ_post_s
end
function calc_gamma_τ_post_mean(a0k_hat_vec_,b0k_hat_vec_;Ns= 1000)
    τ_post_s= sample_gamma_τ_post(a0k_hat_vec_,b0k_hat_vec_;Ns= Ns);
    mean_τ_post = vec.(mean.(τ_post_s,dims=2));

    return mean_τ_post
end
function calc_gamma_τ_post_mean(τ_post_s)
    mean_τ_post = vec.(mean.(τ_post_s,dims=2));
    return mean_τ_post
end



######## same as get_average_posterior_cluster_frequency2 in turingChainProcessing! Do not export until you do multiple dispatch!
function get_average_posterior_cluster_frequency(z_post_s,T,true_z,KMax,KTrue,C_t)
    

    function genData_vec2dict(time_vec)
        T = length(time_vec)
        C_t = [length(c) for c in time_vec]
        final_dict = Dict()
        for t in 1:T
            for c in 1:C_t[t]
                final_dict[(c,t)] = time_vec[t][c]
            end
        end
        # final_dict = Dict((c,t) => time_vec[t][c] for c in 1:C_t[t] for t in 1:T)
        return final_dict
    end

    S = length(z_post_s)
    post_z_dict_s = genData_vec2dict.(z_post_s)
    true_z_dict = genData_vec2dict(true_z)
    counts_mat_s = zeros(Int,KMax,KTrue,T,S)
    for s in 1:S
        for t in 1:T
            for c in 1:C_t[t]
                row = Int(post_z_dict_s[s][c,t])
                col = Int(true_z_dict[c,t])
                counts_mat_s[row,col,t,s] += 1
            end
        end
    end
    avg_counts_mat = mean(counts_mat_s,dims=4)
    return avg_counts_mat
end

#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################

function tidy_vi_make_z_post_s(outputs_dict; S=100)
    rtik_mat = outputs_dict[:rtik_mat_]
    T = length(unique(rtik_mat[:,1]));
    N = length(unique(rtik_mat[:,2]));
    K = length(unique(rtik_mat[:,3]))
    N_t = tidy_get_Nt_from_rtikmat(rtik_mat)
    timepoints = collect(1:T)
    cell_ids = collect(1:N)
    sample_ids = collect(1:S)
    nrows = N * S
    ncols = size(rtik_mat)[2]

    time = innermelt(timepoints, N_t)
    

    rtik_ = [rtik_mat[(i-1)*K+1:i*K,end] for i in 1:N]
    z_samples = [rand(arraydist(Distributions.Categorical.(rtik_))) for s in 1:S]
    flatten_z = recursive_flatten(z_samples)
    z_samples_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    z_samples_mat[:,1] = innermelt(sample_ids,N)
    z_samples_mat[:,2] = outermelt(time,S)
    z_samples_mat[:,3] = outermelt(cell_ids,S)
    z_samples_mat[:,4] = flatten_z

    return z_samples_mat
end
function tidy_vi_make_s_post_s(outputs_dict; S=100)
    c_ttprime_mat = outputs_dict[:c_ttprime_mat_]
    T = length(unique(c_ttprime_mat[:,1]));
    timepoints = collect(1:T)
    sample_ids = collect(1:S)
    nrows = T * S
    ncols = size(c_ttprime_mat)[2]
    time = timepoints
    c_ttprime_ = [c_ttprime_mat[(t-1)*T+1:t*T,end] for t in 1:T]
    s_samples = [rand(arraydist(Distributions.Categorical.(c_ttprime_))) for s in 1:S]
    flatten_s = recursive_flatten(s_samples)
    s_samples_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    s_samples_mat[:,1] = innermelt(sample_ids,T)
    s_samples_mat[:,2] = outermelt(time,S)
    s_samples_mat[:,3] = flatten_s


    
    return s_samples_mat
end
function tidy_vi_make_π_post_s(outputs_dict; S=100)
    θ_hat_mat = outputs_dict[:θ_hat_mat_]
    T = length(unique(θ_hat_mat[:,1]));
    Kplus = length(unique(θ_hat_mat[:,2]));
    timepoints = collect(1:T)
    kplus_id = collect(1:Kplus)
    sample_ids = collect(1:S)
    nrows =  Kplus*T * S
    ncols = size(θ_hat_mat)[2] +1 
    time = innermelt(timepoints,Kplus)
    θ_hat_ = [θ_hat_mat[(t-1)*Kplus+1:t* Kplus,end] for t in 1:T]
    π_samples = [rand.(Distributions.Dirichlet.(θ_hat_)) for s in 1:S]
    flatten_π_samples = recursive_flatten(π_samples)
    π_samples_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    π_samples_mat[:,1] = innermelt(sample_ids,Kplus*T)
    π_samples_mat[:,2] = outermelt(time,S)
    π_samples_mat[:,3] = outermelt(kplus_id,T*S)
    π_samples_mat[:,4] = flatten_π_samples

    
    return π_samples_mat
end
function tidy_calc_π_post_summarization(π_samples_mat;quantile_lo=0.05,quantile_hi=0.95)
    S = length(unique(π_samples_mat[:,1]))
    T = length(unique(π_samples_mat[:,2]))
    Kplus = length(unique(π_samples_mat[:,3]))

    π_samples_values = [π_samples_mat[v:T*Kplus:end,end] for v in 1:T*Kplus]
    π_samples_means = mean.(π_samples_values)
    π_samples_std = std.(π_samples_values)
    π_samples_quantiles = [quantile(v, [quantile_lo, quantile_hi]) for v in π_samples_values]
    π_samples_quantiles_mat =permutedims(reduce(hcat,π_samples_quantiles))
    π_samples_lower_quantile = [q[1] - m for (q, m) in zip(π_samples_quantiles, π_samples_means)]
    π_samples_upper_quantile = [q[2] - m for (q, m) in zip(π_samples_quantiles, π_samples_means)]
    
    π_post_summary = reduce(hcat, [π_samples_mat[1:T*Kplus,1:3],π_samples_means,π_samples_std,π_samples_quantiles_mat,π_samples_lower_quantile,π_samples_upper_quantile]) 
    
    return π_post_summary
end
function tidy_vi_make_normalgamma_μ_τ_post_s(outputs_dict; S=100)
    λ0kmka0kb0k_hat_mat = outputs_dict[:λ0kmka0kb0k_hat_mat_]
    K = length(unique(λ0kmka0kb0k_hat_mat[:,1]));
    G= length(unique(λ0kmka0kb0k_hat_mat[:,2]));
    gene_ids = collect(1:G)
    state_id = collect(1:K)
    sample_ids = collect(1:S)
    nrows =  G*K* S
    ncols = 5 
    states = λ0kmka0kb0k_hat_mat[:,1]
    genes = λ0kmka0kb0k_hat_mat[:,2]
    λ0k_hat_ = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k* G,3] for k in 1:K]
    mk_hat_ = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k* G,4] for k in 1:K]
    a0k_hat_ = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k* G,5] for k in 1:K]
    b0k_hat_ = [λ0kmka0kb0k_hat_mat[(k-1)*G+1:k* G,6] for k in 1:K]

    τ_samples = [[rand(arraydist(Gamma.(a0k_hat_[k], 1 ./ b0k_hat_[k]))) for k in 1:K] for s in 1:S];
    λτ_inv_post_s =[ [(λ0k_hat_[k] .*  τ_samples[s][k]) .^ (-1) for k in 1: K] for s in 1:S];
    μ_samples = [[rand(arraydist(Normal.(mk_hat_[k], sqrt.(λτ_inv_post_s[s][k])))) for k in 1:K] for s in 1:S];
    flatten_τ_samples = recursive_flatten(τ_samples)
    flatten_μ_samples = recursive_flatten(μ_samples)
    μ_τ_samples_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    μ_τ_samples_mat[:,1] = innermelt(sample_ids,G*K)
    μ_τ_samples_mat[:,2] = outermelt(states,S)
    μ_τ_samples_mat[:,3] = outermelt(genes,S)
    μ_τ_samples_mat[:,4] = flatten_μ_samples
    μ_τ_samples_mat[:,5] = flatten_τ_samples

    
    return μ_τ_samples_mat
end
function tidy_calc_μ_τ_post_summarization(μ_τ_samples_mat;quantile_lo=0.05,quantile_hi=0.95)
    S = length(unique(μ_τ_samples_mat[:,1]))
    K = length(unique(μ_τ_samples_mat[:,2]))
    G = length(unique(μ_τ_samples_mat[:,3]))

    τ_samples_values = [μ_τ_samples_mat[v:K*G:end,5] for v in 1:K*G]
    τ_samples_means = mean.(τ_samples_values)
    τ_samples_std = std.(τ_samples_values)
    τ_samples_quantiles = [quantile(v, [quantile_lo, quantile_hi]) for v in τ_samples_values]
    τ_samples_quantiles_mat =permutedims(reduce(hcat,τ_samples_quantiles))
    τ_samples_lower_quantile = [q[1] - m for (q, m) in zip(τ_samples_quantiles, τ_samples_means)]
    τ_samples_upper_quantile = [q[2] - m for (q, m) in zip(τ_samples_quantiles, τ_samples_means)]
    
    τ_post_summary = reduce(hcat, [μ_τ_samples_mat[1:K*G,1:3],τ_samples_means,τ_samples_std,τ_samples_quantiles_mat,τ_samples_lower_quantile,τ_samples_upper_quantile])


    μ_samples_values = [μ_τ_samples_mat[v:K*G:end,4] for v in 1:K*G]
    μ_samples_means = mean.(μ_samples_values)
    μ_samples_std = std.(μ_samples_values)
    μ_samples_quantiles = [quantile(v, [quantile_lo, quantile_hi]) for v in μ_samples_values]
    μ_samples_quantiles_mat =permutedims(reduce(hcat,μ_samples_quantiles))
    μ_samples_lower_quantile = [q[1] - m for (q, m) in zip(μ_samples_quantiles, μ_samples_means)]
    μ_samples_upper_quantile = [q[2] - m for (q, m) in zip(μ_samples_quantiles, μ_samples_means)]
    
    μ_post_summary = reduce(hcat, [μ_τ_samples_mat[1:K*G,1:3],μ_samples_means,μ_samples_std,μ_samples_quantiles_mat,μ_samples_lower_quantile,μ_samples_upper_quantile]) 
    
    return μ_post_summary,τ_post_summary
end
function tidy_get_average_posterior_cluster_frequency(zmat,z_samples_mat,KTrue,KMax)
    S = length(unique(z_samples_mat[:,1]))
    T = length(unique(z_samples_mat[:,2]))
    N = length(unique(z_samples_mat[:,3]))
    sample_ids = collect(1:S)
    timepoints = collect(1:T)
    # KMax = length(unique(z_samples_mat[:,end]))
    # KTrue = length(unique(zmat[:,end]))
    timepoint_freq = countmap(Int.(zmat[:,1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    true_z = [Int.(zmat[st:en,end]) for (st,en) in timeranges]
    z_post_s = [[Int.(z_samples_mat[(s-1)*N+st:(s-1)*N+en,end]) for (st,en) in timeranges] for s in 1:S]
    # z_post_s,T,true_z,KMax,KTrue,C_t
    avg_counts_mat = get_average_posterior_cluster_frequency(z_post_s,T,true_z,KMax,KTrue,N_t)

    nrows = T*KMax*KTrue
    ncols = 4
    tidy_avg_counts_mat = Matrix{Union{Float64,Int}}(undef,nrows,ncols)
    avg_counts_vec = [vec(reshape(avg_counts_mat[:,:,t,1], KMax*KTrue,1)) for t in 1:T]
    tidy_avg_counts_mat[:,1] = innermelt(timepoints,KMax*KTrue)
    tidy_avg_counts_mat[:,2] = outermelt(innermelt(collect(1:KTrue),KMax), T)
    tidy_avg_counts_mat[:,3] = outermelt(collect(1:KMax),T*KTrue)
    tidy_avg_counts_mat[:,4] = recursive_flatten(avg_counts_vec)

    # S = length(z_post_s)
    # post_z_dict_s = genData_vec2dict.(z_post_s)
    # true_z_dict = genData_vec2dict(true_z)
    # counts_mat_s = zeros(Int,KMax,KTrue,T,S)
    # for s in 1:S
    #     for t in 1:T
    #         for c in 1:C_t[t]
    #             row = Int(post_z_dict_s[s][c,t])
    #             col = Int(true_z_dict[c,t])
    #             counts_mat_s[row,col,t,s] += 1
    #         end
    #     end
    # end
    # avg_counts_mat = mean(counts_mat_s,dims=4)
    return tidy_avg_counts_mat, avg_counts_mat
end





