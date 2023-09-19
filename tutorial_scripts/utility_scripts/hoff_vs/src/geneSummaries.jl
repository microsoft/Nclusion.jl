function calc_cell_normal_μ_τ_ll_scores(x,mean_τ_post,mean_μ_post)
    T = length(x);
    C_t = [length(el) for el in x];
    KMax = length(mean_μ_post)
    cell_ll_scores =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T);
    τk = [mean_τ_post[k] for k in 1:KMax];
    μk = [mean_μ_post[k] for k in 1:KMax];
    ll(cell_expression,μk,τk) = [logpdf.(Normal.(μk[k],sqrt.(τk[k] .^ (-1))),cell_expression) for k in 1:length(τk)];
    for t in 1:T
        num_cells = C_t[t]
        cell_ll_scores_t = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            scores = ll(x[t][i],μk,τk)
            cell_ll_scores_t[i] = scores;
        end
        cell_ll_scores[t] = cell_ll_scores_t;
    end
    return cell_ll_scores
end
function calc_cell_normal_μ_τ_l_scores(x,mean_τ_post,mean_μ_post)
    T = length(x);
    C_t = [length(el) for el in x];
    KMax = length(mean_μ_post)
    cell_ll_scores =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T);
    τk = [mean_τ_post[k] for k in 1:KMax];
    μk = [mean_μ_post[k] for k in 1:KMax];
    lh(cell_expression,μk,τk) = [pdf.(Normal.(μk[k],sqrt.(τk[k] .^ (-1))),cell_expression) for k in 1:length(τk)];
    for t in 1:T
        num_cells = C_t[t]
        cell_ll_scores_t = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            scores = lh(x[t][i],μk,τk)
            cell_ll_scores_t[i] = scores;
        end
        cell_ll_scores[t] = cell_ll_scores_t;
    end
    return cell_ll_scores
end
function calc_cell_normal_μ_τ_ll_gene_weights(x,mean_τ_post,mean_μ_post)
    T = length(x);
    C_t = [length(el) for el in x];
    KMax = length(mean_μ_post)
    cell_ll_gene_weights =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    τk = [mean_τ_post[k] for k in 1:KMax];
    μk = [mean_μ_post[k] for k in 1:KMax];
    ll(cell_expression,μk,τk) = [logpdf.(Normal.(μk[k],sqrt.(τk[k] .^ (-1))),cell_expression) for k in 1:length(τk)];
    for t in 1:T
        num_cells = C_t[t]
        cell_ll_gene_weights_t = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            scores = ll(x[t][i],μk,τk)
            gene_weights = norm_weights.(scores)
            cell_ll_gene_weights_t[i]= gene_weights;
        end
        cell_ll_gene_weights[t] = cell_ll_gene_weights_t;
    end
    return cell_ll_gene_weights
end
function calc_cell_normal_μ_τ_ll_gene_weights(cell_ll_scores)
    T = length(cell_ll_scores);
    C_t = [length(el) for el in cell_ll_scores];
    cell_ll_gene_weights =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        num_cells = C_t[t]
        cell_ll_gene_weights_t = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            scores = cell_ll_scores[t][i]
            gene_weights = norm_weights.(scores)
            cell_ll_gene_weights_t[i]= gene_weights;
        end
        cell_ll_gene_weights[t] = cell_ll_gene_weights_t;
    end
    return cell_ll_gene_weights
end
function calc_cell_normal_μ_τ_ll_gene_weightsbygene(cell_ll_scores)
    T = length(cell_ll_scores);
    C_t = [length(el) for el in cell_ll_scores];
    cell_ll_gene_weights =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
    for t in 1:T
        num_cells = C_t[t]
        cell_ll_gene_weights_t = Vector{Vector{Vector{Float64}}}(undef,num_cells)
        for i in 1:num_cells
            scores = cell_ll_scores[t][i]
            scores = collect(eachcol(permutedims(reduce(hcat,scores))))
            gene_weights = norm_weights.(scores)
            cell_ll_gene_weights_t[i]= gene_weights;
        end
        cell_ll_gene_weights[t] = cell_ll_gene_weights_t;
    end
    return cell_ll_gene_weights
end

function calc_gene_importance_weights(x,rtik,mean_τ_post,mean_μ_post,z_post_s;null_precision=10)
    T = length(x)
    C_t = length.(x)
    G = length(x[1][1])
    K = length(rtik[1][1])
    S = length(z_post_s)
    cell_l_scores = calc_cell_normal_μ_τ_l_scores(x,mean_τ_post,mean_μ_post);
    expected_val_cell_cell_l_score = [[[[cell_l_scores[t][i][k][j] .* rtik[t][i][k] for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    expected_val_cell_cell_l_weight = [[[normToProb(expected_val_cell_cell_l_score[t][i][j]) for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    expected_val_cell_cell_l_weight_1 = [[[[expected_val_cell_cell_l_weight[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    # expected_val_cell_cell_l_score = [[[cell_l_scores[t][i][k] .* rtik_[t][i][k]  for k in 1:K] for i in 1:C_t[t]] for t in 1:T]

    # 
    null_cell_l_scores = calc_cell_normal_μ_τ_l_scores(x,[null_precision .* ones(G) for k in 1:K],[zeros(G) for k in 1:K]);
    ration_ = [[[[expected_val_cell_cell_l_score[t][i][j][k] ./ null_cell_l_scores[t][i][k][j]  for k in 1:K] for j in 1:G] for i in 1:C_t[t]] for t in 1:T];
    ration1_ = [[[[ration_[t][i][j][k] for j in 1:G] for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];
    ration1_weight = [[[normToProb(ration1_[t][i][k]) for k in 1:K]  for i in 1:C_t[t]] for t in 1:T];

    gene_significance_weights_mat = Array{Float64}(undef,G,K,S)
    gene_significance_weights = Vector{Vector{Vector{Float64}}}(undef,S)
    for s in 1:S
        clus = Vector{Vector{Vector{Float64}}}(undef,K)
        for k in 1:K
            clus[k] = [zeros(G)]
        end
        for t in 1:T
            for i in 1:C_t[t]
                clus_k = z_post_s[s][t][i]
                # push!(clus[clus_k],expected_val_cell_cell_l_weight_1[t][i][clus_k])
                # push!(clus[clus_k],expected_val_cell_cell_l_score1[t][i][clus_k])
                push!(clus[clus_k],ration1_weight[t][i][clus_k])
                # push!(clus[clus_k],expected_val_cell_cell_l_score1_weight[t][i][clus_k])
    
            end
        end
        N_k = sum(sum.(rtik))
        sum.(clus) ./ N_k
        gene_significance_weights[s] = sum.(clus) ./ N_k
        # println("here")
        gene_significance_weights_mat[:,:,s] = hcat(gene_significance_weights[s]...)
        # println("now here")
    end
    return gene_significance_weights,gene_significance_weights_mat,ration1_weight 
end

function get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
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
function get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,a0_err_hat_vec,b0_err_hat_vec,rtik,v_tikj)
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

function get_bayes_factor(rtik,v_tikj)
    T = length(rtik)
    C_t = length.(rtik)
    K = length(rtik[1][1])
    G= length(v_tikj[1][1][1])
    N_k = sum(sum.(rtik))
    imp_prob = [[[[v_tikj[t][i][k][j][1] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    not_imp_prob = [[[[v_tikj[t][i][k][j][2] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    expected_val_cell_imp_η = [[[imp_prob[t][i][k] .* rtik[t][i][k] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    expected_val_cell_notimp_η = [[[not_imp_prob[t][i][k] .* rtik[t][i][k] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    clus_expected_val_cell_imp_η = [[expected_val_cell_imp_η[t][i][k] for t in 1:T for i in 1:C_t[t]] for k in 1:K];
    clus_expected_val_cell_notimp_η = [[expected_val_cell_notimp_η[t][i][k] for t in 1:T for i in 1:C_t[t]] for k in 1:K];
    clus_sum_expected_val_cell_imp_η = sum.(clus_expected_val_cell_imp_η) ;
    clus_sum_expected_val_cell_notimp_η = sum.(clus_expected_val_cell_notimp_η);
    gene_imp_weights = permutedims(reduce(hcat,[normToProb(col) for col in eachcol(permutedims(reduce(hcat,clus_sum_expected_val_cell_imp_η)))]))
    gene_notimp_weights = permutedims(reduce(hcat,[normToProb(col) for col in eachcol(permutedims(reduce(hcat,clus_sum_expected_val_cell_notimp_η)))]))
    bayes_factor = [v1 ./ v2 for (v1,v2) in zip(clus_sum_expected_val_cell_imp_η,clus_sum_expected_val_cell_notimp_η)]
    return bayes_factor, gene_notimp_weights, gene_imp_weights
end

function dev_calc_inclusion_prob()
        
    clus_expected_val_cell_imp_η  = Vector{Vector{Vector{Float64}}}(undef,K);
    clus_expected_val_cell_notimp_η = Vector{Vector{Vector{Float64}}}(undef,K);
    for k in 1:K
        clus_expected_val_cell_imp_η[k] = [zeros(G)]
        clus_expected_val_cell_notimp_η[k]=[zeros(G)]
    end
    for t in 1:T
        for i in 1:C_t[t]
            clus_k = z_post_s[end][t][i]
            # push!(clus[clus_k],expected_val_cell_cell_l_weight_1[t][i][clus_k])
            # push!(clus[clus_k],expected_val_cell_cell_l_score1[t][i][clus_k])
            push!(clus_expected_val_cell_imp_η[clus_k],expected_val_cell_imp_η[t][i][clus_k])
            push!(clus_expected_val_cell_notimp_η[clus_k],expected_val_cell_notimp_η[t][i][clus_k])

        end
    end
    N_k = sum(sum.(rtik_))
    avg_clus_imp = sum.(clus_expected_val_cell_imp_η) ./ N_k
    avg_clus_notimp = sum.(clus_expected_val_cell_notimp_η) ./ N_k
    avg_clus_imp2 = sum.(clus_expected_val_cell_imp_η) ./ sum(N_k)
    avg_clus_notimp2 = sum.(clus_expected_val_cell_notimp_η) ./ sum(N_k)
    # gene_significance_weights = sum.(clus) ./ N_k
    # gene_significance_weights_mat = hcat(gene_significance_weights...)
    cell_l_scores = calc_cell_normal_μ_τ_l_scores(x_input,mean_τ_post,mean_μ_post);
    imp_prob = [[[[v_tikj_[t][i][k][j][1] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    not_imp_prob = [[[[v_tikj_[t][i][k][j][2] for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    expected_val_cell_imp_η = [[[ imp_prob[t][i][k] .* rtik_[t][i][k] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    expected_val_cell_notimp_η = [[[not_imp_prob[t][i][k] .* rtik_[t][i][k] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
    all_clus_expected_val_cell_imp_η = [[expected_val_cell_imp_η[t][i][k] for t in 1:T for i in 1:C_t[t]] for k in 1:K];
    all_clus_expected_val_cell_notimp_η = [[expected_val_cell_notimp_η[t][i][k] for t in 1:T for i in 1:C_t[t]] for k in 1:K];
    clus_sum_expected_val_cell_imp_η = sum.(all_clus_expected_val_cell_imp_η) ./ sum(sum(sum.(rtik_)));
    clus_sum_expected_val_cell_notimp_η =  sum.(all_clus_expected_val_cell_notimp_η) ./ sum(sum(sum.(rtik_)));
    clus_sum_expected_val_cell_imp_η = [ [va <= 10^-5 ? 0.0 : va  for va in el] for el in clus_sum_expected_val_cell_imp_η]
    clus_sum_expected_val_cell_notimp_η = [ [va <= 10^-5 ? 0.0 : va  for va in el] for el in clus_sum_expected_val_cell_notimp_η]
    gene_imp_weights = permutedims(reduce(hcat,[normToProb(col) for col in eachcol(permutedims(reduce(hcat,clus_sum_expected_val_cell_imp_η)))]))
    gene_notimp_weights = permutedims(reduce(hcat,[normToProb(col) for col in eachcol(permutedims(reduce(hcat,clus_sum_expected_val_cell_notimp_η)))]))
    bayes_factor = [v1 ./ v2 for (v1,v2) in zip(clus_sum_expected_val_cell_imp_η,clus_sum_expected_val_cell_notimp_η)]
 
end


# function calc_cell_normal_μ_τ_ll_gene_weightsbygene(cell_ll_scores)
#     T = length(cell_ll_scores);
#     C_t = [length(el) for el in cell_ll_scores];
#     cell_ll_gene_weights =  Vector{Vector{Vector{Vector{Float64}}}}(undef,T)
#     for t in 1:T
#         num_cells = C_t[t]
#         cell_ll_gene_weights_t = Vector{Vector{Vector{Float64}}}(undef,num_cells)
#         for i in 1:num_cells
#             scores = cell_ll_scores[t][i]
#             scores = collect(eachcol(permutedims(reduce(hcat,scores))))
#             gene_weights = norm_weights.(scores)
#             cell_ll_gene_weights_t[i]= gene_weights;
#         end
#         cell_ll_gene_weights[t] = cell_ll_gene_weights_t;
#     end
#     return cell_ll_gene_weights
# end


function calculate_gene_weight_mean(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=1)
    return mean(gene_weights_clus(cluster_gene_weights,k,genename,post_sample,gene_names,z_post_s))
end
function calculate_gene_weight_variance(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=1)
    return var(gene_weights_clus(cluster_gene_weights,k,genename,post_sample,gene_names,z_post_s))
end
function calculate_gene_weight_stdev(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=1)
    return std(gene_weights_clus(cluster_gene_weights,k,genename,post_sample,gene_names,z_post_s))
end
function calculate_gene_weight_n(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=1)
    return length(gene_weights_clus(cluster_gene_weights,k,genename,post_sample,gene_names,z_post_s))
end
function calculate_gene_weight_posteriorInterval(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=1,pct =0.9)
    up_pct = pct
    lo_pct = 1- up_pct
    n = length(gene_weights_clus(cluster_gene_weights,k,genename,post_sample,gene_names,z_post_s))
    sorted_gene_weights = sort(gene_weights_clus(cluster_gene_weights,k,genename,post_sample,gene_names,z_post_s),rev =true)
    upb_index = ceil(Int,n*up_pct)
    lob_index = ceil(Int,n*lo_pct)
    upb = sorted_gene_weights[upb_index]
    lob = sorted_gene_weights[lob_index]
    return upb,lob
end
function calculate_gene_weight_summary_statistics(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=1,pct =0.9)
    m = calculate_gene_weight_mean(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=post_sample)
    s_sq = calculate_gene_weight_variance(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=post_sample)
    s = calculate_gene_weight_stdev(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=post_sample)
    n = calculate_gene_weight_n(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=post_sample)
    upb,lob = calculate_gene_weight_posteriorInterval(cluster_gene_weights,k,genename,gene_names,z_post_s;post_sample=post_sample,pct =pct)
    return m,s_sq,s,n,upb,lob
end
function calculate_gene_weight_mean(gene_weights_clus_array)
    return mean(gene_weights_clus_array)
end
function calculate_gene_weight_variance(gene_weights_clus_array)
    return var(gene_weights_clus_array)
end
function calculate_gene_weight_stdev(gene_weights_clus_array)
    return std(gene_weights_clus_array)
end
function calculate_gene_weight_n(gene_weights_clus_array)
    return length(gene_weights_clus_array)
end
function calculate_gene_weight_posteriorInterval(gene_weights_clus_array;pct =0.9)
    up_pct = pct
    lo_pct = 1- up_pct
    n = length(gene_weights_clus_array)
    sorted_gene_weights = sort(gene_weights_clus_array,rev =true)
    upb_index = ceil(Int,n*up_pct)
    lob_index = ceil(Int,n*lo_pct)
    upb = sorted_gene_weights[upb_index]
    lob = sorted_gene_weights[lob_index]
    return upb,lob
end
function calculate_gene_weight_summary_statistics(gene_weights_clus_array;pct =0.9)
    m = calculate_gene_weight_mean(gene_weights_clus_array)
    s_sq = calculate_gene_weight_variance(gene_weights_clus_array)
    s = calculate_gene_weight_stdev(gene_weights_clus_array)
    n = calculate_gene_weight_n(gene_weights_clus_array)
    upb,lob = calculate_gene_weight_posteriorInterval(gene_weights_clus_array;pct =pct)
    return m,s_sq,s,n,upb,lob
end
