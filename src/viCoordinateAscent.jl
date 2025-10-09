"""
    log_TrainFeature!(i::Int,training_logger::TrainFeature{U,W},clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
This function logs the training features for a given iteration.
"""
function log_TrainFeature!(i::Int,training_logger::TrainFeature{U,W},clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    T_all = sum(T)
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].m_mu
    end
    training_logger.m_mu[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].s_sq_mu
    end
    training_logger.s_sq_mu[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].y
    end
    training_logger.y[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].h1
    end
    training_logger.h1[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].h2
    end
    training_logger.h2[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].a
    end
    training_logger.a[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].b
    end
    training_logger.b[i] = training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].u
    end
    training_logger.u[i] =  training_logger.JKcache
    training_logger.JKcache .= 0.0
    for k in 1:K
        training_logger.JKcache[:,k] .= clusters[k].v
    end
    training_logger.v[i] = training_logger.JKcache
    training_logger.KplusTcache .= 0.0
    for it in 1:T_all
        training_logger.KplusTcache[:,it] .= conditions[it].d
    end
    training_logger.d[i] = training_logger.KplusTcache
    training_logger.KplusTcache .= 0.0
    training_logger.Tcache .= 0.0
    for it in 1:T_all
        training_logger.Tcache[it] = conditions[it].w1[1]
    end
    training_logger.w1[i] = training_logger.Tcache
    training_logger.Tcache .= 0.0
    for it in 1:T_all
        training_logger.Tcache[it] = conditions[it].w2[1]
    end
    training_logger.w2[i] = training_logger.Tcache
    training_logger.Tcache .= 0.0
    training_logger.Kcache .= 0.0
    for k in 1:K
        training_logger.Kcache[k] = clusters[k].Nk[1]
    end
    training_logger.Nk[i] = training_logger.Kcache
    training_logger.Kcache .= 0.0
    for k in 1:K
        training_logger.Kcache[k] = clusters[k].g1[1]
    end
    training_logger.g1[i] = training_logger.Kcache
    training_logger.Kcache .= 0.0
    for k in 1:K
        training_logger.Kcache[k] = clusters[k].g2[1]
    end
    training_logger.g2[i] = training_logger.Kcache
    training_logger.Kcache .= 0.0
    training_logger.Jcache .= 0.0
    for j in 1:J
        training_logger.Jcache[j] = clusters[1].m_nu[j]
    end
    training_logger.m_nu[i] = training_logger.Jcache
    training_logger.Jcache .= 0.0
    for j in 1:J
        training_logger.Jcache[j] = clusters[1].s_sq_nu[j]
    end
    training_logger.s_sq_nu[i] = training_logger.Jcache
    training_logger.Jcache .= 0.0
    return training_logger
end

"""
    cavi(inputs;delta_rsum_ep = 10^(-6),delta_rsum_lag = 1,delta_occymean_ep = 10^(-6),delta_occymean_lag = 1,logger=nothing,update_clusterwise::Bool = false,elbo_sign_change_max=10.0,check_cluster_interpretability_bool = true,train_h::Bool = true,train_w::Bool = true,train_uv::Bool = true,train_m_nu::Bool = true,train_s_sq_nu::Bool = true,train_ab::Bool = true,multiple_s_sq_updates::Bool = false,burnin::Int = 1,eta_update_mode="Local",sigma_update_mode="Local",lambda_update_mode="Local",remove_small_clusters::Bool = true,use_log::Bool = true)
This is the main variational inference function. It perfroms coordinate ascent to infer the paramters of the NCLUSION Model
"""
function cavi(inputs;delta_rsum_ep = 10^(-6),delta_rsum_lag = 1,delta_occymean_ep = 10^(-6),delta_occymean_lag = 1,logger=nothing,update_clusterwise::Bool = false,elbo_sign_change_max=10.0,check_cluster_interpretability_bool = true,train_h::Bool = true,train_w::Bool = true,train_uv::Bool = true,train_m_nu::Bool = true,train_s_sq_nu::Bool = true,train_ab::Bool = true,multiple_s_sq_updates::Bool = false,burnin::Int = 1,eta_update_mode="Local",sigma_update_mode="Local",lambda_update_mode="Local",remove_small_clusters::Bool = true,use_log::Bool = true)
    # inputs_copy = deepcopy(inputs);
    # inputs = deepcopy(inputs_copy);
    cells,clusters,conditions,matrixconditions,dataparams,modelparams,training_logger  = (; inputs...);
    N = dataparams.N;
    J = dataparams.J;
    K = modelparams.K;
    T = dataparams.T;
    # ##### FOR DEBUGGING Y update only ######
    # r_true = zeros(modelparams.K+1,dataparams.N);
    # for n in 1:dataparams.N
    #   cells[n].r .= zeros(modelparams.K+1);
    #   cells[n].r[cell_cluster_labels[n]] = 1.0;
    #   r_true[cell_cluster_labels[n],n] = 1.0;
    # end
    # empirical_mean_X = hcat([mean(anndata_dict1["X"][:,cell_cluster_labels .== k ], dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_cluster_mean = hcat([mean(used_representation[:,cell_cluster_labels .== k ], dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_cluster_median = hcat([median(used_representation[:,cell_cluster_labels .== k ], dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_gene_mean = mean(used_representation,dims=2)
    # empirical_gene_median = median(used_representation,dims=2)
    # ## empirical_cluster_pip = hcat([mean(used_representation[:,cell_cluster_labels .== k ] .> 0.0, dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_cluster_var = hcat([var(used_representation[:,cell_cluster_labels .== k ], dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_cluster_precision = 1 ./ (empirical_cluster_var .+ 1e-120)
    # empirical_var_notin =  hcat([var(used_representation[:,cell_cluster_labels .!= k ], dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_precision_notin = 1 ./ (empirical_var_notin  .+ 1e-120 )
    # empirical_mean_notin = hcat([mean(used_representation[:,cell_cluster_labels .!= k ], dims=2) for k in sort(unique(cell_cluster_labels))]...)
    # empirical_cluster_pip1 = Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((0 .- empirical_gene_mean) .^2 ./(2 .*1))   .- ((0 .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip2 = Float64.(sigmoid.(log.( (empirical_var_notin .+ 1e-120) ./( empirical_cluster_var .+ 1e-120)) .+ ((empirical_gene_mean .- empirical_mean_notin) .^2 ./(2 .*(empirical_var_notin .+ 1e-120)))   .- ((empirical_gene_mean .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip3 = Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((0.0 .- empirical_mean_notin) .^2 ./(2 .*1.))   .- ((0.0 .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip4 = Float64.(sigmoid.(log.( 1 ./( 1)) .+ ((0.0 .- empirical_mean_notin) .^2 ./(2 .*1.))   .- ((0.0 .- empirical_cluster_mean) .^2 ./(2 .*1))) .>0.5)
    # empirical_cluster_pip5 = Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((empirical_gene_mean .- empirical_mean_notin) .^2 ./(2 .*1.))   .- ((empirical_gene_mean .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip6 = Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((empirical_gene_median .- empirical_mean_notin) .^2 ./(2 .*1.))   .- ((empirical_gene_median .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip7 = Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((0.0 .- empirical_gene_median) .^2 ./(2 .*1.))   .- ((0.0 .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip = Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((0.0 .- empirical_mean_notin) .^2 ./(2 .*1.))   .- ((0.0 .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip =  Float64.(sigmoid.(log.( 1 ./( empirical_cluster_var .+ 1e-120)) .+ ((0.0 .- empirical_gene_median) .^2 ./(2 .*1.))   .- ((0.0 .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    # empirical_cluster_pip =  Float64.(sigmoid.(log.( (empirical_var_notin .+ 1e-120) ./( empirical_cluster_var .+ 1e-120)) .+ ((0.0 .- empirical_gene_median) .^2 ./(2 .*(empirical_var_notin .+ 1e-120)))   .- ((0.0 .- empirical_cluster_mean) .^2 ./(2 .*empirical_cluster_var .+ 1e-120))) .>0.5)
    ## ff =hcat([[mean([el in used_representation_feature_name[sortperm(-empirical_cluster_var[:,k])][1:j] for el in deg[:,"$(k)"]])  for j in 1:length(used_representation_feature_name[sortperm(-empirical_cluster_var[:,k])])] for k in 1:5]...)
    # empirical_global_var = var(used_representation)
    # empirical_cohensd = (empirical_cluster_mean .- empirical_mean_notin) ./ sqrt.((empirical_cluster_var .+ empirical_var_notin) ./2)
    # empirical_ssmd = (empirical_cluster_mean .- empirical_mean_notin) ./ sqrt.(empirical_var_notin .+ 1e-120)
    # empirical_dprime = abs.(empirical_gene_mean .- empirical_cluster_mean) ./ sqrt.((empirical_cluster_var .+ empirical_var_notin) ./2)
    # true_y_all = nothing
    # true_y = nothing
    # if !isnothing(deg)
    #     if all([ "$(el)" in names(deg) for el in unique(cell_cluster_labels)]) &&  all([el in used_representation_feature_name for el in (vec(reshape(Matrix(deg),(size(Matrix(deg))[1]*size(Matrix(deg))[2],1))))])
    #     end
    #     true_y_all = zeros(J,K);
    #     for k in sort(unique(cell_cluster_labels))
    #         for j in 1:J
    #             if used_representation_feature_name[j] in deg[:,"$(k)"]
    #                 true_y_all[j,k] = 1.0;
    #             end
    #         end
    #     end
    # end
    # if !isnothing(true_y_all)
    #     true_y = true_y_all[:,sort(unique(cell_cluster_labels))];
    # end
    ## hcat([empirical_cluster_mean[collect(el) .== 1.0,i] for (i,el) in enumerate(eachcol(true_y))]...)
    ## empirical_cluster_mean[vec(sum(true_y,dims=2).>0.0),:]
    ## hcat([empirical_cluster_var[collect(el) .== 1.0,i] for (i,el) in enumerate(eachcol(true_y))]...)
    ## empirical_cluster_var[vec(sum(true_y,dims=2).>0.0),:]
    ## hcat([empirical_cluster_precision[collect(el) .== 1.0,i] for (i,el) in enumerate(eachcol(true_y))]...)
    ## empirical_cluster_precision[vec(sum(true_y,dims=2).>0.0),:]
    ## empirical_precision_notin[vec(sum(true_y,dims=2).>0.0),:]
    ## marker_genes_j_index = collect(1:J)[vec(sum(true_y,dims=2).>0.0)]
    ## used_representation_feature_name[marker_genes_j_index]
    ## ddd =hcat([[mean([el in used_representation_feature_name[sortperm(-empirical_cohensd[:,k])][1:j] for el in deg[:,"$(k)"]])  for j in 1:length(used_representation_feature_name[sortperm(-empirical_cohensd[:,k])])] for k in 1:5]...)
    ## empirical_cohensd2 = (empirical_cluster_mean .- empirical_gene_mean) ./ sqrt.((empirical_cluster_var .+ empirical_var_notin) ./2)
    ## ddd2 =hcat([[mean([el in used_representation_feature_name[sortperm(-empirical_cohensd2[:,k])][1:j] for el in deg[:,"$(k)"]])  for j in 1:length(used_representation_feature_name[sortperm(-empirical_cohensd2[:,k])])] for k in 1:5]...)
    # for k in sort(unique(cell_cluster_labels))
    #     if !isnothing(deg)
    #         if all([ "$(el)" in names(deg) for el in unique(cell_cluster_labels)]) &&  all([el in used_representation_feature_name for el in (vec(reshape(Matrix(deg),(size(Matrix(deg))[1]*size(Matrix(deg))[2],1))))])
    #             clusters[k].y .= 1/(K*J);#1/100#1.0;#1.0;#
    #             for j in 1:J
    #                 # if used_representation_feature_name[j] in deg[:,"$(k)"]
    #                 #     clusters[k].y[j] = 1.0;
    #                 # end
    #                 clusters[k].y[j] = kappa2*empirical_cluster_pip[j,k];#
    #                 # clusters[k].y[j] = 0.0255*(sigmoid(empirical_ssmd[j,k]) > 0.5);
    #                 # if used_representation_feature_name[j] in reshape(Matrix(deg),(size(Matrix(deg))[1]*size(Matrix(deg))[2],1))
    #                 #     clusters[k].y[j] = 1/J;
    #                 # end
    #                 # clusters[k].y[j] = empirical_cluster_pip[j,k] > median(empirical_cluster_pip);
    #                 # clusters[k].y[j] = empirical_cluster_pip[j,k]
    #                 # if empirical_cluster_mean[j,k] > empirical_gene_mean[j]
    #                 #     clusters[k].y[j] = 0.5;
    #                 # end
    #                 # clusters[k].y[j] = 1 .- sigmoid.((empirical_gene_mean[j] .- empirical_cluster_mean[j,k]).^2 ./ (empirical_cluster_var[j,k]))
    #                 # clusters[k].y[j] = 1 .- sigmoid.((empirical_gene_mean[j] .- empirical_cluster_mean[j,k]).^2 ./ (empirical_cluster_var[j,k]))
    #                 # clusters[k].y[j] = sigmoid.((empirical_gene_mean[j] .- empirical_cluster_mean[j,k]).^2 ./ (empirical_cluster_var[j,k]/empirical_global_var))
    #                 # clusters[k].y[j] = sigmoid.((empirical_gene_mean[j] .- empirical_cluster_mean[j,k]).^1 ./ (empirical_global_var))
    #                 # clusters[k].y[j] = sigmoid.((empirical_gene_mean[j] .- empirical_cluster_mean[j,k]).^2 ./ (empirical_global_var))
    #                 # clusters[k].y[j] = sigmoid.((empirical_gene_median[j] .- empirical_cluster_median[j,k]).^2 ./ (empirical_global_var)) > 0.5
    #             end
    #         end
    #     end
    # end
    # # for k in 1:modelparams.K
    # #  clusters[k].y .= 1/K;
    # # #  clusters[k].a .= 1.0;
    # # #  clusters[k].b .= 1.0;
    # # #  clusters[k].u .= 1.0;
    # # #  clusters[k].v .= 1.0;
    # # end
    # y = hcat([clusters[k].y for k in 1:5]...);
    # [mean([ el in used_representation_feature_name[(y .>= 0.5)[:,k]] for el in deg[:,"$(k)"]]) for k in sort(unique(cell_cluster_labels))]
    # sum(y .>= 0.5,dims=1)
    # sum(sum(y .>= 0.5,dims=2) .> 0)
    # mean([ el in used_representation_feature_name[sum(y .>= 0.5,dims=2) .> 0] for el in vec(reshape(Matrix(deg),(size(Matrix(deg))[1]*size(Matrix(deg))[2],1)))])
    # #m_nu = permutedims(hcat([clusters[1].m_nu[j] for j in 1:J]...))
    # ## all([[el for (j,el) in enumerate(used_representation_feature_name) if clusters[k].y[j] == 1.0] == sort(deg[:,"$(k)"]) for k in sort(unique(cell_cluster_labels))]) ### Just for testing purposes
    #############################
    # for k in 1:modelparams.K
    # #  clusters[k].y .= 1.0;
    # #  clusters[k].a .= 1.0;
    # #  clusters[k].b .= 1.0;
    #  clusters[k].u .= kappa1;
    #  clusters[k].v .= kappa2;
    # end
    # for k in 1:modelparams.K
    #   clusters[k].h1 .= 1.0;
    #   clusters[k].h2 .= 4.0;
    # end
    # for k in 1:modelparams.K
    #   clusters[k].g1[1] = 1/(modelparams.K - (k-1));
    #   clusters[k].g2[1] = 1.0;
    # end 
    # empirical_means = hcat([mean(anndata_dict1["X"][:,cell_cluster_labels .== k ], dims=2) for k in 1:modelparams.K]...)
    # true_means = anndata_dict1["uns"]["mus"] .* anndata_dict1["uns"]["rho"]
    iter = 1;
    elbo_sign_change_counter = 0.0
    elbo_sign_change_max = elbo_sign_change_max
    converged_bool = false;
    is_converged = "false";
    use_log = true;
    cluster_indices = collect(1:modelparams.K);
    noninterpreable_clusters = falses(modelparams.K);
    z_argmax = zeros(Int,length(cells))#[el.z_argmax[1] for el in cells]
    # z_argmax .= [Int(argmax(el.r)) for el in cells]
    # z_argmax .= [Int(argmax(el)) for el in eachcol(r_true)]
    # countmap(z_argmax)
    # Clustering.randindex(cell_cluster_labels,z_argmax)[1]
    # clusters_to_redistribute = falses(modelparams.K);
    num_iter=modelparams.num_iter;
    mean_abs_diff_rsum_lag = zeros(num_iter)
    mean_abs_diff_occymean_lag = zeros(num_iter)
    rsum = zeros(modelparams.K+1,num_iter)
    occymean = zeros(dataparams.J,num_iter)
    rsum_iscoverged = falses(num_iter);
    occymean_iscoverged = falses(num_iter);
    max_continuous_nan_inf_count = 10;
    continuous_nan_inf_counter = 0.0;
    reason_for_non_convergence = "None";
    change_seeds = modelparams.change_seeds;
    current_seed = modelparams.init_seed;
    if change_seeds
        seed_used = "$current_seed, "
    else
        seed_used = "$current_seed"
    end
    log_TrainFeature!(1,training_logger,clusters,conditions,dataparams,modelparams);
    update_d_sum!(conditions,dataparams);
    # d_sum = hcat([conditions[it].d_sum for it in 1:T_all]...);
    #any([any(isnan.(conditions[it].d_sum)) for it in 1:T_all])
    #any([any(isinf.(conditions[it].d_sum)) for it in 1:T_all])   
    #### Basically from mpu_script_v5.jl #######
    update_x_hat!(cells,clusters,dataparams,modelparams);
    #x_hat = hcat([clusters[k].x_hat for k in 1:modelparams.K]...);
    # any([any(isnan.(clusters[k].x_hat)) for k in 1:K])
    # any([any(isinf.(clusters[k].x_hat)) for k in 1:K])
    update_x_hat_sq!(cells,clusters,dataparams,modelparams);
    #x_hat_sq = hcat([clusters[k].x_hat_sq for k in 1:modelparams.K]...);
    #any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])
    #any([any(isinf.(clusters[k].x_hat_sq)) for k in 1:K])
    update_Nk!(cells,clusters,dataparams,modelparams);
    #Nk = permutedims(hcat([clusters[k].Nk for k in 1:modelparams.K]...));
    #any([any(isnan.(clusters[k].Nk)) for k in 1:K])
    #any([any(isinf.(clusters[k].Nk)) for k in 1:K])
    update_Ctt!(cells,conditions,dataparams,modelparams);
    update_CNtk!(cells,matrixconditions,dataparams,modelparams);
    ##################################
    if multiple_s_sq_updates
        update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
        # s_sq_mu = hcat([clusters[k].s_sq_mu for k in 1:modelparams.K]...);
        #any([any(isnan.(clusters[k].s_sq_mu)) for k in 1:K])
        #any([any(isinf.(clusters[k].s_sq_mu)) for k in 1:K])
    end
    while !converged_bool
        _flushed_logger("\t\t\t Starting Iteration $iter...";logger)
        # Local E-STEP
        if train_m_nu
            update_m_nu!(clusters,dataparams,modelparams);
            #m_nu = permutedims(hcat([clusters[1].m_nu[j] for j in 1:J]...));
            #any([any(isnan.(clusters[1].m_nu[j])) for j in 1:J])
            #any([any(isinf.(clusters[1].m_nu[j])) for j in 1:J])
        end
        if train_s_sq_nu
            update_s_sq_nu!(clusters,dataparams,modelparams);
            #s_sq_nu = permutedims(hcat([clusters[1].s_sq_nu[j] for j in 1:J]...));
            #any([any(isnan.(clusters[1].s_sq_nu[j])) for j in 1:J])
            #any([any(isinf.(clusters[1].s_sq_nu[j])) for j in 1:J])
        end
        #### Basically from mpu_script_v5.jl #######
        update_m_mu!(clusters,dataparams,modelparams); #CHANGE?
        # m_mu = hcat([clusters[k].m_mu for k in 1:modelparams.K]...);
        # m_nu .- empirical_cluster_mean
        # mean(empirical_cluster_mean,dims=2 ) .- empirical_cluster_mean
        # m_nu .+ m_mu
        # m_nu .+ y .* m_mu
        #any([any(isnan.(clusters[k].m_mu)) for k in 1:K])
        #any([any(isinf.(clusters[k].m_mu)) for k in 1:K])
        if train_ab ##### not this part but should be the same as in mpu_script_v5.jl when train_ab is true #####
            update_a!(clusters,dataparams,modelparams;sigma_update_mode=sigma_update_mode);
            # a = hcat([clusters[k].a for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].a)) for k in 1:modelparams.K])
            #any([any(isinf.(clusters[k].a)) for k in 1:modelparams.K])
            update_b!(clusters,dataparams,modelparams;sigma_update_mode=sigma_update_mode);
            # b = hcat([clusters[k].b for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].b)) for k in 1:modelparams.K])
            #any([any(isinf.(clusters[k].b)) for k in 1:modelparams.K])
        end
        if train_uv ##### not this part but should be the same as in mpu_script_v5.jl when train_uv is true and update_clusterwise is false #####
            update_u!(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode);
            # u = hcat([clusters[k].u for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].u)) for k in 1:modelparams.K])
            #any([any(isinf.(clusters[k].u)) for k in 1:modelparams.K])
            update_v!(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode);
            # v = hcat([clusters[k].v for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].v)) for k in 1:modelparams.K])
            #any([any(isinf.(clusters[k].v)) for k in 1:modelparams.K])
        end
        if multiple_s_sq_updates
            update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
            # s_sq_mu = hcat([clusters[k].s_sq_mu for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].s_sq_mu)) for k in 1:K])
            #any([any(isinf.(clusters[k].s_sq_mu)) for k in 1:K])
        end
        ##################################
        # hcat([[E_ln_eta(j,clusters[k]) - E_ln_minus_eta(j,clusters[k]) for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * (log(2*pi*clusters[k].s_sq_mu[j]) + 1)  for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * (E_ln_sigma_sq(j,clusters[k]) + E_ln_lambda(j,clusters[k]) - log(clusters[k].s_sq_mu[j]))   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * (E_ln_sigma_sq(j,clusters[k]) + E_ln_lambda(j,clusters[k]) - log(clusters[k].s_sq_mu[j]))   for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * (E_ln_sigma_sq(j,clusters[k]) )   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * (E_ln_lambda(j,clusters[k]))   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * (E_ln_lambda(j,clusters[k]))   for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ - 0.5 * ( - log(clusters[k].s_sq_mu[j]))   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ 0.5 * (clusters[k].m_mu[j]) ^2 /(clusters[k].s_sq_mu[j])   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # used_representation_feature_name[marker_genes_j_index]
        # hcat([[ 0.5 * (clusters[k].m_mu[j]) ^2 /(clusters[k].s_sq_mu[j])   for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # y_tilde = hcat([[E_ln_eta(j,clusters[k]) - E_ln_minus_eta(j,clusters[k]) for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...) .+ hcat([[ - 0.5 * (log(2*pi*clusters[k].s_sq_mu[j]) + 1)  for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...) .+ hcat([[ - 0.5 * (E_ln_sigma_sq(j,clusters[k]) + E_ln_lambda(j,clusters[k]) - log(clusters[k].s_sq_mu[j]))   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...) .+ hcat([[ 0.5 * (clusters[k].m_mu[j]) ^2 /(clusters[k].s_sq_mu[j])   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # y_tilde[marker_genes_j_index,:]
        #  y_int = sigmoid.(hcat([[E_ln_eta(j,clusters[k]) - E_ln_minus_eta(j,clusters[k]) for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...) .+ hcat([[ - 0.5 * (log(2*pi*clusters[k].s_sq_mu[j]) + 1)  for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...) .+ hcat([[ - 0.5 * (E_ln_sigma_sq(j,clusters[k]) + E_ln_lambda(j,clusters[k]) - log(clusters[k].s_sq_mu[j]))   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...) .+ hcat([[ 0.5 * (clusters[k].m_mu[j]) ^2 /(clusters[k].s_sq_mu[j])   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...))
        # y_int[marker_genes_j_index,:]
        # e_ln_lambda(100.,101.)
        # hcat([[ clusters[k].u[j]  for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].v[j]  for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].u[j]  for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].v[j]  for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].m_mu[j]   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].m_mu[j]   for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].m_mu[j] ^2   for j in marker_genes_j_index] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        # hcat([[ clusters[k].s_sq_mu[j]   for j in 1:J] for k in sort(unique([Int(argmax(el.r)) for el in cells]))]...)
        update_y!(clusters,dataparams,modelparams); #CHANGE?
        # y = hcat([clusters[k].y for k in 1:modelparams.K]...);
        # [mean([ el in used_representation_feature_name[(y .>= 0.5)[:,k]] for el in deg[:,"$(k)"]]) for k in sort(unique(cell_cluster_labels))]
        # sum(y[:,1:5],dims=1)
        # sum(sum(y[:,1:5],dims=2) .>= 1)
        # sum(sum(y[:,1:5] .>= 0.5,dims=2) .>= 1)
        #any([any(isnan.(clusters[k].y)) for k in 1:K])
        #any([any(isinf.(clusters[k].y)) for k in 1:K])
        update_d!(clusters, conditions,matrixconditions, dataparams, modelparams; use_log =use_log);
        # d = hcat([conditions[it].d for it in 1:T_all]...);
        #any([any(isnan.(conditions[it].d)) for it in 1:T_all])
        #any([any(isinf.(conditions[it].d)) for it in 1:T_all])   
        update_d_sum!(conditions,dataparams);
        # d_sum = hcat([conditions[it].d_sum for it in 1:T_all]...);
        #any([any(isnan.(conditions[it].d_sum)) for it in 1:T_all])
        #any([any(isinf.(conditions[it].d_sum)) for it in 1:T_all])   
        update_c!(cells,conditions,dataparams,modelparams);
        # c = permutedims(hcat([cells[n].c for n in 1:N]...));
        # any([any(isnan.(cells[n].c)) for n in 1:N])
        # any([any(isinf.(cells[n].c)) for n in 1:N])
        update_r!(cells,clusters,conditions,dataparams,modelparams);
        # r = permutedims(hcat([cells[n].r for n in 1:N]...));
        # sum(r,dims=1)
        # any([any(isnan.(cells[n].r)) for n in 1:N])
        # any([any(isinf.(cells[n].r)) for n in 1:N])
        # z_argmax .= [Int(el.z_argmax[1]) for el in cells]
        z_argmax .= [Int(el.z_argmax[1]) for el in cells]
        _flushed_logger("\t\t\t\t\t\t\t\tNumber of occupied clusters: $(length(unique(z_argmax)))";logger)
        # _flushed_logger("\t\t\t\t\t\t\t\tPips in occupied clusters: $(hcat([clusters[k].y for k in 1:modelparams.K]...)[:,sort(unique(z_argmax))])";logger)
        # unique(z_argmax)
        # countmap(z_argmax)
        # Clustering.randindex(cell_cluster_labels,z_argmax)[1]
        ##################################
        # Global E-STEP
        #### Basically from mpu_script_v5.jl #######
        update_x_hat!(cells,clusters,dataparams,modelparams);
        #x_hat = hcat([clusters[k].x_hat for k in 1:modelparams.K]...);
        # any([any(isnan.(clusters[k].x_hat)) for k in 1:K])
        # any([any(isinf.(clusters[k].x_hat)) for k in 1:K])
        update_x_hat_sq!(cells,clusters,dataparams,modelparams);
        #x_hat_sq = hcat([clusters[k].x_hat_sq for k in 1:modelparams.K]...);
        #any([any(isnan.(clusters[k].x_hat_sq)) for k in 1:K])
        #any([any(isinf.(clusters[k].x_hat_sq)) for k in 1:K])
        update_Nk!(cells,clusters,dataparams,modelparams);
        #Nk = permutedims(hcat([clusters[k].Nk for k in 1:modelparams.K]...));
        #any([any(isnan.(clusters[k].Nk)) for k in 1:K])
        #any([any(isinf.(clusters[k].Nk)) for k in 1:K])
        rsum[:,iter] = vec(sum(permutedims(hcat([cells[n].r for n in 1:N]...)),dims=1));
        occymean[:,iter] = vec(mean(hcat([clusters[k].y for k in 1:modelparams.K]...)[:,unique(z_argmax)],dims=2))
        if iter > delta_rsum_lag
            mean_abs_diff_rsum_lag[iter] = mean(abs.(rsum[:,iter] - rsum[:,iter-delta_rsum_lag]))
            if mean_abs_diff_rsum_lag[iter] < delta_rsum_ep
                rsum_iscoverged[iter] = true
            end
        end
        if iter > delta_occymean_lag
            mean_abs_diff_occymean_lag[iter] = mean(abs.(occymean[:,iter] - occymean[:,iter-delta_occymean_lag]))
            if mean_abs_diff_occymean_lag[iter] < delta_occymean_ep
                occymean_iscoverged[iter] = true
            end
        end
        if multiple_s_sq_updates                 ##### not this part but should be the same as in mpu_script_v5.jl when multiple_s_sq_updates is true #####
            update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
            # s_sq_mu = hcat([clusters[k].s_sq_mu for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].s_sq_mu)) for k in 1:K])
            #any([any(isinf.(clusters[k].s_sq_mu)) for k in 1:K])
        end
        update_Ctt!(cells,conditions,dataparams,modelparams);
        #any([any(isnan.(conditions[it].Ctt)) for it in 1:T_all])
        #any([any(isinf.(conditions[it].Ctt)) for it in 1:T_all])   
        update_CNtk!(cells,matrixconditions,dataparams,modelparams);
        #any([any(isnan.(matrixconditions[it].CNtk)) for it in 1:T_all])
        #any([any(isinf.(matrixconditions[it].CNtk)) for it in 1:T_all])
        update_alpha_Tk_stats!(clusters,conditions,dataparams,modelparams);
        #any([any(isnan.(clusters[k].alpha_Tk)) for k in 1:K])
        #any([any(isinf.(clusters[k].alpha_Tk)) for k in 1:K])
        ##################################
        # Global M-STEP
         #### Basically from mpu_script_v5.jl #######
        # if multiple_s_sq_updates ##### not this part but should be the same as in mpu_script_v5.jl when multiple_s_sq_updates is true #####
        #     update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
        #     # s_sq_mu = hcat([clusters[k].s_sq_mu for k in 1:modelparams.K]...);
        #     #any([any(isnan.(clusters[k].s_sq_mu)) for k in 1:K])
        #     #any([any(isinf.(clusters[k].s_sq_mu)) for k in 1:K])
        # end
        if train_h
            update_h1!(clusters,dataparams,modelparams; eta_update_mode=eta_update_mode);
            # h1 = hcat([clusters[k].h1 for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].h1)) for k in 1:modelparams.K])
            #any([any(isinf.(clusters[k].h1)) for k in 1:modelparams.K])
            update_h2!(clusters,dataparams,modelparams; eta_update_mode=eta_update_mode);
            # h2 = hcat([clusters[k].h2 for k in 1:modelparams.K]...);
            #any([any(isnan.(clusters[k].h2)) for k in 1:modelparams.K])
            #any([any(isinf.(clusters[k].h2)) for k in 1:modelparams.K])
        end
        if train_w
            update_w1!(conditions, dataparams, modelparams);
            # w1 = permutedims(hcat([conditions[it].w1 for it in 1:T_all]...));
            #any([any(isnan.(conditions[it].w1)) for it in 1:T_all])
            #any([any(isinf.(conditions[it].w1)) for it in 1:T_all])
            update_w2!(conditions, dataparams, modelparams);
            # w2 = permutedims(hcat([conditions[it].w2 for it in 1:T_all]...));
            #any([any(isnan.(conditions[it].w2)) for it in 1:T_all])
            #any([any(isinf.(conditions[it].w2)) for it in 1:T_all])
        end
        update_g1g2!(clusters,dataparams,modelparams;use_log =use_log);
        # g1 = permutedims(hcat([clusters[k].g1 for k in 1:modelparams.K]...));
        # g2 = permutedims(hcat([clusters[k].g2 for k in 1:modelparams.K]...));
        #any([any(isnan.(clusters[k].g1)) for k in 1:K])
        #any([any(isinf.(clusters[k].g1)) for k in 1:K])
        #any([any(isnan.(clusters[k].g2)) for k in 1:K])
        #any([any(isinf.(clusters[k].g2)) for k in 1:K])
        # Calculate ELBO
        iter = Int64(iter)
        LB =  ELBO(cells,clusters,conditions,dataparams,modelparams;use_log= true,update_clusterwise=update_clusterwise,eta_update_mode=eta_update_mode);
        if isinf(LB) || isnan(LB)
            _flushed_logger("\t\t\t  WARNING: ELBO IS NAN OR INF...";logger)
            if continuous_nan_inf_counter > 0
                _flushed_logger("\t\t\t  ELBO IS NAN OR INF for $(Int(continuous_nan_inf_counter)) times in a row...";logger)
            end
            continuous_nan_inf_counter += 1
            if continuous_nan_inf_counter > max_continuous_nan_inf_count
                _flushed_logger("\t\t\t  ELBO IS NAN OR INF TOO MANY TIMES... Breaking";logger)
                is_converged = "false"
                reason_for_non_convergence = ""
                reason_for_non_convergence = reason_for_non_convergence * "WARNING: ELBO IS NAN OR INF too many times; "
                converged_bool = true
                # iter += 1
                continue
            end
            if check_cluster_interpretability_bool
                check_cluster_interpretability!(noninterpreable_clusters,clusters,dataparams,modelparams)#;,clusters_to_redistribute
                if change_seeds
                    current_seed = Int(current_seed+1)
                    _flushed_logger("\t\t\t  CHANGING SEED. NEW SEED IS $current_seed...";logger)
                    seed_used = seed_used * "$current_seed, "
                end
                reset_clusters!(cluster_indices[noninterpreable_clusters],cells,clusters,conditions,matrixconditions,dataparams,modelparams,change_seeds,current_seed)#,
                iter += 1
                continue
            end
            if change_seeds
                current_seed = Int(current_seed+1)
                _flushed_logger("\t\t\t  CHANGING SEED. NEW SEED IS $current_seed...";logger)
                seed_used = seed_used * "$current_seed, "
            end
            reset_clusters!(cluster_indices,cells,clusters,conditions,matrixconditions,dataparams,modelparams,change_seeds,current_seed)
            continue
            # break
        else
            if continuous_nan_inf_counter > 0
                _flushed_logger("\t\t\t\t\t Resetting NAN/INF ELBO counter ...";logger)
            end
            continuous_nan_inf_counter = 0
        end
        training_logger.elbo_[iter] = LB
        log_TrainFeature!(Int64(iter+1),training_logger,clusters,conditions,dataparams,modelparams)
        # println("Iteration: $iter, ELBO: $LB")
        _flushed_logger("\t\t\t\t Iteration $iter Completed, Interation ELBO: $LB";logger)
        if iter > burnin+1
            delta_elbo = abs(training_logger.elbo_[iter] - training_logger.elbo_[iter-1])
            if (occymean_iscoverged[iter] && rsum_iscoverged[iter]) || iter>=num_iter || continuous_nan_inf_counter > max_continuous_nan_inf_count
                converged_bool = true
                if iter>=num_iter || continuous_nan_inf_counter > max_continuous_nan_inf_count
                    is_converged = "false"
                    reason_for_non_convergence = ""
                    if iter>=num_iter
                        reason_for_non_convergence = reason_for_non_convergence * "WARNING: Max Iterations Reached; "
                    # elseif delta_elbo > elbo_ep
                    #     reason_for_non_convergence = reason_for_non_convergence * "WARNING: ELBO Change Below Threshold; "
                    elseif continuous_nan_inf_counter >= max_continuous_nan_inf_count
                        reason_for_non_convergence = reason_for_non_convergence * "WARNING: ELBO IS NAN OR INF too many times; "
                    end
                else
                    is_converged = "true"
                    # println("\t Convergence reached at: $iter")
                    _flushed_logger("\t\t\t Convergence reached at: $iter...";logger)
                end
            end
            delta_elbo_sign = training_logger.elbo_[iter] - training_logger.elbo_[iter-1]
            if delta_elbo_sign < 0
                elbo_sign_change_counter += 1
                # println("\t\t  ELBO IS DECREASING!")
                if elbo_sign_change_counter > elbo_sign_change_max
                    # println("\t\t\t  ELBO DECREASED TOO MUCH!")
                    _flushed_logger("\t\t\t\t\t  ELBO DECREASED TOO MUCH! ";logger)
                    is_converged = "false"
                    reason_for_non_convergence = reason_for_non_convergence * "WARNING: ELBO decreased too much; "
                    converged_bool = true
                end
            else
                if elbo_sign_change_counter > 0
                    # println("\t\t RESET ELBO SIGN COUNTER")
                    _flushed_logger("\t\t\t\t\t RESET ELBO SIGN COUNTER";logger)
                    elbo_sign_change_counter = 0
                end
            end
        end
        iter += 1
    end
    # Get final values
    update_x_hat!(cells,clusters,dataparams,modelparams);
    update_x_hat_sq!(cells,clusters,dataparams,modelparams);
    update_Nk!(cells,clusters,dataparams,modelparams);
    update_Ctt!(cells,conditions,dataparams,modelparams);
    update_CNtk!(cells,matrixconditions,dataparams,modelparams);
    if train_m_nu
        update_m_nu!(clusters,dataparams,modelparams);
    end
    if train_s_sq_nu
        update_s_sq_nu!(clusters,dataparams,modelparams);
    end
    if train_ab    
        update_a!(clusters,dataparams,modelparams;sigma_update_mode=sigma_update_mode);
        update_b!(clusters,dataparams,modelparams;sigma_update_mode=sigma_update_mode);
    end
    if train_uv
        update_u!(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode);
        update_v!(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode);
    end
    update_m_mu!(clusters,dataparams,modelparams); #CHANGE?
    update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
    update_y!(clusters,dataparams,modelparams); #CHANGE?
    update_d!(clusters, conditions,matrixconditions, dataparams, modelparams; use_log =use_log);
    update_d_sum!(conditions,dataparams);
    if train_h
        update_h1!(clusters,dataparams,modelparams; eta_update_mode=eta_update_mode);
        update_h2!(clusters,dataparams,modelparams; eta_update_mode=eta_update_mode);
    end
    if train_w
        update_w1!(conditions, dataparams, modelparams);
        update_w2!(conditions, dataparams, modelparams);
    end
    nonemptychain_indx = broadcast(!,ismissing.(training_logger.elbo_) .|| isnan.(training_logger.elbo_)) 
    elbo_ = training_logger.elbo_[nonemptychain_indx]
    truncation_value = length(elbo_) #+ 1
    output_str_list1 = @name elbo_;
    output_key_list1 = Symbol.(naming_vec(output_str_list1));
    output_var_list1 = [elbo_];
    LB = elbo_;
    outputs_dict = OrderedDict{Symbol,Any}();
    addToDict!(outputs_dict,output_key_list1,output_var_list1);
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cells,clusters,conditions,dataparams,modelparams,training_logger);
    output_str_list2 = @name is_converged,truncation_value,LB,seed_used,reason_for_non_convergence;
    output_key_list2 = Symbol.(naming_vec(output_str_list2));
    output_var_list2 = [is_converged,truncation_value,LB,seed_used,reason_for_non_convergence];
    addToDict!(outputs_dict,output_key_list2,output_var_list2);
    training_features_dict = extract_training_features(training_logger,nonemptychain_indx);
    return outputs_dict,training_features_dict
end


"""
    testting_convergence_monitoring(inputs;num_iter=100,use_log=true,eta_update_mode="Global",lambda_update_mode="Local",sigma_update_mode="Local")
A function to test convergence monitoring by tracking mean absolute differences in `rsum` and `occymean` over iterations.
"""
function testting_convergence_monitoring(inputs;num_iter=100,use_log=true,eta_update_mode="Global",lambda_update_mode="Local",sigma_update_mode="Local")
    inputs_copy = deepcopy(inputs);
    cells,clusters,conditions,matrixconditions,dataparams,modelparams,training_logger  = (; inputs_copy...);
    update_d_sum!(conditions,dataparams);
    update_x_hat!(cells,clusters,dataparams,modelparams);
    update_x_hat_sq!(cells,clusters,dataparams,modelparams);
    update_Nk!(cells,clusters,dataparams,modelparams);
    update_Ctt!(cells,conditions,dataparams,modelparams);
    update_CNtk!(cells,matrixconditions,dataparams,modelparams);
    update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
    mean_abs_diff_rsum_lag1 = zeros(num_iter)
    mean_abs_diff_rsum_lag2 = zeros(num_iter)
    mean_abs_diff_rsum_lag3 = zeros(num_iter)
    mean_abs_diff_occymean_lag1 = zeros(num_iter)
    mean_abs_diff_occymean_lag2 = zeros(num_iter)
    mean_abs_diff_occymean_lag3 = zeros(num_iter)
    rsum = zeros(modelparams.K+1,num_iter)
    occymean = zeros(dataparams.J,num_iter)
    for iter in 1:num_iter
        update_m_nu!(clusters,dataparams,modelparams);
        update_s_sq_nu!(clusters,dataparams,modelparams);
        update_m_mu!(clusters,dataparams,modelparams); #CHANGE?
        update_a!(clusters,dataparams,modelparams;sigma_update_mode=sigma_update_mode);
        update_b!(clusters,dataparams,modelparams;sigma_update_mode=sigma_update_mode);
        update_y!(clusters,dataparams,modelparams); #CHANGE?
        update_d!(clusters, conditions,matrixconditions, dataparams, modelparams; use_log =use_log);
        update_d_sum!(conditions,dataparams);
        update_c!(cells,conditions,dataparams,modelparams);
        update_r!(cells,clusters,conditions,dataparams,modelparams);
        z_argmax .= [Int(el.z_argmax[1]) for el in cells]
        _flushed_logger("\t\t\t\t\t\t\t\tNumber of occupied clusters: $(length(unique(z_argmax)))";logger)
        rsum[:,iter] = vec(sum(permutedims(hcat([cells[n].r for n in 1:N]...)),dims=1));
        occymean[:,iter] = vec(mean(hcat([clusters[k].y for k in 1:modelparams.K]...)[:,unique(z_argmax)],dims=2))
        if iter > 1
            mean_abs_diff_rsum_lag1[iter] = mean(abs.(rsum[:,iter] - rsum[:,iter-1]))
            # mean_abs_diff_rsum_lag2[iter] = mean(abs.(rsum_lag1 - rsum_lag2))
            # mean_abs_diff_rsum_lag3[iter] = mean(abs.(rsum_lag2 - rsum_lag3))
            mean_abs_diff_occymean_lag1[iter] = mean(abs.(occymean[:,iter] - occymean[:,iter-1]))
            # mean_abs_diff_occymean_lag2[iter] = mean(abs.(occymean_lag1 - occymean_lag2))
            # mean_abs_diff_occymean_lag3[iter] = mean(abs.(occymean_lag2 - occymean_lag3))
            _flushed_logger("\t\t\t\t\t\t\t\t\t\t mean_abs_diff_rsum_lag1 for iteration $(iter): $(mean_abs_diff_rsum_lag1[iter])";logger)
            _flushed_logger("\t\t\t\t\t\t\t\t\t\t mean_abs_diff_occymean_lag1 for iteration $(iter): $(mean_abs_diff_occymean_lag1[iter])";logger)
        end
        if iter > 2
            mean_abs_diff_rsum_lag2[iter] = mean(abs.(rsum[:,iter] - rsum[:,iter-2]))
            mean_abs_diff_occymean_lag2[iter] = mean(abs.(occymean[:,iter] - occymean[:,iter-2]))
            _flushed_logger("\t\t\t\t\t\t\t\t\t\t mean_abs_diff_rsum_lag2 for iteration $(iter): $(mean_abs_diff_rsum_lag2[iter])";logger)
            _flushed_logger("\t\t\t\t\t\t\t\t\t\t mean_abs_diff_occymean_lag2 for iteration $(iter): $(mean_abs_diff_occymean_lag2[iter])";logger)
        end
        if iter > 3
            mean_abs_diff_rsum_lag3[iter] = mean(abs.(rsum[:,iter] - rsum[:,iter-3]))
            mean_abs_diff_occymean_lag3[iter] = mean(abs.(occymean[:,iter] - occymean[:,iter-3]))
            _flushed_logger("\t\t\t\t\t\t\t\t\t\t mean_abs_diff_rsum_lag3 for iteration $(iter): $(mean_abs_diff_rsum_lag3[iter])";logger)
            _flushed_logger("\t\t\t\t\t\t\t\t\t\t mean_abs_diff_occymean_lag3 for iteration $(iter): $(mean_abs_diff_occymean_lag3[iter])";logger)
        end
        update_x_hat!(cells,clusters,dataparams,modelparams);
        update_x_hat_sq!(cells,clusters,dataparams,modelparams);
        update_Nk!(cells,clusters,dataparams,modelparams);
        update_s_sq_mu!(clusters,dataparams,modelparams); #CHANGE?
        update_Ctt!(cells,conditions,dataparams,modelparams);
        update_CNtk!(cells,matrixconditions,dataparams,modelparams);
        update_alpha_Tk_stats!(clusters,conditions,dataparams,modelparams);
        update_u!(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode);
        update_v!(clusters,dataparams,modelparams;lambda_update_mode=lambda_update_mode);
        update_h1!(clusters,dataparams,modelparams; eta_update_mode=eta_update_mode);
        update_h2!(clusters,dataparams,modelparams; eta_update_mode=eta_update_mode);
        update_w1!(conditions, dataparams, modelparams);
        update_w2!(conditions, dataparams, modelparams);
        update_g1g2!(clusters,dataparams,modelparams;use_log =use_log);
    end

end

"""
    extract_training_features(training_logger,nonemptychain_indx)
A function to extract training features from the training logger after removing any missing or NaN values.
"""
function extract_training_features(training_logger,nonemptychain_indx)
    training_features_dict = OrderedDict{Symbol,Any}()
    elbo_ = training_logger.elbo_
    training_features_dict[:LB_] = elbo_[nonemptychain_indx]
    params_of_interest = [:y, :m_mu, :s_sq_mu, :m_nu, :s_sq_nu, :d, :g1, :g2,:h1, :h2,:Nk,:w1, :w2,:a, :b,:u, :v]
    for param in params_of_interest
        new_key = Symbol(String(param) * "_")
        params = getfield(training_logger,param)
        training_features_dict[new_key] = params[nonemptychain_indx]
    end
    return training_features_dict
end


"""
    check_cluster_interpretability!(noninterpreable_clusters::BitVector,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #clusters_to_redistribute::BitVector,
A function to check the interpretability of clusters based on certain criteria.
"""
function check_cluster_interpretability!(noninterpreable_clusters::BitVector,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #clusters_to_redistribute::BitVector,
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    significance_prop= modelparams.significance_prop
    min_number_cells=modelparams.min_number_cells
    max_number_cells =0.5*N
    min_percent_of_genes =modelparams.min_percent_of_genes
    max_percent_of_genes =modelparams.max_percent_of_genes
    Kplus = K + 1
    max_number_genes = round(Int,J*max_percent_of_genes)
    min_number_genes = round(Int,J*min_percent_of_genes)
    # ;significance_prop=0.5,min_number_cells=1,min_number_genes = 5,max_percent_of_genes = 0.75
    # cluster_indices = collet(1:K)
    # noninterpreable_clusters = Vector{Bool}(undef,K)
    noninterpreable_cluster_counter = 0
    for k in 1:K
        noninterpreable_clusters[k] = false
        # clusters_to_redistribute[k] = false
        cond1 = (clusters[k].Nk[1] > min_number_cells && all(clusters[k].y .>= significance_prop))
        # cond2 = (clusters[k].Nk[1] > min_number_cells && all(clusters[k].y .< significance_prop) &&  sum(clusters[k].y .< significance_prop) > min_number_genes)
        cond2 = (clusters[k].Nk[1] > min_number_cells && all(clusters[k].y .< significance_prop))
        cond3 = (clusters[k].Nk[1] > min_number_cells &&  sum(clusters[k].y .< significance_prop) < min_number_genes)
        cond4 = (clusters[k].Nk[1] > min_number_cells &&  sum(clusters[k].y .>= significance_prop) > max_number_genes)
        # cond5 = ((clusters[k].Nk[1] <= min_number_cells && clusters[k].Nk[1] >= 1.0 )||clusters[k].Nk[1] >= max_number_cells )# && any(clusters[k].y .>= significance_prop) && sum(clusters[k].y .< significance_prop) > min_number_genes && sum(clusters[k].y .>= significance_prop) < max_number_genes
        if cond1 || cond2 || cond4 || cond3 #|| cond5
            # println("\t\t\t\t\t\t  CLUSTER: $k, Nk: $(clusters[k].Nk[1]), number of signifcant pips: $(sum(clusters[k].y .>= significance_prop))")
            # println("\t\t\t\t\t\t  Condition 1: $(cond1), Condition 2: $(cond2), Condition 3: $(cond3), Condition 4: $(cond4), Condition 5: $(cond5)")
            noninterpreable_clusters[k] = true
            noninterpreable_cluster_counter+=1
        else
            noninterpreable_clusters[k] = false
        end
    end
    println("\t\t\t\t\t\t  Number of noninterpretable clusters: $noninterpreable_cluster_counter")
    return noninterpreable_clusters #,clusters_to_redistribute
end
# ::Vector{CellFeature{U,W,J}}

"""
    reset_clusters!(cluster_indices::Vector{Int},cells,clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature,change_seeds::Bool,current_seed::Int) where {U <: AbstractFloat, W <: Int64}
A function to reset specified clusters by reinitializing their parameters and updating related statistics.
"""
function reset_clusters!(cluster_indices::Vector{Int},cells,clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature,change_seeds::Bool,current_seed::Int) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    # J = dataparams.J
    # T = dataparams.T
    # N = dataparams.N
    # K = modelparams.K
    # significance_prop= modelparams.significance_prop
    # min_number_cells=modelparams.min_number_cells
    # max_number_cells =0.5*N
    # min_percent_of_genes =modelparams.min_percent_of_genes
    # max_percent_of_genes =modelparams.max_percent_of_genes
    # Kplus = K + 1
    # max_number_genes = round(Int,J*max_percent_of_genes)
    # min_number_genes = round(Int,J*max_percent_of_genes)
    # max_small_clusters_to_redistribute = 3
    # max_large_clusters_to_redistribute = 1
    # redistribute_small_clusters_counter = 0
    # redistribute_large_clusters_counter = 0
    if change_seeds
        Random.seed!(current_seed)
    end
    println("\t\t\t\t\t\t  RESETTING CLUSTER: $cluster_indices")
    for k in cluster_indices
        # if ((clusters[k].Nk[1] <= min_number_cells && clusters[k].Nk[1] >= 1.0 ) ||clusters[k].Nk[1] >= max_number_cells )
        #     # println("\t\t\t\t\t\t  RESETTING d for CLUSTER: $k")
        #     for t in 1:T
        #         conditions[t].d[k] = 0.000001#exp(randn())
        #     end
        #     update_d_sum!(conditions,dataparams)
        # end
        clusters[k].Nk[1] = 0.0
        clusters[k].x_hat .= randn(J)
        clusters[k].x_hat_sq .= exp.(randn(J))
        clusters[k].m_mu .= randn(J)
        clusters[k].s_sq_mu .= exp.(randn(J))
        clusters[k].g1[1] = logistic(randn())
        clusters[k].g2[1] = exp(randn())
        clusters[k].h1[1] = exp(randn())
        clusters[k].h2[1] = exp(randn())
        clusters[k].y .= logistic.(randn(J))
    end
    # println("\t\t\t\t\t\t  REDISTRIBUTING CLUSTER: $clusters_to_redistribute")


    # collapsed_x_hat = randn(J)
    # collapsed_x_hat_sq = exp.(randn(J))
    # collapsed_m = randn(J)
    # collapsed_s_sq = exp.(randn(J))
    # collapsed_y = logistic.(randn(J))
    # collapsed_g1 = logistic(randn())
    # collapsed_g2 = exp(randn())
    # collapsed_h1 = exp(randn())
    # collapsed_h2 = exp(randn())
    # for k in cluster_indices
    #     clusters[k].Nk[1] = 0.0
    #     clusters[k].x_hat .= collapsed_x_hat
    #     clusters[k].x_hat_sq .= collapsed_x_hat_sq
    #     clusters[k].m .= collapsed_m
    #     clusters[k].s_sq .= collapsed_s_sq
    #     clusters[k].g1[1] = collapsed_g1
    #     clusters[k].g2[1] = collapsed_g2
    #     clusters[k].h1[1] = collapsed_h1
    #     clusters[k].h2[1] = collapsed_h2
    #     clusters[k].y .= collapsed_y
    # end
    update_x_hat!(cells,clusters,dataparams,modelparams);
    update_x_hat_sq!(cells,clusters,dataparams,modelparams);
    update_Nk!(cells,clusters,dataparams,modelparams);
    update_Ctt!(cells,conditions,dataparams,modelparams);
    update_CNtk!(cells,matrixconditions,dataparams,modelparams);
    return cells,clusters,conditions
end
# ::Vector{CellFeature{U,W,J}},
"""
    redistribute_r!(cells,clusters_to_exclude::Vector{Int}, clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #  formerly update_rtik_mpu!
A function to redistribute the responsibility weights `r` for each cell, excluding specified clusters.
"""
function redistribute_r!(cells,clusters_to_exclude::Vector{Int}, clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #  formerly update_rtik_mpu!
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    # @inbounds for n in 1:N
    Threads.@threads for n in 1:N
        # t = first(dataparams.LinearAddress[n])
        # cells[n]._reset!(cells[n].cache,cells[n].BitType)
        cells[n].cache .= 0.0
        adjust_E_ln_pi!(cells[n],conditions)
        # cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
        for k in 1:K
            # cells[n].cache .= 0.0
            E_log_normal_l_j!(cells[n],clusters[k], dataparams)
            @inbounds cells[n].r[k] = cells[n].cache[k]
        end
        zeroing_value = -1.5*maximum([abs(minimum(cells[n].r)), abs(maximum(cells[n].r))])
        for k in clusters_to_exclude
            cells[n].r[k] = zeroing_value
        end
        norm_weights3!(K,cells[n].r)
        if any(isnan.(cells[n].r))
            println("NAN IN R")
        end
        if any(isinf.(cells[n].r))
            println("INF IN R")
        end
    end
    return cells
end

"""
    remove_small_clusters!(z_argmax::Vector{Int},cells,clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
A function to remove small clusters by setting their responsibility weights to zero and updating related statistics.
"""
function remove_small_clusters!(z_argmax::Vector{Int},cells,clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    min_number_cells = modelparams.min_number_cells
    cluster_counts = countmap(z_argmax)
    small_clusters = [key for key in keys(cluster_counts) if cluster_counts[key] < 20]
    if length(small_clusters) == 0
        return cells,clusters,conditions,matrixconditions,dataparams,modelparams
    end
    println("\t\t\t\t\t\t  Removing small clusters: $small_clusters")
    # for k in small_clusters
    #     clusters[k].g1[1] = 1e-12
    # end
    Threads.@threads for n in 1:N
        for k in small_clusters
            cells[n].r[k] = 0.0
        end
        norm_weights3!(K,cells[n].r)
    end
    update_x_hat!(cells,clusters,dataparams,modelparams);
    update_x_hat_sq!(cells,clusters,dataparams,modelparams);
    update_Nk!(cells,clusters,dataparams,modelparams);
    update_Ctt!(cells,conditions,dataparams,modelparams);
    update_CNtk!(cells,matrixconditions,dataparams,modelparams);
    return cells,clusters,conditions,matrixconditions,dataparams,modelparams

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
    extract_scalars_paramter(paramname,scalars,dataparams)
This function extracts the scalar specific parameters from the Scalars object used in inference
"""
function extract_scalars_paramter(paramname,scalars,dataparams)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    T = dataparams.T

    param = isone(length(getfield(scalars,paramname))) ? getfield(scalars,paramname)[1] : getfield(cscalers,paramname) 
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
    J = dataparams.J
    param = [ isone(length(getfield(geneparams[j],paramname))) ? getfield(geneparams[j],paramname)[1] : getfield(geneparams[j],paramname) for j in 1:J]
    return param
end


"""
        extract_elbo_vals_perK(paramname,training_logger)
This function extracts the cluster specific elbo values from the TrainFeature object used in inference
"""
function extract_elbo_vals_perK(paramname,training_logger)
    if typeof(paramname) <: String
        paramname = Symbol(paramname)
    end
    # G = dataparams.G

    vals = getfield(training_logger,paramname)
    return vals
end
"""
    extract_r_paramter(cells,dataparams)
This function extracts the cell-level cluster probability vector for each cell in the CellFeatures object
"""
function extract_r_paramter(cells,dataparams)

    paramname = :_r
    I = dataparams.I
    T = dataparams.T
    N_t = dataparams.N_t
    r = Vector{Vector{Vector{Vector{cells[1].BitType}}}}(undef,I)
    n = 0
    for i in 1:I
        r[i] = Vector{Vector{Vector{cells[1].BitType}}}(undef,T[i])
        for t in 1:T[i]
            r[i][t] = Vector{Vector{cells[1].BitType}}(undef,N_t[i][t])
            for nd in 1:N_t[i][t]
            n += 1
            r[i][t][nd] = cells[n].r  
            end
        end 
    end

    return r
end
"""
    extract_c_paramter(cells,dataparams)
This function extracts the cell-level cluster probability vector for each cell in the CellFeatures object
"""
function extract_c_paramter(cells,dataparams)

    paramname = :_c
    I = dataparams.I
    T = dataparams.T
    N_t = dataparams.N_t
    c =  Vector{Vector{Vector{Vector{cells[1].BitType}}}}(undef,I)
    n = 0
    for i in 1:I
        c[i] = Vector{Vector{Vector{cells[1].BitType}}}(undef,T[i])
        for t in 1:T[i]
            c[i][t] = Vector{Vector{cells[1].BitType}}(undef,N_t[i][t])
            for nd in 1:N_t[i][t]
            n += 1
            c[i][t][nd] = cells[n].c  
            end
        end
    end

    return c
end

"""
    extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,training_logger)
This function extracts all parameters from custom objects and adds them to the previously instantiated output dictionary.
"""
function extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,conditionparams,dataparams,modelparams,training_logger)
    cluster_params_of_interest = [:y, :m_mu, :s_sq_mu, :g1, :g2,:h1, :h2,:u, :v, :x_hat,:x_hat_sq, :Nk,:a, :b, :m_nu, :s_sq_nu]
    condition_params_of_interest = [:d, :w1, :w2]
    gene_params_of_interest = []
    scalars_params_of_interest = []
    perK_elbos = [:elbo_]

    outputs_dict[:r_] = extract_r_paramter(cellpop, dataparams)
    for fn in cluster_params_of_interest
        key = Symbol(String(fn) * "_")
        outputs_dict[key] = extract_cluster_paramter(fn, clusters, modelparams)
    end
    outputs_dict[:c] = extract_c_paramter(cellpop, dataparams)
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
        outputs_dict[key] = extract_elbo_vals_perK(fn, training_logger)
    end
    # for fn in model_params
    #     key = Symbol(String(fn) * "_")
    #     outputs_dict[key] = getfield(modelparams,fn)
    # end
end

# """
#         extract_gene_paramter(paramname,geneparams,dataparams)
#     This function extracts the gene specific parameters from the GeneFeatures object used in inference
# """
# function extract_gene_paramter(paramname,geneparams,dataparams)
#     if typeof(paramname) <: String
#         paramname = Symbol(paramname)
#     end
#     G = dataparams.G

#     param = [ isone(length(getfield(geneparams[j],paramname))) ? getfield(geneparams[j],paramname)[1] : getfield(geneparams[j],paramname) for j in 1:G]
#     return param
# end

# """
#         extract_elbo_vals_perK(paramname,elbolog)
#     This function extracts the cluster specific elbo values from the ElboFeatures object used in inference
# """
# function extract_elbo_vals_perK(paramname,elbolog)
#     if typeof(paramname) <: String
#         paramname = Symbol(paramname)
#     end
#     # G = dataparams.G

#     vals = getfield(elbolog,paramname)
#     return vals
# end

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

# """
#         extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
#     This function extracts all parameters from custom objects and adds them to the previously instantiated output dictionary.
# """
# function extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
#     cluster_params_of_interest = [:yjk_hat, :mk_hat, :v_sq_k_hat, :σ_sq_k_hat, :var_muk,:κk_hat, :Nk, :gk_hat, :hk_hat, :ak_hat, :bk_hat, :x_hat,:x_hat_sq]
#     condition_params_of_interest = [:d_hat_t, :c_tt_prime, :st_hat]
#     gene_params_of_interest = [:λ_sq]

#     outputs_dict[:rtik_] = extract_rtik_paramter(cellpop, dataparams)
#     for fn in cluster_params_of_interest
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_cluster_paramter(fn, clusters, modelparams)
#     end
#     for fn in condition_params_of_interest
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_condition_paramter(fn, conditionparams, dataparams)
#     end
#     for fn in gene_params_of_interest
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_gene_paramter(fn, geneparams, dataparams)
#     end
# end

# """
#         extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog)
#     This function extracts all parameters from custom objects and adds them to the previously instantiated output dictionary.
# """
# function extract_and_add_parameters_to_outputs_dict!(outputs_dict,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams,elbolog)
#     cluster_params_of_interest = [:yjk_hat, :mk_hat, :v_sq_k_hat, :σ_sq_k_hat, :var_muk,:κk_hat, :Nk, :gk_hat, :hk_hat, :ak_hat, :bk_hat, :x_hat,:x_hat_sq]
#     condition_params_of_interest = [:d_hat_t, :c_tt_prime, :st_hat]
#     gene_params_of_interest = [:λ_sq]
#     perK_elbos = [:per_k_elbo]
#     model_params = [:ηk]

#     outputs_dict[:rtik_] = extract_rtik_paramter(cellpop, dataparams)
#     for fn in cluster_params_of_interest
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_cluster_paramter(fn, clusters, modelparams)
#     end
#     for fn in condition_params_of_interest
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_condition_paramter(fn, conditionparams, dataparams)
#     end
#     for fn in gene_params_of_interest
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_gene_paramter(fn, geneparams, dataparams)
#     end
#     for fn in perK_elbos
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = extract_elbo_vals_perK(fn, elbolog)
#     end
#     for fn in model_params
#         key = Symbol(String(fn) * "_")
#         outputs_dict[key] = getfield(modelparams,fn)
#     end
# end
