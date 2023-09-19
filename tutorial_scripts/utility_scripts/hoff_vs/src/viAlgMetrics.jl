


#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################
function tidy_calculate_elapsed_iterations(outputs_dict)
    elbo_ = outputs_dict[:elbo_]
    elsaped_iterations = length(elbo_)
    return elsaped_iterations
end
function tidy_calculate_approx_cluster_utilization(outputs_dict; thresh =1e-5 )
    rtik_mat = outputs_dict[:rtik_mat_]
    Nk_mat = tidy_update_Nk(rtik_mat);
    @views Nk_mat[Nk_mat[:,end] .<=thresh,end] .= 0.0
    return Nk_mat
end
function tidy_calculate_overall_cluster_occupancy_rate(outputs_dict; thresh =1e-5 )
    approx_cluster_utilization_mat = tidy_calculate_approx_cluster_utilization(outputs_dict; thresh =thresh )
    K = length(unique(approx_cluster_utilization_mat[:,1]))
    num_non_empty_clusters = sum(broadcast(!,iszero.(approx_cluster_utilization_mat[:,end])))
    overall_cluster_occupancy_rate = num_non_empty_clusters/K
    return overall_cluster_occupancy_rate,num_non_empty_clusters
end
function tidy_calculate_per_time_approx_cluster_utilization(outputs_dict; thresh =1e-5 )
    rtik_mat = outputs_dict[:rtik_mat_]
    Ntk_mat = tidy_update_Ntk(rtik_mat);
    @views Ntk_mat[Ntk_mat[:,end] .<=thresh,end] .= 0.0
    return Ntk_mat
end
function tidy_calculate_per_time_cluster_occupancy_rate(outputs_dict; thresh =1e-5 )
    per_time_approx_cluster_utilization_mat = tidy_calculate_per_time_approx_cluster_utilization(outputs_dict; thresh =thresh )
    T = length(unique(per_time_approx_cluster_utilization_mat[:,1]))
    K = length(unique(per_time_approx_cluster_utilization_mat[:,2]))
    num_non_empty_clusters = [sum(broadcast(!,iszero.(per_time_approx_cluster_utilization_mat[(t-1)*K+1:t*K,end]))) for t in 1:T]
    per_time_cluster_occupancy_rate = num_non_empty_clusters ./ K
    mean_time_cluster_occupancy_rate = mean(per_time_cluster_occupancy_rate)
    return mean_time_cluster_occupancy_rate,per_time_cluster_occupancy_rate,num_non_empty_clusters
end
