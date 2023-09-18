function getRandIndices(true_z, z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    ari_vov = Vector{Vector}(undef, T)
    ri_vov = Vector{Vector}(undef, T)
    mirkinindx_vov = Vector{Vector}(undef, T)
    hubertindx_vov = Vector{Vector}(undef, T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]), S)
        rand_indices = Clustering.randindex.(true_tp_assgn, pred_tp_assgn)
        ari_vov[t] = [indx[1] for indx in rand_indices]
        ri_vov[t] = [indx[2] for indx in rand_indices]
        mirkinindx_vov[t] = [indx[3] for indx in rand_indices]
        hubertindx_vov[t] = [indx[4] for indx in rand_indices]
    end
    return ari_vov, ri_vov, mirkinindx_vov, hubertindx_vov
end
function getNMI(true_z, z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    nmi_vov = Vector{Vector}(undef, T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]), S)
        nmi_vov[t] = Clustering.mutualinfo.(true_tp_assgn, pred_tp_assgn)
    end
    return nmi_vov
end
function getVmeasure(true_z, z_post_s; beta=1.0)
    T = length(true_z)
    S = length(z_post_s)
    vmeasure_vov = Vector{Vector}(undef, T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]), S)
        vmeasure_vov[t] = Clustering.vmeasure.(true_tp_assgn, pred_tp_assgn, β=beta)
    end
    return vmeasure_vov
end
function getVarInfo(true_z, z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    varinfo_vov = Vector{Vector}(undef, T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]), S)
        varinfo_vov[t] = Clustering.varinfo.(true_tp_assgn, pred_tp_assgn)
    end
    return varinfo_vov
end
function getJaccardSimilarity(true_z, z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    jaccard_vov = Vector{Vector}(undef, T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]), S)
        jaccard_vov[t] = 1 .- Distances.jaccard.(true_tp_assgn, pred_tp_assgn)
    end
    return jaccard_vov
end
function getMetrics_and_Stats(results_dict, algo_str, metric_str, KMax_vec, ylim_lb_vec, ylim_ub_vec; conf_level=0.95)
    T = results_dict["T"]
    algo_metric_mat = Array{Float64,3}(undef, T, 3, length(KMax_vec))
    algo_series_name = Vector{String}(undef, length(KMax_vec))
    for i in 1:length(KMax_vec)
        kmax = KMax_vec[i]
        metric_vov = results_dict["$(algo_str)Kmax$(kmax)Results"]["$(metric_str)_vov"]
        conf_int = t_test.(metric_vov; conf_level=conf_level)
        algo_metric_mat[:, 1, i] = mean.(metric_vov)
        lob = first.(conf_int)
        algo_metric_mat[:, 2, i] = first.(conf_int)
        hib = last.(conf_int)
        algo_metric_mat[:, 3, i] = last.(conf_int)
        # ylim_lb = maximum([0.0 - 0.01 ,minimum(lob)  - 0.1])
        push!(ylim_lb_vec, maximum([0.0 - 0.01, minimum(lob) - 0.1]))
        # ylim_ub = minimum([1.0 + 0.01 , maximum(hib) + 0.1])
        push!(ylim_ub_vec, minimum([1.0 + 0.01, maximum(hib) + 0.1]))
        algo_series_name[i] = "$(chop(algo_str)) (Kmax = $(kmax))"
    end
    return algo_metric_mat, algo_series_name, ylim_lb_vec, ylim_ub_vec
end

function time_invariant_ari(z, z_post_s)
    S = length(z_post_s)
    time_invariant_RandIndices = [[Clustering.randindex(Int.(recursive_flatten(z)), Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_ARI = [el[1] for el in time_invariant_RandIndices]
    return time_invariant_ARI
end
function time_invariant_nmi(z, z_post_s)
    S = length(z_post_s)
    time_invariant_NMI_ = [[Clustering.mutualinfo(Int.(recursive_flatten(z)), Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_NMI = [el[1] for el in time_invariant_NMI_]
    return time_invariant_NMI
end
function time_invariant_vmeasure(z, z_post_s; beta=1.0)
    S = length(z_post_s)
    time_invariant_Vmeasure_ = [[Clustering.vmeasure(Int.(recursive_flatten(z)), Int.(recursive_flatten(z_post_s[s])), β=beta)...] for s in 1:S]
    time_invariant_Vmeasure = [el[1] for el in time_invariant_Vmeasure_]
    return time_invariant_Vmeasure
end
function time_invariant_varinfo(z, z_post_s)
    S = length(z_post_s)
    time_invariant_VarInfo_ = [[Clustering.varinfo(Int.(recursive_flatten(z)), Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_VarInfo = [el[1] for el in time_invariant_VarInfo_]
    return time_invariant_VarInfo
end
function time_invariant_jaccard(z, z_post_s)
    S = length(z_post_s)
    time_invariant_Jaccard_ = [[1 .- Distances.jaccard.(Int.(recursive_flatten(z)), Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_Jaccard = [el[1] for el in time_invariant_Jaccard_]
    return time_invariant_Jaccard
end



function calc_time_invariant_ARI_summarization(ari_invar; conf_level=0.95)



    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    ari_RandIndices_values = ari_invar


    ari_RandIndices_means = mean(ari_RandIndices_values)
    ari_RandIndices_std = std(ari_RandIndices_values)
    if length(ari_RandIndices_values) > 1
        ari_RandIndices_conf_level = t_test(ari_RandIndices_values; conf_level=conf_level)
        ari_RandIndices_quantiles_mat = reduce(hcat, [ari_RandIndices_conf_level...])
    else
        ari_RandIndices_quantiles_mat = [ari_RandIndices_means, ari_RandIndices_means]
    end



    ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    ari_RandIndices_summary[1, 1] = 1
    ari_RandIndices_summary[1, 2] = ari_RandIndices_means
    ari_RandIndices_summary[1, 3] = ari_RandIndices_std
    ari_RandIndices_summary[1, 4] = ari_RandIndices_quantiles_mat[1]
    ari_RandIndices_summary[1, 5] = ari_RandIndices_quantiles_mat[2]


    return ari_RandIndices_summary
end
function calc_time_variant_ARI_summarization(ari_vov, T; conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    ari_RandIndices_values = ari_vov


    ari_RandIndices_means = mean.(ari_RandIndices_values)
    ari_RandIndices_std = std.(ari_RandIndices_values)
    if length(ari_RandIndices_values[1]) > 1
        ari_RandIndices_conf_level = t_test.(ari_RandIndices_values; conf_level=conf_level)
        ari_RandIndices_quantiles_mat = reduce(vcat, reduce.(hcat, [ari_RandIndices_conf_level...]))
    else
        ari_RandIndices_quantiles_mat = ari_RandIndices_means .* ones(Float64, T, 2)
    end

    ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    ari_RandIndices_summary[:, 1] = collect(1:T)
    ari_RandIndices_summary[:, 2] = ari_RandIndices_means
    ari_RandIndices_summary[:, 3] = ari_RandIndices_std
    ari_RandIndices_summary[:, 4] = ari_RandIndices_quantiles_mat[:, 1]
    ari_RandIndices_summary[:, 5] = ari_RandIndices_quantiles_mat[:, 2]


    return ari_RandIndices_summary
end
function calc_time_invariant_NMI_summarization(nmi_invar; conf_level=0.95)
    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    nmi_values = nmi_invar


    nmi_means = mean(nmi_values)
    nmi_std = std(nmi_values)
    if length(nmi_values) > 1
        nmi_conf_level = t_test(nmi_values; conf_level=conf_level)
        nmi_quantiles_mat = reduce(hcat, [nmi_conf_level...])
    else
        nmi_quantiles_mat = [nmi_means, nmi_means]
    end


    nmi_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    nmi_summary[1, 1] = 1
    nmi_summary[1, 2] = nmi_means
    nmi_summary[1, 3] = nmi_std
    nmi_summary[1, 4] = nmi_quantiles_mat[1]
    nmi_summary[1, 5] = nmi_quantiles_mat[2]


    return nmi_summary
end
function calc_time_variant_NMI_summarization(nmi_vov, T; conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    nmi_values = nmi_vov


    nmi_means = mean.(nmi_values)
    nmi_std = std.(nmi_values)
    if length(nmi_values[1]) > 1
        nmi_conf_level = t_test.(nmi_values; conf_level=conf_level)
        nmi_quantiles_mat = reduce(vcat, reduce.(hcat, [nmi_conf_level...]))
    else
        nmi_quantiles_mat = nmi_means .* ones(Float64, T, 2)
    end

    nmi_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    nmi_summary[:, 1] = collect(1:T)
    nmi_summary[:, 2] = nmi_means
    nmi_summary[:, 3] = nmi_std
    nmi_summary[:, 4] = nmi_quantiles_mat[:, 1]
    nmi_summary[:, 5] = nmi_quantiles_mat[:, 2]


    return nmi_summary
end
function calc_time_invariant_Vmeasure_summarization(vmeasure_invar; conf_level=0.95)



    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    vmeasure_values = vmeasure_invar


    vmeasure_means = mean(vmeasure_values)
    vmeasure_std = std(vmeasure_values)
    if length(vmeasure_values) > 1
        vmeasure_conf_level = t_test(vmeasure_values; conf_level=conf_level)
        vmeasure_quantiles_mat = reduce(hcat, [vmeasure_conf_level...])
    else
        vmeasure_quantiles_mat = [vmeasure_means, vmeasure_means]
    end

    vmeasure_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    vmeasure_summary[1, 1] = 1
    vmeasure_summary[1, 2] = vmeasure_means
    vmeasure_summary[1, 3] = vmeasure_std
    vmeasure_summary[1, 4] = vmeasure_quantiles_mat[1]
    vmeasure_summary[1, 5] = vmeasure_quantiles_mat[2]


    return vmeasure_summary
end
function calc_time_variant_Vmeasure_summarization(vmeasure_vov, T; conf_level=0.95)
    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    vmeasure_values = vmeasure_vov


    vmeasure_means = mean.(vmeasure_values)
    vmeasure_std = std.(vmeasure_values)
    if length(vmeasure_values[1]) > 1
        vmeasure_conf_level = t_test.(vmeasure_values; conf_level=conf_level)
        vmeasure_quantiles_mat = reduce(vcat, reduce.(hcat, [vmeasure_conf_level...]))
    else
        vmeasure_quantiles_mat = vmeasure_means .* ones(Float64, T, 2)
    end


    vmeasure_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    vmeasure_summary[:, 1] = collect(1:T)
    vmeasure_summary[:, 2] = vmeasure_means
    vmeasure_summary[:, 3] = vmeasure_std
    vmeasure_summary[:, 4] = vmeasure_quantiles_mat[:, 1]
    vmeasure_summary[:, 5] = vmeasure_quantiles_mat[:, 2]


    return vmeasure_summary
end
function calc_time_invariant_VarInfo_summarization(varinfo_invar; conf_level=0.95)
    varinfo_values = varinfo_invar


    varinfo_means = mean(varinfo_values)
    varinfo_std = std(varinfo_values)
    if length(varinfo_values) > 1
        varinfo_conf_level = t_test(varinfo_values; conf_level=conf_level)
        varinfo_quantiles_mat = reduce(hcat, [varinfo_conf_level...])
    else
        varinfo_quantiles_mat = [varinfo_means, varinfo_means]
    end

    varinfo_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    varinfo_summary[1, 1] = 1
    varinfo_summary[1, 2] = varinfo_means
    varinfo_summary[1, 3] = varinfo_std
    varinfo_summary[1, 4] = varinfo_quantiles_mat[1]
    varinfo_summary[1, 5] = varinfo_quantiles_mat[2]


    return varinfo_summary
end
function calc_time_variant_VarInfo_summarization(varinfo_vov, T; conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    varinfo_values = varinfo_vov


    varinfo_means = mean.(varinfo_values)
    varinfo_std = std.(varinfo_values)
    if length(varinfo_values[1]) > 1
        varinfo_conf_level = t_test.(varinfo_values; conf_level=conf_level)
        varinfo_quantiles_mat = reduce(vcat, reduce.(hcat, [varinfo_conf_level...]))
    else
        varinfo_quantiles_mat = varinfo_means .* ones(Float64, T, 2)
    end

    varinfo_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    varinfo_summary[:, 1] = collect(1:T)
    varinfo_summary[:, 2] = varinfo_means
    varinfo_summary[:, 3] = varinfo_std
    varinfo_summary[:, 4] = varinfo_quantiles_mat[:, 1]
    varinfo_summary[:, 5] = varinfo_quantiles_mat[:, 2]


    return varinfo_summary
end
function calc_time_invariant_Jaccard_summarization(jaccard_invar; conf_level=0.95)
    jaccard_values = jaccard_invar


    jaccard_means = mean(jaccard_values)
    jaccard_std = std(jaccard_values)
    if length(jaccard_values) > 1
        jaccard_conf_level = t_test(jaccard_values; conf_level=conf_level)
        jaccard_quantiles_mat = reduce(hcat, [jaccard_conf_level...])
    else
        jaccard_quantiles_mat = [jaccard_means, jaccard_means]
    end

    jaccard_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    jaccard_summary[1, 1] = 1
    jaccard_summary[1, 2] = jaccard_means
    jaccard_summary[1, 3] = jaccard_std
    jaccard_summary[1, 4] = jaccard_quantiles_mat[1]
    jaccard_summary[1, 5] = jaccard_quantiles_mat[2]


    return jaccard_summary
end
function calc_time_variant_Jaccard_summarization(jaccard_vov, T; conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    jaccard_values = jaccard_vov


    jaccard_means = mean.(jaccard_values)
    jaccard_std = std.(jaccard_values)
    if length(jaccard_values[1]) > 1
        jaccard_conf_level = t_test.(jaccard_values; conf_level=conf_level)
        jaccard_quantiles_mat = reduce(vcat, reduce.(hcat, [jaccard_conf_level...]))
    else
        jaccard_quantiles_mat = jaccard_means .* ones(Float64, T, 2)
    end

    jaccard_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    jaccard_summary[:, 1] = collect(1:T)
    jaccard_summary[:, 2] = jaccard_means
    jaccard_summary[:, 3] = jaccard_std
    jaccard_summary[:, 4] = jaccard_quantiles_mat[:, 1]
    jaccard_summary[:, 5] = jaccard_quantiles_mat[:, 2]


    return jaccard_summary
end
# function calc_time_invariant_ARI_summarization(ari_invar;conf_level=0.95)



#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     ari_RandIndices_values = ari_invar


#     ari_RandIndices_means = mean(ari_RandIndices_values)
#     ari_RandIndices_std = std(ari_RandIndices_values)
#     ari_RandIndices_conf_level = t_test(ari_RandIndices_values; conf_level=conf_level)
#     ari_RandIndices_quantiles_mat =reduce(hcat,[ari_RandIndices_conf_level...])


#     ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
#     ari_RandIndices_summary[1,1] = 1
#     ari_RandIndices_summary[1,2] = ari_RandIndices_means
#     ari_RandIndices_summary[1,3] = ari_RandIndices_std
#     ari_RandIndices_summary[1,4] = ari_RandIndices_quantiles_mat[1]
#     ari_RandIndices_summary[1,5] = ari_RandIndices_quantiles_mat[2]


#     return ari_RandIndices_summary
# end
# function calc_time_variant_ARI_summarization(ari_vov,T;conf_level=0.95)


#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     ari_RandIndices_values = ari_vov


#     ari_RandIndices_means = mean.(ari_RandIndices_values)
#     ari_RandIndices_std = std.(ari_RandIndices_values)
#     ari_RandIndices_conf_level = t_test.(ari_RandIndices_values; conf_level=conf_level)
#     ari_RandIndices_quantiles_mat =reduce(vcat,reduce.(hcat,[ari_RandIndices_conf_level...]))


#     ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
#     ari_RandIndices_summary[:,1] = collect(1:T)
#     ari_RandIndices_summary[:,2] = ari_RandIndices_means
#     ari_RandIndices_summary[:,3] = ari_RandIndices_std
#     ari_RandIndices_summary[:,4] = ari_RandIndices_quantiles_mat[:,1]
#     ari_RandIndices_summary[:,5] = ari_RandIndices_quantiles_mat[:,2]


#     return ari_RandIndices_summary
# end
# function calc_time_invariant_NMI_summarization(nmi_invar;conf_level=0.95)
#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     nmi_values = nmi_invar


#     nmi_means = mean(nmi_values)
#     nmi_std = std(nmi_values)
#     nmi_conf_level = t_test(nmi_values; conf_level=conf_level)
#     nmi_quantiles_mat =reduce(hcat,[nmi_conf_level...])


#     nmi_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
#     nmi_summary[1,1] = 1
#     nmi_summary[1,2] = nmi_means
#     nmi_summary[1,3] = nmi_std
#     nmi_summary[1,4] = nmi_quantiles_mat[1]
#     nmi_summary[1,5] = nmi_quantiles_mat[2]


#     return nmi_summary
# end
# function calc_time_variant_NMI_summarization(nmi_vov,T;conf_level=0.95)


#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     nmi_values = nmi_vov


#     nmi_means = mean.(nmi_values)
#     nmi_std = std.(nmi_values)
#     nmi_conf_level = t_test.(nmi_values; conf_level=conf_level)
#     nmi_quantiles_mat = reduce(vcat,reduce.(hcat,[nmi_conf_level...])) 


#     nmi_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
#     nmi_summary[:,1] = collect(1:T)
#     nmi_summary[:,2] = nmi_means
#     nmi_summary[:,3] = nmi_std
#     nmi_summary[:,4] = nmi_quantiles_mat[:,1]
#     nmi_summary[:,5] = nmi_quantiles_mat[:,2]


#     return nmi_summary
# end
# function calc_time_invariant_Vmeasure_summarization(vmeasure_invar;conf_level=0.95)



#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     vmeasure_values = vmeasure_invar


#     vmeasure_means = mean(vmeasure_values)
#     vmeasure_std = std(vmeasure_values)
#     vmeasure_conf_level = t_test(vmeasure_values; conf_level=conf_level)
#     vmeasure_quantiles_mat =reduce(hcat,[vmeasure_conf_level...])


#     vmeasure_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
#     vmeasure_summary[1,1] = 1
#     vmeasure_summary[1,2] = vmeasure_means
#     vmeasure_summary[1,3] = vmeasure_std
#     vmeasure_summary[1,4] = vmeasure_quantiles_mat[1]
#     vmeasure_summary[1,5] = vmeasure_quantiles_mat[2]


#     return vmeasure_summary
# end
# function calc_time_variant_Vmeasure_summarization(vmeasure_vov,T;conf_level=0.95)
#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     vmeasure_values = vmeasure_vov


#     vmeasure_means = mean.(vmeasure_values)
#     vmeasure_std = std.(vmeasure_values)
#     vmeasure_conf_level = t_test.(vmeasure_values; conf_level=conf_level)
#     vmeasure_quantiles_mat = reduce(vcat,reduce.(hcat,[vmeasure_conf_level...])) 


#     vmeasure_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
#     vmeasure_summary[:,1] = collect(1:T)
#     vmeasure_summary[:,2] = vmeasure_means
#     vmeasure_summary[:,3] = vmeasure_std
#     vmeasure_summary[:,4] = vmeasure_quantiles_mat[:,1]
#     vmeasure_summary[:,5] = vmeasure_quantiles_mat[:,2]


#     return vmeasure_summary
# end
# function calc_time_invariant_VarInfo_summarization(varinfo_invar;conf_level=0.95)
#     varinfo_values = varinfo_invar


#     varinfo_means = mean(varinfo_values)
#     varinfo_std = std(varinfo_values)
#     varinfo_conf_level = t_test(varinfo_values; conf_level=conf_level)
#     varinfo_quantiles_mat =reduce(hcat,[varinfo_conf_level...])


#     varinfo_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
#     varinfo_summary[1,1] = 1
#     varinfo_summary[1,2] = varinfo_means
#     varinfo_summary[1,3] = varinfo_std
#     varinfo_summary[1,4] = varinfo_quantiles_mat[1]
#     varinfo_summary[1,5] = varinfo_quantiles_mat[2]


#     return varinfo_summary
# end
# function calc_time_variant_VarInfo_summarization(varinfo_vov,T;conf_level=0.95)


#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     varinfo_values = varinfo_vov


#     varinfo_means = mean.(varinfo_values)
#     varinfo_std = std.(varinfo_values)
#     varinfo_conf_level = t_test.(varinfo_values; conf_level=conf_level)
#     varinfo_quantiles_mat = reduce(vcat,reduce.(hcat,[varinfo_conf_level...])) 


#     varinfo_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
#     varinfo_summary[:,1] = collect(1:T)
#     varinfo_summary[:,2] = varinfo_means
#     varinfo_summary[:,3] = varinfo_std
#     varinfo_summary[:,4] = varinfo_quantiles_mat[:,1]
#     varinfo_summary[:,5] = varinfo_quantiles_mat[:,2]


#     return varinfo_summary
# end
# function calc_time_invariant_Jaccard_summarization(jaccard_invar;conf_level=0.95)
#     jaccard_values = jaccard_invar


#     jaccard_means = mean(jaccard_values)
#     jaccard_std = std(jaccard_values)
#     jaccard_conf_level = t_test(jaccard_values; conf_level=conf_level)
#     jaccard_quantiles_mat =reduce(hcat,[jaccard_conf_level...])


#     jaccard_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
#     jaccard_summary[1,1] = 1
#     jaccard_summary[1,2] = jaccard_means
#     jaccard_summary[1,3] = jaccard_std
#     jaccard_summary[1,4] = jaccard_quantiles_mat[1]
#     jaccard_summary[1,5] = jaccard_quantiles_mat[2]


#     return jaccard_summary
# end
# function calc_time_variant_Jaccard_summarization(jaccard_vov,T;conf_level=0.95)


#     #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
#     jaccard_values = jaccard_vov


#     jaccard_means = mean.(jaccard_values)
#     jaccard_std = std.(jaccard_values)
#     jaccard_conf_level = t_test.(jaccard_values; conf_level=conf_level)
#     jaccard_quantiles_mat = reduce(vcat,reduce.(hcat,[jaccard_conf_level...])) 


#     jaccard_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
#     jaccard_summary[:,1] = collect(1:T)
#     jaccard_summary[:,2] = jaccard_means
#     jaccard_summary[:,3] = jaccard_std
#     jaccard_summary[:,4] = jaccard_quantiles_mat[:,1]
#     jaccard_summary[:,5] = jaccard_quantiles_mat[:,2]


#     return jaccard_summary
# end
function calc_time_invariant_CVI_summarization(cvi; conf_level=0.95)



    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    cvi_values = cvi


    cvi_means = mean(cvi_values)
    cvi_std = std(cvi_values)
    # cvi_conf_level = t_test(cvi_values; conf_level=conf_level)
    # cvi_quantiles_mat = reduce(hcat, [cvi_conf_level...])
    if length(cvi_values) > 1
        cvi_conf_level = t_test(cvi_values; conf_level=conf_level)
        cvi_quantiles_mat = reduce(hcat, [cvi_conf_level...])
    else
        cvi_quantiles_mat = [cvi_means, cvi_means]
    end

    cvi_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    cvi_summary[1, 1] = 1
    cvi_summary[1, 2] = cvi_means
    cvi_summary[1, 3] = cvi_std
    cvi_summary[1, 4] = cvi_quantiles_mat[1]
    cvi_summary[1, 5] = cvi_quantiles_mat[2]


    return cvi_summary
end

function setup_metrics_list()
    return ["ari", "nmi", "vmeasure", "varinfo", "jaccard", "ch", "csil", "db", "gd43", "gd53", "ps", "rcip", "wb", "xb", "final_elbo"]
end
# function time_invariant_fscore(z,z_post_s)
#     S = length(z_post_s)
#     time_invariant_fscore = []
#     for s in 1:S
#         pred = Int.(recursive_flatten(z_post_s[s]))
#         labels = Int.(recursive_flatten(z))
#         @rput pred;
#         @rput labels;
#         R"""
#             library(PerfMeas)
#             library(yardstick)
#             fmeasure = f_meas_vec(factor(labels), factor(pred))
#         """
#         @rget fmeasure;
#         push!(time_invariant_fscore,fmeasure)
#     end
#     # [[MLJBase.FScore(Int.(recursive_flatten(z_post_s[s])),Int.(recursive_flatten(z)))...] for s in 1:S]
#     time_invariant_fscore = [el[1] for el in time_invariant_fscore]
#     return time_invariant_fscore
# end


#####################################################
#####################################################
################# TIDY FUNCTIONS ####################
#####################################################
#####################################################
function tidy_getRandIndices(zmat, z_samples_mat)
    S = length(unique(z_samples_mat[:, 1]))
    T = length(unique(z_samples_mat[:, 2]))
    N = length(unique(z_samples_mat[:, 3]))

    sample_ids = collect(1:S)
    time_invariant_RandIndices = [[Clustering.randindex(Int.(zmat[:, end]), Int.(z_samples_mat[(s-1)*N+1:s*N, end]))...] for s in 1:S]
    # time_invariant_RandIndices_mat = permutedims(reduce(hcat,time_invariant_RandIndices))
    # time_invariant_RandIndices_mat = reduce(hcat,[sample_ids,time_invariant_RandIndices_mat])
    time_invariant_RandIndices_mat = Matrix{Union{Float64,Int}}(undef, S, 5)
    time_invariant_RandIndices_mat[:, 1] = sample_ids
    time_invariant_RandIndices_mat[:, 2:end] = permutedims(reduce(hcat, time_invariant_RandIndices))


    timepoint_freq = countmap(Int.(zmat[:, 1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
    # println(timeranges)
    # println(N_t)
    if T >= 2
        time_variant_RandIndices = [[Clustering.randindex(Int.(zmat[st:en, end]), Int.(z_samples_mat[(s-1)*N+st:(s-1)*N+en, end]))...] for s in 1:S for (st, en) in timeranges]
    else
        time_variant_RandIndices = [[Clustering.randindex(Int.(zmat[1:end, end]), Int.(z_samples_mat[(s-1)*N+1:(s)*N, end]))...] for s in 1:S]
    end

    # time_variant_RandIndices_mat = permutedims(reduce(hcat,time_variant_RandIndices))
    # time_variant_RandIndices_mat = reduce(hcat,[innermelt(sample_ids,T),time_variant_RandIndices_mat])
    time_variant_RandIndices_mat = Matrix{Union{Float64,Int}}(undef, T * S, 5)
    time_variant_RandIndices_mat[:, 1] = innermelt(sample_ids, T)
    time_variant_RandIndices_mat[:, 2:end] = permutedims(reduce(hcat, time_variant_RandIndices))
    # [(s-1)*N+st:(s-1)*N+en for s in 1:S for (st,en) in timeranges]

    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    return time_invariant_RandIndices_mat, time_variant_RandIndices_mat
end
function tidy_calc_time_invariant_ARI_summarization(time_invariant_RandIndices_mat; conf_level=0.95)
    S = length(unique(time_invariant_RandIndices_mat[:, 1]))


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    ari_RandIndices_values = time_invariant_RandIndices_mat[:, 2]


    ari_RandIndices_means = mean(ari_RandIndices_values)
    ari_RandIndices_std = std(ari_RandIndices_values)
    ari_RandIndices_conf_level = t_test(ari_RandIndices_values; conf_level=conf_level)
    ari_RandIndices_quantiles_mat = reduce(hcat, [ari_RandIndices_conf_level...])


    ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    ari_RandIndices_summary[1, 1] = 1
    ari_RandIndices_summary[1, 2] = ari_RandIndices_means
    ari_RandIndices_summary[1, 3] = ari_RandIndices_std
    ari_RandIndices_summary[1, 4] = ari_RandIndices_quantiles_mat[1]
    ari_RandIndices_summary[1, 5] = ari_RandIndices_quantiles_mat[2]


    return ari_RandIndices_summary
end
function tidy_calc_time_variant_ARI_summarization(time_variant_RandIndices_mat; conf_level=0.95)
    S = length(unique(time_variant_RandIndices_mat[:, 1]))
    sample_freq = countmap(Int.(time_variant_RandIndices_mat[:, 1]))
    T = maximum([sample_freq[key] for key in sort(collect(keys(sample_freq)))])

    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    ari_RandIndices_values = [time_variant_RandIndices_mat[t:T:end, 2] for t in 1:T]


    ari_RandIndices_means = mean.(ari_RandIndices_values)
    ari_RandIndices_std = std.(ari_RandIndices_values)
    ari_RandIndices_conf_level = t_test.(ari_RandIndices_values; conf_level=conf_level)
    ari_RandIndices_quantiles_mat = reduce(vcat, reduce.(hcat, [ari_RandIndices_conf_level...]))


    ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    ari_RandIndices_summary[:, 1] = collect(1:T)
    ari_RandIndices_summary[:, 2] = ari_RandIndices_means
    ari_RandIndices_summary[:, 3] = ari_RandIndices_std
    ari_RandIndices_summary[:, 4] = ari_RandIndices_quantiles_mat[:, 1]
    ari_RandIndices_summary[:, 5] = ari_RandIndices_quantiles_mat[:, 2]


    return ari_RandIndices_summary
end
function tidy_getNMI(zmat, z_samples_mat)
    S = length(unique(z_samples_mat[:, 1]))
    T = length(unique(z_samples_mat[:, 2]))
    N = length(unique(z_samples_mat[:, 3]))
    sample_ids = collect(1:S)
    time_invariant_NMI = [[Clustering.mutualinfo(Int.(zmat[:, end]), Int.(z_samples_mat[(s-1)*N+1:s*N, end]))...] for s in 1:S]
    time_invariant_NMI_mat = Matrix{Union{Float64,Int}}(undef, S, 2)
    time_invariant_NMI_mat[:, 1] = sample_ids
    time_invariant_NMI_mat[:, 2:end] = permutedims(reduce(hcat, time_invariant_NMI))

    timepoint_freq = countmap(Int.(zmat[:, 1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)

    time_variant_NMI = [[Clustering.mutualinfo(Int.(zmat[st:en, end]), Int.(z_samples_mat[(s-1)*N+st:(s-1)*N+en, end]))...] for s in 1:S for (st, en) in timeranges]
    time_variant_NMI_mat = Matrix{Union{Float64,Int}}(undef, T * S, 2)
    time_variant_NMI_mat[:, 1] = innermelt(sample_ids, T)
    time_variant_NMI_mat[:, 2:end] = permutedims(reduce(hcat, time_variant_NMI))

    return time_invariant_NMI_mat, time_variant_NMI_mat
end
function tidy_calc_time_invariant_NMI_summarization(time_invariant_NMI_mat; conf_level=0.95)
    S = length(unique(time_invariant_NMI_mat[:, 1]))


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    nmi_values = time_invariant_NMI_mat[:, 2]


    nmi_means = mean(nmi_values)
    nmi_std = std(nmi_values)
    nmi_conf_level = t_test(nmi_values; conf_level=conf_level)
    nmi_quantiles_mat = reduce(hcat, [nmi_conf_level...])


    nmi_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    nmi_summary[1, 1] = 1
    nmi_summary[1, 2] = nmi_means
    nmi_summary[1, 3] = nmi_std
    nmi_summary[1, 4] = nmi_quantiles_mat[1]
    nmi_summary[1, 5] = nmi_quantiles_mat[2]


    return nmi_summary
end
function tidy_calc_time_variant_NMI_summarization(time_variant_NMI_mat; conf_level=0.95)
    S = length(unique(time_variant_NMI_mat[:, 1]))
    sample_freq = countmap(Int.(time_variant_NMI_mat[:, 1]))
    T = maximum([sample_freq[key] for key in sort(collect(keys(sample_freq)))])


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    nmi_values = [time_variant_NMI_mat[t:T:end, 2] for t in 1:T]


    nmi_means = mean.(nmi_values)
    nmi_std = std.(nmi_values)
    nmi_conf_level = t_test.(nmi_values; conf_level=conf_level)
    nmi_quantiles_mat = reduce(vcat, reduce.(hcat, [nmi_conf_level...]))


    nmi_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    nmi_summary[:, 1] = collect(1:T)
    nmi_summary[:, 2] = nmi_means
    nmi_summary[:, 3] = nmi_std
    nmi_summary[:, 4] = nmi_quantiles_mat[:, 1]
    nmi_summary[:, 5] = nmi_quantiles_mat[:, 2]


    return nmi_summary
end
function tidy_getVmeasure(zmat, z_samples_mat; beta=1.0)
    S = length(unique(z_samples_mat[:, 1]))
    T = length(unique(z_samples_mat[:, 2]))
    N = length(unique(z_samples_mat[:, 3]))
    sample_ids = collect(1:S)
    time_invariant_Vmeasure = [[Clustering.vmeasure(Int.(zmat[:, end]), Int.(z_samples_mat[(s-1)*N+1:s*N, end]), β=beta)...] for s in 1:S]
    time_invariant_Vmeasure_mat = Matrix{Union{Float64,Int}}(undef, S, 2)
    time_invariant_Vmeasure_mat[:, 1] = sample_ids
    time_invariant_Vmeasure_mat[:, 2:end] = permutedims(reduce(hcat, time_invariant_Vmeasure))

    timepoint_freq = countmap(Int.(zmat[:, 1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)

    time_variant_Vmeasure = [[Clustering.vmeasure(Int.(zmat[st:en, end]), Int.(z_samples_mat[(s-1)*N+st:(s-1)*N+en, end]), β=beta)...] for s in 1:S for (st, en) in timeranges]
    time_variant_Vmeasure_mat = Matrix{Union{Float64,Int}}(undef, T * S, 2)
    time_variant_Vmeasure_mat[:, 1] = innermelt(sample_ids, T)
    time_variant_Vmeasure_mat[:, 2:end] = permutedims(reduce(hcat, time_variant_Vmeasure))

    return time_invariant_Vmeasure_mat, time_variant_Vmeasure_mat
end
function tidy_calc_time_invariant_Vmeasure_summarization(time_invariant_Vmeasure_mat; conf_level=0.95)
    S = length(unique(time_invariant_Vmeasure_mat[:, 1]))


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    vmeasure_values = time_invariant_Vmeasure_mat[:, 2]


    vmeasure_means = mean(vmeasure_values)
    vmeasure_std = std(vmeasure_values)
    vmeasure_conf_level = t_test(vmeasure_values; conf_level=conf_level)
    vmeasure_quantiles_mat = reduce(hcat, [vmeasure_conf_level...])


    vmeasure_summary = Matrix{Union{Float64,Int}}(undef, 1, 5)
    vmeasure_summary[1, 1] = 1
    vmeasure_summary[1, 2] = vmeasure_means
    vmeasure_summary[1, 3] = vmeasure_std
    vmeasure_summary[1, 4] = vmeasure_quantiles_mat[1]
    vmeasure_summary[1, 5] = vmeasure_quantiles_mat[2]


    return vmeasure_summary
end
function tidy_calc_time_variant_Vmeasure_summarization(time_variant_Vmeasure_mat; conf_level=0.95)
    S = length(unique(time_variant_Vmeasure_mat[:, 1]))
    sample_freq = countmap(Int.(time_variant_Vmeasure_mat[:, 1]))
    T = maximum([sample_freq[key] for key in sort(collect(keys(sample_freq)))])

    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    vmeasure_values = [time_variant_Vmeasure_mat[t:T:end, 2] for t in 1:T]


    vmeasure_means = mean.(vmeasure_values)
    vmeasure_std = std.(vmeasure_values)
    vmeasure_conf_level = t_test.(vmeasure_values; conf_level=conf_level)
    vmeasure_quantiles_mat = reduce(vcat, reduce.(hcat, [vmeasure_conf_level...]))


    vmeasure_summary = Matrix{Union{Float64,Int}}(undef, T, 5)
    vmeasure_summary[:, 1] = collect(1:T)
    vmeasure_summary[:, 2] = vmeasure_means
    vmeasure_summary[:, 3] = vmeasure_std
    vmeasure_summary[:, 4] = vmeasure_quantiles_mat[:, 1]
    vmeasure_summary[:, 5] = vmeasure_quantiles_mat[:, 2]


    return vmeasure_summary
end
function tidy_getVarInfo(zmat, z_samples_mat)
    S = length(unique(z_samples_mat[:, 1]))
    T = length(unique(z_samples_mat[:, 2]))
    N = length(unique(z_samples_mat[:, 3]))
    sample_ids = collect(1:S)
    time_invariant_VarInfo = [[Clustering.varinfo(Int.(zmat[:, end]), Int.(z_samples_mat[(s-1)*N+1:s*N, end]))...] for s in 1:S]
    time_invariant_VarInfo_mat = Matrix{Union{Float64,Int}}(undef, S, 2)
    time_invariant_VarInfo_mat[:, 1] = sample_ids
    time_invariant_VarInfo_mat[:, 2:end] = permutedims(reduce(hcat, time_invariant_VarInfo))

    timepoint_freq = countmap(Int.(zmat[:, 1]))
    N_t = [timepoint_freq[key] for key in sort(collect(keys(timepoint_freq)))]
    timeranges = zip(collect(0:T-1) .* N_t .+ 1, collect(1:T) .* N_t)

    time_variant_VarInfo = [[Clustering.varinfo(Int.(zmat[st:en, end]), Int.(z_samples_mat[(s-1)*N+st:(s-1)*N+en, end]))...] for s in 1:S for (st, en) in timeranges]
    time_variant_VarInfo_mat = Matrix{Union{Float64,Int}}(undef, T * S, 2)
    time_variant_VarInfo_mat[:, 1] = innermelt(sample_ids, T)
    time_variant_VarInfo_mat[:, 2:end] = permutedims(reduce(hcat, time_variant_VarInfo))

    return time_invariant_VarInfo_mat, time_variant_VarInfo_mat
end
