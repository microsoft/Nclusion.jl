"""
    getRandIndices(true_z, z_post_s)
This function takes in reference labels (true_z) and a vector of vectors that contains inferred labels (z_post_s) and computes the Adjusted Rand Index , Rand index,  Mirkin's index, and the Hubert's index.  
"""

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

"""
    getNMI(true_z, z_post_s)
This function takes in reference labels (true_z) and a vector of vectors that contains inferred labels (z_post_s) and computes the Normalized Mutual information between the reference and each inferred label sample. 
"""

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

"""
    getVmeasure(true_z, z_post_s; beta=1.0)
This function takes in reference labels (true_z) and a vector of vectors that contains inferred labels (z_post_s) and computes the V-measure  between the reference and each inferred label sample. 
"""

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

"""
    getVarInfo(true_z, z_post_s)
This function takes in reference labels (true_z) and a vector of vectors that contains inferred labels (z_post_s) and computes the Variation of information   between the reference and each inferred label sample. 
"""

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

"""
    getJaccardSimilarity(true_z, z_post_s)
This function takes in reference labels (true_z) and a vector of vectors that contains inferred labels (z_post_s) and computes the Jaccard Similarity between the reference and each inferred label sample. 
"""

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

"""
    calc_ARI_summarization(ari_vov, T; conf_level=0.95)
This function takes in calculates the mean Adjusted Rand Index and the upper and lower bounds of the Confidence Interval across inferred samples and for each conditions/timepoints
"""

function calc_ARI_summarization(ari_vov, T; conf_level=0.95)

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


"""
    calc_NMI_summarization(nmi_vov, T; conf_level=0.95)
This function takes in calculates the mean Normalized Mutual Information and the upper and lower bounds of the Confidence Interval across inferred samples and for each conditions/timepoints
"""

function calc_NMI_summarization(nmi_vov, T; conf_level=0.95)

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


"""
    calc_Vmeasure_summarization(vmeasure_vov, T; conf_level=0.95
This function takes in calculates the mean V-measure and the upper and lower bounds of the Confidence Interval across inferred samples and for each conditions/timepoints
"""

function calc_Vmeasure_summarization(vmeasure_vov, T; conf_level=0.95)
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

"""
    calc_VarInfo_summarization(varinfo_vov, T; conf_level=0.95)
This function takes in calculates the mean Variation of information and lower bounds of the Confidence Interval across inferred samples and for each conditions/timepoints
"""

function calc_VarInfo_summarization(varinfo_vov, T; conf_level=0.95)

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


"""
    calc_Jaccard_summarization(jaccard_vov, T; conf_level=0.95)
This function takes in calculates the mean Jaccard Similarity and lower bounds of the Confidence Interval across inferred samples and for each conditions/timepoints
"""

function calc_Jaccard_summarization(jaccard_vov, T; conf_level=0.95)

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

"""
    clustering_quality_metrics(results_df,filepath_;beta_ = 1.0, conf_lvl = 0.95)
This function calculates and saves all extrinsic metrics   
"""

function clustering_quality_metrics(results_df,filepath_;beta_ = 1.0, conf_lvl = 0.95, logger=nothing)
       #Get Called assignments
       called_assignments_vec = results_df[:,end-1]

       #Get Inferred assignments
       inferred_assignments_vec = results_df[:,end]
   
       #Tranform data into format needed for metric functions
       conditionpoints = sort(unique(results_df[:,1]))
       conditionpoints_dict =countmap(results_df[:,1])
       N_t = [conditionpoints_dict[el] for el in conditionpoints]
       ranges_ = get_timeranges(N_t)
       called_assignments = [called_assignments_vec[st:en] for (st,en) in ranges_]
       inferred_assignments = [inferred_assignments_vec[st:en] for (st,en) in ranges_]
       T = length(conditionpoints);
   
    #    @info "Calculating Metrics..."
        _flushed_logger("\t Calculating Metrics...";logger)

       #Get Metrics on a per time/condition basis
       ari_vov,_,_,_ = getRandIndices(called_assignments,[inferred_assignments])
       nmi_vov = getNMI(called_assignments,[inferred_assignments])
       vmeasure_vov = getVmeasure(called_assignments,[inferred_assignments];beta=beta_)
       varinfo_vov = getVarInfo(called_assignments,[inferred_assignments])
       jaccard_vov = getJaccardSimilarity(called_assignments,[inferred_assignments])
   
   
    #    @info "Saving Metrics..."
       _flushed_logger("\t Saving Metrics...";logger)
       #Save ARI as CSV
       ari_var_summary =calc_ARI_summarization(ari_vov,T;conf_level=conf_lvl)
       ari_var_summary_df = DataFrame(ari_var_summary,:auto);
       rename!(ari_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/ARI_Summary.csv", ari_var_summary_df) 
   
       #Save NMI as CSV
       nmi_var_summary = calc_NMI_summarization(nmi_vov,T;conf_level=conf_lvl)
       nmi_var_summary_df = DataFrame(nmi_var_summary,:auto);
       rename!(nmi_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/NMI_Summary.csv", nmi_var_summary_df) 
   
       #Save VMeasure as CSV
       vmeasure_var_summary =calc_Vmeasure_summarization(vmeasure_vov,T;conf_level=conf_lvl)
       vmeasure_var_summary_df = DataFrame(vmeasure_var_summary,:auto);
       rename!(vmeasure_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/VMeasure_Summary.csv", vmeasure_var_summary_df) 
   
       #Save VarInfo as CSV
       varinfo_var_summary = calc_VarInfo_summarization(varinfo_vov,T;conf_level=conf_lvl)
       varinfo_var_summary_df = DataFrame(varinfo_var_summary,:auto);
       rename!(varinfo_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/VarInfo_Summary.csv", varinfo_var_summary_df) 
   
       #Save Jaccard as CSV
       jaccard_var_summary = calc_Jaccard_summarization(jaccard_vov,T;conf_level=conf_lvl)
       jaccard_var_summary_df = DataFrame(jaccard_var_summary,:auto);
       rename!(jaccard_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
       CSV.write(filepath_*"/Jaccard_Summary.csv", jaccard_var_summary_df) 
   
      
end

