println("Loading Packages...")
using Clustering # TSne, MultivariateStats,
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
using Distances
using CSV,DataFrames
using Distributions



function main(ARGS)
    println("Parsing Inputs...")
    #Parsing Arguments
    cluster_results_filename, confidence_level, vmeasure_beta = ARGS


    #If Confidence level is not provided set default to 0.95, otherwise take user input
    if !isempty(confidence_level)
        conf_lvl = parse(Float64, confidence_level)
    else
        conf_lvl = 0.95
    end

    #If VMeasure Beta is not provided set default to 1.0, otherwise take user input
    if !isempty(vmeasure_beta)
        beta = parse(Float64, vmeasure_beta)
    else
        beta = 1.0
    end

    #Get filepath
    filepath = join(split(cluster_results_filename,"/")[1:end-1],"/")
    
    println("Loading Data...")
    #Load Data from filepath
    results_df =  CSV.read(cluster_results_filename,DataFrame);

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

    println("Calculating Metrics...")
    #Get Metrics on a per time/condition basis
    ari_vov,_,_,_ = getRandIndices(called_assignments,[inferred_assignments])
    nmi_vov = getNMI(called_assignments,[inferred_assignments])
    vmeasure_vov = getVmeasure(called_assignments,[inferred_assignments];beta=beta)
    varinfo_vov = getVarInfo(called_assignments,[inferred_assignments])
    jaccard_vov = getJaccardSimilarity(called_assignments,[inferred_assignments])

    #Get overall Metric
    ari_invar = time_invariant_ari(called_assignments,[inferred_assignments])
    nmi_invar =time_invariant_nmi(called_assignments,[inferred_assignments])
    vmeasure_invar =time_invariant_vmeasure(called_assignments,[inferred_assignments],;beta=beta)
    varinfo_invar =time_invariant_varinfo(called_assignments,[inferred_assignments])
    jaccard_invar =time_invariant_jaccard(called_assignments,[inferred_assignments])

    println("Saving Metrics...")
    #Save ARI as CSV
    ari_invar_summary = calc_time_invariant_ARI_summarization(ari_invar;conf_level=conf_lvl);
    ari_invar_summary_df = DataFrame(ari_invar_summary,:auto);
    rename!(ari_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/ARI_Summary_AllTimepoints.csv", ari_invar_summary_df) 
    ari_var_summary =calc_time_variant_ARI_summarization(ari_vov,T;conf_level=conf_lvl)
    ari_var_summary_df = DataFrame(ari_var_summary,:auto);
    rename!(ari_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/ARI_Summary_PerTimepoints.csv", ari_var_summary_df) 

    #Save NMI as CSV
    nmi_invar_summary = calc_time_invariant_NMI_summarization(nmi_invar;conf_level=conf_lvl)
    nmi_invar_summary_df = DataFrame(nmi_invar_summary,:auto);
    rename!(nmi_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/NMI_Summary_AllTimepoints.csv", nmi_invar_summary_df) 
    nmi_var_summary = calc_time_variant_NMI_summarization(nmi_vov,T;conf_level=conf_lvl)
    nmi_var_summary_df = DataFrame(nmi_var_summary,:auto);
    rename!(nmi_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/NMI_Summary_PerTimepoints.csv", nmi_var_summary_df) 

    #Save VMeasure as CSV
    vmeasure_invar_summary = calc_time_invariant_Vmeasure_summarization(vmeasure_invar;conf_level=conf_lvl)
    vmeasure_invar_summary_df = DataFrame(vmeasure_invar_summary,:auto);
    rename!(vmeasure_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/VMeasure_Summary_AllTimepoints.csv", vmeasure_invar_summary_df) 
    vmeasure_var_summary =calc_time_variant_Vmeasure_summarization(vmeasure_vov,T;conf_level=conf_lvl)
    vmeasure_var_summary_df = DataFrame(vmeasure_var_summary,:auto);
    rename!(vmeasure_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/VMeasure_Summary_PerTimepoints.csv", vmeasure_var_summary_df) 

    #Save VarInfo as CSV
    varinfo_invar_summary = calc_time_invariant_VarInfo_summarization(varinfo_invar;conf_level=conf_lvl)
    varinfo_invar_summary_df = DataFrame(varinfo_invar_summary,:auto);
    rename!(varinfo_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/VarInfo_Summary_AllTimepoints.csv", varinfo_invar_summary_df) 
    varinfo_var_summary = calc_time_variant_VarInfo_summarization(varinfo_vov,T;conf_level=conf_lvl)
    varinfo_var_summary_df = DataFrame(varinfo_var_summary,:auto);
    rename!(varinfo_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/VarInfo_Summary_PerTimepoints.csv", varinfo_var_summary_df) 

    #Save Jaccard as CSV
    jaccard_invar_summary = calc_time_invariant_Jaccard_summarization(jaccard_invar;conf_level=conf_lvl)
    jaccard_invar_summary_df = DataFrame(jaccard_invar_summary,:auto);
    rename!(jaccard_invar_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/Jaccard_Summary_AllTimepoints.csv", jaccard_invar_summary_df)
    jaccard_var_summary = calc_time_variant_Jaccard_summarization(jaccard_vov,T;conf_level=conf_lvl)
    jaccard_var_summary_df = DataFrame(jaccard_var_summary,:auto);
    rename!(jaccard_var_summary_df,[:condition,:mean,:sd,:lowerbound,:upperbound])
    CSV.write(filepath*"/Jaccard_Summary_PerTimepoints.csv", jaccard_var_summary_df) 

    println("Finished, Ending Script...")
end

function t_test(x; conf_level=0.95)
    alpha = (1 - conf_level)
    tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
    SE = std(x)./sqrt(length(x))

    lo, hi = mean(x) .+ [-1, 1] .* tstar * SE
    # "($lo, $hi)"
    return lo, hi
end
function recursive_flatten(x::AbstractArray)
    if any(a->typeof(a)<:AbstractArray, x)#eltype(x) <: Vector
        recursive_flatten(vcat(x...))
    else
        return x
    end
end

function getRandIndices(true_z,z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    ari_vov = Vector{Vector}(undef,T)
    ri_vov = Vector{Vector}(undef,T)
    mirkinindx_vov = Vector{Vector}(undef,T)
    hubertindx_vov = Vector{Vector}(undef,T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]),S)
        rand_indices = Clustering.randindex.(true_tp_assgn, pred_tp_assgn)
        ari_vov[t] = [indx[1] for indx in rand_indices]
        ri_vov[t] = [indx[2] for indx in rand_indices]
        mirkinindx_vov[t] = [indx[3] for indx in rand_indices]
        hubertindx_vov[t] =  [indx[4] for indx in rand_indices]
    end
    return ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
end
function getNMI(true_z,z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    nmi_vov = Vector{Vector}(undef,T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]),S)
        nmi_vov[t] = Clustering.mutualinfo.(true_tp_assgn, pred_tp_assgn)
    end
    return nmi_vov
end
function getVmeasure(true_z,z_post_s;beta=1.0)
    T = length(true_z)
    S = length(z_post_s)
    vmeasure_vov = Vector{Vector}(undef,T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]),S)
        vmeasure_vov[t] = Clustering.vmeasure.(true_tp_assgn, pred_tp_assgn, β = beta)
    end
    return vmeasure_vov
end
function getVarInfo(true_z,z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    varinfo_vov = Vector{Vector}(undef,T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]),S)
        varinfo_vov[t] = Clustering.varinfo.(true_tp_assgn, pred_tp_assgn)
    end
    return varinfo_vov
end
function getJaccardSimilarity(true_z,z_post_s)
    T = length(true_z)
    S = length(z_post_s)
    jaccard_vov = Vector{Vector}(undef,T)
    for t in 1:T
        pred_tp_assgn = [el[t] for el in z_post_s]
        true_tp_assgn = fill(Int.(true_z[t]),S)
        jaccard_vov[t] = 1 .- Distances.jaccard.(true_tp_assgn, pred_tp_assgn)
    end
    return jaccard_vov
end

function time_invariant_ari(z,z_post_s)
    S = length(z_post_s)
    time_invariant_RandIndices = [[Clustering.randindex(Int.(recursive_flatten(z)),Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_ARI = [el[1] for el in time_invariant_RandIndices]
    return time_invariant_ARI
end
function time_invariant_nmi(z,z_post_s)
    S = length(z_post_s)
    time_invariant_NMI_ = [[Clustering.mutualinfo(Int.(recursive_flatten(z)),Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_NMI = [el[1] for el in time_invariant_NMI_]
    return time_invariant_NMI
end
function time_invariant_vmeasure(z,z_post_s,;beta=1.0)
    S = length(z_post_s)
    time_invariant_Vmeasure_ = [[Clustering.vmeasure(Int.(recursive_flatten(z)),Int.(recursive_flatten(z_post_s[s])), β = beta)...] for s in 1:S]
    time_invariant_Vmeasure = [el[1] for el in time_invariant_Vmeasure_]
    return time_invariant_Vmeasure
end
function time_invariant_varinfo(z,z_post_s)
    S = length(z_post_s)
    time_invariant_VarInfo_ = [[Clustering.varinfo(Int.(recursive_flatten(z)),Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_VarInfo = [el[1] for el in time_invariant_VarInfo_]
    return time_invariant_VarInfo
end
function time_invariant_jaccard(z,z_post_s)
    S = length(z_post_s)
    time_invariant_Jaccard_ = [[1 .- Distances.jaccard.(Int.(recursive_flatten(z)),Int.(recursive_flatten(z_post_s[s])))...] for s in 1:S]
    time_invariant_Jaccard = [el[1] for el in time_invariant_Jaccard_]
    return time_invariant_Jaccard
end


function calc_time_invariant_ARI_summarization(ari_invar;conf_level=0.95)



    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    ari_RandIndices_values = ari_invar
    
 
    ari_RandIndices_means = mean(ari_RandIndices_values)
    ari_RandIndices_std = std(ari_RandIndices_values)
    if length(ari_RandIndices_values) > 1
        ari_RandIndices_conf_level = t_test(ari_RandIndices_values; conf_level=conf_level)
        ari_RandIndices_quantiles_mat =reduce(hcat,[ari_RandIndices_conf_level...])
    else
        ari_RandIndices_quantiles_mat = [ari_RandIndices_means,ari_RandIndices_means]
    end
    

    
    ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
    ari_RandIndices_summary[1,1] = 1
    ari_RandIndices_summary[1,2] = ari_RandIndices_means
    ari_RandIndices_summary[1,3] = ari_RandIndices_std
    ari_RandIndices_summary[1,4] = ari_RandIndices_quantiles_mat[1]
    ari_RandIndices_summary[1,5] = ari_RandIndices_quantiles_mat[2]

    
    return ari_RandIndices_summary
end
function calc_time_variant_ARI_summarization(ari_vov,T;conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    ari_RandIndices_values = ari_vov
    
 
    ari_RandIndices_means = mean.(ari_RandIndices_values)
    ari_RandIndices_std = std.(ari_RandIndices_values)
    if length(ari_RandIndices_values[1]) > 1
        ari_RandIndices_conf_level = t_test.(ari_RandIndices_values; conf_level=conf_level)
        ari_RandIndices_quantiles_mat =reduce(vcat,reduce.(hcat,[ari_RandIndices_conf_level...]))
    else
        ari_RandIndices_quantiles_mat = ari_RandIndices_means .* ones(Float64,T,2)
    end
    
    ari_RandIndices_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
    ari_RandIndices_summary[:,1] = collect(1:T)
    ari_RandIndices_summary[:,2] = ari_RandIndices_means
    ari_RandIndices_summary[:,3] = ari_RandIndices_std
    ari_RandIndices_summary[:,4] = ari_RandIndices_quantiles_mat[:,1]
    ari_RandIndices_summary[:,5] = ari_RandIndices_quantiles_mat[:,2]

    
    return ari_RandIndices_summary
end
function calc_time_invariant_NMI_summarization(nmi_invar;conf_level=0.95)
    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    nmi_values = nmi_invar
    
 
    nmi_means = mean(nmi_values)
    nmi_std = std(nmi_values)
    if length(nmi_values) > 1
        nmi_conf_level = t_test(nmi_values; conf_level=conf_level)
        nmi_quantiles_mat =reduce(hcat,[nmi_conf_level...])
    else
        nmi_quantiles_mat = [nmi_means,nmi_means]
    end

    
    nmi_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
    nmi_summary[1,1] = 1
    nmi_summary[1,2] = nmi_means
    nmi_summary[1,3] = nmi_std
    nmi_summary[1,4] = nmi_quantiles_mat[1]
    nmi_summary[1,5] = nmi_quantiles_mat[2]

    
    return nmi_summary
end
function calc_time_variant_NMI_summarization(nmi_vov,T;conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    nmi_values = nmi_vov
    
 
    nmi_means = mean.(nmi_values)
    nmi_std = std.(nmi_values)
    if length(nmi_values[1]) > 1
        nmi_conf_level = t_test.(nmi_values; conf_level=conf_level)
        nmi_quantiles_mat = reduce(vcat,reduce.(hcat,[nmi_conf_level...])) 
    else
        nmi_quantiles_mat = nmi_means .* ones(Float64,T,2)
    end
    
    nmi_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
    nmi_summary[:,1] = collect(1:T)
    nmi_summary[:,2] = nmi_means
    nmi_summary[:,3] = nmi_std
    nmi_summary[:,4] = nmi_quantiles_mat[:,1]
    nmi_summary[:,5] = nmi_quantiles_mat[:,2]

    
    return nmi_summary
end
function calc_time_invariant_Vmeasure_summarization(vmeasure_invar;conf_level=0.95)



    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    vmeasure_values = vmeasure_invar
    
 
    vmeasure_means = mean(vmeasure_values)
    vmeasure_std = std(vmeasure_values)
    if length(vmeasure_values) > 1
        vmeasure_conf_level = t_test(vmeasure_values; conf_level=conf_level)
        vmeasure_quantiles_mat =reduce(hcat,[vmeasure_conf_level...])
    else
        vmeasure_quantiles_mat = [vmeasure_means,vmeasure_means]
    end
    
    vmeasure_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
    vmeasure_summary[1,1] = 1
    vmeasure_summary[1,2] = vmeasure_means
    vmeasure_summary[1,3] = vmeasure_std
    vmeasure_summary[1,4] = vmeasure_quantiles_mat[1]
    vmeasure_summary[1,5] = vmeasure_quantiles_mat[2]

    
    return vmeasure_summary
end
function calc_time_variant_Vmeasure_summarization(vmeasure_vov,T;conf_level=0.95)
    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    vmeasure_values = vmeasure_vov
    
 
    vmeasure_means = mean.(vmeasure_values)
    vmeasure_std = std.(vmeasure_values)
    if length(vmeasure_values[1]) > 1
        vmeasure_conf_level = t_test.(vmeasure_values; conf_level=conf_level)
        vmeasure_quantiles_mat = reduce(vcat,reduce.(hcat,[vmeasure_conf_level...]))
    else
        vmeasure_quantiles_mat = vmeasure_means .* ones(Float64,T,2)
    end 

    
    vmeasure_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
    vmeasure_summary[:,1] = collect(1:T)
    vmeasure_summary[:,2] = vmeasure_means
    vmeasure_summary[:,3] = vmeasure_std
    vmeasure_summary[:,4] = vmeasure_quantiles_mat[:,1]
    vmeasure_summary[:,5] = vmeasure_quantiles_mat[:,2]

    
    return vmeasure_summary
end
function calc_time_invariant_VarInfo_summarization(varinfo_invar;conf_level=0.95)
    varinfo_values = varinfo_invar
    
 
    varinfo_means = mean(varinfo_values)
    varinfo_std = std(varinfo_values)
    if length(varinfo_values) > 1
        varinfo_conf_level = t_test(varinfo_values; conf_level=conf_level)
        varinfo_quantiles_mat =reduce(hcat,[varinfo_conf_level...])
    else
        varinfo_quantiles_mat = [varinfo_means,varinfo_means]
    end
    
    varinfo_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
    varinfo_summary[1,1] = 1
    varinfo_summary[1,2] = varinfo_means
    varinfo_summary[1,3] = varinfo_std
    varinfo_summary[1,4] = varinfo_quantiles_mat[1]
    varinfo_summary[1,5] = varinfo_quantiles_mat[2]

    
    return varinfo_summary
end
function calc_time_variant_VarInfo_summarization(varinfo_vov,T;conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    varinfo_values = varinfo_vov
    
 
    varinfo_means = mean.(varinfo_values)
    varinfo_std = std.(varinfo_values)
    if length(varinfo_values[1]) > 1
        varinfo_conf_level = t_test.(varinfo_values; conf_level=conf_level)
        varinfo_quantiles_mat = reduce(vcat,reduce.(hcat,[varinfo_conf_level...])) 
    else
        varinfo_quantiles_mat = varinfo_means .* ones(Float64,T,2)
    end
    
    varinfo_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
    varinfo_summary[:,1] = collect(1:T)
    varinfo_summary[:,2] = varinfo_means
    varinfo_summary[:,3] = varinfo_std
    varinfo_summary[:,4] = varinfo_quantiles_mat[:,1]
    varinfo_summary[:,5] = varinfo_quantiles_mat[:,2]

    
    return varinfo_summary
end
function calc_time_invariant_Jaccard_summarization(jaccard_invar;conf_level=0.95)
    jaccard_values = jaccard_invar
    
 
    jaccard_means = mean(jaccard_values)
    jaccard_std = std(jaccard_values)
    if length(jaccard_values) > 1
        jaccard_conf_level = t_test(jaccard_values; conf_level=conf_level)
        jaccard_quantiles_mat =reduce(hcat,[jaccard_conf_level...])
    else
        jaccard_quantiles_mat = [jaccard_means,jaccard_means]
    end
    
    jaccard_summary = Matrix{Union{Float64,Int}}(undef,1,5) 
    jaccard_summary[1,1] = 1
    jaccard_summary[1,2] = jaccard_means
    jaccard_summary[1,3] = jaccard_std
    jaccard_summary[1,4] = jaccard_quantiles_mat[1]
    jaccard_summary[1,5] = jaccard_quantiles_mat[2]

    
    return jaccard_summary
end
function calc_time_variant_Jaccard_summarization(jaccard_vov,T;conf_level=0.95)


    #ari_vov,ri_vov,mirkinindx_vov,hubertindx_vov
    jaccard_values = jaccard_vov
    
 
    jaccard_means = mean.(jaccard_values)
    jaccard_std = std.(jaccard_values)
    if length(jaccard_values[1]) > 1
        jaccard_conf_level = t_test.(jaccard_values; conf_level=conf_level)
        jaccard_quantiles_mat = reduce(vcat,reduce.(hcat,[jaccard_conf_level...])) 
    else
        jaccard_quantiles_mat = jaccard_means .* ones(Float64,T,2)
    end
    
    jaccard_summary = Matrix{Union{Float64,Int}}(undef,T,5) 
    jaccard_summary[:,1] = collect(1:T)
    jaccard_summary[:,2] = jaccard_means
    jaccard_summary[:,3] = jaccard_std
    jaccard_summary[:,4] = jaccard_quantiles_mat[:,1]
    jaccard_summary[:,5] = jaccard_quantiles_mat[:,2]

    
    return jaccard_summary
end
function get_timeranges(N_t)
    T = length(N_t)
    starts = Vector{Int}(undef,T)
    ends = Vector{Int}(undef,T)
    
    for t in 1:T
        if t == 1
            starts[1] = 0 + 1
            ends[1] = 0 + N_t[1]
            continue
        end
        starts[t] = ends[t-1] + 1
        ends[t] = ends[t-1] + N_t[t]
    end
    return zip(starts,ends)
end




main(ARGS)