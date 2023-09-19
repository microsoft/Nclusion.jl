

function plotNumClusters(tchain, iterations)
    h = map(
        t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value))), 1:iterations
        );
    #Posterior Over Clusters
    clus_plot = histogram(h, xlabel = "Number of clusters", legend = false, bins=10);
    display(clus_plot)
    return clus_plot
end
function plotNumClusters(tchain, burn, iterations)
    h = map(
        t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value))), burn:iterations
        );
    #Posterior Over Clusters
    clus_plot = histogram(h, xlabel = "Number of clusters", legend = false, bins=10);
    display(clus_plot)
    return clus_plot
end

function plotNumClustersOverItr(tchain, burn, iterations)
    ht = map(t -> (t,length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value)))), burn:iterations);
    itr = [el[1] for el in ht]
    clus_num = [el[2] for el in ht]
    clus_plot = plot(itr,clus_num)
    display(clus_plot)
    return clus_plot
end
function plotNumClustersOverItr(tchain, iterations)
    ht = map(t -> (t,length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value)))), 1:iterations);
    itr = [el[1] for el in ht]
    clus_num = [el[2] for el in ht]
    clus_plot = plot(itr,clus_num)
    display(clus_plot)
    return clus_plot
end
function plotlambdaposteriorSamples(tchain,n_samples,K)
    λ_df = getlambdaposteriorDF(tchain,n_samples)
    pv = violin(reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K],legend=false,yticks = 0:maximum(Matrix(λ_df))+5)
    pv = boxplot!(pv, reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K], fillalpha=0.75, linewidth=2,legend=false)
    
    display(pv)
    return pv
end

function plotlambdaposteriorChainAndDensity(tchain,n_samples,K)
    λ_ids = getlambdaIDs(tchain)
    p_λ = plot(tchain[n_samples:end, λ_ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter);
    display(p_λ)
    return p_λ
end

function plotCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
    exclusive_rate, inclusive_rate = getEmpricalCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
    pCCR = plot(collect(1:T),inclusive_rate,labels=reshape(["Inclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:circle,legend = :outerbottom)
    pCCR = plot!(pCCR,collect(1:T),exclusive_rate,labels=reshape(["Exclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:diamond, linestyle=:dashdot, legend = :outerbottom)
    display(pCCR)
    return pCCR
end

function plotAvgPosteriorClusterMembershipAtTime(tchain,tp,T,true_z,KMax,KTrue,C_t; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster")
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    ticklabelStem = ticklabelStem
    ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
    priorlabelStem = priorlabelStem
    priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
    
    chosen_t = tp
    ylim_max = maximum(avg_counts_mat)
    p = groupedbar(avg_counts_mat[:,:,chosen_t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $chosen_t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
    return p
end
function plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,tp,KMax,KTrue; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster")
    ticklabelStem = ticklabelStem
    ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
    priorlabelStem = priorlabelStem
    priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
    
    chosen_t = tp
    ylim_max = maximum(avg_counts_mat)
    p = groupedbar(avg_counts_mat[:,:,chosen_t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $chosen_t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
    return p
end
function plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,tp,KMax,KTrue; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster")
    ticklabelStem = ticklabelStem
    ticklabel =  labelMaker("",KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
    priorlabelStem = priorlabelStem
    priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
    
    chosen_t = tp
    ylim_max = maximum(avg_counts_mat)
    p = groupedbar(avg_counts_mat[:,:,chosen_t], bar_position = :stack, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $chosen_t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency",  label=priorlabel, legend = :outerbottom) #ylims=(0,ylim_max),
    return p
end

function plot_MixtureWeights_post(param, param_name,t,KMax, param_full_name; truth=nothing)
    # plt.figure(figsize=figsize)
    ticklabel =  labelMaker("",KMax)
    p = boxplot(param, whis=[2.5, 97.5], showmeans=true, showfliers=false, xlabel = "Posterior mixture components", ylabel=param_full_name , title= "95% Credible Intervals for $(param_full_name)\n at time T = $(t)",color= reshape([:blue for i in 1:KMax],1,KMax) ,legend = false,xticks=(1:length(ticklabel), ticklabel))
    # plt.
    # plt.ylabel(param_full_name)
    # plt.title("95% Credible Intervals for $(param_full_name)")
    
    # if truth != nothing
    #     for line in truth
    #         plot!(p, line, ls=":")
    #     end
    # end
    return p
end

function makeAvgPosteriorClusterMembershipGIF(name,tchain,T,true_z,KMax,KTrue,C_t; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster", fps =3)
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    ticklabelStem = ticklabelStem
    ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
    priorlabelStem = priorlabelStem
    priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
    ylim_max = maximum(avg_counts_mat)
    anim = @animate for t in 1:T
        p = groupedbar(avg_counts_mat[:,:,t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
    
    end
    gif(anim, name, fps = fps)
end
function makeAvgPosteriorClusterMembershipGIF(name,avg_counts_mat,T,KMax,KTrue; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster", fps =3)
    ticklabelStem = ticklabelStem
    ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
    priorlabelStem = priorlabelStem
    priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
    ylim_max = maximum(avg_counts_mat)
    anim = @animate for t in 1:T
        p = groupedbar(avg_counts_mat[:,:,t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by True Membership at time T = $t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
    
    end
    gif(anim, name, fps = fps)
end


function makePosteriorClusterAssignmentTable_fromSubClusterMax(name,tchain,T,true_z,KMax,KTrue,C_t;latex_file = true)
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    clus_w_maxPop_vec = [get_clus_w_maxPopAtT(avg_counts_mat,t) for t in 1:T]
    clus_w_maxPop_mat = reduce(hcat, clus_w_maxPop_vec) #vs# hcat(clus_w_maxPop_vec...)
    row_lbl = "Posterior Cluster " .* string.(collect(1:KMax))
    c_data = hcat(row_lbl , clus_w_maxPop_mat)
    colHeadings = append!(["Posterior Cluster Index"], ["Cluster identity at Time, t = $t" for t in 1:T])
    colSubHeadings = append!(["[h ∈ {1,...,KMax}]"], ["[Based on max. pop. in cluster]" for i in 1:T])
    header = (colHeadings,colSubHeadings)
    
    if !latex_file
        hl_func = Highlighter((c_data,i,j)-> (c_data[i,j] in FindNonunique_func(c_data[:,j])) && (j>1) , crayon"red bold")
    else
        hl_func = LatexHighlighter((c_data,i,j)-> (c_data[i,j] in FindNonunique_func(c_data[:,j])) && (j>1) , ["color{red}", "textbf"])
    end
    
    
    io = IOBuffer()
    if !latex_file
        pretty_table(io, c_data;header = header, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green", highlighters = hl_func,alignment=:c, autowrap = true, linebreaks = true)
    else
        pretty_table(io, c_data;header = header, highlighters = hl_func,alignment=:c,backend = Val(:latex))
    end
    open(name, "w") do f
        write(f,String(take!(io)))
    end

end




function PosteriorPlotsAndSave!(tchain,thin_IterMCMC,KMax,KTrue,true_z_dict,true_z,T,C_t,relabel_tchain,thin_relabel_filepath,thin_original_filepath,thin_n_burin_samples,fps,plot_name_vec,job_name_stem,param_str,unique_time_id,plot_ext)
    original_prefix = "Orignial"
    PosteriorNumClustersPlotsAndSave!(tchain,thin_IterMCMC,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaSamplesPlotsAndSave!(tchain,thin_n_burin_samples,KMax,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorClusteringRatePlotsAndSave!(tchain,thin_n_burin_samples,T,KTrue,true_z_dict,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaChainPlotsAndSave!(tchain,thin_n_burin_samples,KMax,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaPosteriorValueHistogramPlotsAndSave!(tchain,thin_n_burin_samples,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaRunningMeanPlotsAndSave!(tchain,thin_n_burin_samples,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    # plotAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);

    PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!(tchain,T,true_z,KMax,KTrue,C_t,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id);
    

    relablled_prefix = "Relabelled"
    PosteriorNumClustersPlotsAndSave!(relabel_tchain, length(relabel_tchain) ,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaSamplesPlotsAndSave!(relabel_tchain,1,KMax,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorClusteringRatePlotsAndSave!(relabel_tchain,1,T,KTrue,true_z_dict,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaChainPlotsAndSave!(relabel_tchain,1,KMax,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaPosteriorValueHistogramPlotsAndSave!(relabel_tchain,1,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    PosteriorlambdaRunningMeanPlotsAndSave!(relabel_tchain,1,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
    # plotAvgPosteriorClusterMembershipAndSave!(relabel_tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext); 

    PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!(tchain,T,true_z,KMax,KTrue,C_t,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id);


end
function PosteriorNumClustersPlotsAndSave!(tchain,IterMCMC,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    numclusterplot = plotNumClusters(tchain,IterMCMC)#plotNumClusters(tchain, IterMCMC);
    plot_id = prefix *"NumberOfClusterPost"
    numclusterplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(numclusterplot, numclusterplot_name);
    push!(plot_name_vec, numclusterplot_name)
end
function PosteriorNumClustersPlotsAndSave!(tchain, num_burnin ,IterMCMC,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    numclusterplot = plotNumClusters(tchain, num_burnin, IterMCMC)#plotNumClusters(tchain, IterMCMC);
    plot_id = prefix *"NumberOfClusterPost"
    numclusterplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(numclusterplot, numclusterplot_name);
    push!(plot_name_vec, numclusterplot_name)
end
function PosteriorlambdaSamplesPlotsAndSave!(tchain,n_burin_samples,KMax,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    lambdaposteriorplot = plotlambdaposteriorSamples(tchain,n_burin_samples,KMax);
    plot_id = prefix *"ClusterValuesPost"
    lambdaposteriorplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(lambdaposteriorplot, lambdaposteriorplot_name);
    push!(plot_name_vec, lambdaposteriorplot_name)
end
function PosteriorClusteringRatePlotsAndSave!(tchain,n_burin_samples,T,KTrue,true_z_dict,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    clusteringRatePlot = plotCorrectClusteringRates(tchain,n_burin_samples,T,KTrue,true_z_dict);
    plot_id = prefix *"InclusionExclusionClusteringRates"
    clusteringRateplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(clusteringRatePlot, clusteringRateplot_name);
    push!(plot_name_vec, clusteringRateplot_name)
end
function PosteriorlambdaChainPlotsAndSave!(tchain,n_burin_samples,KMax,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    lambdaChainplot =plotlambdaposteriorChainAndDensity(tchain,n_burin_samples,KMax);
    plot_id = prefix *"LambdaPosteriorChain"
    lambdaChainplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(lambdaChainplot, lambdaChainplot_name);
    push!(plot_name_vec, lambdaChainplot_name)
end
function PosteriorlambdaPosteriorValueHistogramPlotsAndSave!(tchain,n_burin_samples,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    λ_ids = getlambdaIDs(tchain);
    lambdaHistplot = histogram(tchain[n_burin_samples:end, λ_ids, :])
    plot_id = prefix *"LambdaPosteriorValueHistogram"
    lambdaHistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(lambdaHistplot, lambdaHistplot_name);
    push!(plot_name_vec, lambdaHistplot_name)
end
function PosteriorlambdaRunningMeanPlotsAndSave!(tchain,n_burin_samples,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    λ_ids = getlambdaIDs(tchain);
    lambdaRunningMeanplot = meanplot(tchain[n_burin_samples:end, λ_ids, :])
    plot_id = prefix *"LambdaPosteriorRunningMean"
    lambdaRunningMeanplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(lambdaRunningMeanplot, lambdaRunningMeanplot_name);
    push!(plot_name_vec, lambdaRunningMeanplot_name)
end
function plotAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"AvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end
function plotAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"AvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end

function plotAvgPosteriorClusterMembershipAndSave2!(z_post_s,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency2(z_post_s,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"AvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end

function PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!(tchain,T,true_z,KMax,KTrue,C_t,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id)
    txt_table_ext = ".txt"
    latex_file_bool = false
    table_id = prefix *"PosteriorClusterAssignmentTable"
    txt_table_name = filepath * generate_filenameBase(table_id ,job_name_stem,param_str,unique_time_id) * txt_table_ext
    makePosteriorClusterAssignmentTable_fromSubClusterMax(txt_table_name,tchain,T,true_z,KMax,KTrue,C_t;latex_file = latex_file_bool)
    push!(plot_name_vec, txt_table_name)

    tex_table_ext = ".tex"
    latex_file_bool = true
    table_id = prefix *"PosteriorClusterAssignmentTable"
    tex_table_name = filepath * generate_filenameBase(table_id ,job_name_stem,param_str,unique_time_id) * tex_table_ext
    makePosteriorClusterAssignmentTable_fromSubClusterMax(tex_table_name,tchain,T,true_z,KMax,KTrue,C_t;latex_file = latex_file_bool)
    push!(plot_name_vec, tex_table_name)
end


function PosteriorMixtureWeightsAndSave!(tchain,π_,T,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    # datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
    πpost =  extract(tchain, "π_t", is_scalar=false);
    KMax = Int(size(πpost)[2]/T)
    sindx = 1
    eindx = KMax
    for t in 1:T
        param_t = πpost[:,sindx:eindx]
        dataHeatmapPlot = plot_MixtureWeights_post(param_t, :π_t,t,KMax, "mixture weights (π_t)", truth=π_);
        sindx = Int(eindx+1)
        eindx = Int(eindx+KMax)
        static_heat_id = prefix *"PosteriorMixtureWeightsAtTime$t"
        dataHeatmaprPlot_name = filepath * generate_filenameBase(static_heat_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataHeatmapPlot, dataHeatmaprPlot_name);
        push!(plot_name_vec, dataHeatmaprPlot_name)
    end

   
end

function plotStackedAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"StackedAvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end

function plotStackedAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"StackedAvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end

function plotStackedAvgPosteriorClusterMembershipAndSave2!(z_post_s,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency2(z_post_s,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"StackedAvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end
function plotStackedAvgPosteriorClusterMembershipAndSave2!(z_post_s,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,fileNameGeneratorFunc,prefix)
    ticklabelStem = "Posterior Cluster"
    priorlabelStem = "True Cluster"
    avg_counts_mat = get_average_posterior_cluster_frequency2(z_post_s,T,true_z,KMax,KTrue,C_t)
    for t in 1:T
        dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
        static_plot_id = prefix *"StackedAvgPosteriorClusterMembershipAtTime$t"
        dataAvgPosteriorClusterMembership_name = fileNameGeneratorFunc(static_plot_id)
        savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
        push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
    end

    gif_id = prefix *"PosteriorClusterMembership"
    dataAvgPosteriorClusterMembershipGIF_name = fileNameGeneratorFunc(gif_id)
    # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
    push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)

end


# module turingChainPlotting
#     # include("turingChainProcessing.jl")
#     # using .turingChainProcessing
    
#     using Random
#     using Distributions
#     using Turing
#     using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     using Test
#     import Debugger
#     using CSV,DataFrames
#     using LaTeXStrings, TypedTables, PrettyTables

#     curr_dir = ENV["PWD"]
#     src_dir = "/src/"
#     include(curr_dir*src_dir*"IOUtils.jl")
#     using .IOUtils
#     include(curr_dir*src_dir*"InferenceUtils.jl")
#     using .InferenceUtils.turingChainProcessing

#     export plotCorrectClusteringRates,plotlambdaposteriorSamples,plotNumClusters,plotlambdaposteriorChainAndDensity, plotNumClustersOverItr,makeAvgPosteriorClusterMembershipGIF,plotAvgPosteriorClusterMembershipAtTime,makePosteriorClusterAssignmentTable_fromSubClusterMax

#     function plotNumClusters(tchain, iterations)
#         h = map(
#             t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value))), 1:iterations
#             );
#         #Posterior Over Clusters
#         clus_plot = histogram(h, xlabel = "Number of clusters", legend = false, bins=10);
#         display(clus_plot)
#         return clus_plot
#     end
#     function plotNumClusters(tchain, burn, iterations)
#         h = map(
#             t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value))), burn:iterations
#             );
#         #Posterior Over Clusters
#         clus_plot = histogram(h, xlabel = "Number of clusters", legend = false, bins=10);
#         display(clus_plot)
#         return clus_plot
#     end

#     function plotNumClustersOverItr(tchain, burn, iterations)
#         ht = map(t -> (t,length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value)))), burn:iterations);
#         itr = [el[1] for el in ht]
#         clus_num = [el[2] for el in ht]
#         clus_plot = plot(itr,clus_num)
#         display(clus_plot)
#         return clus_plot
#     end
#     function plotNumClustersOverItr(tchain, iterations)
#         ht = map(t -> (t,length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value)))), 1:iterations);
#         itr = [el[1] for el in ht]
#         clus_num = [el[2] for el in ht]
#         clus_plot = plot(itr,clus_num)
#         display(clus_plot)
#         return clus_plot
#     end
#     function plotlambdaposteriorSamples(tchain,n_samples,K)
#         λ_df = getlambdaposteriorDF(tchain,n_samples)
#         pv = violin(reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K],legend=false,yticks = 0:maximum(Matrix(λ_df))+5)
#         pv = boxplot!(pv, reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K], fillalpha=0.75, linewidth=2,legend=false)
        
#         display(pv)
#         return pv
#     end

#     function plotlambdaposteriorChainAndDensity(tchain,n_samples,K)
#         λ_ids = getlambdaIDs(tchain)
#         p_λ = plot(tchain[n_samples:end, λ_ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter);
#         display(p_λ)
#         return p_λ
#     end
    
#     function plotCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
#         exclusive_rate, inclusive_rate = getEmpricalCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
#         pCCR = plot(collect(1:T),inclusive_rate,labels=reshape(["Inclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:circle,legend = :outerbottom)
#         pCCR = plot!(pCCR,collect(1:T),exclusive_rate,labels=reshape(["Exclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:diamond, linestyle=:dashdot, legend = :outerbottom)
#         display(pCCR)
#         return pCCR
#     end

#     function plotAvgPosteriorClusterMembershipAtTime(tchain,tp,T,true_z,KMax,KTrue,C_t; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster")
#         avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         ticklabelStem = ticklabelStem
#         ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
#         priorlabelStem = priorlabelStem
#         priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
        
#         chosen_t = tp
#         ylim_max = maximum(avg_counts_mat)
#         p = groupedbar(avg_counts_mat[:,:,chosen_t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $chosen_t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
#         return p
#     end
#     function plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,tp,KMax,KTrue; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster")
#         ticklabelStem = ticklabelStem
#         ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
#         priorlabelStem = priorlabelStem
#         priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
        
#         chosen_t = tp
#         ylim_max = maximum(avg_counts_mat)
#         p = groupedbar(avg_counts_mat[:,:,chosen_t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $chosen_t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
#         return p
#     end
#     function plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,tp,KMax,KTrue; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster")
#         ticklabelStem = ticklabelStem
#         ticklabel =  labelMaker("",KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
#         priorlabelStem = priorlabelStem
#         priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
        
#         chosen_t = tp
#         ylim_max = maximum(avg_counts_mat)
#         p = groupedbar(avg_counts_mat[:,:,chosen_t], bar_position = :stack, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $chosen_t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency",  label=priorlabel, legend = :outerbottom) #ylims=(0,ylim_max),
#         return p
#     end

#     function makeAvgPosteriorClusterMembershipGIF(name,tchain,T,true_z,KMax,KTrue,C_t; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster", fps =3)
#         avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         ticklabelStem = ticklabelStem
#         ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
#         priorlabelStem = priorlabelStem
#         priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
#         ylim_max = maximum(avg_counts_mat)
#         anim = @animate for t in 1:T
#             p = groupedbar(avg_counts_mat[:,:,t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by  True Membership at time T = $t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
        
#         end
#         gif(anim, name, fps = fps)
#     end
#     function makeAvgPosteriorClusterMembershipGIF(name,avg_counts_mat,T,KMax,KTrue; ticklabelStem = "Posterior Cluster",priorlabelStem = "True Cluster", fps =3)
#         ticklabelStem = ticklabelStem
#         ticklabel =  labelMaker(ticklabelStem,KMax) #["Posterior Cluster 1" "Posterior Cluster 2" "Posterior Cluster 3"] 
#         priorlabelStem = priorlabelStem
#         priorlabel = labelMaker(priorlabelStem,KTrue)#["True Cluster 1" "True Cluster 2" "True Cluster 3"]
#         ylim_max = maximum(avg_counts_mat)
#         anim = @animate for t in 1:T
#             p = groupedbar(avg_counts_mat[:,:,t], bar_position = :dodge, bar_width=0.7, xticks=(1:length(ticklabel), ticklabel), title = "Posterior Cluster Membership \n Stratified by True Membership at time T = $t", xlabel = "Posterior Cluster Index", ylabel = "Average Cluster Frequency", ylims=(0,ylim_max), label=priorlabel, legend = :outerbottom)
        
#         end
#         gif(anim, name, fps = fps)
#     end


#     function makePosteriorClusterAssignmentTable_fromSubClusterMax(name,tchain,T,true_z,KMax,KTrue,C_t;latex_file = true)
#         avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         clus_w_maxPop_vec = [get_clus_w_maxPopAtT(avg_counts_mat,t) for t in 1:T]
#         clus_w_maxPop_mat = reduce(hcat, clus_w_maxPop_vec) #vs# hcat(clus_w_maxPop_vec...)
#         row_lbl = "Posterior Cluster " .* string.(collect(1:KMax))
#         c_data = hcat(row_lbl , clus_w_maxPop_mat)
#         colHeadings = append!(["Posterior Cluster Index"], ["Cluster identity at Time, t = $t" for t in 1:T])
#         colSubHeadings = append!(["[h ∈ {1,...,KMax}]"], ["[Based on max. pop. in cluster]" for i in 1:T])
#         header = (colHeadings,colSubHeadings)
        
#         if !latex_file
#             hl_func = Highlighter((c_data,i,j)-> (c_data[i,j] in FindNonunique_func(c_data[:,j])) && (j>1) , crayon"red bold")
#         else
#             hl_func = LatexHighlighter((c_data,i,j)-> (c_data[i,j] in FindNonunique_func(c_data[:,j])) && (j>1) , ["color{red}", "textbf"])
#         end
        
        
#         io = IOBuffer()
#         if !latex_file
#             pretty_table(io, c_data;header = header, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green", highlighters = hl_func,alignment=:c, autowrap = true, linebreaks = true)
#         else
#             pretty_table(io, c_data;header = header, highlighters = hl_func,alignment=:c,backend = Val(:latex))
#         end
#         open(name, "w") do f
#             write(f,String(take!(io)))
#         end
    
#     end
    

#     export PosteriorPlotsAndSave!,PosteriorNumClustersPlotsAndSave!,PosteriorlambdaSamplesPlotsAndSave!,PosteriorClusteringRatePlotsAndSave!,PosteriorlambdaChainPlotsAndSave!,PosteriorlambdaPosteriorValueHistogramPlotsAndSave!,PosteriorlambdaRunningMeanPlotsAndSave!,plotAvgPosteriorClusterMembershipAndSave!,PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!

#     function PosteriorPlotsAndSave!(tchain,thin_IterMCMC,KMax,KTrue,true_z_dict,true_z,T,C_t,relabel_tchain,thin_relabel_filepath,thin_original_filepath,thin_n_burin_samples,fps,plot_name_vec,job_name_stem,param_str,unique_time_id,plot_ext)
#         original_prefix = "Orignial"
#         PosteriorNumClustersPlotsAndSave!(tchain,thin_IterMCMC,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaSamplesPlotsAndSave!(tchain,thin_n_burin_samples,KMax,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorClusteringRatePlotsAndSave!(tchain,thin_n_burin_samples,T,KTrue,true_z_dict,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaChainPlotsAndSave!(tchain,thin_n_burin_samples,KMax,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaPosteriorValueHistogramPlotsAndSave!(tchain,thin_n_burin_samples,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaRunningMeanPlotsAndSave!(tchain,thin_n_burin_samples,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         # plotAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id,plot_ext);

#         PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!(tchain,T,true_z,KMax,KTrue,C_t,plot_name_vec,thin_original_filepath,original_prefix,job_name_stem,param_str,unique_time_id);
        
    
#         relablled_prefix = "Relabelled"
#         PosteriorNumClustersPlotsAndSave!(relabel_tchain, length(relabel_tchain) ,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaSamplesPlotsAndSave!(relabel_tchain,1,KMax,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorClusteringRatePlotsAndSave!(relabel_tchain,1,T,KTrue,true_z_dict,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaChainPlotsAndSave!(relabel_tchain,1,KMax,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaPosteriorValueHistogramPlotsAndSave!(relabel_tchain,1,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         PosteriorlambdaRunningMeanPlotsAndSave!(relabel_tchain,1,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext);
#         # plotAvgPosteriorClusterMembershipAndSave!(relabel_tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id,plot_ext); 

#         PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!(tchain,T,true_z,KMax,KTrue,C_t,plot_name_vec,thin_relabel_filepath,relablled_prefix,job_name_stem,param_str,unique_time_id);


#     end
#     function PosteriorNumClustersPlotsAndSave!(tchain,IterMCMC,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         numclusterplot = plotNumClusters(tchain,IterMCMC)#plotNumClusters(tchain, IterMCMC);
#         plot_id = prefix *"NumberOfClusterPost"
#         numclusterplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(numclusterplot, numclusterplot_name);
#         push!(plot_name_vec, numclusterplot_name)
#     end
#     function PosteriorNumClustersPlotsAndSave!(tchain, num_burnin ,IterMCMC,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         numclusterplot = plotNumClusters(tchain, num_burnin, IterMCMC)#plotNumClusters(tchain, IterMCMC);
#         plot_id = prefix *"NumberOfClusterPost"
#         numclusterplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(numclusterplot, numclusterplot_name);
#         push!(plot_name_vec, numclusterplot_name)
#     end
#     function PosteriorlambdaSamplesPlotsAndSave!(tchain,n_burin_samples,KMax,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         lambdaposteriorplot = plotlambdaposteriorSamples(tchain,n_burin_samples,KMax);
#         plot_id = prefix *"ClusterValuesPost"
#         lambdaposteriorplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(lambdaposteriorplot, lambdaposteriorplot_name);
#         push!(plot_name_vec, lambdaposteriorplot_name)
#     end
#     function PosteriorClusteringRatePlotsAndSave!(tchain,n_burin_samples,T,KTrue,true_z_dict,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         clusteringRatePlot = plotCorrectClusteringRates(tchain,n_burin_samples,T,KTrue,true_z_dict);
#         plot_id = prefix *"InclusionExclusionClusteringRates"
#         clusteringRateplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(clusteringRatePlot, clusteringRateplot_name);
#         push!(plot_name_vec, clusteringRateplot_name)
#     end
#     function PosteriorlambdaChainPlotsAndSave!(tchain,n_burin_samples,KMax,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         lambdaChainplot =plotlambdaposteriorChainAndDensity(tchain,n_burin_samples,KMax);
#         plot_id = prefix *"LambdaPosteriorChain"
#         lambdaChainplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(lambdaChainplot, lambdaChainplot_name);
#         push!(plot_name_vec, lambdaChainplot_name)
#     end
#     function PosteriorlambdaPosteriorValueHistogramPlotsAndSave!(tchain,n_burin_samples,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         λ_ids = getlambdaIDs(tchain);
#         lambdaHistplot = histogram(tchain[n_burin_samples:end, λ_ids, :])
#         plot_id = prefix *"LambdaPosteriorValueHistogram"
#         lambdaHistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(lambdaHistplot, lambdaHistplot_name);
#         push!(plot_name_vec, lambdaHistplot_name)
#     end
#     function PosteriorlambdaRunningMeanPlotsAndSave!(tchain,n_burin_samples,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         λ_ids = getlambdaIDs(tchain);
#         lambdaRunningMeanplot = meanplot(tchain[n_burin_samples:end, λ_ids, :])
#         plot_id = prefix *"LambdaPosteriorRunningMean"
#         lambdaRunningMeanplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(lambdaRunningMeanplot, lambdaRunningMeanplot_name);
#         push!(plot_name_vec, lambdaRunningMeanplot_name)
#     end
#     function plotAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         ticklabelStem = "Posterior Cluster"
#         priorlabelStem = "True Cluster"
#         avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         for t in 1:T
#             dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
#             static_plot_id = prefix *"AvgPosteriorClusterMembershipAtTime$t"
#             dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#             savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
#             push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
#         end
    
#         gif_id = prefix *"PosteriorClusterMembership"
#         dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
#         push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)
    
#     end
#     function plotAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         ticklabelStem = "Posterior Cluster"
#         priorlabelStem = "True Cluster"
#         avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         for t in 1:T
#             dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTime(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
#             static_plot_id = prefix *"AvgPosteriorClusterMembershipAtTime$t"
#             dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#             savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
#             push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
#         end
    
#         gif_id = prefix *"PosteriorClusterMembership"
#         dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
#         push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)
    
#     end

#     function plotStackedAvgPosteriorClusterMembershipAndSave!(tchain,T,true_z,KMax,KTrue,C_t,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         ticklabelStem = "Posterior Cluster"
#         priorlabelStem = "True Cluster"
#         avg_counts_mat = get_average_posterior_cluster_frequency(tchain,T,true_z,KMax,KTrue,C_t)
#         for t in 1:T
#             dataAvgPosteriorClusterMembership = plotAvgPosteriorClusterMembershipAtTimeStacked(avg_counts_mat,t,KMax,KTrue; ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
#             static_plot_id = prefix *"StackedAvgPosteriorClusterMembershipAtTime$t"
#             dataAvgPosteriorClusterMembership_name = filepath * generate_filenameBase(static_plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#             savefig(dataAvgPosteriorClusterMembership, dataAvgPosteriorClusterMembership_name);
#             push!(plot_name_vec, dataAvgPosteriorClusterMembership_name)
#         end
    
#         gif_id = prefix *"PosteriorClusterMembership"
#         dataAvgPosteriorClusterMembershipGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         # makeAvgPosteriorClusterMembershipGIF(dataAvgPosteriorClusterMembershipGIF_name,avg_counts_mat,T,KMax,KTrue;fps = fps,ticklabelStem = ticklabelStem ,priorlabelStem = priorlabelStem )
#         push!(plot_name_vec, dataAvgPosteriorClusterMembershipGIF_name)
    
#     end
    
#     function PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!(tchain,T,true_z,KMax,KTrue,C_t,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id)
#         txt_table_ext = ".txt"
#         latex_file_bool = false
#         table_id = prefix *"PosteriorClusterAssignmentTable"
#         txt_table_name = filepath * generate_filenameBase(table_id ,job_name_stem,param_str,unique_time_id) * txt_table_ext
#         makePosteriorClusterAssignmentTable_fromSubClusterMax(txt_table_name,tchain,T,true_z,KMax,KTrue,C_t;latex_file = latex_file_bool)
#         push!(plot_name_vec, txt_table_name)
    
#         tex_table_ext = ".tex"
#         latex_file_bool = true
#         table_id = prefix *"PosteriorClusterAssignmentTable"
#         tex_table_name = filepath * generate_filenameBase(table_id ,job_name_stem,param_str,unique_time_id) * tex_table_ext
#         makePosteriorClusterAssignmentTable_fromSubClusterMax(tex_table_name,tchain,T,true_z,KMax,KTrue,C_t;latex_file = latex_file_bool)
#         push!(plot_name_vec, tex_table_name)
#     end

#     # #DO NOT EXPORT HERE!!!! EXPORT IN PlottingUtils!
#     # function labelMaker(item,K)
#     #     return reshape([ item*" "*string(i) for i in 1:K],1,K)
#     # end




# end
