function plotDataGeneratingDist(μ,σ²,mixing_prob)
    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²) 
    data_generating_distribution  = MixtureModel(NegativeBinomial.(r,p),mixing_prob)
    dgd = plot(data_generating_distribution)
    display(dgd)
    return dgd
end
function plotNBMMDataGeneratingDist(μ,σ²,mixing_prob)
    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²) 
    data_generating_distribution  = MixtureModel(NegativeBinomial.(r,p),mixing_prob)
    dgd = plot(data_generating_distribution)
    display(dgd)
    return dgd
end
function plotPMMDataGeneratingDist(μ,mixing_prob)
    data_generating_distribution  = MixtureModel(NegativeBinomial.(μ),mixing_prob)
    dgd = plot(data_generating_distribution)
    display(dgd)
    return dgd
end

function plotDataGeneratingDistAndSave!(μ,σ²,mixing_prob,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
    plot_id = prefix *"DataGeneratingDistribution"
    datagenDistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(datagenDistPlot, datagenDistplot_name);
    push!(plot_name_vec, datagenDistplot_name)
end

function plotNBMMDataGeneratingDistAndSave!(μ,σ²,mixing_prob,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    datagenDistPlot = plotNBMMDataGeneratingDist(μ,σ²,mixing_prob)
    plot_id = prefix *"NBMMDataGeneratingDistribution"
    datagenDistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(datagenDistPlot, datagenDistplot_name);
    push!(plot_name_vec, datagenDistplot_name)
end
function plotPMMDataGeneratingDistAndSave!(μ,mixing_prob,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    datagenDistPlot = plotPMMDataGeneratingDist(μ,mixing_prob)
    plot_id = prefix *"PMMDataGeneratingDistribution"
    datagenDistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
    savefig(datagenDistPlot, datagenDistplot_name);
    push!(plot_name_vec, datagenDistplot_name)
end


function makeDataScatterGIF(name,x,true_z,Ktrue,  C_t,T;fps = 3)
    
    x_array =  reduce(vcat,map(el -> reduce(vcat,el),x))
    max_data_val = maximum(x_array)
    max_indx_val = maximum(vcat(C_t...))

    anim = @animate for t in 1:T
        x_tp = map(el -> reduce(vcat,el),x[t])
        p = scatter(collect(1:C_t[t]),x_tp , group=true_z[t], ylims=(0.0,max_data_val), xlims=(0.0,max_indx_val),title = "Gene Expression of Cells at time T = $t", xlabel = "Cell index", ylabel = "Gene Expression Count",labels=reshape(["Cluster $i" for i in 1:Ktrue],1,Ktrue))
    
    end
    gif(anim, name, fps = fps)
end
function plotDataScatterAtTime(x,true_z,Ktrue,C_t,t)
    x_array =  reduce(vcat,map(el -> reduce(vcat,el),x))
    max_data_val = maximum(x_array)
    max_indx_val = maximum(vcat(C_t...))
    x_tp = map(el -> reduce(vcat,el),x[t])
    p = scatter(collect(1:C_t[t]),x_tp , group=true_z[t], ylims=(0.0,max_data_val), xlims=(0.0,max_indx_val),title = "Gene Expression of Cells at time T = $t", xlabel = "Cell index", ylabel = "Gene Expression Count",labels=reshape(["Cluster $i" for i in 1:Ktrue],1,Ktrue))
    return p
end

function plotDataScatterAndSave!(x,true_z,C_t,T,Ktrue,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    # datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
    for t in 1:T
        dataScatterPlot = plotDataScatterAtTime(x,true_z,Ktrue,C_t,t)
        static_scatter_id = prefix *"DataScatterScatterAtTime$t"
        dataScatterPlot_name = filepath * generate_filenameBase(static_scatter_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataScatterPlot, dataScatterPlot_name);
        push!(plot_name_vec, dataScatterPlot_name)
    end
    gif_id = prefix *"DataScatterGIF"
    datascatterGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # makeDataScatterGIF(datascatterGIF_name,x,true_z,Ktrue,C_t,T;fps = fps)
    push!(plot_name_vec, datascatterGIF_name)
   
end

#DO NOT EXPORT HERE!!!! EXPORT IN PlottingUtils!
# function labelMaker(item,K)
#     return reshape([ item*" "*string(i) for i in 1:K],1,K)
# end

function plotGGenesDataScatterAtTime(x,t,true_z,KTrue;show_plt=false)
    X = permutedims(hcat(x[t]...))
    G = size(X)[2]
    Z = true_z[t]
    data = hcat(Z,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')
    # c_vec = []
    # for c in 1:3
    #     c_val = c * one(eltype(data[:,1]))
    #     push!(c_vec, data_transformed[])
    # end
    c_val = 1 * one(eltype(data[:,1]))
    c_clus = X_transformed[:, Z.==c_val]
    p_data  = scatter(c_clus[1,:], c_clus[2,:],labels="Cluster 1 ")
    plot!(p_data,xlabel="PC1",ylabel="PC2",title = "PCA Projection (Components 1 and 2) of \n $G Genes Expressed in Cells at time T = $t")
    if !isone(length(unique(data[:,1])))
        for c in 2:KTrue
            c_val = c * one(eltype(data[:,1]))
            c_clus = X_transformed[:,Z.==c_val]
            scatter!(p_data,c_clus[1,:], c_clus[2,:],labels="Cluster $c ")
            plot!(p_data,xlabel="PC1",ylabel="PC2",title = "PCA Projection (Components 1 and 2) of \n $G Genes Expressed in Cells at time T = $t")
        end      
    end
    if show_plt
        display(p_data)
    end
    return p_data
end

function plotGGenesDataScatterAndSave!(x,true_z,T,G,KTrue,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext;gene_thres=200)
    # datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
    for t in 1:T
        if G == 1
            dataScatterPlot = plotDataScatterAtTime(x,true_z,KTrue,C_t,t) 
        else
            dataScatterPlot = plotGGenesDataScatterAtTime(x,t,true_z,KTrue)
        end
        static_scatter_id = prefix *"GGenesDataScatterScatterAtTime$t"
        dataScatterPlot_name = filepath * generate_filenameBase(static_scatter_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataScatterPlot, dataScatterPlot_name);
        push!(plot_name_vec, dataScatterPlot_name)
    end
    if G >2
        plotTsne = true
        if G > gene_thres 
            plotTsne = false
        end
        all_time_scatter_PCA = plot_All_Time_Points_PCA(x,true_z,G,KTrue)
        static_scatter_id2_PCA = prefix *"GGenesDataScatterScatterAllTimepointsPCA"
        all_time_scatter_PCA_name = filepath * generate_filenameBase(static_scatter_id2_PCA,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(all_time_scatter_PCA, all_time_scatter_PCA_name);
        push!(plot_name_vec, all_time_scatter_PCA_name)
        if plotTsne
            all_time_scatter_TSNE = plot_All_Time_Points_TSNE(x,true_z,G;tsne_transform = nothing, todisplay=false)
            static_scatter_id2_TSNE = prefix *"GGenesDataScatterScatterAllTimepointsTSNE"
            all_time_scatter_TSNE_name = filepath * generate_filenameBase(static_scatter_id2_TSNE,job_name_stem,param_str,unique_time_id) * plot_ext
            savefig(all_time_scatter_TSNE, all_time_scatter_TSNE_name);
            push!(plot_name_vec, all_time_scatter_TSNE_name)
        end
    end
    # gif_id = prefix *"DataScatterGIF"
    # datascatterGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # #NEED TO MAKE FOR G GENE CASE!! makeDataScatterGIF(datascatterGIF_name,x,true_z,Ktrue,C_t,T;fps = fps)
    # push!(plot_name_vec, datascatterGIF_name)
   
end
function plotGGenesDataScatterAndSave!(x,true_z,T,G,KTrue,fps,plot_name_vec,fileNameGeneratorFunc,prefix;gene_thres=200)
    # datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
    for t in 1:T
        if G == 1
            dataScatterPlot = plotDataScatterAtTime(x,true_z,KTrue,C_t,t) 
        else
            dataScatterPlot = plotGGenesDataScatterAtTime(x,t,true_z,KTrue)
        end
        static_scatter_id = prefix *"GGenesDataScatterScatterAtTime$t"
        dataScatterPlot_name = fileNameGeneratorFunc(static_scatter_id)
        savefig(dataScatterPlot, dataScatterPlot_name);
        push!(plot_name_vec, dataScatterPlot_name)
    end
    if G >2
        plotTsne = true
        if G > gene_thres 
            plotTsne = false
        end
        all_time_scatter_PCA = plot_All_Time_Points_PCA(x,true_z,G,KTrue)
        static_scatter_id2_PCA = prefix *"GGenesDataScatterScatterAllTimepointsPCA"
        all_time_scatter_PCA_name = fileNameGeneratorFunc(static_scatter_id2_PCA)
        savefig(all_time_scatter_PCA, all_time_scatter_PCA_name);
        push!(plot_name_vec, all_time_scatter_PCA_name)
        if plotTsne
            all_time_scatter_TSNE = plot_All_Time_Points_TSNE(x,true_z,G;tsne_transform = nothing, todisplay=false)
            static_scatter_id2_TSNE = prefix *"GGenesDataScatterScatterAllTimepointsTSNE"
            all_time_scatter_TSNE_name= fileNameGeneratorFunc(static_scatter_id2_TSNE)
            savefig(all_time_scatter_TSNE, all_time_scatter_TSNE_name);
            push!(plot_name_vec, all_time_scatter_TSNE_name)
        end
    end
    # gif_id = prefix *"DataScatterGIF"
    # datascatterGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
    # #NEED TO MAKE FOR G GENE CASE!! makeDataScatterGIF(datascatterGIF_name,x,true_z,Ktrue,C_t,T;fps = fps)
    # push!(plot_name_vec, datascatterGIF_name)
   
end

function make_heatmapPlotatT(x,true_z,t,C_t)
    X = x[t]
    Z = true_z[t]
    C = C_t[t]
    cell_indx = collect(1:C)
    tt = [(z,id) for (z , id) in zip(Z,cell_indx)]
    new_order = [el[2] for el in  sort(tt, by= gg -> gg[1])]
    new_X = X[new_order]
    new_Z = Z[new_order]
    h = permutedims(reduce(hcat,new_X))
    p = heatmap(cor(h,dims=2), xlabel="Cell Index",ylabel="Cell Index",title = "Correlation of Cell Expression \n at time $(t)")
    return p
end

function plotGGenesHeatMapAndSave!(x,true_z,T,C_t,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
    # datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
    for t in 1:T
        dataHeatmapPlot = make_heatmapPlotatT(x,true_z,t,C_t)
        static_heat_id = prefix *"GGenesDataHeatmapAtTime$t"
        dataHeatmaprPlot_name = filepath * generate_filenameBase(static_heat_id,job_name_stem,param_str,unique_time_id) * plot_ext
        savefig(dataHeatmapPlot, dataHeatmaprPlot_name);
        push!(plot_name_vec, dataHeatmaprPlot_name)
    end

   
end



function plot_All_Time_Points_PCA(x,z,G,Ktrue)
    T = length(x)

    X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    data = hcat(Z,X)
    M_ = fit(PCA, data[:,2:end]'; maxoutdim=2, pratio=1.0)
    X_transformed =  MultivariateStats.transform(M_, data[:,2:end]')

    c_val = 1 * one(eltype(data[:,1]))
    c_clus = X_transformed[:, Z .==c_val]
    p_data  = scatter(c_clus[1,:], c_clus[2,:],labels="Cluster 1 ")
    for c in 2:Ktrue
        c_val = c * one(eltype(data[:,1]))
        c_clus = X_transformed[:,Z.==c_val]
        scatter!(p_data,c_clus[1,:], c_clus[2,:],labels="Cluster $c ")
    end
    plot!(p_data,xlabel="PC1",ylabel="PC2",title = "PCA Projection (Components 1 and 2) of \n $G Genes Expressed in Cells across all time points")
    return p_data
end
function plot_All_Time_Points_TSNE(x,z,G;tsne_transform = nothing, todisplay=true)


    rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())

    T = length(x)
    
    vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T]))
    X = reduce(vcat,[permutedims(reduce(hcat,x[t])) for t in 1:T])#vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
    Z = vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T])) #vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
    K = Int(maximum(unique(Z)))
    clusterIDs =  sort(unique(Z))
    data = hcat(Z,X)
    
    X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
    if isnothing(tsne_transform)
        Y =  tsne(X_transformed, 2, 50, 1000, 20.0);
    else
        Y = tsne_transform
    end
    plotTitle = "TSne Projection of $G Genes \n Expressed in Cells across $T time points"
    theplot = scatter(Y[:,1], Y[:,2], marker=(4,4,:circle,stroke(0)),group = Z,  colors=Int.(Z), labels=labelMaker("Cluster ", K),legend=:outertopright)
    plot!(theplot,title = plotTitle,xlabel="TSne 1",ylabel="TSne 2")
    if todisplay
        display(theplot)
    end
    return theplot
end


# module syntheticDataPlotting
#     include("MathUtils.jl")
#     using .MathUtils

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

#     # Want to plot the distribution over cells being clustered correctly together
#     export plotDataGeneratingDist,plotNBMMDataGeneratingDist,plotPMMDataGeneratingDist, plotDataGeneratingDistAndSave!,plotNBMMDataGeneratingDistAndSave!,plotPMMDataGeneratingDistAndSave!,plotDataScatterAtTime,makeDataScatterGIF,plotDataScatterAndSave!
#     function plotDataGeneratingDist(μ,σ²,mixing_prob)
#         r = calc_r(μ,σ²)
#         p = calc_p(μ,σ²) 
#         data_generating_distribution  = MixtureModel(NegativeBinomial.(r,p),mixing_prob)
#         dgd = plot(data_generating_distribution)
#         display(dgd)
#         return dgd
#     end
#     function plotNBMMDataGeneratingDist(μ,σ²,mixing_prob)
#         r = calc_r(μ,σ²)
#         p = calc_p(μ,σ²) 
#         data_generating_distribution  = MixtureModel(NegativeBinomial.(r,p),mixing_prob)
#         dgd = plot(data_generating_distribution)
#         display(dgd)
#         return dgd
#     end
#     function plotPMMDataGeneratingDist(μ,mixing_prob)
#         data_generating_distribution  = MixtureModel(NegativeBinomial.(μ),mixing_prob)
#         dgd = plot(data_generating_distribution)
#         display(dgd)
#         return dgd
#     end

#     function plotDataGeneratingDistAndSave!(μ,σ²,mixing_prob,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
#         plot_id = prefix *"DataGeneratingDistribution"
#         datagenDistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(datagenDistPlot, datagenDistplot_name);
#         push!(plot_name_vec, datagenDistplot_name)
#     end

#     function plotNBMMDataGeneratingDistAndSave!(μ,σ²,mixing_prob,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         datagenDistPlot = plotNBMMDataGeneratingDist(μ,σ²,mixing_prob)
#         plot_id = prefix *"NBMMDataGeneratingDistribution"
#         datagenDistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(datagenDistPlot, datagenDistplot_name);
#         push!(plot_name_vec, datagenDistplot_name)
#     end
#     function plotPMMDataGeneratingDistAndSave!(μ,mixing_prob,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         datagenDistPlot = plotPMMDataGeneratingDist(μ,mixing_prob)
#         plot_id = prefix *"PMMDataGeneratingDistribution"
#         datagenDistplot_name = filepath * generate_filenameBase(plot_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         savefig(datagenDistPlot, datagenDistplot_name);
#         push!(plot_name_vec, datagenDistplot_name)
#     end

    
#     function makeDataScatterGIF(name,x,true_z,Ktrue,  C_t,T;fps = 3)
        
#         x_array =  reduce(vcat,map(el -> reduce(vcat,el),x))
#         max_data_val = maximum(x_array)
#         max_indx_val = maximum(vcat(C_t...))

#         anim = @animate for t in 1:T
#             x_tp = map(el -> reduce(vcat,el),x[t])
#             p = scatter(collect(1:C_t[t]),x_tp , group=true_z[t], ylims=(0.0,max_data_val), xlims=(0.0,max_indx_val),title = "Gene Expression of Cells at time T = $t", xlabel = "Cell index", ylabel = "Gene Expression Count",labels=reshape(["Cluster $i" for i in 1:Ktrue],1,Ktrue))
        
#         end
#         gif(anim, name, fps = fps)
#     end
#     function plotDataScatterAtTime(x,true_z,Ktrue,C_t,t)
#         x_array =  reduce(vcat,map(el -> reduce(vcat,el),x))
#         max_data_val = maximum(x_array)
#         max_indx_val = maximum(vcat(C_t...))
#         x_tp = map(el -> reduce(vcat,el),x[t])
#         p = scatter(collect(1:C_t[t]),x_tp , group=true_z[t], ylims=(0.0,max_data_val), xlims=(0.0,max_indx_val),title = "Gene Expression of Cells at time T = $t", xlabel = "Cell index", ylabel = "Gene Expression Count",labels=reshape(["Cluster $i" for i in 1:Ktrue],1,Ktrue))
#         return p
#     end

#     function plotDataScatterAndSave!(x,true_z,C_t,T,Ktrue,fps,plot_name_vec,filepath,prefix,job_name_stem,param_str,unique_time_id,plot_ext)
#         # datagenDistPlot = plotDataGeneratingDist(μ,σ²,mixing_prob)
#         for t in 1:T
#             dataScatterPlot = plotDataScatterAtTime(x,true_z,Ktrue,C_t,t)
#             static_scatter_id = prefix *"DataScatterScatterAtTime$t"
#             dataScatterPlot_name = filepath * generate_filenameBase(static_scatter_id,job_name_stem,param_str,unique_time_id) * plot_ext
#             savefig(dataScatterPlot, dataScatterPlot_name);
#             push!(plot_name_vec, dataScatterPlot_name)
#         end
#         gif_id = prefix *"DataScatterGIF"
#         datascatterGIF_name = filepath * generate_filenameBase(gif_id,job_name_stem,param_str,unique_time_id) * plot_ext
#         # makeDataScatterGIF(datascatterGIF_name,x,true_z,Ktrue,C_t,T;fps = fps)
#         push!(plot_name_vec, datascatterGIF_name)
       
#     end
    
#     #DO NOT EXPORT HERE!!!! EXPORT IN PlottingUtils!
#     function labelMaker(item,K)
#         return reshape([ item*" "*string(i) for i in 1:K],1,K)
#     end

# end

