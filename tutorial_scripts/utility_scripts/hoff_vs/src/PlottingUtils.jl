module PlottingUtils
    

    using Random
    using Distributions
    using Turing
    using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
    using TSne, MultivariateStats, Clustering
    using Test
    import Debugger
    using CSV,DataFrames
    using LaTeXStrings, TypedTables, PrettyTables,ColorSchemes,Gnuplot, Colors, ColorSchemes
    using Gnuplot

    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"
    include(curr_dir*src_dir*"IOUtils.jl")
    using .IOUtils
    include(curr_dir*src_dir*"MathUtils.jl")
    using .MathUtils
    include(curr_dir*src_dir*"MCMCInferenceUtils.jl")
    using .MCMCInferenceUtils
    include(curr_dir*src_dir*"MetricsUtils.jl")
    using .MetricsUtils
    include(curr_dir*src_dir*"DataPreprocessingUtils.jl")
    using .DataPreprocessingUtils
    include(curr_dir*src_dir*"TidyUtils.jl")
    using .TidyUtils
    # using .InferenceUtils.turingChainProcessing

            #From "synDataPlotting.jl"
    export plotDataGeneratingDist,
            plotNBMMDataGeneratingDist,plotPMMDataGeneratingDist,
            plotDataGeneratingDistAndSave!,
            plotNBMMDataGeneratingDistAndSave!,
            plotPMMDataGeneratingDistAndSave!,
            plotDataScatterAtTime,
            makeDataScatterGIF,
            plotDataScatterAndSave!,
            plot_All_Time_Points_PCA,
            plot_All_Time_Points_TSNE,
            plotGGenesDataScatterAtTime,
            plotGGenesDataScatterAndSave!,
            make_heatmapPlotatT,
            plotGGenesHeatMapAndSave!,

            #from This file 
            labelMaker,
            
            
            #From "turingChainPlotting.jl"
            plotCorrectClusteringRates,
            plotlambdaposteriorSamples,
            plotNumClusters,
            plotlambdaposteriorChainAndDensity,
            plotNumClustersOverItr,
            makeAvgPosteriorClusterMembershipGIF,
            plotAvgPosteriorClusterMembershipAtTime,
            plotAvgPosteriorClusterMembershipAtTimeStacked,
            makePosteriorClusterAssignmentTable_fromSubClusterMax,
            PosteriorPlotsAndSave!,PosteriorNumClustersPlotsAndSave!,
            PosteriorlambdaSamplesPlotsAndSave!,
            PosteriorClusteringRatePlotsAndSave!,
            PosteriorlambdaChainPlotsAndSave!,
            PosteriorlambdaPosteriorValueHistogramPlotsAndSave!,
            PosteriorlambdaRunningMeanPlotsAndSave!,
            plotAvgPosteriorClusterMembershipAndSave!,
            plotAvgPosteriorClusterMembershipAndSave2!,
            PosteriorClusterAssignmentTable_fromSubClusterMaxAndSave!,
            PosteriorMixtureWeightsAndSave!,
            plotStackedAvgPosteriorClusterMembershipAndSave!,
            plotStackedAvgPosteriorClusterMembershipAndSave2!,
            plot_MixtureWeights_post,

            #From "plotModelMetrics.jl"
            plotNMI,
            plotPosteriorClusteringNMIAndSave!,
            plotAdjustedRandIndx,
            plotNMI,plotPosteriorClusteringNMIAndSave!,
            plotAdjustedRandIndx,
            plotPosteriorClusteringARIAndSave!,
            plotVMeasure,
            plotPosteriorClusteringVMeasureAndSave!,
            plotVarInfo,
            plotPosteriorClusteringVarInfoAndSave!,
            plotTimingsSubplot,
            plotMetricsSubplot,

            gnu_stackedhistogram_plot,
            gnu_allcellshistogram_plot,
            gnuplotcluster_geneimportanceweights,
            gnu_stackedhistogram_plot,
            gnu_multiplot_1feature_type_binhistogram_plot,
            gnu_multiplot_1feature_type_binhistogram_withmeans_plot,
            gnu_1feature_all_binhistogam_plot,
            gnu_plot_1gene_data_hist,
            gnu_multiplot_1gene_data_hist,
            gnu_multiplot_stackedhistogram_plot,
            gnu_multiplot_stackedhistogram_plot2,
            gnu_multiplot_allcellshistogram_plot,
            gnu_plot_All_Time_Points_PCA,
            gnu_plot_All_Time_Points_TSNE,
            gnu_multiplot_plot_All_Time_Points_TSNE,
            gnu_multiplot_plot_All_Time_Points_PCA,
            gnu_multiplot_plot_All_Time_Points_PCA_circles,
            gnu_multiplot_plot_All_Time_Points_TSNE_circles,
            gnu_elbo_plot,
            gnu_tidy_multiplot_plot_All_Time_Points_PCA_ToyData_circles,
            gnu_multiplot_plot_All_Time_Points_PCA_4Plots_circles,
            gnu_multiplot_plot_All_Time_Points_TSNE_4Plots_circles,
            gnu_elbovsK_plot,
            gnu_ARIvsK_plot,
            gnu_tidy_multiplot_plot_All_Time_Points_TSNE_ToyData_circles,
            gnu_tidy_multiplot_plot_All_Time_Points_PCA_4Plots,
            gnu_tidy_multiplot_plot_All_Time_Points_TSNE_4Plots,
            gnu_multiplot_plot_inclusion_probabilities,
            gnu_plot_dotplot,
            gnu_plot_allocations_barchart,
            gnu_plot_memory_barchart,
            gnu_plot_avg_time_barchart,
            gnu_plot_median_time_barchart,
            gnu_plot_num_evals_barchart,
            gnu_plot_cumulative_avg_time,
            gnu_plot_cumulative_alloc,
            gnu_plot_cumulative_memory,
            gnu_plot_varying_allocations,
            gnu_plot_varying_memory,
            gnu_plot_varying_avg_time,
            gnu_plot_varying_median_time,
            gnu_plot_varying_num_evals,
            gnu_plot_inferred_signal_null_dist,
            gnu_multiplot_all_inferred_signal_null_dist,
            gnu_plot_bayesfactors,
            gnu_multiplot_grand_tour,
            gnu_ARIvsMaxIter_plot, 
            gnu_FinalElbovsMaxIter_plot,
            gnu_multiplot_elbo_vs_metrics_correlation,
            gnu_ParamVsMetric_scatterplot,
            gnu_multiplot_ParamVsMetric_scatterplot,
            gnu_plot_speedups,
            gnu_multiplot_metrics_speedups,
            gnu_multiplot_functions_speedups,


            get_tidy_tsne_transform,
            get_tsne_transform,
            get_tidy_PCA_transform,
            get_PCA_transform,
            get_conversion_units,
            reorder_colors

    function get_tidy_tsne_transform(xmat; ndims = 2,reduce_dims=50,max_iter = 1000,perplexity =20.0)
        rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())
        T = length(unique(xmat[:,1]));
        N = length(unique(xmat[:,2]));
        G = length(unique(xmat[:,3]));
        # timepoints = collect(1:T);
        # states_id = collect(1:K);
        # cell_ids = collect(1:N);
        # gene_ids = collect(1:G);
        N_t = tidy_get_Nt_from_xmat(xmat);
        x_vec = [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
        timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
        x = [x_vec[st:en] for (st,en) in timeranges]
        T = length(x)
        X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
        X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
        Y =  tsne(X_transformed, ndims, reduce_dims, max_iter, perplexity);
        return Y
    end
    function get_tsne_transform(x; ndims = 2,reduce_dims=50,max_iter = 1000,perplexity =20.0)
        rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())
        T = length(x)
        # vec(reduce(vcat,[permutedims(reduce(hcat,z[t])) for t in 1:T]))
        X = reduce(vcat,[permutedims(reduce(hcat,x[t])) for t in 1:T])#vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)    
        X_transformed =  rescale(X, dims=1);#MultivariateStats.transform(M_, data[:,2:end]')
        Y =  tsne(X_transformed, ndims, reduce_dims, max_iter, perplexity);
        return Y
    end
    
    function get_tidy_PCA_transform(xmat)
        T = length(unique(xmat[:,1]));
        N = length(unique(xmat[:,2]));
        G = length(unique(xmat[:,3]));
        # timepoints = collect(1:T);
        # states_id = collect(1:K);
        # cell_ids = collect(1:N);
        # gene_ids = collect(1:G);
        N_t = tidy_get_Nt_from_xmat(xmat);
    
    
        x_vec = [xmat[(i-1)*G+1:(i)*G,end] for i in 1:N];
        timeranges = tidy_get_timeranges(N_t)#zip(collect(0:T-1).*  N_t  .+ 1, collect(1:T).*  N_t)
        x = [x_vec[st:en] for (st,en) in timeranges]
    
    
        T = length(x)
    
        X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
        # Z = vec(vcat([permutedims(hcat(z[t]...)) for t in 1:T]...))
        data = X
        M_ = fit(PCA, data'; maxoutdim=2, pratio=1.0)
        X_transformed =  MultivariateStats.transform(M_, data')
        return X_transformed
    end
    
    function get_PCA_transform(x)
        T = length(x)
        X = vcat([permutedims(hcat(x[t]...)) for t in 1:T]...)
        data = X 
        M_ = fit(PCA, data'; maxoutdim=2, pratio=1.0)
        X_transformed =  MultivariateStats.transform(M_, data')
        return X_transformed
    end

    function get_conversion_units(units;unit_type="")
        conversion_const= 1.0
        if !isnothing(units)
            if unit_type == "time"
                if units == "sec"
                    units_str = "(s)"
                    conversion_const= 1.0*10^(-9.)
                elseif units == "min"
                    units_str = "(min)"
                    conversion_const= 1.0*10^(-9.)/60
                elseif units == "hrs"
                    units_str = "(hrs)"
                    conversion_const= 1.0*10^(-9.)/(60*60)
                elseif units == "ms"
                    units_str = "(ms)"
                    conversion_const= 1.0*10^(-6.)
                elseif units == "micros"
                    units_str = "(µs)"
                    conversion_const= 1.0*10^(-3.)
                else
                    units_str = "(ns)"
                end
            elseif unit_type == "mem"
                if units == "GiB"
                    units_str = "(GiB)"
                    conversion_const= 9.313225746154785*10^(-10.)
                elseif units == "MiB"
                    units_str = "(MiB)"
                    conversion_const= 9.5367431640625*10^(-7.)
                elseif units == "KiB"
                    units_str = "(KiB)"
                    conversion_const= 0.00097656
                else
                    units_str = "(Bytes)"
                end
            else
                units_str = ""
                conversion_const=1.0
            end
        else
            if unit_type == "mem"
                units_str = "(Bytes)"
            elseif unit_type == "time"
                units_str = "(ns)"
            else
                units_str = ""
            end
        end
    
        return units_str,conversion_const
    end

    function reorder_colors(linCol_,sorted_unique_clusters)
        num_colors = length(linCol_)
        num_clusters = length(sorted_unique_clusters)
        color_order = Vector{Int}(undef,num_clusters)
        not_used_color_order = Vector{Int}(undef,num_colors-num_clusters)
        # pca_used_color_order = collect(1:num_colors)
        curr_index_counter = IndexCounter(1)
        curr_index = getCount(curr_index_counter)
        n_curr_index_counter = IndexCounter(1)
        n_curr_index = getCount(n_curr_index_counter)
        for i in 1:num_colors
            if curr_index <= num_clusters
                if i in sorted_unique_clusters
                    location_index = findall(x->x==i, sorted_unique_clusters)[1]
                    color_order[location_index] = i
                    increment!(curr_index_counter)
                    curr_index = getCount(curr_index_counter)
                elseif n_curr_index <= num_colors-num_clusters
                    not_used_color_order[n_curr_index] = i
                    increment!(n_curr_index_counter)
                    n_curr_index = getCount(n_curr_index_counter)
                end
            elseif n_curr_index <= num_colors-num_clusters
                not_used_color_order[n_curr_index] = i
                increment!(n_curr_index_counter)
                n_curr_index = getCount(n_curr_index_counter)
            end
        end
        not_used_color_order = vcat(sort(not_used_color_order[1:num_clusters],rev=true),not_used_color_order[num_clusters+1:end])
        color_order = vcat(color_order,not_used_color_order)
        # println(color_order)
        return color_order
    end
    include(curr_dir*src_dir*"synDataPlotting.jl")
    include(curr_dir*src_dir*"turingChainPlotting.jl")
    include(curr_dir*src_dir*"plotModelMetrics.jl")
    include(curr_dir*src_dir*"plotGNUPlots.jl")
    function labelMaker(item,K)
        return reshape([ item*" "*string(i) for i in 1:K],1,K)
    end
end