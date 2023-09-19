module MetricsUtils


    using  Clustering # TSne, MultivariateStats,
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
    using Distances
#     using Revise


    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"
    include(curr_dir*src_dir*"MathUtils.jl")
    using .MathUtils
    include(curr_dir*src_dir*"TidyUtils.jl")
    using .TidyUtils

    
           #From "modelMetrics.jl"
    export getRandIndices,
           getNMI,
           getVmeasure,
           getVarInfo,
           getJaccardSimilarity,
           getMetrics_and_Stats,
           time_invariant_ari,
           time_invariant_nmi,
           time_invariant_vmeasure,
           time_invariant_jaccard,
           time_invariant_varinfo,
           calc_time_invariant_ARI_summarization,
           calc_time_variant_ARI_summarization,
           calc_time_invariant_NMI_summarization,
           calc_time_variant_NMI_summarization,
           calc_time_invariant_Vmeasure_summarization,
           calc_time_variant_Vmeasure_summarization,
           calc_time_invariant_VarInfo_summarization,
           calc_time_variant_VarInfo_summarization,
           calc_time_invariant_Jaccard_summarization,
           calc_time_variant_Jaccard_summarization,
           calc_time_invariant_CVI_summarization,
           setup_metrics_list,

           
           tidy_getRandIndices,
           tidy_calc_time_invariant_ARI_summarization,
           tidy_calc_time_variant_ARI_summarization,
           tidy_getNMI,
           tidy_calc_time_invariant_NMI_summarization,
           tidy_calc_time_variant_NMI_summarization,
           tidy_getVmeasure,
           tidy_calc_time_invariant_Vmeasure_summarization,
           tidy_calc_time_variant_Vmeasure_summarization,
           tidy_getVarInfo,


           tidy_calculate_elapsed_iterations,
           tidy_calculate_approx_cluster_utilization,
           tidy_calculate_overall_cluster_occupancy_rate,
           tidy_calculate_per_time_approx_cluster_utilization,
           tidy_calculate_per_time_cluster_occupancy_rate

    include(curr_dir*src_dir*"modelMetrics.jl")
    include(curr_dir*src_dir*"viAlgMetrics.jl")
end