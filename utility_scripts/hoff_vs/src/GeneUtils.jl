module GeneUtils
    using Random
    using Distributions
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
    using CSV,DataFrames

    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"
    include(curr_dir*src_dir*"IOUtils.jl")
    using .IOUtils
    include(curr_dir*src_dir*"MathUtils.jl")
    using .MathUtils

    export  create_vec_cluster_scores_mat, 
            create_vec_cluster_gene_weights_mat, 
            topGenesClusters_byCalledCluster, 
            topGenesClusters_byAllCells,
            cells_in_k_clus_index,
            cells_not_in_k_clus_index,
            gene_weights_clus,
            gene_weights_not_clus,
            gene_weights_all,
            get_geneID,
            generate_topNGeneWeightValues_DF,
            generate_topNGeneWeightStats_DF,
            create_vec_cluster_gene_weightsbygene_mat,
            topGenesClustersbygene_byCalledCluster,
            topGenesClustersbygene_byAllCells,
            gene_weights_clus_tensor,
            gene_weights_not_clus_tensor,
            generate_topNGeneWeightStats_DF_tensor,

            calc_cell_normal_μ_τ_ll_scores,
            calc_cell_normal_μ_τ_l_scores,
            calc_cell_normal_μ_τ_ll_gene_weights,
            calc_cell_normal_μ_τ_ll_gene_weightsbygene,
            calculate_gene_weight_mean,
            calculate_gene_weight_variance,
            calculate_gene_weight_stdev,
            calculate_gene_weight_n,
            calc_gene_importance_weights,
            get_gene_PIP,
            get_bayes_factor,
            calculate_gene_weight_posteriorInterval,
            calculate_gene_weight_summary_statistics





    include(curr_dir*src_dir*"geneSummaries.jl")
    include(curr_dir*src_dir*"geneUtilityFunctions.jl")
end