module MCMCInferenceUtils

    using  Clustering # TSne, MultivariateStats,
    using Random
    using Distributions
    using Turing
    using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
    using Test
    import Debugger
    using CSV,DataFrames
    using LinearAlgebra
    using Hungarian

    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"
    include(curr_dir*src_dir*"DataPreprocessingUtils.jl")
    using .DataPreprocessingUtils
    # include(curr_dir*src_dir*"synDataPreprocess.jl")
    # using  .syntheticDataPreprocessing
    
            #From "customInference.jl"
    export addNewTimepoint!,
           addClusterFromGammaBaseDist!,
           addToExistingPoisGammaCluster!,
           generate_V_h_init,
           generate_θ_t,
           generate_π_h_init,
           _2ndBlockedGibbs_AllTimeStep,
           _BlockedGibbs_AllTimeStepInit,
           _3rdBlockedGibbs_AllTimeStep,
           _4thBlockedGibbs_AllTimeStep,
           _1stBlockedGibbs_AllTimeStep,
           _1stBlockedGibbs_1stTimeStep,
           _firstAttempt,
           
           #From "modelMetrics.jl"
        #    getRandIndices,
        #    getNMI,
        #    getVmeasure,
        #    getVarInfo,


           #From "turingInferenceModels.jl"
           timeseries_indep_dp_pmm1,
           timeseries_indep_dp_pmm2,
           timeseries_indep_dp_pmm2FactorialDesign,
           timeseries_indep_dp_pmm3,
           timeseries_indep_dp_pmm_MixtureModelLiklihood1,
           timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood1,
           timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood2,
           timeseries_indep_dp_pmm_MultivariateMixtureModelLiklihood22,
           timeseries_hdp_pmm_MultivariateMixtureModelLiklihood1,
           timeseries_hdp_pmm_MultivariateMixtureModelLiklihood2,
           timeseries_hdp_pmm_MultivariateMixtureModelLiklihood3,
           timeseries_hdp_gmm_MultivariateMixtureModelLiklihood1,
           timeseries_indep_dp_gmm1,_indep_dp_gmm1,_indep_dp_pmm1,
           timeseries_indep_dp_pmm1_fixedλ,infiniteTimeSeriesPMM,
           dp_pmm_sb2,
           dp_pmm_sb1,
           timeseries_Hdpsb_dp_pmm1,
           timeseries_Hdpsb_dp_pmm2, 
           timeseries_ddp_pmm, 

           #From "turingChainProcessing.jl"
           getEmpricalCorrectClusteringRates, 
           getzIDs, 
           getzDF, 
           getlambdaIDs, 
           getlambdaposteriorAvg, 
           getlambdaposteriorDF,
           getmuposteriorAvg,
           getmuposteriorDF,
           getmuIDs,
           getsigmaposteriorAvg,
           getsigmaposteriorDF,
           getsigmaIDs,
           generateDataDF_from_TuringChain,
           get_var_ids,
           get_chain_var,
           get_static_chn_param_univariate,
           get_dynamic_chn_param_univariate,
           unpack_data_chn_param,
           get_πposterior,
           get_Poisson_cell_lhprob,
           get_Normal_cell_lhprob,
           get_P_tensor,
           stephens_relabelling,
           update_nu, 
           init_nu,
           calculate_Q_hat,
           calculate_Cost_per_t,
           calculate_Cost_across_t,
           create_CostMat_across_t,
           get_average_posterior_cluster_frequency,
           get_chn_values_dict,
           relabel_chain_variable_dict,
           make_θ_t_var_names,
           make_π_t_var_names,
           make_v_var_names,
           make_z_var_names,
           make_λ_var_names,
           make_m_var_names,
           make_s_var_names,
           make_θ_t_vov2Mat,
           make_π_t_vov2Mat,
           make_v_vov2Mat,
           make_z_vov2Mat,
           make_λ_vov2Mat,
           make_m_vov2Mat,
           make_s_vov2Mat,
           make_vectorParam_vov2Mat,
           partition_matrix_cluster_df_names,
           partition_matrix_cluster_df_names2,
           _position1_Matrix_re_func,
           _position2_Matrix_re_func,
           _position2_Matrix_re_func2,
           π_t_post_vov2Mat,
           _get_πposterior,
           _get_HDP_πposterior,
           partition_cluster_df_names,
           get_cell_lhprob_mat,
           get_post_cluster_assgn,
           get_post_cluster_assgnMat2,
           get_post_cluster_assgnMat3,
           _get_chnMatrix_values_dict,
           _relabel_chain_variable_dict,
           make_PIP_t_var_names,
           make_z_t_vov2Mat,
           make_PIP_t_vov2Mat,
           make_λmat_vov2Mat,
           make_λmat_var_names,
           getLambdaInferenceModel2_newNames,
           extract,
           re_func, 
           time_re_func,
           cellID_re_func, 
           _position1_re_func,
           _position2_re_func,
           _position3_re_func,
           general_re_func,
           get_average_posterior_cluster_frequency2,
           get_clus_w_maxPopAtT,
           FindNonunique_func


    include(curr_dir*src_dir*"customInference.jl")
    include(curr_dir*src_dir*"turingInferenceModels.jl")
    include(curr_dir*src_dir*"turingChainProcessing.jl")
#     include("modelMetrics.jl")
end