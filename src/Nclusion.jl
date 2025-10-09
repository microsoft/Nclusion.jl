module Nclusion

    # Write your package code here.

    using Distributions
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, LinearAlgebra, HypothesisTests,SpecialFunctions
    using Optim
    using Flux
    using DataFrames
    using JSON, JSON3
    using Dates
    using MultivariateStats, Clustering
    using OrderedCollections
    using CSV
    using LaTeXStrings, TypedTables, PrettyTables
    using JLD2,FileIO
    using Random
    using Test
    using Colors, ColorSchemes
    using BenchmarkTools
    using Profile
    using HDF5
    using ClusterValidityIndices
    using StaticArrays
    using Pkg
    using Distributed
    using Distances
    using Logging
    using Base.Threads
    using Logging,LoggingExtras
    import Debugger

    curr_dir = ENV["PWD"]
    src_dir = "/src/"

    export @name, 
        naming_vec,
    #    set_current_value,
        addToDict!,
        addToOrderedDict!,
        get_unique_time_id,
        load_data,
        preparing_data,
        select_cells_hvgs,
        initialize_model_parameters,
        run_nclusion,
        run_cavi,
        make_ids,
        summarize_parameters,
        mk_outputs_filepath,
        mk_outputs_pathname,
        saving_summary_file,
        save_embeddings,
        save_pips,
        save_Nk,
        make_nclusion_inputs,
        make_labels,
        save_labels,
        _flushed_logger,
        append_summary,
        create_results_dict,
        benchmark_nclusion,
        setup_experiment_tag,
        sort_projection_matrix_key_names,select_data_representation,
        verify_initialization_type,
        size_concentration_contractions,get_highly_variable_genes_bool,
        subset_on_highly_variable_genes_bool,center_and_scale_data_cols,
        center_and_scale_matrix_cols,size_reorder_clusters,
        calculate_s_values,
        return_occupied_clusters,
        subset_on_occupied_clusters,
        posterior_summaries,
        posterior_summaries_on_projection,
        make_new_embeddings,
        generate_report,
        submit_slurm_job_to_generate_report


    export recursive_flatten,
           outermelt, 
           innermelt,
           initialize_dict,
           initialize_ordered_dict,
           add_variables_to_ordered_dict!,
           @add_variables_to_ordered_dict_as_string!,
           @add_variables_to_ordered_dict_as_string!,get_m_init_from_initialization_approach


    export   t_test,
             norm_weights,
             normToProb,
             norm_weights3,
             norm_weights3!, 
             normToProb3!,
             sigmoidNorm!,
             ln_Gamma_distribution_normalizer,
             ln_Beta_distribution_normalizer,
             ln_Dirichlet_distribution_nomralizer

    # export getRandIndices,
    #        getNMI,
    #        getVmeasure,
    #        getVarInfo,
    #        getJaccardSimilarity,
    #        time_invariant_ari,
    #        time_invariant_nmi,
    #        time_invariant_vmeasure,
    #        time_invariant_jaccard,
    #        time_invariant_varinfo,
    #        calc_time_invariant_ARI_summarization,
    #        calc_time_variant_ARI_summarization,
    #        calc_time_invariant_NMI_summarization,
    #        calc_time_variant_NMI_summarization,
    #        calc_time_invariant_Vmeasure_summarization,
    #        calc_time_variant_Vmeasure_summarization,
    #        calc_time_invariant_VarInfo_summarization,
    #        calc_time_variant_VarInfo_summarization,
    #        calc_time_invariant_Jaccard_summarization,
    #        calc_time_variant_Jaccard_summarization,
    #        calc_time_invariant_CVI_summarization


    export cavi,
           extract_cluster_paramter,
           extract_scalars_paramter,
           extract_condition_paramter,
           extract_gene_paramter,
           extract_rtik_paramter,
           extract_elbo_vals_perK,
           extract_r_paramter,
           extract_c_paramter,
           extract_and_add_parameters_to_outputs_dict!,
   
           #Functions that calculate parts of the elbo
           calc_DataElbo,
           calc_Hz,
           calc_SurragateLowerBound_unconstrained,
           calc_SurragateLowerBound,
           calc_Hs,
           calc_wAllocationsLowerBound,
           calc_GammaElbo,
           calc_alphaElbo,
           calc_ImportanceElbo,
           calc_Hv,
           calc_DataElbo7,
           calc_DataElbo12,
           log_TrainFeature!,
           testting_convergence_monitoring,
           extract_training_features,
           check_cluster_interpretability!,
           reset_clusters!,
           redistribute_r!,
           remove_small_clusters!,
   
   
           calc_DataElbo25_fast3,
           calc_Hz_fast3,
           calc_HyjkSurragateLowerBound_unconstrained,
           calc_wAllocationsLowerBound_fast3,
           calc_alphaElbo_fast3,
           calc_HsGammaAlphaElbo_fast3,
           calc_HsElbo,
           calculate_elbo,
           calculate_elbo_perK,
           calc_SurragateLowerBound_unconstrained_elbo,
           calc_DataElbo_perK,
           calculate_elbo_mpu,
           calc_DataElbo_mpu,
           Lp_data,
           Lp_pi_tau_phi,
           Lp_tau_omega,
           Lp_pi_chi,
           Lp_surrogate,
           Lp_chi,
           Lp_omega,
           Lp_sigma_sq,
           Lp_mu,
           Lp_nu,
           Lp_lambda,
           Lp_rho,
           Lp_eta,
           Lq_r,
           Lq_c,
           Lq_d,
           Lq_w1w2,
           Lq_g1g2,
           Lq_ab,
           Lq_ms_sq_mu,
           Lq_ms_sq_nu,
           Lq_uv,
           Lq_y,
           Lq_h1h2,
           ELBO,
   
           #Math Util function
           c_Ga,
           c_Beta,
   
           #Variational Distribution Updates
           update_rtik_mpu!,
           update_Ntk_mpu!,
           update_c_ttprime_mpu!,
           update_d_hat_mpu!,
           update_d_hat_sum_mpu!,
           update_Nk_mpu!,
           update_x_hat_k_mpu!,
           update_x_hat_sq_k_mpu!,
           update_mk_hat_mpu!,
           update_v_sq_k_hat_mpu!,
           update_yjk_mpu!,
           update_var_muk_hat_mpu!,
           update_κk_hat_mpu!,
           update_σ_sq_k_hat_mpu!,
           update_λ_sq_hat_mpu!,
           update_st_hat_mpu!,
           update_Tk_mpu!,
           update_gh_hat_mpu!,
           update_ηk!,
           update_w1!,
           update_w2!,
           update_d!,
           update_d_sum!,
           update_r!,
           adjust_E_ln_pi!,
           E_log_normal_l_j!,
           update_c!,
           update_y!,
           update_h1!,
           update_h2!,
           update_u!,
           update_v!,
           update_a!,
           update_b!,
           update_m_mu!,
           update_s_sq_mu!,
           update_m_nu!,
           update_s_sq_nu!,
           g_transfromed_surrogate_lb,
           update_g1g2!,
           update_x_hat!,
           update_Ctt!,
           update_Nk!,
           update_x_hat_sq!,
           update_CNtk!,
           update_alpha_Tk_stats!,
   
   
           #Expected value functions
           βk_expected_value,
           uk_expected_value,
           logUk_expected_value,
           log1minusUk_expected_value,
           log_π_expected_value,
           log_τ_kj_expected_value,
           log_τ_k_expected_value,
           τ_μ_expected_value,
           e_ln_sigma_sq,
           e_ll_sq_diff_mu,
           e_ln_pi,
           e_ln_omega,
           e_ln_minusomega,
           e_omega,
           e_minusomega,
           e_chi,
           e_minus_chi,
           e_ln_chi,
           e_ln_minus_chi,
           e_ln_lambda,
           e_one_over_lambda,
           e_mu_sq,
           e_mu,
           e_prior_sq_diff_mu,
           e_nu_sq,
           e_nu_diff_sq,
           e_nu,
           e_prior_sq_diff_nu,
           e_ln_eta,
           e_ln_minus_eta,
           e_eta,
           e_minus_eta,
           recursive_minus_e_chi_cumprod,
           log_of_recursive_minus_e_chi_cumprod,
           expectation_sbk,
           e_one_over_sigma_sq,
           E_ln_sigma_sq,
           E_one_over_sigma_sq,
           E_ll_sq_diff_mu,
           E_ln_pi,
           E_ln_omega,
           E_omega,
           E_minusomega,
           E_ln_minusomega,
           E_chi,
           E_minus_chi,
           E_ln_chi,
           E_ln_minus_chi,
           E_ln_lambda,
           E_one_over_lambda,
           E_mu_sq,
           E_mu,
           E_prior_sq_diff_mu,
           E_nu_sq,
           E_nu,
           E_prior_sq_diff_nu,
           E_nu_diff_sq,
           E_ln_eta,
           E_ln_minus_eta,
           E_eta,
           E_minus_eta,
           recursive_minus_E_chi_cumprod,
           log_of_recursive_minus_E_chi_cumprod,
           expectation_SBk,
           recursive_cumsum_E_ln_minusomega,

           αt_expected_value,
           γ_expected_value,
           log_γ_expected_value,
           log_αt_expected_value,
           log_w_ttprime_expected_value,
           log_tilde_wt_expected_value,
           log_minus_tilde_wt_expected_value,
   
           log_τ_kj_error_expected_value,
           log_τ_k_error_expected_value,
           τ_μ_error_expected_value,
   
           adjust_e_log_π_tk3!,
           log_π_expected_value_fast3!,
           expectation_βk,
           recursive_minus_e_uk_cumprod,
           expectation_log_normal_l_j!,
           expectation_log_π_tk,
           expectation_log_tilde_wtt,
           expectation_log_minus_tilde_wtt,
           expectation_uk,
           expectation_αt,
           expectation_log_αt,
           expectation_logUk,
           expectation_log1minusUk,
   
   
           #Intializations functions
           init_params_states,
           init_θ_hat_tk,
           init_ηtkj_prior,
           
           init_mk_hat!,
           init_λ_sq_vec!,
           init_σ_sq_k_vec!,
           init_v_sq_k_hat_vec!,
           init_ghk_hat_vec!,
           init_c_ttprime_hat_vec!,
           init_d_hat_vec!,
           init_st_hat_vec!,
           init_yjk_vec!,
           init_rtik_vec!,
           initialize_VariationalInference_types!,
           low_memory_initialization,
           init_c!,
           init_r!,
           init_w1!,
           init_w2!,
           init_d!,
           init_d_k,
           init_g1!,
           init_m_mu!,
           init_s_sq_mu!,
           init_y!,
           init_h1!,
           init_h2!,
           init_a!,
           init_b!,
           init_m_nu!,
           init_s_sq_nu!,
           init_u!,
           init_v!,
           generate_fake_cells,
           generate_fake_time_conditions,
           generate_fake_clusters,
           generate_fake_model_params,
           generate_fake_data_params,
           generate_fake_inputs,

   
           #Surragate Lower Bound optimzation functions
           #from viSurragateHdpUpdates.jls
           SurragateLowerBound_util,
           SurragateLowerBound_unconstrained_util,
           g_constrained!,
           g_unconstrained!,
           genterate_Delta_mk,
   
           #Closures for Surragate Lower Bound optimzation
           #from viSurragateHdpUpdates.jl
           g_constrained_closure!,
           g_unconstrained_closure!,
           SurragateLowerBound_closure,
           SurragateLowerBound_unconstrained_closure,
   
           Features,
           check_nothing_type,
           CellFeature,
           MatrixConditionFeature,
           DataFeature,
           GeneFeatures,
           TrainFeature,
           ClusterFeature,
           ConditionFeature,
           ModelParameterFeature,
           ElboFeatures,
           get_timeranges,
           get_linear_index_as_ragged_array,
           get_linear_time_condition_update_neighbors,
           get_linear_time_condition_network_neighbors,
           _reset!


    
    include("processing.jl")
    include("math.jl")
    include("modelMetrics.jl")
    include("viCustomType.jl")
    include("viCoordinateAscent.jl")
    include("viElboCalculations.jl")
    include("viExpectations.jl")
    include("viInitializations.jl")
    include("viSurragateHdpUpdates.jl")
    include("viVariationalUpdates.jl")
end
