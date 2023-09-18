module DataGenerationUtils
    export fake_mvGausssian_data_for_testing,
    fake_mvGausssian_data_largeG_indepG_for_testing,
    fake_mvGausssian_data_indepG_for_testing,
    fake_Poisson_data_for_testing,
    fake_generateGPhenNormal_RoundedCensored_indepG_for_testing,
    fake_mvGausssian_corrdata_for_testing,
    fake_mvCountsViaCopulas_corrdata_for_testing,
    generate_number_of_DEG,
    generate_DEG_locations,
    generate_DEG_values!,
    generate_nonDEG_values!,
    generate_DEG_means,
    generate_rand_covariance_matrix,

    # From "customSynDataGen.jl"
    generate1Phen_SimpleCase1,
    generateGMM,
    generateGPhen_SimpleCase1,
    generateGPhenNormal_SimpleCase1,
    check_if_truth_dict_valid,
    

    #From "turingSynDataGen.jl"
    gen1Phen_SimpleCase1_NB_TimeInvar,
    gen4gene_Corr_PGMvLN,gen1Phen_SimpleCase1_Normal_TimeInvar,
    NBMixtureModel,NormalMixtureModel,
    genGPhen_SimpleCase1_NB_TimeInvar,
    gen1Phen_SimpleCase1_Poisson_TimeInvar,
    gen1Phen_SimpleCase1_PoissonGamma_TimeInvar

    curr_dir = ENV["PWD"]
    src_dir = "/hoff_vs/src/"

    using Random
    using Distributions
    using Turing
    using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
    using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics,LinearAlgebra, Combinatorics
    using Test
    import Debugger
    using CSV,DataFrames


    include(curr_dir*src_dir*"MathUtils.jl")
    using .MathUtils
    include(curr_dir*src_dir*"TidyUtils.jl")
    using .TidyUtils
    
    include(curr_dir*src_dir*"customSynDataGen.jl")
    include(curr_dir*src_dir*"turingSynDataGen.jl")
    # include("simGenerateGaussianMixtureData.jl")
end