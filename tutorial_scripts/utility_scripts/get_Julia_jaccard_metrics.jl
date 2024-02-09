ENV["GKSwstype"] = "100"

using Logging,LoggingExtras


logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;

function flushed_logger(logger,msg)
    with_logger(logger) do
        @info msg
    end
end


# @info "Using $(Threads.nthreads()) thread(s)...."
# @info "Loading Packages...."
flushed_logger(logger,"Using $(Threads.nthreads()) thread(s)....")
flushed_logger(logger,"Loading Packages....")
using Random
using Distributions
using Flux
# using Flux, Turing
# using Turing.Variational
# using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, LinearAlgebra, HypothesisTests
using Test
# import Debugger
using CSV,DataFrames
using JSON, JSON3
using Dates
using TSne, MultivariateStats, Clustering,Distances
using LaTeXStrings, TypedTables, PrettyTables
# using RCall
using Gnuplot, Colors, ColorSchemes
using SpecialFunctions
using Optim
using BenchmarkTools
using DataFramesMeta
using Profile
using JLD2,FileIO
using Hyperopt
using OrderedCollections
using HDF5
using ClusterValidityIndices
using StaticArrays
using BlackBoxOptim

using Pkg
using Distributed

curr_dir = ENV["PWD"]
src_dir = "/hoff_vs/src/"
env_location = curr_dir*"/nclsn/bin/python"
ENV["PYTHON"] = env_location
# Pkg.build("PyCall")
# Pkg.build("RCall")
# using RCall
# using PyCall

# # @info "Loading NCLUSION Modules...."
# flushed_logger(logger,"Loading NCLUSION Modules....")

# include(curr_dir*src_dir*"PlottingUtils.jl")
# using .PlottingUtils
# include(curr_dir*src_dir*"MathUtils.jl")
# using .MathUtils
# include(curr_dir*src_dir*"MCMCInferenceUtils.jl")
# using .MCMCInferenceUtils
# include(curr_dir*src_dir*"DataGenUtils.jl")
# using .DataGenerationUtils
# include(curr_dir*src_dir*"DataPreprocessingUtils.jl")
# using .DataPreprocessingUtils
# include(curr_dir*src_dir*"MetricsUtils.jl")
# using .MetricsUtils
# include(curr_dir*src_dir*"IOUtils.jl")
# using .IOUtils
# include(curr_dir*src_dir*"VariationalInferenceUtils.jl")
# using .VariationalInferenceUtils
# include(curr_dir*src_dir*"TidyUtils.jl")
# using .TidyUtils


# include(curr_dir*src_dir*"test_viFunctions.jl")
# include(curr_dir*src_dir*"baselines_RCall.jl")
# include(curr_dir*src_dir*"viDevelopment_ModelVersions.jl")

# @everywhere using Distributed
# addprocs(1)

"""
        get_timeranges(N_t)
    This function returns the linear indices that contain cells from the same condition.
    ```math

    ```
"""
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
    return [(st,en) for (st,en) in zip(starts,ends)]
end


flushed_logger(logger,"Setting Plotting Settings....")
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        # @info "Enabling Headless Plotting...."
        flushed_logger(logger,"\t Enabling Headless Plotting....")
        
        Gnuplot.options.gpviewer = false
    else
        Gnuplot.options.gpviewer = true
    end
else
    Gnuplot.options.gpviewer = true
end

if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        # @info "Setting Plotting enviroment variables...."
        flushed_logger(logger,"\t Setting Plotting enviroment variables....")
        to_display=false
    else
        to_display=true
    end
else
    to_display=true
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
    calc_Jaccard_summarization(jaccard_vov, T; conf_level=0.95)
This function takes in calculates the mean Jaccard Similarity and lower bounds of the Confidence Interval across inferred samples and for each conditions/timepoints
"""

function calc_Jaccard_summarization(jaccard_vov, condition_idx,method,hvgs,dataset,data_partition,eval_metric; conf_level=0.95)

    jaccard_values = jaccard_vov

    if data_partition == ""
        data_partition = "ALL"
    else
        data_partition = parse(Int64, data_partition)
    end

    hvgs = parse(Int64, hvgs)

    score = mean(jaccard_values)[1]
    condition_idx = parse(Int64, condition_idx)
    # condition_idx,data_partition,method,dataset,hvgs,eval_metric,score
    jaccard_summary = Matrix{Any}(undef, condition_idx, 7)
    jaccard_summary[1, 1] = condition_idx
    jaccard_summary[1, 2] = data_partition
    jaccard_summary[1, 3] = method
    jaccard_summary[1, 4] = dataset
    jaccard_summary[1, 5] = hvgs
    jaccard_summary[1, 6] = eval_metric
    jaccard_summary[1, 7] = score


    return jaccard_summary
end

function _clustering_quality_metrics(results_df,filepath_,label_column_indx,query_column_indx,condition_idx,method,hvgs,dataset,data_partition,eval_metric;beta_ = 1.0, conf_lvl = 0.95, logger=nothing)

    #Get Called assignments (col #4)
    if isnothing(label_column_indx)
        called_assignments_vec = results_df[:,end-2]
    else
        called_assignments_vec = results_df[:,label_column_indx]
    end
    

    #Get Inferred assignments (col #6 or #5)
    if isnothing(query_column_indx)
        inferred_assignments_vec = results_df[:,end]
    else
        inferred_assignments_vec = results_df[:,query_column_indx]
    end

    #Tranform data into format needed for metric functions
    conditionpoints = sort(unique(results_df[:,1]))
    conditionpoints_dict =countmap(results_df[:,1])
    N_t = [conditionpoints_dict[el] for el in conditionpoints]
    ranges_ = get_timeranges(N_t)
    called_assignments = [called_assignments_vec[st:en] for (st,en) in ranges_]
    inferred_assignments = [inferred_assignments_vec[st:en] for (st,en) in ranges_]
    T = length(conditionpoints);




    #Get Metrics on a per time/condition basis
    jaccard_vov = getJaccardSimilarity(called_assignments,[inferred_assignments])

    #Save Jaccard as CSV
    # condition_idx,data_partition,method,dataset,hvgs,eval_metric,score
    jaccard_var_summary = calc_Jaccard_summarization(jaccard_vov,condition_idx,method,hvgs,dataset,data_partition,eval_metric;conf_level=conf_lvl)
    jaccard_var_summary_df = DataFrame(jaccard_var_summary,:auto);
    output_suffix = ""
    if data_partition == ""
        data_partition = "ALL"
        output_suffix = ""
    else
        data_partition = data_partition
        output_suffix = data_partition
    end
    rename!(jaccard_var_summary_df,[:condition_idx,:data_partition,:method,:dataset,:hvgs,:eval_metric,:score])
    CSV.write(filepath_*"/"*eval_metric*"_"*method*"_"*dataset*"_"*hvgs*"hvgs"*output_suffix*".csv", jaccard_var_summary_df) 

   
end



function main(ARGS)
    filename,label_column_indx,query_column_indx,condition_idx,method,hvgs,dataset,data_partition,eval_metric,outfiliepath = ARGS

    label_column_indx = parse(Int64, label_column_indx)
    query_column_indx = parse(Int64, query_column_indx)
    


    flushed_logger(logger, "Loading Cluster Memberships")
    cluster_results_df =CSV.read(filename,DataFrame,header=true)
    _clustering_quality_metrics(cluster_results_df,outfiliepath,label_column_indx,query_column_indx,condition_idx,method,hvgs,dataset,data_partition,eval_metric)
    flushed_logger(logger, "Finished script")
end



main(ARGS)
