ENV["GKSwstype"] = "100"

using Logging,LoggingExtras
logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;

function flushed_logger(msg;logger=nothing)
    if !isnothing(logger)
        with_logger(logger) do
            @info msg
        end
    end
end
flushed_logger("Using $(Threads.nthreads()) thread(s)....";logger)
flushed_logger("Loading Packages....";logger)
using Random
using Distributions
using Flux
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, LinearAlgebra, HypothesisTests
using Test
using CSV,DataFrames
using JSON, JSON3
using Dates
using TSne, MultivariateStats, Clustering
using LaTeXStrings, TypedTables, PrettyTables
using Gnuplot, Colors, ColorSchemes
using SpecialFunctions
using Optim
using BenchmarkTools
using Profile
using JLD2,FileIO
using OrderedCollections
using HDF5
using ClusterValidityIndices
using StaticArrays
using Pkg
using Distributed
curr_dir = ENV["PWD"]
src_dir = "/src/"


flushed_logger("Loading NCLUSION Modules....";logger)
include(curr_dir*src_dir*"nclusion.jl")
using .nclusion


flushed_logger("Setting Plotting Settings....";logger)
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        flushed_logger("\t Enabling Headless Plotting....";logger)
        
        Gnuplot.options.gpviewer = false
    else
        Gnuplot.options.gpviewer = true
    end
else
    Gnuplot.options.gpviewer = true
end
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        flushed_logger("\t Setting Plotting enviroment variables....";logger)
        to_display=false
    else
        to_display=true
    end
else
    to_display=true
end


function main(ARGS)
    datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir = ARGS


    if !isempty(alpha1)
        alpha1 = parse(Float64, alpha1)
    else
        alpha1 = 1.0
    end
    if !isempty(gamma1)
        gamma1 = parse(Float64, gamma1)
    else
        gamma1 = 1.0
    end
    if !isempty(KMax)
        KMax = parse(Int64, KMax)
    else
        KMax = 10
    end

    if !isempty(seed)
        seed = parse(Int64, seed)
    else
        seed = 12345
    end

    if !isempty(elbo_ep)
        elbo_ep = parse(Float64, elbo_ep)
    else
        elbo_ep = 10^(-6)
    end
    if isempty(dataset)
        dataset = ""
    end
    if isempty(outdir)
        outdir = ""
    end

    outputs_dict = run_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = logger)
    filepath = outputs_dict[:filepath]
    filename = "$filepath/output.jld2"
    flushed_logger("Saving Outputs...";logger)
    jldsave(filename,true;outputs_dict=outputs_dict)


    flushed_logger("Finishing Script...";logger)
end
main(ARGS)
