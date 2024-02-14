ENV["GKSwstype"] = "100"

using Logging,LoggingExtras
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

include(curr_dir*src_dir*"Nclusion.jl")
using .Nclusion


logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;

datafilename1 = "/users/cnwizu/data/cnwizu/SCoOP-sc/data/pdac-biopsy/5000hvgs_pdac_biopsy_preprocessed2.h5ad" # 
alpha1 = 1 * 10^(-7.0)
gamma1 = 1 * 10^(-7.0)
KMax = 25
seed = 12345
elbo_ep = 10^(-0.0)
num_iter = 500
dataset = "pdac_biopsy"
outdir = "$curr_dir"
save_metrics=true

outputs_dict = run_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = logger,num_iter = num_iter,save_metrics=save_metrics)
filepath = outputs_dict[:filepath]
filename = "$filepath/output.jld2"

jldsave(filename,true;outputs_dict=outputs_dict)


