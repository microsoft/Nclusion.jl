

using Logging,LoggingExtras
using Random
using Distributions
using Flux
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, LinearAlgebra, HypothesisTests
using Test
using CSV,DataFrames
using JSON, JSON3
using Dates
using MultivariateStats, Clustering
using LaTeXStrings, TypedTables, PrettyTables
using Colors, ColorSchemes
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

# include(curr_dir*src_dir*"Nclusion.jl")
using Nclusion


logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;

datafilename1 = curr_dir*"/data/write/pbmc3k.h5ad" # 
alpha0 = 1 * 10^(-0.0)
gamma0 = 1 * 10^(-0.0)
KMax = 25
seed = 12345
elbo_ep = 10^(-0.0)
num_iter = 500
dataset_name = "example"
outdir = "$curr_dir"

outputs_dict = run_nclusion(datafilename1;
 logger = logger,
 num_iter = num_iter,
 outdir=outdir,
 dataset_name = dataset_name,
 seed = seed,
 elbo_ep = elbo_ep,
 alpha0 = alpha0,
 gamma0 = gamma0,
 KMax = KMax,
 )
filepath = outputs_dict[:filepath]
filename = "$filepath/output.jld2"

jldsave(filename,true;outputs_dict=outputs_dict)


