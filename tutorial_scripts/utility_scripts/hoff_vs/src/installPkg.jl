using Pkg
Pkg.update()



isinstalled(pkg::AbstractString) = pkg != "METADATA" && pkg != "REQUIRE" && pkg[1] != '.' && Pkg.cd(isdir, pkg)
# 
pkgList =  ["PProf", "ForwardDiff", "DistributionsAD", "EvalMetrics", "Measurements", "MLDatasets", "OrderedCollections", "CategoricalArrays", "Clustering", "Permutations", "StatsFuns", "JLD2", "Zygote", "Expectations", "DataFrames", "SpecialFunctions", "Turing", "MCMCChains", "StatsModels", "MLJ", "Distributions", "JSON3", "AdvancedVI", "TimerOutputs", "MLJBase", "Combinatorics", "TSne", "KernelDensity", "HDF5", "Muon", "ProgressMeter", "Gadfly", "FileIO", "TypedTables", "TimeSeries", "Printf", "PrettyTables", "CSV", "Stan", "Revise", "Test", "Random", "Debugger", "BenchmarkTools", "Tulip", "Optim", "SortingAlgorithms", "StatsPlots", "ColorSchemes", "ReverseDiff", "PyCall", "DataFramesMeta", "JSON", "StatsBase", "Distances", "Flux", "Plots", "Profile", "NLSolversBase", "PyPlot", "StaticArrays", "XLSX", "RCall", "Chain", "GLM", "Hyperopt", "MLBase", "Lasso", "HypothesisTests", "JuMP", "OptimalTransport", "LaTeXStrings", "StochasticOptimalTransport", "ZygoteRules", "JuliaFormatter", "MultivariateStats", "RDatasets", "Colors", "NMF", "Gnuplot", "ClusterValidityIndices", "Hungarian", "BlackBoxOptim", "Optimisers"]
for pkg in pkgList
    if !(pkg in keys(Pkg.project().dependencies))
        if !isinstalled(pkg)
            Pkg.add(pkg)
        end
    end  
end
curr_dir  = ENV["PWD"] # = "/mnt/e/cnwizu/Playground/SCoOP-sc" #
src_dir = "/hoff_vs/src/"
env_location = curr_dir*"/nclsn/bin/python"
ENV["PYTHON"] = env_location
Pkg.build("PyCall")
Pkg.build("RCall")



using Random
using Distributions
using Flux, Turing
using Turing.Variational
using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics, LinearAlgebra, HypothesisTests,Distances
using Test
using Debugger
using CSV,DataFrames
using JSON, JSON3
using Dates
using TSne, MultivariateStats, Clustering
using LaTeXStrings, TypedTables, PrettyTables
using RCall
using Gnuplot, Colors, ColorSchemes
using SpecialFunctions
using Optim
using BenchmarkTools
using DataFramesMeta
using Profile
using JLD2,FileIO,HDF5
using Hyperopt
using OrderedCollections
using StaticArrays
using ClusterValidityIndices
using ProgressMeter
using PyCall
using Hungarian

