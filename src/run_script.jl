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


flushed_logger("Loading NCLUSION Modules....";logger)
# include(curr_dir*src_dir*"Nclusion.jl")
using Nclusion



function main(ARGS)
    datafilename1,KMax,alpha0,gamma0,phi1,phi2, kappa1, kappa2, xi1, xi2, varphi1, varphi2,nu0,sigma_sq_nu,seed,elbo_ep,num_iter,dataset_name,outdir,delta_rsum_ep, delta_rsum_lag, delta_occymean_ep, delta_occymean_lag, elbo_sign_decrease_max_tolerance, num_samples_posterior_samples, s_value_thresh, burnin = ARGS


    if !isempty(alpha0)
        alpha0 = parse(Float64, alpha0)
    else
        alpha0 = 1.0
    end
    if !isempty(gamma0)
        gamma0 = parse(Float64, gamma0)
    else
        gamma0 = 1.0
    end

    if !isempty(phi1)
        phi1 = parse(Float64, phi1)
    else
        phi1 = 1.0
    end

    if !isempty(phi2)
        phi2 = parse(Float64, phi2)
    else
        phi2 = 1.0
    end

    if !isempty(kappa1)
        kappa1 = parse(Float64, kappa1)
    else
        kappa1 = 1.0
    end

    if !isempty(kappa2)
        kappa2 = parse(Float64, kappa2)
    else
        kappa2 = 1.0
    end

    if !isempty(xi1)
        xi1 = parse(Float64, xi1)
    else
        xi1 = 1.0
    end
    if !isempty(xi2)
        xi2 = parse(Float64, xi2)
    else
        xi2 = 1.0
    end
    if !isempty(varphi1)
        varphi1 = parse(Float64, varphi1)
    else
        varphi1 = 1.0
    end
    if !isempty(varphi2)
        varphi2 = parse(Float64, varphi2)
    else
        varphi2 = 1.0
    end

    if !isempty(nu0)
        nu0 = parse(Float64, nu0)
    else
        nu0 = 0.0
    end
    if !isempty(sigma_sq_nu)
        sigma_sq_nu = parse(Float64, sigma_sq_nu)
    else
        sigma_sq_nu = 1e-12
    end
    if !isempty(KMax)
        KMax = parse(Int64, KMax)
    else
        KMax = 25
    end

    if !isempty(seed)
        seed = parse(Int64, seed)
    else
        seed = 12345
    end

    if !isempty(elbo_ep)
        elbo_ep = parse(Float64, elbo_ep)
    else
        elbo_ep = 1.0e-6
    end
    if isempty(num_iter)
        num_iter = parse(Float64, num_iter)
    else
        num_iter=500
    end
    if isempty(dataset_name)
        dataset_name="foo"
    end
    if isempty(outdir)
        outdir = "./"
    end
    if !isempty(rand_init)
        rand_init = parse(Bool, rand_init)
    else
        rand_init = true
    end
    if !isempty(delta_rsum_ep)
        delta_rsum_ep = parse(Float64, delta_rsum_ep)
    else
        delta_rsum_ep = 1.0
    end
    if !isempty(delta_rsum_lag)
        delta_rsum_lag = parse(Int64, delta_rsum_lag)
    else
        delta_rsum_lag = 2
    end
    if !isempty(delta_occymean_ep)
        delta_occymean_ep = parse(Float64, delta_occymean_ep)
    else
        delta_occymean_ep = 1.0
    end
    if !isempty(delta_occymean_lag)
        delta_occymean_lag = parse(Int64, delta_occymean_lag)
    else
        delta_occymean_lag = 2
    end
    if !isempty(elbo_sign_decrease_max_tolerance)
        elbo_sign_decrease_max_tolerance = parse(Int64, elbo_sign_decrease_max_tolerance)
    else
        elbo_sign_decrease_max_tolerance = 10
    end
    if !isempty(num_samples_posterior_samples)
        num_samples_posterior_samples = parse(Int64, num_samples_posterior_samples)
    else
        num_samples_posterior_samples = 1000
    end
    if !isempty(s_value_thresh)
        s_value_thresh = parse(Float64, s_value_thresh)
    else
        s_value_thresh=0.05
    end
    if !isempty(burnin)
        burnin = parse(Int64, burnin)
    else
        burnin = 10
    end


    outputs_dict = run_nclusion(datafilename1;
     logger = logger,
     num_iter = num_iter,
     outdir=outdir,
     dataset_name = dataset_name,
     seed = seed,
     elbo_ep = elbo_ep,
     alpha0 = alpha0,
     gamma0 = gamma0,
     phi1 = phi1,
     phi2 = phi2,
     kappa1 = kappa1,
     kappa2 = kappa2,
     xi1 = xi1,
     xi2 = xi2,
     varphi1 = varphi1,
     varphi2 = varphi2,
     nu0 = nu0,
     sigma_sq_nu = sigma_sq_nu,
     KMax = KMax,
     rand_init=rand_init,
     delta_rsum_ep=delta_rsum_ep,
     delta_rsum_lag=delta_rsum_lag,
     delta_occymean_ep=delta_occymean_ep,
     delta_occymean_lag=delta_occymean_lag,
     elbo_sign_decrease_max_tolerance=elbo_sign_decrease_max_tolerance,
     num_samples_posterior_samples=num_samples_posterior_samples,
     s_value_thresh=s_value_thresh,
     burnin=burnin,
     )
    filepath = outputs_dict[:filepath]
    filename = "$filepath/output.jld2"
    flushed_logger("Saving Outputs...";logger)
    jldsave(filename,true;outputs_dict=outputs_dict)


    flushed_logger("Finishing Script...";logger)
end
main(ARGS)
