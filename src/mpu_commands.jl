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
flushed_logger(logger,"Using $(Threads.nthreads()) thread(s)....")
flushed_logger(logger,"Loading Packages....")
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
using DataFramesMeta
using Profile
using JLD2,FileIO
using Hyperopt
using OrderedCollections
using HDF5
using ClusterValidityIndices
using StaticArrays
using Pkg
using Distributed
curr_dir = ENV["PWD"]
src_dir = "/src/"


flushed_logger(logger,"Loading NCLUSION Modules....")
include(curr_dir*src_dir*"nclusion.jl")
using .nclusion


flushed_logger(logger,"Setting Plotting Settings....")
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
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
        flushed_logger(logger,"\t Setting Plotting enviroment variables....")
        to_display=false
    else
        to_display=true
    end
else
    to_display=true
end


flushed_logger(logger,"Setting User Defined paramaters...")
datafilename1 = "anadata.h5ad"
alpha1 = 1.0
gamma1 = 1.0
KMax = 10
num_var_feat = 2500
scale_factor = 10000
seed = 12345
elbo_ep = 10^(-6)




flushed_logger(logger,"Loading Data and Metadata...")
Random.seed!(seed)
fid1 = h5open(datafilename1,"r")
anndata_dict1 = read(fid1)
x_mat = anndata_dict1["raw"]["X"]
N = size(x_mat)[2]



flushed_logger(logger,"Preparing Data...")
gene_ids = anndata_dict1["var"]["_index"]
cell_ids = anndata_dict1["obs"]["_index"]
cell_cluster_labels = anndata_dict1["obs"]["cell_type"]["codes"]
unique_clusters_labels =  [join(uppercasefirst.(split(el,"_")), " " ) for el in anndata_dict1["obs"]["cell_type"]["categories"]]
KCalled = length(unique_clusters_labels)
cluster_label_dict = OrderedDict(unique_clusters_labels[i] => i for i in 1:KCalled)
rev_cluster_label_dict = OrderedDict(i =>  unique_clusters_labels[i] for i in 1:KCalled)
cell_labels = [rev_cluster_label_dict[el+1] for el in cell_cluster_labels]
if all(isa.(cell_labels,Number))
    if iszero(minimum(cell_labels))
        cell_labels = cell_labels .+ 1
    end
end
cell_cluster_labels = cell_cluster_labels .+ 1
num_cnts = anndata_dict1["obs"]["total_counts"]
x = nothing
z_true = nothing
C_t = nothing
numi = nothing
top_genes_bool = nothing
top_genes = nothing

flushed_logger(logger,"Selecting HVGs...")
T=1
x_temp = Vector{Vector{Vector{Float64}}}(undef,T)
x_temp[1] = [Float64.(collect(col)) for col in eachcol(x_mat)]
numi_temp = Vector{Vector{Float64}}(undef,T)
numi_temp[1] = Float64.(collect(num_cnts))
log_norm_x = raghavan2021_lognormalization(x_temp;scaling_factor=scale_factor,pseudocount=1.0,numi=numi_temp)
log_norm_xmat =  hcat(vcat(log_norm_x...)...)
gene_std_vec = vec(std(log_norm_xmat, dims=2))
sorted_indx  = sortperm(gene_std_vec,rev=true)
topGGenes = num_var_feat
top_genes = gene_ids[sorted_indx][1:topGGenes]
top_genes_bool = [in(el,Set(top_genes)) for  el in gene_ids]
T = 1
x = Vector{Vector{Vector{Float64}}}(undef,T)
z_true = Vector{Vector{Int}}(undef,T)
C_t = Vector{Float64}(undef,T)
numi = Vector{Vector{Float64}}(undef,T)
num_cells = N
x_mat_to_use = nothing
x_mat_to_use = x_mat
chosen_cells = randperm(N)[1:num_cells]
cell_intersect = Set(chosen_cells)
cell_intersect_bool = [in(el,cell_intersect) for  el in collect(1:N)]
C_t[1] = sum(cell_intersect_bool)
x[1] = [Float64.(collect(col)) for col in eachcol(x_mat_to_use[top_genes_bool,cell_intersect_bool])]
z_true[1] = Int.(collect(cell_cluster_labels))[cell_intersect_bool]
numi[1] = Float64.(collect(num_cnts))[cell_intersect_bool]
G = length(x[1][1])


flushed_logger(logger,"Creating NCLUSION Input...")
x_input = raghavan2021_lognormalization(x;scaling_factor=scale_factor,pseudocount=1.0,numi=numi);



flushed_logger(logger,"Setting Up Model parameters...")
rand_init = false   
uniform_theta_init = true;
ep = 0.001;
L_init = 20#1#
ηk_L = LinRange(1/(G),1/10,L_init)
ηk_L = collect(ηk_L)
L = length(ηk_L)
final_elbo_vec = Vector{Float64}(undef,L)
elbo_vec = Vector{Any}(undef,L)
rtik_vec = Vector{Any}(undef,L)
yjk_vec = Vector{Any}(undef,L)
perK_elbo_vec = Vector{Any}(undef,L)
outputs_dict_vec = Vector{Any}(undef,L)
delta_t_vec=Vector{Float64}(undef,L)
nk_perL_vec = Vector{Any}(undef,L)
ηk_vec = Vector{Any}(undef,L)
num_iter, num_local_iter = 1000,1;
elbologger = [ElboFeatures(l,KMax,num_iter) for l in 1:L]
per_K_bool = true
update_η_bool = false
ηk_trend_vec_ = Vector{Any}(undef,L)

flushed_logger(logger,"Maximum K initialized at $KMax")
flushed_logger(logger,"Starting Variational Inference")
elsaped_time = @elapsed begin
    for l in 1:L
        flushed_logger(logger, "\t Model $l running now...")
        
        ηk,α0,γ0,ϕ0 = ηk_L[l],alpha1,gamma1,1.
        input_str_list = @name x_input, KMax, ηk,α0,γ0,ϕ0, elbolog num_iter, num_local_iter;
        input_key_list = Symbol.(naming_vec(input_str_list));
        input_var_list = [x_input, KMax, ηk,α0,γ0,ϕ0,elbologger[l], num_iter, num_local_iter];
        inputs_dict = OrderedDict()
        addToDict!(inputs_dict,input_key_list,input_var_list);
        st = time()
        outputs_dict = variational_inference_dynamicHDP_vshoff_mpu(inputs_dict;uniform_theta_init = true, rand_init = rand_init, update_η_bool= update_η_bool,mt_mode="optimal",elbo_ep = elbo_ep);
        dt = time() - st
        flushed_logger(logger, "\t \t Finished Training Model $l. Model $l took $dt seconds to run...")
        elbo_, rtik_, yjk_hat_, mk_hat_, v_sq_k_hat_, σ_sq_k_hat_, var_muk_, κk_hat_, Nk_, gk_hat_, hk_hat_, ak_hat_, bk_hat_,x_hat_, x_hat_sq_, d_hat_t_, c_tt_prime_, st_hat_, λ_sq_,per_k_elbo_,ηk_, Tk_, chain_dict, is_converged, truncation_value, ηk_trend_ = (; outputs_dict...);
        final_elbo_vec[l] = elbo_[end]
        perK_elbo_vec[l] = per_k_elbo_
        ηk_vec[l] = ηk_
        nk_perL_vec[l] = Nk_
        ηk_trend_vec_[l] = ηk_trend_
        elbo_vec[l] = elbo_
        rtik_vec[l] = rtik_
        yjk_vec[l] = yjk_hat_
        delta_t_vec[l] = dt
        outputs_dict_vec[l] = outputs_dict

    end
end

flushed_logger(logger, "$L Models Took a total of $elsaped_time seconds to run")
flushed_logger(logger,"\t Averaging Parameters....")
importance_weights = norm_weights(final_elbo_vec)
weighted_yjk_vec = importance_weights .* yjk_vec;
weighted_rtik_vec = importance_weights .* rtik_vec;
weighted_elbo_vec = importance_weights .* elbo_vec;
weighted_elapsed_time = importance_weights .* delta_t_vec;
mk_hat_L = [el[:mk_hat_] for el in outputs_dict_vec];
v_sq_k_hat_L = [el[:v_sq_k_hat_] for el in outputs_dict_vec];
σ_sq_k_hat_L = [el[:σ_sq_k_hat_] for el in outputs_dict_vec];
Nk_L = [el[:Nk_] for el in outputs_dict_vec];
d_hat_t_L = [el[:d_hat_t_] for el in outputs_dict_vec];
c_tt_prime_L = [el[:c_tt_prime_] for el in outputs_dict_vec];
st_hat_L = [el[:st_hat_] for el in outputs_dict_vec];
λ_sq_L = [el[:λ_sq_] for el in outputs_dict_vec];
weighted_mk_vec = importance_weights .* mk_hat_L;
weighted_v_sq_k_vec = importance_weights .* v_sq_k_hat_L;
weighted_σ_sq_k_vec = importance_weights .* σ_sq_k_hat_L;
weighted_Nk_vec = importance_weights .* Nk_L;
weighted_d_vec = importance_weights .* d_hat_t_L;
weighted_c_tt_prime_vec = importance_weights .* c_tt_prime_L;
weighted_st_vec = importance_weights .* st_hat_L;
weighted_λ_sq_vec = importance_weights .* λ_sq_L;
old_weighted_elbo_vec = deepcopy(weighted_elbo_vec)
maxLen = maximum(length.(old_weighted_elbo_vec))
for l in 1:L
    curr_len = length(old_weighted_elbo_vec[l])
    diff_ = maxLen - curr_len
    if !iszero(diff_)
        weighted_elbo_vec[l] = vcat(old_weighted_elbo_vec[l],zeros(diff_))
    end
end
mean_elbo = sum(weighted_elbo_vec)
pip = sum(weighted_yjk_vec)
mean_rtik = sum(weighted_rtik_vec)
mean_mk= sum(weighted_mk_vec)
mean_v_sq_k = sum(weighted_v_sq_k_vec)
mean_σ_sq_k = sum(weighted_σ_sq_k_vec)
mean_Nk = sum(weighted_Nk_vec)
mean_d = sum(weighted_d_vec)
mean_c_tt_prime = sum(weighted_c_tt_prime_vec)
mean_st = sum(weighted_st_vec)
mean_λ_sq = sum(weighted_λ_sq_vec) 

flushed_logger(logger, "Calculating Model Fit")
z_argmax = [argmax.(el) for el in  mean_rtik];
z_infer = vcat(z_argmax...);
num_posterior_samples = 10;
z_post_s = [[z_infer]]
ari_vov, _, _, _= getRandIndices(z_true, z_post_s)
flushed_logger(logger, "ARI $(mean(mean(ari_vov)))")
flushed_logger(logger, "Finishing Script...")
