using Random
using Distributions
using Turing
using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
using Test
import Debugger
using CSV,DataFrames

Random.seed!(1234)

# First, Generate data 

# 1-gene NB model, 3 clusters, 3 different means, same variances, parameters do not change over time
# 4 time point, same number of cells per time point, equal proportion of cells come from each cluster

k = 3
T = 4
Cₜ =  2*ones(Int,T) #[2,2,2,2]
μ = [2.,4.,6.]
σ² = 10.0*ones(k)

@model function gen1Phen_SimpleCase1_NB_TimeInvar(x,T,k,C_t,μ,σ²,α,mixing_prob)
    if x === missing
        x = Vector{Vector}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = tzeros(Int, cells)
        end
    end
    if μ >= σ²
        throw(DomainError(μ, "μ must be less than σ²"))
    end
    if k !== maximum([length(μ), length(σ²)])
        error("must ha μ and σ² for each of the k components")
    end
    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²)
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end
    # mixing_prob ~ Dirichlet(α .* ones(k) ./ k)
    H = NegativeBinomial.(r,p)
    z = Vector{Vector}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z
    # data_dist = MixtureModel(NegativeBinomial.(r,p),mixing_prob)

    
end
C_t = 6 .* ones(Int,T)
α = 2
model = gen1Phen_SimpleCase1_NB_TimeInvar(missing,T,k,C_t,μ,σ²,α,ones(k) ./ k)
c = sample(model, Prior(), 1)
DataFrame(c)


function generate1Phen_SimpleCase1(;k = 3, T = 4, Cₜ = 2*ones(Int,T),μ = [2.,4.,6.], σ² = 10.0*ones(k)  )
    truth_dict = Dict{Tuple,Float64}()
    data_dict = Dict{Tuple,VecOrMat}()
    r = μ.^2 ./ (σ² .- μ)
    p = r ./ (r .+ μ)
    for t in 1:T
        c = Cₜ[t]
        phen = rand.(NegativeBinomial.(r,p),c)
        phen = vcat(phen...)
        cluster_assgn = [i for i in 1:k for j in 1:c]
        concat_data = hcat(cluster_assgn,phen)
        concat_data = concat_data[shuffle(1:end),:]
        cluster_assgn, phen = concat_data[:,1], concat_data[:,2:end]
        cell_ids = [ i for i in 1:c*k]
        for i in 1:length(cell_ids)
            id = cell_ids[i]
            truth_dict[(id,t)] = cluster_assgn[id]
            data_dict[(id,t)] = phen[id,:]
        end
    end
    return truth_dict, data_dict
end

function check_if_truth_dict_valid(truth_dict,k,T,cellsInCluster)
    unique_clusters = unique([truth_dict[key] for key in keys(truth_dict)])
    unique_clusters = sort(unique_clusters)
    true_clusters =  collect(1:k .* ones(eltype(unique_clusters)))
    unique_clusters == true_clusters # Checks if number of clusters is correct
    if unique_clusters == true_clusters
        counts_chckr = Vector{Bool}()
        for cluster in 1:k
            true_number_of_cells_in_cluster = sum([cellsInCluster[t] for t in 1:T]) 
            gen_number_of_cells_in_clusters = sum([truth_dict[key]==k for key in keys(truth_dict)])
            push!(counts_chckr, true_number_of_cells_in_cluster == gen_number_of_cells_in_clusters)
        end
        if  all(counts_chckr)
            return true
        else
            return "counts within clusters is NOT correct"
        end
    else
        return "number of clusters is NOT correct"
    end
end

# TESTING THE ABOVE FUNCTION
cellsInCluster = Cₜ
test_truth_dict, test_data_dict = generate1Phen_SimpleCase1()
check_if_truth_dict_valid(test_truth_dict,k,T,cellsInCluster)

#GENERATING SYNTHETIC DATA
truth_dict, data_dict = generate1Phen_SimpleCase1()

function StickBreaking(β_prime)
    one_minus_β_prime = 1. .- β_prime
    cumprod_one_minus_β_prime = cumprod(one_minus_β_prime, dims = 1)
    β_prime_copy = deepcopy(β_prime)
    push!(β_prime_copy, 1.0)
    pushfirst!(cumprod_one_minus_β_prime, 1.0)
    β = β_prime_copy .* cumprod_one_minus_β_prime
    return β
end


function _firstAttempt(data_dict)
    DP_α = 1
    a,b = 1,1
    lambda_prior_dist = Gamma(a,b) # alpha = a = shape, b = rate = 1/scale = 1/beta
    # λ = rand(lambda_prior_dist)
    likelihood_vec = Vector()
    H = Poisson
    # push!(likelihood_vec,H(λ))

    # For Time Point 1!
    tp = 1

    key_array_tp1 = collect(keys(data_dict))
    key_array_tp1 = key_array_tp1[last.(key_array_tp1) .== tp]
    key_array_tp1 = sort(key_array_tp1, by = first)
    _tp1 = [(first(key),data_dict[key]) for key in key_array_tp1]

    tp1_cluster_assgn = Dict()
    n_ht = Dict(t => Dict() for t in 1:T) # Number of cells in cluster h in this current (they use iteration to refer to simulation iterations (like in Gibbs) but I have to consider time point too!)
    itr = 1
    curr_k_max = 0
    for cell in _tp1
        if itr == 1
            # println("here")
            tp1_cluster_assgn[1] = [cell]
            n_ht[tp][1] = 1
            curr_k_max += 1
            # println("Current K_max = $curr_k_max")
            λ = rand(Gamma(lambda_prior_dist.α + last(cell)[1], 1 + 1/(lambda_prior_dist.θ))) # Only works for One gene case!
            push!(likelihood_vec,H(λ))
        else
            # println("now here")
            # println("Current K_max = $curr_k_max")
            phen_val = last(cell)[1]
            L_vals = [pdf(L,phen_val) for L in likelihood_vec]
            n_h = [n_ht[tp][key] for key in sort(collect(keys(n_ht[tp])))]
            cluster_prob = L_vals .* n_h
            ###### NOW IM GUESSING ######
            new_cluster_prob = DP_α / (itr - 1 + DP_α)
            new_cluster_indicator = rand(Bernoulli(new_cluster_prob))
            if isone(new_cluster_indicator)
                curr_k_max += 1 # Increase total Number
                # println("WE BIG NOW")
                # println("Current K_max = $curr_k_max")
                
                tp1_cluster_assgn[curr_k_max] = [cell] # Add to Cluster
                n_ht[tp][curr_k_max] = 1 # Increase Cluster Count
                λ = rand(Gamma(lambda_prior_dist.α + last(cell)[1], 1 + 1/(lambda_prior_dist.θ))) #rand(lambda_prior_dist) ## Sample new parameters
                push!(likelihood_vec,H(λ))
            else
                norm_cluster_prob =  softmax(cluster_prob)#cluster_prob ./sum(cluster_prob)
                cluster_val = rand(Categorical(norm_cluster_prob)) # Determine Cluster
                push!(tp1_cluster_assgn[cluster_val],cell)# Add to Cluster
                n_ht[tp][cluster_val] += 1 # Increase Cluster Count
                n_ = n_ht[tp][cluster_val] 
                suff_stat = sum(last.(tp1_cluster_assgn[cluster_val]))[1]
                λ = rand(Gamma(lambda_prior_dist.α + suff_stat, n_ + 1/(lambda_prior_dist.θ)))  ## Sample new parameters
                likelihood_vec[cluster_val] = H(λ)
            end
        end
        # println("We outside")
        # println("Current K_max = $curr_k_max")
        itr += 1
    end
    # println("After Time 1, The current K is $curr_k_max")
    tp_cluster_assgn = [tp1_cluster_assgn]
    _tp_vec = [_tp1]
    ka_= collect(keys(data_dict))
    cells_in_t = [length(_tp1)]
    for tp in 2:T
        key_array_ = ka_[last.(ka_) .== tp]
        key_array_ = sort(key_array_, by = first)
        _tp = [(first(key),data_dict[key]) for key in key_array_]
        push!(_tp_vec,_tp) #Just keeping track for dubuging purposes
        push!(cells_in_t,length(_tp))
        new_tp_cluster_assgn = Dict()
        for i in 1:curr_k_max
            n_ht[tp][i] = 0
            new_tp_cluster_assgn[i] = []
        end
        
        push!(tp_cluster_assgn,new_tp_cluster_assgn)
        itr = 1
        for cell in _tp
            phen_val = last(cell)[1]
            L_vals = [pdf(L,phen_val) for L in likelihood_vec]
            n_h = [n_ht[tp][key] for key in sort(collect(keys(n_ht[tp])))] .+ [n_ht[tp-1][key] for key in sort(collect(keys(n_ht[tp-1])))]
            cluster_prob = L_vals .* n_h
            ###### NOW IM GUESSING ######
            new_cluster_prob = DP_α / (cells_in_t[tp-1]  +itr - 1 + DP_α)
            new_cluster_indicator = rand(Bernoulli(new_cluster_prob))
            if isone(new_cluster_indicator)
                curr_k_max += 1 # Increase total Number
                tp_cluster_assgn[tp][curr_k_max]= [cell]
                n_ht[tp][curr_k_max] = 1 # Increase Cluster Count
                for i in 1:tp-1
                    n_ht[i][curr_k_max] = 0 # Update Kmax for past time points
                    tp_cluster_assgn[i][curr_k_max] = []  # Update Kmax for past time points
                end
                λ = rand(Gamma(lambda_prior_dist.α + last(cell)[1], 1 + 1/(lambda_prior_dist.θ))) #rand(lambda_prior_dist) ## Sample new parameters
                push!(likelihood_vec,H(λ))
            else
                norm_cluster_prob = softmax(cluster_prob)#cluster_prob ./sum(cluster_prob)
                cluster_val = rand(Categorical(norm_cluster_prob)) # Determine Cluster
                push!(tp_cluster_assgn[tp][cluster_val],cell)# Add to Cluster
                n_ht[tp][cluster_val] += 1 # Increase Cluster Count
                n_ = n_ht[tp][cluster_val] + n_ht[tp-1][cluster_val]
                old_suff_stat_bool = isempty(last.(tp_cluster_assgn[tp-1][cluster_val]))
                if old_suff_stat_bool
                    old_suff_stat = 0 
                else
                    old_suff_stat = sum(last.(tp_cluster_assgn[tp-1][cluster_val]))[1]
                end
                new_suff_stat = sum(last.(tp_cluster_assgn[tp][cluster_val]))[1] 
                suff_stat = new_suff_stat + old_suff_stat
                λ = rand(Gamma(lambda_prior_dist.α + suff_stat, n_ + 1/(lambda_prior_dist.θ)))  ## Sample new parameters
                likelihood_vec[cluster_val] = H(λ)
            end


            itr +=1
        end
        # println("After Time $tp, The current K is $curr_k_max")
    end
    return curr_k_max,n_ht,tp_cluster_assgn,likelihood_vec
end
# Debugger.@enter f(data_dict)
sim_time = 1000
k_max_simtime = Vector{Int}(undef,sim_time)
for s in 1:sim_time
    k_max_simtime[s],_,_,_ =_firstAttempt(data_dict)
end
# f(data_dict)
histogram(k_max_simtime)

function getCustomersAtTimepoint(data_dict,tp)
    ka_ = collect(keys(data_dict))
    key_array_ = ka_[last.(ka_) .== tp]
    key_array_ = sort(key_array_, by = first)
    _tp = [(first(key),data_dict[key]) for key in key_array_]
    return _tp
end
getCustomersAtTimepoint(data_dict,2)
getCustomersAtTimepoint(data_dict,1)



## Cluster CHOSEN Does not EXISTS
# First TP and itr = 1


# First TP and itr =/= 1


# NOT first TP and itr shouldnt matter (right?) ...



function addNewTimepoint!(prev_tp,
                        tp_cluster_assgn,
                        n_ht)
    curr_tp = prev_tp + 1
    prev_keys1 = collect(keys(tp_cluster_assgn[prev_tp]))
    prev_keys2 = collect(keys(n_ht[prev_tp]))
    for  key in prev_keys1
        tp_cluster_assgn[curr_tp][key] = []
    end
    for key in prev_keys2
        n_ht[curr_tp][key] =  0
    end
    # tp_cluster_assgn[curr_tp] = Dict(tp_cluster_assgn[curr_tp][key] => [] for key in prev_keys1) 
    # n_ht[curr_tp] = Dict(n_ht[curr_tp][key] => 0 for key in prev_keys2)
end


function addClusterFromGammaBaseDist!(customer, 
                                            tp, 
                                            curr_k_max, 
                                            tp_cluster_assgn, 
                                            n_ht,
                                            likelihood_vec,
                                            lambda_prior_dist;
                                            ft=false, H=Poisson )

    curr_k_max += 1 # Increase total Number
    tp_cluster_assgn[tp][curr_k_max] = [customer] # Add to Cluster
    n_ht[tp][curr_k_max] = 1   # Increase Cluster Count
    if !ft
        for i in 1:tp-1
            n_ht[i][curr_k_max] = 0 # Update Kmax for past time points
            tp_cluster_assgn[i][curr_k_max] = []  # Update Kmax for past time points
        end
    end
    suff_stat = last(customer)[1]
    n_ = 1
    shape_update = lambda_prior_dist.α + suff_stat
    scale_update = 1/(lambda_prior_dist.θ) +  n_
    λ = rand(Gamma(shape_update,scale_update))
    push!(likelihood_vec,H(λ))
    return curr_k_max
end



function test_addClusterFromGammaBaseDist!_init(data_dict)
    test_val_dict = Dict()
    test_val_dict["customers_tp1"] = getCustomersAtTimepoint(data_dict,1)
    test_val_dict["customers_tp2"] = getCustomersAtTimepoint(data_dict,2)
    test_val_dict["tp_cluster_assgn"] = [Dict() for i in 1:4]
    test_val_dict["n_ht"] = [Dict() for i in 1:4]
    test_val_dict["likelihood_vec"] = []
    test_val_dict["lambda_prior_dist"] = Gamma(1,1)
    return  test_val_dict #test_customers1, test_customers2,test_curr_k_max,test_tp_cluster_assgn,test_n_ht,test_likelihood_vec,test_lambda_prior_dist
end


function test_addClusterFromGammaBaseDist!_(data_dict)
    # test_customers1, test_customers2,test_curr_k_max,test_tp_cluster_assgn,test_n_ht,test_likelihood_vec,test_lambda_prior_dist  =  test_addClusterFromGammaBaseDist!_init(data_dict)
    addCluster_test_val = test_addClusterFromGammaBaseDist!_init(data_dict)
    
    # MAKE SURE FIRST CASE IS OK: First TP and itr = 1
    tp = 1
    curr_k_max = 0
    itr = 1
    curr_k_max = addClusterFromGammaBaseDist!(addCluster_test_val["customers_tp1"][itr],
                            tp,
                            curr_k_max, 
                            addCluster_test_val["tp_cluster_assgn"], 
                            addCluster_test_val["n_ht"],
                            addCluster_test_val["likelihood_vec"],
                            addCluster_test_val["lambda_prior_dist"];
                            ft=true, H=Poisson )

    #TODO: Make tests
    @testset "MAKE SURE FIRST CASE IS OK: First TP and itr = 1" begin
        @test curr_k_max == 1
        @test haskey(addCluster_test_val["tp_cluster_assgn"][1],1)
        @test !haskey(addCluster_test_val["tp_cluster_assgn"][1], 2)
        @test addCluster_test_val["tp_cluster_assgn"][1][1] == [(1,data_dict[(1,1)])]
        @test length(collect(keys(addCluster_test_val["tp_cluster_assgn"][1]))) ==   1
        @test haskey(addCluster_test_val["n_ht"][1],1)
        @test !haskey(addCluster_test_val["n_ht"][1],2)
        @test addCluster_test_val["n_ht"][1][1] == 1
        @test length(collect(keys(addCluster_test_val["n_ht"][1]))) ==   1
        @test length(addCluster_test_val["likelihood_vec"]) == 1
    end

    # MAKE SURE SECOND  CASE IS OK: First TP and itr =/= 1
     
    tp = 1
    itr = 2
    curr_k_max = addClusterFromGammaBaseDist!(addCluster_test_val["customers_tp1"][itr],
                            tp,
                            curr_k_max, 
                            addCluster_test_val["tp_cluster_assgn"], 
                            addCluster_test_val["n_ht"],
                            addCluster_test_val["likelihood_vec"],
                            addCluster_test_val["lambda_prior_dist"];
                            ft=true, H=Poisson )
    
    #TODO: Make tests
    @testset "MAKE SURE SECOND  CASE IS OK: First TP and itr =/= 1" begin
        @test curr_k_max == 2
        @test haskey(addCluster_test_val["tp_cluster_assgn"][1],2)
        @test !haskey(addCluster_test_val["tp_cluster_assgn"][1], 3)
        @test addCluster_test_val["tp_cluster_assgn"][1][2] == [(2,data_dict[(2,1)])]
        @test length(collect(keys(addCluster_test_val["tp_cluster_assgn"][1]))) ==   2
        @test haskey(addCluster_test_val["n_ht"][1],2)
        @test !haskey(addCluster_test_val["n_ht"][1],3)
        @test addCluster_test_val["n_ht"][1][2] == 1
        @test length(collect(keys(addCluster_test_val["n_ht"][1]))) ==   2
        @test length(addCluster_test_val["likelihood_vec"]) == 2
    end
    # 
    addNewTimepoint!(tp,addCluster_test_val["tp_cluster_assgn"],addCluster_test_val["n_ht"])
    @testset "MAKE SURE NEW TIME POINT IS ADDED:" begin
        @test curr_k_max == 2
        @test haskey(addCluster_test_val["tp_cluster_assgn"][2],1)
        @test haskey(addCluster_test_val["tp_cluster_assgn"][2],2)
        @test isempty(addCluster_test_val["tp_cluster_assgn"][2][1])
        @test isempty(addCluster_test_val["tp_cluster_assgn"][2][2])
        @test length(collect(keys(addCluster_test_val["tp_cluster_assgn"][2]))) ==   2
        @test haskey(addCluster_test_val["n_ht"][2],1)
        @test haskey(addCluster_test_val["n_ht"][2],2)
        @test addCluster_test_val["n_ht"][2][1] == 0
        @test addCluster_test_val["n_ht"][2][2] == 0
        @test length(collect(keys(addCluster_test_val["n_ht"][2]))) ==   2
    end
    # MAKE SURE THIRD  CASE IS OK: NOT first TP and itr shouldnt matter (right?) ...
    tp = 2
    itr = 1
    curr_k_max =  addClusterFromGammaBaseDist!(addCluster_test_val["customers_tp2"][itr],
                            tp,
                            curr_k_max, 
                            addCluster_test_val["tp_cluster_assgn"], 
                            addCluster_test_val["n_ht"],
                            addCluster_test_val["likelihood_vec"],
                            addCluster_test_val["lambda_prior_dist"];
                            ft=false, H=Poisson )

    #TODO: Make tests
    @testset "MAKE SURE THIRD  CASE IS OK: NOT first TP and itr shouldnt matter (right?)" begin
        @test curr_k_max == 3
        @test haskey(addCluster_test_val["tp_cluster_assgn"][2],3)
        @test haskey(addCluster_test_val["tp_cluster_assgn"][1],3)
        @test isempty(addCluster_test_val["tp_cluster_assgn"][1][3])
        @test addCluster_test_val["tp_cluster_assgn"][2][3] == [(1,data_dict[(1,2)])]
        @test length(collect(keys(addCluster_test_val["tp_cluster_assgn"][1]))) ==   3
        @test haskey(addCluster_test_val["n_ht"][2],3)
        @test haskey(addCluster_test_val["n_ht"][1],3)
        @test addCluster_test_val["n_ht"][2][3] == 1
        @test addCluster_test_val["n_ht"][1][3] == 0
        @test length(collect(keys(addCluster_test_val["n_ht"][1]))) ==   3
        @test length(addCluster_test_val["likelihood_vec"]) == 3
    end
    return
end


test_addClusterFromGammaBaseDist!_(data_dict)


# addNewTimepoint!(tp,updateCluster_test_val["tp_cluster_assgn"],updateCluster_test_val["n_ht"])

## Cluster CHOSEN EXISTS
# First TP 

# NOT first TP 


function addToExistingPoisGammaCluster!(customer,
                        cluster_prob,
                        tp,
                        tp_cluster_assgn,
                        n_ht,
                        likelihood_vec,
                        lambda_prior_dist;
                        ft=false, H=Poisson, is_cluster_prob_norm = false 
                        )
    
    if is_cluster_prob_norm
        norm_cluster_prob = cluster_prob
    else
        norm_cluster_prob = softmax(cluster_prob)#cluster_prob ./sum(cluster_prob)
    end
    cluster_val = rand(Categorical(norm_cluster_prob)) # Determine Cluster
    push!(tp_cluster_assgn[tp][cluster_val],customer)# Add to Cluster
    n_ht[tp][cluster_val] += 1 # Increase Cluster Count
    if  !ft
        prev_t_suff_stat_bool = isempty(last.(tp_cluster_assgn[tp-1][cluster_val]))
        if prev_t_suff_stat_bool
            prev_t_suff_stat = 0 
        else
            prev_t_suff_stat = sum(last.(tp_cluster_assgn[tp-1][cluster_val]))[1]
        end
        prev_t_n_ = n_ht[tp-1][cluster_val]
    else
        prev_t_suff_stat = 0 
        prev_t_n_ = 0 
    end
    curr_t_suff_stat = sum(last.(tp_cluster_assgn[tp][cluster_val]))[1] 
    curr_t_n_ = n_ht[tp][cluster_val]
    suff_stat = curr_t_suff_stat + prev_t_suff_stat
    n_ = curr_t_n_ + prev_t_n_
    λ = rand(Gamma(lambda_prior_dist.α + suff_stat, n_ + 1/(lambda_prior_dist.θ)))  ## Sample new parameters
    likelihood_vec[cluster_val] = H(λ)
end

function test_addToExistingPoisGammaCluster!_(data_dict)
    updateCluster_test_val = test_addClusterFromGammaBaseDist!_init(data_dict)
    tp = 1
    curr_k_max = 0
    itr = 1
    curr_k_max = addClusterFromGammaBaseDist!(updateCluster_test_val["customers_tp1"][itr],
                            tp,
                            curr_k_max, 
                            updateCluster_test_val["tp_cluster_assgn"], 
                            updateCluster_test_val["n_ht"],
                            updateCluster_test_val["likelihood_vec"],
                            updateCluster_test_val["lambda_prior_dist"];
                            ft=true, H=Poisson )
    itr += 1
    curr_k_max = addClusterFromGammaBaseDist!(updateCluster_test_val["customers_tp1"][itr],
                            tp,
                            curr_k_max, 
                            updateCluster_test_val["tp_cluster_assgn"], 
                            updateCluster_test_val["n_ht"],
                            updateCluster_test_val["likelihood_vec"],
                            updateCluster_test_val["lambda_prior_dist"];
                            ft=true, H=Poisson )
    cluster_prob = [1.0,0.0]
    itr += 1
    addToExistingPoisGammaCluster!(updateCluster_test_val["customers_tp1"][itr],
                                    cluster_prob,
                                    tp,
                                    updateCluster_test_val["tp_cluster_assgn"], 
                                    updateCluster_test_val["n_ht"],
                                    updateCluster_test_val["likelihood_vec"],
                                    updateCluster_test_val["lambda_prior_dist"];
                                    ft=true, H=Poisson,is_cluster_prob_norm = true )
    #TODO: Make tests
    @testset "MAKE SURE FIRST CASE IS OK: First TP" begin
        @test curr_k_max == 2
        @test itr == 3
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][1],1)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][1],2)
        @test !haskey(updateCluster_test_val["tp_cluster_assgn"][1], 3)
        @test length(updateCluster_test_val["tp_cluster_assgn"][1][1]) ==   2
        @test length(collect(keys(updateCluster_test_val["tp_cluster_assgn"][1]))) ==   2
        @test updateCluster_test_val["tp_cluster_assgn"][1][1] == [(1,data_dict[(1,1)]), (3,data_dict[(3,1)])]
        @test updateCluster_test_val["tp_cluster_assgn"][1][2] == [(2,data_dict[(2,1)])]
        @test haskey(updateCluster_test_val["n_ht"][1],1)
        @test haskey(updateCluster_test_val["n_ht"][1],2)
        @test !haskey(updateCluster_test_val["n_ht"][1],3)
        @test updateCluster_test_val["n_ht"][1][1] == 2
        @test updateCluster_test_val["n_ht"][1][2] == 1
        @test length(collect(keys(updateCluster_test_val["n_ht"][1]))) ==   2
        @test length(updateCluster_test_val["likelihood_vec"]) == 2
    end
    addNewTimepoint!(tp,updateCluster_test_val["tp_cluster_assgn"],updateCluster_test_val["n_ht"])
    @testset "MAKE SURE NEW TIME POINT IS ADDED:" begin
        @test curr_k_max == 2
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][2],1)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][2],2)
        @test isempty(updateCluster_test_val["tp_cluster_assgn"][2][1])
        @test isempty(updateCluster_test_val["tp_cluster_assgn"][2][2])
        @test length(collect(keys(updateCluster_test_val["tp_cluster_assgn"][2]))) ==   2
        @test haskey(updateCluster_test_val["n_ht"][2],1)
        @test haskey(updateCluster_test_val["n_ht"][2],2)
        @test updateCluster_test_val["n_ht"][2][1] == 0
        @test updateCluster_test_val["n_ht"][2][2] == 0
        @test length(collect(keys(updateCluster_test_val["n_ht"][2]))) ==   2
    end


    tp = 2
    itr = 1
    curr_k_max =  addClusterFromGammaBaseDist!(updateCluster_test_val["customers_tp2"][itr],
                            tp,
                            curr_k_max, 
                            updateCluster_test_val["tp_cluster_assgn"], 
                            updateCluster_test_val["n_ht"],
                            updateCluster_test_val["likelihood_vec"],
                            updateCluster_test_val["lambda_prior_dist"];
                            ft=false, H=Poisson )
    itr += 1
    curr_k_max =  addClusterFromGammaBaseDist!(updateCluster_test_val["customers_tp2"][itr],
                            tp,
                            curr_k_max, 
                            updateCluster_test_val["tp_cluster_assgn"], 
                            updateCluster_test_val["n_ht"],
                            updateCluster_test_val["likelihood_vec"],
                            updateCluster_test_val["lambda_prior_dist"];
                            ft=false, H=Poisson )
    
    itr += 1
    cluster_prob = [1.0,0.0,0.0,0.0]
    addToExistingPoisGammaCluster!(updateCluster_test_val["customers_tp2"][itr], #3
                                    cluster_prob,
                                    tp,
                                    updateCluster_test_val["tp_cluster_assgn"], 
                                    updateCluster_test_val["n_ht"],
                                    updateCluster_test_val["likelihood_vec"],
                                    updateCluster_test_val["lambda_prior_dist"];
                                    ft=false, H=Poisson, is_cluster_prob_norm = true  )
    cluster_prob = [0.0,0.0,1.0,0.0]
    itr += 1
    addToExistingPoisGammaCluster!(updateCluster_test_val["customers_tp2"][itr], #4
                                    cluster_prob,
                                    tp,
                                    updateCluster_test_val["tp_cluster_assgn"], 
                                    updateCluster_test_val["n_ht"],
                                    updateCluster_test_val["likelihood_vec"],
                                    updateCluster_test_val["lambda_prior_dist"];
                                    ft=false, H=Poisson, is_cluster_prob_norm = true  )
    @testset "MAKE SURE SECOND CASE IS OK: NOT first TP" begin
        @test curr_k_max == 4
        @test itr == 4
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][2],3)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][2],4)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][2],1)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][2],2)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][1],3)
        @test haskey(updateCluster_test_val["tp_cluster_assgn"][1],4)
        @test isempty(updateCluster_test_val["tp_cluster_assgn"][1][3])
        @test isempty(updateCluster_test_val["tp_cluster_assgn"][1][4])
        @test isempty(updateCluster_test_val["tp_cluster_assgn"][2][2])
        @test length(updateCluster_test_val["tp_cluster_assgn"][2][1]) ==   1
        @test length(updateCluster_test_val["tp_cluster_assgn"][2][3]) ==   2
        @test length(updateCluster_test_val["tp_cluster_assgn"][2][4]) ==   1
        @test length(collect(keys(updateCluster_test_val["tp_cluster_assgn"][1]))) ==  4 
        @test length(collect(keys(updateCluster_test_val["tp_cluster_assgn"][2]))) ==  4 
        @test updateCluster_test_val["tp_cluster_assgn"][2][1] == [(3,data_dict[(3,2)])]
        @test updateCluster_test_val["tp_cluster_assgn"][2][3] == [(1,data_dict[(1,2)]), (4,data_dict[(4,2)])]
        @test updateCluster_test_val["tp_cluster_assgn"][2][4] == [(2,data_dict[(2,2)])]
        @test haskey(updateCluster_test_val["n_ht"][2],3)
        @test haskey(updateCluster_test_val["n_ht"][2],4)
        @test haskey(updateCluster_test_val["n_ht"][2],1)
        @test haskey(updateCluster_test_val["n_ht"][2],2)
        @test haskey(updateCluster_test_val["n_ht"][1],3)
        @test haskey(updateCluster_test_val["n_ht"][1],4)
        @test updateCluster_test_val["n_ht"][2][1] == 1
        @test updateCluster_test_val["n_ht"][2][2] == 0
        @test updateCluster_test_val["n_ht"][2][3] == 2
        @test updateCluster_test_val["n_ht"][2][4] == 1
        @test length(collect(keys(updateCluster_test_val["n_ht"][1]))) ==   4
        @test length(collect(keys(updateCluster_test_val["n_ht"][2]))) ==   4
        @test length(updateCluster_test_val["likelihood_vec"]) == 4

        # @test haskey(addCluster_test_val["tp_cluster_assgn"][2],3)
        # @test haskey(addCluster_test_val["tp_cluster_assgn"][1],3)
        # @test isempty(addCluster_test_val["tp_cluster_assgn"][1][3])
        # @test addCluster_test_val["tp_cluster_assgn"][2][3] == [(1,data_dict[(1,2)])]
        # @test length(collect(keys(addCluster_test_val["tp_cluster_assgn"][1]))) ==   3
        # @test haskey(addCluster_test_val["n_ht"][2],3)
        # @test haskey(addCluster_test_val["n_ht"][1],3)
        # @test addCluster_test_val["n_ht"][2][3] == 1
        # @test length(collect(keys(addCluster_test_val["n_ht"][1]))) ==   3
        # @test length(addCluster_test_val["likelihood_vec"]) == 3
    end
end
test_addToExistingPoisGammaCluster!_(data_dict)


# for tp in 1:T
#     if tp == 1
#         ft_bool = true
#     else
#         ft_bool = false
#     end
#     customers = getCustomersAtTimepoint(data_dict,tp)
#     itr = 1
#     for customer in customers
#         if itr == 1 && ft_bool
#             addClusterFromGammaBaseDist!(customer,
#                                         tp,
#                                         curr_k_max,
#                                         )
#         else

#         end
#     end

# end

# addClusterFromGammaBaseDist!(addCluster_test_val["customers_tp1"][itr],tp,curr_k_max, addCluster_test_val["tp_cluster_assgn"], addCluster_test_val["n_ht"],addCluster_test_val["likelihood_vec"],addCluster_test_val["lambda_prior_dist"];ft=ft_bool)


function seqCRP(customers,tp_cluster_assgn,n_ht, curr_k_max;ft=false)
    itr = 1

end
function _secondAttempt()
    
end




# Just the first Time step!
function _1stBlockedGibbs_1stTimeStep(data_dict)
    K_max = 25
    gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
    gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
    rand_gibbs_init_cluster = [ rand(1:K_max,c) for c in 3 .* Cₜ]
    H = Poisson
    a,b = 1,1
    a_θ, b_θ = 1,0.1
    theta_prior = Gamma(a_θ,1/b_θ)
    θ_t =  [rand(theta_prior) for t in 1:T]
    likelihood_prior = Gamma(a,1/b)
    likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
    likelihood_vec_init = deepcopy(likelihood_vec)
    likelihood_vec = [likelihood_vec]
    V_h_init = Vector{Vector}(undef,T)
    π_h_init = Vector{Vector}(undef,T)
    for t in 1:T
        V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
        π_h_init[t] = stickbreak(V_h_init[t])
    end
    tp = 1
    cells_ = getCustomersAtTimepoint(data_dict,tp)
    phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
    cluster_prob = [π_h_init[tp] .* phen for phen in phen_liklihood]
    norm_cluster_prob = softmax.(cluster_prob)
    cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
    for c in 1:length(cluster_assignment)
        assgn = cluster_assignment[c]
        cell = cells_[c]
        push!(gibbs_tp_cluster_assgn[tp][assgn],cell)
        gibbs_n_ht[tp][assgn]+=1
    end

    #Update π
    n_ = [gibbs_n_ht[tp][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
    V_h_1 =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
    π_h_1 = stickbreak(V_h_1[1:K_max-1])

    #Update λ
    curr_n_h = [gibbs_n_ht[tp][key] for key in 1:K_max]
    prev_n_h = zeros(Int,K_max) # [gibbs_n_ht[tp-1][key] for key in 1:K_max]
    n_h = curr_n_h + prev_n_h
    get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h]) 
    curr_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
    prev_suff_stats = zeros(eltype(last(gibbs_tp_cluster_assgn[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
    # [cells_[1] for h in 1:K_max]
    suff_stats = curr_suff_stats + prev_suff_stats
    upd_a = likelihood_prior.α .+ suff_stats
    upd_b = likelihood_prior.θ .+ n_h
    upd_λ =  rand.(Gamma.(upd_a, upd_b))
    push!(likelihood_vec, H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ])

    # Update θ 
    max_K_at_t = K_max #maximum(cluster_assignment)
    new_θ_t = deepcopy(θ_t)
    new_θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_h_1[1:max_K_at_t]))))
end

_1stBlockedGibbs_1stTimeStep(data_dict)


function _1stBlockedGibbs_AllTimeStep(data_dict)
    ###### IN LOOP FORM ###########
    K_max = 25
    gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
    gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
    rand_gibbs_init_cluster = [ rand(1:K_max,c) for c in 3 .* Cₜ]
    H = Poisson
    a,b = 1,1
    a_θ, b_θ = 1,0.1
    theta_prior = Gamma(a_θ,1/b_θ)
    θ_t =  [rand(theta_prior) for t in 1:T]
    likelihood_prior = Gamma(a,1/b)
    likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
    likelihood_vec_init = deepcopy(likelihood_vec)
    likelihood_vec = [likelihood_vec]
    V_h_init = Vector{Vector}(undef,T)
    π_h_init = Vector{Vector}(undef,T)
    for t in 1:T
        V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
        π_h_init[t] = stickbreak(V_h_init[t])
    end

    for tp in 1:T
        cells_ = getCustomersAtTimepoint(data_dict,tp)

        # Sample States
        phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]

        cluster_prob = [π_h_init[tp] .* phen for phen in phen_liklihood]
        norm_cluster_prob = softmax.(cluster_prob)
        cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
        for c in 1:length(cluster_assignment)
            assgn = cluster_assignment[c]
            cell = cells_[c]
            push!(gibbs_tp_cluster_assgn[tp][assgn],cell)
            gibbs_n_ht[tp][assgn]+=1
        end
        
        #Update π
        if tp == 1
            n_ = [gibbs_n_ht[tp][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
            V_ht =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht = stickbreak(V_ht[1:K_max-1])
        else
            n_ = [gibbs_n_ht[tp][key] for key in 1:K_max] .+ [gibbs_n_ht[tp-1][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
            V_ht =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht = stickbreak(V_ht[1:K_max-1])
        end

        
        
        #Update λ
        get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
        curr_n_h = [gibbs_n_ht[tp][key] for key in 1:K_max]
        curr_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
        if tp == 1
            prev_n_h = zeros(Int,K_max) # [gibbs_n_ht[tp-1][key] for key in 1:K_max]
            prev_suff_stats = zeros(eltype(last(gibbs_tp_cluster_assgn[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
        else
            prev_n_h =  [gibbs_n_ht[tp-1][key] for key in 1:K_max]
            prev_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
        end
        n_h = curr_n_h + prev_n_h 
        suff_stats = curr_suff_stats + prev_suff_stats
        upd_a = likelihood_prior.α .+ suff_stats
        upd_b = likelihood_prior.θ .+ n_h
        upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
        push!(likelihood_vec, H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
        
        # Update θ_t 
        max_K_at_t = K_max #maximum(cluster_assignment)
        # new_θ_t = deepcopy(θ_t)
        θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[1:max_K_at_t]))))
    end
    return θ_t
end
vv = _1stBlockedGibbs_AllTimeStep(data_dict) 

function _2ndBlockedGibbs_AllTimeStep(itr,data_dict,V_ht,π_ht,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t;K_max = 5,H = Poisson,a = 1,b=1,a_θ =1, b_θ =1/5)
    ###### IN LOOP FORM ###########
    if maximum(unique([last(key) for key in keys(data_dict)])) == length(unique([last(key) for key in keys(data_dict)]))
        T = length(unique([last(key) for key in keys(data_dict)]))
    else
        T = maximum(unique([last(key) for key in keys(data_dict)])) 
    end

    tp_cluster_assgn_t = [Dict(h => [] for h in 1:K_max) for i in 1:T]
    n_ht_t = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
    # rand_gibbs_init_cluster = [ rand(1:K_max,c) for c in 3 .* Cₜ]
    

    theta_prior = Gamma(a_θ,1/b_θ)
    θ_t =  [rand(theta_prior) for t in 1:T]
    likelihood_prior = Gamma(a,1/b)
    likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
    # likelihood_vec_init = deepcopy(likelihood_vec)
    


    likelihood_vec_t = [likelihood_vec]

    cluster_assignment_t = Vector{Vector}(undef,T)
    norm_cluster_prob_t = Vector{Vector}(undef,T)
    
    for tp in 1:T
        cells_ = getCustomersAtTimepoint(data_dict,tp)

        # Sample States
        phen_liklihood = [pdf.(likelihood_vec_t[tp],last(cell)[1]) for cell in cells_]
        old_π_ht = π_ht[tp]
        cluster_prob = [old_π_ht .* phen for phen in phen_liklihood]
        norm_cluster_prob = softmax.(cluster_prob)
        norm_cluster_prob_t[tp] = norm_cluster_prob
        cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
        
        cluster_assignment_t[tp] = cluster_assignment 
        for c in 1:length(cluster_assignment)
            assgn = cluster_assignment[c]
            cell = cells_[c]
            push!(tp_cluster_assgn_t[tp][assgn],cell) # tp_cluster_assgn_t
            n_ht_t[tp][assgn]+=1 # n_ht_t
        end
        
        #Update π
        if tp == 1
            n_ = [n_ht_t[tp][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        else
            n_ = [n_ht_t[tp][key] for key in 1:K_max] .+ [n_ht_t[tp-1][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        end

        
        
        #Update λ
        get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
        curr_n_h = [n_ht_t[tp][key] for key in 1:K_max]
        curr_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        if tp == 1
            prev_n_h = zeros(Int,K_max) # [n_ht_t[tp-1][key] for key in 1:K_max]
            prev_suff_stats = zeros(eltype(last(tp_cluster_assgn_t[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        else
            prev_n_h =  [n_ht_t[tp-1][key] for key in 1:K_max]
            prev_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        end
        n_h = curr_n_h + prev_n_h 
        suff_stats = curr_suff_stats + prev_suff_stats
        upd_a = likelihood_prior.α .+ suff_stats
        upd_b = likelihood_prior.θ .+ n_h
        upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
        # likelihood_vec_t[tp] = H.(upd_λ)
        likelihood_vec = H.(upd_λ)
        push!(likelihood_vec_t, likelihood_vec)#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
        
        # Update θ_t 
        ####################
        #TODO: THIS Update yeilds a θ value in the  
        # max_K_at_t = K_max #maximum(cluster_assignment)
        # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
        θ_t[tp] = θ_t[tp]
    end
    gibbs_tp_cluster_assgn[itr] = tp_cluster_assgn_t
    gibbs_n_ht[itr] = n_ht_t
    gibbs_likelihood_vec_t[itr] = likelihood_vec_t
    gibbs_V_ht[itr] = V_ht
    gibbs_π_ht[itr] = π_ht
    gibbs_cluster_assignment[itr] = cluster_assignment_t
    gibbs_norm_cluster_prob[itr] = norm_cluster_prob_t
    gibbs_θ_t[itr] = θ_t
    return θ_t
end



K_max = 5
function generate_θ_t(theta_prior,T)
    # theta_prior = Gamma(a_θ,1/b_θ)
    θ_t =  [rand(theta_prior) for t in 1:T]
    return θ_t
end
function generate_V_h_init(T,K_max,θ_t)
    V_h_init = Vector{Vector}(undef,T)
    for t in 1:T
        V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
    end
    return V_h_init
end
function generate_π_h_init(T,K_max,θ_t)
    π_h_init = Vector{Vector}(undef,T)
    V_h_init = generate_V_h_init(T,K_max,θ_t)
    for t in 1:T
        π_h_init[t] = stickbreak(V_h_init[t])
    end
    return π_h_init, V_h_init
end

a_θ,b_θ  = 0.1,0.1
theta_prior = Gamma(a_θ,1/b_θ)
θ_t = generate_θ_t(theta_prior,T)
gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,1)
gibbs_n_ht = Vector{Vector{Dict}}(undef,1)
gibbs_likelihood_vec_t = Vector{Vector}(undef,1)
gibbs_V_ht= Vector{Vector}(undef,1)
gibbs_π_ht = Vector{Vector}(undef,1)
gibbs_cluster_assignment= Vector{Vector}(undef,1)
gibbs_norm_cluster_prob = Vector{Vector}(undef,1)
gibbs_θ_t = Vector{Vector}(undef,1)
π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
ss = _2ndBlockedGibbs_AllTimeStep(1,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t) 
ss2 = _2ndBlockedGibbs_AllTimeStep(1,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht, gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t, b_θ =1/25)


test_itr = 10000
gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,test_itr)
gibbs_n_ht = Vector{Vector{Dict}}(undef,test_itr) 
gibbs_likelihood_vec_t = Vector{Vector}(undef,test_itr)
gibbs_V_ht= Vector{Vector}(undef,test_itr)
gibbs_π_ht = Vector{Vector}(undef,test_itr)
gibbs_cluster_assignment= Vector{Vector}(undef,test_itr)
gibbs_norm_cluster_prob = Vector{Vector}(undef,test_itr)
gibbs_θ_t = Vector{Vector}(undef,test_itr)
π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
for k in 1:test_itr
    _2ndBlockedGibbs_AllTimeStep(k,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t)
end

# likelihood_vec
# V_h
# π_h

# cluster_assignment
# norm_cluster_prob
# θ_t

# gibbs_tp_cluster_assgn
# gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
# gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
# n_ht_dict
K_max = 25
numBurnin = 10000
maxNumItr = 10000


# V_h_init = Vector{Vector}(undef,T)
# π_h_init = Vector{Vector}(undef,T)

# for t in 1:T
#     V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
#     π_h_init[t] = stickbreak(V_h_init[t])
# end

function _BlockedGibbs_AllTimeStepInit(data_dict;K_max = 5,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = 1,b=1,a_θ =0.1, b_θ =1/10)
    gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
    gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
    theta_prior = theta_prior_dist(a_θ,1/b_θ)
    θ_t =  [rand(theta_prior) for t in 1:T]
    likelihood_prior =  likelihood_prior_dist(a,1/b)
    likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
    likelihood_vec = [likelihood_vec]
    # V_h_init = Vector{Vector}(undef,T)
    # π_h_init = Vector{Vector}(undef,T)
    # for t in 1:T
    #     V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
    #     π_h_init[t] = stickbreak(V_h_init[t])
    # end
    π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
    V_ht = Vector{Vector}(undef,T)
    π_ht = Vector{Vector}(undef,T)

    for tp in 1:T
        cells_ = getCustomersAtTimepoint(data_dict,tp)

        # Sample States
        phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]

        cluster_prob = [π_h_init[tp] .* phen for phen in phen_liklihood]
        norm_cluster_prob = softmax.(cluster_prob)
        cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
        for c in 1:length(cluster_assignment)
            assgn = cluster_assignment[c]
            cell = cells_[c]
            push!(gibbs_tp_cluster_assgn[tp][assgn],cell)
            gibbs_n_ht[tp][assgn]+=1
        end
        
        #Update π
        if tp == 1
            n_ = [gibbs_n_ht[tp][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        else
            n_ = [gibbs_n_ht[tp][key] for key in 1:K_max] .+ [gibbs_n_ht[tp-1][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        end

        
        
        #Update λ
        get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
        curr_n_h = [gibbs_n_ht[tp][key] for key in 1:K_max]
        curr_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
        if tp == 1
            prev_n_h = zeros(Int,K_max) # [gibbs_n_ht[tp-1][key] for key in 1:K_max]
            prev_suff_stats = zeros(eltype(last(gibbs_tp_cluster_assgn[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
        else
            prev_n_h =  [gibbs_n_ht[tp-1][key] for key in 1:K_max]
            prev_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
        end
        n_h = curr_n_h + prev_n_h 
        suff_stats = curr_suff_stats + prev_suff_stats
        upd_a = likelihood_prior.α .+ suff_stats
        upd_b = likelihood_prior.θ .+ n_h
        upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
        push!(likelihood_vec, H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
        
        # Update θ_t 
        max_K_at_t = K_max #maximum(cluster_assignment)
        # new_θ_t = deepcopy(θ_t)
        # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
        θ_t[tp] = θ_t[tp]
    end
    return θ_t, likelihood_vec, π_ht, V_ht,likelihood_prior,theta_prior 
end

_BlockedGibbs_AllTimeStepInit(data_dict;K_max = 5,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = 1,b=1,a_θ =1, b_θ =1/5)

function _3rdBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t;K_max = 5,H = Poisson,a = 1,b=1,a_θ =1, b_θ =1/5)
    ###### IN LOOP FORM ###########
    if maximum(unique([last(key) for key in keys(data_dict)])) == length(unique([last(key) for key in keys(data_dict)]))
        T = length(unique([last(key) for key in keys(data_dict)]))
    else
        T = maximum(unique([last(key) for key in keys(data_dict)])) 
    end

    #Updates from Previous iteration
    θ_t =  deepcopy(gibbs_θ_t[itr-1])
    likelihood_vec = deepcopy(gibbs_likelihood_vec_t[itr-1])
    π_ht = deepcopy(gibbs_π_ht[itr-1])
    V_ht = deepcopy(gibbs_V_ht[itr-1])
    
    #iteration specific values initialization
    tp_cluster_assgn_t = [Dict(h => [] for h in 1:K_max) for i in 1:T]
    n_ht_t = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
    likelihood_vec_t = []
    cluster_assignment_t = Vector{Vector}(undef,T)
    norm_cluster_prob_t = Vector{Vector}(undef,T)
    
    for tp in 1:T
        cells_ = getCustomersAtTimepoint(data_dict,tp)

        # Sample States
        phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
        old_π_ht = π_ht[tp]
        cluster_prob = [old_π_ht .* phen for phen in phen_liklihood]
        norm_cluster_prob = softmax.(cluster_prob)
        norm_cluster_prob_t[tp] = norm_cluster_prob
        cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
        
        cluster_assignment_t[tp] = cluster_assignment 
        for c in 1:length(cluster_assignment)
            assgn = cluster_assignment[c]
            cell = cells_[c]
            push!(tp_cluster_assgn_t[tp][assgn],cell) # tp_cluster_assgn_t
            n_ht_t[tp][assgn]+=1 # n_ht_t
        end
        
        #Update π
        if tp == 1
            n_ = [n_ht_t[tp][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        else
            n_ = [n_ht_t[tp][key] for key in 1:K_max] .+ [n_ht_t[tp-1][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        end

        
        
        #Update λ
        get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
        curr_n_h = [n_ht_t[tp][key] for key in 1:K_max]
        curr_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        if tp == 1
            prev_n_h = zeros(Int,K_max) # [n_ht_t[tp-1][key] for key in 1:K_max]
            prev_suff_stats = zeros(eltype(last(tp_cluster_assgn_t[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        else
            prev_n_h =  [n_ht_t[tp-1][key] for key in 1:K_max]
            prev_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        end
        n_h = curr_n_h + prev_n_h 
        suff_stats = curr_suff_stats + prev_suff_stats
        upd_a = likelihood_prior.α .+ suff_stats
        upd_b = likelihood_prior.θ .+ n_h
        upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
        # likelihood_vec_t[tp] = H.(upd_λ)
        # likelihood_vec = H.(upd_λ)
        push!(likelihood_vec_t,H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
        
        # Update θ_t 
        ####################
        #TODO: THIS Update yeilds a θ value in the  
        # max_K_at_t = K_max #maximum(cluster_assignment)
        # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
        θ_t[tp] = θ_t[tp]
    end
    gibbs_tp_cluster_assgn[itr] = tp_cluster_assgn_t
    gibbs_n_ht[itr] = n_ht_t
    gibbs_likelihood_vec_t[itr] = likelihood_vec_t
    gibbs_V_ht[itr] = V_ht
    gibbs_π_ht[itr] = π_ht
    gibbs_cluster_assignment[itr] = cluster_assignment_t
    gibbs_norm_cluster_prob[itr] = norm_cluster_prob_t
    gibbs_θ_t[itr] = θ_t
end

# V_h only depends of the prev. time step; λ only depends on the current time step
function _4thBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t;K_max = 5,H = Poisson,a = 1,b=1,a_θ =1, b_θ =1/5)
    ###### IN LOOP FORM ###########
    if maximum(unique([last(key) for key in keys(data_dict)])) == length(unique([last(key) for key in keys(data_dict)]))
        T = length(unique([last(key) for key in keys(data_dict)]))
    else
        T = maximum(unique([last(key) for key in keys(data_dict)])) 
    end

    #Updates from Previous iteration
    θ_t =  deepcopy(gibbs_θ_t[itr-1])
    likelihood_vec = deepcopy(gibbs_likelihood_vec_t[itr-1])
    π_ht = deepcopy(gibbs_π_ht[itr-1])
    V_ht = deepcopy(gibbs_V_ht[itr-1])
    
    #iteration specific values initialization
    tp_cluster_assgn_t = [Dict(h => [] for h in 1:K_max) for i in 1:T]
    n_ht_t = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
    likelihood_vec_t = []
    cluster_assignment_t = Vector{Vector}(undef,T)
    norm_cluster_prob_t = Vector{Vector}(undef,T)
    
    for tp in 1:T
        cells_ = getCustomersAtTimepoint(data_dict,tp)

        # Sample States
        phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
        old_π_ht = π_ht[tp]
        cluster_prob = [old_π_ht .* phen for phen in phen_liklihood]
        norm_cluster_prob = softmax.(cluster_prob)
        norm_cluster_prob_t[tp] = norm_cluster_prob
        cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
        
        cluster_assignment_t[tp] = cluster_assignment 
        for c in 1:length(cluster_assignment)
            assgn = cluster_assignment[c]
            cell = cells_[c]
            push!(tp_cluster_assgn_t[tp][assgn],cell) # tp_cluster_assgn_t
            n_ht_t[tp][assgn]+=1 # n_ht_t
        end
        
        #Update π
        if tp == 1
            n_ = [n_ht_t[tp][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        else
            n_ = [n_ht_t[tp-1][key] for key in 1:K_max] #.+ [n_ht_t[tp][key] for key in 1:K_max] #sort(collect(keys(n_ht_t[tp])))
            V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
            π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
        end

        
        
        #Update λ
        get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
        curr_n_h = [n_ht_t[tp][key] for key in 1:K_max]
        curr_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        if tp == 1
            prev_n_h = zeros(Int,K_max) # [n_ht_t[tp-1][key] for key in 1:K_max]
            prev_suff_stats = zeros(eltype(last(tp_cluster_assgn_t[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        else
            prev_n_h =  [n_ht_t[tp-1][key] for key in 1:K_max]
            prev_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
        end
        n_h = curr_n_h + prev_n_h 
        suff_stats = curr_suff_stats #+ prev_suff_stats
        upd_a = likelihood_prior.α .+ suff_stats
        upd_b = likelihood_prior.θ .+ n_h
        upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
        # likelihood_vec_t[tp] = H.(upd_λ)
        # likelihood_vec = H.(upd_λ)
        push!(likelihood_vec_t,H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
        
        # Update θ_t 
        ####################
        #TODO: THIS Update yeilds a θ value in the  
        # max_K_at_t = K_max #maximum(cluster_assignment)
        # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
        θ_t[tp] = θ_t[tp]
    end
    gibbs_tp_cluster_assgn[itr] = tp_cluster_assgn_t
    gibbs_n_ht[itr] = n_ht_t
    gibbs_likelihood_vec_t[itr] = likelihood_vec_t
    gibbs_V_ht[itr] = V_ht
    gibbs_π_ht[itr] = π_ht
    gibbs_cluster_assignment[itr] = cluster_assignment_t
    gibbs_norm_cluster_prob[itr] = norm_cluster_prob_t
    gibbs_θ_t[itr] = θ_t
    return θ_t
end

K_max = 5
a = 6
b= 1
a_θ = 1
b_θ =1/5
numBurnin = 10000
maxNumItr = 10000
maxNumItr = maxNumItr + numBurnin

gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,maxNumItr+1)
gibbs_n_ht = Vector{Vector{Dict}}(undef,maxNumItr+1) 
gibbs_likelihood_vec_t = Vector{Vector}(undef,maxNumItr+1)
gibbs_V_ht= Vector{Vector}(undef,maxNumItr+1)
gibbs_π_ht = Vector{Vector}(undef,maxNumItr+1)
gibbs_cluster_assignment= Vector{Vector}(undef,maxNumItr+1)
gibbs_norm_cluster_prob = Vector{Vector}(undef,maxNumItr+1)
gibbs_θ_t = Vector{Vector}(undef,maxNumItr+1)



rand_init_θ_t, rand_init_likelihood_vec, rand_init_π_ht, rand_init_V_ht,likelihood_prior,theta_prior  = _BlockedGibbs_AllTimeStepInit(data_dict;K_max = K_max,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = a,b=b,a_θ =a_θ, b_θ =b_θ)

gibbs_likelihood_vec_t[1] = rand_init_likelihood_vec
gibbs_V_ht[1] =  rand_init_V_ht
gibbs_π_ht[1] = rand_init_π_ht
gibbs_θ_t[1] = 0.001ones(T)#rand_init_θ_t

# _3rdBlockedGibbs_AllTimeStep(2,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,K_max = K_max)

for step in 2:maxNumItr+1
    itr = step
    _3rdBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,a = a,b=b,K_max = K_max)
    # _4thBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,a = a,b=b,K_max = K_max)
end

function _1stBlockedGibbs()
    
end


# DP PMM model under stick-breaking construction
@model dp_pmm_sb1(x, T, K) = begin
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=6
    b=1
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, K-1)
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            end
        end
    end

    
    
end



x = [[data_dict[(i,t)][1] for i in 1:6] for t in 1:4]
true_clusters = [[truth_dict[(i,t)] for i in 1:6] for t in 1:4]
coordinates = [[(i,t) for i in 1:6] for t in 1:4]
K = 5
T = 4
m=dp_pmm_sb1(x, T, K)

# Set random seed for reproducibility
Random.seed!(0);

# Compile time approx. 32s.
# Run time approx. 70s.

@time hmc_chain = begin
    burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
    n_samples = 500
    iterations = burn + n_samples
    n_components = 10
    stepsize = 0.01
    nleapfrog = floor(Int, 1 / stepsize)
 
    chain = sample(dp_pmm_sb1(x, T, K), 
           HMC(stepsize, nleapfrog),
           iterations)
end
;

###### REPEATED BELOW##########
@model dp_pmm_sb2(x, T, K) = begin
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=6
    b=1
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, K-1)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                z[t][c] ~ Categorical(eta)
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                z[t][c] ~ Categorical(eta)
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
            end
        end
        # k[t] = k_t
    end

    
    
end
burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
n_samples = 500
iterations = burn + n_samples
K = 25
T = 4
m2=dp_pmm_sb2(x, T, K)
# pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))
tchain = sample(m2, PG(100),iterations);
k = map(
    t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :k), :].value))),
    1:iterations
);
ids = findall(map(name -> occursin("λ", string(name)), names(tchain)));
p = plot(tchain[:, ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter)
###############################################
###############################################
###############################################


function jitter!(a::Array, factor=1.0)
    @assert eltype(a) <: AbstractFloat
    a .+= rand(size(a)...) .* factor
end

@model function infiniteTimeSeriesPMM(x,T)
    # Hyper-parameters, i.e. concentration parameter and parameters of H.

    # α = 1.0
    # μ0 = 0.0
    # σ0 = 1.0
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)

    a=6
    b=1
    
    


    rpm = DirichletProcess.(θ_t)
    
    # Define the base distribution, i.e. expected value of the Dirichlet process.
    H = Gamma(a,1/b)
    
    # Latent assignment.
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
        
    # Locations of the infinitely many clusters.
    λ = tzeros(Float64,0)
    nk = Vector{Vector{Int}}(undef, T)
    for t in 1:T
        for c in 1:nobs_t[t]
            # Number of clusters.
            K = maximum(vcat(z[1:t]...))
            nk[t] = Vector{Int}(map(k -> sum(z[t] .== k), 1:K))

            # Draw the latent assignment.
            if t==1
                n__ = nk[t]#z[t][c] ~ ChineseRestaurantProcess(rpm[t], nk[t])
            else
                n__ =  nk[t] .+ nk[t-1]#z[t][c] ~ ChineseRestaurantProcess(rpm[t], nk[t] .+ nk[t-1])
            end
            z[t][c] ~ ChineseRestaurantProcess(rpm[t], n__)
            
        
            # Create a new cluster?
            if z[t][c] > K
                if t != 1
                    for i in 1:t-1
                        push!(nk[i],zero(eltype(nk[i])))
                    end
                end
                push!(λ, 0.0)

                # Draw location of new cluster.
                # println("Here at $K")
                λ[z[t][c]] ~ H
            end
                    
            # Draw observation.
            x[t][c] ~ Poisson(λ[z[t][c]])
        end
    end
end

m3 = infiniteTimeSeriesPMM(x,T)
iterations = 1000
@time smc_chain = begin
    burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
    n_samples = 1000
    iterations = burn + n_samples
 
    chain = sample(m3, SMC(), iterations);
end

@model function timeseries_ddp_pmm(x,T,K,W)
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_0 ~ Gamma(a_θ,1/b_θ)
    # θ_t ~ filldist(θ_0, T)
    θ_t = Vector{Float64}(undef, T)

    a=6
    b=1
    H_0 = Gamma(a,1/b)
    λ = Vector{Vector}(undef, T)
    β = Vector{Vector}(undef, T)
    η = Vector{Vector}(undef, T) # or Vector{Float64}(undef, T)
    # λ ~ filldist(H_0, K)


    v = Vector{Vector}(undef, K-1)
    # crm = DirichletProcess.(θ_t)
    crm = Vector{DirichletProcess{Float64}}()
    
    τ = Vector{Float64}(undef, T)
    ϕ = Vector{Float64}(undef, T)
    W[1] = 0.0


   
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    
    
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                x[t][c] ~ MixtureModel(Poisson,λ, eta)
            end
        end
    end
end





###### REPEATED ABOVE################
@model dp_pmm_sb2(x, T, K) = begin
    nobs_t = length.(x)
    a_θ = 1
    b_θ = 1
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=6
    b=1
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, K-1)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                eta = stickbreak(v[t])
                z[t][c] ~ Categorical(eta)
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
            else
                eta = stickbreak(v[t])
                z[t][c] ~ Categorical(eta)
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, eta)
            end
        end
        # k[t] = k_t
    end

    
    
end
burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
n_samples = 500
iterations = burn + n_samples
K = 25
T = 4
m2=dp_pmm_sb2(x, T, K)
# pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))
tchain = sample(m2, PG(100),iterations);
k = map(
    t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :k), :].value))),
    1:iterations
);
ids = findall(map(name -> occursin("λ", string(name)), names(tchain)));
p = plot(tchain[:, ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter)

###############################################
###############################################
###############################################


@model function timeseries_indep_dp_pmm1(x, T, K;a_θ = 1,b_θ = 1,a=6,b = 1) # basically the same as dp_pmm_sb2 but now we can change Hyperparameters a_θ,b_θ ,a,b
    nobs_t = length.(x)
    a_θ = a_θ
    b_θ = b_θ
    θ_t ~ filldist(Gamma(a_θ,1/b_θ), T)
    crm = DirichletProcess.(θ_t)
    # v ~ filldist.(StickBreakingProcess.(crm), K - 1)
    a=a
    b=b
    λ ~ filldist(Gamma(a,1/b), K)
    v = Vector{Vector}(undef, K-1)
    π_ = Vector{Vector}(undef, K)
    z = Vector(map(c-> tzeros(Int, nobs_t[c]), 1:T))
    for t in 1:T
        v[t] ~ filldist(StickBreakingProcess(crm[t]), K - 1)
        # k_t = Vector{Int}(undef, nobs_t[t])
        for c in 1:nobs_t[t]
            if t ==1 
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
            else
                π_[t] = stickbreak(v[t])
                z[t][c] ~ Categorical(π_[t])
                x[t][c] ~ Poisson(λ[z[t][c]])#MixtureModel(Poisson,λ, π_[t])
            end
        end
        # k[t] = k_t
    end
end

#Posterior Over Clusters Numbers
function plotNumClusters(tchain, iterations)
    h = map(
        t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value))), 1:iterations
        );
    #Posterior Over Clusters
    clus_plot = histogram(h, xlabel = "Number of clusters", legend = false);
    display(clus_plot)
    return clus_plot
end


function getlambdaIDs(tchain)
    λ_ids = findall(map(name -> occursin("λ", string(name)), names(tchain)))
    return λ_ids
end

function plotlambdaposteriorChainAndDensity(tchain,n_samples,K)
    λ_ids = getlambdaIDs(tchain)
    p_λ = plot(tchain[n_samples:end, λ_ids, :]; legend=true, labels=reshape(["λ $i" for i in 1:K],1,K), colordim=:parameter);
    display(p_λ)
    return p_λ
end

function getlambdaposteriorDF(tchain,n_samples)
    λ_ids = getlambdaIDs(tchain)
    λ_df = DataFrame(tchain[n_samples:end, λ_ids, :])[:,3:end]
    return λ_df
end
function getlambdaposteriorAvg(tchain,n_samples)
    λ_df = getlambdaposteriorDF(tchain,n_samples)
    cluster_means = [mean(col) for col in  eachcol(λ_df)]
    cluster_std = [std(col) for col in  eachcol(λ_df)]
    cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
    return cluster_mat
end
function plotlambdaposteriorSamples(tchain,n_samples,K)
    λ_df = getlambdaposteriorDF(tchain,n_samples)
    pv = violin(reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K],legend=false,yticks = 0:maximum(Matrix(λ_df))+5)
    pv = boxplot!(pv, reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K], fillalpha=0.75, linewidth=2,legend=false)
    
    display(pv)
    return pv
end

function re_func(t,var)
    return Regex(string(var)*"\\["*string(t)*"\\]\\[[0-9]+\\]")
end
function time_re_func(t,var)
    return Regex(string(var)*"\\["*string(t)*"\\]\\[[0-9]+\\]")
end
function cellID_re_func(id,var)
    return Regex(string(var)*"\\[[0-9]+\\]\\["*string(id)*"\\]")
end


function calc_r(μ,σ²)
    r = μ.^2 ./ (σ² .- μ)
    return r
end

function calc_p(μ,σ²)
    r = calc_r(μ,σ²)
    p = r ./ (r .+ μ)
    return p
end
function calc_dispersion(μ,a_)
    # r = calc_r(μ,σ²)
    σ² = μ .+ a_ .* μ .^ 2 
    return σ²
end


# Want to the distribution over cells being clustered correctly together
function plotDataGeneratingDist(μ,σ²,mixing_prob)
    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²) 
    data_generating_distribution  = MixtureModel(NegativeBinomial.(r,p),mixing_prob)
    dgd = plot(data_generating_distribution)
    display(dgd)
    return data_generating_distribution
end


function getzIDs(tchain)
    z_ids = findall(map(name -> occursin("z", string(name)), names(tchain)));
    return z_ids
end
function getzDF(tchain,n_samples)
    z_ids = getzIDs(tchain)
    z_df = DataFrame(tchain[n_samples:end, z_ids, :])[:,3:end]
    return z_df
end
function getTrueClusterMembershipDict(T,k,truth_dict)
    true_cluster_membership_dict = Dict(t => Dict(i => Set([]) for i in 1:k ) for t in 1:T)
    for t in 1:T
        cells = getCustomersAtTimepoint(truth_dict,t)
        for c in cells
            true_cluster_id = last(c)
            push!(true_cluster_membership_dict[t][true_cluster_id], first(c))
        end 
    end
    return true_cluster_membership_dict
end



function getEmpricalCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
    empirical_cluster_membership_exclusive = Dict(t => Dict() for t in 1:T) #Exactly correct, no other cells included in the cluster
    empirical_cluster_membership_inclusive = Dict(t => Dict() for t in 1:T) #Correct, but other cells are in the cluster
    numCellsPerT = [length(getCustomersAtTimepoint(truth_dict,t)) for t in 1:T]
    true_cluster_membership_dict = getTrueClusterMembershipDict(T,k,truth_dict)
    z_df = getzDF(tchain,n_samples)
    for r in eachrow(z_df)
        for t in 1:T
            numCellsT = numCellsPerT[t]
            membership_dict_t = true_cluster_membership_dict[t]
            tp_col = [m.match for m in match.(time_re_func(t,"z"),names(r)) if !isnothing(m) ]
            unique_emp_cluster_assgin = unique(collect(r[tp_col]))
            
            temp_dict = Dict(clus => Set([]) for clus in unique_emp_cluster_assgin)
            
            tp_cells = r[tp_col]
            #Assume linear indexing of cell ids!!!!
            for c in 1:numCellsT
                cell_col = [m.match for m in match.(cellID_re_func(c,"z"),names(tp_cells)) if !isnothing(m) ]
                cell_cluster_assgn = collect(tp_cells[cell_col])[1]
                push!(temp_dict[cell_cluster_assgn],c)
            end
            all_empirical_clusters = [temp_dict[key] for key in keys(temp_dict)]
            for key in keys(membership_dict_t)
                true_cluster_membership = membership_dict_t[key]
                inclusive_check = [issubset(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
                exlusive_check = [issetequal(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
                if any(inclusive_check)
                    if haskey(empirical_cluster_membership_inclusive[t],true_cluster_membership)
                        empirical_cluster_membership_inclusive[t][true_cluster_membership] += 1
                    else
                        empirical_cluster_membership_inclusive[t][true_cluster_membership] = 1
                    end 
                end
                if any(exlusive_check)
                    if haskey(empirical_cluster_membership_exclusive[t],true_cluster_membership)
                        empirical_cluster_membership_exclusive[t][true_cluster_membership] +=1
                    else
                        empirical_cluster_membership_exclusive[t][true_cluster_membership] = 1
                    end
                end
            end
            
        end
    end
    emp_clus_mem_counts_inclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
    emp_clus_mem_counts_exclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
    for t in 1:T
        for clus in 1:k
            key = true_cluster_membership_dict[t][clus]
            if haskey(empirical_cluster_membership_exclusive[t],key)
                emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
            end
            # emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
            if haskey(empirical_cluster_membership_inclusive[t],key)
                emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
            end
            # emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
        end
    end
    exclusive_rate = [[emp_clus_mem_counts_exclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
    inclusive_rate = [[emp_clus_mem_counts_inclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
    return exclusive_rate, inclusive_rate
end


function plotCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
    exclusive_rate, inclusive_rate = getEmpricalCorrectClusteringRates(tchain,n_samples,T,k,truth_dict)
    pCCR = plot(collect(1:T),inclusive_rate,labels=reshape(["Inclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:circle,legend = :outerbottom)
    pCCR = plot!(pCCR,collect(1:T),exclusive_rate,labels=reshape(["Exclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:diamond, linestyle=:dashdot, legend = :outerbottom)
    display(pCCR)
    return pCCR
end

# Set random seed for reproducibility
Random.seed!(0);

burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
n_samples = 500
PG_param = 100
# iterations = burn + n_samples
x = [[data_dict[(i,t)][1] for i in 1:6] for t in 1:4]
true_clusters = [[truth_dict[(i,t)] for i in 1:6] for t in 1:4]
coordinates = [[(i,t) for i in 1:6] for t in 1:4]
K = 12
true_k = 3
T = 4
m2=timeseries_indep_dp_pmm1(x, T, K)
# pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))

# tchain = sample(m2, PG(100),iterations);
@time PG_chain = begin
    burnin = burn  # NOTE: The burn in is also returned. Can't be discarded.
    post_samples = n_samples
    iterations = burnin + post_samples
    PG_param = PG_param
 
    tchain = sample(m2, 
            PG(PG_param),
           iterations)
end;




# h = map(
#     t -> length(unique(vec(tchain[t, MCMCChains.namesingroup(tchain, :z), :].value))),
#     1:iterations
# );
# #Posterior Over Clusters
# histogram(h, xlabel = "Number of clusters", legend = false)


# Plots Posterior Over Cluster Parameters λ
# λ_ids = findall(map(name -> occursin("λ", string(name)), names(tchain)));




# λ_df = DataFrame(tchain[n_samples:end, λ_ids, :])[:,3:end]
# cluster_means = [mean(col) for col in  eachcol(λ_df)]
# cluster_std = [std(col) for col in  eachcol(λ_df)]
# cluster_mat = hcat(collect(1:length(cluster_means)), cluster_means, cluster_std)
# violin(reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K],legend=false,yticks = 0:maximum(Matrix(λ_df))+5)
# boxplot!(reshape(["λ $i" for i in 1:K],1,K),[λ_df[:,i] for i in 1:K], fillalpha=0.75, linewidth=2,legend=false) #[λ_df[:,1] λ_df[:,2]])




# Plots Posterior Over Cluster Parameters π_ at last time point
# v_df = DataFrame(tchain[1000, MCMCChains.namesingroup(tchain, :v), :])[!,3:end]
# v_t_vec = [[v_df[!,j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(v_df)) if !isnothing(m) ]] for t in 1:T]
# # π_ids = findall(map(name -> occursin("π_", string(name)), names(tchain)));
# π_t_post = stickbreak.(v_t_vec)
# histogram(π_t_post[1]) #Not as helpful!

# k = 3
# T = 4
# Cₜ =  2*ones(Int,T) #[2,2,2,2]
# μ = [2.,4.,6.]
# σ² = 10.0*ones(k)
# r = μ.^2 ./ (σ² .- μ)
# p = r ./ (r .+ μ)
# data_generating_distribution  = MixtureModel(NegativeBinomial.(r,p),[1/3,1/3,1/3])
# plot(data_generating_distribution)
# z_ids = findall(map(name -> occursin("z", string(name)), names(tchain)));
# z_df = DataFrame(tchain[n_samples:end, z_ids, :])[:,3:end]
# true_cluster_membership_dict = Dict(t => Dict(i => Set([]) for i in 1:k ) for t in 1:T)
# for t in 1:T
#     cells = getCustomersAtTimepoint(truth_dict,t)
#     for c in cells
#         true_cluster_id = last(c)
#         push!(true_cluster_membership_dict[t][true_cluster_id], first(c))
#     end 
# end
# empirical_cluster_membership_exclusive = Dict(t => Dict() for t in 1:T) #Exactly correct, no other cells included in the cluster
# empirical_cluster_membership_inclusive = Dict(t => Dict() for t in 1:T) #Correct, but other cells are in the cluster
# numCellsPerT = [length(getCustomersAtTimepoint(truth_dict,t)) for t in 1:T]

# for r in eachrow(z_df)
#     for t in 1:T
#         numCellsT = numCellsPerT[t]
#         membership_dict_t = true_cluster_membership_dict[t]
#         tp_col = [m.match for m in match.(time_re_func(t,"z"),names(r)) if !isnothing(m) ]
#         unique_emp_cluster_assgin = unique(collect(r[tp_col]))
        
#         temp_dict = Dict(clus => Set([]) for clus in unique_emp_cluster_assgin)
        
#         tp_cells = r[tp_col]
#         #Assume linear indexing of cell ids!!!!
#         for c in 1:numCellsT
#             cell_col = [m.match for m in match.(cellID_re_func(c,"z"),names(tp_cells)) if !isnothing(m) ]
#             cell_cluster_assgn = collect(tp_cells[cell_col])[1]
#             push!(temp_dict[cell_cluster_assgn],c)
#         end
#         all_empirical_clusters = [temp_dict[key] for key in keys(temp_dict)]
#         for key in keys(membership_dict_t)
#             true_cluster_membership = membership_dict_t[key]
#             inclusive_check = [issubset(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
#             exlusive_check = [issetequal(true_cluster_membership, emp_cluster) for emp_cluster in all_empirical_clusters]
#             if any(inclusive_check)
#                 if haskey(empirical_cluster_membership_inclusive[t],true_cluster_membership)
#                     empirical_cluster_membership_inclusive[t][true_cluster_membership] += 1
#                 else
#                     empirical_cluster_membership_inclusive[t][true_cluster_membership] = 1
#                 end 
#             end
#             if any(exlusive_check)
#                 if haskey(empirical_cluster_membership_exclusive[t],true_cluster_membership)
#                     empirical_cluster_membership_exclusive[t][true_cluster_membership] +=1
#                 else
#                     empirical_cluster_membership_exclusive[t][true_cluster_membership] = 1
#                 end
#             end
#         end
        
#     end
# end
# emp_clus_mem_counts_inclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
# emp_clus_mem_counts_exclusive = Dict(t => Dict(i => 0 for i in 1:k ) for t in 1:T)
# for t in 1:T
#     for clus in 1:k
#         key = true_cluster_membership_dict[t][clus]
#         if haskey(empirical_cluster_membership_exclusive[t],key)
#             emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
#         end
#         # emp_clus_mem_counts_exclusive[t][clus] = empirical_cluster_membership_exclusive[t][key]
#         if haskey(empirical_cluster_membership_inclusive[t],key)
#             emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
#         end
#         # emp_clus_mem_counts_inclusive[t][clus] = empirical_cluster_membership_inclusive[t][key]
#     end
# end
# exclusive_rate = [[emp_clus_mem_counts_exclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]
# inclusive_rate = [[emp_clus_mem_counts_inclusive[t][clus] for t in 1:T] for clus in 1:k ] ./ size(z_df)[1]


# plot(collect(1:T),inclusive_rate,labels=reshape(["Inclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:circle)
# plot!(collect(1:T),exclusive_rate,labels=reshape(["Exclusive Rate for Cluster $i" for i in 1:k],1,k),marker=:diamond, linestyle=:dashdot)


plotNumClusters(tchain, iterations)
plotlambdaposteriorChainAndDensity(tchain,n_samples,K)
λ_ids = getlambdaIDs(tchain)
meanplot(tchain[n_samples:end, λ_ids, :])
histogram(tchain[n_samples:end, λ_ids, :])
plotlambdaposteriorSamples(tchain,n_samples,K)
v_df = DataFrame(tchain[1000, MCMCChains.namesingroup(tchain, :v), :])[!,3:end]
v_t_vec = [[v_df[!,j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(v_df)) if !isnothing(m) ]] for t in 1:T]
π_t_post = stickbreak.(v_t_vec)
histogram(π_t_post[1]) #Not as helpful!
plotDataGeneratingDist(μ,σ²,[1/3,1/3,1/3])
plotCorrectClusteringRates(tchain,n_samples,T,true_k,truth_dict)


####################################################################################################################

# NEW DATA GENERATED: This was a mixture with very different components. 

####################################################################################################################

# Want to the distribution over cells being clustered correctly together
k = 3
T = 4
Cₜ =  2*ones(Int,T) #[2,2,2,2]
μ = [1.,20.,50.]
a_σ = .5
σ² = μ .+ a_σ .* μ .^ 2 
r = μ.^2 ./ (σ² .- μ)
p = r ./ (r .+ μ)
data_generating_distribution2  = MixtureModel(NegativeBinomial.(r,p),[1/3,1/3,1/3])
plot(data_generating_distribution2)

_truth_dict, _data_dict = generate1Phen_SimpleCase1(μ = μ, σ² = σ²  )

# Set random seed for reproducibility
Random.seed!(0);

burn = 500  # NOTE: The burn in is also returned. Can't be discarded.
n_samples = 1000
PG_param = 100
# iterations = burn + n_samples
x = [[_data_dict[(i,t)][1] for i in 1:6] for t in 1:4]
true_clusters = [[_truth_dict[(i,t)] for i in 1:6] for t in 1:4]
coordinates = [[(i,t) for i in 1:6] for t in 1:4]
K = 12
true_k = 3
T = 4
m3=timeseries_indep_dp_pmm1(x, T, K)
# pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))

# tchain = sample(m2, PG(100),iterations);
@time PG_chain = begin
    burnin = burn  # NOTE: The burn in is also returned. Can't be discarded.
    post_samples = n_samples
    iterations = burnin + post_samples
    PG_param = PG_param
 
    tchain = sample(m3, 
            PG(PG_param),
           iterations)
end;


plotNumClusters(tchain, iterations)
plotlambdaposteriorChainAndDensity(tchain,n_samples,K)
λ_ids = getlambdaIDs(tchain)
meanplot(tchain[n_samples:end, λ_ids, :])
histogram(tchain[n_samples:end, λ_ids, :])
plotlambdaposteriorSamples(tchain,n_samples,K)
v_df = DataFrame(tchain[1000, MCMCChains.namesingroup(tchain, :v), :])[!,3:end]
v_t_vec = [[v_df[!,j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(v_df)) if !isnothing(m) ]] for t in 1:T]
π_t_post = stickbreak.(v_t_vec)
histogram(π_t_post[1]) #Not as helpful!
plotDataGeneratingDist(μ,σ²,[1/3,1/3,1/3])
plotCorrectClusteringRates(tchain,n_samples,T,true_k,truth_dict)

k = 3
T = 4
a_σ = 0.5
Cₜ =  2*ones(Int,T) #[2,2,2,2]
μ = [2.,4.,6.]
σ² = calc_dispersion(μ,a_σ)

@model function gen1Phen_SimpleCase1_NB_TimeInvar(x,T,k,C_t,μ,σ²,α,mixing_prob)
    if x === missing
        x = Vector{Vector}(undef,T)
        for t in 1:T
            cells = C_t[t]
            x[t] = tzeros(Int, cells)
        end
    end
    if μ >= σ²
        throw(DomainError(μ, "μ must be less than σ²"))
    end
    if k !== maximum([length(μ), length(σ²)])
        error("must ha μ and σ² for each of the k components")
    end
    r = calc_r(μ,σ²)
    p = calc_p(μ,σ²)
    if mixing_prob === missing
        _prob ~ Dirichlet(α .* ones(k) ./ k)
    else
        _prob ~ arraydist(Dirac.(mixing_prob))
    end
    # mixing_prob ~ Dirichlet(α .* ones(k) ./ k)
    H = NegativeBinomial.(r,p)
    z = Vector{Vector}(undef,T)
    for t in 1:T
        cells = C_t[t]
        z[t] = tzeros(Int, cells)
        for c in 1:cells
            z[t][c] ~ Categorical(_prob)
            x[t][c] ~ H[z[t][c]]
        end 
        
    end
    return x,z
    # data_dist = MixtureModel(NegativeBinomial.(r,p),mixing_prob)

    
end
C_t = 6 .* ones(Int,T)
α = 2
model = gen1Phen_SimpleCase1_NB_TimeInvar(missing,T,k,C_t,μ,σ²,α,ones(k) ./ k)
c = sample(model, Prior(), 1)
data_df = DataFrame(c)[!,3:end]
col_names_z = [[m.match for m in match.(time_re_func(t,"z"),names(data_df)) if !isnothing(m) ] for t in 1:T]
z_true = [[data_df[:,c][1] for c in col_names_z[t]] for t in 1:T]
col_names_x = [[m.match for m in match.(time_re_func(t,"x"),names(data_df)) if !isnothing(m) ] for t in 1:T]
x_data =  [[data_df[:,c][1:end]   for c in col_names_x[t]] for t in 1:T]
mixing_prob_ids = findall(map(name -> occursin("_prob", string(name)), names(data_df)))
mixing_prob = vec(collect(data_df[1,mixing_prob_ids]))
plotDataGeneratingDist(μ,σ²,mixing_prob)
function genData_vec2dict(time_vec)
    T = length(time_vec)
    C_t = [length(c) for c in time_vec]
    final_dict = Dict()
    for t in 1:T
        for c in 1:C_t[t]
            final_dict[(c,t)] = time_vec[t][c]
        end
    end
    # final_dict = Dict((c,t) => time_vec[t][c] for c in 1:C_t[t] for t in 1:T)
    return final_dict
end
data_dict2 = genData_vec2dict(x_data)
truth_dict2 = genData_vec2dict(z_true)
x = [[data_dict2[(i,t)][1] for i in 1:C_t[t]] for t in 1:T]
true_clusters = [[truth_dict2[(i,t)] for i in 1:C_t[t]] for t in 1:T]
coordinates = [[(i,t) for i in 1:C_t[t]] for t in 1:T]
K = 12
true_k = 3
T = 4
m4=timeseries_indep_dp_pmm1(x, T, K)
# pmm_sampler = Gibbs(PG(100, :k_t), HMC(0.05, 10))

# tchain = sample(m2, PG(100),iterations);
@time PG_chain = begin
    burnin = burn  # NOTE: The burn in is also returned. Can't be discarded.
    post_samples = n_samples
    iterations = burnin + post_samples
    PG_param = PG_param
 
    tchain = sample(m4, 
            PG(PG_param),
           iterations)
end;

plotNumClusters(tchain, iterations)
plotlambdaposteriorChainAndDensity(tchain,n_samples,K)
λ_ids = getlambdaIDs(tchain)
meanplot(tchain[n_samples:end, λ_ids, :])
histogram(tchain[n_samples:end, λ_ids, :])
plotlambdaposteriorSamples(tchain,n_samples,K)
v_df = DataFrame(tchain[1000, MCMCChains.namesingroup(tchain, :v), :])[!,3:end]
v_t_vec = [[v_df[!,j][1] for j in [m.match for m in match.(time_re_func(t,"v"),names(v_df)) if !isnothing(m) ]] for t in 1:T]
π_t_post = stickbreak.(v_t_vec)
histogram(π_t_post[1]) #Not as helpful!
plotDataGeneratingDist(μ,σ²,mixing_prob)
plotCorrectClusteringRates(tchain,n_samples,T,true_k,truth_dict)