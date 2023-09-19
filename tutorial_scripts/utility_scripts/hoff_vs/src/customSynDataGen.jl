

function generate1Phen_SimpleCase1(;k = 3, T = 4, CperK = 2*ones(Int,T),μ = [2.,4.,6.], σ² = 10.0*ones(k)  )
    truth_dict = Dict{Tuple,Float64}()
    data_dict = Dict{Tuple,VecOrMat}()
    r = μ.^2 ./ (σ² .- μ)
    p = r ./ (r .+ μ)
    for t in 1:T
        c = CperK[t]
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
function generateGPhen_SimpleCase1(k, T, CperK ,μ, σ²)
    truth_dict = Dict{Tuple,Float64}()
    data_dict = Dict{Tuple,VecOrMat}()
    r = vec(map(clus -> calc_r(μ[clus],σ²[clus]), 1:k))
    p = vec(map(clus -> calc_p(μ[clus],σ²[clus]), 1:k))#calc_p(μ,σ²)
    for t in 1:T
        c = CperK[t]
        dist_ = vec(map(clus -> arraydist(NegativeBinomial.(r[clus],p[clus])), 1:k))
        phen = permutedims.(rand.(dist_,c))
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
function generateGPhenNormal_SimpleCase1(k, T, CperK ,μ, σ²)
    truth_dict = Dict{Tuple,Float64}()
    data_dict = Dict{Tuple,VecOrMat}()
    # r = vec(map(clus -> calc_r(μ[clus],σ²[clus]), 1:k))
    # p = vec(map(clus -> calc_p(μ[clus],σ²[clus]), 1:k))#calc_p(μ,σ²)
    function censor_and_round(v)
        v[v .< 0.0] .= 0.0
        return round.(v) 
    end
    for t in 1:T
        c = CperK[t]
        dist_ = vec(map(clus -> arraydist(Normal.(μ[clus],σ²[clus])), 1:k))
        phen = permutedims.(rand.(dist_,c))
        phen = map(phn -> censor_and_round(phen[phn]), 1:k)
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
function generateGMM(;k = 3, T = 4, CperK = 2*ones(Int,T),μ = [2.,4.,6.], σ² = 10.0*ones(k)  )
    truth_dict = Dict{Tuple,Float64}()
    data_dict = Dict{Tuple,VecOrMat}()
    for t in 1:T
        c = CperK[t]
        phen = rand.(Normal.(μ,σ²),c)
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



function fake_mvGausssian_data_for_testing(G,C_t,Ktrue;μ =nothing,τ = nothing, mix_prob =nothing,same_prob_t = true,dynamic = false)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* 10.0 .^(power[k]) .* ones(G)#rand(G)
        end
    end
    if isnothing(τ)
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            x[t][c] = rand(MultivariateNormal(μ[z[t][c]],diagm(sqrt.((τ[z[t][c]]) .^(-1))) ))
        end
    end
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob,μ,τ 
end

function fake_mvGausssian_corrdata_for_testing(C_t,KTrue,μk ,Σk; mix_prob =nothing,same_prob_t = true,dynamic = false)
    T = length(C_t)
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(KTrue)./KTrue 
        else
            mix_prob = rand(Dirichlet(KTrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:KTrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            x[t][c] = rand(MultivariateNormal(μk[z[t][c]],Σk[z[t][c]] ))
        end
    end
    
    return x,z,mix_prob
end

function fake_mvCountsViaCopulas_corrdata_for_testing(C_t,KTrue,μk ,Σk; mix_prob =nothing,same_prob_t = true,dynamic = false,epsilon=0.001)
    T = length(C_t)
    G = length(μk[1])
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(KTrue)./KTrue 
        else
            mix_prob = rand(Dirichlet(KTrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:KTrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Int64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Int64}}(undef,C_t[t])
        for c in 1:C_t[t]
            corr_normal_data = rand(MultivariateNormal(zeros(G),Σk[z[t][c]] ))
            corr_uniform_data = cdf.(Normal(0,1),corr_normal_data)
            corr_gamma_data = quantile.(Gamma.(μk[z[t][c]] .+ epsilon,1),corr_uniform_data)
            x[t][c] = rand.(Poisson.(corr_gamma_data))
        end
    end
    
    return x,z,mix_prob
end

function fake_mvGausssian_data_largeG_indepG_for_testing(G,C_t,Ktrue;μ =nothing,τ = nothing, mix_prob =nothing,same_prob_t = true,dynamic = false)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* 100.0 .^(power[k]) .* rand(G)
        end
    end
    if isnothing(τ)
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            x[t][c] = rand(arraydist(Normal.(μ[z[t][c]],sqrt.((τ[z[t][c]]) .^(-1)) )))
        end
    end
        # τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob,μ,τ 
end



function fake_mvGausssian_data_indepG_for_testing(G,C_t,Ktrue;μ =nothing,τ = nothing, mix_prob =nothing,same_prob_t = true,dynamic = false,μ_magnitude = 10.0,largeG = true)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* μ_magnitude .^(power[k]) .* rand(G)
        end
    elseif typeof(μ) <: VecOrMat && eltype(μ) <: Number
        μ = [μ[k] .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(τ)
        τ = [1.0 .* ones(G) for k in 1:Ktrue]
    elseif typeof(τ) <: VecOrMat && eltype(τ) <: Number
        τ = [τ[k] .* ones(G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[t]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            if largeG 
                z[t][c]
                μ[z[t][c]]
                τ[z[t][c]]
                x[t][c] = rand(arraydist(Normal.(μ[z[t][c]],sqrt.((τ[z[t][c]]) .^(-1)) )))
            else
                x[t][c] = rand(MultivariateNormal(μ[z[t][c]],diagm(sqrt.((τ[z[t][c]]) .^(-1))) ))
            end
            
        end
    end
        # τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob,μ,τ 
end


function fake_Poisson_data_for_testing(G,causalG,housekeepingG,C_t,Ktrue;μ_causal =nothing,τ = 10, mix_prob =nothing,same_prob_t = true,dynamic = false,a_gamma_causal =2 ,b_gamma_causal = 3,a_gamma_noise =1 ,b_gamma_noise = 2,causal_influence = 0.2,noise_influence = 0.2)
    function censor_and_round(v)
        v[v .< 0.0] .= 0.0
        return round.(v) 
    end
    noiseG = G - causalG - housekeepingG
    # noise_influence = 1-causal_influence
    num_noise_genes_to_choose = Int(ceil(noise_influence * noiseG))
    num_causal_genes_to_choose = Int(ceil(causal_influence * causalG))
    causalG_combination_indices = Vector{Vector{Int}}(undef,housekeepingG)
    noiseG_combination_indices = Vector{Vector{Int}}(undef,housekeepingG)
    for g in 1:housekeepingG
        causalG_combination_indices[g] = sample(1:causalG, num_causal_genes_to_choose,replace=false)
        noiseG_combination_indices[g] = sample(1:noiseG, num_noise_genes_to_choose,replace=false)
    end
    T = length(C_t)
    if isnothing(μ_causal)
        μ_causal = Vector{Vector{Float64}}(undef,Ktrue)
        μ_noise = rand(filldist(Gamma(a_gamma_noise,b_gamma_noise),noiseG))
        for k in 1:Ktrue
            μ_causal[k] = rand(filldist(Gamma(a_gamma_causal,b_gamma_causal),causalG))
        end
    # elseif typeof(μ) <: VecOrMat && eltype(μ) <: Number
    #     μ = [μ[k] .* ones(causalG) for k in 1:Ktrue]
    end
    # if isnothing(τ)
    #     τ = [1.0 .* ones(G) for k in 1:Ktrue]
    # # elseif typeof(τ) <: VecOrMat && eltype(τ) <: Number
    # #     τ = [τ[k] .* ones(G) for k in 1:Ktrue]
    # end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = ones(Ktrue)./Ktrue 
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end

        if dynamic
            α_0 = 1.0
            mix_prob = [rand(Dirichlet( α_0 .* mix_prob)) for t in 1:T]
        else
            mix_prob = [mix_prob for t in 1:T]
        end
    end
    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            causal_expression = rand(arraydist(Normal.(μ_causal[z[t][c]],τ^(-1))))
            noise_expression = rand(arraydist(Normal.(μ_noise,τ^(-1))))
            causal_housekeeping_contribution = [causal_expression[el] for el in causalG_combination_indices]
            noise_housekeeping_contribution = [noise_expression[el] for el in noiseG_combination_indices]
            μ_housekeeping = 1 ./(num_noise_genes_to_choose + num_causal_genes_to_choose) .* (sum.(causal_housekeeping_contribution) .+ sum.(noise_housekeeping_contribution))
            houskeeping_expression = rand(arraydist(Normal.(μ_housekeeping,τ^(-1))))
            x[t][c] = censor_and_round(reduce(vcat,[causal_expression,houskeeping_expression,noise_expression]))
            # if largeG 
            #     x[t][c] = rand(arraydist(Normal.(μ[z[t][c]],(τ[z[t][c]]) .^(-1) )))
            # else
            #     x[t][c] = rand(MultivariateNormal(μ[z[t][c]],diagm((τ[z[t][c]]) .^(-1)) ))
            # end
            
        end
    end
        # τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z,mix_prob, μ_causal, μ_noise
end

function fake_generateGPhenNormal_RoundedCensored_indepG_for_testing(C_t,Ktrue,μ,τ, mix_prob;largeG = true)
    function censor_and_round(v)
        v[v .< 0.0] .= 0.0
        return round.(Int,v) 
    end
    T = length(C_t)

    
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[t]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Int}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Int}}(undef,C_t[t])
        for c in 1:C_t[t]
            if largeG 
                x[t][c] =censor_and_round( rand(arraydist(Normal.(μ[z[t][c]],sqrt.((τ[z[t][c]]) .^(-1)) ))))
            else
                x[t][c] = censor_and_round(rand(MultivariateNormal(μ[z[t][c]],diagm(sqrt.((τ[z[t][c]]) .^(-1))) )))
            end
            
        end
    end
        # τ = [1.0 .* ones(G) for k in 1:Ktrue]
    return x,z
end

function generate_number_of_DEG(G,percent_important)
    num_important_features = round(Int,G*percent_important)
    if iszero(num_important_features)
        num_important_features=1
    end
    return Int(num_important_features)
end
function generate_DEG_locations(G,K,num_important_features;random_locations=false)
    deg_locations = Vector{Vector{Bool}}(undef,K)
    if !random_locations
        for k in 1:K
            deg_locations[k] = falses(G)
            true_indices = mod1.((k-1)*num_important_features+1:k*num_important_features, G)
            deg_locations[k][true_indices] .= true
        end
    else
        s = Set()
        println("here")
        while length(s) < K
            i = randperm(G)[1:num_important_features]
            push!(s,i)
        end
        # used_locations = shuffle(collect(combinations(BigInt.(1:G),BigInt(num_important_features))))[1:K]
        used_locations = collect(collect.(s))
        # println(used_locations)
        for k in 1:K
            deg_locations[k] = falses(G)
            true_indices = used_locations[k]
            deg_locations[k][true_indices] .= true
        end
    end
    return deg_locations
end
function generate_DEG_values!(μk,deg_locations,high_de;low_de=0,random_deg_values=false)
    K = length(deg_locations)
    G = length(deg_locations[1])
    num_important_features = sum(deg_locations[1])
    for k in 1:K
        cluster_deg_locations = deg_locations[k]
        if random_deg_values
            μk[k][cluster_deg_locations] .= rand(Uniform(low_de,high_de),num_important_features)
        else
            μk[k][cluster_deg_locations] .= high_de
        end
    end
    return μk
end
function generate_nonDEG_values!(μk,deg_locations,high_nonde;low_nonde=0,random_nondeg_values=false)
    K = length(deg_locations)
    G = length(deg_locations[1])
    num_non_important_features = G - sum(deg_locations[1])
    for k in 1:K
        cluster_nondeg_locations = broadcast(!,deg_locations[k])
        if random_nondeg_values
            μk[k][cluster_nondeg_locations] .= rand(Uniform(low_nonde,high_nonde),num_non_important_features)
        else
            μk[k][cluster_nondeg_locations] .= high_nonde
        end
    end
    return μk
end
function generate_DEG_means(G,indices,low_de,high_de)
    means_vec = zeros(G)
    num_deg = length(indices)
    means_vec[indices] =  rand(Uniform(low_de,high_de),num_deg)
    return means_vec
end
function generate_rand_covariance_matrix(G, KTrue;corr_feat = true,one_corr_mat = true)
    if corr_feat && one_corr_mat
        common_cov_mat = rand(LKJ(G,1));
        Σk= [common_cov_mat for k in 1:KTrue];
    elseif corr_feat && !one_corr_mat
        Σk= [rand(LKJ(G,1)) for k in 1:KTrue];
    else
        common_cov_mat = Matrix(1.0I, G, G);Matrix(1.0I, G, G);
        Σk= [common_cov_mat for k in 1:KTrue];
    end  
    return Σk  
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





# module customSyntheticDataGeneration

#     # curr_dir = ENV["PWD"]
#     # src_dir = "/src/"

#     # using Random
#     # using Distributions
#     # using Turing
#     # using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     # using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     # using Test
#     # import Debugger
#     # using CSV,DataFrames


#     # include(curr_dir*src_dir*"MathUtils.jl")
#     # using .MathUtils
#     #### Data Generation Functions ######
#     export generate1Phen_SimpleCase1,generateGMM
#     function generate1Phen_SimpleCase1(;k = 3, T = 4, CperK = 2*ones(Int,T),μ = [2.,4.,6.], σ² = 10.0*ones(k)  )
#         truth_dict = Dict{Tuple,Float64}()
#         data_dict = Dict{Tuple,VecOrMat}()
#         r = μ.^2 ./ (σ² .- μ)
#         p = r ./ (r .+ μ)
#         for t in 1:T
#             c = CperK[t]
#             phen = rand.(NegativeBinomial.(r,p),c)
#             phen = vcat(phen...)
#             cluster_assgn = [i for i in 1:k for j in 1:c]
#             concat_data = hcat(cluster_assgn,phen)
#             concat_data = concat_data[shuffle(1:end),:]
#             cluster_assgn, phen = concat_data[:,1], concat_data[:,2:end]
#             cell_ids = [ i for i in 1:c*k]
#             for i in 1:length(cell_ids)
#                 id = cell_ids[i]
#                 truth_dict[(id,t)] = cluster_assgn[id]
#                 data_dict[(id,t)] = phen[id,:]
#             end
#         end
#         return truth_dict, data_dict
#     end
#     function generateGPhen_SimpleCase1(k, T, CperK ,μ, σ²)
#         truth_dict = Dict{Tuple,Float64}()
#         data_dict = Dict{Tuple,VecOrMat}()
#         r = vec(map(clus -> calc_r(μ[clus],σ²[clus]), 1:k))
#         p = vec(map(clus -> calc_p(μ[clus],σ²[clus]), 1:k))#calc_p(μ,σ²)
#         for t in 1:T
#             c = CperK[t]
#             dist_ = vec(map(clus -> arraydist(NegativeBinomial.(r[clus],p[clus])), 1:k))
#             phen = permutedims.(rand.(dist_,c))
#             phen = vcat(phen...)
#             cluster_assgn = [i for i in 1:k for j in 1:c]
#             concat_data = hcat(cluster_assgn,phen)
#             concat_data = concat_data[shuffle(1:end),:]
#             cluster_assgn, phen = concat_data[:,1], concat_data[:,2:end]
#             cell_ids = [ i for i in 1:c*k]
#             for i in 1:length(cell_ids)
#                 id = cell_ids[i]
#                 truth_dict[(id,t)] = cluster_assgn[id]
#                 data_dict[(id,t)] = phen[id,:]
#             end
#         end
#         return truth_dict, data_dict
#     end
#     function generateGPhenNormal_SimpleCase1(k, T, CperK ,μ, σ²)
#         truth_dict = Dict{Tuple,Float64}()
#         data_dict = Dict{Tuple,VecOrMat}()
#         # r = vec(map(clus -> calc_r(μ[clus],σ²[clus]), 1:k))
#         # p = vec(map(clus -> calc_p(μ[clus],σ²[clus]), 1:k))#calc_p(μ,σ²)
#         function censor_and_round(v)
#             v[v .< 0.0] .= 0.0
#             return round.(v) 
#         end
#         for t in 1:T
#             c = CperK[t]
#             dist_ = vec(map(clus -> arraydist(Normal.(μ[clus],σ²[clus])), 1:k))
#             phen = permutedims.(rand.(dist_,c))
#             phen = map(phn -> censor_and_round(phen[phn]), 1:k)
#             phen = vcat(phen...)
#             cluster_assgn = [i for i in 1:k for j in 1:c]
#             concat_data = hcat(cluster_assgn,phen)
#             concat_data = concat_data[shuffle(1:end),:]
#             cluster_assgn, phen = concat_data[:,1], concat_data[:,2:end]
#             cell_ids = [ i for i in 1:c*k]
#             for i in 1:length(cell_ids)
#                 id = cell_ids[i]
#                 truth_dict[(id,t)] = cluster_assgn[id]
#                 data_dict[(id,t)] = phen[id,:]
#             end
#         end
#         return truth_dict, data_dict
#     end
#     function generateGMM(;k = 3, T = 4, CperK = 2*ones(Int,T),μ = [2.,4.,6.], σ² = 10.0*ones(k)  )
#         truth_dict = Dict{Tuple,Float64}()
#         data_dict = Dict{Tuple,VecOrMat}()
#         for t in 1:T
#             c = CperK[t]
#             phen = rand.(Normal.(μ,σ²),c)
#             phen = vcat(phen...)
#             cluster_assgn = [i for i in 1:k for j in 1:c]
#             concat_data = hcat(cluster_assgn,phen)
#             concat_data = concat_data[shuffle(1:end),:]
#             cluster_assgn, phen = concat_data[:,1], concat_data[:,2:end]
#             cell_ids = [ i for i in 1:c*k]
#             for i in 1:length(cell_ids)
#                 id = cell_ids[i]
#                 truth_dict[(id,t)] = cluster_assgn[id]
#                 data_dict[(id,t)] = phen[id,:]
#             end
#         end
#         return truth_dict, data_dict
#     end



#     #### Validation Functions ######
#     export check_if_truth_dict_valid
#     function check_if_truth_dict_valid(truth_dict,k,T,cellsInCluster)
#         unique_clusters = unique([truth_dict[key] for key in keys(truth_dict)])
#         unique_clusters = sort(unique_clusters)
#         true_clusters =  collect(1:k .* ones(eltype(unique_clusters)))
#         unique_clusters == true_clusters # Checks if number of clusters is correct
#         if unique_clusters == true_clusters
#             counts_chckr = Vector{Bool}()
#             for cluster in 1:k
#                 true_number_of_cells_in_cluster = sum([cellsInCluster[t] for t in 1:T]) 
#                 gen_number_of_cells_in_clusters = sum([truth_dict[key]==k for key in keys(truth_dict)])
#                 push!(counts_chckr, true_number_of_cells_in_cluster == gen_number_of_cells_in_clusters)
#             end
#             if  all(counts_chckr)
#                 return true
#             else
#                 return "counts within clusters is NOT correct"
#             end
#         else
#             return "number of clusters is NOT correct"
#         end
#     end


# end