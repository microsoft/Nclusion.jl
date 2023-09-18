
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
# a_θ,b_θ  = 0.1,0.1
# theta_prior = Gamma(a_θ,1/b_θ)
# θ_t = generate_θ_t(theta_prior,T)
# gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,1)
# gibbs_n_ht = Vector{Vector{Dict}}(undef,1)
# gibbs_likelihood_vec_t = Vector{Vector}(undef,1)
# gibbs_V_ht= Vector{Vector}(undef,1)
# gibbs_π_ht = Vector{Vector}(undef,1)
# gibbs_cluster_assignment= Vector{Vector}(undef,1)
# gibbs_norm_cluster_prob = Vector{Vector}(undef,1)
# gibbs_θ_t = Vector{Vector}(undef,1)
# π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
# ss = _2ndBlockedGibbs_AllTimeStep(1,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t) 
# ss2 = _2ndBlockedGibbs_AllTimeStep(1,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht, gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t, b_θ =1/25)


# test_itr = 10000
# gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,test_itr)
# gibbs_n_ht = Vector{Vector{Dict}}(undef,test_itr) 
# gibbs_likelihood_vec_t = Vector{Vector}(undef,test_itr)
# gibbs_V_ht= Vector{Vector}(undef,test_itr)
# gibbs_π_ht = Vector{Vector}(undef,test_itr)
# gibbs_cluster_assignment= Vector{Vector}(undef,test_itr)
# gibbs_norm_cluster_prob = Vector{Vector}(undef,test_itr)
# gibbs_θ_t = Vector{Vector}(undef,test_itr)
# π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
# for k in 1:test_itr
#     _2ndBlockedGibbs_AllTimeStep(k,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t)
# end

# K_max = 25
# numBurnin = 10000
# maxNumItr = 10000
 
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

# _BlockedGibbs_AllTimeStepInit(data_dict;K_max = 5,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = 1,b=1,a_θ =1, b_θ =1/5)

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




# K_max = 5
# a = 6
# b= 1
# a_θ = 1
# b_θ =1/5
# numBurnin = 10000
# maxNumItr = 10000
# maxNumItr = maxNumItr + numBurnin

# gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,maxNumItr+1)
# gibbs_n_ht = Vector{Vector{Dict}}(undef,maxNumItr+1) 
# gibbs_likelihood_vec_t = Vector{Vector}(undef,maxNumItr+1)
# gibbs_V_ht= Vector{Vector}(undef,maxNumItr+1)
# gibbs_π_ht = Vector{Vector}(undef,maxNumItr+1)
# gibbs_cluster_assignment= Vector{Vector}(undef,maxNumItr+1)
# gibbs_norm_cluster_prob = Vector{Vector}(undef,maxNumItr+1)
# gibbs_θ_t = Vector{Vector}(undef,maxNumItr+1)



# rand_init_θ_t, rand_init_likelihood_vec, rand_init_π_ht, rand_init_V_ht,likelihood_prior,theta_prior  = _BlockedGibbs_AllTimeStepInit(data_dict;K_max = K_max,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = a,b=b,a_θ =a_θ, b_θ =b_θ)

# gibbs_likelihood_vec_t[1] = rand_init_likelihood_vec
# gibbs_V_ht[1] =  rand_init_V_ht
# gibbs_π_ht[1] = rand_init_π_ht
# gibbs_θ_t[1] = 0.001ones(T)#rand_init_θ_t

# # _3rdBlockedGibbs_AllTimeStep(2,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,K_max = K_max)

# for step in 2:maxNumItr+1
#     itr = step
#     _3rdBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,a = a,b=b,K_max = K_max)
#     # _4thBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,a = a,b=b,K_max = K_max)
# end



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
#vv = _1stBlockedGibbs_AllTimeStep(data_dict) 

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
# _1stBlockedGibbs_1stTimeStep(data_dict)

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
# TODO: ADD THESE TESTS
# sim_time = 1000
# k_max_simtime = Vector{Int}(undef,sim_time)
# for s in 1:sim_time
#     k_max_simtime[s],_,_,_ =_firstAttempt(data_dict)
# end
# histogram(k_max_simtime)

##### NOT FINISHED ######################
function seqCRP(customers,tp_cluster_assgn,n_ht, curr_k_max;ft=false)
    itr = 1

end
function _secondAttempt()
   #TODO 
end
function _1stBlockedGibbs()
    #TODO
end


# module customInfernce
#     include("synDataPreprocess.jl")
#     using  .syntheticDataPreprocessing

#     using Random
#     using Distributions
#     using Turing
#     using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
#     using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
#     using Test
#     import Debugger
#     using CSV,DataFrames

#     export addNewTimepoint!
#     function addNewTimepoint!(prev_tp,
#                     tp_cluster_assgn,
#                                 n_ht)
#         curr_tp = prev_tp + 1
#         prev_keys1 = collect(keys(tp_cluster_assgn[prev_tp]))
#         prev_keys2 = collect(keys(n_ht[prev_tp]))
#         for  key in prev_keys1
#             tp_cluster_assgn[curr_tp][key] = []
#         end
#         for key in prev_keys2
#             n_ht[curr_tp][key] =  0
#         end
#     end

#     export addClusterFromGammaBaseDist!
#     function addClusterFromGammaBaseDist!(customer, 
#                                                 tp, 
#                                                 curr_k_max, 
#                                                 tp_cluster_assgn, 
#                                                 n_ht,
#                                                 likelihood_vec,
#                                                 lambda_prior_dist;
#                                                 ft=false, H=Poisson )

#         curr_k_max += 1 # Increase total Number
#         tp_cluster_assgn[tp][curr_k_max] = [customer] # Add to Cluster
#         n_ht[tp][curr_k_max] = 1   # Increase Cluster Count
#         if !ft
#             for i in 1:tp-1
#             n_ht[i][curr_k_max] = 0 # Update Kmax for past time points
#             tp_cluster_assgn[i][curr_k_max] = []  # Update Kmax for past time points
#             end
#         end
#         suff_stat = last(customer)[1]
#         n_ = 1
#         shape_update = lambda_prior_dist.α + suff_stat
#         scale_update = 1/(lambda_prior_dist.θ) +  n_
#         λ = rand(Gamma(shape_update,scale_update))
#         push!(likelihood_vec,H(λ))
#         return curr_k_max
#     end

#     export addToExistingPoisGammaCluster!
#     function addToExistingPoisGammaCluster!(customer,
#                                             cluster_prob,
#                                             tp,
#                                             tp_cluster_assgn,
#                                             n_ht,
#                                             likelihood_vec,
#                                             lambda_prior_dist;
#                                             ft=false, H=Poisson, is_cluster_prob_norm = false 
#                                             )

#         if is_cluster_prob_norm
#             norm_cluster_prob = cluster_prob
#         else
#             norm_cluster_prob = softmax(cluster_prob)#cluster_prob ./sum(cluster_prob)
#         end
#         cluster_val = rand(Categorical(norm_cluster_prob)) # Determine Cluster
#         push!(tp_cluster_assgn[tp][cluster_val],customer)# Add to Cluster
#         n_ht[tp][cluster_val] += 1 # Increase Cluster Count
#         if  !ft
#             prev_t_suff_stat_bool = isempty(last.(tp_cluster_assgn[tp-1][cluster_val]))
#             if prev_t_suff_stat_bool
#                 prev_t_suff_stat = 0 
#             else
#                 prev_t_suff_stat = sum(last.(tp_cluster_assgn[tp-1][cluster_val]))[1]
#             end
#             prev_t_n_ = n_ht[tp-1][cluster_val]
#         else
#             prev_t_suff_stat = 0 
#             prev_t_n_ = 0 
#         end
#         curr_t_suff_stat = sum(last.(tp_cluster_assgn[tp][cluster_val]))[1] 
#         curr_t_n_ = n_ht[tp][cluster_val]
#         suff_stat = curr_t_suff_stat + prev_t_suff_stat
#         n_ = curr_t_n_ + prev_t_n_
#         λ = rand(Gamma(lambda_prior_dist.α + suff_stat, n_ + 1/(lambda_prior_dist.θ)))  ## Sample new parameters
#         likelihood_vec[cluster_val] = H(λ)
#     end

#     export generate_V_h_init, generate_θ_t, generate_π_h_init
#     function generate_θ_t(theta_prior,T)
#         # theta_prior = Gamma(a_θ,1/b_θ)
#         θ_t =  [rand(theta_prior) for t in 1:T]
#         return θ_t
#     end
#     function generate_V_h_init(T,K_max,θ_t)
#         V_h_init = Vector{Vector}(undef,T)
#         for t in 1:T
#             V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
#         end
#         return V_h_init
#     end
#     function generate_π_h_init(T,K_max,θ_t)
#         π_h_init = Vector{Vector}(undef,T)
#         V_h_init = generate_V_h_init(T,K_max,θ_t)
#         for t in 1:T
#             π_h_init[t] = stickbreak(V_h_init[t])
#         end
#         return π_h_init, V_h_init
#     end

#     export _2ndBlockedGibbs_AllTimeStep
#     function _2ndBlockedGibbs_AllTimeStep(itr,data_dict,V_ht,π_ht,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t;K_max = 5,H = Poisson,a = 1,b=1,a_θ =1, b_θ =1/5)
#         ###### IN LOOP FORM ###########
#         if maximum(unique([last(key) for key in keys(data_dict)])) == length(unique([last(key) for key in keys(data_dict)]))
#             T = length(unique([last(key) for key in keys(data_dict)]))
#         else
#             T = maximum(unique([last(key) for key in keys(data_dict)])) 
#         end
    
#         tp_cluster_assgn_t = [Dict(h => [] for h in 1:K_max) for i in 1:T]
#         n_ht_t = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
#         # rand_gibbs_init_cluster = [ rand(1:K_max,c) for c in 3 .* Cₜ]
        
    
#         theta_prior = Gamma(a_θ,1/b_θ)
#         θ_t =  [rand(theta_prior) for t in 1:T]
#         likelihood_prior = Gamma(a,1/b)
#         likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
#         # likelihood_vec_init = deepcopy(likelihood_vec)
        
    
    
#         likelihood_vec_t = [likelihood_vec]
    
#         cluster_assignment_t = Vector{Vector}(undef,T)
#         norm_cluster_prob_t = Vector{Vector}(undef,T)
        
#         for tp in 1:T
#             cells_ = getCustomersAtTimepoint(data_dict,tp)
    
#             # Sample States
#             phen_liklihood = [pdf.(likelihood_vec_t[tp],last(cell)[1]) for cell in cells_]
#             old_π_ht = π_ht[tp]
#             cluster_prob = [old_π_ht .* phen for phen in phen_liklihood]
#             norm_cluster_prob = softmax.(cluster_prob)
#             norm_cluster_prob_t[tp] = norm_cluster_prob
#             cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
            
#             cluster_assignment_t[tp] = cluster_assignment 
#             for c in 1:length(cluster_assignment)
#                 assgn = cluster_assignment[c]
#                 cell = cells_[c]
#                 push!(tp_cluster_assgn_t[tp][assgn],cell) # tp_cluster_assgn_t
#                 n_ht_t[tp][assgn]+=1 # n_ht_t
#             end
            
#             #Update π
#             if tp == 1
#                 n_ = [n_ht_t[tp][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             else
#                 n_ = [n_ht_t[tp][key] for key in 1:K_max] .+ [n_ht_t[tp-1][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             end
    
            
            
#             #Update λ
#             get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
#             curr_n_h = [n_ht_t[tp][key] for key in 1:K_max]
#             curr_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             if tp == 1
#                 prev_n_h = zeros(Int,K_max) # [n_ht_t[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = zeros(eltype(last(tp_cluster_assgn_t[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             else
#                 prev_n_h =  [n_ht_t[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             end
#             n_h = curr_n_h + prev_n_h 
#             suff_stats = curr_suff_stats + prev_suff_stats
#             upd_a = likelihood_prior.α .+ suff_stats
#             upd_b = likelihood_prior.θ .+ n_h
#             upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
#             # likelihood_vec_t[tp] = H.(upd_λ)
#             likelihood_vec = H.(upd_λ)
#             push!(likelihood_vec_t, likelihood_vec)#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
            
#             # Update θ_t 
#             ####################
#             #TODO: THIS Update yeilds a θ value in the  
#             # max_K_at_t = K_max #maximum(cluster_assignment)
#             # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
#             θ_t[tp] = θ_t[tp]
#         end
#         gibbs_tp_cluster_assgn[itr] = tp_cluster_assgn_t
#         gibbs_n_ht[itr] = n_ht_t
#         gibbs_likelihood_vec_t[itr] = likelihood_vec_t
#         gibbs_V_ht[itr] = V_ht
#         gibbs_π_ht[itr] = π_ht
#         gibbs_cluster_assignment[itr] = cluster_assignment_t
#         gibbs_norm_cluster_prob[itr] = norm_cluster_prob_t
#         gibbs_θ_t[itr] = θ_t
#         return θ_t
#     end
#     # a_θ,b_θ  = 0.1,0.1
#     # theta_prior = Gamma(a_θ,1/b_θ)
#     # θ_t = generate_θ_t(theta_prior,T)
#     # gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,1)
#     # gibbs_n_ht = Vector{Vector{Dict}}(undef,1)
#     # gibbs_likelihood_vec_t = Vector{Vector}(undef,1)
#     # gibbs_V_ht= Vector{Vector}(undef,1)
#     # gibbs_π_ht = Vector{Vector}(undef,1)
#     # gibbs_cluster_assignment= Vector{Vector}(undef,1)
#     # gibbs_norm_cluster_prob = Vector{Vector}(undef,1)
#     # gibbs_θ_t = Vector{Vector}(undef,1)
#     # π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
#     # ss = _2ndBlockedGibbs_AllTimeStep(1,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t) 
#     # ss2 = _2ndBlockedGibbs_AllTimeStep(1,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht, gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t, b_θ =1/25)
    
    
#     # test_itr = 10000
#     # gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,test_itr)
#     # gibbs_n_ht = Vector{Vector{Dict}}(undef,test_itr) 
#     # gibbs_likelihood_vec_t = Vector{Vector}(undef,test_itr)
#     # gibbs_V_ht= Vector{Vector}(undef,test_itr)
#     # gibbs_π_ht = Vector{Vector}(undef,test_itr)
#     # gibbs_cluster_assignment= Vector{Vector}(undef,test_itr)
#     # gibbs_norm_cluster_prob = Vector{Vector}(undef,test_itr)
#     # gibbs_θ_t = Vector{Vector}(undef,test_itr)
#     # π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
#     # for k in 1:test_itr
#     #     _2ndBlockedGibbs_AllTimeStep(k,data_dict,V_h_init,π_h_init,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht,gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t)
#     # end

#     # K_max = 25
#     # numBurnin = 10000
#     # maxNumItr = 10000
#      export _BlockedGibbs_AllTimeStepInit, _3rdBlockedGibbs_AllTimeStep, _4thBlockedGibbs_AllTimeStep
#     function _BlockedGibbs_AllTimeStepInit(data_dict;K_max = 5,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = 1,b=1,a_θ =0.1, b_θ =1/10)
#         gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
#         gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
#         theta_prior = theta_prior_dist(a_θ,1/b_θ)
#         θ_t =  [rand(theta_prior) for t in 1:T]
#         likelihood_prior =  likelihood_prior_dist(a,1/b)
#         likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
#         likelihood_vec = [likelihood_vec]
#         # V_h_init = Vector{Vector}(undef,T)
#         # π_h_init = Vector{Vector}(undef,T)
#         # for t in 1:T
#         #     V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
#         #     π_h_init[t] = stickbreak(V_h_init[t])
#         # end
#         π_h_init, V_h_init= generate_π_h_init(T,K_max,θ_t)
#         V_ht = Vector{Vector}(undef,T)
#         π_ht = Vector{Vector}(undef,T)
    
#         for tp in 1:T
#             cells_ = getCustomersAtTimepoint(data_dict,tp)
    
#             # Sample States
#             phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
    
#             cluster_prob = [π_h_init[tp] .* phen for phen in phen_liklihood]
#             norm_cluster_prob = softmax.(cluster_prob)
#             cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
#             for c in 1:length(cluster_assignment)
#                 assgn = cluster_assignment[c]
#                 cell = cells_[c]
#                 push!(gibbs_tp_cluster_assgn[tp][assgn],cell)
#                 gibbs_n_ht[tp][assgn]+=1
#             end
            
#             #Update π
#             if tp == 1
#                 n_ = [gibbs_n_ht[tp][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             else
#                 n_ = [gibbs_n_ht[tp][key] for key in 1:K_max] .+ [gibbs_n_ht[tp-1][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             end
    
            
            
#             #Update λ
#             get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
#             curr_n_h = [gibbs_n_ht[tp][key] for key in 1:K_max]
#             curr_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#             if tp == 1
#                 prev_n_h = zeros(Int,K_max) # [gibbs_n_ht[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = zeros(eltype(last(gibbs_tp_cluster_assgn[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#             else
#                 prev_n_h =  [gibbs_n_ht[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#             end
#             n_h = curr_n_h + prev_n_h 
#             suff_stats = curr_suff_stats + prev_suff_stats
#             upd_a = likelihood_prior.α .+ suff_stats
#             upd_b = likelihood_prior.θ .+ n_h
#             upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
#             push!(likelihood_vec, H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
            
#             # Update θ_t 
#             max_K_at_t = K_max #maximum(cluster_assignment)
#             # new_θ_t = deepcopy(θ_t)
#             # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
#             θ_t[tp] = θ_t[tp]
#         end
#         return θ_t, likelihood_vec, π_ht, V_ht,likelihood_prior,theta_prior 
#     end
    
#     # _BlockedGibbs_AllTimeStepInit(data_dict;K_max = 5,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = 1,b=1,a_θ =1, b_θ =1/5)
    
#     function _3rdBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t;K_max = 5,H = Poisson,a = 1,b=1,a_θ =1, b_θ =1/5)
#         ###### IN LOOP FORM ###########
#         if maximum(unique([last(key) for key in keys(data_dict)])) == length(unique([last(key) for key in keys(data_dict)]))
#             T = length(unique([last(key) for key in keys(data_dict)]))
#         else
#             T = maximum(unique([last(key) for key in keys(data_dict)])) 
#         end
    
#         #Updates from Previous iteration
#         θ_t =  deepcopy(gibbs_θ_t[itr-1])
#         likelihood_vec = deepcopy(gibbs_likelihood_vec_t[itr-1])
#         π_ht = deepcopy(gibbs_π_ht[itr-1])
#         V_ht = deepcopy(gibbs_V_ht[itr-1])
        
#         #iteration specific values initialization
#         tp_cluster_assgn_t = [Dict(h => [] for h in 1:K_max) for i in 1:T]
#         n_ht_t = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
#         likelihood_vec_t = []
#         cluster_assignment_t = Vector{Vector}(undef,T)
#         norm_cluster_prob_t = Vector{Vector}(undef,T)
        
#         for tp in 1:T
#             cells_ = getCustomersAtTimepoint(data_dict,tp)
    
#             # Sample States
#             phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
#             old_π_ht = π_ht[tp]
#             cluster_prob = [old_π_ht .* phen for phen in phen_liklihood]
#             norm_cluster_prob = softmax.(cluster_prob)
#             norm_cluster_prob_t[tp] = norm_cluster_prob
#             cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
            
#             cluster_assignment_t[tp] = cluster_assignment 
#             for c in 1:length(cluster_assignment)
#                 assgn = cluster_assignment[c]
#                 cell = cells_[c]
#                 push!(tp_cluster_assgn_t[tp][assgn],cell) # tp_cluster_assgn_t
#                 n_ht_t[tp][assgn]+=1 # n_ht_t
#             end
            
#             #Update π
#             if tp == 1
#                 n_ = [n_ht_t[tp][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             else
#                 n_ = [n_ht_t[tp][key] for key in 1:K_max] .+ [n_ht_t[tp-1][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             end
    
            
            
#             #Update λ
#             get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
#             curr_n_h = [n_ht_t[tp][key] for key in 1:K_max]
#             curr_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             if tp == 1
#                 prev_n_h = zeros(Int,K_max) # [n_ht_t[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = zeros(eltype(last(tp_cluster_assgn_t[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             else
#                 prev_n_h =  [n_ht_t[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             end
#             n_h = curr_n_h + prev_n_h 
#             suff_stats = curr_suff_stats + prev_suff_stats
#             upd_a = likelihood_prior.α .+ suff_stats
#             upd_b = likelihood_prior.θ .+ n_h
#             upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
#             # likelihood_vec_t[tp] = H.(upd_λ)
#             # likelihood_vec = H.(upd_λ)
#             push!(likelihood_vec_t,H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
            
#             # Update θ_t 
#             ####################
#             #TODO: THIS Update yeilds a θ value in the  
#             # max_K_at_t = K_max #maximum(cluster_assignment)
#             # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
#             θ_t[tp] = θ_t[tp]
#         end
#         gibbs_tp_cluster_assgn[itr] = tp_cluster_assgn_t
#         gibbs_n_ht[itr] = n_ht_t
#         gibbs_likelihood_vec_t[itr] = likelihood_vec_t
#         gibbs_V_ht[itr] = V_ht
#         gibbs_π_ht[itr] = π_ht
#         gibbs_cluster_assignment[itr] = cluster_assignment_t
#         gibbs_norm_cluster_prob[itr] = norm_cluster_prob_t
#         gibbs_θ_t[itr] = θ_t
#     end
    
#     # V_h only depends of the prev. time step; λ only depends on the current time step
#     function _4thBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t;K_max = 5,H = Poisson,a = 1,b=1,a_θ =1, b_θ =1/5)
#         ###### IN LOOP FORM ###########
#         if maximum(unique([last(key) for key in keys(data_dict)])) == length(unique([last(key) for key in keys(data_dict)]))
#             T = length(unique([last(key) for key in keys(data_dict)]))
#         else
#             T = maximum(unique([last(key) for key in keys(data_dict)])) 
#         end
    
#         #Updates from Previous iteration
#         θ_t =  deepcopy(gibbs_θ_t[itr-1])
#         likelihood_vec = deepcopy(gibbs_likelihood_vec_t[itr-1])
#         π_ht = deepcopy(gibbs_π_ht[itr-1])
#         V_ht = deepcopy(gibbs_V_ht[itr-1])
        
#         #iteration specific values initialization
#         tp_cluster_assgn_t = [Dict(h => [] for h in 1:K_max) for i in 1:T]
#         n_ht_t = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
#         likelihood_vec_t = []
#         cluster_assignment_t = Vector{Vector}(undef,T)
#         norm_cluster_prob_t = Vector{Vector}(undef,T)
        
#         for tp in 1:T
#             cells_ = getCustomersAtTimepoint(data_dict,tp)
    
#             # Sample States
#             phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
#             old_π_ht = π_ht[tp]
#             cluster_prob = [old_π_ht .* phen for phen in phen_liklihood]
#             norm_cluster_prob = softmax.(cluster_prob)
#             norm_cluster_prob_t[tp] = norm_cluster_prob
#             cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
            
#             cluster_assignment_t[tp] = cluster_assignment 
#             for c in 1:length(cluster_assignment)
#                 assgn = cluster_assignment[c]
#                 cell = cells_[c]
#                 push!(tp_cluster_assgn_t[tp][assgn],cell) # tp_cluster_assgn_t
#                 n_ht_t[tp][assgn]+=1 # n_ht_t
#             end
            
#             #Update π
#             if tp == 1
#                 n_ = [n_ht_t[tp][key] for key in 1:K_max]#sort(collect(keys(n_ht_t[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             else
#                 n_ = [n_ht_t[tp-1][key] for key in 1:K_max] #.+ [n_ht_t[tp][key] for key in 1:K_max] #sort(collect(keys(n_ht_t[tp])))
#                 V_ht[tp] =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht[tp] = stickbreak(V_ht[tp][1:K_max-1])
#             end
    
            
            
#             #Update λ
#             get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
#             curr_n_h = [n_ht_t[tp][key] for key in 1:K_max]
#             curr_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             if tp == 1
#                 prev_n_h = zeros(Int,K_max) # [n_ht_t[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = zeros(eltype(last(tp_cluster_assgn_t[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             else
#                 prev_n_h =  [n_ht_t[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = [isempty(get_cluster_stats(tp_cluster_assgn_t,tp,h)) ? 0 : sum(get_cluster_stats(tp_cluster_assgn_t,tp,h))[1] for h in 1:K_max]
#             end
#             n_h = curr_n_h + prev_n_h 
#             suff_stats = curr_suff_stats #+ prev_suff_stats
#             upd_a = likelihood_prior.α .+ suff_stats
#             upd_b = likelihood_prior.θ .+ n_h
#             upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
#             # likelihood_vec_t[tp] = H.(upd_λ)
#             # likelihood_vec = H.(upd_λ)
#             push!(likelihood_vec_t,H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
            
#             # Update θ_t 
#             ####################
#             #TODO: THIS Update yeilds a θ value in the  
#             # max_K_at_t = K_max #maximum(cluster_assignment)
#             # θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[tp][1:max_K_at_t]))))
#             θ_t[tp] = θ_t[tp]
#         end
#         gibbs_tp_cluster_assgn[itr] = tp_cluster_assgn_t
#         gibbs_n_ht[itr] = n_ht_t
#         gibbs_likelihood_vec_t[itr] = likelihood_vec_t
#         gibbs_V_ht[itr] = V_ht
#         gibbs_π_ht[itr] = π_ht
#         gibbs_cluster_assignment[itr] = cluster_assignment_t
#         gibbs_norm_cluster_prob[itr] = norm_cluster_prob_t
#         gibbs_θ_t[itr] = θ_t
#         return θ_t
#     end


#     # K_max = 5
#     # a = 6
#     # b= 1
#     # a_θ = 1
#     # b_θ =1/5
#     # numBurnin = 10000
#     # maxNumItr = 10000
#     # maxNumItr = maxNumItr + numBurnin
    
#     # gibbs_tp_cluster_assgn = Vector{Vector{Dict}}(undef,maxNumItr+1)
#     # gibbs_n_ht = Vector{Vector{Dict}}(undef,maxNumItr+1) 
#     # gibbs_likelihood_vec_t = Vector{Vector}(undef,maxNumItr+1)
#     # gibbs_V_ht= Vector{Vector}(undef,maxNumItr+1)
#     # gibbs_π_ht = Vector{Vector}(undef,maxNumItr+1)
#     # gibbs_cluster_assignment= Vector{Vector}(undef,maxNumItr+1)
#     # gibbs_norm_cluster_prob = Vector{Vector}(undef,maxNumItr+1)
#     # gibbs_θ_t = Vector{Vector}(undef,maxNumItr+1)
    
    
    
#     # rand_init_θ_t, rand_init_likelihood_vec, rand_init_π_ht, rand_init_V_ht,likelihood_prior,theta_prior  = _BlockedGibbs_AllTimeStepInit(data_dict;K_max = K_max,H = Poisson,theta_prior_dist = Gamma, likelihood_prior_dist=Gamma, a = a,b=b,a_θ =a_θ, b_θ =b_θ)
    
#     # gibbs_likelihood_vec_t[1] = rand_init_likelihood_vec
#     # gibbs_V_ht[1] =  rand_init_V_ht
#     # gibbs_π_ht[1] = rand_init_π_ht
#     # gibbs_θ_t[1] = 0.001ones(T)#rand_init_θ_t
    
#     # # _3rdBlockedGibbs_AllTimeStep(2,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,K_max = K_max)
    
#     # for step in 2:maxNumItr+1
#     #     itr = step
#     #     _3rdBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,a = a,b=b,K_max = K_max)
#     #     # _4thBlockedGibbs_AllTimeStep(itr,data_dict,likelihood_prior,gibbs_tp_cluster_assgn,gibbs_n_ht,gibbs_likelihood_vec_t,gibbs_V_ht, gibbs_π_ht,gibbs_cluster_assignment,gibbs_norm_cluster_prob,gibbs_θ_t,a = a,b=b,K_max = K_max)
#     # end
    

#     export _1stBlockedGibbs_AllTimeStep
#     function _1stBlockedGibbs_AllTimeStep(data_dict)
#         ###### IN LOOP FORM ###########
#         K_max = 25
#         gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
#         gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
#         rand_gibbs_init_cluster = [ rand(1:K_max,c) for c in 3 .* Cₜ]
#         H = Poisson
#         a,b = 1,1
#         a_θ, b_θ = 1,0.1
#         theta_prior = Gamma(a_θ,1/b_θ)
#         θ_t =  [rand(theta_prior) for t in 1:T]
#         likelihood_prior = Gamma(a,1/b)
#         likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
#         likelihood_vec_init = deepcopy(likelihood_vec)
#         likelihood_vec = [likelihood_vec]
#         V_h_init = Vector{Vector}(undef,T)
#         π_h_init = Vector{Vector}(undef,T)
#         for t in 1:T
#             V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
#             π_h_init[t] = stickbreak(V_h_init[t])
#         end
    
#         for tp in 1:T
#             cells_ = getCustomersAtTimepoint(data_dict,tp)
    
#             # Sample States
#             phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
    
#             cluster_prob = [π_h_init[tp] .* phen for phen in phen_liklihood]
#             norm_cluster_prob = softmax.(cluster_prob)
#             cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
#             for c in 1:length(cluster_assignment)
#                 assgn = cluster_assignment[c]
#                 cell = cells_[c]
#                 push!(gibbs_tp_cluster_assgn[tp][assgn],cell)
#                 gibbs_n_ht[tp][assgn]+=1
#             end
            
#             #Update π
#             if tp == 1
#                 n_ = [gibbs_n_ht[tp][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
#                 V_ht =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht = stickbreak(V_ht[1:K_max-1])
#             else
#                 n_ = [gibbs_n_ht[tp][key] for key in 1:K_max] .+ [gibbs_n_ht[tp-1][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
#                 V_ht =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#                 π_ht = stickbreak(V_ht[1:K_max-1])
#             end
    
            
            
#             #Update λ
#             get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h])
#             curr_n_h = [gibbs_n_ht[tp][key] for key in 1:K_max]
#             curr_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#             if tp == 1
#                 prev_n_h = zeros(Int,K_max) # [gibbs_n_ht[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = zeros(eltype(last(gibbs_tp_cluster_assgn[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#             else
#                 prev_n_h =  [gibbs_n_ht[tp-1][key] for key in 1:K_max]
#                 prev_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#             end
#             n_h = curr_n_h + prev_n_h 
#             suff_stats = curr_suff_stats + prev_suff_stats
#             upd_a = likelihood_prior.α .+ suff_stats
#             upd_b = likelihood_prior.θ .+ n_h
#             upd_λ =  rand.(Gamma.(upd_a, upd_b  ))
#             push!(likelihood_vec, H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ]) 
            
#             # Update θ_t 
#             max_K_at_t = K_max #maximum(cluster_assignment)
#             # new_θ_t = deepcopy(θ_t)
#             θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_ht[1:max_K_at_t]))))
#         end
#         return θ_t
#     end
#     #vv = _1stBlockedGibbs_AllTimeStep(data_dict) 
    
#     # Just the first Time step!
#     export  _1stBlockedGibbs_1stTimeStep
#     function _1stBlockedGibbs_1stTimeStep(data_dict)
#         K_max = 25
#         gibbs_tp_cluster_assgn = [Dict(h => [] for h in 1:K_max) for i in 1:T]
#         gibbs_n_ht = [Dict(h => 0 for h in 1:K_max) for i in 1:T]
#         rand_gibbs_init_cluster = [ rand(1:K_max,c) for c in 3 .* Cₜ]
#         H = Poisson
#         a,b = 1,1
#         a_θ, b_θ = 1,0.1
#         theta_prior = Gamma(a_θ,1/b_θ)
#         θ_t =  [rand(theta_prior) for t in 1:T]
#         likelihood_prior = Gamma(a,1/b)
#         likelihood_vec = [H(rand(likelihood_prior)) for h in 1:K_max]
#         likelihood_vec_init = deepcopy(likelihood_vec)
#         likelihood_vec = [likelihood_vec]
#         V_h_init = Vector{Vector}(undef,T)
#         π_h_init = Vector{Vector}(undef,T)
#         for t in 1:T
#             V_h_init[t] = rand(Beta(1,θ_t[t]),K_max-1)
#             π_h_init[t] = stickbreak(V_h_init[t])
#         end
#         tp = 1
#         cells_ = getCustomersAtTimepoint(data_dict,tp)
#         phen_liklihood = [pdf.(likelihood_vec[tp],last(cell)[1]) for cell in cells_]
#         cluster_prob = [π_h_init[tp] .* phen for phen in phen_liklihood]
#         norm_cluster_prob = softmax.(cluster_prob)
#         cluster_assignment =  rand.(Categorical.(norm_cluster_prob))
#         for c in 1:length(cluster_assignment)
#             assgn = cluster_assignment[c]
#             cell = cells_[c]
#             push!(gibbs_tp_cluster_assgn[tp][assgn],cell)
#             gibbs_n_ht[tp][assgn]+=1
#         end

#         #Update π
#         n_ = [gibbs_n_ht[tp][key] for key in 1:K_max]#sort(collect(keys(gibbs_n_ht[tp])))
#         V_h_1 =rand.(Beta.(1 .+ n_,θ_t[tp] .+ cumsum(n_) ))
#         π_h_1 = stickbreak(V_h_1[1:K_max-1])

#         #Update λ
#         curr_n_h = [gibbs_n_ht[tp][key] for key in 1:K_max]
#         prev_n_h = zeros(Int,K_max) # [gibbs_n_ht[tp-1][key] for key in 1:K_max]
#         n_h = curr_n_h + prev_n_h
#         get_cluster_stats(cluster_assgn,t,h) = last.(cluster_assgn[t][h]) 
#         curr_suff_stats = [isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#         prev_suff_stats = zeros(eltype(last(gibbs_tp_cluster_assgn[tp][cluster_assignment[1]][1])),K_max) #[isempty(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h)) ? 0 : sum(get_cluster_stats(gibbs_tp_cluster_assgn,tp,h))[1] for h in 1:K_max]
#         # [cells_[1] for h in 1:K_max]
#         suff_stats = curr_suff_stats + prev_suff_stats
#         upd_a = likelihood_prior.α .+ suff_stats
#         upd_b = likelihood_prior.θ .+ n_h
#         upd_λ =  rand.(Gamma.(upd_a, upd_b))
#         push!(likelihood_vec, H.(upd_λ))#push!(likelihood_vec, [H(new_lam) for new_lam in upd_λ ])

#         # Update θ 
#         max_K_at_t = K_max #maximum(cluster_assignment)
#         new_θ_t = deepcopy(θ_t)
#         new_θ_t[tp] = rand(Gamma(a_θ + max_K_at_t - 1, 1/(b_θ) + sum(log.(1 .- V_h_1[1:max_K_at_t]))))
#     end
#     # _1stBlockedGibbs_1stTimeStep(data_dict)

#     export  _firstAttempt
#     function _firstAttempt(data_dict)
#         DP_α = 1
#         a,b = 1,1
#         lambda_prior_dist = Gamma(a,b) # alpha = a = shape, b = rate = 1/scale = 1/beta
#         # λ = rand(lambda_prior_dist)
#         likelihood_vec = Vector()
#         H = Poisson
#         # push!(likelihood_vec,H(λ))
    
#         # For Time Point 1!
#         tp = 1
    
#         key_array_tp1 = collect(keys(data_dict))
#         key_array_tp1 = key_array_tp1[last.(key_array_tp1) .== tp]
#         key_array_tp1 = sort(key_array_tp1, by = first)
#         _tp1 = [(first(key),data_dict[key]) for key in key_array_tp1]
    
#         tp1_cluster_assgn = Dict()
#         n_ht = Dict(t => Dict() for t in 1:T) # Number of cells in cluster h in this current (they use iteration to refer to simulation iterations (like in Gibbs) but I have to consider time point too!)
#         itr = 1
#         curr_k_max = 0
#         for cell in _tp1
#             if itr == 1
#                 # println("here")
#                 tp1_cluster_assgn[1] = [cell]
#                 n_ht[tp][1] = 1
#                 curr_k_max += 1
#                 # println("Current K_max = $curr_k_max")
#                 λ = rand(Gamma(lambda_prior_dist.α + last(cell)[1], 1 + 1/(lambda_prior_dist.θ))) # Only works for One gene case!
#                 push!(likelihood_vec,H(λ))
#             else
#                 # println("now here")
#                 # println("Current K_max = $curr_k_max")
#                 phen_val = last(cell)[1]
#                 L_vals = [pdf(L,phen_val) for L in likelihood_vec]
#                 n_h = [n_ht[tp][key] for key in sort(collect(keys(n_ht[tp])))]
#                 cluster_prob = L_vals .* n_h
#                 ###### NOW IM GUESSING ######
#                 new_cluster_prob = DP_α / (itr - 1 + DP_α)
#                 new_cluster_indicator = rand(Bernoulli(new_cluster_prob))
#                 if isone(new_cluster_indicator)
#                     curr_k_max += 1 # Increase total Number
#                     # println("WE BIG NOW")
#                     # println("Current K_max = $curr_k_max")
                    
#                     tp1_cluster_assgn[curr_k_max] = [cell] # Add to Cluster
#                     n_ht[tp][curr_k_max] = 1 # Increase Cluster Count
#                     λ = rand(Gamma(lambda_prior_dist.α + last(cell)[1], 1 + 1/(lambda_prior_dist.θ))) #rand(lambda_prior_dist) ## Sample new parameters
#                     push!(likelihood_vec,H(λ))
#                 else
#                     norm_cluster_prob =  softmax(cluster_prob)#cluster_prob ./sum(cluster_prob)
#                     cluster_val = rand(Categorical(norm_cluster_prob)) # Determine Cluster
#                     push!(tp1_cluster_assgn[cluster_val],cell)# Add to Cluster
#                     n_ht[tp][cluster_val] += 1 # Increase Cluster Count
#                     n_ = n_ht[tp][cluster_val] 
#                     suff_stat = sum(last.(tp1_cluster_assgn[cluster_val]))[1]
#                     λ = rand(Gamma(lambda_prior_dist.α + suff_stat, n_ + 1/(lambda_prior_dist.θ)))  ## Sample new parameters
#                     likelihood_vec[cluster_val] = H(λ)
#                 end
#             end
#             # println("We outside")
#             # println("Current K_max = $curr_k_max")
#             itr += 1
#         end
#         # println("After Time 1, The current K is $curr_k_max")
#         tp_cluster_assgn = [tp1_cluster_assgn]
#         _tp_vec = [_tp1]
#         ka_= collect(keys(data_dict))
#         cells_in_t = [length(_tp1)]
#         for tp in 2:T
#             key_array_ = ka_[last.(ka_) .== tp]
#             key_array_ = sort(key_array_, by = first)
#             _tp = [(first(key),data_dict[key]) for key in key_array_]
#             push!(_tp_vec,_tp) #Just keeping track for dubuging purposes
#             push!(cells_in_t,length(_tp))
#             new_tp_cluster_assgn = Dict()
#             for i in 1:curr_k_max
#                 n_ht[tp][i] = 0
#                 new_tp_cluster_assgn[i] = []
#             end
            
#             push!(tp_cluster_assgn,new_tp_cluster_assgn)
#             itr = 1
#             for cell in _tp
#                 phen_val = last(cell)[1]
#                 L_vals = [pdf(L,phen_val) for L in likelihood_vec]
#                 n_h = [n_ht[tp][key] for key in sort(collect(keys(n_ht[tp])))] .+ [n_ht[tp-1][key] for key in sort(collect(keys(n_ht[tp-1])))]
#                 cluster_prob = L_vals .* n_h
#                 ###### NOW IM GUESSING ######
#                 new_cluster_prob = DP_α / (cells_in_t[tp-1]  +itr - 1 + DP_α)
#                 new_cluster_indicator = rand(Bernoulli(new_cluster_prob))
#                 if isone(new_cluster_indicator)
#                     curr_k_max += 1 # Increase total Number
#                     tp_cluster_assgn[tp][curr_k_max]= [cell]
#                     n_ht[tp][curr_k_max] = 1 # Increase Cluster Count
#                     for i in 1:tp-1
#                         n_ht[i][curr_k_max] = 0 # Update Kmax for past time points
#                         tp_cluster_assgn[i][curr_k_max] = []  # Update Kmax for past time points
#                     end
#                     λ = rand(Gamma(lambda_prior_dist.α + last(cell)[1], 1 + 1/(lambda_prior_dist.θ))) #rand(lambda_prior_dist) ## Sample new parameters
#                     push!(likelihood_vec,H(λ))
#                 else
#                     norm_cluster_prob = softmax(cluster_prob)#cluster_prob ./sum(cluster_prob)
#                     cluster_val = rand(Categorical(norm_cluster_prob)) # Determine Cluster
#                     push!(tp_cluster_assgn[tp][cluster_val],cell)# Add to Cluster
#                     n_ht[tp][cluster_val] += 1 # Increase Cluster Count
#                     n_ = n_ht[tp][cluster_val] + n_ht[tp-1][cluster_val]
#                     old_suff_stat_bool = isempty(last.(tp_cluster_assgn[tp-1][cluster_val]))
#                     if old_suff_stat_bool
#                         old_suff_stat = 0 
#                     else
#                         old_suff_stat = sum(last.(tp_cluster_assgn[tp-1][cluster_val]))[1]
#                     end
#                     new_suff_stat = sum(last.(tp_cluster_assgn[tp][cluster_val]))[1] 
#                     suff_stat = new_suff_stat + old_suff_stat
#                     λ = rand(Gamma(lambda_prior_dist.α + suff_stat, n_ + 1/(lambda_prior_dist.θ)))  ## Sample new parameters
#                     likelihood_vec[cluster_val] = H(λ)
#                 end
    
    
#                 itr +=1
#             end
#             # println("After Time $tp, The current K is $curr_k_max")
#         end
#         return curr_k_max,n_ht,tp_cluster_assgn,likelihood_vec
#     end
#     # TODO: ADD THESE TESTS
#     # sim_time = 1000
#     # k_max_simtime = Vector{Int}(undef,sim_time)
#     # for s in 1:sim_time
#     #     k_max_simtime[s],_,_,_ =_firstAttempt(data_dict)
#     # end
#     # histogram(k_max_simtime)
    
    
#     function seqCRP(customers,tp_cluster_assgn,n_ht, curr_k_max;ft=false)
#         itr = 1
    
#     end
#     function _secondAttempt()
        
#     end

#     function _1stBlockedGibbs()
    
#     end

# end


# include("synDataPreprocess.jl")
# using  .syntheticDataPreprocessing

# using Random
# using Distributions
# using Turing
# using Turing.RandomMeasures: stickbreak, DirichletProcess, StickBreakingProcess, ChineseRestaurantProcess
# using StatsBase, StatsFuns, StatsModels, StatsPlots, Statistics
# using Test
# import Debugger
# using CSV,DataFrames

