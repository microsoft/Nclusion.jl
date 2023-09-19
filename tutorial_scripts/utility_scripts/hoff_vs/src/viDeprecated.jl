
function variational_inference_MvGaussian(x, G,K,γ,α0,m0,mu0,num_iter, num_local_iter; rand_init = false,ep = 0.001)
    T = length(x)
    v0 = G + 1
    if rand_init
        rtik,Ntk,θ_hat,m0k_vec,Mk_vec,W0,Wk_vec,vk_vec, rho_hat, omega_hat = rand_init_params_MvGaussian(T,G,K,γ,v0,m0,mu0;ep= ep)
    else
        rtik,Ntk,θ_hat,m0k_vec,Mk_vec,W0,Wk_vec,vk_vec, rho_hat, omega_hat = init_params_MvGaussian(T,G,K,γ,v0,m0,mu0)
    end
    # @debug W0
    # updateVariationalZdist(x,G,m0k_vec,Mk_vec,Wk_vec,vk_vec,θ_hat)
    num_iter = num_iter
    elbo = Vector{Float64}(undef,num_iter)
    for iter in 1:num_iter
        # Local updates

        # @debug "Iteration $iter"
        # if iter == 3
        #     @debug "start $iter" 
        #     @debug Wk_vec
        # end
        # if iter == 2
        #     @debug "start $iter" 
        #     @debug Wk_vec

        # end
        for local_iter in 1:num_local_iter
            # if iter == 3
            #     @debug rtik 
            # end
            rtik = updateVariationalZdist_MvGaussian(x,G,m0k_vec,Mk_vec,Wk_vec,vk_vec,θ_hat)
            Ntk = calc_Ntk_MvGaussian(rtik)
            θ_hat = updateVariationalπdist_MvGaussian(α0,γ,K,Ntk)
            
        end
        # if iter == 4
        #     # @debug "end $iter"
        #     # @debug Wk_vec
        #     return
        # end

        # if iter == 2
        #     @debug rtik 
        #     @debug Ntk
        # end
        xbar_k = update_xbar_k_MvGaussian(rtik,Ntk,x)
        Sk_vec = update_Sk_MvGaussian(x,xbar_k,Ntk,rtik)

        
        
        # Global updates
            # m0k_vec
        m0k_vec = update_m0k_MvGaussian(m0k_vec, Ntk)
            # vk_vec
        vk_vec = update_vk_MvGaussian(vk_vec, Ntk)
            # Mk_vec
        Mk_vec = update_Mk_MvGaussian(Mk_vec,m0k_vec,mu0,m0, Ntk,xbar_k)
            # Wk_vec
        Wk_vec = update_Wk_MvGaussian(xbar_k,mu0,m0,Ntk,Sk_vec,W0)
            # Tk
        T_k = calc_Tk(θ_hat)
        # LowerBound_LG_unconstrained(c,d,G,γ,α_0,T_k)
        # c_hat = StatsFuns.logit.(rho_hat)
        # d_hat = log.(omega_hat)
        # cd_vec_args = [c_hat , d_hat]
        # cd_vec_args0  = permutedims(reduce(hcat,cd_vec_args))
        # LB_LG_unconstrained = LowerBound_LG2__unconstrained_closure(G,γ,α0,T_k)
        # gg_uncon! = g_unconstrained_closure!(G,γ,α0,T_k)
        # lb_lg_results = optimize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = 1_000))
        # # @debug Optim.converged(lb_lg_results)
        # rho_hat = sigmoid.(lb_lg_results.minimizer[1,:])
        # omega_hat = exp.(lb_lg_results.minimizer[2,:])

        # c_hat = StatsFuns.logit.(rho_hat)
        # d_hat = log.(omega_hat)
        c_hat,d_hat = update_rho_omega_hat_MvGaussian(rho_hat,omega_hat,G,γ,α0,T_k;optim_max_iter=1000)
        #Calculate Elbo
        #calc_Hz(rtik) +
        # Have to that the negative of LowerBound_LG_unconstrained to get the max value and undo the steps needed for the optimization. Optim.jl only finds minimum so I that find the min of the neg. LowerBound_LG_unconstrained to get the max
        # @debug Ntk
        # @debug -1 * LowerBound_LG_unconstrained(c_hat,d_hat,G,γ,α0,T_k)
        elbo[iter] =   calc_Hz(rtik) + -1 * LowerBound_LG_unconstrained_MvGaussian(c_hat,d_hat,G,γ,α0,T_k)
        # if iter == 2
        #     @debug "end $iter"
        #     @debug Wk_vec
        # end
    end
    return elbo,rtik,Ntk
end
function update_rho_omega_hat_MvGaussian(rho_hat,omega_hat,G,γ,α0,T_k;optim_max_iter=1000)
    c_hat = StatsFuns.logit.(rho_hat)
    d_hat = log.(omega_hat)
    cd_vec_args = [c_hat , d_hat]
    cd_vec_args0  = permutedims(reduce(hcat,cd_vec_args))
    LB_LG_unconstrained = LowerBound_LG2__unconstrained_closure_MvGaussian(G,γ,α0,T_k)
    gg_uncon! = g_unconstrained_closure_MvGaussian!(G,γ,α0,T_k)
    lb_lg_results = optimize(LB_LG_unconstrained,gg_uncon!,cd_vec_args0, LBFGS(linesearch=Optim.LineSearches.BackTracking()), Optim.Options(iterations = optim_max_iter))
    # @debug Optim.converged(lb_lg_results)
    new_rho_hat = sigmoid.(lb_lg_results.minimizer[1,:])
    new_omega_hat = exp.(lb_lg_results.minimizer[2,:])

    new_c_hat = StatsFuns.logit.(new_rho_hat)
    new_d_hat = log.(new_omega_hat)
    return new_c_hat,new_d_hat
end
function init_params_MvGaussian(T,G,K,γ,v0,m0,mu0)
    θ_hat =  [ones(K+1) ./ K for t in 1:T]#[βk_expected_value(γ,K) for t in 1:T]
    m0k_vec =  m0 .* ones(K)
    Mk_vec = [ mu0 .* ones(G) for k in 1:K]
    W0 = Matrix{Float64}(I,G,G)#rand_pos_def_scale_mat(G)
    Wk_vec =  [W0 for k in 1:K] # [Matrix{Float64}(I,G,G) for k in 1:Ktrue] # [Symmetric(rand(G,G)) for k in 1:Ktrue]
    vk_vec = v0 .* ones(K) .+ 1.0 
    rho_hat = 0.75 .* ones(K)
    omega_hat = 2.0 .* ones(K)
    rtik= nothing
    Ntk = nothing
    return rtik,Ntk,θ_hat,m0k_vec,Mk_vec,W0,Wk_vec,vk_vec, rho_hat, omega_hat
end
function rand_init_params_MvGaussian(T,G,K,γ,v0,m0,mu0;ep = 0.001)
    θ_hat =  [βk_expected_value(γ,K) for t in 1:T]
    m0k_vec =  m0 .* rand(K) .+ ep
    Mk_vec = [ mu0 .* rand(G) .+ ep for k in 1:K]
    W0 = rand_pos_def_scale_mat(G)
    Wk_vec =  [W0 for k in 1:K] # [Matrix{Float64}(I,G,G) for k in 1:Ktrue] # [Symmetric(rand(G,G)) for k in 1:Ktrue]
    vk_vec = v0 .* ones(K) .+ 1.0 
    rho_hat = rand(K)#0.75 .* ones(K)
    omega_hat = 10.0 .* rand(K)
    rtik= nothing
    Ntk = nothing
    return rtik,Ntk,θ_hat,m0k_vec,Mk_vec,W0,Wk_vec,vk_vec, rho_hat, omega_hat
end

function fake_mvGausssian_data_for_testing(G,C_t,Ktrue;μ =nothing,Λ = nothing, mix_prob =nothing,same_prob_t = true)
    T = length(C_t)
    if isnothing(μ)
        μ = Vector{Vector{Float64}}(undef,Ktrue)
        power = [floor(k/2) for k in 1:Ktrue]
        neg = [(-1.0)^(k) for k in 1:Ktrue]
        neg[1] = 1.0
        for k in 1:Ktrue
            μ[k] = neg[k] .* 10.0 .^(power[k]) .* rand(G)
        end
    end
    if isnothing(Λ)
        Λ = [Matrix{Float64}(I, G,G) for k in 1:Ktrue]
    end
    if isnothing(mix_prob)
        if same_prob_t
            mix_prob = [ones(Ktrue)./Ktrue for t in 1:T]
        else
            mix_prob = rand(Dirichlet(Ktrue, 1.0))
        end
    end
    # assgn_mix_model = Vector{Vector}(undef,T)
    assgn_mix_model = [MixtureModel(Dirac.(collect(1:Ktrue)),mix_prob[T]) for t in 1:T ]
    z = [rand(assgn_mix_model[t],C_t[t]) for t in 1:T]
    x = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        x[t] = Vector{Vector{Float64}}(undef,C_t[t])
        for c in 1:C_t[t]
            x[t][c] = rand(MultivariateNormal(μ[z[t][c]],Λ[z[t][c]]))
        end
    end
    return x,z,mix_prob 
end

function fake_mvGausssian_inferenceVariables_for_testing(G,Ktrue)
    m0k_vec = ones(Ktrue)
    Mk_vec = [ones(G) for k in 1:Ktrue]
    # want to make Positive definite Matrix
    
    Wk_vec =  [rand_pos_def_scale_mat(G) for k in 1:Ktrue] # [Matrix{Float64}(I,G,G) for k in 1:Ktrue] # [Symmetric(rand(G,G)) for k in 1:Ktrue]
    vk_vec = G .* ones(Ktrue) .+ 1 
    θ_hat = [append!(ones(Ktrue) ./ Ktrue, 0.0) for t in 1:T]
    return m0k_vec,Mk_vec,Wk_vec,vk_vec,θ_hat
end
function rand_pos_def_scale_mat(G)
    A = randn(G,G);
    A = A'*A;
    A = (A + A')/2
    return A
end
# Update r
function update_rtik_MvGaussian(G, e_μ_Λ, e_log_Λ,e_log_π)
    Gln2π = G * log(2 * π)
    
    p_tilde_k = e_log_π .+ 1/2 .* (e_log_Λ .- Gln2π  .-  e_μ_Λ)
    # @debug "e_log_π: $e_log_π"
    # @debug "e_log_Λ: $e_log_Λ"
    # @debug "e_μ_Λ: $e_μ_Λ"
    # @debug "p_tilde_k: $p_tilde_k"
    # @debug "Before p_tilde_k: $p_tilde_k"
    # p_tilde_k = p_tilde_k .- maximum(p_tilde_k)
    # @debug "After p_tilde_k: $p_tilde_k"
    val_sum = StatsFuns.logsumexp(p_tilde_k)
    rtik = exp.(p_tilde_k .- val_sum)
    # @debug "$rtik"
    # @debug "Before ritk: $rtik"
    
    shifted_rtik = rtik .+ eps(1.0)
    rtik = shifted_rtik ./ sum(shifted_rtik)
    
    # @debug "After ritk: $rtik"
    return rtik
end

function μk_Λk_expected_value_MvGaussian(G,m0k, vk,x_ti,Mk,Wk)
    diffence_vector = x_ti .- Mk
    AtBA = transpose(diffence_vector) * Wk * diffence_vector
    e_μk_Λk = G * 1/m0k + vk * AtBA
    return e_μk_Λk 
end

function μ_Λ_expected_value_MvGaussian(G,m0k_vec, vk_vec,x_ti,Mk_vec,Wk_vec)
    K = length(m0k_vec)
    e_μ_Λ = Vector{Float64}(undef,K)
    for k in 1:K
        m0k = m0k_vec[k]
        vk = vk_vec[k]
        Mk = Mk_vec[k]
        Wk = Wk_vec[k]
        e_μ_Λ[k] =  μk_Λk_expected_value_MvGaussian(G,m0k, vk,x_ti,Mk,Wk)
    end
    return e_μ_Λ
end

function log_Λk_expected_value_MvGaussian(G,K,vk,Wk)
    e_log_Λk = Vector{Float64}(undef,K)
    Gln2 = G * log(2)
    det_Wk = det(Wk)
    e_log_Λk = 0.0
    for g in 1:G
        e_log_Λk += digamma(1\2 * (vk + 1 - g) ) + Gln2 + log(det_Wk) # Make sure we need to add Gln2  and  log(det_Wk) in each iteration
    end
    # param_val_dict["log_Λk_expected_value"] = e_log_Λk
    return e_log_Λk
end
function log_Λ_expected_value_MvGaussian(G,vk_vec,Wk_vec)
    K = length(vk_vec)
    e_log_Λ = Vector{Float64}(undef,K)
    for k in 1:K
        vk = vk_vec[k]
        Wk = Wk_vec[k]
        e_log_Λ[k] = log_Λk_expected_value_MvGaussian(G,K,vk,Wk)
    end
    return e_log_Λ
end

function log_π_expected_value_MvGaussian(θ_hat)
    K = Int64(length(θ_hat)-1)
    digamma_sum = digamma(sum(θ_hat[1:K]))
    e_log_π = digamma.(θ_hat[1:K]) .- digamma_sum
    # param_val_dict["log_π_expected_value"] = e_log_π
    return e_log_π
end
function updateVariationalZdist_MvGaussian(x,G,m0k_vec,Mk_vec,Wk_vec,vk_vec,θ_hat)
    T = length(x)
    rtik = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        # @debug "updateVariationalZdist Time $t"
        cells = C_t[t]
        e_log_Λ = log_Λ_expected_value_MvGaussian(G,vk_vec,Wk_vec)
        e_log_π = log_π_expected_value_MvGaussian(θ_hat[t])
        rtik[t] = Vector{Vector{Float64}}(undef,cells)
        for i in 1:cells
            e_μ_Λ = μ_Λ_expected_value_MvGaussian(G,m0k_vec, vk_vec,x[t][i],Mk_vec,Wk_vec)
            
            rtik[t][i] = update_rtik_MvGaussian(G, e_μ_Λ, e_log_Λ,e_log_π)
        end
        # println(t)
        # println(size(e_log_Λ))
        # println(size(e_log_π))
        # println(size(e_μ_Λ))
        # println("&&&&&")
        
    end
    # println("TTTTTT")
    # println(size(rtik[1][1]))
    return rtik
end


function calc_Ntk_MvGaussian(rtik)
    T = length(rtik)
    C_t = [length(el) for el in rtik]
    K = length(rtik[1][1])
    Ntk =  Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        Nt_vec = zeros(Float64,K+1) # K+1 is supposed to be the >K which aggregates the infinity extra unused topics
        Nt_vec[1:K] =  sum(rtik[t]) #1 .+ sum(rtik[t]) ##### BIGGGGGG BAUSE
        # for k in 1:K
        #     Nt_sum = 0.0
        #     for i in C_t[t]
        #         Nt_sum += rtik[t][i][k]
        #     end
        #     Nt_vec[k] = Nt_sum
        # end
        Ntk[t] = Nt_vec
    end
    return Ntk
end



function βk_expected_value_MvGaussian(γ,K) # Check this!
    e_uk = 1/(1+γ)
    e_βk = Vector{Float64}(undef,K+1)
    for k in 1:K
        e_βk[k] = 1/γ * e_uk^(k-1) - e_uk^k
    end
    gamma_plus =Float64( 1.0+γ)
    inv_gamma = Float64(1/γ)
    e_βk[K+1] = (inv_gamma*gamma_plus^(-K)  - gamma_plus^(-K-1)) / (1 -  gamma_plus^(-1))
    # e_βk = [1/γ * e_uk^(k-1) - e_uk^k  for k in 1:K]
    return e_βk
end
function update_theta_hat_MvGaussian(α0,γ,K,N_k)
    e_βk =  βk_expected_value_MvGaussian(γ,K)
    θ_hat = α0 * e_βk  .+ N_k
    return θ_hat
end

function updateVariationalπdist_MvGaussian(α0,γ,K,Ntk)
    T = length(Ntk)
    θ_hat = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        N_k = Ntk[t]
        θ_hat[t] = update_theta_hat_MvGaussian(α0,γ,K,N_k)
    end
    return θ_hat
end


function update_m0k_MvGaussian(m0k_vec, Ntk)
    K = length(m0k_vec)
    Nk = sum(Ntk)
    m0k_vec = m0k_vec .+ Nk[1:K]
    return m0k_vec
end

function update_vk_MvGaussian(vk_vec, Ntk)
    K = length(vk_vec)
    Nk = sum(Ntk)
    vk_vec = vk_vec .+ Nk[1:K]
    return vk_vec
end


function update_xbar_k_MvGaussian(rtik,Ntk,x)
    K = length(Ntk[1]) - 1
    C_t = [length(el) for el in x]
    T = length(Ntk)
    Nk = sum(Ntk)[1:K]
    xbar_k = Vector{Vector{Float64}}(undef,K)
    
    for k in 1:K
        Nk_scale = 1/Nk[k]
        xbar_k[k] = Nk_scale * sum([rtik[t][c][k] .* x[t][c] for t in 1:T for c in 1:C_t[t] ]) # Check this one too!
    end
    return xbar_k
end


function update_Mk_MvGaussian(Mk_vec,m0k_vec,mu0,m0, Ntk,xbar_k)
    K = length(Mk_vec)
    Nk = sum(Ntk)[1:K]
    Nk_prod_xbar_k =  Nk .* xbar_k
    mu0_prood_m0 = mu0 * m0 .* ones(Float64,K) 
    Mk_vec = 1 ./ m0k_vec .* map(g -> g[1] .+ g[2] , zip(mu0_prood_m0, Nk_prod_xbar_k))
    return m0k_vec
end

function update_Sk_MvGaussian(x,xbar_k,Ntk,rtik)
    K = length(xbar_k)
    T = length(x)
    C_t = [length(el) for el in x]
    Nk = sum(Ntk)[1:K]
    Sk_vec = Vector{VecOrMat{Float64}}(undef,K)
    for k in 1:K
        Sk = []
        for t in 1:T
            for i in 1:C_t[t]
                diff_vect = x[t][i] - xbar_k[k]
                push!(Sk, rtik[t][i][k] .*   diff_vect   * transpose(diff_vect))
            end
        end
        Nk_scale = 1 ./ Nk[k]
        Sk_vec[k] = Nk_scale .* sum(Sk)
    end
    return Sk_vec
end

function update_Wk_MvGaussian(xbar_k,mu0,m0,Ntk,Sk_vec,W0)
    K = length(xbar_k)
    Wk_vec = Vector{VecOrMat{Float64}}(undef,K)
    Nk = sum(Ntk)[1:K]
    W0_inv = inv(W0)
    for k in 1:K
        N_k = Nk[k]
        S_k = Sk_vec[k]
        xbar_ = xbar_k[k]
        diff_vect = xbar_ .- mu0
        m0_to_N_k_ratio =  m0*N_k/(N_k + m0)
        # @debug "m0_to_N_k_ratio for state $k: $m0_to_N_k_ratio"
        Wk =  W0_inv + N_k .* S_k + m0_to_N_k_ratio .* diff_vect   * transpose(diff_vect)
        # @debug "S_k for state $k: $(S_k)"
        Wk_vec[k] = inv(Wk)
    end
    return Wk_vec
end


function logUk_expected_value_MvGaussian(rho_hat,omega_hat)
    return digamma.(rho_hat .* omega_hat) .- digamma.(omega_hat)
end
function log1minusUk_expected_value_MvGaussian(rho_hat,omega_hat)
    return digamma.((1.0 .- rho_hat) .* omega_hat) .- digamma.(omega_hat)
end
function calc_Tk(e_log_π)
    K = length(sum(e_log_π)) -1
    return sum(e_log_π)[1:K]
end
function LowerBound_LG_MvGaussian(rho_hat,omega_hat,G,γ,α_0,T_k)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value_MvGaussian(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value_MvGaussian(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value_MvGaussian(γ,K)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    # lb_lg_k = [-cB[k] + (G + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (G*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    lb_lg_k = -1.0 .* c_B .+  (G .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (G .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α_0 .* e_βk .* T_k
    #take Negative to find max  
    lb_lg = -sum(lb_lg_k)
    return lb_lg
end
function LowerBound_LG_unconstrained_MvGaussian(c,d,G,γ,α_0,T_k)
    rho_hat = sigmoid.(c)
    omega_hat = exp.(d)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value_MvGaussian(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value_MvGaussian(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value_MvGaussian(γ,K)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    # lb_lg_k = [-cB[k] + (G + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (G*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    lb_lg_k = -1.0 .* c_B .+  (G .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (G .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α_0 .* e_βk .* T_k
    #take Negative to find max  
    lb_lg = -sum(lb_lg_k)
    return lb_lg
end
# function LowerBound_LG_test(x,G,γ,α_0,T_k)
#     #rho_hat ->  x[1,:]
#     #omega_hat -> x[2,:]
#     c_B = beta.(x[1,:] .* x[2,:] , (1.0 .- x[1,:]) .* x[2,:])
#     e_logUk = logUk_expected_value(x[1,:],x[2,:])
#     e_log1minusUk =  log1minusUk_expected_value(x[1,:],x[2,:])
#     K = length(x[1,:])
#     e_βk = βk_expected_value(γ,K)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
#     k_vec = collect(1:K)
#     lb_lg_k = -1.0 .* c_B .+  (G .+ 1. .- x[1,:] .* x[2,:]) .*e_logUk  .+  (G .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- x[1,:]) .* x[2,:]) .* e_log1minusUk  .+  α_0 .* e_βk .* T_k
#     lb_lg = sum(lb_lg_k)
#     return lb_lg
# end
function calc_Hz(rtik)
    z_entropy = 0.0
    T = length(rtik)
    C_t = [length(el) for el in rtik]
    K = length(rtik[1][1])
    # @debug K, C_t, T
    
    # for k in 1:K
    #     for t in 1:T
    #         for i in 1:C_t[t]
    #             # if isnan(rtik[t][i][k]*log(rtik[t][i][k]))
    #             #     @debug rtik[t][i][k], log(rtik[t][i][k]), t, i, k
    #             # end
    #             log_rtik = log(rtik[t][i][k])
    #             if isinf(log_rtik)
    #                 log_rtik = 0.0
    #             end
    #             z_entropy += rtik[t][i][k]*log_rtik
    #         end
    #     end
    # end
    
    for t in 1:T
        for i in 1:C_t[t]
                # if isnan(rtik[t][i][k]*log(rtik[t][i][k]))
                #     @debug rtik[t][i][k], log(rtik[t][i][k]), t, i, k
                # end
            z_entropy += entropy(rtik[t][i][1:K])
        end
    end
    return  -z_entropy
end
# function LowerBound_LG_k(rho_hat,omega_hat,G,γ,α_0,T_k,k)
#     cB = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
#     e_logUk = logUk_expected_value(rho_hat,omega_hat)
#     e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
#     K = length(rho_hat)
#     e_βk = βk_expected_value(γ,K) # e_βk = βk_expected_value(γ,K)[1:K]
#     k_vec = collect(1:K)[k]
#     # lb_lg_k = [-cB[k] + (G + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (G*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
#     lb_lg_k = -cB[k] .+  (G .+ 1. .- rho_hat[k] .* omega_hat[k]) .*e_logUk[k]  +  (G .* ( K .+ 1. .- k_vec) .+ γ - (1.0 .- rho_hat[k]) .* omega_hat[k]) .* e_log1minusUk[k]  .+  α_0 .* e_βk[k] .* T_k[k]
#     lb_lg = sum(lb_lg_k)
#     return lb_lg
# end
# function LowerBound_LG2_k(vec_args,G,γ,α_0,T_k,k)
#     rho_hat = vec_args[1]
#     omega_hat = vec_args[2]
#     LowerBound_LG_k(rho_hat,omega_hat,G,γ,α_0,T_k,k)
# end
# LowerBound_LG2_k_closure(G,γ,α_0,T_k,k) = vec_args ->  LowerBound_LG2_k(vec_args,G,γ,α_0,T_k,k)

function LowerBound_LG2_MvGaussian(vec_args,G,γ,α_0,T_k)
    rho_hat = vec_args[1, :]
    omega_hat = vec_args[2, :]
    LowerBound_LG_MvGaussian(rho_hat,omega_hat,G,γ,α_0,T_k)
end
function LowerBound_LG2_unconstrained_MvGaussian(vec_args,G,γ,α_0,T_k)
    c = vec_args[1, :] # rho_hat
    d = vec_args[2, :] # omega_hat
    LowerBound_LG_unconstrained_MvGaussian(c,d,G,γ,α_0,T_k)
end

function genterate_Delta_mk_MvGaussian(rho_hat, γ)
    K = length(rho_hat)
    e_βk = βk_expected_value_MvGaussian(γ,K)
    Δ_mk =  Matrix{Float64}(undef,K,K+1)
    for k in 1:K+1
        for m in 1:K
            if m < k
                Δ_mk[m,k] = -1.0 /(1 - rho_hat[m]) * e_βk[k]
            elseif m == k
                Δ_mk[m,k] = 1.0 / rho_hat[m] * e_βk[k]
            elseif m > 0
                Δ_mk[m,k] = 0.0
            end
        end
    end
    return Δ_mk
end

function g_constrained_MvGaussian!(GD, x, G,γ,α_0,T_k)
    rho_hat =  x[1,:]
    omega_hat = x[2,:]

    K = length(x[1,:])
    k_vec = collect(1:K)
    rho_omega_hat =  rho_hat.*omega_hat 
    minusrho_omega_hat =  (1 .- rho_hat) .* omega_hat 
    minusrho_hat = (1 .- rho_hat)
    Δ_mk = genterate_Delta_mk_MvGaussian(rho_hat, γ)
    Δ_prod_T =  Δ_mk[:,1:K] .* T_k
    a_Δ_prod_T = α_0 .* permutedims(sum(Δ_prod_T,dims=1))
    #take Negative to find max  

    GD[1,:] = -1.0 .* ( omega_hat .* (G .+ 1 .-  rho_omega_hat) .* polygamma.(1,rho_omega_hat) .- omega_hat .*( G .*(K .+ 1.0 .- k_vec)  .+ γ .- minusrho_omega_hat) .* polygamma.(1,minusrho_omega_hat) .+ a_Δ_prod_T) #  rho
    #-2.0 .* (1.0 .- x[1,:]) .- 400.0 .* (x[2,:] .- x[1,:].^2) .* x[1,:]
    GD[2,:] = -1.0 .* ((G .+ 1.0 .- rho_omega_hat) .* (rho_hat .* polygamma.(1,rho_omega_hat) .- polygamma.(1,omega_hat)) .+ ( G .* (K .+ 1.0 .- k_vec) .+ γ .-  minusrho_omega_hat) .* ( minusrho_hat .* polygamma.(1, minusrho_omega_hat)  .- polygamma.(1,omega_hat) )) #  omega
     #(x[2,:] .- x[1,:].^2)
end
function g_unconstrained_MvGaussian!(GD, x, G,γ,α_0,T_k)
    rho_hat =  x[1,:] # c
    omega_hat = x[2,:] # d

    K = length(x[1,:])
    k_vec = collect(1:K)
    rho_omega_hat =  rho_hat.*omega_hat 
    minusrho_omega_hat =  (1 .- rho_hat) .* omega_hat 
    minusrho_hat = (1 .- rho_hat)
    Δ_mk = genterate_Delta_mk_MvGaussian(rho_hat, γ)
    Δ_prod_T =  Δ_mk[:,1:K] .* T_k
    a_Δ_prod_T = α_0 .* permutedims(sum(Δ_prod_T,dims=1))
    #take Negative to find max  

    GD[1,:] = -1.0 .*(rho_hat .* minusrho_hat .* (omega_hat .* (G .+ 1 .-  rho_omega_hat) .* polygamma.(1,rho_omega_hat) .- omega_hat .*( G .*(K .+ 1.0 .- k_vec)  .+ γ .- minusrho_omega_hat) .* polygamma.(1,minusrho_omega_hat) .+ a_Δ_prod_T)) # c and rho
    #-2.0 .* (1.0 .- x[1,:]) .- 400.0 .* (x[2,:] .- x[1,:].^2) .* x[1,:] 
    GD[2,:] = -1.0 .* ( omega_hat .* ((G .+ 1.0 .- rho_omega_hat) .* (rho_hat .* polygamma.(1,rho_omega_hat) .- polygamma.(1,omega_hat)) .+ ( G .* (K .+ 1.0 .- k_vec) .+ γ .-  minusrho_omega_hat) .* ( minusrho_hat .* polygamma.(1, minusrho_omega_hat)  .- polygamma.(1,omega_hat) ))) # d and omega
     #(x[2,:] .- x[1,:].^2)
end
# g_constrained!(Matrix{Float64}(undef,2,K), vec_args0, G,γ,α_0,T_k)
# g_unconstrained!(Matrix{Float64}(undef,2,K), cd_vec_args0, G,γ,α_0,T_k)
g_constrained_closure_MvGaussian!(G,γ,α_0,T_k) = (GD, x)  ->  g_constrained_MvGaussian!(GD, x, G,γ,α_0,T_k)
g_unconstrained_closure_MvGaussian!(G,γ,α_0,T_k) = (GD, x)  ->  g_unconstrained_MvGaussian!(GD, x, G,γ,α_0,T_k)




LowerBound_LG2_closure_MvGaussian(G,γ,α_0,T_k) = vec_args ->  LowerBound_LG2_MvGaussian(vec_args,G,γ,α_0,T_k)
LowerBound_LG2__unconstrained_closure_MvGaussian(G,γ,α_0,T_k) = vec_args ->  LowerBound_LG2_unconstrained_MvGaussian(vec_args,G,γ,α_0,T_k)
LowerBound_LG2_test_closure_MvGaussian(G,γ,α_0,T_k) = vec_args ->  LowerBound_LG_test_MvGaussian(vec_args,G,γ,α_0,T_k)

