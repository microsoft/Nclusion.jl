function generate_a0k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    a0k_hat_vec = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G,K))
    else
        values = ones(float_type,G,K)
    end
    
    for k in 1:K
        @views a0k_hat_vec[k] = values[:,k]
    end
    return a0k_hat_vec
end
function generate_b0k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    b0k_hat_vec = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = covert.(float_type,rand(Uniform(0,1),G,K))
    else
        values = ones(float_type,G,K)
    end
    
    for k in 1:K
        @views b0k_hat_vec[k] = values[:,k]
    end
    return b0k_hat_vec
end
function generate_v_sq_k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    v_sq_k_hat_vec = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G,K))
    else
        values = ones(float_type,G,K)
    end
    
    for k in 1:K
        @views v_sq_k_hat_vec[k] = values[:,k]
    end
    return v_sq_k_hat_vec
end
function generate_σ_sq_k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    σ_sq_k_hat_vec = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G,K))
    else
        values = ones(float_type,G,K)
    end
    
    for k in 1:K
        @views σ_sq_k_hat_vec[k] = values[:,k]
    end
    return σ_sq_k_hat_vec
end
function generate_λ0k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    λ0k_hat_vec = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G,K))
    else
        values = ones(float_type,G,K)
    end
    
    for k in 1:K
        @views λ0k_hat_vec[k] = values[:,k]
    end
    return λ0k_hat_vec
end

function generate_λj_testing(;G=20000,rand_init=false, float_type=Float64)
    λj_hat_vec = Vector{float_type}(undef,G)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G))
    else
        values = ones(float_type,G)
    end
    
    for j in 1:G
        @views λj_hat_vec[j] = values[j]
    end
    return λj_hat_vec
end

function generate_mk_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    m0k_hat_vec = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G,K))
    else
        values = ones(float_type,G,K)
    end
    
    for k in 1:K
        @views m0k_hat_vec[k] = values[:,k]
    end
    return m0k_hat_vec
end
function generate_a0_err_testing(;G=20000,rand_init=false, float_type=Float64)
    a0_err_hat_vec = Vector{float_type}(undef,G)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G))
    else
        values = ones(float_type,G)
    end
    
    @views a0_err_hat_vec = values
    return a0_err_hat_vec
end
function generate_b0_err_testing(;G=20000,rand_init=false, float_type=Float64)
    b0_err_hat_vec = Vector{float_type}(undef,G)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G))
    else
        values = ones(float_type,G)
    end
    
    @views b0_err_hat_vec = values
    return b0_err_hat_vec
end
function generate_λ0_err_testing(;G=20000,rand_init=false, float_type=Float64)
    λ0_err_hat_vec = Vector{float_type}(undef,G)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G))
    else
        values = ones(float_type,G)
    end
    
    @views λ0_err_hat_vec = values
    return λ0_err_hat_vec
end
function generate_m_err_testing(;G=20000,rand_init=false, float_type=Float64)
    m_err_hat_vec = Vector{float_type}(undef,G)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),G))
    else
        values = zeros(float_type,G)
    end
    
    @views m_err_hat_vec = values
    return m_err_hat_vec
end
function generate_rhok_testing(;K=50,rand_init=false, float_type=Float64)
    rhok_hat_vec = Vector{float_type}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),K))
    else
        values = 0.5 .* ones(float_type,K)
    end
    @views rhok_hat_vec = values
    return rhok_hat_vec
end
function generate_gk_testing(;K=50,rand_init=false, float_type=Float64)
    gk_hat_vec = Vector{float_type}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),K))
    else
        values = 0.5 .* ones(float_type,K)
    end
    @views gk_hat_vec = values
    return gk_hat_vec
end
function generate_omegak_testing(;K=50,rand_init=false, float_type=Float64)
    omegak_hat_vec = Vector{float_type}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,2),K))
    else
        values = 2.0 .* ones(float_type,K)
    end
    @views omegak_hat_vec = values
    return omegak_hat_vec
end
function generate_hk_testing(;K=50,rand_init=false, float_type=Float64)
    hk_hat_vec = Vector{float_type}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,2),K))
    else
        values = 2.0 .* ones(float_type,K)
    end
    @views hk_hat_vec = values
    return hk_hat_vec
end
function generate_a_γ_testing(;rand_init=false, float_type=Float64)
    a_γ_hat = 0.0
    if rand_init
        a_γ_hat = convert(float_type,rand(Uniform(0,10)))
    else
        a_γ_hat = convert(float_type,1.0)
    end
    return a_γ_hat
end
function generate_b_γ_testing(;rand_init=false, float_type=Float64)
    b_γ_hat = 0.0
    if rand_init
        b_γ_hat = convert(float_type,rand(Uniform(0,10)))
    else
        b_γ_hat = convert(float_type,1.0) 
    end
    return b_γ_hat
end
function generate_awt_testing(;T=5,rand_init=false, float_type=Float64)
    awt_hat_vec = Vector{float_type}(undef,T)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),T))
    else
        values = ones(float_type,T)
    end
    @views awt_hat_vec = values
    return awt_hat_vec
end
function generate_bwt_testing(;T=5,rand_init=false, float_type=Float64)
    bwt_hat_vec = Vector{float_type}(undef,T)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,1),T))
    else
        values = ones(float_type,T)
    end
    @views bwt_hat_vec = values
    return bwt_hat_vec
end
function generate_a_αt_testing(;T=5,rand_init=false, float_type=Float64)
    a_αt_hat_vec = Vector{float_type}(undef,T)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,10),T))
    else
        values = ones(float_type,T)
    end
    @views a_αt_hat_vec = values
    return a_αt_hat_vec
end
function generate_b_αt_testing(;T=5,rand_init=false, float_type=Float64)
    b_αt_hat_vec = Vector{float_type}(undef,T)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,10),T))
    else
        values = ones(float_type,T)
    end
    @views b_αt_hat_vec = values
    return b_αt_hat_vec
end
function generate_st_testing(;T=5,rand_init=false, float_type=Float64)
    st_hat_vec = Vector{float_type}(undef,T)
    if rand_init
        values = convert.(float_type,rand(Uniform(0,10),T))
    else
        values = ones(float_type,T)
    end
    @views st_hat_vec = values
    return st_hat_vec
end
function generate_c_ttprime_testing(;T=5,rand_init=false, float_type=Float64)
    c_ttprime_vec = Vector{Vector{float_type}}(undef,T)
    if rand_init
        values = convert.(float_type,rand(Dirichlet(ones(T) ./T),T))
    else
        values = 1/T ./ ones(float_type,T,T)
    end

    for t in 1:T
        @views c_ttprime_vec[t] = values[:,t]
    end
    return c_ttprime_vec
end
function generate_θ_hat_testing(;T=5,K=50,rand_init=false, float_type=Float64)
    θ_hat_vec = Vector{Vector{float_type}}(undef,T)
    Kplus = K+1 
    if rand_init
        values = convert.(float_type,rand(Uniform(0,10),Kplus,T))
    else
        values = ones(float_type,Kplus,T)
    end
    for t in 1:T
        @views θ_hat_vec[t] = values[:,t]
    end
    return θ_hat_vec
end
function generate_d_hat_testing(;T=5,K=50,rand_init=false, float_type=Float64)
    d_hat_vec = Vector{Vector{float_type}}(undef,T)
    Kplus = K+1 
    if rand_init
        values = convert.(float_type,rand(Uniform(0,10),Kplus,T))
    else
        values = ones(float_type,Kplus,T)
    end
    for t in 1:T
        @views d_hat_vec[t] = values[:,t]
    end
    return d_hat_vec
end

function generate_rtik_testing(;T=5,K=50,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    if rand_init
        rtik = [[convert.(float_type,rand(Dirichlet(ones(K) ./K))) for i in 1:C_t[t]] for t in 1:T]
    else
        rtik = [[ones(float_type,K) ./K for i in 1:C_t[t]] for t in 1:T]
    end
    return rtik
end

function generate_v_tikj_testing(;G=20000,K=50,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    v_tikj = Vector{Vector{Vector{Vector{Vector{float_type}}}}}(undef,T)
    if rand_init
        for t in 1:T
            v_tikj[t]=Vector{Vector{Vector{Vector{float_type}}}}(undef,C_t[t])
            for i in 1:C_t[t] 
                v_tikj[t][i] = Vector{Vector{Vector{float_type}}}(undef,K)
                for k in 1:K
                    v_tikj[t][i][k] = Vector{Vector{float_type}}(undef,G)
                    for j in 1:G
                        v_tikj[t][i][k][j] = Vector{float_type}(undef,2)
                        value = convert.(float_type,rand(Dirichlet(ones(2) ./2)))
                        v_tikj[t][i][k][j][1] = value[1]
                        v_tikj[t][i][k][j][2] = value[2]
                    end
                end
            end
        end
        # v_tikj = [[[[rand(Dirichlet(ones(2) ./2))  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    else
        # v_tikj = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
        for t in 1:T
            v_tikj[t]=Vector{Vector{Vector{Vector{float_type}}}}(undef,C_t[t])
            for i in 1:C_t[t] 
                v_tikj[t][i] = Vector{Vector{Vector{float_type}}}(undef,K)
                for k in 1:K
                    v_tikj[t][i][k] = Vector{Vector{float_type}}(undef,G)
                    for j in 1:G
                        v_tikj[t][i][k][j] = Vector{float_type}(undef,2) 
                        value =  ones(float_type,2) ./2
                        v_tikj[t][i][k][j][1] = value[1]
                        v_tikj[t][i][k][j][2] = value[2]
                    end
                end
            end
        end
    end
    return v_tikj
end

function generate_pip_kj_imp_weights_testing(;G=20000,K=50,rand_init=false, float_type=Float64)

    pip_kj = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Dirichlet(ones(G) ./G),K))
    else
        values = 1/G ./ ones(float_type,G,K)
    end

    for k in 1:K
        @views pip_kj[k] = values[:,k]
    end
    return pip_kj
end

function generate_y_kj_imp_weights_testing(;G=20000,K=50,rand_init=false, float_type=Float64)

    y_kj = Vector{Vector{float_type}}(undef,K)
    if rand_init
        values = convert.(float_type,rand(Dirichlet(ones(G) ./G),K))
    else
        values = 1/G ./ ones(float_type,G,K)
    end

    for k in 1:K
        @views y_kj[k] = values[:,k]
    end
    return y_kj
end
function generate_η_kj_prior_testing(;G=20000,K=50,pct_important=0.5,rand_init=false, float_type=Float64)
    η_prior = Vector{Vector{float_type}}(undef,K)
    η_prior = [[pct_important for j in 1:G] for k in 1:K]

    for k in 1:K
        @views η_prior[k] =  pct_important .* ones(float_type,G)
    end
    return η_prior
end

function generate_x_testing(;G=20000,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    x = Vector{Vector{Vector{float_type}}}(undef,T)
    if rand_init
        for t in 1:T
            x[t]=Vector{Vector{float_type}}(undef,C_t[t])
            for i in 1:C_t[t] 
                x[t][i] = Vector{float_type}(undef,G)
                x[t][i]=convert.(float_type,rand(G))
            end
        end
    else
        for t in 1:T
            x[t]=Vector{Vector{float_type}}(undef,C_t[t])
            for i in 1:C_t[t] 
                x[t][i] = Vector{float_type}(undef,G)
                x[t][i]= ones(float_type,G)
            end
        end
    end
    return x
end

function generate_e_log_π_testing(;T=5,K=50, float_type=Float64)
    e_log_π = [ones(float_type,K) for t in 1:T]
    return e_log_π
end

function generate_e_log_τ_testing(;K=50, float_type=Float64)
    e_log_τ = ones(float_type,K)
    return e_log_τ
end
function generate_e_log_τkj_testing(;G=20000,K=50, float_type=Float64)
    e_log_τkj = [ones(float_type,G) for k in 1:K]
    return e_log_τkj
end
function generate_e_log_τkj_err_testing(;G=20000,K=50, float_type=Float64)
    e_log_τkj = ones(float_type,G)
    return e_log_τkj
end
function generate_e_τ_μ_testing(;K=50,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    e_τ_μ = Vector{Vector{Vector{float_type}}}(undef,T)
    for t in 1:T
        e_τ_μ[t]=Vector{Vector{float_type}}(undef,C_t[t])
        for i in 1:C_t[t] 
            e_τ_μ[t][i] = Vector{float_type}(undef,K)
            e_τ_μ[t][i]=ones(float_type,K)
        end
    end

    return e_τ_μ
end
function generate_e_τ_μ_tikj_testing(;G=20000,K=50,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    e_τ_μ_tikj = Vector{Vector{Vector{Vector{float_type}}}}(undef,T)
    for t in 1:T
        e_τ_μ_tikj[t]=Vector{Vector{Vector{float_type}}}(undef,C_t[t])
        for i in 1:C_t[t] 
            e_τ_μ_tikj[t][i] = Vector{Vector{float_type}}(undef,K)
            for k in 1:K
                e_τ_μ_tikj[t][i][k] = Vector{float_type}(undef,G)
                e_τ_μ_tikj[t][i][k] = ones(float_type,G)
            end
        end
    end

    return e_τ_μ_tikj
end

function generate_e_τ_0j_err_testing(;G=20000,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    e_τ_0j_err = Vector{Vector{Vector{float_type}}}(undef,T)
    for t in 1:T
        e_τ_0j_err[t]=Vector{Vector{float_type}}(undef,C_t[t])
        for i in 1:C_t[t] 
            e_τ_0j_err[t][i] = Vector{float_type}(undef,G)
            e_τ_0j_err[t][i] = ones(float_type,G)
        end
    end

    return e_τ_0j_err
end
function generate_e_τ_0_err_testing(;T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    e_τ_0_err = Vector{Vector{float_type}}(undef,T)
    for t in 1:T
        e_τ_0_err[t]=Vector{float_type}(undef,C_t[t])
        e_τ_0_err[t]=ones(float_type,C_t[t])
    end

    return e_τ_0_err
end
function generate_Ntk_testing(;T=5,K=50, float_type=Float64)
    Kplus = K+1
    Ntk = [ones(float_type,Kplus) for t in 1:T]
    return Ntk
end

function generate_N_signal_testing(;G=20000,K=50,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    N_signal = [[[ones(float_type,G) for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return N_signal
end
function generate_N_error_testing(;G=20000,K=50,T=5,N=1_000_000,C_t=nothing,rand_init=false, float_type=Float64)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end

    N_error = [[[ones(float_type,G) for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    return N_error
end
function generate_Nj_error_testing(;G=20000,rand_init=false, float_type=Float64)
    Nj_error = ones(float_type,G)
    return Nj_error
end
function generate_Nkj_signal_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    Nkj_signal = [ones(float_type,G) for k in 1:K]
    return Nkj_signal
end
function generate_Nk_testing(;K=50,rand_init=false, float_type=Float64)
    Nk = ones(float_type,K)
    return Nk
end
function generate_x_hat_k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    x_hat_k = [ones(float_type,G) for k in 1:K]
    return x_hat_k
end
function generate_x_hat_sq_k_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    x_hat_sq_k = [ones(float_type,G) for k in 1:K]
    return x_hat_sq_k
end
function generate_x_hat_error_testing(;G=20000,rand_init=false, float_type=Float64)
    x_hat_err = ones(float_type,G)
    return x_hat_err
end
function generate_x_hatk_signal_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    x_hatk_signal = [ones(float_type,G) for k in 1:K]
    return x_hatk_signal
end
function generate_x_hat_sq_error_testing(;G=20000,rand_init=false, float_type=Float64)
    x_hat_sq_err = ones(float_type,G)
    return x_hat_sq_err
end
function generate_x_hatk_sq_signal_testing(;G=20000,K=50,rand_init=false, float_type=Float64)
    x_hatk_sq_signal = [ones(float_type,G) for k in 1:K]
    return x_hatk_sq_signal
end
function generate_e_γ_testing(;rand_init=false, float_type=Float64)
    return convert(float_type,1.0)
end
function generate_Tαk_testing(;K=50,rand_init=false, float_type=Float64)
    Kplus = K+1
    Tαk = ones(float_type,Kplus) 
    return Tαk
end
function generate_Tk_testing(;K=50,rand_init=false, float_type=Float64)
    Kplus = K+1
    Tk = ones(float_type,Kplus) 
    return Tk
end
function generate_hyperparamters_testing(;rand_init=false, float_type=Float64)
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    return convert.(float_type,(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w))
end
function generate_init_vec_testing(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=20000,K=50,T=5,N=1_000_000,C_t=nothing, nothing_init=false)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0);

    if nothing_init
        λ0k_hat_vec_init = nothing; # 
        mk_hat_vec_init = nothing;
        a0k_hat_vec_init = nothing; #
        b0k_hat_vec_init =  nothing; #
        rhok_hat_vec_init, omegak_hat_vec_init = nothing,nothing;
        a_γ_hat_init = nothing;
        b_γ_hat_init = nothing;
        awt_hat_vec_init = nothing;
        bwt_hat_vec_init =  nothing;
        a_αt_hat_vec_init = nothing;
        b_αt_hat_vec_init = nothing;
        c_ttprime_vec_init =  nothing;
        θ_hat_vec_init = nothing;
        rtik_init = nothing;
    else
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]; # 
        mk_hat_vec_init = [μ0_vec for k in 1:K];
        a0k_hat_vec_init = [a0_vec for k in 1:K]; #
        b0k_hat_vec_init =  [b0_vec for k in 1:K]; #
        rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K);
        a_γ_hat_init = a_γ;
        b_γ_hat_init = b_γ;
        awt_hat_vec_init = [adot_w for t in 1:T];
        bwt_hat_vec_init =  [bdot_w for t in 1:T];
        a_αt_hat_vec_init = [a_α for t in 1:T];
        b_αt_hat_vec_init = [b_α for t in 1:T];
        c_ttprime_vec_init =  [ones(T) ./T  for t in 1:T];
        θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T];
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T];
    end

    return λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init, omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,rtik_init
end
function generate_hyperparamters_vs12_testing(;rand_init=false, float_type=Float64)
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    return convert.(float_type(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err))
end
function generate_init_vec_vs12_testing(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err;G=20000,K=50,T=5,N=1_000_000,C_t=nothing, nothing_init=false)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0);
    λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);
    if nothing_init
        λ0k_hat_vec_init = nothing; # 
        mk_hat_vec_init = nothing;
        a0k_hat_vec_init = nothing; #
        b0k_hat_vec_init =  nothing; #
        rhok_hat_vec_init, omegak_hat_vec_init = nothing,nothing;
        a_γ_hat_init = nothing;
        b_γ_hat_init = nothing;
        awt_hat_vec_init = nothing;
        bwt_hat_vec_init =  nothing;
        a_αt_hat_vec_init = nothing;
        b_αt_hat_vec_init = nothing;
        c_ttprime_vec_init =  nothing;
        θ_hat_vec_init = nothing;
        rtik_init = nothing;
        v_tikj_vec_init = nothing;
        λ0_err_hat_vec_init = nothing;
        m_err_hat_vec_init =nothing;
        a0_err_hat_vec_init = nothing;
        b0_err_hat_vec_init = nothing;
    else
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]; # 
        mk_hat_vec_init = [μ0_vec for k in 1:K];
        a0k_hat_vec_init = [a0_vec for k in 1:K]; #
        b0k_hat_vec_init =  [b0_vec for k in 1:K]; #
        rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K);
        a_γ_hat_init = a_γ;
        b_γ_hat_init = b_γ;
        awt_hat_vec_init = [adot_w for t in 1:T];
        bwt_hat_vec_init =  [bdot_w for t in 1:T];
        a_αt_hat_vec_init = [a_α for t in 1:T];
        b_αt_hat_vec_init = [b_α for t in 1:T];
        c_ttprime_vec_init =  [ones(T) ./T  for t in 1:T];
        θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T];
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T];
        v_tikj_vec_init = [[[[ones(2) ./2  for j in 1:G] for k in 1:K] for i in 1:C_t[t]] for t in 1:T];
        λ0_err_hat_vec_init = λ0_err_vec;
        m_err_hat_vec_init =μ0_err_vec;
        a0_err_hat_vec_init = a0_err_vec;
        b0_err_hat_vec_init = b0_err_vec;
    end

    return λ0_vec, μ0_vec, a0_vec, b0_vec,λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init, omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,rtik_init,v_tikj_vec_init,λ0_err_hat_vec_init,m_err_hat_vec_init,a0_err_hat_vec_init,b0_err_hat_vec_init
end
function generate_hyperparamters_testing18(;rand_init=false)
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    return convert.(float_type,(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w))
end
function generate_init_vec_testing18(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=20000,K=50,T=5,N=1_000_000,C_t=nothing, nothing_init=false)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0);

    if nothing_init
        λ0k_hat_vec_init = nothing; # 
        mk_hat_vec_init = nothing;
        a0k_hat_vec_init = nothing; #
        b0k_hat_vec_init =  nothing; #
        rhok_hat_vec_init, omegak_hat_vec_init = nothing,nothing;
        a_γ_hat_init = nothing;
        b_γ_hat_init = nothing;
        awt_hat_vec_init = nothing;
        bwt_hat_vec_init =  nothing;
        a_αt_hat_vec_init = nothing;
        b_αt_hat_vec_init = nothing;
        c_ttprime_vec_init =  nothing;
        θ_hat_vec_init = nothing;
        rtik_init = nothing;
        pip_kj_init = nothing;
    else
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]; # 
        mk_hat_vec_init = [μ0_vec for k in 1:K];
        a0k_hat_vec_init = [a0_vec for k in 1:K]; #
        b0k_hat_vec_init =  [b0_vec for k in 1:K]; #
        rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K);
        a_γ_hat_init = a_γ;
        b_γ_hat_init = b_γ;
        awt_hat_vec_init = [adot_w for t in 1:T];
        bwt_hat_vec_init =  [bdot_w for t in 1:T];
        a_αt_hat_vec_init = [a_α for t in 1:T];
        b_αt_hat_vec_init = [b_α for t in 1:T];
        c_ttprime_vec_init =  [ones(T) ./T  for t in 1:T];
        θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T];
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T];
        pip_kj_init =  [ones(G) ./G  for j in 1:G]
    end

    return λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init, omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,rtik_init,pip_kj_init
end
function generate_hyperparamters_testing25(;rand_init=false, float_type=Float64)
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    return convert.(float_type,(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w))
end
function generate_hyperparamters_testing(;rand_init=false, float_type=Float64)
    ηk,α0,γ0,ϕ0 = 1.0, 1.0, 1.0, 1.0
    return convert.(float_type,(ηk,α0,γ0,ϕ0 ))
end
function generate_init_vec_testing25(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=20000,K=50,T=5,N=1_000_000,C_t=nothing, nothing_init=false)
    if isnothing(C_t)
        C_t = Int.(N/T .* ones(Int,T))
    else
        N=sum(C_t)
    end
    λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0);

    if nothing_init
        λ0k_hat_vec_init = nothing; # 
        mk_hat_vec_init = nothing;
        a0k_hat_vec_init = nothing; #
        b0k_hat_vec_init =  nothing; #
        rhok_hat_vec_init, omegak_hat_vec_init = nothing,nothing;
        a_γ_hat_init = nothing;
        b_γ_hat_init = nothing;
        awt_hat_vec_init = nothing;
        bwt_hat_vec_init =  nothing;
        a_αt_hat_vec_init = nothing;
        b_αt_hat_vec_init = nothing;
        c_ttprime_vec_init =  nothing;
        θ_hat_vec_init = nothing;
        rtik_init = nothing;
        pip_kj_init = nothing;
    else
        λ0k_hat_vec_init = [λ0_vec for k in 1:K]; # 
        mk_hat_vec_init = [μ0_vec for k in 1:K];
        a0k_hat_vec_init = [a0_vec for k in 1:K]; #
        b0k_hat_vec_init =  [b0_vec for k in 1:K]; #
        rhok_hat_vec_init, omegak_hat_vec_init = init_params_states(K);
        a_γ_hat_init = a_γ;
        b_γ_hat_init = b_γ;
        awt_hat_vec_init = [adot_w for t in 1:T];
        bwt_hat_vec_init =  [bdot_w for t in 1:T];
        a_αt_hat_vec_init = [a_α for t in 1:T];
        b_αt_hat_vec_init = [b_α for t in 1:T];
        c_ttprime_vec_init =  [ones(T) ./T  for t in 1:T];
        θ_hat_vec_init = [ones(K+1) ./(K+1)  for t in 1:T];
        rtik_init = [[ones(K) ./K for i in 1:C_t[t]] for t in 1:T];
        pip_kj_init =  [ones(G) ./G  for j in 1:G]
    end

    return λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init, omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,rtik_init,pip_kj_init
end
function generate_Glog_testing(;G=20000,rand_init=false, float_type=Float64)
    Glog = G*log(2π)
    return convert(float_type,Glog)
end
function generate_logpi_testing(;rand_init=false, float_type=Float64)
    logpi = log(2π)
    return convert(float_type,logpi)
end


#λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
function test_generation_functions(;G=20000,K=50,N=1_000_000,T=5, float_type=Float64)
    C_t = Int.(N/T .* ones(Int,T))

    rtik = generate_rtik_testing(;T=T,N=N,C_t=C_t,rand_init = false, float_type=float_type);
    r_rtik = generate_rtik_testing(;T=T,N=N,C_t=C_t,rand_init = true, float_type=float_type);
    @test length(rtik) == T
    @test all(length.(rtik) .== C_t)
    @test all(isprobvec.(vcat(rtik...)))
    @test length(r_rtik) == T
    @test all(length.(r_rtik) .== C_t)
    @test all(isprobvec.(vcat(r_rtik...)))
    
    a0_err_hat_vec = generate_a0_err_testing(;G=G,rand_init=false, float_type=float_type);
    r_a0_err_hat_vec = generate_a0_err_testing(;G=G,rand_init=true, float_type=float_type);
    @test length(a0_err_hat_vec) == G
    @test length(r_a0_err_hat_vec) == G

    b0_err_hat_vec = generate_b0_err_testing(;G=G,rand_init=false, float_type=float_type);
    r_b0_err_hat_vec = generate_b0_err_testing(;G=G,rand_init=true, float_type=float_type);
    @test length(b0_err_hat_vec) == G
    @test length(r_b0_err_hat_vec) == G

    λ0_err_hat_vec = generate_λ0_err_testing(;G=G,rand_init=false, float_type=float_type);
    r_λ0_err_hat_vec = generate_λ0_err_testing(;G=G,rand_init=true, float_type=float_type);
    @test length(λ0_err_hat_vec) == G
    @test length(r_λ0_err_hat_vec) == G

    m_err_hat_vec = generate_m_err_testing(;G=G,rand_init=false, float_type=float_type);
    r_m_err_hat_vec = generate_m_err_testing(;G=G,rand_init=true, float_type=float_type);
    @test length(m_err_hat_vec) == G
    @test length(r_m_err_hat_vec) == G

    Kplus = K + 1
    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=false, float_type=float_type);
    r_θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=false, float_type=float_type);
    @test length(θ_hat_vec) == T
    @test all(length.(θ_hat_vec) .== Kplus)
    @test length(r_θ_hat_vec) == T
    @test all(length.(r_θ_hat_vec) .== Kplus)
  


    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=false, float_type=float_type);
    r_c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=true, float_type=float_type);
    @test length(c_ttprime_vec) == T
    @test all(length.(c_ttprime_vec) .== T)
    @test all(isprobvec.(c_ttprime_vec))
    @test length(r_c_ttprime_vec) == T
    @test all(length.(r_c_ttprime_vec) .== T)
    @test all(isprobvec.(r_c_ttprime_vec))
    

    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=false, float_type=float_type);
    r_a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=true, float_type=float_type);
    @test length(a_αt_hat_vec) == T
    @test length(r_a_αt_hat_vec) == T

    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=false, float_type=float_type);
    r_b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=true, float_type=float_type);
    @test length(b_αt_hat_vec) == T
    @test length(r_b_αt_hat_vec) == T

    awt_hat_vec = generate_awt_testing(;T=T,rand_init=false, float_type=float_type);
    r_awt_hat_vec = generate_awt_testing(;T=T,rand_init=true, float_type=float_type);
    @test length(awt_hat_vec) == T
    @test length(r_awt_hat_vec) == T

    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=false, float_type=float_type);
    r_bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=true, float_type=float_type);
    @test length(bwt_hat_vec) == T
    @test length(r_bwt_hat_vec) == T


    a_γ_hat = generate_a_γ_testing(;rand_init=false, float_type=float_type);
    r_a_γ_hat = generate_a_γ_testing(;rand_init=true, float_type=float_type);
    @test typeof(a_γ_hat) <: float_type
    @test typeof(r_a_γ_hat) <: float_type

    b_γ_hat = generate_b_γ_testing(;rand_init=false, float_type=float_type);
    r_b_γ_hat = generate_b_γ_testing(;rand_init=true, float_type=float_type);
    @test typeof(b_γ_hat) <: float_type
    @test typeof(r_b_γ_hat) <: float_type

    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=false, float_type=float_type);
    r_omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=true, float_type=float_type);
    @test length(omegak_hat_vec) == K
    @test length(r_omegak_hat_vec) == K

    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=false, float_type=float_type);
    r_rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=true, float_type=float_type);
    @test length(rhok_hat_vec) == K
    @test length(r_rhok_hat_vec) == K


    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=false, float_type=float_type);
    r_mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=true, float_type=float_type);
    @test length(mk_hat_vec) == K
    @test all([length(mk_hat_vec[k]) == G for k in 1:K])
    @test all([all(mk_hat_vec[k] .== mk_hat_vec[1]) for k in 1:K])
    @test length(r_mk_hat_vec) == K
    @test all([length(r_mk_hat_vec[k]) == G for k in 1:K])

    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=false, float_type=float_type);
    r_λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=true, float_type=float_type);
    @test length(λ0k_hat_vec) == K
    @test all([length(λ0k_hat_vec[k]) == G for k in 1:K])
    @test all([all(λ0k_hat_vec[k] .== λ0k_hat_vec[1]) for k in 1:K])
    @test length(r_λ0k_hat_vec) == K
    @test all([length(r_λ0k_hat_vec[k]) == G for k in 1:K])

    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=false, float_type=float_type);
    r_b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=true, float_type=float_type);
    @test length(b0k_hat_vec) == K
    @test all([length(b0k_hat_vec[k]) == G for k in 1:K])
    @test all([all(b0k_hat_vec[k] .== b0k_hat_vec[1]) for k in 1:K])
    @test length(r_b0k_hat_vec) == K
    @test all([length(r_b0k_hat_vec[k]) == G for k in 1:K])

    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=false, float_type=float_type);
    r_a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=true, float_type=float_type);
    @test length(a0k_hat_vec) == K
    @test all([length(a0k_hat_vec[k]) == G for k in 1:K])
    @test all([all(a0k_hat_vec[k] .== a0k_hat_vec[1]) for k in 1:K])
    @test length(r_a0k_hat_vec) == K
    @test all([length(r_a0k_hat_vec[k]) == G for k in 1:K])
end
# using Random, Distributions, Test
# K = 50;T=5;N=1_000_000;G=20_000;C_t = Int.(N/T .* ones(Int,T));
# K = 6;T=2;N=10;G=4;C_t = Int.(N/T .* ones(Int,T));

# test_generation_functions()
#
# e_log_π = generate_e_log_π_testing(;T=T,K=K)
# e_log_τ = generate_e_log_τ_testing(;K=K)
# e_log_τkj = generate_e_log_τkj_testing(;G=G,K=K)
# e_log_τkj_err = generate_e_log_τkj_err_testing(;G=G,K=K)
# e_τ_μ = generate_e_τ_μ_testing(;K=K,T=T,N=N,C_t=nothing,rand_init=false)
# e_τ_μ_tikj = generate_e_τ_μ_tikj_testing(;G=G,K=K,T=T,N=N,C_t=nothing,rand_init=false)
# e_τ_0j_err = generate_e_τ_0j_err_testing(;G=G,T=T,N=N,C_t=nothing,rand_init=false)
# e_τ_0_err = generate_e_τ_0_err_testing(;T=T,N=N,C_t=nothing,rand_init=false)
# Ntk = generate_Ntk_testing(;T=T,K=K)
# N_signal = generate_N_signal_testing(;G=G,K=K,T=T,N=N,C_t=nothing,rand_init=false)
# N_error = generate_N_error_testing(;G=G,K=K,T=T, N=N ,C_t=nothing,rand_init=false)
# Nj_error = generate_Nj_error_testing(;G=G,rand_init=false)
# Nkj_signal = generate_Nkj_signal_testing(;G=G,K=K,rand_init=false)
# Nk = generate_Nk_testing(;K=K,rand_init=false)
# x_hat_k = generate_x_hat_k_testing(;G=G,K=K,rand_init=false)
# x_hat_sq_k = generate_x_hat_sq_k_testing(;G=G,K=K,rand_init=false)
# x_hat_error = generate_x_hat_error_testing(;G=G,rand_init=false)
# x_hatk_signal = generate_x_hatk_signal_testing(;G=G,K=K,rand_init=false)
# x_hat_sq_error = generate_x_hat_sq_error_testing(;G=G,rand_init=false)
# x_hatk_sq_signal = generate_x_hatk_sq_signal_testing(;G=G,K=K,rand_init=false)
# e_γ = generate_e_γ_testing(;rand_init=false)
# Tαk = generate_Tαk_testing(;K=K,rand_init=false)
# λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = generate_hyperparamters_testing(;rand_init=false)
# λ0_vec, μ0_vec, a0_vec, b0_vec,λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init, omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,rtik_init = generate_init_vec_testing(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=G,K=K,T=T,N=N,C_t=nothing, nothing_init=false)
# λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err = generate_hyperparamters_vs12_testing(;rand_init=false)
# λ0_vec, μ0_vec, a0_vec, b0_vec,λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec,λ0k_hat_vec_init,mk_hat_vec_init,a0k_hat_vec_init,b0k_hat_vec_init,rhok_hat_vec_init, omegak_hat_vec_init,a_γ_hat_init,b_γ_hat_init,awt_hat_vec_init,bwt_hat_vec_init,a_αt_hat_vec_init,b_αt_hat_vec_init,c_ttprime_vec_init,θ_hat_vec_init,rtik_init,v_tikj_vec_init,λ0_err_hat_vec_init,m_err_hat_vec_init,a0_err_hat_vec_init,b0_err_hat_vec_init = generate_init_vec_vs12_testing(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err;G=G,K=K,T=T,N=N,C_t=nothing, nothing_init=false)
# ηkj_prior = initialize_η_kj_prior(x,K)
# Glog = generate_Glog_testing(;G=G,rand_init=false)
# logpi = generate_logpi_testing(;rand_init=false)


function test_βk_expected_value()
    rho_hat_vec1 = [0.25,0.25,0.25]
    omega_hat_vec1 = [2,2,2]
    e_βk_vec1 = βk_expected_value(rho_hat_vec1, omega_hat_vec1)
    @testset "Check when rho_hat is 0.25 for all K and omega_hat is 2 for all K; K = 3" begin
        @test length(e_βk_vec1) == length(rho_hat_vec1)+1 
        @test all(e_βk_vec1 .≈ [0.25, 0.1875, 0.140625 , 0.421875])
    end
end
function test_uk_expected_value()
    rho_hat_vec1 = [1,1]
    omega_hat_vec1 = [1,1]
    e_uk_vec1 = uk_expected_value(rho_hat_vec1, omega_hat_vec1)
    @testset "Check when rho_hat and omega_hat are vectors of 1's" begin
        @test length(e_uk_vec1) == length(rho_hat_vec1)
        @test all(e_uk_vec1 .== 1)
    end

    rho_hat_vec2 = [0.5,0.5]
    omega_hat_vec2 = [2,2]
    e_uk_vec2 = uk_expected_value(rho_hat_vec2, omega_hat_vec2)
    @testset "Check when rho_hat is 0.5 for all K and omega_hat is 2 for all K" begin
        @test length(e_uk_vec2) == length(rho_hat_vec2)
        @test all(e_uk_vec2 .== 0.5)
    end
end


function init_unittest_params_HDP()
    C_t = [2,2]
    T = length(C_t)
    K = 3
    G = 4
    α0,γ = 1.,1.
    x = [[[1,2,3,4],[4,3,2,1]], [[4,3,2,1],[3,2,1,4]]]   
    rtik = [[[0.75,0.125,0.125],[0.125,0.125,0.75]],[[0.125,0.125,0.75],[0.125,0.75,0.125]]] 
    λ0_vec =  [1.,1.5,1.,1.5]
    μ0_vec =  [2.,2.,2.,2.]
    a0_vec =  [2.0,2.5,2.0,2.5]
    b0_vec =  [0.5,0.5,0.5,0.5]
    rhok_hat_vec, omegak_hat_vec = init_params_states(K)
    θ_hat = [ones(K+1) ./(K+1)  for t in 1:T]#init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)#
    Ntk_test = update_Ntk(rtik)
    Ntk_true = [[0.875, 0.25,0.875,0.],[0.25,0.875,0.875,0.]]
    @test all(Ntk_test .≈ Ntk_true)
    Nk_test = update_Nk(rtik)
    Nk_true = [1.125, 1.125,1.75]
    @test all(Nk_test .≈ Nk_true)
    cell11_rtikxti = [rtik[1][1][1] .* x[1][1], rtik[1][1][2] .* x[1][1],rtik[1][1][3] .* x[1][1]]
    cell12_rtikxti = [rtik[1][2][1] .* x[1][2], rtik[1][2][2] .* x[1][2],rtik[1][2][3] .* x[1][2]]
    cell21_rtikxti = [rtik[2][1][1] .* x[2][1], rtik[2][1][2] .* x[2][1],rtik[2][1][3] .* x[2][1]]
    cell22_rtikxti = [rtik[2][2][1] .* x[2][2], rtik[2][2][2] .* x[2][2],rtik[2][2][3] .* x[2][2]]
    rtikxti = [[cell11_rtikxti,cell12_rtikxti],[cell21_rtikxti,cell22_rtikxti]]
    xbar_k_true = []
    for k in 1:K
        state_sum = zeros(G)
        for t in 1:T
            for i in 1:C_t[t]
                state_sum +=rtikxti[t][i][k]
            end
        end 
        push!(xbar_k_true,state_sum)   
    end
    xbar_k_test = update_xbar_k(x,rtik)
    @test all(xbar_k_test[1] .≈ xbar_k_true[1])
    @test all(xbar_k_test[2] .≈ xbar_k_true[2])
    @test all(xbar_k_test[3] .≈ xbar_k_true[3])

    cell11_rtikxti_minusXbar = [rtik[1][1][1] .* (x[1][1] .-  xbar_k_test[1]) .^2 , rtik[1][1][2] .* (x[1][1] .-  xbar_k_test[2]) .^2, rtik[1][1][3] .* (x[1][1] .- xbar_k_test[3]) .^2 ]

    cell12_rtikxti_minusXbar  = [rtik[1][2][1] .* (x[1][2] .- xbar_k_test[1]) .^2 , rtik[1][2][2] .*( x[1][2] .- xbar_k_test[2] ) .^ 2 ,rtik[1][2][3] .* (x[1][2] .- xbar_k_test[3]) .^2]

    cell21_rtikxti_minusXbar  = [rtik[2][1][1] .* (x[2][1] .- xbar_k_test[1]) .^ 2, rtik[2][1][2] .* (x[2][1] .- xbar_k_test[2]) .^ 2 ,rtik[2][1][3] .*( x[2][1] .- xbar_k_test[3]) .^2]

    cell22_rtikxti_minusXbar  = [rtik[2][2][1] .* (x[2][2] .- xbar_k_test[1]) .^2, rtik[2][2][2] .* (x[2][2] .- xbar_k_test[2]) .^2  ,rtik[2][2][3] .* (x[2][2] .- xbar_k_test[3]) .^2]

    rtikxti_minusXbar  = [[cell11_rtikxti_minusXbar,cell12_rtikxti_minusXbar],[cell21_rtikxti_minusXbar,cell22_rtikxti_minusXbar]]
    
    
    sk_true = []
    for k in 1:K
        state_sum = zeros(G)
        for t in 1:T
            for i in 1:C_t[t]
                state_sum += rtikxti_minusXbar[t][i][k]
            end
        end
        push!(sk_true, state_sum)
    end
    sk_test = update_sk_GIndepGaussian(x,xbar_k_test,rtik)
    @test all(sk_test[1] .≈ sk_true[1])
    @test all(sk_test[2] .≈ sk_true[2])
    @test all(sk_test[3] .≈ sk_true[3])

    
    λ0k_hat_vec_test = update_λ0k_hat(λ0_vec,Nk_test)
    λ0k_hat_vec_true = []
    for  k in 1:K
        λ0k_hat = Nk_test[k] .+ λ0_vec
        push!(λ0k_hat_vec_true, λ0k_hat)
    end
    @test all(λ0k_hat_vec_test[1] .≈ λ0k_hat_vec_true[1])
    @test all(λ0k_hat_vec_test[2] .≈ λ0k_hat_vec_true[2])
    @test all(λ0k_hat_vec_test[3] .≈ λ0k_hat_vec_true[3])

    a0k_hat_vec_test =  update_a0k_hat(a0_vec,Nk_test)
    a0k_hat_vec_true = []
    for k in 1:K
        a0k_hat = 1/2 .* Nk_test[k] .+ a0_vec
        push!(a0k_hat_vec_true, a0k_hat)
    end
    @test all(a0k_hat_vec_test[1] .≈ a0k_hat_vec_true[1])
    @test all(a0k_hat_vec_test[2] .≈ a0k_hat_vec_true[2])
    @test all(a0k_hat_vec_test[3] .≈ a0k_hat_vec_true[3])

    mk_hat_vec_test= update_mk_hat(λ0_vec,μ0_vec, Nk_test,xbar_k_test)
    mk_hat_vec_true = []
    for k in 1:K
        numer = λ0_vec .* μ0_vec .+ Nk_test[k] .* xbar_k_test[k]
        denom = λ0_vec .+ Nk_test[k]
        mk_hat = numer ./ denom
        push!(mk_hat_vec_true,mk_hat)
    end
    @test all(mk_hat_vec_test[1] .≈ mk_hat_vec_true[1])
    @test all(mk_hat_vec_test[2] .≈ mk_hat_vec_true[2])
    @test all(mk_hat_vec_test[3] .≈ mk_hat_vec_true[3])

    b0k_hat_vec_test = update_b0k_hat(b0_vec,λ0_vec,μ0_vec, Nk_test,xbar_k_test,sk_test)
    b0k_hat_vec_true = []
    for k in 1:K
        numer = λ0_vec * Nk_test[k] .* (xbar_k_test[k] .- μ0_vec) .^2
        denom = λ0_vec .+ Nk_test[k]
        b0k_hat = b0_vec .+ 1/2 .* sk_test[k] .+ 1/2 .* numer ./ denom
        push!(b0k_hat_vec_true,b0k_hat)
    end
    @test all(b0k_hat_vec_test[1] .≈ b0k_hat_vec_true[1])
    @test all(b0k_hat_vec_test[2] .≈ b0k_hat_vec_true[2])
    @test all(b0k_hat_vec_test[3] .≈ b0k_hat_vec_true[3])
    

    GlogTwoπ = G*log(2π)
    θ_hat = [ones(K+1) ./(K+1)  for t in 1:T]
    e_log_π_test = log_π_expected_value(θ_hat) # T by K
    e_log_π_true = []
    for t in 1:T
        e_log_π_t = []
        for k in 1:K
        e_log_π_tk = digamma.(θ_hat[t][k]) .- digamma(sum(θ_hat[t][1:K]))
        push!(e_log_π_t, e_log_π_tk )
        end
        push!(e_log_π_true,e_log_π_t)
    end
    @test all(e_log_π_test[1] .≈ e_log_π_true[1])
    @test all(e_log_π_test[2] .≈ e_log_π_true[2])


    e_log_τ_test1 = sum.(log_τ_expected_value_old.(a0k_hat_vec_test, a0k_hat_vec_test))
    e_log_τ_true = []
    for k in 1:K
        e_log_τ_kj = digamma.(a0k_hat_vec_test[k]) .- log.(b0k_hat_vec_test[k])
        e_log_τ_k = sum(e_log_τ_kj)
        push!(e_log_τ_true,e_log_τ_k)
    end
    @test all(e_log_τ_test1 .≈ e_log_τ_true)

    e_log_τ_test2 = log_τ_k_expected_value_new1(a0k_hat_vec_test, b0k_hat_vec_test)
    @test all(e_log_τ_test2 .≈ e_log_τ_true)
    ###############################################################
    e_τ_μ_kj_test = τ_μ_expected_value_old(x,λ0k_hat_vec_test,mk_hat_vec_test,a0k_hat_vec_test, b0k_hat_vec_test)
    e_τ_μ_test = [[[ sum(e_τ_μ_kj_test[t][i][k]) for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    e_τ_μ_kj_true = Vector{Vector{Vector{Vector}}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt =  Vector{Vector{Vector}}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                e_τ_μ_kjti_cell  =  a0k_hat_vec_test[k] ./  b0k_hat_vec_test[k] .*  (x[t][i] .- mk_hat_vec_test[k]) .^2 .+ 1 ./λ0k_hat_vec_test[k]
                e_τ_μ_kjti[k] = e_τ_μ_kjti_cell# push!(e_τ_μ_kjti,e_τ_μ_kjti_cell)
            end
            e_τ_μ_kjt[i] = e_τ_μ_kjti#push!(e_τ_μ_kjt,e_τ_μ_kjti)
        end
        e_τ_μ_kj_true[t] = e_τ_μ_kjt #push!(e_τ_μ_kj_true, e_τ_μ_kjt)
    end
    test_t = rand(collect(1:T))
    test_i = rand(collect(1:C_t[test_t]))
    test_k = rand(collect(1:K))
    println("for t: $test_t, i:$test_i, k:$test_k")
    @test all(e_τ_μ_kj_test[test_t][test_i][test_k] .=== e_τ_μ_kj_true[test_t][test_i][test_k])

    e_τ_μ_true = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_1 =  Vector{Vector{Float64}}(undef,cells)
        for i in 1:cells
            e_τ_μ_2 =  Vector{Float64}(undef,K)
            for k in 1:K
                e_τ_μ_3_cell  = sum(e_τ_μ_kj_true[t][i][k])
                e_τ_μ_2[k] = e_τ_μ_3_cell#push!(e_τ_μ_2,e_τ_μ_3_cell)
            end
            e_τ_μ_1[i] = e_τ_μ_2 #push!(e_τ_μ_1,e_τ_μ_2)
        end
        e_τ_μ_true[t] =e_τ_μ_1#push!(e_τ_μ_true, e_τ_μ_1)
    end
    test_t = rand(collect(1:T))
    test_i = rand(collect(1:C_t[test_t]))
    println("for t: $test_t, i:$test_i")
    @test all( e_τ_μ_test[test_t][test_i]  .== e_τ_μ_true[test_t][test_i])

    e_τ_μ_kj_test2,e_τ_μ_test2 =  τ_μ_expected_value_new1(x,λ0k_hat_vec_test,mk_hat_vec_test,a0k_hat_vec_test, b0k_hat_vec_test)
    test_t = rand(collect(1:T))
    test_i = rand(collect(1:C_t[test_t]))
    test_k = rand(collect(1:K))
    println("for t: $test_t, i:$test_i, k:$test_k")
    @test all(e_τ_μ_kj_test2[test_t][test_i][test_k] .=== e_τ_μ_kj_true[test_t][test_i][test_k])
    @test all(e_τ_μ_test2[test_t][test_i]  .== e_τ_μ_true[test_t][test_i])
    ########################################################################
    rtik_test = update_rtik(GlogTwoπ,e_log_π_test,e_log_τ_test2,e_τ_μ_test2)
    rtik_true = Vector{Vector{Vector{Float64}}}(undef, T)
    for t in 1:T
        cells = C_t[t]
        rti_true = Vector{Vector{Float64}}(undef, cells)
        for i in 1:cells
            p_tik_true = Vector{Float64}(undef, K)
            for k in 1:K
                p_tik_true[k] = e_log_π_test[t][k] .- 1/2 .*GlogTwoπ .+ 1/2 .*e_log_τ_test2[k] .- 1/2 .* e_τ_μ_test2[t][i][k]
            end
            valsum = StatsFuns.logsumexp(p_tik_true)
            valdiff = p_tik_true .- valsum
            val = exp.(valdiff)
            rti_true[i] = val
        end
        rtik_true[t] = rti_true
    end
    @test all(isprobvec(rtik_test[1][1]) && isprobvec(rtik_true[1][1]))
    @test all(isprobvec(rtik_test[1][2]) && isprobvec(rtik_true[1][2]))
    @test all(isprobvec(rtik_test[2][1]) && isprobvec(rtik_true[2][1]))
    @test all(isprobvec(rtik_test[2][2]) && isprobvec(rtik_true[2][2]))

    @test all(rtik_test[1][1] .== rtik_true[1][1])
    @test all(rtik_test[1][2] .== rtik_true[1][2])
    @test all(rtik_test[2][1] .== rtik_true[2][1])
    @test all(rtik_test[2][2] .== rtik_true[2][2])

    #######################################################################
    θ_hat_test = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk_test,α0)
    θ_hat_true = Vector{Vector{Float64}}(undef, T)
    e_βk_true = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    for t in 1:T
        θ_hat_t_true = Vector{Float64}(undef, K)
        for k in 1:K
            θ_hat_t_true[k] = Ntk_test[t][k] .+ α0 .* e_βk_true[k] 
        end
        θ_hat_true[t] = θ_hat_t_true
    end
    @test all(θ_hat_test[1][1:K] .== θ_hat_true[1])
    @test all(θ_hat_test[2][1:K] .== θ_hat_true[2])
    

    #######################################################################
    Tk_test = update_Tk(θ_hat)

    e_log_π_Tk = log_π_expected_value(θ_hat)
    Tk_true = Vector{Float64}(undef,K)
    for k in 1:K
        Tk_true_t = 0.0
        for t in 1:T
            Tk_true_t += e_log_π_Tk[t][k]
        end
        Tk_true[k] = Tk_true_t
    end

    @test all(Tk_test .== Tk_true)


    #######################################################################

end

function init_unittest_params_dHDP()
    C_t = [2,2]
    T = length(C_t)
    K = 3
    G = 4
    α0,γ = 1.,1.
    x = [[[1,2,3,4],[4,3,2,1]], [[4,3,2,1],[3,2,1,4]]]   
    rtik = [[[0.75,0.125,0.125],[0.125,0.125,0.75]],[[0.125,0.125,0.75],[0.125,0.75,0.125]]] 
    c_ttprime = [[1.0,0.0],[0.75,0.25]] #each element is of length T
    λ0_vec =  [1.,1.5,1.,1.5]
    μ0_vec =  [2.,2.,2.,2.]
    a0_vec =  [2.0,2.5,2.0,2.5]
    b0_vec =  [0.5,0.5,0.5,0.5]
    awt_hat = [0.5]
    bwt_hat = [1.25]
    rhok_hat_vec, omegak_hat_vec = init_params_states(K)
    θ_hat = [ones(K+1) ./(K+1)  for t in 1:T]#init_θ_hat_tk(T,K,rhok_hat_vec, omegak_hat_vec;rand_init = false)#
    Ntk_test = update_Ntk(rtik)
    Ntk_true = [[0.875, 0.25,0.875,0.],[0.25,0.875,0.875,0.]]
    @test all(Ntk_test .≈ Ntk_true)
    Nk_test = update_Nk(rtik)
    Nk_true = [1.125, 1.125,1.75]
    @test all(Nk_test .≈ Nk_true)
    cell11_rtikxti = [rtik[1][1][1] .* x[1][1], rtik[1][1][2] .* x[1][1],rtik[1][1][3] .* x[1][1]]
    cell12_rtikxti = [rtik[1][2][1] .* x[1][2], rtik[1][2][2] .* x[1][2],rtik[1][2][3] .* x[1][2]]
    cell21_rtikxti = [rtik[2][1][1] .* x[2][1], rtik[2][1][2] .* x[2][1],rtik[2][1][3] .* x[2][1]]
    cell22_rtikxti = [rtik[2][2][1] .* x[2][2], rtik[2][2][2] .* x[2][2],rtik[2][2][3] .* x[2][2]]
    rtikxti = [[cell11_rtikxti,cell12_rtikxti],[cell21_rtikxti,cell22_rtikxti]]
    xbar_k_true = []
    for k in 1:K
        state_sum = zeros(G)
        for t in 1:T
            for i in 1:C_t[t]
                state_sum +=rtikxti[t][i][k]
            end
        end 
        push!(xbar_k_true,state_sum)   
    end
    xbar_k_test = update_xbar_k(x,rtik)
    @test all(xbar_k_test[1] .≈ xbar_k_true[1])
    @test all(xbar_k_test[2] .≈ xbar_k_true[2])
    @test all(xbar_k_test[3] .≈ xbar_k_true[3])

    cell11_rtikxti_minusXbar = [rtik[1][1][1] .* (x[1][1] .-  xbar_k_test[1]) .^2 , rtik[1][1][2] .* (x[1][1] .-  xbar_k_test[2]) .^2, rtik[1][1][3] .* (x[1][1] .- xbar_k_test[3]) .^2 ]

    cell12_rtikxti_minusXbar  = [rtik[1][2][1] .* (x[1][2] .- xbar_k_test[1]) .^2 , rtik[1][2][2] .*( x[1][2] .- xbar_k_test[2] ) .^ 2 ,rtik[1][2][3] .* (x[1][2] .- xbar_k_test[3]) .^2]

    cell21_rtikxti_minusXbar  = [rtik[2][1][1] .* (x[2][1] .- xbar_k_test[1]) .^ 2, rtik[2][1][2] .* (x[2][1] .- xbar_k_test[2]) .^ 2 ,rtik[2][1][3] .*( x[2][1] .- xbar_k_test[3]) .^2]

    cell22_rtikxti_minusXbar  = [rtik[2][2][1] .* (x[2][2] .- xbar_k_test[1]) .^2, rtik[2][2][2] .* (x[2][2] .- xbar_k_test[2]) .^2  ,rtik[2][2][3] .* (x[2][2] .- xbar_k_test[3]) .^2]

    rtikxti_minusXbar  = [[cell11_rtikxti_minusXbar,cell12_rtikxti_minusXbar],[cell21_rtikxti_minusXbar,cell22_rtikxti_minusXbar]]
    
    
    sk_true = []
    for k in 1:K
        state_sum = zeros(G)
        for t in 1:T
            for i in 1:C_t[t]
                state_sum += rtikxti_minusXbar[t][i][k]
            end
        end
        push!(sk_true, state_sum)
    end
    sk_test = update_sk_GIndepGaussian(x,xbar_k_test,rtik)
    @test all(sk_test[1] .≈ sk_true[1])
    @test all(sk_test[2] .≈ sk_true[2])
    @test all(sk_test[3] .≈ sk_true[3])

    
    λ0k_hat_vec_test = update_λ0k_hat(λ0_vec,Nk_test)
    λ0k_hat_vec_true = []
    for  k in 1:K
        λ0k_hat = Nk_test[k] .+ λ0_vec
        push!(λ0k_hat_vec_true, λ0k_hat)
    end
    @test all(λ0k_hat_vec_test[1] .≈ λ0k_hat_vec_true[1])
    @test all(λ0k_hat_vec_test[2] .≈ λ0k_hat_vec_true[2])
    @test all(λ0k_hat_vec_test[3] .≈ λ0k_hat_vec_true[3])

    a0k_hat_vec_test =  update_a0k_hat(a0_vec,Nk_test)
    a0k_hat_vec_true = []
    for k in 1:K
        a0k_hat = 1/2 .* Nk_test[k] .+ a0_vec
        push!(a0k_hat_vec_true, a0k_hat)
    end
    @test all(a0k_hat_vec_test[1] .≈ a0k_hat_vec_true[1])
    @test all(a0k_hat_vec_test[2] .≈ a0k_hat_vec_true[2])
    @test all(a0k_hat_vec_test[3] .≈ a0k_hat_vec_true[3])

    mk_hat_vec_test= update_mk_hat(λ0_vec,μ0_vec, Nk_test,xbar_k_test)
    mk_hat_vec_true = []
    for k in 1:K
        numer = λ0_vec .* μ0_vec .+ Nk_test[k] .* xbar_k_test[k]
        denom = λ0_vec .+ Nk_test[k]
        mk_hat = numer ./ denom
        push!(mk_hat_vec_true,mk_hat)
    end
    @test all(mk_hat_vec_test[1] .≈ mk_hat_vec_true[1])
    @test all(mk_hat_vec_test[2] .≈ mk_hat_vec_true[2])
    @test all(mk_hat_vec_test[3] .≈ mk_hat_vec_true[3])

    b0k_hat_vec_test = update_b0k_hat(b0_vec,λ0_vec,μ0_vec, Nk_test,xbar_k_test,sk_test)
    b0k_hat_vec_true = []
    for k in 1:K
        numer = λ0_vec * Nk_test[k] .* (xbar_k_test[k] .- μ0_vec) .^2
        denom = λ0_vec .+ Nk_test[k]
        b0k_hat = b0_vec .+ 1/2 .* sk_test[k] .+ 1/2 .* numer ./ denom
        push!(b0k_hat_vec_true,b0k_hat)
    end
    @test all(b0k_hat_vec_test[1] .≈ b0k_hat_vec_true[1])
    @test all(b0k_hat_vec_test[2] .≈ b0k_hat_vec_true[2])
    @test all(b0k_hat_vec_test[3] .≈ b0k_hat_vec_true[3])
    

    GlogTwoπ = G*log(2π)
    θ_hat = [ones(K+1) ./(K+1)  for t in 1:T]
    e_log_π_test = log_π_expected_value(θ_hat) # T by K
    e_log_π_true = []
    for t in 1:T
        e_log_π_t = []
        for k in 1:K
        e_log_π_tk = digamma.(θ_hat[t][k]) .- digamma(sum(θ_hat[t][1:K]))
        push!(e_log_π_t, e_log_π_tk )
        end
        push!(e_log_π_true,e_log_π_t)
    end
    @test all(e_log_π_test[1] .≈ e_log_π_true[1])
    @test all(e_log_π_test[2] .≈ e_log_π_true[2])


    e_log_τ_test1 = sum.(log_τ_expected_value_old.(a0k_hat_vec_test, a0k_hat_vec_test))
    e_log_τ_true = []
    for k in 1:K
        e_log_τ_kj = digamma.(a0k_hat_vec_test[k]) .- log.(b0k_hat_vec_test[k])
        e_log_τ_k = sum(e_log_τ_kj)
        push!(e_log_τ_true,e_log_τ_k)
    end
    @test all(e_log_τ_test1 .≈ e_log_τ_true)

    e_log_τ_test2 = log_τ_k_expected_value_new1(a0k_hat_vec_test, b0k_hat_vec_test)
    @test all(e_log_τ_test2 .≈ e_log_τ_true)
    ###############################################################
    e_τ_μ_kj_test = τ_μ_expected_value_old(x,λ0k_hat_vec_test,mk_hat_vec_test,a0k_hat_vec_test, b0k_hat_vec_test)
    e_τ_μ_test = [[[ sum(e_τ_μ_kj_test[t][i][k]) for k in 1:K] for i in 1:C_t[t]] for t in 1:T]
    e_τ_μ_kj_true = Vector{Vector{Vector{Vector}}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_kjt =  Vector{Vector{Vector}}(undef,cells)
        for i in 1:cells
            e_τ_μ_kjti = Vector{Vector{Float64}}(undef,K)
            for k in 1:K
                e_τ_μ_kjti_cell  =  a0k_hat_vec_test[k] ./  b0k_hat_vec_test[k] .*  (x[t][i] .- mk_hat_vec_test[k]) .^2 .+ 1 ./λ0k_hat_vec_test[k]
                e_τ_μ_kjti[k] = e_τ_μ_kjti_cell# push!(e_τ_μ_kjti,e_τ_μ_kjti_cell)
            end
            e_τ_μ_kjt[i] = e_τ_μ_kjti#push!(e_τ_μ_kjt,e_τ_μ_kjti)
        end
        e_τ_μ_kj_true[t] = e_τ_μ_kjt #push!(e_τ_μ_kj_true, e_τ_μ_kjt)
    end
    test_t = rand(collect(1:T))
    test_i = rand(collect(1:C_t[test_t]))
    test_k = rand(collect(1:K))
    println("for t: $test_t, i:$test_i, k:$test_k")
    @test all(e_τ_μ_kj_test[test_t][test_i][test_k] .=== e_τ_μ_kj_true[test_t][test_i][test_k])

    e_τ_μ_true = Vector{Vector{Vector{Float64}}}(undef,T)
    for t in 1:T
        cells = C_t[t]
        e_τ_μ_1 =  Vector{Vector{Float64}}(undef,cells)
        for i in 1:cells
            e_τ_μ_2 =  Vector{Float64}(undef,K)
            for k in 1:K
                e_τ_μ_3_cell  = sum(e_τ_μ_kj_true[t][i][k])
                e_τ_μ_2[k] = e_τ_μ_3_cell#push!(e_τ_μ_2,e_τ_μ_3_cell)
            end
            e_τ_μ_1[i] = e_τ_μ_2 #push!(e_τ_μ_1,e_τ_μ_2)
        end
        e_τ_μ_true[t] =e_τ_μ_1#push!(e_τ_μ_true, e_τ_μ_1)
    end
    test_t = rand(collect(1:T))
    test_i = rand(collect(1:C_t[test_t]))
    println("for t: $test_t, i:$test_i")
    @test all( e_τ_μ_test[test_t][test_i]  .== e_τ_μ_true[test_t][test_i])

    e_τ_μ_kj_test2,e_τ_μ_test2 =  τ_μ_expected_value_new1(x,λ0k_hat_vec_test,mk_hat_vec_test,a0k_hat_vec_test, b0k_hat_vec_test)
    test_t = rand(collect(1:T))
    test_i = rand(collect(1:C_t[test_t]))
    test_k = rand(collect(1:K))
    println("for t: $test_t, i:$test_i, k:$test_k")
    @test all(e_τ_μ_kj_test2[test_t][test_i][test_k] .=== e_τ_μ_kj_true[test_t][test_i][test_k])
    @test all(e_τ_μ_test2[test_t][test_i]  .== e_τ_μ_true[test_t][test_i])
    ########################################################################
    rtik_test = update_rtik(GlogTwoπ,e_log_π_test,e_log_τ_test2,e_τ_μ_test2)
    rtik_true = Vector{Vector{Vector{Float64}}}(undef, T)
    for t in 1:T
        cells = C_t[t]
        rti_true = Vector{Vector{Float64}}(undef, cells)
        for i in 1:cells
            p_tik_true = Vector{Float64}(undef, K)
            for k in 1:K
                p_tik_true[k] = e_log_π_test[t][k] .- 1/2 .*GlogTwoπ .+ 1/2 .*e_log_τ_test2[k] .- 1/2 .* e_τ_μ_test2[t][i][k]
            end
            valsum = StatsFuns.logsumexp(p_tik_true)
            valdiff = p_tik_true .- valsum
            val = exp.(valdiff)
            rti_true[i] = val
        end
        rtik_true[t] = rti_true
    end
    @test all(isprobvec(rtik_test[1][1]) && isprobvec(rtik_true[1][1]))
    @test all(isprobvec(rtik_test[1][2]) && isprobvec(rtik_true[1][2]))
    @test all(isprobvec(rtik_test[2][1]) && isprobvec(rtik_true[2][1]))
    @test all(isprobvec(rtik_test[2][2]) && isprobvec(rtik_true[2][2]))

    @test all(rtik_test[1][1] .== rtik_true[1][1])
    @test all(rtik_test[1][2] .== rtik_true[1][2])
    @test all(rtik_test[2][1] .== rtik_true[2][1])
    @test all(rtik_test[2][2] .== rtik_true[2][2])

    #######################################################################
    θ_hat_test = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk_test,α0)
    θ_hat_true = Vector{Vector{Float64}}(undef, T)
    e_βk_true = βk_expected_value(rhok_hat_vec, omegak_hat_vec)
    for t in 1:T
        θ_hat_t_true = Vector{Float64}(undef, K)
        for k in 1:K
            θ_hat_t_true[k] = Ntk_test[t][k] .+ α0 .* e_βk_true[k] 
        end
        θ_hat_true[t] = θ_hat_t_true
    end
    @test all(θ_hat_test[1][1:K] .== θ_hat_true[1])
    @test all(θ_hat_test[2][1:K] .== θ_hat_true[2])
    

    #######################################################################
    Tk_test = update_Tk(θ_hat)

    e_log_π_Tk = log_π_expected_value(θ_hat)
    Tk_true = Vector{Float64}(undef,K)
    for k in 1:K
        Tk_true_t = 0.0
        for t in 1:T
            Tk_true_t += e_log_π_Tk[t][k]
        end
        Tk_true[k] = Tk_true_t
    end

    @test all(Tk_test .== Tk_true)


    #######################################################################
    e_log_tilde_wt_vec = log_tilde_wt_expected_value.(awt_hat,bwt_hat)
    e_log_minus_tilde_wt_vec = log_minus_tilde_wt_expected_value.(awt_hat,bwt_hat)
    e_log_w2_true = [0+e_log_minus_tilde_wt_vec[1], e_log_tilde_wt_vec[1]]
    e_log_w2_test = log_w_ttprime_expected_value(awt_hat,bwt_hat)[2]
    @test all(e_log_w2_true .== e_log_w2_test)

    
    awt_hat2,bwt_hat2 = [0.5,1.],[1.,0.5]
    e_log_tilde_wt_vec2 = log_tilde_wt_expected_value.(awt_hat2,bwt_hat2)
    e_log_minus_tilde_wt_vec2 = log_minus_tilde_wt_expected_value.(awt_hat2,bwt_hat2)
    e_log_w2_true2 = [0+e_log_minus_tilde_wt_vec2[1], e_log_tilde_wt_vec2[1]]
    e_log_w2_test2 = log_w_ttprime_expected_value(awt_hat2,bwt_hat2)[2]
    @test all(e_log_w2_true2 .== e_log_w2_test2)

    e_log_w3_true2 = [0+e_log_minus_tilde_wt_vec2[1] + e_log_minus_tilde_wt_vec2[2], e_log_tilde_wt_vec2[1] + e_log_minus_tilde_wt_vec2[2], e_log_tilde_wt_vec2[2]]
    e_log_w3_test2 = log_w_ttprime_expected_value(awt_hat2,bwt_hat2)[3]
    @test all(e_log_w3_true2 .== e_log_w3_test2)


end