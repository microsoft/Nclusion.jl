"""
    init_c!(c_init,T,N_t;rand_init = false)
This function initializes the c variable for the variational inference algorithm. If c_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_c!(c_init,T,N_t;rand_init = false)
    if isnothing(c_init) && rand_init
        c_init = [[[rand(Dirichlet(ones(t) ./t));zeros(T-t)] for n in 1:N_t[t] ] for t in 1:T]
    elseif isnothing(c_init) && !rand_init
        c_init = [[[ones(t) ./t;zeros(T-t)] for n in 1:N_t[t]]  for t in 1:T]
    end
    return c_init
end
"""
    init_r!(r_init,K,T,N_t;rand_init = false)
This function initializes the r variable for the variational inference algorithm. If r_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
""" 
function init_r!(r_init,K,T,N_t;rand_init = false)
    if isnothing(r_init) && rand_init
        r_init = [[[rand(Dirichlet(ones(K) ./K));zeros(1) ] for i in 1:N_t[t]] for t in 1:T]
    elseif  isnothing(r_init) && !rand_init
        r_init = [[[ones(K) ./K ;zeros(1)  ] for i in 1:N_t[t]] for t in 1:T]
    end
    return r_init
end

"""
    init_w1!(w1_init,T;rand_init = false)
This function initializes the w1 variable for the variational inference algorithm. If w1_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
""" 
function init_w1!(w1_init,T;rand_init = false)
    if isnothing(w1_init) && rand_init
        w1_init = exp.(randn(T))
    elseif isnothing(w1_init) && !rand_init
        w1_init = ones(T)
    end
    w1_init[1] = 1.0
    return w1_init
end

"""
    init_w2!(w2_init,T;rand_init = false)
This function initializes the w2 variable for the variational inference algorithm. If w2_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_w2!(w2_init,T;rand_init = false)
    if isnothing(w2_init) && rand_init
        w2_init = exp.(randn(T))
    elseif isnothing(w2_init) && !rand_init
        w2_init = ones(T)
    end
    w2_init[1] = 1.0
    return w2_init
end

"""
    init_d!(d_init,K,T;rand_init = false,uniform_theta_init=true, g1_init = nothing, g2_init= nothing)
This function initializes the d variable for the variational inference algorithm. If d_init is not provided, it can be initialized randomly, uniformly, or based on the g1 and g2 parameters depending on the flags provided.
"""
function init_d!(d_init,K,T;rand_init = false,uniform_theta_init=true, g1_init = nothing, g2_init= nothing)
    if isnothing(d_init)
        if uniform_theta_init
            d_init = [ones(K+1) ./(K+1)  for t in 1:T]
        else
            if rand_init
                d_init = [rand(K+1) for t in 1:T]
            else
                d_init = init_d_k(T,g1_init, g2_init);
            end
        end
    end
    return d_init
end

"""
    init_d_k(T,g1_init, g2_init)
This function initializes the d variable for the variational inference algorithm based on the g1 and g2 parameters.
"""
function init_d_k(T,g1_init, g2_init)
    d_vec = [[expectation_sbk(k,K+1,g1_init,g2_init; use_log = true) for i in 1:K+1] for t in 1:T]
    return d_vec
end

"""
    init_g1!(g1_init,K;rand_init = false)
This function initializes the g1 variable for the variational inference algorithm. If g1_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_g1!(g1_init,K;rand_init = false)
    if isnothing(g1_init) && rand_init
        g1_init = logistic.(randn(K))
    elseif isnothing(g1_init) && !rand_init
        g1_init = ones(K)
    end
    return  g1_init
end

"""
    init_g2!(g2_init,K;rand_init = false)
This function initializes the g2 variable for the variational inference algorithm. If g2_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_g2!(g2_init,K;rand_init = false)
    if isnothing(g2_init) && rand_init
        g2_init = exp.(randn(K))
    elseif isnothing(g2_init) && !rand_init
        g2_init = ones(K)
    end
    return  g2_init
end

"""
    init_m_mu!(m_mu_init,K,J;rand_init = false)
This function initializes the m_mu variable for the variational inference algorithm. If m_mu_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_m_mu!(m_mu_init,K,J;rand_init = false)
    mu0_vec = zeros(J)
    if isnothing(m_mu_init) && rand_init
        m_mu_init = [randn(J) for k in 1:K]
    elseif isnothing(m_mu_init) && !rand_init
        m_mu_init = [mu0_vec for k in 1:K]
    end
    return m_mu_init
end

"""
    init_s_sq_mu!(s_sq_mu_init,K,J;rand_init = false)
This function initializes the s_sq_mu variable for the variational inference algorithm. If s_sq_mu_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_s_sq_mu!(s_sq_mu_init,K,J;rand_init = false)
    s_sq0_vec = ones(J)
    if isnothing(s_sq_mu_init) && rand_init
        s_sq_mu_init = [exp.(randn(J)) for k in 1:K]
    elseif isnothing(s_sq_mu_init) && !rand_init
        s_sq_mu_init = [s_sq0_vec for k in 1:K]
    end
    return s_sq_mu_init
end
"""
    init_y!(y_init,K,J;rand_init = false)   
This function initializes the y variable for the variational inference algorithm. If y_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_y!(y_init,K,J;rand_init = false)
    y0_vec = 0.5*ones(J)
    if isnothing(y_init) && rand_init
        y_init = [logistic.(randn(J)) for k in 1:K]
    elseif isnothing(y_init) && !rand_init
        y_init = [y0_vec for k in 1:K]
    end
    return y_init
end

"""
    init_h1!(h1_init,K,J;rand_init = false)
This function initializes the h1 variable for the variational inference algorithm. If h1_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
""" 
function init_h1!(h1_init,K,J;rand_init = false)
    a, b = 0.0, 1.0
    if isnothing(h1_init) && rand_init
        h1_init = rand(K) #exp.(randn(K))
    elseif isnothing(h1_init) && !rand_init
        h1_init = ones(K)
    end
    return h1_init
end

"""
    init_h2!(h2_init,K,J;rand_init = false)
This function initializes the h2 variable for the variational inference algorithm. If h2_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_h2!(h2_init,K,J;rand_init = false)
    if isnothing(h2_init) && rand_init
        h2_init = rand(K) #exp.(randn(K))
    elseif isnothing(h2_init) && !rand_init
        h2_init = ones(K)
    end
    return h2_init
end

"""
    init_a!(a_init,K,J;rand_init = false)
This function initializes the a variable for the variational inference algorithm. If a_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_a!(a_init,K,J;rand_init = false)
    if isnothing(a_init) && rand_init
        a_init = exp.(randn(J))
    elseif isnothing(a_init) && !rand_init
        a_init = ones(J)
    end
    return a_init
end

"""
    init_b!(b_init,K,J;rand_init = false)
This function initializes the b variable for the variational inference algorithm. If b_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_b!(b_init,K,J;rand_init = false)
    if isnothing(b_init) && rand_init
        b_init = exp.(randn(J))
    elseif isnothing(b_init) && !rand_init
        b_init = ones(J)
    end
    return b_init
end
"""
    init_m_nu!(m_nu_init,K,J;rand_init = false)
This function initializes the m_nu variable for the variational inference algorithm. If m_nu_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_m_nu!(m_nu_init,K,J;rand_init = false)
    nu0_vec = zeros(J)
    if isnothing(m_nu_init) && rand_init
        m_nu_init = randn(J)
    elseif isnothing(m_nu_init) && !rand_init
        m_nu_init = nu0_vec
    end
    return m_nu_init
end

"""
    init_s_sq_nu!(s_sq_nu_init,K,J;rand_init = false)
This function initializes the s_sq_nu variable for the variational inference algorithm. If s_sq_nu_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_s_sq_nu!(s_sq_nu_init,K,J;rand_init = false)
    s_sq0_nu_vec = ones(J)
    if isnothing(s_sq_nu_init) && rand_init
        s_sq_nu_init = exp.(randn(J))
    elseif isnothing(s_sq_nu_init) && !rand_init
        s_sq_nu_init = s_sq0_nu_vec
    end
    return s_sq_nu_init
end

"""
    init_u!(u_init,K,J;rand_init = false)
This function initializes the u variable for the variational inference algorithm. If u_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_u!(u_init,K,J;rand_init = false)
    if isnothing(u_init) && rand_init
        u_init = exp.(randn(float_type,1))
    elseif isnothing(u_init) && !rand_init
        u_init = 1.0
    end
    return u_init
end

"""
    init_v!(v_init,K,J;rand_init = false)
This function initializes the v variable for the variational inference algorithm. If v_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_v!(v_init,K,J;rand_init = false)
    if isnothing(v_init) && rand_init
        v_init = exp(randn())
    elseif isnothing(v_init) && !rand_init
        v_init = 1.0
    end
    return v_init
end

"""
    generate_fake_cells(data_input,dataparams,modelparams,SEED=2020)
This function generates fake cell data for testing purposes. It creates a list of CellFeature objects with random values for their attributes based on the provided data input, data parameters, and model parameters.
"""
function generate_fake_cells(data_input,dataparams,modelparams,SEED=2020)
    Random.seed!(SEED)
    I = dataparams.I
    T = dataparams.T
    N_t = dataparams.N_t
    K = modelparams.K
    N = dataparams.N
    cells = [CellFeature(i,t,n,K,T,data_input[i][t][n]) for i in 1:I for t in 1:T[i] for n in 1:N_t[i][t]]
    a, b = 0.0, 10.0
    for n in 1:N
        i = first(cells[n].i) # Alternatively, i = dataparams.LinearAddress[n][1]
        t = first(cells[n].t) # Alternatively, t = dataparams.LinearAddress[n][2]
        unnormalized_r = rand(K) .* (b - a) .+ a  # Generates a 1D array of K random numbers between a and b
        normalized_r = normToProb(unnormalized_r)  # Normalizes the array to sum to 1
        unnormalized_c = rand(t) .* (b - a) .+ a  # Generates a 1D array of T random numbers between a and b
        normalized_c = normToProb(unnormalized_c)  # Normalizes the array to sum to 1
        cells[n].r .= [normalized_r;zeros(1)]  # Generates a 1D array of K+1 random integers between 1 and 10
        cells[n].c .= [normalized_c;zeros(T[i]-t)]  # Generates a 1D array of T random integers between 1 and 10
        cells[n].cache .= 0.0  # Generates a 1D array of J random integers between 1 and 10
    end
    return cells
end

"""
    generate_fake_time_conditions(T,K,N_t,SEED=2020,float_type=Float64,condition_update_neighbors=nothing,condition_network_neighbors=nothing)
This function generates fake time condition data for testing purposes. It creates a list of ConditionFeature and MatrixConditionFeature objects with random values for their attributes based on the provided time points, number of clusters, and number of cells at each time point. Optional parameters allow for custom neighbor structures.
"""
function generate_fake_time_conditions(T,K,N_t,SEED=2020,float_type=Float64,condition_update_neighbors=nothing,condition_network_neighbors=nothing)
    Random.seed!(SEED)
    Kplus = K+1;
    I = length(T)
    T_all = sum(T)
    if isnothing(condition_update_neighbors)
        condition_update_neighbors=get_linear_time_condition_update_neighbors(N_t;get_ragged_array=true)
    end
    if isnothing(condition_network_neighbors)
        condition_network_neighbors=get_linear_time_condition_network_neighbors(N_t;get_ragged_array=true)
    end
    conditions = [ConditionFeature(i,t,K,T[i],condition_update_neighbors[i][t],condition_network_neighbors[i][t];float_type=float_type) for i in 1:I for t in 1:T[i]];
    matrixconditions = [MatrixConditionFeature(i,t,K,T,condition_update_neighbors[i][t],condition_network_neighbors[i][t];float_type=float_type) for i in 1:I for t in 1:T[i]];
    it=1
    for i in 1:I
        for t in 1:T[i]
            conditions[it].d .= rand(1:10, Kplus)  # Generates a 1D array of K+1 random integers between 1 and 10
            conditions[it].d_sum[1] = sum(conditions[t].d)
            matrixconditions[it].cache .= 0.0
            conditions[it].Ctt .= rand(1:10, T[i])  # Generates a 1D array of T random integers between 1 and 10
            # conditions[t].time_cache .= 0.0
            matrixconditions[it].CNtk .= rand(1:10,T[i],Kplus)
            matrixconditions[it].suffstats_cache .= zeros(T[i],Kplus)
            conditions[it].w1[1] = 1.0
            conditions[it].w2[1] = 1.0
            it += 1
        end
    end
    return conditions,matrixconditions
end

"""
    generate_fake_clusters(K,modelparams,SEED=2020,float_type=Float64)
This function generates fake cluster data for testing purposes. It creates a list of ClusterFeature objects with random values for their attributes based on the provided number of clusters, model parameters, and seed.
"""
function generate_fake_clusters(K,modelparams,SEED=2020,float_type=Float64)
    Random.seed!(SEED)
    Kplus = K+1;
    clusters = [ClusterFeature(k,J;float_type=float_type) for k in 1:Kplus];
    a, b = 0.0, 10.0
    for k in 1:K
        clusters[k].m_mu .= randn(J)  # Generates a 1D array of J random N(0,1) distributed numbers
        clusters[k].s_sq_mu .= rand(J) .* (b - a) .+ a  # Generates a random number between a and b
        clusters[k].y .= logistic.(randn(J))  # Generates a 1D array of J random probabilities between 0 and 1
        clusters[k].h1 .= rand(J) .* (b - a) .+ a # Generates a random number between a and b
        clusters[k].h2 .= rand(J) .* (b - a) .+ a # Generates a random number between a and b
        clusters[k].Nk[1] = rand() * (b - a) + a 
        clusters[k].x_hat .= randn(J)  # Generates a 1D array of J random N(0,1) distributed numbers
        clusters[k].x_hat_sq .=  rand(J) .* (b - a) .+ a  # Generates a random number between a and b
        clusters[k].g1[1] = rand() # Generates a random probability between 0 and 1
        clusters[k].g2[1] = rand() * (b - a) + a # Generates a random number between a and b
        clusters[k].u .= rand(J) .* (b - a) .+ a  # Generates a random number between a and b
        clusters[k].v .= rand(J) .* (b - a) .+ a  # Generates a random number between a and b
        clusters[k].a .= rand(J) .* (b - a) .+ a  # Generates a random number between a and b
        clusters[k].b .= rand(J) .* (b - a) .+ a  # Generates a random number between a and b
        clusters[k].alpha_Tk[1] = rand() * (b - a) + a # Generates a random number between a and b
        clusters[k].cache .= 0.0  # Generates a 1D array of J random integers between 1 and 10
    end
    clusters[Kplus].m_mu .= 0.0  # Generates a 1D array of J zeros
    clusters[Kplus].s_sq_mu .= 1.0  # Generates a 1D array of J ones
    clusters[Kplus].y .= e_eta(modelparams.varphi1[1],modelparams.varphi2[1])
    clusters[Kplus].h1 .= modelparams.varphi1[1] .* ones(J)
    clusters[Kplus].h2 .= modelparams.varphi2[1] .* ones(J)
    clusters[Kplus].u .= modelparams.kappa1 .* ones(J)
    clusters[Kplus].v .= modelparams.kappa2 .* ones(J)
    clusters[Kplus].a .= modelparams.xi1 .* ones(J)
    clusters[Kplus].b .= modelparams.xi2 .* ones(J)
    clusters[Kplus].Nk[1] = 0.0 # Generates a 1D array of J zeros
    clusters[Kplus].x_hat .= 0.0  # Generates a 1D array of J zeros
    clusters[Kplus].x_hat_sq .= 1.0  # Generates a 1D array of J ones
    clusters[Kplus].g1[1] = 1.0
    clusters[Kplus].g2[1] = modelparams.gamma0[1]
    clusters[Kplus].alpha_Tk[1] = 1.0
    clusters[Kplus].cache .= 0.0  # Generates a 1D array of J zeros
    return clusters
end


"""
    generate_fake_model_params(data_input,K,SEED=2020)  
This function generates fake model parameters for testing purposes. It creates a ModelParameterFeature object with random values for its attributes based on the provided data input, number of clusters, and seed.
"""
function generate_fake_model_params(data_input,K,SEED=2020)
    Random.seed!(SEED)
    a, b = 0.0, 10.0
    V = eltype(data_input[1][1][1])
    I = length(data_input)
    T = [length(data_input[i]) for i in 1:I]
    J = length(data_input[1][1][1])
    alpha0 = [[rand() * (b - a) + a for t in 1:T[i]] for i in 1:I]
    gamma0 = rand() * (b - a) + a 
    phi1 = rand() * (b - a) + a 
    phi2 = rand() * (b - a) + a 
    kappa1 = rand() * (b - a) + a 
    kappa2 = rand() * (b - a) + a 
    xi1 = rand() * (b - a) + a 
    xi2 = rand() * (b - a) + a 
    varphi1 = rand() * (b - a) + a 
    varphi2 = rand() * (b - a) + a 
    nu0 = [randn() for j in 1:J]
    sigma_sq_nu = [rand() * (b - a) + a for j in 1:J]
    num_iter = rand(1:10) 
    uniform_theta_init = false
    rand_init = false
    change_seeds = false
    significance_prop = 0.5
    min_number_cells=100
    min_percent_cells=.100
    min_percent_of_genes=0.05
    max_percent_of_genes=0.75
    init_seed = 2020
    return ModelParameterFeature(data_input,K,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,num_iter,uniform_theta_init,rand_init,change_seeds,init_seed)
end

"""
    generate_fake_data_params(data_input,SEED=2020)
This function generates fake data parameters for testing purposes. It creates a DataFeature object based on the provided data input and seed.
"""
function generate_fake_data_params(data_input,SEED=2020)
    Random.seed!(SEED)
    return DataFeature(data_input)
end


"""
    generate_fake_inputs(I,TMax,J,K;SEED=2020,float_type=Float64,NMax=10,guarantee_an_idividual_with_singleton_timepoint=true)
This function generates fake inputs for testing purposes. It creates a list of fake data, data parameters, model parameters, time conditions, and cluster data based on the provided dimensions and seed.
"""
function generate_fake_inputs(I,TMax,J,K;SEED=2020,float_type=Float64,NMax=10,guarantee_an_idividual_with_singleton_timepoint=true)
    T = generate_fake_T(I,TMax,SEED,guarantee_an_idividual_with_singleton_timepoint);
    data_input = generate_fake_dataset(I,T,J,SEED;NMax=NMax);
    dataparams = generate_fake_data_params(data_input,SEED);
    N = dataparams.N
    N_t = dataparams.N_t
    LinearAddress = dataparams.LinearAddress
    TimeRanges = dataparams.TimeRanges
    modelparams = generate_fake_model_params(data_input,K,SEED);
    conditions,matrixconditions = generate_fake_time_conditions(T,K,N_t,SEED);
    cells = generate_fake_cells(data_input,dataparams,modelparams,SEED);
    clusters = generate_fake_clusters(K,modelparams,SEED,float_type);
    return I,T,J,K,N,N_t,LinearAddress,TimeRanges,data_input,dataparams,modelparams,conditions,matrixconditions,cells,clusters# = generate_fake_inputs(I,TMax,J,K,;SEED=SEED,float_type=Float64,NMax=NMax);
end

"""
 init_params_states(K)
This function initializes the parameters for the states in the variational inference algorithm. It returns two vectors, rho_hat_vec and omega_hat_vec, both of length K, with predefined values.
"""
function init_params_states(K)
    rho_hat_vec = 0.25 .* ones(Float64,K)
    omega_hat_vec = 2 .* ones(Float64,K)
    return rho_hat_vec, omega_hat_vec
end

######################################################


"""
    init_mk_hat!(mk_hat_init,x,K,G;rand_init = false)
This function initializes the mean vector for each cluster in the variational inference algorithm. It returns a vector of mean vectors, mk_hat_init, based on the provided data input, number of clusters, number of genes, and initialization type.
"""
function init_mk_hat!(mk_hat_init,x,K,G;rand_init = false)
    μ0_vec = ones(G)
    if isnothing(mk_hat_init) && rand_init
        mk_hat_init = [rand(Uniform( minimum(reduce(vcat,reduce(vcat,x)))-1,maximum(reduce(vcat,reduce(vcat,x)))+1),length(μ0_vec)) for k in 1:K]
    elseif isnothing(mk_hat_init) && !rand_init
        mk_hat_init = [μ0_vec for k in 1:K]
    end
    return mk_hat_init
end

"""
    init_λ_sq_vec!(λ_sq_init,G;rand_init = false, lo=0,hi=1)
This function initializes the lambda squared vector for each gene in the variational inference algorithm. It returns a vector of lambda squared values, λ_sq_init, based on the provided number of genes and initialization type.
"""
function init_λ_sq_vec!(λ_sq_init,G;rand_init = false, lo=0,hi=1)
    λ_sq_vec = ones(G)
    if isnothing(λ_sq_init) && rand_init
        λ_sq_init = rand(Uniform(lo,hi),length(λ_sq_vec))
    elseif isnothing(λ_sq_init) && !rand_init
        λ_sq_init = λ_sq_vec
    end
    return λ_sq_init
end

"""
    init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = false, lo=0,hi=1)
This function initializes the sigma squared vector for each cluster in the variational inference algorithm. It returns a vector of sigma squared values, σ_sq_k_init, based on the provided number of clusters, number of genes, and initialization type.
"""
function init_σ_sq_k_vec!(σ_sq_k_init,K,G;rand_init = false, lo=0,hi=1)
    σ_sq_k_vec = ones(G)
    if isnothing(σ_sq_k_init) && rand_init
        σ_sq_k_init = [rand(Uniform(lo,hi),length(σ_sq_k_vec)) for k in 1:K]
    elseif isnothing(σ_sq_k_init) && !rand_init
        σ_sq_k_init = [σ_sq_k_vec for k in 1:K] #
    end
    return σ_sq_k_init
end

"""
    init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = false, lo=0,hi=1)
This function initializes the variance vector for each cluster in the variational inference algorithm. It returns a vector of variance vectors, v_sq_k_hat_init, based on the provided number of clusters, number of genes, and initialization type.
"""
function init_v_sq_k_hat_vec!(v_sq_k_hat_init,K,G;rand_init = false, lo=0,hi=1)
    v_sq_k_vec = ones(G)
    if isnothing(v_sq_k_hat_init) && rand_init
        v_sq_k_hat_init =  [rand(Uniform(lo,hi),length(v_sq_k_vec)) for k in 1:K]
    elseif isnothing(v_sq_k_hat_init) && !rand_init
        v_sq_k_hat_init =  [v_sq_k_vec for k in 1:K] #
    end 
    return v_sq_k_hat_init
end

"""
    init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = false, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2)
This function initializes the gamma and eta vectors for each cluster in the variational inference algorithm. It returns a vector of gamma values, gk_hat_init, and a vector of eta values, hk_hat_init, based on the provided number of clusters and initialization type.
"""
function init_ghk_hat_vec!(gk_hat_init,hk_hat_init,K;rand_init = false, g_lo=0,g_hi=1, h_lo= 0,h_hi = 2)
    if isnothing(gk_hat_init) || isnothing(hk_hat_init)
        if rand_init
            gk_hat_init = rand(Uniform(g_lo,g_hi), (K,));
            hk_hat_init = rand(Uniform(h_lo,h_hi), (K,));
        else
            gk_hat_init, hk_hat_init = init_params_states(K)
        end
    end
    return gk_hat_init,hk_hat_init
end

"""
    init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = false)
This function initializes the conditional probability vector for each time point in the variational inference algorithm. It returns a vector of conditional probability vectors, c_ttprime_init, based on the provided number of time points and initialization type.
"""
function init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = false)
    if isnothing(c_ttprime_init) && rand_init
        c_ttprime_init = [rand(Dirichlet(ones(T) ./T)) for t in 1:T]
    elseif isnothing(c_ttprime_init) && !rand_init
        c_ttprime_init = [ones(T) ./T  for t in 1:T]
    end
    
    return c_ttprime_init
end

"""
    init_d_hat_vec!(d_hat_init,K,T;rand_init = false,uniform_theta_init=false, gk_hat_init = nothing, hk_hat_init= nothing)
This function initializes the d_hat variable for the variational inference algorithm. If d_hat_init is not provided, it can be initialized randomly, uniformly, or based on the gk_hat and hk_hat parameters depending on the flags provided.
"""
function init_d_hat_vec!(d_hat_init,K,T;rand_init = false,uniform_theta_init=false, gk_hat_init = nothing, hk_hat_init= nothing)
    if isnothing(d_hat_init)
        if uniform_theta_init
            d_hat_init = [ones(K+1) ./(K+1)  for t in 1:T]#
        else
            if rand_init
                d_hat_init = [rand(K+1) for t in 1:T]
            else
                d_hat_init = init_d_hat_tk(T,gk_hat_init, hk_hat_init);
            end
        end
    end

    return d_hat_init
end

function init_d_hat_tk(T,g_hat_vec, h_hat_vec)
    d_hat_vec = [βk_expected_value(g_hat_vec, h_hat_vec) for t in 1:T]
    return d_hat_vec
end

"""
    init_yjk_vec!(yjk_init,G,K;rand_init = false)
This function initializes the yjk variable for the variational inference algorithm. If yjk_init is not provided, it can be initialized randomly or uniformly depending on the flags provided.
"""
function init_yjk_vec!(yjk_init,G,K;rand_init = false)
    if isnothing(yjk_init) && rand_init
        yjk_init = [[rand(Beta(1.,1.)) for j in 1:G] for k in 1:K] 
    elseif isnothing(yjk_init) && !rand_init
        yjk_init =[[0.5 for j in 1:G] for k in 1:K] 
    end

    return yjk_init
end

"""
    init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
This function initializes the st_hat variable for the variational inference algorithm. If st_hat_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_st_hat_vec!(st_hat_init,T,ϕ0;rand_init = false, lo=0,hi=1)
    if isnothing(st_hat_init) && rand_init
        st_hat_init = [rand(Uniform(lo,hi)) for t in 1:T]
    elseif isnothing(st_hat_init) && !rand_init
        st_hat_init = [ϕ0 for t in 1:T]
    end
    return st_hat_init
end

"""
    init_rtik_vec!(rtik_init,K,T,N_t;rand_init = false)
This function initializes the rtik variable for the variational inference algorithm. If rtik_init is not provided, it can be initialized randomly or uniformly based on the rand_init flag.
"""
function init_rtik_vec!(rtik_init,K,T,N_t;rand_init = false)
    if isnothing(rtik_init) && rand_init
        rtik_init = [[rand(Dirichlet(ones(K) ./K)) for i in 1:N_t[t]] for t in 1:T]
    elseif  isnothing(rtik_init) && !rand_init
        rtik_init = [[ones(K) ./K for i in 1:N_t[t]] for t in 1:T]
    end
    

    return rtik_init
end

#####################################################
#####################################################
################# FAST FUNCTIONS ####################
#####################################################
#####################################################
#####################
#####################

"""
    initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init)
This function initializes the types for the variational inference algorithm. It sets the initial values for various parameters in the cell population, clusters, condition parameters, data parameters, model parameters, and gene parameters based on the provided initial values.
"""
function initialize_VariationalInference_types!(cellpop,clusters,conditionparams,dataparams,modelparams,geneparams,mk_hat_init,v_sq_k_hat_init,λ_sq_init,σ_sq_k_init,gk_hat_init,hk_hat_init,d_hat_init,rtik_init,yjk_init,c_ttprime_init,st_hat_init)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    N_t = dataparams.N_t
    K = modelparams.K
    for j in 1:G
        geneparams[j].λ_sq[1] = λ_sq_init[j]
    end
    for k in 1:K
        clusters[k].mk_hat .= copy(mk_hat_init[k])
        clusters[k].v_sq_k_hat .= copy(v_sq_k_hat_init[k])
        clusters[k].σ_sq_k_hat .= copy(σ_sq_k_init[k])
        clusters[k].var_muk .= copy(yjk_init[k] .* (mk_hat_init[k].^2 .+ v_sq_k_hat_init[k]) .- yjk_init[k] .* mk_hat_init[k].^2)
        clusters[k].κk_hat .= copy(yjk_init[k] .* mk_hat_init[k])
        clusters[k].yjk_hat .= copy(yjk_init[k])
        clusters[k].gk_hat[1] = gk_hat_init[k] 
        clusters[k].hk_hat[1] = hk_hat_init[k]
        clusters[k].ak_hat[1] = StatsFuns.logit(gk_hat_init[k])
        clusters[k].bk_hat[1] = log(hk_hat_init[k])
    end
    n = 0
    for t in 1:T
        for i in 1:N_t[t]
        n += 1
        cellpop[n].rtik .= copy(rtik_init[t][i])
        end
        conditionparams[t].c_tt_prime .= copy(c_ttprime_init[t])
        conditionparams[t].d_hat_t .= copy(d_hat_init[t])
        conditionparams[t].d_hat_t_sum[1] = sum(d_hat_init[t])
        if t ==T
            conditionparams[t].st_hat[1] = 0.0
        else
            conditionparams[t].st_hat[1] = st_hat_init[t]
        end
    end

    return cellpop,clusters,conditionparams,dataparams,modelparams,geneparams
end
