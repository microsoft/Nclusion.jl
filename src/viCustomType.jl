"""
    check_nothing_type(val)
This function checks if a value is `nothing` or an array of `Nothing` type.
"""
function check_nothing_type(val)
    is_nothing_bool = isnothing(val)
    is_nothing_array_bool  = val isa AbstractArray && eltype(val) <: Nothing
    if is_nothing_bool || is_nothing_array_bool
        return true
    else
        return false
    end
end

"""
    Features
This is an abstract type for all Features tracked by NCLUSION
"""
abstract type Features end 


"""
    CellFeature
This is type that allows NCLUSION to track all of the cell-specific features during inference
"""
struct CellFeature{U <: AbstractFloat, W <:Int64, J} <: Features 
    i::W
    t::W
    n::W
    r::Vector{U}
    c::Vector{U}
    z_argmax::Vector{W}
    cache::Vector{U}
    x::NTuple{J, U}
    xsq::NTuple{J, U}
    BitType::DataType
    function CellFeature(i,t,n,K,T,data;rand_init = false,c_init = nothing,r_init = nothing)
        float_type = eltype(data)
        U = eltype(data)
        J = length(data)
        W = typeof(t)
        Kplus = K+1
        data_sq = data .^2
        numel = 0
        if typeof(T) <: AbstractArray
            numel = maximum([numel,Kplus,J,recursive_flatten(T)...])
            if check_nothing_type(c_init) && rand_init
                c_init = [rand(Dirichlet(ones(t) ./t));zeros(T[i]-t)]
            elseif check_nothing_type(c_init) && !rand_init
                c_init = [ones(t) ./t;zeros(T[i]-t)]
            end
        else
            numel = maximum([numel,Kplus,J,T])
            if check_nothing_type(c_init) && rand_init
                c_init = [rand(Dirichlet(ones(t) ./t));zeros(T-t)]
            elseif check_nothing_type(c_init) && !rand_init
                c_init = [ones(t) ./t;zeros(T-t)]
            end
        end
        if check_nothing_type(r_init) && rand_init
            r_init = [rand(Dirichlet(ones(K) ./K));zeros(1) ]
        elseif  check_nothing_type(r_init) && !rand_init
            r_init = [ones(K) ./K ;zeros(1)  ]
        end
        z_argmax = [Int(argmax(r_init))]
        new{U,W,J}(i,t,n,r_init,c_init,z_argmax,zeros(U,numel),Tuple(data),Tuple(data_sq),U)
    end
end

"""
    MatrixConditionFeature
This is type that allows NCLUSION to track all of the matrix condition-specific features during inference
"""
struct MatrixConditionFeature{U <: AbstractFloat, W <:Int64} <: Features#cluster_features
    i::W
    t::W
    CNtk::Matrix{U}
    suffstats_cache::Matrix{U}
    cache::Vector{U}
    condition_update_neighbors::Vector{W}
    condition_network_neighbors::Vector{W}
    BitType::DataType
    function MatrixConditionFeature(i,t,K,T,condition_update_neighbors,condition_network_neighbors;float_type=Float64)
        U = float_type
        W = typeof(t)
        Kplus = K+1
        # numel = maximum([Kplus,T])
        numel = 0
        if typeof(T) <: AbstractArray
            numel = maximum([numel,Kplus,recursive_flatten(T)...])
            T_used = T[i]
        else
            numel = maximum([numel,Kplus,T])
            T_used = T
        end
        new{U,W}(i,t,zeros(U,T_used,Kplus),zeros(U,T_used,Kplus),zeros(U,numel),condition_update_neighbors,condition_network_neighbors,U)
    end
end
"""
    ClusterFeature
This is type that allows NCLUSION to track all of the cluster-specific features during inference
"""
struct ClusterFeature{U <: AbstractFloat, W <:Int64} <: Features #cluster_features
    k::W
    m_mu::Vector{U}
    s_sq_mu::Vector{U}
    m_nu::Vector{U}
    s_sq_nu::Vector{U}
    y::Vector{U}
    Nk::Vector{U}#size 1 so we can mutate it
    x_hat::Vector{U}
    x_hat_sq::Vector{U}
    a::Vector{U}
    b::Vector{U}
    u::Vector{U}
    v::Vector{U}
    g1::Vector{U}#size 1 so we can mutate it
    g2::Vector{U}#size 1 so we can mutate it
    h1::Vector{U}#size 1 so we can mutate it
    h2::Vector{U}#size 1 so we can mutate it
    alpha_Tk::Vector{U}#size 1 so we can mutate it
    cache::Vector{U}# accumulation::Vector{U}
    BitType::DataType
    function ClusterFeature(k,J;float_type=Float64,m_mu_init = nothing,s_sq_mu_init = nothing,m_nu_init = nothing,s_sq_nu_init = nothing,y_init = nothing,g1_init = nothing,g2_init = nothing,h1_init = nothing,h2_init = nothing,a_init = nothing, b_init = nothing,u_init = nothing, v_init = nothing,rand_init = false,update_clusterwise=false,eta_update_mode="Local",lambda_update_mode="Local",sigma_update_mode="Local")
        U = float_type
        W = typeof(k)
        if check_nothing_type(g1_init) && rand_init
            g1_init = logistic.(randn(float_type,1))
        elseif check_nothing_type(g1_init) && !rand_init
            g1_init = ones(float_type,1)
        end
        if check_nothing_type(g2_init) && rand_init
            g2_init = exp.(randn(float_type,1))
        elseif check_nothing_type(g2_init) && !rand_init
            g2_init = ones(float_type,1)
        end
        if check_nothing_type(m_mu_init) && rand_init
            m_mu_init = randn(float_type,J) 
        elseif check_nothing_type(m_mu_init) && !rand_init
            m_mu_init = zeros(float_type,J)
        end
        if check_nothing_type(s_sq_mu_init) && rand_init
            s_sq_mu_init = exp.(randn(float_type,J)) 
        elseif check_nothing_type(s_sq_mu_init) && !rand_init
            s_sq_mu_init = ones(float_type,J)
        end
        if check_nothing_type(m_nu_init) 
            m_nu_init = zeros(float_type,J)
        end
        if check_nothing_type(s_sq_nu_init)
            s_sq_nu_init = ones(float_type,J)
        end
        if check_nothing_type(y_init) && rand_init
            y_init = logistic.(randn(float_type,J))
        elseif check_nothing_type(y_init) && !rand_init
            y_init = 0.5*ones(float_type,J)
        end
        a, b = 0.0, 1.0
        eta_update_mode_bool = (eta_update_mode == "Local") || (eta_update_mode == "Genewise")
        if check_nothing_type(h1_init) && rand_init && eta_update_mode_bool
            h1_init = rand(float_type,J) #exp.(randn(K))
        elseif check_nothing_type(h1_init) && rand_init && !eta_update_mode_bool
            h1_init = rand(float_type,1) .* ones(float_type,J) 
        elseif check_nothing_type(h1_init) && !rand_init
            h1_init = ones(float_type,J)
        end
        if check_nothing_type(h2_init) && rand_init && eta_update_mode_bool
            h2_init = rand(float_type,J) #exp.(randn(K))
        elseif check_nothing_type(h2_init) && rand_init && !eta_update_mode_bool
            h2_init = rand(float_type,1) .* ones(float_type,J)
        elseif check_nothing_type(h2_init) && !rand_init
            h2_init = ones(float_type,J)
        end
        lambda_update_mode_bool = (lambda_update_mode == "Local") || (lambda_update_mode == "Genewise")
        if check_nothing_type(u_init) && rand_init && lambda_update_mode_bool
            u_init = exp.(randn(float_type,J))
        elseif check_nothing_type(u_init) && rand_init && !lambda_update_mode_bool
            u_init = rand(float_type,1) .* ones(float_type,J)
        elseif (check_nothing_type(u_init) && !rand_init) #|| !update_clusterwise
            u_init = ones(float_type,J)
        end
        if check_nothing_type(v_init) && rand_init && lambda_update_mode_bool
            v_init = exp.(randn(float_type,J))
        elseif check_nothing_type(v_init) && rand_init && !lambda_update_mode_bool
            v_init = rand(float_type,1) .* ones(float_type,J) 
        elseif (check_nothing_type(v_init) && !rand_init) #|| !update_clusterwise
            v_init = ones(float_type,J)
        end
        sigma_update_mode_bool = (sigma_update_mode == "Local") || (sigma_update_mode == "Genewise")
        if check_nothing_type(a_init) && rand_init && sigma_update_mode_bool
            a_init = exp.(randn(float_type,J)) 
        elseif check_nothing_type(a_init) && rand_init && !sigma_update_mode_bool
            a_init = rand(float_type,1) .* ones(float_type,J) 
        elseif (check_nothing_type(a_init) && !rand_init) #|| !update_clusterwise
            a_init = ones(float_type,J)
        end
        if check_nothing_type(b_init) && rand_init && sigma_update_mode_bool
            b_init = exp.(randn(float_type,J)) 
        elseif check_nothing_type(b_init) && rand_init && !sigma_update_mode_bool
            b_init = rand(float_type,1) .* ones(float_type,J) 
        elseif (check_nothing_type(b_init) && !rand_init) #|| !update_clusterwise
            b_init = ones(float_type,J)
        end
        new{U,W}(k,m_mu_init,s_sq_mu_init,m_nu_init,s_sq_nu_init,y_init,zeros(U,1),zeros(U,J),zeros(U,J),a_init, b_init,u_init,v_init,g1_init,g2_init,h1_init,h2_init,zeros(U,1),zeros(U,J),U)
    end
end

"""
    GeneFeatures
This is type that allows NCLUSION to track all of the gene specific features during inference
"""
struct GeneFeatures{U <: AbstractFloat, W <:Int64,P <: Function} <: Features #cluster_features
    j::W
    λ_sq::Vector{U}#size 1 so we can mutate it
    cache::Vector{U}# accumulation::Vector{U}
    _reset!::P
    BitType::DataType
    function GeneFeatures(j;float_type=Float64)
        U = float_type
        W = typeof(j)
        function _reset!(val_vec::Vector{U},BitType::DataType) where U <: AbstractFloat
            for indx in eachindex(val_vec)
                val_vec[indx] = zero(BitType)
            end
            return val_vec
        end
        P = typeof(_reset!)
        new{U,W,P}(j,Vector{U}(undef,1),zeros(U,1),_reset!,U)
    end
end

"""
    ConditionFeature
This is type that allows NCLUSION to track all of the condition specific features during inference
"""
struct ConditionFeature{U <: AbstractFloat, W <:Int64} <: Features#cluster_features
    i::W
    t::W
    d_sum::Vector{U}#size 1 so we can mutate it
    d::Vector{U}
    Ctt::Vector{U}
    w1::Vector{U}#size 1 so we can mutate it
    w2::Vector{U}#size 1 so we can mutate it
    condition_update_neighbors::Vector{W}
    condition_network_neighbors::Vector{W}
    BitType::DataType
    function ConditionFeature(i,t,K,T,condition_update_neighbors,condition_network_neighbors;float_type=Float64,d_init = nothing,w1_init = nothing,w2_init = nothing,rand_init = false)
        U = float_type
        W = typeof(t)
        Kplus = K+1
        if check_nothing_type(d_init) && rand_init
            d_init = rand(Float64,Kplus) 
        elseif check_nothing_type(d_init) && !rand_init
            d_init =ones(Float64,Kplus) ./(Kplus) 
        end
        if t ==1
            w1_init = ones(float_type,1)
            w2_init = ones(float_type,1)
        else
            if check_nothing_type(w1_init) && rand_init
                w1_init = exp.(randn(float_type,1))
            elseif check_nothing_type(w1_init) && !rand_init
                w1_init = ones(float_type,1)
            end
            if check_nothing_type(w2_init) && rand_init
                w2_init = exp.(randn(float_type,1))
            elseif check_nothing_type(w2_init) && !rand_init
                w2_init = ones(float_type,1)
            end
        end
        numel = 0
        if typeof(T) <: AbstractArray
            numel = T[i]
        else
            numel = T
        end
        new{U,W}(i,t,zeros(U,1),d_init,zeros(U,numel),w1_init,w2_init,condition_update_neighbors,condition_network_neighbors,U)
    end
end

"""
    DataFeature
This is type that allows NCLUSION to track all of other dataset-specific features during inference
"""
struct DataFeature{U <: AbstractFloat,W <: Int64} <: Features#cluster_features
    I::W
    T::Vector{W}
    J::W
    N_t::Vector{Vector{W}}
    N::W
    LinearAddress::Vector{Tuple{W, W,W}}
    TimeRanges::Vector{Tuple{W, W}}
    Jlog::U
    logpi::U
    BitType::DataType
    function DataFeature(data_input)
        U = eltype(data_input[1][1][1])
        I = length(data_input)
        T = [length(data_input[i]) for i in 1:I]
        N_t = [[length(data_input[i][t]) for t in  1:T[i]] for  i in 1:I]
        N = convert(U,sum([sum(el) for el in N_t]))
        J =length(data_input[1][1][1])
        Jlog = convert(U,J*log(2π))
        logpi = convert(U,log(2π))
        W = typeof(I)
        LinearAddress = [(i,t,n) for i in 1:I for t in 1:T[i] for n in 1:N_t[i][t]]
        TimeRanges = get_timeranges(N_t)
        new{U,W}(I,T,J,N_t,N,LinearAddress,TimeRanges,Jlog,logpi,U)
    end
end

"""
    ElboFeatures
This is type that allows NCLUSION to track all the elbo during inference
"""
struct ElboFeatures{U <: AbstractFloat, W <:Int64} <: Features#cluster_features
    l::W
    elbo_::Vector{Union{Missing,U}}#size 1 so we can mutate it
    per_k_elbo::Matrix{Union{Missing,U}}
    BitType::DataType
    function ElboFeatures(l,K,num_iter;float_type=Float64)
        U = float_type
        W = typeof(l)
        new{U,W}(l,Vector{Union{Missing,Float64}}(undef,num_iter),Matrix{Union{Missing,Float64}}(undef,K,num_iter),U)
    end
end


"""
        ModelParameterFeature
    This is type that allows NCLUSION to track all of the user defined model parameters during inference
"""
struct ModelParameterFeature{U <: AbstractFloat,V <:AbstractFloat,W <: Int64} <: Features #cluster_features
    K::W
    alpha0::Vector{Vector{U}}#size is equal to the number of individuals and the number of time points within each individual
    gamma0::Vector{U}#size 1 so we can mutate it
    phi1::Vector{U}#size 1 so we can mutate it
    phi2::Vector{U}#size 1 so we can mutate it
    kappa1::Vector{U}#size 1 so we can mutate it
    kappa2::Vector{U}#size 1 so we can mutate it
    xi1::Vector{U}#size 1 so we can mutate it
    xi2::Vector{U}#size 1 so we can mutate it
    varphi1::Vector{U}#size 1 so we can mutate it
    varphi2::Vector{U}#size 1 so we can mutate it
    nu0::Vector{U}#size 1 so we can mutate it
    sigma_sq_nu::Vector{U}#size 1 so we can mutate it
    significance_prop::U
    min_number_cells::U
    min_percent_cells::U
    min_percent_of_genes::U
    max_percent_of_genes::U
    num_iter::W
    uniform_theta_init::Bool
    rand_init::Bool
    change_seeds::Bool
    init_seed::W
    BitType::DataType
    # ep::U
    # elbo_ep::U
    function ModelParameterFeature(data_input,K,alpha0,gamma0,phi1,phi2,kappa1,kappa2,xi1,xi2,varphi1,varphi2,nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,num_iter,uniform_theta_init,rand_init,change_seeds,seed)
        U = eltype(data_input[1][1][1])
        V = eltype(data_input[1][1][1])
        W = typeof(K)
        I = length(data_input)
        T = [length(data_input[i]) for i in 1:I]
        J = length(data_input[1][1][1])
        if typeof(alpha0) <: Number
            alpha0 = [[alpha0 for t in 1:T[i]] for i in 1:I]
        end
        if typeof(nu0) <: Number
            nu0 = [nu0 for j in 1:J]
        end
        if typeof(sigma_sq_nu) <: Number
            sigma_sq_nu = [sigma_sq_nu for j in 1:J]
        end
        new{U,V,W}(K,alpha0,[gamma0],[phi1],[phi2],[kappa1],[kappa2],[xi1],[xi2],[varphi1],[varphi2],nu0,sigma_sq_nu,significance_prop,min_number_cells,min_percent_cells,min_percent_of_genes,max_percent_of_genes,num_iter,uniform_theta_init,rand_init,change_seeds,seed,V)
    end
end


"""
    TrainFeature
This is type that allows NCLUSION to track all the feature changes during inference
"""
struct TrainFeature{U <: AbstractFloat, W <:Int64} <: Features#cluster_features
    l::W
    elbo_::Vector{Union{Missing,U}}#size 1 so we can mutate it
    d::Vector{Union{Missing,Matrix{U}}}
    y::Vector{Union{Missing,Matrix{U}}}
    m_mu::Vector{Union{Missing,Matrix{U}}}
    s_sq_mu::Vector{Union{Missing,Matrix{U}}}
    h1::Vector{Union{Missing,Matrix{U}}}
    h2::Vector{Union{Missing,Matrix{U}}}
    Nk::Vector{Union{Missing,Vector{U}}}
    w1::Vector{Union{Missing,Vector{U}}}
    w2::Vector{Union{Missing,Vector{U}}}
    g1::Vector{Union{Missing,Vector{U}}}
    g2::Vector{Union{Missing,Vector{U}}}
    a::Vector{Union{Missing,Matrix{U}}}
    b::Vector{Union{Missing,Matrix{U}}}
    m_nu::Vector{Union{Missing,Vector{U}}}
    s_sq_nu::Vector{Union{Missing,Vector{U}}}
    u::Vector{Union{Missing,Matrix{U}}}
    v::Vector{Union{Missing,Matrix{U}}}
    Tcache::Vector{U}# accumulation::Vector{U}
    Kcache::Vector{U}# accumulation::Vector{U}
    Jcache::Vector{U}# accumulation::Vector{U}
    JKcache::Matrix{U}# accumulation::Vector{U}
    KplusTcache::Matrix{U}# accumulation::Vector{U}
    BitType::DataType
    function TrainFeature(l,T_all,K,J,num_iter;float_type=Float64)
        Kplus = K+1
        U = float_type
        W = typeof(l)
        elbo_init = Vector{Union{Missing,float_type}}(missing,num_iter)
        d_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        y_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        m_mu_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        s_sq_mu_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        h1_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        h2_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        Nk_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        w1_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        w2_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        g1_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        g2_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        a_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        b_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        m_nu_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        s_sq_nu_init = Vector{Union{Missing,Vector{float_type}}}(missing,num_iter)
        u_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        v_init = Vector{Union{Missing,Matrix{float_type}}}(missing,num_iter)
        Tcache = zeros(float_type,T_all)
        Kcache = zeros(float_type,K)
        Jcache = zeros(float_type,J)
        JKcache = zeros(float_type,J,K)
        KplusTcache = zeros(float_type,Kplus,T_all)
        new{U,W}(l,elbo_init,d_init,y_init,m_mu_init,s_sq_mu_init,h1_init,h2_init,Nk_init,w1_init,w2_init,g1_init,g2_init,a_init,b_init,m_nu_init,s_sq_nu_init,u_init,v_init,Tcache,Kcache,Jcache,JKcache,KplusTcache,U)
    end
end

""""
    get_timeranges(N_t)
This function returns the linear indices that contain cells from the same condition.
"""
function get_timeranges(N_t)
    I = length(N_t)
    T = [length(el) for el in N_t]
    T_all = 0
    for i in 1:I
        for t in 1:T[i]
            T_all += 1
        end
    end
    starts = Vector{Int}(undef,T_all)
    ends = Vector{Int}(undef,T_all)
    it = 1
    for i in 1:I    
        for t in 1:T[i]
            if it == 1
                starts[1] = 0 + 1
                ends[1] = 0 + N_t[1][1]
                it += 1
                continue
            end
            starts[it] = ends[it-1] + 1
            ends[it] = ends[it-1] + N_t[i][t]
            it += 1
        end
    end
    return [(st,en) for (st,en) in zip(starts,ends)]
end

"""
    get_linear_index_as_ragged_array(N_t)
This function returns the linear indices that contain cells from the same condition as a ragged array.
"""
function get_linear_index_as_ragged_array(N_t)
    I = length(N_t)
    T = [length(N_t[i]) for i in 1:I]
    N = convert(typeof(I),sum([sum(el) for el in N_t]))
    T_all = sum(T)
    linear_sample_index = collect(1:N)
    linear_time_index = collect(1:T_all)
    linear_sample_index_as_ragged_array = Vector{Vector{Vector{Int}}}(undef,I)
    linear_time_index_as_ragged_array = Vector{Vector{Int}}(undef,I)
    sample_counter=1
    time_counter=1
    for i in 1:I
        linear_sample_index_as_ragged_array[i] = Vector{Vector{Int}}(undef,T[i])
        linear_time_index_as_ragged_array[i] = Vector{Int}(undef,T[i])
        for t in 1:T[i]
            linear_sample_index_as_ragged_array[i][t] = Vector{Int}(undef,N_t[i][t])
            linear_time_index_as_ragged_array[i][t] = linear_time_index[time_counter]
            time_counter+=1
            for n in 1:N_t[i][t]
                linear_sample_index_as_ragged_array[i][t][n] = linear_sample_index[sample_counter]
                sample_counter+=1
            end
        end
    end
    return linear_sample_index_as_ragged_array, linear_time_index_as_ragged_array
end

"""
    get_linear_time_condition_update_neighbors(input;get_ragged_array=false)
This function returns the linear indices that contain cells from the same condition as a ragged array.
"""
function get_linear_time_condition_update_neighbors(input;get_ragged_array=false)
    I = length(input)
    T = [length(el) for el in input]
    if get_ragged_array
        _ , linear_time_index_as_ragged_array = get_linear_index_as_ragged_array(input)
    else
        linear_time_index_as_ragged_array = input
    end
    condition_update_neighbors = Vector{Vector{Vector{Int}}}(undef,I)
    for i in 1:I
        condition_update_neighbors[i] = Vector{Vector{Int}}(undef,T[i])
        for t in 1:T[i]
            if t != 1 
                condition_update_neighbors[i][t] = linear_time_index_as_ragged_array[i][t:T[i]]
            else 
                condition_update_neighbors[i][t] = [0]
            end
        end
    end
    return condition_update_neighbors
end

"""
    get_linear_time_condition_network_neighbors(input;get_ragged_array=false)
This function returns the linear indices that contain cells from the same condition as a ragged array.
"""
function get_linear_time_condition_network_neighbors(input;get_ragged_array=false)
    I = length(input)
    T = [length(el) for el in input]
    if get_ragged_array
        _ , linear_time_index_as_ragged_array = get_linear_index_as_ragged_array(input)
    else
        linear_time_index_as_ragged_array = input
    end
    condition_network_neighbors = Vector{Vector{Vector{Int}}}(undef,I)
    for i in 1:I
        condition_network_neighbors[i] = Vector{Vector{Int}}(undef,T[i])
        for t in 1:T[i]
            condition_network_neighbors[i][t] = linear_time_index_as_ragged_array[i][t:T[i]]
        end
    end
    return condition_network_neighbors
end


"""
    _reset!(val_vec::Vector{U},BitType::DataType)
Performs an inplace setting of values in a vector to 0. Maintains the type of the variable prior to reset.
"""
function _reset!(val_vec::Vector{U},BitType::DataType) where U <: AbstractFloat
    for indx in eachindex(val_vec)
        val_vec[indx] = zero(BitType)
    end
    return val_vec
end