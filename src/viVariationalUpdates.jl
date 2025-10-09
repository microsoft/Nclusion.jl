
############################################
############# UPDATE FUNCTIONS #############
############################################

"""
    update_w1!(conditions,dataparams,modelparams)
This function updates the w1 parameter in the ConditionFeature type. This is one of the two a varational parameters used to to appoximate the omega parameter in the model. The omega parameter follows a Beta distribution in the true posterior. The omega parameter follows a Beta distribution in the variational approximation.  The update is done in place.
```math
    p(omega_{t_{i}} | phi_1, phi_2) = text{Beta}(omega_{t} | phi_1, phi_2) 

    q(omega_{t_{i}} | hat{w}_{1t_{i}}, hat{w}_{2t_{i}}) = text{Beta}(omega_{t_{i}} | hat{w}_{1t_{i}}, hat{w}_{2t_{i}})

    hat{w}_{1t_{i}} rightarrow  phi_1 +  left[ sum_{t_{i}'=t_{i}}^{T_{i}}C^{ left(t_{i}' right)}_{t_{i}} right]
```
"""
function update_w1!(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <:Int64} #  formerly part of update_st_hat_mpu!
    float_type = dataparams.BitType
    T = dataparams.T
    I = dataparams.I
    T_all = sum(T)
    for it in 1:T_all
        condition_i = conditions[it].i
        condition_t = conditions[it].t
        if condition_t == 1
            conditions[it].w1[1]=1.0
            continue
        end
        w1 = 0.0
        w1 += modelparams.phi1[1]
        for t_prime in conditions[it].condition_update_neighbors
            w1 += conditions[t_prime].Ctt[condition_t]
        end
        conditions[it].w1[1] = w1
    end
    return conditions
end

"""
    update_w2!(conditions,dataparams,modelparams)
This function updates the w2 parameter in the ConditionFeature type. This is one of the two a varational parameters used to to appoximate the omega parameter in the model. The omega parameter follows a Beta distribution in the true posterior. The omega parameter follows a Beta distribution in the variational approximation.  The update is done in place.
```math
    p(omega_{t_{i}} | phi_1, phi_2) = text{Beta}(omega_{t} | phi_1, phi_2) 

    q(omega_{t_{i}} | hat{w}_{1t_{i}}, hat{w}_{2t_{i}}) = text{Beta}(omega_{t_{i}} | hat{w}_{1t_{i}}, hat{w}_{2t_{i}})

    hat{w}_{2t_{i}} rightarrow phi_2 +  left[ sum_{t_{i}'=t_{i}}^{T_{i}} sum_{m=1}^{t_{i}-1}C^{ left(t_{i}' right)}_{m} right]
```
    
"""
function update_w2!(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <:Int64} #  formerly part of update_st_hat_mpu!
    float_type = dataparams.BitType
    T = dataparams.T
    I = dataparams.I
    T_all = sum(T)
    full_condition_network_neighbors = nothing
    for it in 1:T_all
        condition_i = conditions[it].i
        condition_t = conditions[it].t
        if condition_t == 1
            conditions[it].w2[1]=1.0
            # full_condition_network_neighbors = conditions[it].condition_network_neighbors
            continue
        end
        w2 = 0.0
        w2 += modelparams.phi2[1]
        for t_prime in conditions[it].condition_update_neighbors
            for m in 1:condition_t-1
                w2 += conditions[t_prime].Ctt[m]
            end
        end
        conditions[it].w2[1] = w2
    end
    return conditions
    # conditions[1].w2[1]=1.0
    # for t in 2:T
    #     w2 = 0.0
    #     w2 += modelparams.phi2[1]
    #     for t_prime in t:T
    #         for m in 1:t-1
    #             w2 += conditions[t_prime].Ctt[m]
    #         end
    #     end
    #     conditions[t].w2[1] = w2
    # end
    # return conditions
end


"""
    update_d!(conditions,dataparams,modelparams)
This function updates the d parameter in the ConditionFeature type. This is a varational parameter used to to appoximate the pi parameter in the model. The pi parameter follows a Dirichlet distribution in the true posterior. The pi parameter follows a Dirichlet distribution in the variational approximation.  The update is done in place.
```math
    p(left { hat{pi}^{(t_{i})}_{k} right }_{k=1}^K  | left { alpha_{t_{i}} * textbf{SB}(chi_k)  right }_{k=1}^{K} ) = text{Dirichlet} left(left { hat{pi}^{(t_{i})}_{k} right }_{k=1}^K  | left { alpha_{t_{i}} * textbf{SB}(chi_k)  right }_{k=1}^{K} right)

    q(left { hat{pi}^{(t_{i})}_{k} right }_{k=1}^K  | left { hat{d}^{(t_{i})}_{k} right }_{k=1}^K ) = text{Dirichlet} left(left { hat{pi}^{(t_{i})}_{k} right }_{k=1}^K  | left { hat{d}^{(t_{i})}_{k} right }_{k=1}^K right)

    prod_{i=1}^{I} prod_{t_{i}=1}^{T_{i}}  text{Dirichlet} left( left { hat{d}^{(t_{i})}_{k} right }_{k=1}^K right) quad  text{s.t.}

    hat{d}^{(t_{i})}_{k} rightarrow  alpha_{t_{i}} times hat{g}_{1k} times prod_{ ell=1}^{k-1}(1- hat{g}_{1 ell}) + sum_{t_{i}'=t_{i}}^{T_{i}} sum_{n=1}^{N_{t_{i}'}} hat{r}_{nk}^{ left(t_{i}' right)} hat{c}_{nt_{i}}^{ left(t_{i}' right)}
```
    
"""
function update_d!(clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;use_log=false) where {U <: AbstractFloat, W <: Int64} #  formerly update_d_hat_mpu!
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    T_all = sum(T)
    @inbounds for it in 1:T_all
        # _reset!(conditionparams[t].d_hat_t, float_type)
        conditions_i = conditions[it].i
        conditions_t = conditions[it].t
        conditions[it].d .= 0.0
        matrixconditions[it].suffstats_cache .= 0.0
        # @inbounds for tt in t:T[i]
        @inbounds for tt in conditions[it].condition_network_neighbors
            matrixconditions[it].suffstats_cache .+= matrixconditions[tt].CNtk
        end
        @inbounds for k in 1:Kplus
            E_SBk = expectation_SBk(k,clusters,modelparams;use_log = use_log)
            conditions[it].d[k] = modelparams.alpha0[conditions_i][conditions_t] * E_SBk + matrixconditions[it].suffstats_cache[conditions_t,k]
        end
        # conditionparams[t].d_hat_t[k] = modelparams.α0 * e_βk + updated_Ntk_sum
    end
    return conditions
end


"""
     update_d_sum!(conditions::Vector{ConditionFeature{U,W}}, dataparams::DataFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_d_sum!(conditions::Vector{ConditionFeature{U,W}}, dataparams::DataFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    T = dataparams.T
    T_all = sum(T)
    for it in 1:T_all
        conditions[it].d_sum[1] = sum(conditions[it].d)
    end
    return conditions
end


"""
    update_r!(cells,clusters,conditions,dataparams,modelparams)
This function updates the r parameter in the CellFeature type. This is a varational parameter used to to appoximate the psi parameter in the model. The psi parameter follows a Categorical distribution in the true posterior. The psi parameter follows a Categorical distribution in the variational approximation.  The update is done in place.
```math
    prod_{k=1}^{K} p( [psi_{t_{i}} = k] | { pi^{(t_{i})}_k }_{k=1}^{K}, tau_{n}^{(t_{i})}) = text{Categorical}([psi_{t_{i}} = k] | { pi^{tau_{n}^{(t_{i})}}_k }_{k=1}^{K}, )

    prod_{k=1}^{K} q( [psi_{t_{i}} = k] | { hat{r}^{(t_{i})}_{nk} }_{k=1}^{K}) = text{Categorical}([psi_{t_{i}} = k] | { hat{r}^{(t_{i})}_{nk} }_{k=1}^{K})

    =prod_{i=1}^{I} prod_{t_{i}=1}^{T_{i}} prod_{n=1}^{N_{t_{i}}} text{Categorical} left( left { hat{r}^{(t_{i})}_{nk} right }_{k=1}^K right) quad  text{s.t.}

    hat{r}^{(t_{i})}_{nk} rightarrow frac{ exp left( ln tilde{r}^{(t_{i})}_{nk} right)}{ sum_{k'=1}^{K} exp left( ln tilde{r}^{(t_{i})}_{nk'} right)}
    
    ln tilde{r}^{(t_{i})}_{nk} rightarrow sum_{t_{i}'=1}^{t_{i}} hat{c}_{nt_{i}'}^{(t_{i})} left[ boldsymbol{ Psi} left( hat{d}^{ left(t_{i}' right)}_k right) -  boldsymbol{ Psi} left( sum_{k'=1}^{K} hat{d}^{ left(t_{i}' right)}_{k'} right) right] - frac{J}{2} ln 2 pi  -  frac{1}{2} sum_{j=1}^{J} left( ln hat{b}_j -  boldsymbol{ Psi}( hat{a}_j)  right) - frac{1}{2} sum_{j=1}^{J} left[ frac{ hat{a}_j}{ hat{b}_j } right] left[x^{(t_{i}) 2}_{nj}-2x^{(t_{i})}_{nj} hat{y}_{kj} hat{m}_{kj}+ hat{y}_{kj} hat{m}_{kj}^2 + hat{y}_{kj} hat{s}_{kj}^2  right]
```
"""
function update_r!(cells, clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #  formerly update_rtik_mpu!
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    # @inbounds for n in 1:N
    Threads.@threads for n in 1:N
        # t = first(dataparams.LinearAddress[n])
        # cells[n]._reset!(cells[n].cache,cells[n].BitType)
        cells[n].cache .= 0.0
        adjust_E_ln_pi!(cells[n],conditions,dataparams)
        # cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
        for k in 1:K
            # cells[n].cache .= 0.0
            E_log_normal_l_j!(cells[n],clusters[k], dataparams)
            @inbounds cells[n].r[k] = cells[n].cache[k]
        end
        norm_weights3!(K,cells[n].r)
        cells[n].z_argmax[1] = argmax(cells[n].r)
    end
    return cells
end
"""
    adjust_E_ln_pi!(cell::CellFeature,conditions::Vector{ConditionFeature{U,W}}, dataparams::DataFeature) where {U <: AbstractFloat, W <: Int64}
"""
function adjust_E_ln_pi!(cell::CellFeature,conditions::Vector{ConditionFeature{U,W}}, dataparams::DataFeature) where {U <: AbstractFloat, W <: Int64}# formerly adjust_e_log_π_tk3!
    Kplus = length(conditions[1].d)
    t = cell.t
    n = cell.n
    i = cell.i
    T = dataparams.T
    # conditions[t].e_log_pi_t_cache .= 0.0
    for k in 1:Kplus
        conditions_pis_sums = 0.0
        @simd for tt in 1:t
            @fastmath @inbounds conditions_pis_sums += cell.c[tt] * E_ln_pi(k,conditions[sum(T[1:i-1])+tt])#(digamma(conditionparams[tt].d_hat_t[k]) - digamma(conditionparams[tt].d_hat_t_sum[1])) 
        end
        # conditions[t].e_log_pi_t_cache[k] = conditions_pis_sums
        cell.cache[k] = conditions_pis_sums
    end
    return cell
end
"""
    E_log_normal_l_j!(cellfeature::CellFeature,clusterfeature::ClusterFeature, dataparams::DataFeature)
"""
function E_log_normal_l_j!(cellfeature::CellFeature,clusterfeature::ClusterFeature, dataparams::DataFeature)# where {U <: AbstractFloat, W <: Int64} # formerly expectation_log_normal_l_j
    J = dataparams.J
    k = clusterfeature.k
    for j in 1:J
        @inbounds cellfeature.cache[k]  += - 0.5 * dataparams.logpi +  -0.5 * E_ln_sigma_sq(j,clusterfeature) - 0.5 * E_one_over_sigma_sq(j,clusterfeature) * E_ll_sq_diff_mu(j,cellfeature,clusterfeature) # 1/J * ()
    end
    return cellfeature
end


"""
    update_c!(cells,clusters,conditions,dataparams,modelparams)
This function updates the c parameter in the CellFeature type. This is a varational parameter used to to appoximate the tau parameter in the model. The tau parameter follows a Categorical distribution in the true posterior. The tau parameter follows a Categorical distribution in the variational approximation.  The update is done in place.
```math
    prod_{t_{i}'=1}^{t_{i}}  p( [tau_{n}^{(t_{i})} = t_{i}'] | textbf{SB}(omega_{t_{i}'}) ) = text{Categorical}([tau_{n}^{(t_{i})} = t_{i}'] | { textbf{SB}(omega_{t_{i}'}) }_{t_{i}'=1}^{t_{i}}, )

    prod_{t_{i}'=1}^{t_{i}} q( [tau_{n}^{(t_{i})} = t_{i}'] | hat{c}^{t_{i}}_{nt_{i}'} ) = text{Categorical}([tau_{n}^{(t_{i})} = t_{i}'] | { hat{c}^{t_{i}}_{nt_{i}'} }_{t_{i}'=1}^{t_{i}}, )

    =prod_{i=1}^{I} prod_{t_{i}=1}^{T_{i}} prod_{n=1}^{N_{t_{i}}} text{Categorical} left( left{hat{c}^{(t_{i})}_{nt_{i}'}right}_{t_{i}'=1}^{t_{i}}right)  quad   text{s.t.}
    hat{c}^{(t_{i})}_{nt_{i}'} rightarrow frac{ exp left( ln tilde{c}^{(t_{i})}_{nt_{i}'} right)}{ sum_{ ell=1}^{t_{i}} exp left( ln tilde{c}^{(t_{i})}_{n ell} right)}
    ln tilde{c}^{(t_{i})}_{nt_{i}'} rightarrow sum_{k=1}^{K} hat{r}_{nk}^{(t_{i})} left[ boldsymbol{ Psi} left( hat{d}^{ left(t_{i}' right)}_k right) -  boldsymbol{ Psi} left( sum_{k'=1}^{K} hat{d}^{ left(t_{i}' right)}_{k'} right) right]+ left( left[ boldsymbol{ Psi} left( hat{w}_{1t_{i}'} right) -  boldsymbol{ Psi} left( hat{w}_{1t_{i}'} +  hat{w}_{2t_{i}'} right) right] +  sum_{m=t_{i}'+1}^{t_{i}} left[ boldsymbol{ Psi} left( hat{w}_{2m} right) -  boldsymbol{ Psi} left( hat{w}_{1m} +  hat{w}_{2m} right) right]  right)
```
"""
function update_c!(cells, conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #  formerly update_c_ttprime_mpu!
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = modelparams.K + 1
    # @inbounds for n in 1:N Threads.@threads 
    for n in 1:N
        i = cells[n].i#dataparams.LinearAddress[n][1]
        t = cells[n].t#dataparams.LinearAddress[n][2]
        cells[n].c .= 0.0
        for tt in 1:t
            rd_sum = 0.0
            for k in 1:Kplus
                rd_sum += cells[n].r[k]*E_ln_pi(k,conditions[sum(T[1:i-1])+tt])
            end
            rd_sum += E_ln_omega(conditions[sum(T[1:i-1])+tt])
            ttplus = tt+1
            for m in ttplus:t
                rd_sum +=  E_ln_minusomega(conditions[sum(T[1:i-1])+m])
            end
            cells[n].c[tt] = rd_sum
        end
        norm_weights3!(t,cells[n].c)
    end
    return cells,conditions,dataparams,modelparams
end


"""
    update_y!(clusters,dataparams,modelparams)
This function updates the r parameter in the CellFeature type. This is a varational parameter used to to appoximate the psi parameter in the model. The psi parameter follows a Bernoulli distribution in the true posterior. The psi parameter follows a Bernoulli distribution in the variational approximation.  The update is done in place.
```math
    p( rho_{kj}| eta_k)= text{Bernoulli}( rho_{kj}| eta_k)   
    q( rho_{kj}|  hat{y}_{kj})= text{Bernoulli}( rho_{kj}|  hat{y}_{kj})   
    hat{y}_{kj} rightarrow  frac{1}{1+e^{- ln  tilde{y}_{kj} }}   
        ln  tilde{y}_{kj} rightarrow  boldsymbol{ Psi} ( hat{h}_{1k} ) -  boldsymbol{ Psi} ( hat{h}_{2k} ) -frac{1}{2}(ln 2pi hat{s}^2_{kj} + 1) -frac{1}{2}( ln hat{v} -  boldsymbol{ Psi}(hat{u}) + ln hat{b}_j -  boldsymbol{ Psi} (hat{a}_j) - ln hat{s}^2_{kj} )  +  frac{1}{2} frac{ hat{m}^2_{kj}}{ hat{s}^2_{kj}}
```
"""
function update_y!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} #  formerly update_yjk_mpu!
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].y .= 0.0
        for j in 1:J
            clusters[k].y[j] = E_ln_eta(j,clusters[k]) - E_ln_minus_eta(j,clusters[k]) - 0.5 * (log(2*pi*clusters[k].s_sq_mu[j]) + 1)  - 0.5 * (E_ln_sigma_sq(j,clusters[k]) + E_ln_lambda(j,clusters[k]) - log(clusters[k].s_sq_mu[j])) +  0.5 * (clusters[k].m_mu[j]) ^2 /(clusters[k].s_sq_mu[j])#*E_one_over_sigma_sq(j,clusters[k]) * E_one_over_lambda(j,clusters[k]) #*clusters[k].Nk[1] 
        end
        # sigmoidNorm!(clusters[k].y) 
        clusters[k].y .= sigmoid.(clusters[k].y) 
    end
    return clusters
end


"""
    update_h1!(clusters,dataparams,modelparams)
This function updates the h1 parameter in the ClusterFeature type. This is a varational parameter used to to appoximate the eta parameter in the model. The eta parameter follows a Beta distribution in the true posterior. The Beta parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p ( eta_{k}| varphi_1, varphi_2 )= text{Beta} ( eta_{k}| varphi_1, varphi_2 )   
    q( eta_{k}| hat{h}_1, hat{h}_2)= text{Beta}( eta_{k}| hat{h}_{1k}, hat{h}_{2k})   
        hat{h}_{1k} rightarrow  varphi_1 +   sum_{j=1}^{J}  hat{y}_{kj} 
```
"""
function update_h1!( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; eta_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K

    if eta_update_mode == "Clusterwise"
        for k in 1:K
            clusters[k].cache .= 0.0
            clusters[k].cache[1] +=  modelparams.varphi1[1]
            for j in 1:J
                clusters[k].cache[1] += clusters[k].y[j]
            end
            clusters[k].h1 .= clusters[k].cache[1] .* ones(J)
        end
    elseif  eta_update_mode == "Genewise"
        clusters[1].cache .= 0.0
        clusters[1].cache .+=  modelparams.varphi1[1]*ones(J)
        for k in 1:K
            clusters[1].cache .+= clusters[k].y
        end
        for k in 1:K
            clusters[k].h1 .= clusters[1].cache
        end
    elseif  eta_update_mode == "Global"
        clusters[1].cache .= 0.0
        clusters[1].cache[1] +=  modelparams.varphi1[1]
        for k in 1:K
            for j in 1:J
                clusters[1].cache[1] += clusters[k].y[j]
            end
        end
        for k in 1:K
            clusters[k].h1 .= clusters[1].cache[1] .* ones(J)
        end
    elseif  eta_update_mode == "Local"
        for k in 1:K
            clusters[k].h1 .= 0.0
            # clusters[k].h1 .+= modelparams.varphi1[1]*ones(J)
            for j in 1:J
                clusters[k].h1[j] +=  modelparams.varphi1[1] + clusters[k].y[j]
            end
        end
    else
        error("Invalid eta_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return clusters
end



"""
    update_h2!(clusters,dataparams,modelparams)
This function updates the h2 parameter in the ClusterFeature type. This is a varational parameter used to to appoximate the eta parameter in the model. The eta parameter follows a Beta distribution in the true posterior. The Beta parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p ( eta_{k}| varphi_1, varphi_2 )= text{Beta} ( eta_{k}| varphi_1, varphi_2 )   
    q( eta_{k}| hat{h}_1, hat{h}_2)= text{Beta}( eta_{k}| hat{h}_{1k}, hat{h}_{2k})   
        hat{h}_{2k} rightarrow  varphi_2 +   sum_{j=1}^{J} (1- hat{y}_{kj}) 
```
"""
function update_h2!( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; eta_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    if eta_update_mode == "Clusterwise"
        for k in 1:K
            clusters[k].cache .= 0.0
            clusters[k].cache[1] +=  modelparams.varphi2[1]
            for j in 1:J
                clusters[k].cache[1] += (1 - clusters[k].y[j])
            end
            clusters[k].h2 .= clusters[k].cache[1] .* ones(J)
        end
    elseif  eta_update_mode == "Genewise"
        clusters[1].cache .= 0.0
        clusters[1].cache .+=  modelparams.varphi2[1]*ones(J)
        for k in 1:K
            clusters[1].cache .+= (1 .- clusters[k].y)
        end
        for k in 1:K
            clusters[k].h2 .= clusters[1].cache
        end
    elseif  eta_update_mode == "Global" 
        clusters[1].cache .= 0.0
        clusters[1].cache[1] +=  modelparams.varphi2[1]
        for k in 1:K
            for j in 1:J
                clusters[1].cache[1] += (1 - clusters[k].y[j])
            end
        end
        for k in 1:K
            clusters[k].h2 .= clusters[1].cache[1] .* ones(J)
        end
    elseif  eta_update_mode == "Local"
        for k in 1:K
            clusters[k].h2 .= 0.0
            # clusters[k].h2 .+= *ones(J)
            for j in 1:J
                clusters[k].h2[j] +=  modelparams.varphi2[1]+(1 - clusters[k].y[j])
            end
        end
    else
        error("Invalid eta_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return clusters
end


"""
    update_u!(clusters,dataparams,modelparams)
This function updates the u parameter in the ScalarFeature type. This is a varational parameter used to to appoximate the lambda parameter in the model. The lambda parameter follows a Inverse Gamma distribution in the true posterior. The lambda parameter follows a Inverse Gamma distribution in the variational approximation.  The update is done in place.
```math
    p( lambda| kappa_1, kappa_2)= text{InverseGamma}( lambda| kappa_1, kappa_2)   
    q( lambda| hat{u}, hat{v})= text{InverseGamma}( lambda| hat{u}, hat{v})   
    hat{u} rightarrow  kappa_1 +  frac{1}{2} sum_{k=1}^{K}  sum_{j=1}^{J}   hat{y}_{kj}
```
"""
function update_u!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; lambda_update_mode="Local") where {U <: AbstractFloat, W <: Int64}  #  formerly update_ujk_mpu!
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    if lambda_update_mode == "Clusterwise"
        for k in 1:K
            clusters[k].cache .= 0.0
            clusters[k].cache[1] +=  modelparams.kappa1[1]
            for j in 1:J
                clusters[k].cache[1] += 0.5 * clusters[k].y[j]
            end
            clusters[k].u .= clusters[k].cache[1] .* ones(J)
        end
    elseif lambda_update_mode == "Genewise"
        clusters[1].cache .= 0.0
        clusters[1].cache .+=  modelparams.kappa1[1]*ones(J)
        for k in 1:K
            clusters[1].cache .+= 0.5 .*clusters[k].y
        end
        for k in 1:K
            clusters[k].u .= clusters[1].cache
        end
    elseif lambda_update_mode == "Global"
        clusters[1].cache .= 0.0
        clusters[1].cache[1] +=  modelparams.kappa1[1]
        for k in 1:K
            for j in 1:J
                clusters[1].cache[1] += 0.5 *clusters[k].y[j]
            end
        end
        for k in 1:K
            clusters[k].u .= clusters[1].cache[1] .* ones(J)
        end
    elseif lambda_update_mode == "Local"
        for k in 1:K
            clusters[k].u .= 0.0
            for j in 1:J
                clusters[k].u[j] += modelparams.kappa1[1] + 0.5*clusters[k].y[j]
            end
        end
    else
        error("Invalid lambda_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return clusters
end


"""
    update_v!(scalars,dataparams,modelparam)
This function updates the v parameter in the ScalarFeature type. This is a varational parameter used to to appoximate the lambda parameter in the model. The lambda parameter follows a Inverse Gamma distribution in the true posterior. The lambda parameter follows a Inverse Gamma distribution in the variational approximation.  The update is done in place.
```math
    p( lambda| kappa_1, kappa_2)= text{InverseGamma}( lambda| kappa_1, kappa_2)   
    q( lambda| hat{u}, hat{v})= text{InverseGamma}( lambda| hat{u}, hat{v})   
    hat{v} rightarrow kappa_2 +  frac{1}{2} sum_{k=1}^{K}  sum_{j=1}^{J}   hat{y}_{kj} * frac{ hat{a}_j}{ hat{b}_j}  * ( hat{s}_{kj}^2+  hat{m}_{kj}^2 )
```
"""
function update_v!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; lambda_update_mode="Local") where {U <: AbstractFloat, W <: Int64} #  formerly update_vjk_mpu!
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    if lambda_update_mode == "Clusterwise"
        for k in 1:K
            clusters[k].cache .= 0.0
            clusters[k].cache[1] +=  modelparams.kappa2[1]
            for j in 1:J
                clusters[k].cache[1] += 0.5 * clusters[k].y[j] * E_one_over_sigma_sq(j,clusters[k]) * E_mu_sq(j,clusters[k])
            end
            clusters[k].v .= clusters[k].cache[1] .* ones(J)
        end
    elseif lambda_update_mode == "Genewise"
        clusters[1].cache .= 0.0
        clusters[1].cache .+=  modelparams.kappa2[1]*ones(J)
        for k in 1:K
            for j in 1:J
                clusters[1].cache[j] .+= 0.5 .*clusters[k].y[j] .* E_one_over_sigma_sq.(j,clusters[k]) .* E_mu_sq.(j,clusters[k])
            end
        end
        for k in 1:K
            clusters[k].v .= clusters[1].cache
        end
    elseif lambda_update_mode == "Global"
        clusters[1].cache .= 0.0
        clusters[1].cache[1] +=  modelparams.kappa2[1]
        for k in 1:K
            for j in 1:J
                clusters[1].cache[1] += 0.5 *clusters[k].y[j] * E_one_over_sigma_sq(j,clusters[k]) * E_mu_sq(j,clusters[k])
            end
        end
        for k in 1:K
            clusters[k].v .= clusters[1].cache[1] .* ones(J)
        end
    elseif lambda_update_mode == "Local"
        for k in 1:K
            clusters[k].v .= 0.0
            for j in 1:J
                clusters[k].v[j] +=  modelparams.kappa2[1] + 0.5*clusters[k].y[j] * E_one_over_sigma_sq(j,clusters[k]) * E_mu_sq(j,clusters[k])
            end
        end
    else
        error("Invalid lambda_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return clusters
end

"""
    update_a!(clusters,dataparams,modelparams)
This function updates the a parameter in the GeneFeature type. This is a varational parameter used to to appoximate the sigma^2 parameter in the model. The sigma^2 parameter follows a Inverse Gamma distribution in the true posterior. The sigma^2 parameter follows a Inverse Gamma distribution in the variational approximation.  The update is done in place.
```math
    p( sigma^2_j| xi_1, xi_2)= text{InverseGamma}( sigma^2_j| xi_1, xi_2)   
    q( sigma^2_{j}| hat{a}_j, hat{b}_j)= text{InverseGamma}( sigma^2_{j}| hat{a}_j, hat{b}_j)   
    hat{a}_{j} rightarrow  xi_1   +  frac{1}{2}(N+sum_{k=1}^{K}hat{y}_{kj})
```
"""
function update_a!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; sigma_update_mode="Local") where {U <: AbstractFloat, W <: Int64} 
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    if sigma_update_mode == "Clusterwise"
        error("sigma_update_mode not implemented for 'Clusterwise'")
    elseif sigma_update_mode == "Genewise"
        clusters[1].cache .= 0.0
        # clusters[1].cache .+=  *ones(J)
        for k in 1:K
            clusters[1].cache .+= clusters[k].y
        end
        for k in 1:K
            clusters[k].a .= modelparams.xi1[1] .+ 0.5 .* (clusters[1].cache .+ N)
        end
    elseif sigma_update_mode == "Global"
        error("sigma_update_mode not implemented for 'Global'")
    elseif sigma_update_mode == "Local"
        for k in 1:K
            clusters[k].a .= 0.0
            for j in 1:J
                clusters[k].a[j] += modelparams.xi1[1] + 0.5 * (clusters[k].y[j] + clusters[k].Nk[1])
            end
        end
    else
        error("Invalid sigma_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return clusters
end

"""
    update_b!(clusters,clusters,dataparams,modelparams)
This function updates the b parameter in the GeneFeature type. This is a varational parameter used to to appoximate the sigma^2 parameter in the model. The sigma^2 parameter follows a Inverse Gamma distribution in the true posterior. The sigma^2 parameter follows a Inverse Gamma distribution in the variational approximation.  The update is done in place.
```math
    p( sigma^2_j| xi_1, xi_2)= text{InverseGamma}( sigma^2_j| xi_1, xi_2)   
    q( sigma^2_{j}| hat{a}_j, hat{b}_j)= text{InverseGamma}( sigma^2_{j}| hat{a}_j, hat{b}_j)   
    hat{b}_{j} rightarrow xi_2 + frac{1}{2}sum_{k=1}^{K} frac{hat{u}_{k}}{hat{v}_{k}}hat{y}_{kj}(hat{s}^2_{kj} + hat{m}^2_{kj}) +sum_{k=1}^{K}left[ frac{1}{2}left[hat{x}_{kj}^{2}-2hat{x}_{kj}left(hat{y}_{kj}hat{m}_{kj} +hat{m}_{nu j} right)+N_khat{y}_{kj}(hat{m}_{kj}^2 +hat{s}_{kj}^2 )+N_k( hat{m}_{nu j}^2 +hat{s}_{nu j}^2 ) + 2N_khat{y}_{kj}hat{m}_{mu kj}hat{m}_{nu j} right]right]
```
"""
function update_b!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; sigma_update_mode="Local") where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    if sigma_update_mode == "Clusterwise"
        error("sigma_update_mode not implemented for 'Clusterwise'")
    elseif sigma_update_mode == "Genewise"
        clusters[1].cache .= 0.0
        clusters[1].cache .+= modelparams.xi2[1] .*ones(J)
        for k in 1:K
            for j in 1:J
                clusters[1].cache[j] += 0.5 * clusters[k].y[j] * E_one_over_lambda(j,clusters[k]) * E_mu_sq(j,clusters[k])
                clusters[1].cache[j] += 0.5 * (clusters[k].x_hat_sq[j] - 2 * clusters[k].x_hat[j] * (clusters[k].y[j] * clusters[k].m_mu[j]  +  clusters[k].m_nu[j]) + clusters[k].Nk[1] * clusters[k].y[j] * E_mu_sq(j,clusters[k]) + clusters[k].Nk[1] * E_nu_sq(j,clusters[k]) + 2 * clusters[k].Nk[1] * clusters[k].m_nu[j]*clusters[k].y[j]*clusters[k].m_mu[j])
            end
        end
        for k in 1:K
            clusters[k].b .= clusters[1].cache
        end
    elseif sigma_update_mode == "Global"
        error("sigma_update_mode not implemented for 'Global'")
    elseif sigma_update_mode == "Local"
        for k in 1:K
            clusters[k].b .= 0.0
            for j in 1:J
                clusters[k].b[j] += modelparams.xi2[1]
                clusters[k].b[j] += 0.5 * clusters[k].y[j] * E_one_over_lambda(j,clusters[k]) * E_mu_sq(j,clusters[k])
                clusters[k].b[j] += 0.5 * (clusters[k].x_hat_sq[j] - 2 * clusters[k].x_hat[j] * (clusters[k].y[j] * clusters[k].m_mu[j]  +  clusters[k].m_nu[j]) + clusters[k].Nk[1] * clusters[k].y[j] * E_mu_sq(j,clusters[k]) + clusters[k].Nk[1] * E_nu_sq(j,clusters[k]) + 2 * clusters[k].Nk[1] * clusters[k].m_nu[j]*clusters[k].y[j]*clusters[k].m_mu[j])
            end
        end
    else
        error("Invalid sigma_update_mode. Please choose from 'Clusterwise', 'Genewise', 'Global', 'Local'")
    end
    return clusters
end


"""
    update_m_mu!(clusters,dataparams,modelparams)
This function updates the m parameter in the ClusterFeature type. This is a varational parameter used to to appoximate the mu parameter in the model. The mu parameter follows a Normal distribution in the true posterior. The mu parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p( mu_{kj}| lambda, sigma^2_j) = text{Normal}( mu_{kj}| lambda, sigma^2_j)   
    q( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})= hat{y}_{kj} text{Normal}( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})   
    hat{m}_{kj} rightarrow  frac{hat{x}_{kj} - hat{m}_{nu j}N_{k}}{left(frac{hat{u}_{k}}{hat{v}_{k}} + N_{k} right)}
```
"""
function update_m_mu!( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} 
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        for j in 1:J
            clusters[k].m_mu[j] = (clusters[k].x_hat[j] - clusters[k].m_nu[j]*clusters[k].Nk[1])/ (E_one_over_lambda(j,clusters[k]) + clusters[k].Nk[1])
        end
    end
    return clusters
end

"""
    update_s_sq_mu!(clusters,dataparams,modelparams)
This function updates the s^2 parameter in the ClusterFeature type. This is a varational parameter used to to appoximate the mu parameter in the model. The mu parameter follows a Normal distribution in the true posterior. The mu parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p( mu_{kj}| lambda, sigma^2_j) = text{Normal}( mu_{kj}| lambda, sigma^2_j)   
    q( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})= text{Normal}( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})   
    hat{s}^{2}_{kj} rightarrow left[frac{hat{a}_{j}}{hat{b}_{j}}left(frac{hat{u}_{k}}{hat{v}_{k}} + N_{k}right) right]^{-1}
```
"""
function update_s_sq_mu!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        for j in 1:J
            clusters[k].s_sq_mu[j] = 1 / (E_one_over_sigma_sq(j,clusters[k]) * (E_one_over_lambda(j,clusters[k]) + clusters[k].Nk[1]))
        end
    end
    return clusters
end

 

"""
    update_m_nu!(clusters,dataparams,modelparams)
This function updates the m parameter in the ClusterFeature type. This is a varational parameter used to to appoximate the mu parameter in the model. The mu parameter follows a Normal distribution in the true posterior. The mu parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p( mu_{kj}| lambda, sigma^2_j) = text{Normal}( mu_{kj}| lambda, sigma^2_j)   
    q( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})= hat{y}_{kj} text{Normal}( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})   
    hat{m}_{kj} rightarrow  frac{hat{x}_{kj} - hat{m}_{nu j}N_{k}}{left(frac{hat{u}_{k}}{hat{v}_{k}} + N_{k} right)}
```
"""
function update_m_nu!( clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64} 
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    clusters[1].cache .= 0.0
    for j in 1:J
        numerator = 0.0
        denominator = 0.0
        for k in 1:K
            numerator += E_one_over_sigma_sq(j,clusters[k])*(clusters[k].x_hat[j] - clusters[k].m_mu[j]*clusters[k].y[j]*clusters[k].Nk[1])
            denominator += E_one_over_sigma_sq(j,clusters[k])*clusters[k].Nk[1]
        end
        clusters[1].cache[j] = (numerator + (modelparams.nu0[j]/modelparams.sigma_sq_nu[j])) / (denominator + 1/modelparams.sigma_sq_nu[j])
    end
    for k in 1:K
        clusters[k].m_nu .= clusters[1].cache
    end
    return clusters
end

"""
    update_s_sq_nu!(clusters,dataparams,modelparams)
This function updates the s^2 parameter in the ClusterFeature type. This is a varational parameter used to to appoximate the mu parameter in the model. The mu parameter follows a Normal distribution in the true posterior. The mu parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p( mu_{kj}| lambda, sigma^2_j) = text{Normal}( mu_{kj}| lambda, sigma^2_j)   
    q( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})= text{Normal}( mu_{kj}|  hat{m}_{kj}, hat{s}^2_{kj})   
    hat{s}^{2}_{kj} rightarrow left[frac{hat{a}_{j}}{hat{b}_{j}}left(frac{hat{u}_{k}}{hat{v}_{k}} + N_{k}right) right]^{-1}
```
"""
function update_s_sq_nu!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    clusters[1].cache .= 0.0
    for j in 1:J
        denominator = 0.0
        for k in 1:K
            denominator += E_one_over_sigma_sq(j,clusters[k])*clusters[k].Nk[1]
        end
        clusters[1].cache[j] = 1 /(denominator  + 1/modelparams.sigma_sq_nu[j])
    end
    for k in 1:K
        clusters[k].s_sq_nu .= clusters[1].cache
    end
    return clusters
end

"""
    g_transfromed_surrogate_lb(x::Vector{U}, T::Vector{Int64}, K::Int, alpha0::Vector{U}, gamma0::U,alpha_Tks::Vector{U};use_log=true)  where {U <: AbstractFloat}
This function computes the surrogate lower bound for the g1 and g2 parameters in the variational approximation of the eta parameter in the model. The eta parameter follows a Beta distribution in the true posterior. The Beta parameter follows a Normal distribution in the variational approximation. The surrogate lower bound is used to optimize the g1 and g2 parameters using gradient-based optimization methods. The function takes as input the transformed variables x, the T vector, the number of clusters K, the alpha0 vector, the gamma0 scalar, and the alpha_Tks vector. The function returns the negative of the surrogate lower bound since Optim.jl minimizes functions.
"""
function g_transfromed_surrogate_lb(x::Vector{U}, T::Vector{Int64}, K::Int, alpha0::Vector{U}, gamma0::U,alpha_Tks::Vector{U};use_log=true)  where {U <: AbstractFloat}
    # Transform variables
    I = length(T)
    T_all = sum(T)
    g1 = logistic.(x[1:K]) # Logistic transform to keep 0 < g1 < 1
    g2 = exp.(x[K+1:2K])  # Exponential transform to keep g2 > 0
    # Now define the objective function with g1 and g2
    L_G = 0.0
    for it in 1:T_all
        L_G += K*log(alpha0[it])
    end
    for k in 1:K
        alpha_part = g1[k]*g2[k]
        beta_part = (1 - g1[k])*g2[k]
        loggamma_term1 = loggamma(alpha_part)
        loggamma_term2 = loggamma(beta_part)
        loggamma_term3 = loggamma(alpha_part + beta_part)
        e_ln_chi_term = e_ln_chi(g1[k],g2[k])
        e_ln_minus_chi_term = e_ln_minus_chi(g1[k],g2[k])
        e_sbk = expectation_sbk(k,K,g1,g2; use_log=use_log)
        L_G += loggamma_term1 + loggamma_term2 - loggamma_term3 + (T_all+1-alpha_part)*e_ln_chi_term + (T_all*(K+1-k)+gamma0-beta_part)*e_ln_minus_chi_term +e_sbk*alpha_Tks[k]
    end
    return -L_G  # Since Optim.jl minimizes, return negative of the objective
end

"""
    update_g1g2!(clusters,dataparams,modelparams,use_log=true)
This function updates the g1 and g2 parameters in the ClusterFeature type. This is a varational parameter used to to appoximate the eta parameter in the model. The eta parameter follows a Beta distribution in the true posterior. The Beta parameter follows a Normal distribution in the variational approximation.  The update is done in place.
```math
    p( chi_k|1, gamma_0)= text{Beta}( chi_k|1, gamma_0)   
    q( chi_k| hat{g}_{1k}, hat{g}_{2k})= text{Beta}( chi_k| hat{g}_{1k} hat{g}_{2k},(1- hat{g}_{1k}) hat{g}_{2k})   
        hat{g}_{1k}, hat{g}_{2k} &=  text{argmax}_{ hat{g}_{1k}, hat{g}_{2k}}  mathcal{L}_{G}( cdot)  
        text{s.t. } 0 <  hat{g}_{1k} < 1&,  hat{g}_{2k} > 0  text{ for } k = 1,...,K
```
"""
function update_g1g2!(clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature; use_log=true) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    transformed_g1 = zeros(K)
    transformed_g2 = zeros(K)
    alpha_Tks = zeros(K)
    transformed_x = zeros(2*K)
    alpha0 = recursive_flatten(modelparams.alpha0)
    gamma0 = modelparams.gamma0[1]
    for k in 1:K
        transformed_g1[k] = logit(clusters[k].g1[1]) 
        transformed_g2[k] = log(clusters[k].g2[1])
        alpha_Tks[k] = clusters[k].alpha_Tk[1]
    end
    transformed_x[1:K] = transformed_g1
    transformed_x[K+1:2K] = transformed_g2
    # result = optimize(x -> g_transfromed_surrogate_lb(x, T, K, alpha0, gamma0, alpha_Tks; use_log=use_log), transformed_x, LBFGS(;m=2))#Optim.Adam(;  alpha=0.1)
    result = optimize(x -> g_transfromed_surrogate_lb(x, T, K, alpha0, gamma0, alpha_Tks; use_log=use_log), transformed_x, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = 100))#Optim.Adam(;  alpha=0.1)
    # result = optimize(x -> g_transfromed_surrogate_lb(x, T, K, alpha0, gamma0, alpha_Tks; use_log=use_log), transformed_x, LBFGS(),Optim.Options(iterations = 100))#Optim.Adam(;  alpha=0.1)
    # Extract the optimized values and transform back
    optimal_x = Optim.minimizer(result)
    optimal_g1 = logistic.(optimal_x[1:K])#1 ./ (1 .+ exp.(-optimal_x[1:K]))  # Transform back to 0 < g1 < 1
    optimal_g2 = exp.(optimal_x[K+1:2K])#            # Transform back to g2 > 0
    for k in 1:K
        clusters[k].g1[1] = optimal_g1[k]
        if optimal_g1[k] == 0.0
            optimal_g1[k] = 1e-179
        else
            clusters[k].g1[1] = optimal_g1[k]
        end
        if optimal_g1[k] == 1.0
            optimal_g1[k] = 1-1e-12
        else
            clusters[k].g1[1] = optimal_g1[k]
        end
        if optimal_g2[k] == 0.0
            optimal_g2[k] = 1e-6
        else
            clusters[k].g2[1] = optimal_g2[k]
        end
        # clusters[k].g2[1] = optimal_g2[k]
    end
    return clusters
end


"""
    update_x_hat!(cells,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_x_hat!(cells,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    J = dataparams.J
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        @. clusters[k].x_hat = 0.0
        @inbounds @fastmath for n in 1:N
            @. clusters[k].x_hat +=   cells[n].x * cells[n].r[k]
        end

    end
    return clusters
    # return x_hat_k
end
# ::Vector{CellFeature{U,W,J}}

"""
    update_Ctt!(cells,conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_Ctt!(cells,conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    J = dataparams.J
    N = dataparams.N
    K = modelparams.K
    T = dataparams.T
    T_all = sum(T)
    # Threads.@threads for t in 1:T
    #     @. conditions[t].Ctt = 0.0
    #     @inbounds @fastmath for n in 1:N
    #         @. conditions[t].Ctt +=  cells[n].c[t]
    #     end
    # end
    Threads.@threads for it in 1:T_all
        @. conditions[it].Ctt = 0.0
        # condition_t = conditions[it].t
        @inbounds @fastmath for n in dataparams.TimeRanges[it][1]:dataparams.TimeRanges[it][2]
            @. conditions[it].Ctt +=  cells[n].c
        end
    end
    return conditions
    # return x_hat_k
end
# ::Vector{CellFeature{U,W,J}}

"""
    update_Nk!(cells,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_Nk!(cells,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].Nk[1] = 0.0
        Nk_sum = 0.0
        @inbounds @fastmath for n in 1:N #
            Nk_sum += cells[n].r[k]
        end
        clusters[k].Nk[1] = Nk_sum
    end
    return clusters
end
# ::Vector{CellFeature{U,W,J}}

"""
    update_x_hat_sq!(cells,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_x_hat_sq!(cells,clusters::Vector{ClusterFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    J = dataparams.J
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        @. clusters[k].x_hat_sq = 0.0
        @inbounds @fastmath for n in 1:N
            @. clusters[k].x_hat_sq +=   cells[n].xsq * cells[n].r[k]
        end
    end
    return clusters
end

# ::Vector{CellFeature{U,W,J}}
"""
    update_CNtk!(cells,matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_CNtk!(cells,matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    N = dataparams.N
    T = dataparams.T
    K = modelparams.K
    T_all = sum(T)
    # Threads.@threads for t in 1:T
    #     @. matrixconditions[t].CNtk = 0.0
    #     @inbounds @fastmath for n in 1:N
    #         @. matrixconditions[t].CNtk +=  cells[n].c * cells[n].r'
    #     end
    # end
    Threads.@threads for it in 1:T_all
        @. matrixconditions[it].CNtk = 0.0
        @inbounds @fastmath for n in dataparams.TimeRanges[it][1]:dataparams.TimeRanges[it][2]
            @. matrixconditions[it].CNtk +=  cells[n].c * cells[n].r'
        end
    end
    return matrixconditions
end

# ::Vector{CellFeature{U,W,J}}

"""
    update_alpha_Tk_stats!(clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
"""
function update_alpha_Tk_stats!(clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    float_type = dataparams.BitType
    J = dataparams.J
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1
    T_all = sum(T)
    for k in 1:K
        clusters[k].alpha_Tk[1] = 0.0
        for it in 1:T_all
            condition_i = conditions[it].i
            condition_t = conditions[it].t
            clusters[k].alpha_Tk[1] += modelparams.alpha0[condition_i][condition_t] * E_ln_pi(k,conditions[it])
        end
    end
    return clusters
end

############################################
############################################
############################################


"""
    update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams)
"""
function update_λ_sq_hat_mpu!(geneparams,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    K = modelparams.K

    for j in 1:G
        geneparams[j].cache[1] = 0.0
        yjk_sum = 0.0 
        for k in 1:K
            yjk_sum += clusters[k].yjk_hat[j]
        end
        for k in 1:K
            geneparams[j].cache[1] += 10. + clusters[k].var_muk[j] + clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end
        geneparams[j].λ_sq[1] = geneparams[j].cache[1] ./ (10. + yjk_sum)
    end
    return geneparams
end

"""
    update_yjk_mpu!(clusters,geneparams,dataparams,modelparams)
"""
function update_yjk_mpu!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    for k in 1:K

        clusters[k].yjk_hat .= 0.0
        logodds = log10((modelparams.ηk[1])  / (1 - (modelparams.ηk[1])))
        logp10 = log(10)*logodds 
        for j in 1:G
            clusters[k].yjk_hat[j] =  logp10 + log(sqrt(clusters[k].v_sq_k_hat[j]) / (sqrt(geneparams[j].λ_sq[1]))) + 0.5 * (clusters[k].mk_hat[j]) ^2 /clusters[k].v_sq_k_hat[j] 
        end
        sigmoidNorm!(clusters[k].yjk_hat) 
    end


end

"""
    update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams)
"""
function update_mk_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].mk_hat .= 0.0
        for j in 1:G
            clusters[k].mk_hat[j] +=  (geneparams[j].λ_sq[1] * clusters[k].x_hat[j]) /  (clusters[k].Nk[1] * geneparams[j].λ_sq[1]  + clusters[k].σ_sq_k_hat[j]  )
        end
    end

    return clusters
end

"""
    update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
"""
function update_v_sq_k_hat_mpu!(clusters,geneparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].v_sq_k_hat .= 0.0
        for j in 1:G
            clusters[k].v_sq_k_hat[j] +=   (geneparams[j].λ_sq[1] * clusters[k].σ_sq_k_hat[j] ) / (clusters[k].Nk[1] * geneparams[j].λ_sq[1] + clusters[k].σ_sq_k_hat[j])
        end  
    end
    return clusters
end

"""
    update_var_muk_hat_mpu!(clusters, dataparams,modelparams)
"""
function update_var_muk_hat_mpu!(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    for k in 1:K
        clusters[k].var_muk .= 0.0
        for j in 1:G
            clusters[k].var_muk[j] +=  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j] ^2  + clusters[k].v_sq_k_hat[j])   -  clusters[k].yjk_hat[j] * (clusters[k].mk_hat[j]) ^2
        end

    end
    return clusters
end
"""
    update_κk_hat_mpu!(clusters, dataparams,modelparams)
"""
function update_κk_hat_mpu!(clusters, dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    for k in 1:K
        clusters[k].κk_hat .= 0.0

        clusters[k].κk_hat .+=  clusters[k].yjk_hat .* clusters[k].mk_hat
    end  

    return clusters
end

"""
    update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams)
"""
function update_σ_sq_k_hat_mpu!(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K

    for k in 1:K
        clusters[k].σ_sq_k_hat .= 0.0

        clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 10.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* clusters[k].κk_hat .+  clusters[k].Nk  .* (clusters[k].var_muk .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ 10.0)
        # clusters[k].σ_sq_k_hat .+=   1 ./(clusters[k].Nk .+ 10.0) .* (clusters[k].x_hat_sq .- 2.0 .*  clusters[k].x_hat .* (clusters[k].yjk_hat .* clusters[k].mk_hat) .+  clusters[k].Nk  .* ((clusters[k].yjk_hat .* (clusters[k].mk_hat .^2  .+ clusters[k].v_sq_k_hat)   .-  clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ clusters[k].yjk_hat .* (clusters[k].mk_hat) .^2) .+ 10.0)
        
    end

    return clusters
end

"""
    update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams)
"""
function update_rtik_mpu!(cellpop,clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for i in 1:N
        t = first(dataparams.LinearAddress[i])
        adjust_e_log_π_tk3!(t,conditionparams)
        # cellpop[i]._reset!(cellpop[i].rtik,cellpop[i].BitType)
        cellpop[i].rtik .= 0.0
        for k in 1:K
            expectation_log_normal_l_j!(clusters[k],cellpop[i], dataparams)
            cell_gene_sums = 0.0
            for el in cellpop[i].cache
                cell_gene_sums+=el
            end
            # cell_gene_sums  = sum(cellpop[i].cache)
            cellpop[i].rtik[k] = conditionparams[t].e_log_π_t_cache[k] + cell_gene_sums
        end
        norm_weights3!(cellpop[i].rtik)
    end
    return cellpop
end

"""
    update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams)
"""
function update_Ntk_mpu!(cellpop,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    @inbounds for t in 1:T
        _reset!(conditionparams[t].Ntk,float_type)
        st,en = dataparams.TimeRanges[t]
        @inbounds for i in st:en
            @inbounds for k in 1:K #@fastmath 
                conditionparams[t].Ntk[k] += cellpop[i].rtik[k]
            end
        end
    end

    return conditionparams
end

"""
    update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams)
"""
function update_d_hat_mpu!(clusters,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K + 1

    @inbounds for t in 1:T
        _reset!(conditionparams[t].d_hat_t, float_type)
        @inbounds for k in 1:Kplus
            e_βk = expectation_βk(k, clusters, modelparams)
            updated_Ntk_sum = 0.0
            @inbounds for tt in 1:T
                updated_Ntk_sum += conditionparams[t].c_tt_prime[tt] * conditionparams[tt].Ntk[k]
            end
            conditionparams[t].d_hat_t[k] = modelparams.α0 * e_βk + updated_Ntk_sum
        end
    end
    return conditionparams
end

"""
    update_d_hat_sum_mpu!(conditionparams,dataparams)
"""
function update_d_hat_sum_mpu!(conditionparams,dataparams)
    float_type = dataparams.BitType
    T = dataparams.T
    for t in 1:T
        conditionparams[t].d_hat_t_sum[1] = sum(conditionparams[t].d_hat_t)
    end
    return conditionparams
end

"""
    update_c_ttprime_mpu!(conditionparams,dataparams,modelparams)
"""
function update_c_ttprime_mpu!(conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    log_π_expected_value_fast3!(conditionparams,dataparams,modelparams)
    T = dataparams.T
    Kplus = modelparams.K + 1
    
    for t in 1:T
        _reset!(conditionparams[t].c_tt_prime, float_type)
        for tt in 1:t
            ctt_k_accum = 0.0
            for k in 1:Kplus
                ctt_k_accum += conditionparams[t].Ntk[k] * conditionparams[tt].e_log_π_t_cache[k]
            end
            ctt_w_accum = 0.0
            if !isone(t)
                if isone(tt)
                    e_log_tilde_wt = 0.0
                else
                    e_log_tilde_wt = log_tilde_wt_expected_value(1.0, conditionparams[tt-1].st_hat[1])#e_log_tilde_wt_vec[tt-1]
                end
                sum_e_log_minus_tilde_wt = 0.0
                if tt != t
                    for tprime in tt:t-1
                        sum_e_log_minus_tilde_wt += expectation_log_minus_tilde_wtt(1.0, conditionparams[tprime].st_hat[1])
                    end
                end
                ctt_w_accum += sum_e_log_minus_tilde_wt + e_log_tilde_wt
            end
            conditionparams[t].c_tt_prime[tt] = ctt_k_accum + ctt_w_accum
        end
        norm_weights3!(t, conditionparams[t].c_tt_prime)
    end
    return conditionparams
end

"""
    update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams)
"""
function update_Tk_mpu!(Tk,conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    Kplus = K+1
    Tk .= 0.0
    @inbounds for k in 1:Kplus
        @inbounds @fastmath for t in 1:T
            e_log_π = expectation_log_π_tk(conditionparams[t].d_hat_t[k],conditionparams[t].d_hat_t_sum[1])
            Tk[k] +=  e_log_π
        end
    end
end

"""
    update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=100000)
"""
function update_gh_hat_mpu!(clusters,dataparams,modelparams,Tk;optim_max_iter=100000)
    float_type = dataparams.BitType
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0
    a_hat = [clusters[k].ak_hat[1] for k in 1:K ]
    b_hat =  [clusters[k].bk_hat[1] for k in 1:K ]
    ab_vec_args = [a_hat , b_hat]
    ab_vec_args0  = permutedims(reduce(hcat,ab_vec_args))
    LB_LG_unconstrained = SurragateLowerBound_unconstrained_closure(T,γ0,α0,Tk)
    gg_uncon! = g_unconstrained_closure!(T,γ0,α0,Tk)
    lb_lg_results = Optim.maximize(LB_LG_unconstrained,gg_uncon!,ab_vec_args0, GradientDescent(linesearch=Optim.LineSearches.BackTracking()),Optim.Options(iterations = optim_max_iter))
    for k in 1:K
        clusters[k].gk_hat[1] = sigmoid(lb_lg_results.res.minimizer[1,k])#new_rho_hat[k]
        clusters[k].hk_hat[1] = exp(lb_lg_results.res.minimizer[2,k])#new_omega_hat[k]
        clusters[k].ak_hat[1] = StatsFuns.logit(clusters[k].gk_hat[1]) #new_c_hat[k]
        clusters[k].bk_hat[1] = log(clusters[k].hk_hat[1])#new_d_hat[k]
    end
    return clusters
end

"""
    update_Nk_mpu!(cellpop,clusters,dataparams,modelparams)
"""
function update_Nk_mpu!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    N = dataparams.N
    K = modelparams.K
    for k in 1:K
        clusters[k].Nk .= 0.0
        Nk_sum = 0.0
        for i in 1:N
            Nk_sum += cellpop[i].rtik[k]
        end
        clusters[k].Nk[1] = Nk_sum
    end
    return clusters
end

"""
    update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams)
"""
function update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].x_hat .= 0.0
        @inbounds @fastmath for i in 1:N
            clusters[k].x_hat .+=   cellpop[i].x .* cellpop[i].rtik[k]
        end

    end

    return clusters
    # return x_hat_k
end

"""
    update_x_hat_k_mpu!(cellpop,clusters,dataparams,modelparams)
"""
function update_x_hat_sq_k_mpu!(cellpop,clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    N = dataparams.N
    K = modelparams.K
    Threads.@threads for k in 1:K
        clusters[k].x_hat_sq .= 0.0
        @inbounds @fastmath for i in 1:N
            clusters[k].x_hat_sq .+=   cellpop[i].xsq .* cellpop[i].rtik[k]
        end

    end
    return clusters
end

"""
    update_st_hat_mpu!(conditionparams,dataparams,modelparams)
"""
function update_st_hat_mpu!(conditionparams,dataparams,modelparams) 
    float_type = dataparams.BitType
    T = dataparams.T

    for t in 2:T
        st_hat = 0.0
        st_hat += modelparams.ϕ0
        for t_prime in t:T
            for l in 1:t-1
                st_hat += conditionparams[t_prime].c_tt_prime[l]
            end
        end
        conditionparams[t-1].st_hat[1]=st_hat
    end
    conditionparams[T].st_hat[1]=0.0
    return conditionparams
end

"""
    update_ηk!(clusters,dataparams,modelparams)
"""
function update_ηk!(clusters,dataparams,modelparams)
    float_type = dataparams.BitType
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K


    for k in 1:K

        modelparams.ηk[1]  = 0.0
        modelparams.ηk[1]  = sum(clusters[k].yjk_hat) / sum(1 .- clusters[k].yjk_hat)
        sigmoidNorm!(modelparams.ηk) 
    end

end
