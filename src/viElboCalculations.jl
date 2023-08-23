"""
        calc_SurragateLowerBound(rho_hat,omega_hat,T,γ,α0,Tk)
    Calculates the Surrogate Bound
    ```math

    ```
"""
function calc_SurragateLowerBound(rho_hat,omega_hat,T,γ,α0,Tk)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K]
    k_vec = collect(1:K)
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α0 .* e_βk .* Tk
    lb_lg = sum(lb_lg_k)
    return lb_lg
end

"""
        calc_SurragateLowerBound_unconstrained(c,d,T,γ,α0,Tk)
    Calculates the unconstrained Surrogate Bound 
    ```math

    ```
"""
function calc_SurragateLowerBound_unconstrained(c,d,T,γ,α0,Tk)
    rho_hat = sigmoid.(c)
    omega_hat = exp.(d)
    c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    e_logUk = logUk_expected_value(rho_hat,omega_hat)
    e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    K = length(rho_hat)
    e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    k_vec = collect(1:K)
    lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  α0 .* e_βk .* Tk[1:K]
    lb_lg = sum(lb_lg_k)
    return lb_lg
end


"""
        c_Ga(a0, b0)
    Calculates the log of the Gamma function
    ```math

    ```
"""
function c_Ga(a0, b0)
    a0 .* log.(b0) .- loggamma.(a0)
end

"""
       c_Beta(a0, b0)
    Calculates the log of the Beta function
    ```math

    ```
"""
function c_Beta(a0, b0)
    - logbeta.(a0,b0) 
end

# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calculate_elbo(Tk,cellpop,clusters,geneparams,conditionparams,dataparams,modelparams)
#     elbo_val = 0.0
#     elbo_val += calc_DataElbo(cellpop,clusters,geneparams,dataparams,modelparams)
#     elbo_val += calc_Hz_fast3(cellpop,clusters,dataparams)
#     elbo_val += calc_HyjkSurragateLowerBound_unconstrained(Tk,clusters,dataparams,modelparams)
#     elbo_val += calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
#     # elbo_val += calc_alphaElbo_fast3(conditionparams,dataparams,modelparams)
#     # elbo_val += calc_HsGammaAlphaElbo_fast3(a_γ,b_γ,a_γ_hat,b_γ_hat,conditionparams,dataparams,modelparams)
#     elbo_val +=  calc_HsElbo(conditionparams, dataparams, modelparams)
#     return elbo_val
# end


# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calculate_elbo_perK(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
#     K = modelparams.K
#     for k in 1:K
#         elbolog.per_k_elbo[k,iter] = 0.0
#     end
#     elbo_val = 0.0
#     dataelbo,elbolog = calc_DataElbo_perK(cellpop,clusters,geneparams,elbolog,dataparams,modelparams,iter)
#     zentropy = calc_Hz_fast3(cellpop,clusters,dataparams)
#     lg_elbo,elbolog =calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
#     w_elbo = calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
#     sentropy = calc_HsElbo(conditionparams, dataparams, modelparams)
#     elbo_val +=  dataelbo + zentropy + lg_elbo + w_elbo + sentropy
#     for k in 1:K
#         elbolog.per_k_elbo[k,iter] += zentropy + w_elbo + sentropy
#     end
#     return elbo_val,elbolog
# end

"""
        calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
    Calculates the current iterations elbo. Wrapper function that calls the component elbo calculation functions
"""
function calculate_elbo_mpu(Tk,cellpop,clusters,geneparams,conditionparams,elbolog,dataparams,modelparams,iter)
    K = modelparams.K
    for k in 1:K
        elbolog.per_k_elbo[k,iter] = 0.0
    end
    elbo_val = 0.0
    dataelbo,elbolog = calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    # dataelbo = 0.0
    zentropy = calc_Hz_fast3(cellpop,clusters,dataparams)
    lg_elbo,elbolog =calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    # lg_elbo = 0.0
    w_elbo = calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    sentropy = calc_HsElbo(conditionparams, dataparams, modelparams)
    elbo_val +=  dataelbo + zentropy + lg_elbo + w_elbo + sentropy
    for k in 1:K
        elbolog.per_k_elbo[k,iter] += zentropy + w_elbo + sentropy
    end
    return elbo_val,elbolog
end


# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calc_DataElbo(cellpop,clusters,geneparams,dataparams,modelparams)
#     float_type = dataparams.BitType
#     # if isnothing(float_type)
#     #     float_type =eltype(x[1][1])
#     # end
#     G = dataparams.G
#     T = dataparams.T
#     N = dataparams.N
#     K = modelparams.K
#     one_half_const = 1/2
#     λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
#     data_elbo = 0.0
#     @fastmath @inbounds @simd for k in 1:K
#         clusters[k].cache .= 0.0
#         # @fastmath @inbounds @simd for i in 1:N
#         #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
#         # end
#         # try
#         #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#         # catch
#         #     println("K =  $k")
#         #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
#         #     println("λ_sq_vec =  $(λ_sq_vec)")
#         #     println("ηk =  $(modelparams.ηk)")
#         #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
#         #     println("var_muk =  $(clusters[k].var_muk )")
#         #     println("mk_hat =  $(clusters[k].mk_hat)")

#         # end

        
#     #    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#        data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1] ) .+ (1 .- clusters[k].yjk_hat) .* log((1-(modelparams.ηk[1]) )) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#     end
#     #     
#     return data_elbo
# end

# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calc_DataElbo_perK(cellpop,clusters,geneparams,elbolog,dataparams,modelparams,iter)
#     float_type = dataparams.BitType
#     # if isnothing(float_type)
#     #     float_type =eltype(x[1][1])
#     # end
#     G = dataparams.G
#     T = dataparams.T
#     N = dataparams.N
#     K = modelparams.K
#     one_half_const = 1/2
#     λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
#     data_elbo = 0.0
    
#     @fastmath @inbounds @simd for k in 1:K
#         clusters[k].cache .= 0.0
#         yjk_entropy_perK = 0.0
#         # @fastmath @inbounds @simd for i in 1:N
#         #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
#         # end
#         # try
#         #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#         # catch
#         #     println("K =  $k")
#         #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
#         #     println("λ_sq_vec =  $(λ_sq_vec)")
#         #     println("ηk =  $(modelparams.ηk)")
#         #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
#         #     println("var_muk =  $(clusters[k].var_muk )")
#         #     println("mk_hat =  $(clusters[k].mk_hat)")

#         # end

        
#     #    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#     perK_data_elbo = sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1]) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk[1])) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#     yjk_entropy_perK += entropy(clusters[k].yjk_hat)
#     yjk_entropy_perK = -yjk_entropy_perK
#     perK_ebloval =   perK_data_elbo + yjk_entropy_perK
#     elbolog.per_k_elbo[k,iter] += perK_ebloval
#     data_elbo += perK_ebloval
#     end
#     #     
#     return data_elbo,elbolog
# end

"""
        calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    Calculates the current iterations data elbo.
"""
function calc_DataElbo_mpu(clusters,geneparams,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    G = dataparams.G
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    one_half_const = 1/2
    # λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
    data_elbo = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        clusters[k].cache .= 0.0
        yjk_entropy_perK = 0.0
        # @fastmath @inbounds @simd for i in 1:N
        #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
        # end
        # try
        #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
        # catch
        #     println("K =  $k")
        #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
        #     println("λ_sq_vec =  $(λ_sq_vec)")
        #     println("ηk =  $(modelparams.ηk)")
        #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
        #     println("var_muk =  $(clusters[k].var_muk )")
        #     println("mk_hat =  $(clusters[k].mk_hat)")

        # end

        
    ####    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


    # perK_data_elbo = sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1]) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk[1])) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


    # perK_data_elbo = sum( .+ . .+  .+)

    perK_data_elbo = 0.0
    for j in 1:G
        perK_data_elbo += -1*one_half_const * clusters[k].Nk[1] * log(2π)
        perK_data_elbo += -1*one_half_const * clusters[k].Nk[1]  * log(clusters[k].σ_sq_k_hat[j])
        perK_data_elbo += -1*one_half_const * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat_sq[j] 
        perK_data_elbo +=  1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat[j] * clusters[k].κk_hat[j]
        perK_data_elbo +=  -1*one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].var_muk[j] 
        perK_data_elbo += -1* one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 - one_half_const * clusters[k].yjk_hat[j] * log(geneparams[j].λ_sq[1]) 
        perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].var_muk[j]
        perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 
        perK_data_elbo += clusters[k].yjk_hat[j] * log(modelparams.ηk[1]) + (1 - clusters[k].yjk_hat[j]) * log((1-modelparams.ηk[1]))
        perK_data_elbo += one_half_const * clusters[k].yjk_hat[j] * log(clusters[k].v_sq_k_hat[j])
        perK_data_elbo +=  one_half_const * clusters[k].yjk_hat[j]
    end
    

    yjk_entropy_perK += entropy(clusters[k].yjk_hat)
    yjk_entropy_perK = -yjk_entropy_perK
    perK_ebloval =   perK_data_elbo + yjk_entropy_perK
    elbolog.per_k_elbo[k,iter] += perK_ebloval
    data_elbo += perK_ebloval
    end
    #     
    return data_elbo,elbolog
end

# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calc_DataElbo_mpu2(clusters,geneparams,elbolog,dataparams,modelparams,iter)
#     float_type = dataparams.BitType
#     # if isnothing(float_type)
#     #     float_type =eltype(x[1][1])
#     # end
#     G = dataparams.G
#     T = dataparams.T
#     N = dataparams.N
#     K = modelparams.K
#     one_half_const = 1/2
#     λ_sq_vec = [geneparams[j].λ_sq[1] for j in 1:G]
#     data_elbo = 0.0
    
#     @fastmath @inbounds @simd for k in 1:K
#         clusters[k].cache .= 0.0
#         yjk_entropy_perK = 0.0
#         # @fastmath @inbounds @simd for i in 1:N
#         #     clusters[k].cache .+= cellpop[i].rtik[k] .* clusters[k].pip_k .*(cellpop[i].x .-   clusters[k].mk_hat) .^ 2
#         # end
#         # try
#         #     data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)
#         # catch
#         #     println("K =  $k")
#         #     println("σ_sq_k_hat =  $(clusters[k].σ_sq_k_hat)")
#         #     println("λ_sq_vec =  $(λ_sq_vec)")
#         #     println("ηk =  $(modelparams.ηk)")
#         #     println("v_sq_k_hat =  $(clusters[k].v_sq_k_hat)")
#         #     println("var_muk =  $(clusters[k].var_muk )")
#         #     println("mk_hat =  $(clusters[k].mk_hat)")

#         # end

        
#     ####    data_elbo += sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk)) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


#     perK_data_elbo = sum(-one_half_const .* clusters[k].Nk .* log(2π) .- one_half_const .* clusters[k].Nk .* log.(clusters[k].σ_sq_k_hat) .- one_half_const .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat_sq .+  1 ./ clusters[k].σ_sq_k_hat .* clusters[k].x_hat .* clusters[k].κk_hat .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].var_muk .- one_half_const .* clusters[k].Nk .* 1 ./ clusters[k].σ_sq_k_hat .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .- one_half_const .* clusters[k].yjk_hat .* log.(λ_sq_vec) .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].var_muk .-  one_half_const .* 1 ./λ_sq_vec .* clusters[k].yjk_hat .*  (clusters[k].mk_hat) .^2 .+ clusters[k].yjk_hat .* log(modelparams.ηk[1]) .+ (1 .- clusters[k].yjk_hat) .* log((1-modelparams.ηk[1])) .+ one_half_const .* clusters[k].yjk_hat .* log.(clusters[k].v_sq_k_hat ) .+ one_half_const .* clusters[k].yjk_hat)


#     # perK_data_elbo = sum( .+ . .+  .+)

#     # perK_data_elbo = 0.0
#     # for j in 1:G
#     #     perK_data_elbo += -1*one_half_const * clusters[k].Nk[1] * log(2π)
#     #     perK_data_elbo += -1*one_half_const * clusters[k].Nk[1]  * log(clusters[k].σ_sq_k_hat[j])
#     #     perK_data_elbo += -1*one_half_const * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat_sq[j] 
#     #     perK_data_elbo +=  1 / clusters[k].σ_sq_k_hat[j] * clusters[k].x_hat[j] * clusters[k].κk_hat[j]
#     #     perK_data_elbo +=  -1*one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].var_muk[j] 
#     #     perK_data_elbo += -1* one_half_const * clusters[k].Nk[1]  * 1 / clusters[k].σ_sq_k_hat[j] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 - one_half_const * clusters[k].yjk_hat[j] * log(geneparams[j].λ_sq[1]) 
#     #     perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].var_muk[j]
#     #     perK_data_elbo += -1*one_half_const * 1 /geneparams[j].λ_sq[1] * clusters[k].yjk_hat[j] *  (clusters[k].mk_hat[j]) ^2 
#     #     perK_data_elbo += clusters[k].yjk_hat[j] * log(modelparams.ηk[1]) + (1 - clusters[k].yjk_hat[j]) * log((1-modelparams.ηk[1]))
#     #     perK_data_elbo += one_half_const * clusters[k].yjk_hat[j] * log(clusters[k].v_sq_k_hat[j])
#     #     perK_data_elbo +=  one_half_const * clusters[k].yjk_hat[j]
#     # end
    

#     yjk_entropy_perK += entropy(clusters[k].yjk_hat)
#     yjk_entropy_perK = -yjk_entropy_perK
#     perK_ebloval =   perK_data_elbo + yjk_entropy_perK
#     elbolog.per_k_elbo[k,iter] += perK_ebloval
#     data_elbo += perK_ebloval
#     end
#     #     
#     return data_elbo,elbolog
# end

"""
        calc_Hz_fast3(cellpop,clusters,dataparams)
    Fast calculation of the current iterations cell cluster assignment entropy.
"""
function calc_Hz_fast3(cellpop,clusters,dataparams)
    float_type = dataparams.BitType
    N = dataparams.N
    z_entropy = 0.0
    @fastmath @inbounds @simd for i in 1:N
        z_entropy += entropy(cellpop[i].rtik)
    end
    return  -z_entropy
end

# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calc_HyjkSurragateLowerBound_unconstrained(Tk,clusters,dataparams,modelparams)
#     float_type = dataparams.BitType
#     T = dataparams.T
#     K = modelparams.K
#     α0 = modelparams.α0
#     γ0 = modelparams.γ0
#     # rho_hat = sigmoid.(c)
#     # omega_hat = exp.(d)
#     # c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
#     # e_logUk = logUk_expected_value(rho_hat,omega_hat)
#     # e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
#     # K = length(rho_hat)
#     # e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
#     # k_vec = collect(1:K)
#     # # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
#     # lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ e_γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  e_βk .* Tαk[1:K]
#     #take Negative to find max  
#     # lb_lg = sum(lb_lg_k)

#     lb_lg = 0.0
#     yjk_entropy = 0.0
    
#     @fastmath @inbounds @simd for k in 1:K
#         c_B = beta(clusters[k].gk_hat[1] * clusters[k].hk_hat[1] , (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1])
#         e_βk = expectation_βk(k,clusters,modelparams)
#         e_logUk = expectation_logUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
#         e_log1minusUk = expectation_log1minusUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
#         lb_lg += -1.0 * c_B +  (T + 1. - clusters[k].gk_hat[1] * clusters[k].hk_hat[1]) *e_logUk  +  (T * ( K + 1. - k) + γ0 - (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1]) * e_log1minusUk  +  e_βk * α0 * Tk[k]
#         yjk_entropy += entropy(clusters[k].yjk_hat)
#     end
#     yjk_entropy = -yjk_entropy 


#     return lb_lg + yjk_entropy
# end

"""
        calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    Calculates the current iterations surrogate priors' elbo.
"""
function calc_SurragateLowerBound_unconstrained_elbo(Tk,clusters,elbolog,dataparams,modelparams,iter)
    float_type = dataparams.BitType
    T = dataparams.T
    K = modelparams.K
    α0 = modelparams.α0
    γ0 = modelparams.γ0
    # rho_hat = sigmoid.(c)
    # omega_hat = exp.(d)
    # c_B = beta.(rho_hat .* omega_hat , (1.0 .- rho_hat) .* omega_hat)
    # e_logUk = logUk_expected_value(rho_hat,omega_hat)
    # e_log1minusUk =  log1minusUk_expected_value(rho_hat,omega_hat)
    # K = length(rho_hat)
    # e_βk = βk_expected_value(rho_hat,omega_hat)[1:K] # e_βk = βk_expected_value(γ,K)[1:K]
    # k_vec = collect(1:K)
    # # lb_lg_k = [-cB[k] + (T + 1 - rho_hat[k] * omega_hat[k])*e_logUk[k] +  (T*(K+1-k) + γ -(1.0 .- rho_hat[k]) * omega_hat[k])*e_log1minusUk[k] + α_0*e_βk[k]*T_k[k]  for k in 1:K]
    # lb_lg_k = -1.0 .* c_B .+  (T .+ 1. .- rho_hat .* omega_hat) .*e_logUk  .+  (T .* ( K .+ 1. .- k_vec) .+ e_γ .- (1.0 .- rho_hat) .* omega_hat) .* e_log1minusUk  .+  e_βk .* Tαk[1:K]
    #take Negative to find max  
    # lb_lg = sum(lb_lg_k)

    lb_lg = 0.0
    
    @fastmath @inbounds @simd for k in 1:K
        perK_lg_ebloval = 0.0
        c_B = beta(clusters[k].gk_hat[1] * clusters[k].hk_hat[1] , (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1])
        e_βk = expectation_βk(k,clusters,modelparams)
        e_logUk = expectation_logUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
        e_log1minusUk = expectation_log1minusUk(clusters[k].gk_hat[1] , clusters[k].hk_hat[1])
        perK_lg_ebloval += -1.0 * c_B +  (T + 1. - clusters[k].gk_hat[1] * clusters[k].hk_hat[1]) *e_logUk  +  (T * ( K + 1. - k) + γ0 - (1.0 - clusters[k].gk_hat[1]) * clusters[k].hk_hat[1]) * e_log1minusUk  +  e_βk * α0 * Tk[k]
        elbolog.per_k_elbo[k,iter] += perK_lg_ebloval
        lb_lg += perK_lg_ebloval
    end


    return lb_lg,elbolog
end

"""
        calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    Calculates the current iterations dynamic priors' elbo.
"""
function calc_wAllocationsLowerBound(conditionparams,dataparams,modelparams)
    float_type = dataparams.BitType
    # if isnothing(float_type)
    #     float_type =eltype(x[1][1])
    # end
    T = dataparams.T
    # T = length(c_ttprime)
    wAlloc_elbo = 0.0



    @fastmath @inbounds for t in 2:T
        # a_cttprime = 0.0
        # @fastmath @inbounds @simd for t_prime_a in t:T
        #     a_cttprime += conditionparams[t_prime_a].ctt_prime[t]
        # end

        b_cttprime = 0.0
        @fastmath @inbounds for t_prime_b in t:T
            @fastmath @inbounds @simd for l in 1:t-1
                # c_string = "+ c$(t_prime)$(l) "
                b_cttprime += conditionparams[t_prime_b].c_tt_prime[l]
                # sum_string *= c_string
            end
        end
        c_Beta_p = -logbeta(1, modelparams.ϕ0)#c_Beta(adot_w,bdot_w)
        c_Beta_q = -(-logbeta(1, conditionparams[t-1].st_hat[1]))#c_Beta(conditionparams[t-1].awt_hat[1],conditionparams[t-1].bwt_hat[1])
        c_Beta_pq = c_Beta_p + c_Beta_q
        # adot_w_awt_hat_vec = modelparams.adot_w - conditionparams[t-1].awt_hat[1]
        ϕ0_st_hat_vec = modelparams.ϕ0- conditionparams[t-1].st_hat[1]
        e_log_tilde_wt = expectation_log_tilde_wtt(1, conditionparams[t-1].st_hat[1])
        e_log_minus_tilde_wt = expectation_log_minus_tilde_wtt(1, conditionparams[t-1].st_hat[1])


        wAlloc_elbo_t = c_Beta_pq + (1) * e_log_tilde_wt + (ϕ0_st_hat_vec + b_cttprime) * e_log_minus_tilde_wt
        wAlloc_elbo += wAlloc_elbo_t
    end
    return wAlloc_elbo
end

# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calc_alphaElbo_fast3(conditionparams,dataparams,modelparams)
#     float_type = dataparams.BitType
#     T = dataparams.T
#     α_elbo = 0.0
#     @simd for t in 1:T
#         αt_elbo = c_Ga(modelparams.a_α,modelparams.b_α) - c_Ga(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) + (modelparams.a_α - conditionparams[t].a_αt_hat[1]) *expectation_log_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) - (modelparams.b_α -conditionparams[t].b_αt_hat[1]) * expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
#         α_elbo += αt_elbo
#     end
#     return α_elbo
# end


# """
#         get_unique_time_id()
#     This is an example of Docstring. This function receives two 
#     numbers x and y and returns the sum of the squares.
#     ```math

#     ```
# """
# function calc_HsGammaAlphaElbo_fast3(a_γ,b_γ,a_γ_hat,b_γ_hat,conditionparams,dataparams,modelparams)
#     float_type = dataparams.BitType
#     T = dataparams.T
#     α_elbo = 0.0
#     s_entropy = 0.0
#     @fastmath @inbounds @simd for t in 1:T
#         αt_elbo = c_Ga(modelparams.a_α,modelparams.b_α) - c_Ga(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) + (modelparams.a_α - conditionparams[t].a_αt_hat[1]) *expectation_log_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1]) - (modelparams.b_α -conditionparams[t].b_αt_hat[1]) * expectation_αt(conditionparams[t].a_αt_hat[1],conditionparams[t].b_αt_hat[1])
#         α_elbo += αt_elbo
#         s_entropy += entropy(conditionparams[t].c_tt_prime)
#     end

        
#     s_entropy = -s_entropy

#     c_Gamma_γ_p = c_Ga(a_γ,b_γ)
#     c_Gamma_γ_q = -c_Ga(a_γ_hat,b_γ_hat)
#     e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
#     e_log_γ  = log_γ_expected_value(a_γ_hat,b_γ_hat)
#     γ_elbo = c_Gamma_γ_p + c_Gamma_γ_q+(a_γ- a_γ_hat)*e_log_γ - (b_γ-b_γ_hat) * e_γ


#     return α_elbo + γ_elbo + s_entropy
# end

"""
        calc_HsElbo(conditionparams, dataparams, modelparams)
    Calculation of the current iterations conditions assignment entropy.
"""
function calc_HsElbo(conditionparams, dataparams, modelparams)
    float_type = dataparams.BitType
    T = dataparams.T


    s_entropy = 0.0
    @fastmath @inbounds @simd for t in 1:T
        s_entropy += entropy(conditionparams[t].c_tt_prime)
    end


    s_entropy = -s_entropy

    return s_entropy
end


