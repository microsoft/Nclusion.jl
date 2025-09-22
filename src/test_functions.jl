
# Define a test function to verify update_w1! function
function test_update_w1!(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <:Int64} 
    # Call the update_w1! function
    updated_conditions = update_w1!(conditions, dataparams, modelparams)
    # Extract the necessary parameters
    phi1 = modelparams.phi1[1]
    float_type = dataparams.BitType
    T = dataparams.T
    I = dataparams.I
    T_all = sum(T)
    # Loop through each time step t and compare to the expected value from the math
    it = 1
    for i in 1:I
        @test updated_conditions[it].w1[1] ≈ 1.0
        it += 1
        for t in 2:T[i]
            expected_w1 = phi1
            index_deltas = collect(1:length(collect(t:T[i]))) .- 1
            for dt_prime in index_deltas
                expected_w1 += conditions[it+dt_prime].Ctt[t]
            end
            # Check that the calculated w1 matches the expected value
            @test updated_conditions[it].w1[1] ≈ expected_w1
            it += 1
        end
    end
    println("All tests passed!")
end
# Define a test function to verify update_w1! function
# ::Vector{CellFeature{U,W,J}} 
function test_update_w1!(cells,conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <:Int64} 
    # Extract the necessary parameters
    phi1 = modelparams.phi1[1]
    float_type = dataparams.BitType
    T = dataparams.T
    I = dataparams.I
    T_all = sum(T)
    # Loop through each time step t and compare to the expected value from the math
    it = 1
    for i in 1:I
        for t in 1:T[i]
            conditions[it].Ctt .= 0.0
            for n in dataparams.TimeRanges[it][1]:dataparams.TimeRanges[it][2]
                conditions[it].Ctt .+= cells[n].c
            end
            it += 1
        end
    end
    # Call the update_w1! function
    updated_conditions = update_w1!(conditions, dataparams, modelparams)
    # Loop through each time step t and compare to the expected value from the math
    it = 1
    for i in 1:I
        @test updated_conditions[it].w1[1] ≈ 1.0
        it += 1
        for t in 2:T[i]
            expected_w1 = phi1
            index_deltas = collect(1:length(collect(t:T[i]))) .- 1
            for dt_prime in index_deltas
                expected_w1 += conditions[it+dt_prime].Ctt[t]
            end
            # Check that the calculated w1 matches the expected value
            @test updated_conditions[it].w1[1] ≈ expected_w1
            it += 1
        end
    end
    println("All tests passed!")
end


# Define a test function to verify update_w2! function
function test_update_w2!(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <:Int64} 
    # Call the update_w2! function
    updated_conditions = update_w2!(conditions, dataparams, modelparams)
    # Extract the necessary parameters
    phi2 = modelparams.phi2[1]
    float_type = dataparams.BitType
    T = dataparams.T
    I = dataparams.I
    T_all = sum(T)
    # Loop through each time step t and compare to the expected value from the math
    it = 1
    for i in 1:I
        @test updated_conditions[it].w2[1] ≈ 1.0
        it += 1
        for t in 2:T[i]
            expected_w2 = phi2
            index_deltas = collect(1:length(collect(t:T[i]))) .- 1
            for dt_prime in index_deltas
                for m in 1:t-1
                    expected_w2 += conditions[it+dt_prime].Ctt[m]
                end
            end
            # Check that the calculated w2 matches the expected value
            @test updated_conditions[it].w2[1] ≈ expected_w2
            it += 1
        end
    end
    println("All tests passed!")
end
# Define a test function to verify update_w2! function
# ::Vector{CellFeature{U,W,J}} 
function test_update_w2!(cells,conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <:Int64}
    # Extract the necessary parameters
    phi2 = modelparams.phi2[1]
    float_type = dataparams.BitType
    T = dataparams.T
    I = dataparams.I
    T_all = sum(T)
    # Loop through each time step t and compare to the expected value from the math
    it = 1
    for i in 1:I
        for t in 1:T[i]
            conditions[it].Ctt .= 0.0
            for n in dataparams.TimeRanges[it][1]:dataparams.TimeRanges[it][2]
                conditions[it].Ctt .+= cells[n].c
            end
            it += 1
        end
    end
    # Call the update_w2! function
    updated_conditions = update_w2!(conditions, dataparams, modelparams)
    # Loop through each time step t and compare to the expected value from the math
    it = 1
    for i in 1:I
        @test updated_conditions[it].w2[1] ≈ 1.0
        it += 1
        for t in 2:T[i]
            expected_w2 = phi2
            index_deltas = collect(1:length(collect(t:T[i]))) .- 1
            for dt_prime in index_deltas
                for m in 1:t-1
                    expected_w2 += conditions[it+dt_prime].Ctt[m]
                end
            end
            # Check that the calculated w2 matches the expected value
            @test updated_conditions[it].w2[1] ≈ expected_w2
            it += 1
        end
    end
    # for t in 2:T
    #     expected_w2 = phi2
    #     for t_prime in t:T
    #         for m in 1:t-1
    #             expected_w2 += conditions[t_prime].Ctt[m]
    #         end
    #     end
    #     # Check that the calculated w2 matches the expected value
    #     @test updated_conditions[t].w2[1] ≈ expected_w2
    # end
    println("All tests passed!")
end



# Define a test function to verify update_d! function
# ::Vector{CellFeature{U,W,J}},
function test_update_d!(cells, clusters::Vector{ClusterFeature{U,W}}, conditions::Vector{ConditionFeature{U,W}}, matrixconditions::Vector{MatrixConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature;from_summary=false,use_log=false) where {U <: AbstractFloat, W <: Int64}
    # Extract the necessary parameters
    T = dataparams.T
    K = modelparams.K
    T_all = sum(T)
    Kplus = K + 1
    # Loop through each time step t and compare to the expected value from the math
    if from_summary
        for it in 1:T_all
            matrixconditions[it].CNtk .= 0.0
            for n in dataparams.TimeRanges[it][1]:dataparams.TimeRanges[it][2]
                matrixconditions[it].CNtk .+= cells[n].c * transpose(cells[n].r)
            end
        end
    else
        for it in 1:T_all
            matrixconditions[it].CNtk .= 0.0
            for n in dataparams.TimeRanges[it][1]:dataparams.TimeRanges[it][2]
                matrixconditions[it].CNtk .+= cells[n].c * transpose(cells[n].r)
            end
        end
    end
    # Call the update_d! function
    updated_conditions = update_d!(clusters, conditions,matrixconditions, dataparams, modelparams;use_log = use_log)
    # Loop through each time step t and compare to the expected value from the math
    it  = 1 
    for i in 1:I
        for t in 1:T[i]
            alpha0 = modelparams.alpha0[i][t]
            expected_d = zeros(Kplus)
            ss_cache = zeros(T[i],Kplus)
            index_deltas = collect(1:length(collect(t:T[i]))) .- 1
            for dt_prime in index_deltas
                ss_cache .+= matrixconditions[it + dt_prime].CNtk
            end
            for k in 1:Kplus
                E_SBk = expectation_SBk(k,clusters,modelparams;use_log = use_log)
                expected_d[k] = alpha0 * E_SBk + ss_cache[t,k]
            end
            # Check that the calculated d matches the expected value
            @test all(updated_conditions[it].d .≈ expected_d)
            it += 1
        end
    end
    println("All tests passed!")
end


# Define a test function to verify update_r! function
# ::Vector{CellFeature{U,W,J}}
function test_update_r!(cells, clusters::Vector{ClusterFeature{U,W}},conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    # Extract the necessary parameters
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    J = dataparams.J
    Kplus = K + 1
    # Call the update_d! function
    updated_cells = update_r!(cells,clusters,conditions,dataparams,modelparams)
    # Loop through each cell i and compare to the expected value from the math
    @inbounds for n in 1:N
        i = cells[n].i#dataparams.LinearAddress[n][1]
        t = cells[n].t#dataparams.LinearAddress[n][2]
        e_log_pi_t_cache = zeros(Kplus)
        for k in 1:Kplus
            conditions_pis_sums = 0.0
            for tt in 1:t
                conditions_pis_sums += cells[n].c[tt] * E_ln_pi(k,conditions[sum(T[1:i-1])+tt])#(digamma(conditionparams[tt].d_hat_t[k]) - digamma(conditionparams[tt].d_hat_t_sum[1])) 
            end
            e_log_pi_t_cache[k] = conditions_pis_sums
        end
        expected_r = zeros(Kplus)
        for k in 1:K
            cells_cache = zeros(J)
            for j in 1:J
                cells_cache[j]  =  -0.5 * E_ln_sigma_sq(j,clusters[k]) - 0.5 * E_one_over_sigma_sq(j,clusters[k]) * E_ll_sq_diff_mu(j,cells[n],clusters[k])
            end
            cell_gene_sums = - 0.5 * dataparams.Jlog
            for el in cells_cache
                cell_gene_sums+=el
            end
            expected_r[k] = e_log_pi_t_cache[k] + cell_gene_sums
        end
        norm_weights3!(K,expected_r)
        if all(updated_cells[n].r .≈ expected_r) == false
            println("updated_cells[n].r .≈ expected_r: ",updated_cells[n].r .≈ expected_r)
            for k in 1:K
                if updated_cells[n].r[k] != expected_r[k]
                    println("updated_cells[n].r[k]: ",updated_cells[n].r[k])
                    println("expected_r[k]: ",expected_r[k])
                end
            end
        end
        @test all(updated_cells[n].r .≈ expected_r)
    end
    println("All tests passed!")
end



# Define a test function to verify update_r! function
# ::Vector{CellFeature{U,W,J}}
function test_update_c!(cells, conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature,modelparams::ModelParameterFeature) where {U <: AbstractFloat, W <: Int64}
    # Extract the necessary parameters
    float_type = dataparams.BitType
    T = dataparams.T
    N = dataparams.N
    K = modelparams.K
    J = dataparams.J
    Kplus = K + 1
    update_c!(cells,conditions,dataparams,modelparams);
    # Call the update_d! function
    # updated_cells = update_c!(cells,conditions,dataparams,modelparams)
    updated_cells,conditions,dataparams,modelparams = update_c!(cells,conditions,dataparams,modelparams)
    # Loop through each cell i and compare to the expected value from the math
    @inbounds for n in 1:N
        i = cells[n].i#dataparams.LinearAddress[n][1]
        t = cells[n].t#dataparams.LinearAddress[n][2]
        available_timepoints = []
        expected_c = zeros(T[i])
        for tt in 1:t
            number_of_summed_terms = []
            for k in 1:Kplus
                append!(number_of_summed_terms,1.0)
                expected_c[tt] += cells[n].r[k]*E_ln_pi(k,conditions[sum(T[1:i-1])+tt])
            end
            append!(available_timepoints,tt)
            available_future_timepoints = []
            append!(number_of_summed_terms,1.0)
            expected_c[tt] += E_ln_omega(conditions[sum(T[1:i-1])+tt])
            for m in tt+1:t
                append!(available_future_timepoints,m)
                append!(number_of_summed_terms,1.0)
                expected_c[tt] +=  E_ln_minusomega(conditions[sum(T[1:i-1])+m])
            end
            if t-(tt+1)<0
                available_future_timepoints_length = 0
            else
                available_future_timepoints_length = t-(tt+1) + 1
            end
            if length(available_future_timepoints) != available_future_timepoints_length
                println("length(available_future_timepoints): ",length(available_future_timepoints))
                println("available_future_timepoints_length: ",available_future_timepoints_length)
                println(t)
                println(tt)
            end
            @test length(available_future_timepoints) == available_future_timepoints_length
            @test length(number_of_summed_terms) == Kplus + available_future_timepoints_length + 1
            # updated_cells[n].c[tt] ≈ sum(number_of_summed_terms)
        end
        @test length(available_timepoints) == t
        norm_weights3!(t,expected_c)
        @test all(updated_cells[n].c .≈ expected_c)
    end
    println("All tests passed!")
end


function test_recursive_cumsum_E_ln_minusomega(conditions::Vector{ConditionFeature{U,W}},dataparams::DataFeature) where {U <: AbstractFloat, W <: Int64}
    T = dataparams.T
    I = length(T)
    for i in 1:I
        for t in 1:T[i]
            for tt in 1:t
                expected = 0.0
                for m in tt+1:t
                    expected += E_ln_minusomega(conditions[sum(T[1:i-1])+m])
                end
                if !isapprox(recursive_cumsum_E_ln_minusomega(i,t,tt+1,T,conditions), expected)
                    println("t: ",t," tt: ",tt)
                end
                @test recursive_cumsum_E_ln_minusomega(i,t,tt+1,T,conditions) ≈ expected
            end
        end
    end
    println("All tests passed!")
end