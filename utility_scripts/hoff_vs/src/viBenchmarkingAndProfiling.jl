function create_results_dict(run_,function_name)
    results_dict = OrderedDict{Symbol,Vector{Union{String,Int,Float64}}}()
    results_dict[:name] = [function_name]
    results_dict[:num_alloc] = [run_.allocs]
    results_dict[:memory] = [run_.memory]
    results_dict[:times] = run_.times
    results_dict[:avg_time] = [mean(run_.times)]
    results_dict[:med_time] = [median(run_.times)]
    results_dict[:max_time] = [maximum(run_.times)]
    results_dict[:min_time] = [minimum(run_.times)]
    results_dict[:num_samples] =[ run_.params.samples]
    results_dict[:num_evals] = [run_.params.evals]
    return results_dict
end
function get_function_name(results_dict)
    return results_dict[:name][1]
end
function run_variational_inference_dynamicHDP_benchmarks(;G=20000,K=50,T=5,N=1_000_000,rand_init=false,nothing_init=false,seed =12345,to_display= true)
    Random.seed!(seed)
    unique_time_id = get_unique_time_id()
    C_t = Int.(N/T .* ones(Int,T));
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = generate_hyperparamters_testing(;rand_init=rand_init)
    λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_init,mk_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init, omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,c_ttprime_init,θ_hat_init,rtik_init = generate_init_vec_testing(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=G,K=K,T=T,N=N,C_t=C_t, nothing_init=nothing_init)
    Glog = generate_Glog_testing(;G=G,rand_init=rand_init)
    x = generate_x_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)
    rtik = generate_rtik_testing(;T=T,K=K,N=N,C_t=C_t,rand_init=rand_init);

    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=rand_init);
    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=rand_init);
    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=rand_init);
    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=rand_init);
    awt_hat_vec = generate_awt_testing(;T=T,rand_init=rand_init);
    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=rand_init);
    a_γ_hat = generate_a_γ_testing(;rand_init=rand_init);
    b_γ_hat = generate_b_γ_testing(;rand_init=rand_init);
    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=rand_init);
    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=rand_init);
    c_hat_vec,d_hat_vec =  StatsFuns.logit.(rhok_hat_vec),log.(omegak_hat_vec)
    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=rand_init);
    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=rand_init);
    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    e_log_π = generate_e_log_π_testing(;T=T,K=K)
    e_log_τ = generate_e_log_τ_testing(;K=K)
    e_τ_μ = generate_e_τ_μ_testing(;K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    Ntk = generate_Ntk_testing(;T=T,K=K)
    Nk = generate_Nk_testing(;K=K,rand_init=rand_init)
    x_hat_k = generate_x_hat_k_testing(;G=G,K=K,rand_init=rand_init)
    x_hat_sq_k = generate_x_hat_sq_k_testing(;G=G,K=K,rand_init=rand_init)
    e_γ = generate_e_γ_testing(;rand_init=rand_init)
    Tαk = generate_Tαk_testing(;K=K,rand_init=rand_init)

    
    
    

    b_init_params_genes = benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    b_init_mk_hat = benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    b_init_λ0k_hat_vec = benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    b_init_a0k_hat_vec = benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    b_init_b0k_hat_vec = benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    b_init_ρωk_hat_vec = benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    b_init_a_γ_hat_vec =benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    b_init_b_γ_hat_vec=benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    b_init_awt_hat_vec = benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    b_init_bwt_hat_vec = benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    b_init_a_αt_hat_vec = benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    b_init_b_αt_hat_vec = benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    b_init_c_ttprime_hat_vec = benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    b_init_θ_hat_vec=benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    b_init_rtik_vec=benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    b_log_π_expected_value = benchmark_log_π_expected_value(θ_hat_vec)
    b_log_τ_k_expected_value = benchmark_log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec)
    b_τ_μ_expected_value = benchmark_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_update_rtik = benchmark_update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
    b_update_Ntk = benchmark_update_Ntk(rtik)
    b_update_c_ttprime = benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    b_update_θ_hat = benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    b_update_Nk = benchmark_update_Nk(rtik)
    b_update_x_hat_k = benchmark_update_x_hat_k(x,rtik)
    b_update_x_hat_sq_k = benchmark_update_x_hat_sq_k(x,rtik)
    b_update_λ0k_hat= benchmark_update_λ0k_hat(λ0_vec,Nk)
    b_update_mk_hat_usingXhat= benchmark_update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
    b_update_a0k_hat_usingXhat =benchmark_update_a0k_hat_usingXhat(a0_vec,Nk)
    b_update_b0k_hat_usingXhat =benchmark_update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)
    b_update_αt = benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    b_update_awt_hat = benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    b_update_bwt_hat = benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    b_update_γ = benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    b_γ_expected_value = benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    b_update_Tαk =benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    b_update_rho_omega_hat = benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    b_calc_DataElbo = benchmark_calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    b_calc_Hz = benchmark_calc_Hz(rtik)
    b_calc_SurragateLowerBound_unconstrained = benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    b_calc_Hs = benchmark_calc_Hs(c_ttprime_vec)
    b_calc_wAllocationsLowerBound = benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    b_calc_GammaElbo = benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    b_calc_alphaElbo = benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)

    parameters = OrderedDict(Symbol("G") => G, Symbol("K") => K,Symbol("N") => N,Symbol("T") => T,Symbol("C_t") => C_t,Symbol("rand_init") => rand_init,Symbol("nothing_init") => nothing_init)
    
    key_list_ = Symbol.(vcat("parameters" ,get_function_name.([b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_calc_DataElbo,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo])))
    var_list_ = [parameters,b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_calc_DataElbo,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo]
    results_dict = OrderedDict{Symbol,OrderedDict}()
    addToDict!(results_dict,key_list_,var_list_);
    

    println("Preparing saving Directory...")
    filepathSummaries = "outputs/Benchmarking-ModelOutputs/"
    model_used = "variational_inference_dynamicHDP_dev"
    modeltype = "Benchmarking $model_used [$unique_time_id]"
    verbal_summary= ". Benchmarking the functions from $model_used, using $G genes, $K maximum clusters, $N cells, and $T condition(s)  "
    inferenceModel_dir = ""
    topDirName = "Benchmarking-ModelOutputs/"
    model_run_id = "dev-$unique_time_id"
    dir_path,topDirPath,dir_name,subdir_name,filepath = makeTidyOutputDirectories(curr_dir, topDirName,unique_time_id)
    println("Plotting Results from Run...")
    #filenamebase=filepath*unique_time_id*"_dataPlotWithLabelsPCA_CalledvsInferred"*".png"
    gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    

    # gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".svg",fig_size=(1100,900));
    # gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".svg",fig_size=(1100,900));
    # gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".svg",fig_size=(1100,900));
    # gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".svg",fig_size=(1100,900));
    # gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".svg",fig_size=(1100,900));

    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)

    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=false)

    println("Saving Summarization File(s)...")
    trunc_filepath = join(split(filepath,"/")[1:end-3],"/") *"/"
    create_benchmark_summarization_file(trunc_filepath,unique_time_id,modeltype,results_dict,verbal_summary)
    println("Saving Results(s)...")
    jldsave(filepath*"results_"*unique_time_id*".jld2"; results_dict=results_dict);
    # b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_m_err_hat,b_init_λ0_err_hat_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_x_hat_error_vs_forloops,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_a0_err_hat_usingXhat12,b_update_λ0_err_hat,b_update_m_err_hat_usingXhat,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_update_v_tikj12,b_calc_DataElbo,b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv
    return results_dict,unique_time_id

end
# gg= run_variational_inference_dynamicHDP_benchmarks(;G=G,K=K,T=T,N=N,rand_init=false,nothing_init=false);

function run_variational_inference_dynamicHDP_vs12_benchmarks(;G=20000,K=50,T=5,N=1_000_000,rand_init=false,nothing_init=false,seed=12345)
    Random.seed!(seed)
    unique_time_id = get_unique_time_id()
    C_t = Int.(N/T .* ones(Int,T));
    
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err = generate_hyperparamters_vs12_testing(;rand_init=rand_init)
    λ0_vec, μ0_vec, a0_vec, b0_vec,λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec,λ0k_hat_init,mk_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init, omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,c_ttprime_init,θ_hat_init,rtik_init,v_tikj_init,λ0_err_hat_vec_init,m_err_hat_init,a0_err_hat_init,b0_err_hat_init = generate_init_vec_vs12_testing(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w,μ0_err,λ0_err,a0_err,b0_err;G=G,K=K,T=T,N=N,C_t=C_t, nothing_init=nothing_init)

    Glog = generate_Glog_testing(;G=G,rand_init=rand_init)
    logpi = generate_logpi_testing(;rand_init=rand_init)
    x = generate_x_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)
    rtik = generate_rtik_testing(;T=T,K=K,N=N,C_t=C_t,rand_init=rand_init);
    v_tikj = generate_v_tikj_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    a0_err_hat_vec = generate_a0_err_testing(;G=G,rand_init=rand_init);
    b0_err_hat_vec = generate_b0_err_testing(;G=G,rand_init=rand_init);
    λ0_err_hat_vec = generate_λ0_err_testing(;G=G,rand_init=rand_init);
    m_err_hat_vec = generate_m_err_testing(;G=G,rand_init=rand_init);
    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=rand_init);
    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=rand_init);
    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=rand_init);
    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=rand_init);
    awt_hat_vec = generate_awt_testing(;T=T,rand_init=rand_init);
    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=rand_init);
    a_γ_hat = generate_a_γ_testing(;rand_init=rand_init);
    b_γ_hat = generate_b_γ_testing(;rand_init=rand_init);
    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=rand_init);
    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=rand_init);
    c_hat_vec,d_hat_vec =  StatsFuns.logit.(rhok_hat_vec),log.(omegak_hat_vec)
    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=rand_init);
    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=rand_init);
    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);

    pct_important=0.5
    ηkj_prior = generate_η_kj_prior_testing(;G=G,K=K,pct_important=pct_important,rand_init=rand_init)

    e_log_π = generate_e_log_π_testing(;T=T,K=K)
    e_log_τkj = generate_e_log_τkj_testing(;G=G,K=K)
    e_log_τj_err = generate_e_log_τkj_err_testing(;G=G,K=K)
    e_τ_μ_tikj = generate_e_τ_μ_tikj_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    e_τ_0j_err = generate_e_τ_0j_err_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)

    Ntk = generate_Ntk_testing(;T=T,K=K)
    N_signal = generate_N_signal_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    N_error = generate_N_error_testing(;G=G,K=K,T=T, N=N ,C_t=C_t,rand_init=rand_init)
    Nj_error = generate_Nj_error_testing(;G=G,rand_init=rand_init)
    Nkj_signal = generate_Nkj_signal_testing(;G=G,K=K,rand_init=rand_init)
    x_hatk_signal = generate_x_hatk_signal_testing(;G=G,K=K,rand_init=rand_init)
    x_hat_sq_err = generate_x_hat_sq_error_testing(;G=G,rand_init=rand_init)
    x_hatk_sq_signal = generate_x_hatk_sq_signal_testing(;G=G,K=K,rand_init=rand_init)
    e_γ = generate_e_γ_testing(;rand_init=rand_init)
    Tαk = generate_Tαk_testing(;K=K,rand_init=rand_init)


    
    

    b_init_params_genes = benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    b_init_params_genes_err = benchmark_init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err)
    b_init_mk_hat = benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    b_init_λ0k_hat_vec = benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    b_init_a0k_hat_vec = benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    b_init_b0k_hat_vec = benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    b_init_ρωk_hat_vec = benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    b_init_a_γ_hat_vec =benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    b_init_b_γ_hat_vec=benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    b_init_awt_hat_vec = benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    b_init_bwt_hat_vec = benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    b_init_a_αt_hat_vec = benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    b_init_b_αt_hat_vec = benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    b_init_c_ttprime_hat_vec = benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    b_init_θ_hat_vec=benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    b_init_rtik_vec=benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    b_init_v_tikj_vec=benchmark_init_v_tikj_vec(v_tikj_init,G,K,T,C_t;rand_init=false)
    b_init_a0_err_hat_vec=benchmark_init_a0_err_hat_vec(a0_err_hat_init,a0_err_vec;rand_init=false)
    b_init_b0_err_hat_vec =benchmark_init_b0_err_hat_vec(b0_err_hat_init,b0_err_vec;rand_init=false)
    

 
    b_log_π_expected_value = benchmark_log_π_expected_value(θ_hat_vec)
    b_τ_μ_expected_value = benchmark_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_log_τ_kj_expected_value = benchmark_log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec)
    b_log_τ_kj_error_expected_value = benchmark_log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
    b_τ_μ_error_expected_value12 = benchmark_τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec)
    b_update_rtik_vs12 = benchmark_update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec)
    b_update_Ntk = benchmark_update_Ntk(rtik)
    b_update_c_ttprime = benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    b_update_θ_hat = benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    b_update_N = benchmark_update_N(rtik,v_tikj)
    b_errorNj12 = benchmark_errorNj12(N_error)
    b_signalNkj = benchmark_signalNkj(N_signal)
    b_update_x_hatk_signal_vs_forloops = benchmark_update_x_hatk_signal_vs_forloops(x,N_signal)
    b_update_x_hat_sq_error_vs_forloops12 = benchmark_update_x_hat_sq_error_vs_forloops12(x,N_error)
    b_update_x_hatk_sq_signal_vs_forloops = benchmark_update_x_hatk_sq_signal_vs_forloops(x,N_signal)

    b_update_a0_err_hat_usingXhat12 = benchmark_update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
    b_update_b0_err_hat_usingXhat12 = benchmark_update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)
    b_update_λ0k_signal_hat= benchmark_update_λ0k_signal_hat(λ0_vec,Nkj_signal)
    b_update_a0k_signal_hat_usingXhat =benchmark_update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
    b_update_mk_signal_hat_usingXhat =benchmark_update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
    b_update_b0k_signal_hat_usingXhat = benchmark_update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
    b_n_τ_μ_expected_value = benchmark_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_n_log_τ_kj_expected_value = benchmark_log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec)
    b_n_log_τ_kj_error_expected_value = benchmark_log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
    b_n_τ_μ_error_expected_value12 = benchmark_τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec)
    b_update_v_tikj12 = benchmark_update_v_tikj12(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,ηkj_prior)
    b_update_αt = benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    b_update_awt_hat = benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    b_update_bwt_hat = benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    b_update_γ = benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    b_γ_expected_value = benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    b_update_Tαk =benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    b_update_rho_omega_hat = benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    b_calc_DataElbo12 = benchmark_calc_DataElbo12(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
    b_calc_Hz = benchmark_calc_Hz(rtik)
    b_calc_SurragateLowerBound_unconstrained = benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    b_calc_Hs = benchmark_calc_Hs(c_ttprime_vec)
    b_calc_wAllocationsLowerBound = benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    b_calc_GammaElbo = benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    b_calc_alphaElbo = benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    b_calc_ImportanceElbo = benchmark_calc_ImportanceElbo(v_tikj,ηkj_prior)
    b_calc_Hv = benchmark_calc_Hv(v_tikj)

    parameters = OrderedDict(Symbol("G") => G, Symbol("K") => K,Symbol("N") => N,Symbol("T") => T,Symbol("C_t") => C_t,Symbol("rand_init") => rand_init,Symbol("nothing_init") => nothing_init)

    # get_function_name.([
    key_list_ = Symbol.(vcat("parameters",get_function_name.([b_init_params_genes]),"init_err_params_genes",get_function_name.([b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_a0_err_hat_usingXhat12,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat]),"new_τ_μ_expected_value","new_og_τ_kj_expected_value","new_log_τ_kj_error_expected_value","new_τ_μ_error_expected_value12",get_function_name.([b_update_v_tikj12,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv])))
    var_list_ =  [parameters, b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_a0_err_hat_usingXhat12,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_n_τ_μ_expected_value,b_n_log_τ_kj_expected_value,b_n_log_τ_kj_error_expected_value,b_n_τ_μ_error_expected_value12,b_update_v_tikj12,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv]
    results_dict = OrderedDict{Symbol,OrderedDict}()
    addToDict!(results_dict,key_list_,var_list_);


    # 
    
    println("Preparing saving Directory...")
    filepathSummaries = "outputs/Benchmarking-ModelOutputs/"
    model_used = "variational_inference_dynamicHDP_vs12"
    modeltype = "Benchmarking $model_used [$unique_time_id]"
    verbal_summary= ". Benchmarking the functions from $model_used, using $G genes, $K maximum clusters, $N cells, and $T condition(s)  "
    inferenceModel_dir = ""
    topDirName = "Benchmarking-ModelOutputs/"
    model_run_id = "vs12-$unique_time_id"
    dir_path,topDirPath,dir_name,subdir_name,filepath = makeTidyOutputDirectories(curr_dir, topDirName,unique_time_id)
    println("Plotting Results from Run...")
    #filenamebase=filepath*unique_time_id*"_dataPlotWithLabelsPCA_CalledvsInferred"*".png"
    gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    

    # gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".svg",fig_size=(1100,900));
    # gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".svg",fig_size=(1100,900));
    # gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".svg",fig_size=(1100,900));
    # gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".svg",fig_size=(1100,900));
    # gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".svg",fig_size=(1100,900));

    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)

    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=false)

    println("Saving Summarization File(s)...")
    trunc_filepath = join(split(filepath,"/")[1:end-3],"/") *"/"
    create_benchmark_summarization_file(trunc_filepath,unique_time_id,modeltype,results_dict,verbal_summary)
    println("Saving Results(s)...")
    jldsave(filepath*"results_"*unique_time_id*".jld2"; unique_time_id = unique_time_id ,results_dict=results_dict);
    
    # b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_m_err_hat,b_init_λ0_err_hat_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_x_hat_error_vs_forloops,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_a0_err_hat_usingXhat12,b_update_λ0_err_hat,b_update_m_err_hat_usingXhat,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_update_v_tikj12,b_calc_DataElbo,b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv
    return results_dict,unique_time_id
end
# hh = run_variational_inference_dynamicHDP_vs12_benchmarks(;G=G,K=K,T=T,N=N,rand_init=false,nothing_init=false);

function run_variational_inference_dynamicHDP_vs18_benchmarks(;G=20000,K=50,T=5,N=1_000_000,rand_init=false,nothing_init=false,seed =12345,to_display= true)
    Random.seed!(seed)
    unique_time_id = get_unique_time_id()
    C_t = Int.(N/T .* ones(Int,T));
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = generate_hyperparamters_testing18(;rand_init=rand_init)
    λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_init,mk_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init, omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,c_ttprime_init,θ_hat_init,rtik_init,pip_kj_init = generate_init_vec_testing18(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=G,K=K,T=T,N=N,C_t=C_t, nothing_init=nothing_init)
    Glog = generate_Glog_testing(;G=G,rand_init=rand_init)
    x = generate_x_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)
    rtik = generate_rtik_testing(;T=T,K=K,N=N,C_t=C_t,rand_init=rand_init);
    pip_kj = generate_pip_kj_imp_weights_testing(;G=G,K=K,rand_init=rand_init)
    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=rand_init);
    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=rand_init);
    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=rand_init);
    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=rand_init);
    awt_hat_vec = generate_awt_testing(;T=T,rand_init=rand_init);
    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=rand_init);
    a_γ_hat = generate_a_γ_testing(;rand_init=rand_init);
    b_γ_hat = generate_b_γ_testing(;rand_init=rand_init);
    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=rand_init);
    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=rand_init);
    c_hat_vec,d_hat_vec =  StatsFuns.logit.(rhok_hat_vec),log.(omegak_hat_vec)
    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=rand_init);
    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=rand_init);
    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    e_log_π = generate_e_log_π_testing(;T=T,K=K)
    e_log_τ = generate_e_log_τ_testing(;K=K)
    e_τ_μ = generate_e_τ_μ_testing(;K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    e_log_τkj = generate_e_log_τkj_testing(;G=G,K=K)
    e_log_τj_err = generate_e_log_τkj_err_testing(;G=G,K=K)
    e_τ_μ_tikj = generate_e_τ_μ_tikj_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    e_τ_0j_err = generate_e_τ_0j_err_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)
    Ntk = generate_Ntk_testing(;T=T,K=K)
    rpip = generate_N_signal_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    Nkj = generate_Nkj_signal_testing(;G=G,K=K,rand_init=rand_init)
    # Nk = generate_Nk_testing(;K=K,rand_init=rand_init)
    x_hat_k = generate_x_hat_k_testing(;G=G,K=K,rand_init=rand_init)
    x_hat_sq_k = generate_x_hat_sq_k_testing(;G=G,K=K,rand_init=rand_init)
    e_γ = generate_e_γ_testing(;rand_init=rand_init)
    Tαk = generate_Tαk_testing(;K=K,rand_init=rand_init)

    
    
    

    b_init_params_genes = benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    b_init_mk_hat = benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    b_init_λ0k_hat_vec = benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    b_init_a0k_hat_vec = benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    b_init_b0k_hat_vec = benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    b_init_ρωk_hat_vec = benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    b_init_a_γ_hat_vec =benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    b_init_b_γ_hat_vec=benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    b_init_awt_hat_vec = benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    b_init_bwt_hat_vec = benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    b_init_a_αt_hat_vec = benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    b_init_b_αt_hat_vec = benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    b_init_c_ttprime_hat_vec = benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    b_init_θ_hat_vec=benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    b_init_rtik_vec=benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    b_init_pip_kj_vec= benchmark_init_pip_kj_vec(pip_kj_init,G,K;rand_init=false)
    b_log_π_expected_value = benchmark_log_π_expected_value(θ_hat_vec)
    b_log_τ_k_expected_value = benchmark_log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec)
    b_τ_μ_expected_value = benchmark_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_update_rtik = benchmark_update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec)
    b_update_Ntk = benchmark_update_Ntk(rtik)
    b_update_c_ttprime = benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    b_update_θ_hat = benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    # b_update_Nk = benchmark_update_Nk(rtik)
    b_update_rpip = benchmark_update_N_rpip18(rtik,pip_kj)
    b_Nkj = benchmark_update_Nkj18(rpip)
    b_update_x_hat_k = benchmark_update_x_hat_k(x,rtik)
    b_update_x_hat_sq_k = benchmark_update_x_hat_sq_k(x,rtik)
    b_update_λ0k_hat= benchmark_update_λ0k_hat(λ0_vec,Nkj)
    b_update_mk_hat_usingXhat= benchmark_update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    b_update_a0k_hat_usingXhat =benchmark_update_a0k_hat_usingXhat18(a0_vec,Nkj)
    b_update_b0k_hat_usingXhat =benchmark_update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    b_update_αt = benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    b_update_awt_hat = benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    b_update_bwt_hat = benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    b_update_γ = benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    b_γ_expected_value = benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    b_update_Tαk =benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    b_update_rho_omega_hat = benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    b_get_gene_PIP= benchmark_get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
    b_calc_DataElbo18 = benchmark_calc_DataElbo18(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    b_calc_Hz = benchmark_calc_Hz(rtik)
    b_calc_SurragateLowerBound_unconstrained = benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    b_calc_Hs = benchmark_calc_Hs(c_ttprime_vec)
    b_calc_wAllocationsLowerBound = benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    b_calc_GammaElbo = benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    b_calc_alphaElbo = benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    b_calc_Hpip = benchmark_calc_Hpip(pip_kj)

    parameters = OrderedDict(Symbol("G") => G, Symbol("K") => K,Symbol("N") => N,Symbol("T") => T,Symbol("C_t") => C_t,Symbol("rand_init") => rand_init,Symbol("nothing_init") => nothing_init)

    key_list_ = Symbol.(vcat("parameters" ,get_function_name.([b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_rpip,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo18,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip])))
    var_list_ = [parameters,b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_rpip,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo18,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip]
    results_dict = OrderedDict{Symbol,OrderedDict}()
    addToDict!(results_dict,key_list_,var_list_);
    

    println("Preparing saving Directory...")
    filepathSummaries = "outputs/Benchmarking-ModelOutputs/"
    model_used = "variational_inference_dynamicHDP_vs18"
    modeltype = "Benchmarking $model_used [$unique_time_id]"
    verbal_summary= ". Benchmarking the functions from $model_used, using $G genes, $K maximum clusters, $N cells, and $T condition(s)  "
    inferenceModel_dir = ""
    topDirName = "Benchmarking-ModelOutputs/"
    model_run_id = "vs18-$unique_time_id"
    dir_path,topDirPath,dir_name,subdir_name,filepath = makeTidyOutputDirectories(curr_dir, topDirName,unique_time_id)
    println("Plotting Results from Run...")
    #filenamebase=filepath*unique_time_id*"_dataPlotWithLabelsPCA_CalledvsInferred"*".png"
    gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".png",fig_size=(1100,900),save_svg_copy=true ,to_display= to_display);
    

    # gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".svg",fig_size=(1100,900));
    # gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".svg",fig_size=(1100,900));
    # gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".svg",fig_size=(1100,900));
    # gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".svg",fig_size=(1100,900));
    # gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".svg",fig_size=(1100,900));

    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true ,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true ,to_display= to_display)

    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=false)

    println("Saving Summarization File(s)...")
    trunc_filepath = join(split(filepath,"/")[1:end-3],"/") *"/"
    create_benchmark_summarization_file(trunc_filepath,unique_time_id,modeltype,results_dict,verbal_summary)
    println("Saving Results(s)...")
    jldsave(filepath*"results_"*unique_time_id*".jld2"; results_dict=results_dict);
    # b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_m_err_hat,b_init_λ0_err_hat_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_x_hat_error_vs_forloops,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_a0_err_hat_usingXhat12,b_update_λ0_err_hat,b_update_m_err_hat_usingXhat,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_update_v_tikj12,b_calc_DataElbo,b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv
    return results_dict,unique_time_id

end


function run_variational_inference_dynamicHDP_vs25_benchmarks(;G=20000,K=50,T=5,N=1_000_000,rand_init=false,nothing_init=false,seed =12345,to_display= true)
    Random.seed!(seed)
    unique_time_id = get_unique_time_id()
    C_t = Int.(N/T .* ones(Int,T));
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = generate_hyperparamters_testing25(;rand_init=rand_init)
    λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_init,mk_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init, omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,c_ttprime_init,θ_hat_init,rtik_init,pip_kj_init = generate_init_vec_testing25(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=G,K=K,T=T,N=N,C_t=C_t, nothing_init=nothing_init)
    Glog = generate_Glog_testing(;G=G,rand_init=rand_init)
    x = generate_x_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)
    rtik = generate_rtik_testing(;T=T,K=K,N=N,C_t=C_t,rand_init=rand_init);
    pip_kj = generate_pip_kj_imp_weights_testing(;G=G,K=K,rand_init=rand_init)
    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=rand_init);
    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=rand_init);
    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=rand_init);
    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=rand_init);
    awt_hat_vec = generate_awt_testing(;T=T,rand_init=rand_init);
    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=rand_init);
    a_γ_hat = generate_a_γ_testing(;rand_init=rand_init);
    b_γ_hat = generate_b_γ_testing(;rand_init=rand_init);
    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=rand_init);
    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=rand_init);
    c_hat_vec,d_hat_vec =  StatsFuns.logit.(rhok_hat_vec),log.(omegak_hat_vec)
    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=rand_init);
    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=rand_init);
    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    e_log_π = generate_e_log_π_testing(;T=T,K=K)
    e_log_τ = generate_e_log_τ_testing(;K=K)
    e_τ_μ = generate_e_τ_μ_testing(;K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    e_log_τkj = generate_e_log_τkj_testing(;G=G,K=K)
    e_log_τj_err = generate_e_log_τkj_err_testing(;G=G,K=K)
    e_τ_μ_tikj = generate_e_τ_μ_tikj_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    e_τ_0j_err = generate_e_τ_0j_err_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init)
    Ntk = generate_Ntk_testing(;T=T,K=K)
    rpip = generate_N_signal_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init)
    Nkj = generate_Nkj_signal_testing(;G=G,K=K,rand_init=rand_init)
    # Nk = generate_Nk_testing(;K=K,rand_init=rand_init)
    x_hat_k = generate_x_hat_k_testing(;G=G,K=K,rand_init=rand_init)
    x_hat_sq_k = generate_x_hat_sq_k_testing(;G=G,K=K,rand_init=rand_init)
    e_γ = generate_e_γ_testing(;rand_init=rand_init)
    Tαk = generate_Tαk_testing(;K=K,rand_init=rand_init)

    
    
    

    b_init_params_genes = benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    b_init_mk_hat = benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    b_init_λ0k_hat_vec = benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    b_init_a0k_hat_vec = benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    b_init_b0k_hat_vec = benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    b_init_ρωk_hat_vec = benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    b_init_a_γ_hat_vec =benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    b_init_b_γ_hat_vec=benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    b_init_awt_hat_vec = benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    b_init_bwt_hat_vec = benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    b_init_a_αt_hat_vec = benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    b_init_b_αt_hat_vec = benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    b_init_c_ttprime_hat_vec = benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    b_init_θ_hat_vec=benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    b_init_rtik_vec=benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    b_init_pip_kj_vec= benchmark_init_pip_kj_vec(pip_kj_init,G,K;rand_init=false)
    b_log_π_expected_value = benchmark_log_π_expected_value25(θ_hat_vec)
    b_log_τ_k_expected_value = benchmark_log_τ_k_expected_value25(a0k_hat_vec, b0k_hat_vec)
    b_τ_μ_expected_value = benchmark_τ_μ_expected_value25(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_update_rtik = benchmark_update_rtik_vs25(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec)
    b_update_Ntk = benchmark_update_Ntk(rtik)
    b_update_c_ttprime = benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    b_update_θ_hat = benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    # b_update_Nk = benchmark_update_Nk(rtik)
    b_update_rpip = benchmark_update_N_rpip25(rtik,pip_kj)
    b_Nkj = benchmark_update_Nkj25(rpip)
    b_update_x_hat_k = benchmark_update_x_hat_k25(x,rtik)
    b_update_x_hat_sq_k = benchmark_update_x_hat_sq_k25(x,rtik)
    b_update_λ0k_hat= benchmark_update_λ0k_hat(λ0_vec,Nkj)
    b_update_mk_hat_usingXhat= benchmark_update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    b_update_a0k_hat_usingXhat =benchmark_update_a0k_hat_usingXhat25(a0_vec,Nkj)
    b_update_b0k_hat_usingXhat =benchmark_update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    b_update_αt = benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    b_update_awt_hat = benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    b_update_bwt_hat = benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    b_update_γ = benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    b_γ_expected_value = benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    b_update_Tαk =benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    b_update_rho_omega_hat = benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    b_get_gene_PIP= benchmark_get_gene_PIP25(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
    b_calc_DataElbo25 = benchmark_calc_DataElbo25(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    b_calc_Hz = benchmark_calc_Hz(rtik)
    b_calc_SurragateLowerBound_unconstrained = benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    b_calc_Hs = benchmark_calc_Hs(c_ttprime_vec)
    b_calc_wAllocationsLowerBound = benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    b_calc_GammaElbo = benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    b_calc_alphaElbo = benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    b_calc_Hpip = benchmark_calc_Hpip(pip_kj)

    parameters = OrderedDict(Symbol("G") => G, Symbol("K") => K,Symbol("N") => N,Symbol("T") => T,Symbol("C_t") => C_t,Symbol("rand_init") => rand_init,Symbol("nothing_init") => nothing_init)

    key_list_ = Symbol.(vcat("parameters" ,get_function_name.([b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_rpip,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo25,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip])))
    var_list_ = [parameters,b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_rpip,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo25,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip]
    results_dict = OrderedDict{Symbol,OrderedDict}()
    addToDict!(results_dict,key_list_,var_list_);
    

    println("Preparing saving Directory...")
    filepathSummaries = "outputs/Benchmarking-ModelOutputs/"
    model_used = "variational_inference_dynamicHDP_vs25"
    modeltype = "Benchmarking $model_used [$unique_time_id]"
    verbal_summary= ". Benchmarking the functions from $model_used, using $G genes, $K maximum clusters, $N cells, and $T condition(s)  "
    inferenceModel_dir = ""
    topDirName = "Benchmarking-ModelOutputs/"
    model_run_id = "vs25-$unique_time_id"
    dir_path,topDirPath,dir_name,subdir_name,filepath = makeTidyOutputDirectories(curr_dir, topDirName,model_run_id)
    println("Plotting Results from Run...")
    #filenamebase=filepath*unique_time_id*"_dataPlotWithLabelsPCA_CalledvsInferred"*".png"
    gnu_plot_memory_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_MemoryUsage"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_AverageTime"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_MedianTime"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_Allocations"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_NumEvals"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    

    # gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".svg",fig_size=(1100,900));
    # gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".svg",fig_size=(1100,900));
    # gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".svg",fig_size=(1100,900));
    # gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".svg",fig_size=(1100,900));
    # gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".svg",fig_size=(1100,900));

    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)

    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=false)

    println("Saving Summarization File(s)...")
    trunc_filepath = join(split(filepath,"/")[1:end-3],"/") *"/"
    create_benchmark_summarization_file(trunc_filepath,model_run_id,modeltype,results_dict,verbal_summary)
    println("Saving Results(s)...")
    jldsave(filepath*"results_"*model_run_id*".jld2"; results_dict=results_dict);
    # b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_m_err_hat,b_init_λ0_err_hat_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_x_hat_error_vs_forloops,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_a0_err_hat_usingXhat12,b_update_λ0_err_hat,b_update_m_err_hat_usingXhat,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_update_v_tikj12,b_calc_DataElbo,b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv
    return results_dict,unique_time_id

end

function run_variational_inference_dynamicHDP_vs25_fast_benchmarks(;G=20000,K=50,T=5,N=1_000_000,rand_init=false,nothing_init=false,seed = 12345,to_display= true)
    Random.seed!(seed)
    unique_time_id = get_unique_time_id()
    C_t = Int.(N/T .* ones(Int,T));
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = generate_hyperparamters_testing25(;rand_init=rand_init)
    λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_init,mk_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init, omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,c_ttprime_init,θ_hat_init,rtik_init,pip_kj_init = generate_init_vec_testing25(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=G,K=K,T=T,N=N,C_t=C_t, nothing_init=nothing_init);
    Glog = generate_Glog_testing(;G=G,rand_init=rand_init);
    x = generate_x_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init);
    rtik = generate_rtik_testing(;T=T,K=K,N=N,C_t=C_t,rand_init=rand_init);
    pip_kj = generate_pip_kj_imp_weights_testing(;G=G,K=K,rand_init=rand_init)
    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=rand_init);
    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=rand_init);
    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=rand_init);
    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=rand_init);
    awt_hat_vec = generate_awt_testing(;T=T,rand_init=rand_init);
    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=rand_init);
    a_γ_hat = generate_a_γ_testing(;rand_init=rand_init);
    b_γ_hat = generate_b_γ_testing(;rand_init=rand_init);
    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=rand_init);
    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=rand_init);
    c_hat_vec,d_hat_vec =  StatsFuns.logit.(rhok_hat_vec),log.(omegak_hat_vec)
    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=rand_init);
    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=rand_init);
    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    e_log_π = generate_e_log_π_testing(;T=T,K=K);
    e_log_τ = generate_e_log_τ_testing(;K=K);
    e_τ_μ = generate_e_τ_μ_testing(;K=K,T=T,N=N,C_t=C_t,rand_init=rand_init);
    e_log_τkj = generate_e_log_τkj_testing(;G=G,K=K);
    e_log_τj_err = generate_e_log_τkj_err_testing(;G=G,K=K);
    e_τ_μ_tikj = generate_e_τ_μ_tikj_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init);
    e_τ_0j_err = generate_e_τ_0j_err_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init);
    Ntk = generate_Ntk_testing(;T=T,K=K);
    rpip = generate_N_signal_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init);
    Nkj = generate_Nkj_signal_testing(;G=G,K=K,rand_init=rand_init);
    # Nk = generate_Nk_testing(;K=K,rand_init=rand_init)
    x_hat_k = generate_x_hat_k_testing(;G=G,K=K,rand_init=rand_init);
    x_hat_sq_k = generate_x_hat_sq_k_testing(;G=G,K=K,rand_init=rand_init);
    e_γ = generate_e_γ_testing(;rand_init=rand_init);
    Tαk = generate_Tαk_testing(;K=K,rand_init=rand_init);

    
    
    

    b_init_params_genes = benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    b_init_mk_hat = benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    b_init_λ0k_hat_vec = benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    b_init_a0k_hat_vec = benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    b_init_b0k_hat_vec = benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    b_init_ρωk_hat_vec = benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    b_init_a_γ_hat_vec =benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    b_init_b_γ_hat_vec=benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    b_init_awt_hat_vec = benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    b_init_bwt_hat_vec = benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    b_init_a_αt_hat_vec = benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    b_init_b_αt_hat_vec = benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    b_init_c_ttprime_hat_vec = benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    b_init_θ_hat_vec=benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    b_init_rtik_vec=benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    b_init_pip_kj_vec= benchmark_init_pip_kj_vec(pip_kj_init,G,K;rand_init=false)
    # b_log_π_expected_value = benchmark_log_π_expected_value25_fast(θ_hat_vec)
    # b_log_τ_k_expected_value = benchmark_log_τ_k_expected_value25_fast(a0k_hat_vec, b0k_hat_vec)
    # b_τ_μ_expected_value = benchmark_τ_μ_expected_value25_fast(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_update_rtik = benchmark_update_rtik_vs25_fast(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec)
    b_update_Ntk = benchmark_update_Ntk(rtik)
    b_update_c_ttprime = benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    b_update_θ_hat = benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    # b_update_Nk = benchmark_update_Nk(rtik)
    # b_update_rpip = benchmark_update_N_rpip25_fast(rtik,pip_kj)
    b_Nkj = benchmark_update_Nkj25_fast(rtik, pip_kj);
    b_update_x_hat_k = benchmark_update_x_hat_k25_fast(x,rtik,pip_kj)
    b_update_x_hat_sq_k = benchmark_update_x_hat_sq_k25_fast(x,rtik,pip_kj)
    b_update_λ0k_hat= benchmark_update_λ0k_hat(λ0_vec,Nkj)
    b_update_mk_hat_usingXhat= benchmark_update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    b_update_a0k_hat_usingXhat =benchmark_update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    b_update_b0k_hat_usingXhat =benchmark_update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    b_update_αt = benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    b_update_awt_hat = benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    b_update_bwt_hat = benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    b_update_γ = benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    b_γ_expected_value = benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    b_update_Tαk =benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    b_update_rho_omega_hat = benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    b_get_gene_PIP= benchmark_get_gene_PIP25_fast(x,rtik,mk_hat_vec,a0k_hat_vec,b0k_hat_vec;null_precision=10)
    b_calc_DataElbo25 = benchmark_calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    b_calc_Hz = benchmark_calc_Hz(rtik)
    b_calc_SurragateLowerBound_unconstrained = benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    b_calc_Hs = benchmark_calc_Hs(c_ttprime_vec)
    b_calc_wAllocationsLowerBound = benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    b_calc_GammaElbo = benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    b_calc_alphaElbo = benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    b_calc_Hpip = benchmark_calc_Hpip(pip_kj)

    parameters = OrderedDict(Symbol("G") => G, Symbol("K") => K,Symbol("N") => N,Symbol("T") => T,Symbol("C_t") => C_t,Symbol("rand_init") => rand_init,Symbol("nothing_init") => nothing_init)

    key_list_ = Symbol.(vcat("parameters" ,get_function_name.([b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo25,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip])))
    var_list_ = [parameters,b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo25,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip]
    results_dict = OrderedDict{Symbol,OrderedDict}()
    addToDict!(results_dict,key_list_,var_list_);
    

    println("Preparing saving Directory...")
    filepathSummaries = "outputs/Benchmarking-ModelOutputs/"
    model_used = "variational_inference_dynamicHDP_vs25_fast"
    modeltype = "Benchmarking $model_used [$unique_time_id]"
    model_run_id = "vs25_fast-$unique_time_id"
    verbal_summary= ". Benchmarking the functions from $model_used, using $G genes, $K maximum clusters, $N cells, and $T condition(s)  "
    inferenceModel_dir = ""
    topDirName = "Benchmarking-ModelOutputs/"
    dir_path,topDirPath,dir_name,subdir_name,filepath = makeTidyOutputDirectories(curr_dir, topDirName,unique_time_id)
    println("Plotting Results from Run...")
    #filenamebase=filepath*unique_time_id*"_dataPlotWithLabelsPCA_CalledvsInferred"*".png"
    gnu_plot_memory_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_MemoryUsage"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_AverageTime"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_MedianTime"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_Allocations"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_NumEvals"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    

    # gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".svg",fig_size=(1100,900));
    # gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".svg",fig_size=(1100,900));
    # gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".svg",fig_size=(1100,900));
    # gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".svg",fig_size=(1100,900));
    # gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".svg",fig_size=(1100,900));

    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)

    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=false)

    println("Saving Summarization File(s)...")
    trunc_filepath = join(split(filepath,"/")[1:end-3],"/") *"/"
    create_benchmark_summarization_file(trunc_filepath,model_run_id,modeltype,results_dict,verbal_summary)
    println("Saving Results(s)...")
    jldsave(filepath*"results_"*model_run_id*".jld2"; results_dict=results_dict);
    # b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_m_err_hat,b_init_λ0_err_hat_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_x_hat_error_vs_forloops,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_a0_err_hat_usingXhat12,b_update_λ0_err_hat,b_update_m_err_hat_usingXhat,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_update_v_tikj12,b_calc_DataElbo,b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv
    return results_dict,unique_time_id

end

function run_variational_inference_dynamicHDP_vs25_fast2_benchmarks(;G=20000,K=50,T=5,N=1_000_000,rand_init=false,nothing_init=false,seed = 12345,to_display= true)
    Random.seed!(seed)
    unique_time_id = get_unique_time_id()
    C_t = Int.(N/T .* ones(Int,T));
    λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w = generate_hyperparamters_testing25(;rand_init=rand_init)
    λ0_vec, μ0_vec, a0_vec, b0_vec,λ0k_hat_init,mk_hat_init,a0k_hat_init,b0k_hat_init,rhok_hat_init, omegak_hat_init,a_γ_hat_init,b_γ_hat_init,awt_hat_init,bwt_hat_init,a_αt_hat_init,b_αt_hat_init,c_ttprime_init,θ_hat_init,rtik_init,pip_kj_init = generate_init_vec_testing25(λ0,μ0,a0,b0,a_α,b_α,a_γ,b_γ,adot_w,bdot_w;G=G,K=K,T=T,N=N,C_t=C_t, nothing_init=nothing_init);
    Glog = generate_Glog_testing(;G=G,rand_init=rand_init);
    x = generate_x_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init);
    rtik = generate_rtik_testing(;T=T,K=K,N=N,C_t=C_t,rand_init=rand_init);
    pip_kj = generate_pip_kj_imp_weights_testing(;G=G,K=K,rand_init=rand_init)
    θ_hat_vec = generate_θ_hat_testing(;T=T,K=K,rand_init=rand_init);
    c_ttprime_vec = generate_c_ttprime_testing(;T=T,rand_init=rand_init);
    a_αt_hat_vec = generate_a_αt_testing(;T=T,rand_init=rand_init);
    b_αt_hat_vec = generate_b_αt_testing(;T=T,rand_init=rand_init);
    awt_hat_vec = generate_awt_testing(;T=T,rand_init=rand_init);
    bwt_hat_vec = generate_bwt_testing(;T=T,rand_init=rand_init);
    a_γ_hat = generate_a_γ_testing(;rand_init=rand_init);
    b_γ_hat = generate_b_γ_testing(;rand_init=rand_init);
    omegak_hat_vec = generate_omegak_testing(;K=K,rand_init=rand_init);
    rhok_hat_vec = generate_rhok_testing(;K=K,rand_init=rand_init);
    c_hat_vec,d_hat_vec =  StatsFuns.logit.(rhok_hat_vec),log.(omegak_hat_vec)
    mk_hat_vec = generate_mk_testing(;G=G,K=K,rand_init=rand_init);
    λ0k_hat_vec = generate_λ0k_testing(;G=G,K=K,rand_init=rand_init);
    b0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    a0k_hat_vec = generate_a0k_testing(;G=G,K=K,rand_init=rand_init);
    e_log_π = generate_e_log_π_testing(;T=T,K=K);
    e_log_τ = generate_e_log_τ_testing(;K=K);
    e_τ_μ = generate_e_τ_μ_testing(;K=K,T=T,N=N,C_t=C_t,rand_init=rand_init);
    e_log_τkj = generate_e_log_τkj_testing(;G=G,K=K);
    e_log_τj_err = generate_e_log_τkj_err_testing(;G=G,K=K);
    e_τ_μ_tikj = generate_e_τ_μ_tikj_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init);
    e_τ_0j_err = generate_e_τ_0j_err_testing(;G=G,T=T,N=N,C_t=C_t,rand_init=rand_init);
    Ntk = generate_Ntk_testing(;T=T,K=K);
    rpip = generate_N_signal_testing(;G=G,K=K,T=T,N=N,C_t=C_t,rand_init=rand_init);
    Nkj = generate_Nkj_signal_testing(;G=G,K=K,rand_init=rand_init);
    # Nk = generate_Nk_testing(;K=K,rand_init=rand_init)
    x_hat_k = generate_x_hat_k_testing(;G=G,K=K,rand_init=rand_init);
    x_hat_sq_k = generate_x_hat_sq_k_testing(;G=G,K=K,rand_init=rand_init);
    e_γ = generate_e_γ_testing(;rand_init=rand_init);
    Tαk = generate_Tαk_testing(;K=K,rand_init=rand_init);

    
    
    

    b_init_params_genes = benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    b_init_mk_hat = benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    b_init_λ0k_hat_vec = benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    b_init_a0k_hat_vec = benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    b_init_b0k_hat_vec = benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    b_init_ρωk_hat_vec = benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    b_init_a_γ_hat_vec =benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    b_init_b_γ_hat_vec=benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    b_init_awt_hat_vec = benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    b_init_bwt_hat_vec = benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    b_init_a_αt_hat_vec = benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    b_init_b_αt_hat_vec = benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    b_init_c_ttprime_hat_vec = benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    b_init_θ_hat_vec=benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    b_init_rtik_vec=benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    b_init_pip_kj_vec= benchmark_init_pip_kj_vec(pip_kj_init,G,K;rand_init=false)
    # b_log_π_expected_value = benchmark_log_π_expected_value25_fast(θ_hat_vec)
    # b_log_τ_k_expected_value = benchmark_log_τ_k_expected_value25_fast(a0k_hat_vec, b0k_hat_vec)
    # b_τ_μ_expected_value = benchmark_τ_μ_expected_value25_fast(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    b_update_rtik = benchmark_update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec)
    b_update_Ntk = benchmark_update_Ntk(rtik)
    b_update_c_ttprime = benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    b_update_θ_hat = benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    # b_update_Nk = benchmark_update_Nk(rtik)
    # b_update_rpip = benchmark_update_N_rpip25_fast(rtik,pip_kj)
    b_Nkj = benchmark_update_Nkj25_fast(rtik, pip_kj);
    b_update_x_hat_k = benchmark_update_x_hat_k25_fast(x,rtik,pip_kj)
    b_update_x_hat_sq_k = benchmark_update_x_hat_sq_k25_fast(x,rtik,pip_kj)
    b_update_λ0k_hat= benchmark_update_λ0k_hat(λ0_vec,Nkj)
    b_update_mk_hat_usingXhat= benchmark_update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nkj,x_hat_k)
    b_update_a0k_hat_usingXhat =benchmark_update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    b_update_b0k_hat_usingXhat =benchmark_update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj,x_hat_k,x_hat_sq_k)
    b_update_αt = benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    b_update_awt_hat = benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    b_update_bwt_hat = benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    b_update_γ = benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    b_γ_expected_value = benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    b_update_Tαk =benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    b_update_rho_omega_hat = benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    b_get_gene_PIP= benchmark_get_gene_PIP25_fast2(x,rtik,mk_hat_vec,a0k_hat_vec,b0k_hat_vec;null_precision=10)
    b_calc_DataElbo25 = benchmark_calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    b_calc_Hz = benchmark_calc_Hz(rtik)
    b_calc_SurragateLowerBound_unconstrained = benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    b_calc_Hs = benchmark_calc_Hs(c_ttprime_vec)
    b_calc_wAllocationsLowerBound = benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    b_calc_GammaElbo = benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    b_calc_alphaElbo = benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    b_calc_Hpip = benchmark_calc_Hpip(pip_kj)

    parameters = OrderedDict(Symbol("G") => G, Symbol("K") => K,Symbol("N") => N,Symbol("T") => T,Symbol("C_t") => C_t,Symbol("rand_init") => rand_init,Symbol("nothing_init") => nothing_init)

    key_list_ = Symbol.(vcat("parameters" ,get_function_name.([b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo25,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip])))
    var_list_ = [parameters,b_init_params_genes,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_pip_kj_vec,b_update_rtik,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_Nkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat,b_get_gene_PIP,b_calc_DataElbo25,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_Hpip]
    results_dict = OrderedDict{Symbol,OrderedDict}()
    addToDict!(results_dict,key_list_,var_list_);
    

    println("Preparing saving Directory...")
    filepathSummaries = "outputs/Benchmarking-ModelOutputs/"
    model_used = "variational_inference_dynamicHDP_vs25_fast2"
    modeltype = "Benchmarking $model_used [$unique_time_id]"
    model_run_id = "vs25_fast2-$unique_time_id"
    verbal_summary= ". Benchmarking the functions from $model_used, using $G genes, $K maximum clusters, $N cells, and $T condition(s)  "
    inferenceModel_dir = ""
    topDirName = "Benchmarking-ModelOutputs/"
    dir_path,topDirPath,dir_name,subdir_name,filepath = makeTidyOutputDirectories(curr_dir, topDirName,unique_time_id)
    println("Plotting Results from Run...")
    #filenamebase=filepath*unique_time_id*"_dataPlotWithLabelsPCA_CalledvsInferred"*".png"
    gnu_plot_memory_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_MemoryUsage"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_AverageTime"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_MedianTime"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_Allocations"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_BarChart_NumEvals"*".png",fig_size=(1100,900),save_svg_copy=true,to_display= to_display);
    

    # gnu_plot_memory_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MemoryUsage"*".svg",fig_size=(1100,900));
    # gnu_plot_avg_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_AverageTime"*".svg",fig_size=(1100,900));
    # gnu_plot_median_time_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_MedianTime"*".svg",fig_size=(1100,900));
    # gnu_plot_allocations_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_Allocations"*".svg",fig_size=(1100,900));
    # gnu_plot_num_evals_barchart(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_BarChart_NumEvals"*".svg",fig_size=(1100,900));

    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_AverageTime"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_Allocations"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=true,save_svg_copy=true,to_display= to_display)
    gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*model_run_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".png",fig_size=(1100,900),plot_prop=false,save_svg_copy=true,to_display= to_display)

    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_avg_time(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_AverageTime"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_alloc(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_Allocations"*".svg",fig_size=(1100,900),plot_prop=false)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=true)
    # gnu_plot_cumulative_memory(results_dict;filenamebase=filepath*unique_time_id*"_BenchmarkResults_CumulPlot_MemoryUsage"*".svg",fig_size=(1100,900),plot_prop=false)

    println("Saving Summarization File(s)...")
    trunc_filepath = join(split(filepath,"/")[1:end-3],"/") *"/"
    create_benchmark_summarization_file(trunc_filepath,model_run_id,modeltype,results_dict,verbal_summary)
    println("Saving Results(s)...")
    jldsave(filepath*"results_"*model_run_id*".jld2"; results_dict=results_dict);
    # b_init_params_genes,b_init_params_genes_err,b_init_mk_hat,b_init_λ0k_hat_vec,b_init_a0k_hat_vec,b_init_b0k_hat_vec,b_init_ρωk_hat_vec,b_init_a_γ_hat_vec,b_init_b_γ_hat_vec,b_init_awt_hat_vec,b_init_bwt_hat_vec,b_init_a_αt_hat_vec,b_init_b_αt_hat_vec,b_init_c_ttprime_hat_vec,b_init_θ_hat_vec,b_init_rtik_vec,b_init_v_tikj_vec,b_init_m_err_hat,b_init_λ0_err_hat_vec,b_init_a0_err_hat_vec,b_init_b0_err_hat_vec,b_log_π_expected_value,b_log_τ_k_expected_value,b_τ_μ_expected_value,b_log_τ_kj_expected_value,b_log_τ_kj_error_expected_value,b_τ_μ_error_expected_value12,b_update_rtik,b_update_rtik_vs12,b_update_Ntk,b_update_c_ttprime,b_update_θ_hat,b_update_Nk,b_update_N,b_errorNj12,b_signalNkj,b_update_x_hat_k,b_update_x_hat_sq_k,b_update_x_hat_error_vs_forloops,b_update_x_hatk_signal_vs_forloops,b_update_x_hat_sq_error_vs_forloops12,b_update_x_hatk_sq_signal_vs_forloops,b_update_λ0k_hat,b_update_mk_hat_usingXhat, b_update_a0k_hat_usingXhat,b_update_b0k_hat_usingXhat,b_update_a0_err_hat_usingXhat12,b_update_λ0_err_hat,b_update_m_err_hat_usingXhat,b_update_b0_err_hat_usingXhat12,b_update_λ0k_signal_hat,b_update_a0k_signal_hat_usingXhat,b_update_mk_signal_hat_usingXhat,b_update_b0k_signal_hat_usingXhat,b_update_αt,b_update_awt_hat,b_update_bwt_hat,b_update_γ,b_γ_expected_value,b_update_Tαk, b_update_rho_omega_hat, b_update_v_tikj12,b_calc_DataElbo,b_calc_DataElbo12,b_calc_Hz,b_calc_SurragateLowerBound_unconstrained,b_calc_Hs,b_calc_wAllocationsLowerBound,b_calc_GammaElbo,b_calc_alphaElbo,b_calc_ImportanceElbo,b_calc_Hv
    return results_dict,unique_time_id

end

# mk_hat_init = init_mk_hat!(mk_hat_init,x,K,μ0_vec;rand_init = rand_init)
# mk_hat_init = init_mk_hat!(mk_hat_init,x,K,μ0_vec;rand_init = rand_init)
function benchmark_init_mk_hat(mk_hat_init,x,K,μ0_vec;rand_init=false)
    function_name = "init_mk_hat!"
    run_ = @benchmark mk_hat_init = init_mk_hat!($mk_hat_init,$x,$K,$μ0_vec;rand_init = $rand_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1)
# λ0k_hat_init = init_λ0k_hat_vec!(λ0k_hat_init,K,λ0_vec;rand_init = rand_init, lo=0,hi=1)
function benchmark_init_λ0k_hat_vec(λ0k_hat_init,K,λ0_vec;rand_init=false)
    function_name = "init_λ0k_hat_vec!"
    run_ = @benchmark λ0k_hat_init = init_λ0k_hat_vec!($λ0k_hat_init,$K,$λ0_vec;rand_init = $rand_init, lo=0,hi=1)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
# a0k_hat_init = init_a0k_hat_vec!(a0k_hat_init,K,a0_vec;rand_init = rand_init, lo=0,hi=1);
function benchmark_init_a0k_hat_vec(a0k_hat_init,K,a0_vec;rand_init=false)
    function_name = "init_a0k_hat_vec!"
    run_ = @benchmark a0k_hat_init = init_a0k_hat_vec!($a0k_hat_init,$K,$a0_vec;rand_init = $rand_init, lo=0,hi=1);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
# b0k_hat_init = init_b0k_hat_vec!(b0k_hat_init,K,b0_vec;rand_init = rand_init, lo=0,hi=1);
function benchmark_init_b0k_hat_vec(b0k_hat_init,K,b0_vec;rand_init=false)
    function_name = "init_b0k_hat_vec!"
    run_ = @benchmark b0k_hat_init = init_b0k_hat_vec!($b0k_hat_init,$K,$b0_vec;rand_init = $rand_init, lo=0,hi=1);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
# rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!(rhok_hat_init,omegak_hat_init,K;rand_init = rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
function benchmark_init_ρωk_hat_vec(rhok_hat_init,omegak_hat_init,K;rand_init=false)
    function_name = "init_ρωk_hat_vec!"
    run_ = @benchmark rhok_hat_init,omegak_hat_init = init_ρωk_hat_vec!($rhok_hat_init,$omegak_hat_init,$K;rand_init = $rand_init, ρ_lo=0,ρ_hi=1, ω_lo= 0,ω_hi = 2);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
# a_γ_hat_init = init_a_γ_hat_vec!(a_γ_hat_init,a_γ;rand_init = rand_init, lo=0,hi=10);
function benchmark_init_a_γ_hat_vec(a_γ_hat_init,a_γ;rand_init=false)
    function_name = "init_a_γ_hat_vec!"
    run_ = @benchmark a_γ_hat_init = init_a_γ_hat_vec!($a_γ_hat_init,$a_γ;rand_init = $rand_init, lo=0,hi=10);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
# b_γ_hat_init = init_b_γ_hat_vec!(b_γ_hat_init,b_γ;rand_init = rand_init, lo=0,hi=10);
function benchmark_init_b_γ_hat_vec(b_γ_hat_init,b_γ;rand_init=false)
    function_name = "init_b_γ_hat_vec!"
    run_ = @benchmark b_γ_hat_init = init_b_γ_hat_vec!($b_γ_hat_init,$b_γ;rand_init = $rand_init, lo=0,hi=10);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
# awt_hat_init = init_awt_hat_vec!(awt_hat_init,T,adot_w;rand_init = rand_init, lo=0,hi=1);
function benchmark_init_awt_hat_vec(awt_hat_init,T,adot_w;rand_init=false)
    function_name = "init_awt_hat_vec!"
    run_ = @benchmark awt_hat_init = init_awt_hat_vec!($awt_hat_init,$T,$adot_w;rand_init = $rand_init, lo=0,hi=1)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
# bwt_hat_init =init_bwt_hat_vec!(bwt_hat_init,T,bdot_w;rand_init = rand_init, lo=0,hi=1);
function benchmark_init_bwt_hat_vec(bwt_hat_init,T,bdot_w;rand_init=false)
    function_name = "init_bwt_hat_vec!"
    run_ = @benchmark bwt_hat_init = init_bwt_hat_vec!($bwt_hat_init,$T,$bdot_w;rand_init = $rand_init, lo=0,hi=1)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
# a_αt_hat_init = init_a_αt_hat_vec!(a_αt_hat_init,T,a_α;rand_init = rand_init, lo=0,hi=10);
function benchmark_init_a_αt_hat_vec(a_αt_hat_init,T,a_α;rand_init=false)
    function_name = "init_a_αt_hat_vec!"
    run_ = @benchmark a_αt_hat_init = init_a_αt_hat_vec!($a_αt_hat_init,$T,$a_α;rand_init = $rand_init, lo=0,hi=10)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
# b_αt_hat_init =  init_b_αt_hat_vec!(b_αt_hat_init,T,b_α;rand_init = rand_init, lo=0,hi=10);
function benchmark_init_b_αt_hat_vec(b_αt_hat_init,T,b_α;rand_init=false)
    function_name = "init_b_αt_hat_vec!"
    run_ = @benchmark b_αt_hat_init = init_b_αt_hat_vec!($b_αt_hat_init,$T,$b_α;rand_init = $rand_init, lo=0,hi=10)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
# c_ttprime_init = init_c_ttprime_hat_vec!(c_ttprime_init,T;rand_init = rand_init);
function benchmark_init_c_ttprime_hat_vec(c_ttprime_init,T;rand_init=false)
    function_name = "init_c_ttprime_hat_vec!"
    run_ = @benchmark c_ttprime_init = init_c_ttprime_hat_vec!($c_ttprime_init,$T;rand_init = $rand_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= omegak_hat_init)
# θ_hat_init = init_θ_hat_vec!(θ_hat_init,K,T;rand_init = rand_init,uniform_theta_init=uniform_theta_init, rhok_hat_init = rhok_hat_init, omegak_hat_init= omegak_hat_init)
function benchmark_init_θ_hat_vec(θ_hat_init,K,T,rhok_hat_init,omegak_hat_init;rand_init=false,uniform_theta_init=true)
    function_name = "init_θ_hat_vec!"
    run_ = @benchmark θ_hat_init = init_θ_hat_vec!($θ_hat_init,$K,$T;rand_init = $rand_init,uniform_theta_init=$uniform_theta_init, rhok_hat_init = $rhok_hat_init, omegak_hat_init= $omegak_hat_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
# rtik_init = init_rtik_vec!(rtik_init,K,T,C_t;rand_init = rand_init)
function benchmark_init_rtik_vec(rtik_init,K,T,C_t;rand_init=false)
    function_name = "init_rtik_vec!"
    run_ = @benchmark rtik_init = init_rtik_vec!($rtik_init,$K,$T,$C_t;rand_init = $rand_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

function benchmark_init_pip_kj_vec(pip_kj_init,G,K;rand_init=false)
    function_name = "init_pip_kj_vec!"
    run_ = @benchmark pip_kj_init = init_pip_kj_vec!($pip_kj_init,$G,$K;rand_init = $rand_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# v_tikj_init = init_v_tikj_vec!(v_tikj_init,G,K,T,C_t;rand_init = rand_init)
function benchmark_init_v_tikj_vec(v_tikj_init,G,K,T,C_t;rand_init=false)
    function_name = "init_v_tikj_vec!"
    run_ = @benchmark v_tikj_init = init_v_tikj_vec!($v_tikj_init,$G,$K,$T,$C_t;rand_init = $rand_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# m_err_hat_init = init_m_err_hat!(m_err_hat_init,x,μ0_err_vec;rand_init = rand_init)
function benchmark_init_m_err_hat(m_err_hat_init,x,μ0_err_vec;rand_init=false)
    function_name = "init_m_err_hat!"
    run_ = @benchmark m_err_hat_init = init_m_err_hat!($m_err_hat_init,$x,$μ0_err_vec;rand_init = $rand_init)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# λ0_err_hat_init = init_λ0_err_hat_vec!(λ0_err_hat_init,λ0_err_vec;rand_init = rand_init, lo=0,hi=1)
function benchmark_init_λ0_err_hat_vec(λ0_err_hat_init,λ0_err_vec;rand_init=false)
    function_name = "init_λ0_err_hat_vec!"
    run_ = @benchmark λ0_err_hat_init = init_λ0_err_hat_vec!($λ0_err_hat_init,$λ0_err_vec;rand_init = $rand_init, lo=0,hi=1)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# a0_err_hat_init = init_a0_err_hat_vec!(a0_err_hat_init,a0_err_vec;rand_init = rand_init, lo=0,hi=1)
function benchmark_init_a0_err_hat_vec(a0_err_hat_init,a0_err_vec;rand_init=false)
    function_name = "init_a0_err_hat_vec!"
    run_ = @benchmark a0_err_hat_init = init_a0_err_hat_vec!($a0_err_hat_init,$a0_err_vec;rand_init = $rand_init, lo=0,hi=1);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# b0_err_hat_init = init_b0_err_hat_vec!(b0_err_hat_init,b0_err_vec;rand_init = rand_init, lo=0,hi=1)
function benchmark_init_b0_err_hat_vec(b0_err_hat_init,b0_err_vec;rand_init=false)
    function_name = "init_b0k_hat_vec!"
    run_ = @benchmark b0_err_hat_init = init_b0_err_hat_vec!($b0_err_hat_init,$b0_err_vec;rand_init = $rand_init, lo=0,hi=1);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes(G,λ0,μ0,a0,b0)
# λ0_err_vec, μ0_err_vec, a0_err_vec, b0_err_vec = init_params_genes(G,λ0_err,μ0_err,a0_err,b0_err);
function benchmark_init_params_genes(G,λ0,μ0,a0,b0)
    function_name = "init_params_genes"
    run_ = @benchmark λ0_vec, μ0_vec, a0_vec, b0_vec = init_params_genes($G,$λ0,$μ0,$a0,$b0);
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_log_π = log_π_expected_value(θ_hat_vec) # T by K
# e_log_π = log_π_expected_value(θ_hat_vec) # T by K
# e_log_π = log_π_expected_value(θ_hat_vec) # T by K
function benchmark_log_π_expected_value(θ_hat_vec)
    function_name = "log_π_expected_value"
    run_ = @benchmark e_log_π = log_π_expected_value($θ_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
# e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
# e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
function benchmark_log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec)
    function_name = "log_τ_k_expected_value"
    run_ = @benchmark e_log_τ = log_τ_k_expected_value($a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
# e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K
# e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
# n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K
# n_e_τ_μ_tikj,_ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec); # T by C_t by K by G and T by C_t by K
function benchmark_τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    function_name = "τ_μ_expected_value"
    run_ = @benchmark e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value($x,$λ0k_hat_vec,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
# e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
# n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
# n_e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
function benchmark_log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec)
    function_name = "log_τ_kj_expected_value"
    run_ = @benchmark e_log_τkj = log_τ_kj_expected_value($a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
# e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
# n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
# n_e_log_τj_err = log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec);
function benchmark_log_τ_kj_error_expected_value(a0_err_hat_vec, b0_err_hat_vec)
    function_name = "log_τ_kj_error_expected_value"
    run_ = @benchmark e_log_τj_err = log_τ_kj_error_expected_value($a0_err_hat_vec, $b0_err_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
# e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
# n_e_τ_0j_err,_  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
# n_e_τ_0j_err,_  = τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec);
function benchmark_τ_μ_error_expected_value12(x, a0_err_hat_vec, b0_err_hat_vec)
    function_name = "τ_μ_error_expected_value12"
    run_ = @benchmark e_τ_0j_err,e_τ_0_err  = τ_μ_error_expected_value12($x, $a0_err_hat_vec, $b0_err_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end




# rtik = update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
function benchmark_update_rtik(Glog,e_log_π,e_log_τ,e_τ_μ,c_ttprime_vec)
    function_name = "update_rtik"
    run_ = @benchmark rtik = update_rtik($Glog,$e_log_π,$e_log_τ,$e_τ_μ,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# rtik = update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);
# rtik = update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec);
function benchmark_update_rtik_vs12(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,v_tikj,c_ttprime_vec)
    function_name = "update_rtik_vs12"
    run_ = @benchmark rtik = update_rtik_vs12($Glog,$e_log_π,$e_log_τkj,$e_τ_μ_tikj,$e_log_τj_err,$e_τ_0j_err,$v_tikj,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# Ntk = update_Ntk(rtik)
# Ntk = update_Ntk(rtik)
# Ntk = update_Ntk(rtik)
function benchmark_update_Ntk(rtik)
    function_name = "update_Ntk"
    run_ = @benchmark Ntk = update_Ntk($rtik)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
# c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
# c_ttprime_vec = update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
function benchmark_update_c_ttprime(awt_hat_vec,bwt_hat_vec,rtik,θ_hat_vec)
    function_name = "update_c_ttprime"
    run_ = @benchmark c_ttprime_vec = update_c_ttprime($awt_hat_vec,$bwt_hat_vec,$rtik,$θ_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
# θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
# θ_hat_vec = update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec) 
function benchmark_update_θ_hat(rhok_hat_vec, omegak_hat_vec,Ntk,a_αt_hat_vec,b_αt_hat_vec,c_ttprime_vec)
    function_name = "update_θ_hat"
    run_ = @benchmark θ_hat_vec = update_θ_hat($rhok_hat_vec, $omegak_hat_vec,$Ntk,$a_αt_hat_vec,$b_αt_hat_vec,$c_ttprime_vec) 
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# Nk = update_Nk(rtik)
function benchmark_update_Nk(rtik)
    function_name = "update_Nk"
    run_ = @benchmark Nk = update_Nk($rtik) 
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# N_signal,N_error = update_N(rtik,v_tikj);
# N_signal,N_error = update_N(rtik,v_tikj);
function benchmark_update_N(rtik,v_tikj)
    function_name = "update_N"
    run_ = @benchmark N_signal,N_error = update_N($rtik,$v_tikj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# Nj_error = update_errorNj12(N_error)
# Nj_error = update_errorNj12(N_error)
function benchmark_errorNj12(N_error)
    function_name = "update_errorNj12"
    run_ = @benchmark Nj_error = update_errorNj12($N_error) 
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# Nkj_signal = update_signalNkj(N_signal)
# Nkj_signal = update_signalNkj(N_signal)
function benchmark_signalNkj(N_signal)
    function_name = "update_signalNkj"
    run_ = @benchmark Nkj_signal = update_signalNkj($N_signal) 
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hat_k = update_x_hat_k(x,rtik)
function benchmark_update_x_hat_k(x,rtik)
    function_name = "update_x_hat_k"
    run_ = @benchmark x_hat_k = update_x_hat_k($x, $rtik)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hat_sq_k = update_x_hat_sq_k(x,rtik)
function benchmark_update_x_hat_sq_k(x,rtik)
    function_name = "update_x_hat_sq_k"
    run_ = @benchmark x_hat_sq_k = update_x_hat_sq_k($x, $rtik)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
# x_hat_err = update_x_hat_error_vs_forloops(x,N_error)
function benchmark_update_x_hat_error_vs_forloops(x,N_error)
    function_name = "update_x_hat_error_vs_forloops"
    run_ = @benchmark x_hat_err = update_x_hat_error_vs_forloops($x, $N_error)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hatk_signal = update_x_hatk_signal_vs_forloops(x,N_signal)
# x_hat_k = update_x_hatk_signal_vs_forloops(x,N_signal)
function benchmark_update_x_hatk_signal_vs_forloops(x,N_signal)
    function_name = "update_x_hatk_signal_vs_forloops"
    run_ = @benchmark x_hatk_signal = update_x_hatk_signal_vs_forloops($x, $N_signal)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hat_sq_err = update_x_hat_sq_error_vs_forloops12(x,N_error)
# x_hat_sq_err = update_x_hat_sq_error_vs_forloops12(x,N_error)
function benchmark_update_x_hat_sq_error_vs_forloops12(x,N_error)
    function_name = "update_x_hat_sq_error_vs_forloops12"
    run_ = @benchmark x_hat_sq_err = update_x_hat_sq_error_vs_forloops12($x, $N_error)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
# x_hat_sq_k = update_x_hatk_sq_signal_vs_forloops(x,N_signal)
function benchmark_update_x_hatk_sq_signal_vs_forloops(x,N_signal)
    function_name = "update_x_hatk_sq_signal_vs_forloops"
    run_ = @benchmark x_hatk_sq_signal = update_x_hatk_sq_signal_vs_forloops($x, $N_signal)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# λ0k_hat_vec = update_λ0k_hat(λ0_vec,Nk)
function benchmark_update_λ0k_hat(λ0_vec,Nk)
    function_name = "update_λ0k_hat"
    run_ = @benchmark λ0k_hat_vec = update_λ0k_hat($λ0_vec,$Nk)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# mk_hat_vec= update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
function benchmark_update_mk_hat_usingXhat(λ0_vec,μ0_vec, Nk,x_hat_k)
    function_name = "update_mk_hat_usingXhat"
    run_ = @benchmark mk_hat_vec= update_mk_hat_usingXhat($λ0_vec,$μ0_vec,$Nk,$x_hat_k)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# a0k_hat_vec = update_a0k_hat_usingXhat(a0_vec,Nk)
function benchmark_update_a0k_hat_usingXhat(a0_vec,Nk)
    function_name = "update_a0k_hat_usingXhat"
    run_ = @benchmark a0k_hat_vec = update_a0k_hat_usingXhat($a0_vec,$Nk)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# b0k_hat_vec = update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)
function benchmark_update_b0k_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nk,x_hat_k,x_hat_sq_k)
    function_name = "update_b0k_hat_usingXhat"
    run_ = @benchmark b0k_hat_vec = update_b0k_hat_usingXhat($b0_vec,$λ0_vec,$μ0_vec,$Nk,$x_hat_k,$x_hat_sq_k)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# a0_err_hat_vec = update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
# a0_err_hat_vec = update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
function benchmark_update_a0_err_hat_usingXhat12(a0_err_vec,Nj_error)
    function_name = "update_a0_err_hat_usingXhat12"
    run_ = @benchmark a0_err_hat_vec = update_a0_err_hat_usingXhat12($a0_err_vec,$Nj_error)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
# λ0_err_hat_vec = update_λ0_err_hat(λ0_err_vec,Nj_error)
function benchmark_update_λ0_err_hat(λ0_err_vec,Nj_error)
    function_name = "update_λ0_err_hat"
    run_ = @benchmark λ0_err_hat_vec = update_λ0_err_hat($λ0_err_vec,$Nj_error)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
# m_err_hat_vec = update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
function benchmark_update_m_err_hat_usingXhat(λ0_err_vec,μ0_err_vec, Nj_error,x_hat_err)
    function_name = "update_m_err_hat_usingXhat"
    run_ = @benchmark m_err_hat_vec = update_m_err_hat_usingXhat($λ0_err_vec,$μ0_err_vec,$Nj_error,$x_hat_err)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# b0_err_hat_vec = update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)
# b0_err_hat_vec = update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)
function benchmark_update_b0_err_hat_usingXhat12(b0_err_vec,x_hat_sq_err)
    function_name = "update_b0_err_hat_usingXhat12"
    run_ = @benchmark b0_err_hat_vec = update_b0_err_hat_usingXhat12($b0_err_vec,$x_hat_sq_err)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
# λ0k_hat_vec = update_λ0k_signal_hat(λ0_vec,Nkj_signal)
function benchmark_update_λ0k_signal_hat(λ0_vec,Nkj_signal)
    function_name = "update_λ0k_signal_hat"
    run_ = @benchmark λ0k_hat_vec = update_λ0k_signal_hat($λ0_vec,$Nkj_signal)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
# a0k_hat_vec = update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
function benchmark_update_a0k_signal_hat_usingXhat(a0_vec,Nkj_signal)
    function_name = "update_a0k_signal_hat_usingXhat"
    run_ = @benchmark a0k_hat_vec = update_a0k_signal_hat_usingXhat($a0_vec,$Nkj_signal)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
# mk_hat_vec= update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hat_k)
function benchmark_update_mk_signal_hat_usingXhat(λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal)
    function_name = "update_mk_signal_hat_usingXhat"
    run_ = @benchmark mk_hat_vec= update_mk_signal_hat_usingXhat($λ0_vec,$μ0_vec,$Nkj_signal,$x_hatk_signal)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
# b0k_hat_vec = update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hat_k,x_hat_sq_k)
function benchmark_update_b0k_signal_hat_usingXhat(b0_vec,λ0_vec,μ0_vec, Nkj_signal,x_hatk_signal,x_hatk_sq_signal)
    function_name = "update_b0k_signal_hat_usingXhat"
    run_ = @benchmark b0k_hat_vec = update_b0k_signal_hat_usingXhat($b0_vec,$λ0_vec,$μ0_vec,$Nkj_signal,$x_hatk_signal,$x_hatk_sq_signal)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
# a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
# a_αt_hat_vec,b_αt_hat_vec = update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
function benchmark_update_αt(a_α,b_α,rhok_hat_vec, omegak_hat_vec,θ_hat_vec)
    function_name = "update_αt"
    run_ = @benchmark a_αt_hat_vec,b_αt_hat_vec = update_αt($a_α,$b_α,$rhok_hat_vec, $omegak_hat_vec,$θ_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
# awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
# awt_hat_vec = update_awt_hat(adot_w, c_ttprime_vec)
function benchmark_update_awt_hat(adot_w, c_ttprime_vec)
    function_name = "update_awt_hat"
    run_ = @benchmark awt_hat_vec = update_awt_hat($adot_w,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
# bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
# bwt_hat_vec = update_bwt_hat(bdot_w, c_ttprime_vec)
function benchmark_update_bwt_hat(bdot_w, c_ttprime_vec)
    function_name = "update_bwt_hat"
    run_ = @benchmark bwt_hat_vec = update_bwt_hat($bdot_w,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
# a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
# a_γ_hat,b_γ_hat = update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
function benchmark_update_γ(a_γ,b_γ,rhok_hat_vec, omegak_hat_vec)
    function_name = "update_γ"
    run_ = @benchmark a_γ_hat,b_γ_hat = update_γ($a_γ,$b_γ,$rhok_hat_vec, $omegak_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_γ = γ_expected_value(a_γ_hat,b_γ_hat) # γ_expected_value(a_γ,b_γ)
# e_γ = γ_expected_value(a_γ,b_γ)
# e_γ = γ_expected_value(a_γ_hat,b_γ_hat)
function benchmark_γ_expected_value(a_γ_hat,b_γ_hat)
    function_name = "γ_expected_value"
    run_ = @benchmark e_γ = γ_expected_value($a_γ_hat,$b_γ_hat)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
# Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
# Tαk = update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)#update_Tk(θ_hat)
function benchmark_update_Tαk(θ_hat_vec,a_αt_hat_vec,b_αt_hat_vec)
    function_name = "update_Tαk"
    run_ = @benchmark Tαk = update_Tαk($θ_hat_vec,$a_αt_hat_vec,$b_αt_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
# rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
# rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
function benchmark_update_rho_omega_hat(rhok_hat_vec, omegak_hat_vec,T,e_γ,Tαk;optim_max_iter=1000)
    function_name = "update_rho_omega_hat"
    run_ = @benchmark rhok_hat_vec, omegak_hat_vec, c_hat_vec,d_hat_vec = update_rho_omega_hat($rhok_hat_vec, $omegak_hat_vec,$T,$e_γ,$Tαk;optim_max_iter=$optim_max_iter)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end



# v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
# v_tikj,_  = update_v_tikj12(Glog,rtik,n_e_log_τkj,n_e_τ_μ_tikj,n_e_log_τj_err,n_e_τ_0j_err,ηkj_prior);
function benchmark_update_v_tikj12(Glog,rtik,e_log_τkj,e_τ_μ_tikj,e_log_τj_err,e_τ_0j_err,ηkj_prior)
    function_name = "update_v_tikj12"
    run_ = @benchmark v_tikj,_  = update_v_tikj12($Glog,$rtik,$e_log_τkj,$e_τ_μ_tikj,$e_log_τj_err,$e_τ_0j_err,$ηkj_prior)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# data_elbo = calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
function benchmark_calc_DataElbo(x,rtik,Nk,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    function_name = "calc_DataElbo"
    run_ = @benchmark data_elbo = calc_DataElbo($x,$rtik,$Nk,$mk_hat_vec,$μ0_vec,$λ0k_hat_vec,$λ0_vec,$a0k_hat_vec,$a0_vec, $b0k_hat_vec,$b0_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# rtik = update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
function benchmark_update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec)
    function_name = "update_rtik_vs18"
    run_ = @benchmark rtik = update_rtik_vs18($Glog,$e_log_π,$e_log_τkj,$e_τ_μ_tikj,$pip_kj,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rpip = update_N_rpip18(rtik,pip_kj)
function benchmark_update_N_rpip18(rtik,pip_kj)
    function_name = "update_N_rpip18"
    run_ = @benchmark rpip = update_N_rpip18($rtik,$pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# Nkj = update_Nkj18(rpip)
function benchmark_update_Nkj18(rpip)
    function_name = "update_Nkj18"
    run_ = @benchmark Nkj = update_Nkj18($rpip) 
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# a0k_hat_vec = update_a0k_hat_usingXhat18(a0_vec,Nkj)
function benchmark_update_a0k_hat_usingXhat18(a0_vec,Nkj)
    function_name = "update_a0k_hat_usingXhat18"
    run_ = @benchmark a0k_hat_vec = update_a0k_hat_usingXhat18($a0_vec,$Nkj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# pip_kj = get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
function benchmark_get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
    function_name = "get_gene_PIP"
    run_ = @benchmark pip_kj = get_gene_PIP($x,$mk_hat_vec,$a0k_hat_vec,$b0k_hat_vec,$rtik;null_precision=$null_precision)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# data_elbo = calc_DataElbo18(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
function benchmark_calc_DataElbo18(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    function_name = "calc_DataElbo18"
    run_ = @benchmark data_elbo = calc_DataElbo18($x,$rpip,$Nkj,$mk_hat_vec,$μ0_vec,$λ0k_hat_vec,$λ0_vec,$a0k_hat_vec,$a0_vec, $b0k_hat_vec,$b0_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# pip_entropy = calc_Hpip(pip_kj);
function benchmark_calc_Hpip(pip_kj)
    function_name = "calc_Hpip"
    run_ = @benchmark assgn_entropy =  calc_Hpip($pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# e_log_π = log_π_expected_value(θ_hat_vec) # T by K
function benchmark_log_π_expected_value25(θ_hat_vec)
    function_name = "log_π_expected_value25"
    run_ = @benchmark e_log_π = log_π_expected_value25($θ_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
function benchmark_log_τ_k_expected_value25(a0k_hat_vec, b0k_hat_vec)
    function_name = "log_τ_k_expected_value25"
    run_ = @benchmark e_log_τ = log_τ_k_expected_value25($a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
function benchmark_τ_μ_expected_value25(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    function_name = "τ_μ_expected_value25"
    run_ = @benchmark e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value25($x,$λ0k_hat_vec,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
function benchmark_log_τ_kj_expected_value25(a0k_hat_vec, b0k_hat_vec)
    function_name = "log_τ_kj_expected_value25"
    run_ = @benchmark e_log_τkj = log_τ_kj_expected_value25($a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# rtik = update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
function benchmark_update_rtik_vs25(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj,pip_kj,c_ttprime_vec)
    function_name = "update_rtik_vs25"
    run_ = @benchmark rtik = update_rtik_vs25($Glog,$e_log_π,$e_log_τkj,$e_τ_μ_tikj,$pip_kj,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# x_hat_k = update_x_hat_k(x,rtik)
function benchmark_update_x_hat_k25(x,rtik)
    function_name = "update_x_hat_k25"
    run_ = @benchmark x_hat_k = update_x_hat_k25($x, $rtik)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hat_sq_k = update_x_hat_sq_k(x,rtik)
function benchmark_update_x_hat_sq_k25(x,rtik)
    function_name = "update_x_hat_sq_k25"
    run_ = @benchmark x_hat_sq_k = update_x_hat_sq_k25($x, $rtik)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rpip = update_N_rpip18(rtik,pip_kj)
function benchmark_update_N_rpip25(rtik,pip_kj)
    function_name = "update_N_rpip25"
    run_ = @benchmark rpip = update_N_rpip25($rtik,$pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# Nkj = update_Nkj18(rpip)
function benchmark_update_Nkj25(rpip)
    function_name = "update_Nkj25"
    run_ = @benchmark Nkj = update_Nkj25($rpip) 
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# a0k_hat_vec = update_a0k_hat_usingXhat18(a0_vec,Nkj)
function benchmark_update_a0k_hat_usingXhat25(a0_vec,Nkj)
    function_name = "update_a0k_hat_usingXhat25"
    run_ = @benchmark a0k_hat_vec = update_a0k_hat_usingXhat25($a0_vec,$Nkj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# pip_kj = get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
function benchmark_get_gene_PIP25(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
    function_name = "get_gene_PIP25"
    run_ = @benchmark pip_kj = get_gene_PIP25($x,$mk_hat_vec,$a0k_hat_vec,$b0k_hat_vec,$rtik;null_precision=$null_precision)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# data_elbo = calc_DataElbo18(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
function benchmark_calc_DataElbo25(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    function_name = "calc_DataElbo25"
    run_ = @benchmark data_elbo = calc_DataElbo18($x,$rpip,$Nkj,$mk_hat_vec,$μ0_vec,$λ0k_hat_vec,$λ0_vec,$a0k_hat_vec,$a0_vec, $b0k_hat_vec,$b0_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end





# e_log_π = log_π_expected_value(θ_hat_vec) # T by K
function benchmark_log_π_expected_value25_fast(θ_hat_vec)
    function_name = "log_π_expected_value25_fast"
    run_ = @benchmark e_log_π = log_π_expected_value25_fast($θ_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# e_log_τ = log_τ_k_expected_value(a0k_hat_vec, b0k_hat_vec) # K by 1
function benchmark_log_τ_k_expected_value25_fast(a0k_hat_vec, b0k_hat_vec)
    function_name = "log_τ_k_expected_value25_fast"
    run_ = @benchmark e_log_τ = log_τ_k_expected_value25_fast($a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec) # T by C_t by K by G and T by C_t by K
function benchmark_τ_μ_expected_value25_fast(x,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec)
    function_name = "τ_μ_expected_value25_fast"
    run_ = @benchmark e_τ_μ_tikj,e_τ_μ = τ_μ_expected_value25_fast($x,$λ0k_hat_vec,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# e_log_τkj = log_τ_kj_expected_value(a0k_hat_vec, b0k_hat_vec);
function benchmark_log_τ_kj_expected_value25_fast(a0k_hat_vec, b0k_hat_vec)
    function_name = "log_τ_kj_expected_value25_fast"
    run_ = @benchmark e_log_τkj = log_τ_kj_expected_value25_fast($a0k_hat_vec, $b0k_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# rtik = update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
function benchmark_update_rtik_vs25_fast(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec)
    function_name = "update_rtik_vs25_fast"
    run_ = @benchmark rtik = update_rtik_vs25_fast($x,$Glog,$θ_hat_vec,$λ0k_hat_vec,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec, $pip_kj,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# x_hat_k = update_x_hat_k(x,rtik)
function benchmark_update_x_hat_k25_fast(x,rtik,pip_kj)
    function_name = "update_x_hat_k25_fast"
    run_ = @benchmark x_hat_k = update_x_hat_k25_fast($x,$rtik, $pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# x_hat_sq_k = update_x_hat_sq_k(x,rtik)
function benchmark_update_x_hat_sq_k25_fast(x,rtik,pip_kj)
    function_name = "update_x_hat_sq_k25_fast"
    run_ = @benchmark x_hat_sq_k = update_x_hat_sq_k25_fast($x,$rtik, $pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rpip = update_N_rpip18(rtik,pip_kj)
function benchmark_update_N_rpip25_fast(rtik,pip_kj)
    function_name = "update_N_rpip25_fast"
    run_ = @benchmark rpip = update_N_rpip25_fast($rtik,$pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# Nkj = update_Nkj18(rpip)
function benchmark_update_Nkj25_fast(rtik,pip_kj)
    function_name = "update_Nkj25_fast"
    run_ = @benchmark Nkj = update_Nkj25_fast($rtik, $pip_kj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# a0k_hat_vec = update_a0k_hat_usingXhat18(a0_vec,Nkj)
function benchmark_update_a0k_hat_usingXhat25_fast(a0_vec,Nkj)
    function_name = "update_a0k_hat_usingXhat25_fast"
    run_ = @benchmark a0k_hat_vec = update_a0k_hat_usingXhat25_fast($a0_vec,$Nkj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# pip_kj = get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
function benchmark_get_gene_PIP25_fast(x,rtik,mk_hat_vec,a0k_hat_vec,b0k_hat_vec;null_precision=10)
    function_name = "get_gene_PIP25_fast"
    run_ = @benchmark pip_kj = get_gene_PIP25_fast($x,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec,$rtik;null_precision=$null_precision)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# data_elbo = calc_DataElbo18(x,rpip,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
function benchmark_calc_DataElbo25_fast(x,rtik,pip_kj,Nkj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec, b0k_hat_vec,b0_vec)
    function_name = "calc_DataElbo25_fast"
    run_ = @benchmark data_elbo = calc_DataElbo25_fast($x,$rtik,$pip_kj,$Nkj,$mk_hat_vec,$μ0_vec,$λ0k_hat_vec,$λ0_vec,$a0k_hat_vec,$a0_vec, $b0k_hat_vec,$b0_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# rtik = update_rtik_vs18(Glog,e_log_π,e_log_τkj,e_τ_μ_tikj, pip_kj,c_ttprime_vec)
function benchmark_update_rtik_vs25_fast2(x,Glog,θ_hat_vec,λ0k_hat_vec,mk_hat_vec,a0k_hat_vec, b0k_hat_vec, pip_kj,c_ttprime_vec)
    function_name = "update_rtik_vs25_fast2"
    run_ = @benchmark rtik = update_rtik_vs25_fast($x,$Glog,$θ_hat_vec,$λ0k_hat_vec,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec, $pip_kj,$c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# pip_kj = get_gene_PIP(x,mk_hat_vec,a0k_hat_vec,b0k_hat_vec,rtik;null_precision=10)
function benchmark_get_gene_PIP25_fast2(x,rtik,mk_hat_vec,a0k_hat_vec,b0k_hat_vec;null_precision=10)
    function_name = "get_gene_PIP25_fast2"
    run_ = @benchmark pip_kj = get_gene_PIP25_fast2($x,$mk_hat_vec,$a0k_hat_vec, $b0k_hat_vec,$rtik;null_precision=$null_precision)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end
# data_elbo = calc_DataElbo12(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
function benchmark_calc_DataElbo12(x,rtik,v_tikj,mk_hat_vec,μ0_vec,λ0k_hat_vec,λ0_vec,a0k_hat_vec,a0_vec,a0_err_vec, a0_err_hat_vec,b0k_hat_vec,b0_vec, b0_err_vec,b0_err_hat_vec)
    function_name = "calc_DataElbo12"
    run_ = @benchmark data_elbo = calc_DataElbo12($x,$rtik,$v_tikj,$mk_hat_vec,$μ0_vec,$λ0k_hat_vec,$λ0_vec,$a0k_hat_vec,$a0_vec,$a0_err_vec, $a0_err_hat_vec,$b0k_hat_vec,$b0_vec, $b0_err_vec,$b0_err_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# assgn_entropy =  calc_Hz(rtik) 
# assgn_entropy =  calc_Hz(rtik) 
function benchmark_calc_Hz(rtik)
    function_name = "calc_Hz"
    run_ = @benchmark assgn_entropy =  calc_Hz($rtik)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
# dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
function benchmark_calc_SurragateLowerBound_unconstrained(c_hat_vec,d_hat_vec,T,e_γ,Tαk)
    function_name = "calc_SurragateLowerBound_unconstrained"
    run_ = @benchmark dHDP_surragate_elbo = calc_SurragateLowerBound_unconstrained($c_hat_vec,$d_hat_vec,$T,$e_γ,$Tαk)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# s_entropy = calc_Hs(c_ttprime_vec)
# s_entropy = calc_Hs(c_ttprime_vec)
function benchmark_calc_Hs(c_ttprime_vec)
    function_name = "calc_Hs"
    run_ = @benchmark s_entropy = calc_Hs($c_ttprime_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
# wAlloc_elbo = calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
function benchmark_calc_wAllocationsLowerBound(c_ttprime_vec, adot_w,bdot_w,awt_hat_vec,bwt_hat_vec)
    function_name = "calc_wAllocationsLowerBound"
    run_ = @benchmark wAlloc_elbo = calc_wAllocationsLowerBound($c_ttprime_vec, $adot_w,$bdot_w,$awt_hat_vec,$bwt_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
# γ_elbo = calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
function benchmark_calc_GammaElbo(a_γ,b_γ,a_γ_hat,b_γ_hat)
    function_name = "calc_GammaElbo"
    run_ = @benchmark γ_elbo = calc_GammaElbo($a_γ,$b_γ,$a_γ_hat,$b_γ_hat)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
# α_elbo = calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
function benchmark_calc_alphaElbo(a_α,b_α,a_αt_hat_vec,b_αt_hat_vec)
    function_name = "calc_alphaElbo"
    run_ = @benchmark α_elbo = calc_alphaElbo($a_α,$b_α,$a_αt_hat_vec,$b_αt_hat_vec)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


# imp_elbo = calc_ImportanceElbo(v_tikj,ηkj_prior)#calc_ImportanceElbo(a_η,b_η,a_ηkj_hat,b_ηkj_hat,v_tikj,e_log_ηkj,e_log_minus_ηkj)
function benchmark_calc_ImportanceElbo(v_tikj,ηkj_prior)
    function_name = "calc_ImportanceElbo"
    run_ = @benchmark imp_elbo = calc_ImportanceElbo($v_tikj,$ηkj_prior)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end

# v_entropy = calc_Hv(v_tikj)
function benchmark_calc_Hv(v_tikj)
    function_name = "calc_Hv"
    run_ = @benchmark v_entropy = calc_Hv($v_tikj)
    results_dict = create_results_dict(run_,function_name)
    return results_dict
end


