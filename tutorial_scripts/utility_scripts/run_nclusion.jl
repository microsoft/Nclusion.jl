using nclusion
ENV["GKSwstype"] = "100"
curr_dir = ENV["PWD"]

logger = FormatLogger() do io, args
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end;

function flushed_logger(msg;logger=nothing)
    if !isnothing(logger)
        with_logger(logger) do
            @info msg
        end
    end
end
flushed_logger("Using $(Threads.nthreads()) thread(s)....";logger)
flushed_logger("Loading Packages....";logger)


flushed_logger("Setting Plotting Settings....";logger)
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        flushed_logger("\t Enabling Headless Plotting....";logger)
        
        Gnuplot.options.gpviewer = false
    else
        Gnuplot.options.gpviewer = true
    end
else
    Gnuplot.options.gpviewer = true
end
if "GKSwstype" in collect(keys(ENV)) 
    if ENV["GKSwstype"] == "100"
        flushed_logger("\t Setting Plotting enviroment variables....";logger)
        to_display=false
    else
        to_display=true
    end
else
    to_display=true
end

function main(ARGS)
    datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,num_iter,dataset,save_metrics,outdir = ARGS


    if !isempty(alpha1)
        alpha1 = parse(Float64, alpha1)
    else
        alpha1 = 1.0
    end
    if !isempty(gamma1)
        gamma1 = parse(Float64, gamma1)
    else
        gamma1 = 1.0
    end
    if !isempty(KMax)
        KMax = parse(Int64, KMax)
    else
        KMax = 25
    end

    if !isempty(seed)
        seed = parse(Int64, seed)
    else
        seed = 12345
    end

    if !isempty(elbo_ep)
        elbo_ep = parse(Float64, elbo_ep)
    else
        elbo_ep = 10^(-0)
    end
    if isempty(num_iter)
        num_iter = parse(Float64, num_iter)
    else
        num_iter = 500
    end
    if isempty(dataset)
        dataset = ""
    end
    if isempty(outdir)
        outdir = ""
    end
    if !isempty(save_metrics)
        save_metrics = parse(Bool, save_metrics)
    else
        save_metrics = false
    end


    outputs_dict = run_nclusion(datafilename1,KMax,alpha1,gamma1,seed,elbo_ep,dataset,outdir; logger = logger,num_iter=num_iter,save_metrics=save_metrics)
    filepath = outputs_dict[:filepath]
    filename = "$filepath/output.jld2"
    flushed_logger("Saving Outputs...";logger)
    jldsave(filename,true;outputs_dict=outputs_dict)


    flushed_logger("Finishing Script...";logger)
end
main(ARGS)