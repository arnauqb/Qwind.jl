using ArgParse, Printf

function parse_cl()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--model", "-m"
        help = "Model number to run"
        arg_type = Int
        required = true
    end
    return parsed_args = parse_args(ARGS, s)
end

function cosma_script(;
    qwind_path,
    n_cpus = 1,
    job_name = "qwind",
    script_number = 1,
    partition = "cosma",
    account="durham",
    stdout_path = nothing,
    max_time = "72:00:00",
    save_path = nothing,
    run_script_path = nothing,
)
    script = """
    #!/bin/bash -l
    #SBATCH --ntasks $n_cpus
    #SBATCH -J $(job_name[1:4])_$(@sprintf "%03d" script_number)
    #SBATCH -p $partition
    #SBATCH -A $account
    #SBATCH -o $stdout_path.out
    #SBATCH -e $stdout_path.err
    #SBATCH -t $max_time
    
    export JULIA_PROJECT=$qwind_path


    julia --project=$qwind_path -e 'push!(LOAD_PATH, "@pkglock"); using PkgLock; PkgLock.instantiate_precompile()'
    julia $run_script_path -m $script_number
    """
    open(save_path, "w") do io
        write(io, script)
    end
end

function make_cosma_scripts(
    n_models;
    qwind_path,
    n_cpus = 1,
    path = nothing,
    job_name = "qwind",
    partition = "cosma",
    account = "durham",
    max_time = 72,
)
    run_script_path = path * "/run_model.jl"
    mkpath(path * "/stdout")
    submit_paths = []
    for i = 1:n_models
        model_name = @sprintf "model_%03d" i
        stdout_path = path * "/stdout/$model_name"
        submit_script_path = path * "/$model_name/submit.sh"
        cosma_script(
            qwind_path=qwind_path,
            n_cpus = n_cpus,
            job_name = job_name,
            partition = partition,
            account=account,
            max_time = max_time,
            run_script_path = run_script_path,
            save_path = submit_script_path,
            stdout_path = stdout_path,
            script_number=i,
        )
        push!(submit_paths, submit_script_path)
    end
    submit_all_script = join(["sbatch $submit_path\nsleep 10m" for submit_path in submit_paths], "\n")
    open(path * "/submit_all.sh", "w") do io
        write(io, submit_all_script)
    end
    return
end
