using ArgParse, Printf
export create_models_folders

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
    sysimage_path,
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

    julia --sysimage "$sysimage_path" $run_script_path -m $script_number
    """
    open(save_path, "w") do io
        write(io, script)
    end
end

function make_cosma_scripts(
    n_models;
    qwind_path,
    sysimage_path,
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
            sysimage_path=sysimage_path,
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
    submit_all_script = join(["sbatch $submit_path\nsleep 1m" for submit_path in submit_paths], "\n")
    open(path * "/submit_all.sh", "w") do io
        write(io, submit_all_script)
    end
    return
end

function create_running_script(save_path; n_cpus, max_time, account, partition)
    text = """using Distributed, ClusterManagers
    pids = addprocs_slurm($(n_cpus-1),
                          topology=:master_worker,
                          p=\"$partition\",
                          A=\"$account\",
                          t=\"$max_time\",
                          job_file_loc=\"$save_path/cpu_logs\")
    @everywhere pushfirst!(Base.DEPOT_PATH, \"/tmp/julia.cache\")
    println(\"Running on \$(nprocs()) cores.\")
    @everywhere using LinearAlgebra
    @everywhere BLAS.set_num_threads(1)
    println(\"Qwind single node\")
    using Qwind, Printf
    println(\"Done\")
    @info \"Compiling Qwind...\"
    flush(stdout)
    flush(stderr)

    @info \"Done\"
    flush(stdout)
    flush(stderr)

    args = parse_cl()
    model_num = args[\"model\"]
    model_name = @sprintf(\"model_%03d\", model_num)
    model_path = \"$save_path\" * \"/\$model_name/config.yaml\"

    model = Model(model_path)
    run!(model)
    """
    open(save_path * "/run_model.jl", "w") do io
        write(io, text)
    end
end


function parse_variation(value)
    value_split = split(value, " ")[2:end]
    if value_split[1] == "linear"
        return range(
            parse(Float64, value_split[2]),
            parse(Float64, value_split[3]),
            length = parse(Int, value_split[4]),
        )
    elseif value_split[1] == "log"
        return 10 .^ range(
            log10(parse(Float64, value_split[2])),
            log10(parse(Float64, value_split[3])),
            length = parse(Int, value_split[4]),
        )
    elseif value_split[1] == "grid"
        values = nothing
        try
            return parse.(Float64, split(value_split[2], ","))
        catch
            return split(value_split[2], ",")
        end
    else
        error("Type of variation not supported")
    end
end

function parse_configs(config::Dict)
    paths = []
    values = []
    for (path, value) in iter_paths(config)
        if typeof(value) == String
            if occursin("@vary", value)
                push!(paths, path)
                push!(values, parse_variation(value))
            end
        end
    end
    ret = []
    for value_perm in Iterators.product(values...)
        cc = deepcopy(config)
        for (i, value) in enumerate(value_perm)
            path = paths[i]
            set_value_in_path!(cc, path, value)
        end
        push!(ret, cc)
    end
    ret
end

function create_models_folders(config::Dict)
    save_folder = config[:integrator][:save_path]
    save_folder_base = save_folder
    i = 0
    while isdir(save_folder)
        save_folder = save_folder_base * "_$i"
        i += 1
    end
    @info "Saving results to $save_folder"
    flush()
    mkpath(save_folder)
    configs = parse_configs(config)
    model_dict = Dict()
    for (i, config) in enumerate(configs)
        model_name = @sprintf("model_%03d", i)
        model_dict[model_name] = config
        model_folder = save_folder * "/" * model_name
        mkdir(model_folder)
        config[:integrator][:save_path] = model_folder
        YAML.write_file(model_folder * "/config.yaml", config)
    end
    sys_config = config[:system]
    create_running_script(
        save_folder,
        n_cpus = sys_config[:n_cpus],
        partition = sys_config[:partition],
        account = sys_config[:account],
        max_time = sys_config[:max_time],
    )
    YAML.write_file(save_folder * "/all_configs.yaml", model_dict)
    make_cosma_scripts(length(configs), path = save_folder; configs[1][:system]...)
end

create_models_folders(config::String) =
    create_models_folders(YAML.load_file(config, dicttype = Dict{Symbol,Any}))
