using ArgParse, Printf

export searchsorted_nearest,
    searchsorted_first,
    countsignchanges,
    remove_close_elements,
    iter_paths,
    get_value_in_path,
    set_value_in_path!,
    parse_cl,
    make_cosma_scripts,
    d_euclidean

d_euclidean(r0, r1, z0, z1) = sqrt((r0-r1)^2 + (z0-z1)^2)

function searchsorted_nearest(a, x)
    idx = searchsortedfirst(a, x)
    if (idx == 1)
        return idx
    end
    if (idx > length(a))
        return length(a)
    end
    if (a[idx] == x)
        return idx
    end
    if (abs(a[idx] - x) < abs(a[idx - 1] - x))
        return idx
    else
        return idx - 1
    end
end

function searchsorted_first(a, x, direction = 1)
    idx = searchsortedfirst(a, x)
    if (idx > length(a))
        return length(a)
    end
    if a[idx] == x
        return idx
    end
    if (idx == 1)
        return 1
    end
    (direction == 0) && (direction = 1)
    return Int(idx - (direction + 1) / 2)
end

function countsignchanges(array::Vector{Float64}, reference = 0)
    counter = 0
    if length(array) == 0
        return 0
    end
    current_sign = sign(array[1] + reference)
    for elem in array
        if sign(elem + reference) != current_sign
            current_sign = sign(elem + reference)
            counter += 1
        end
    end
    return counter
end

function remove_close_elements(args...; digits = 4)
    ret = [[array[1] for array in args]]
    ret = hcat(ret...)
    for i = 2:length(args[1])
        noinsert = 0
        toinsert = zeros(length(args))
        for (j, array) in enumerate(args)
            element = trunc(array[i], digits = digits)
            if element in ret[j, :]
                noinsert += 1
            end
            toinsert[j] = element
        end
        if noinsert < length(args)
            ret = hcat(ret, toinsert)
        end
    end
    return [ret[i, :] for i = 1:length(args)]
end

function iter_paths(dict)
    function iter(d, path)
        paths = []
        for (k, v) in pairs(d)
            if typeof(v) <: Dict
                paths = vcat(paths, iter(v, vcat(path, k)))
            end
            push!(paths, (vcat(path, k), v))
        end
        return paths
    end
    return iter(dict, [])
end

function get_value_in_path(d, path)
    i = 1
    while i < length(path)
        d = d[path[i]]
        i += 1
    end
    return d[path[end]]
end

function set_value_in_path!(d, path, value)
    i = 1
    while i < length(path)
        if !(path[i] in keys(d))
            d[path[i]] = Dict()
        end
        d = d[path[i]]
        i += 1
    end
    d[path[end]] = value
end


function parse_cl()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--model", "-m"
        help = "Model number to run"
        arg_type = Int
        required = true
    end
    return parsed_args = parse_args(ARGS, s)
end

function cosma_script(;
    n_cpus = 1,
    job_name = "qwind",
    script_number = 1,
    queue = "cosma",
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
    #SBATCH -p $queue
    #SBATCH -A $account
    #SBATCH -o $stdout_path.out
    #SBATCH -e $stdout_path.err
    #SBATCH -t $max_time

    stdbuf -o0 -e0 julia -p $n_cpus $run_script_path -m $script_number
    """
    open(save_path, "w") do io
        write(io, script)
    end
end

function make_cosma_scripts(
    n_models;
    n_cpus = 1,
    path = nothing,
    job_name = "qwind",
    queue = "cosma",
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
            n_cpus = n_cpus,
            job_name = job_name,
            queue = queue,
            account=account,
            max_time = max_time,
            run_script_path = run_script_path,
            save_path = submit_script_path,
            stdout_path = stdout_path,
            script_number=i,
        )
        push!(submit_paths, submit_script_path)
    end
    submit_all_script = join(["sbatch $submit_path" for submit_path in submit_paths], "\n")
    open(path * "/submit_all.sh", "w") do io
        write(io, submit_all_script)
    end
    return
end
