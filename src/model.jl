export Model, run!, run_parallel!, run_iteration!, create_models_folders
import Qwind: create_and_run_integrator

using YAML, Printf, ProgressMeter

mutable struct Model{T<:AbstractFloat}
    config::Dict
    wind_grid::Grid{T}
    bh::BlackHole{T}
    rad::Radiation{T}
    rt::RadiativeTransfer{T}
    ic::InitialConditions{T}
end

function Model(config::Dict)
    bh = BlackHole(config)
    rad = getfield(Qwind, Symbol(config[:radiation][:mode]))(bh, config)
    wind_grid = Rectangular(config)
    ic = getfield(Qwind, Symbol(config[:initial_conditions][:mode]))(rad, bh, config)
    rt = getfield(Qwind, Symbol(config[:radiative_transfer][:mode]))(rad, config)
    save_path = config[:integrator][:save_path]
    #if isdir(save_path)
    #    mv(save_path, save_path * "_backup", force=true)
    #end
    return Model(config, wind_grid, bh, rad, rt, ic)
end

Model(config_path::String) = Model(YAML.load_file(config_path, dicttype = Dict{Symbol,Any}))
update_model!(model::Model, rt::RadiativeTransfer) = model.rt = rt

create_and_run_integrator(model::Model; r0, linewidth, line_id, atol, rtol) =
    create_and_run_integrator(
        model.rt,
        model.wind_grid,
        model.ic,
        r0 = r0,
        linewidth = linewidth,
        line_id = line_id,
        atol = atol,
        rtol = rtol,
        #save_path = save_path,
    )

run_parallel!(config::String, iterations_dict) =
    run_parallel!(YAML.load_file(config, dicttype = Dict{Symbol,Any}), iterations_dict)
run!(config::String, iterations_dict = nothing) = run!(
    YAML.load_file(config, dicttype = Dict{Symbol,Any}),
    iterations_dict = iterations_dict,
)
run!(config::Dict, iterations_dict = nothing) =
    run!(Model(config), iterations_dict = iterations_dict)

initialize_integrators(model::Model) = initialize_integrators(
    model.rt,
    model.wind_grid,
    model.ic,
    atol = model.config[:integrator][:atol],
    rtol = model.config[:integrator][:rtol],
)

function run_integrators!(model::Model, iterations_dict::Dict; it_num, parallel=true)
    lines_range, lines_widths = get_initial_radii_and_linewidths(model)
    f(i) = create_and_run_integrator(
        model,
        linewidth = lines_widths[i],
        r0 = lines_range[i],
        atol = model.config[:integrator][:atol],
        rtol = model.config[:integrator][:rtol],
        line_id = i,
    )
    @info "Starting iteration $it_num ..."
    @info "Integrating $(length(lines_range)) trajectories ..."
    flush()
    if parallel
        if is_logging(stderr)
            # in an hpc, don't show progress bar
            integrators = pmap(f, 1:length(lines_range), batch_size = 10)
        else
            integrators = @showprogress pmap(f, 1:length(lines_range), batch_size = 10)
        end
    else
        integrators = f.(1:length(lines_range))
    end
    iterations_dict[it_num]["integrators"] = integrators
    return integrators
end

function run_iteration!(model::Model, iterations_dict::Dict; it_num, parallel=true)
    if it_num ∉ keys(iterations_dict)
        iterations_dict[1] = Dict()
    end
    save_path = model.config[:integrator][:save_path]
    integrators = run_integrators!(model, iterations_dict, it_num=it_num, parallel=parallel)
    @info "Integration of iteration $it_num ended!"
    @info "Saving results..."
    flush()
    wind_properties = save_wind(integrators, model, save_path, it_num)
    @info "Done"
    flush()
    @info "Wind properties"
    @info "Mass loss fraction $(wind_properties["mass_loss_fraction"])"
    flush()
    radiative_transfer = update_radiative_transfer(model.rt, integrators)
    update_model!(model, radiative_transfer)
    @info "Done"
    flush()
    iterations_dict[it_num + 1] = Dict()
    iterations_dict[it_num + 1]["radiative_transfer"] = model.rt
    return
end

function run!(model::Model, iterations_dict = nothing; start_it=1, n_iterations=nothing, parallel=true)
    if iterations_dict === nothing
        iterations_dict = Dict()
    end
    save_path = model.config[:integrator][:save_path]
    @info "Saving results to $save_path"
    flush()
    # iterations
    if n_iterations === nothing
        n_iterations = model.config[:integrator][:n_iterations]
    end
    iterations_dict[1] = Dict()
    iterations_dict[1]["radiative_transfer"] = model.rt
    for it = start_it:(start_it+n_iterations-1)
        run_iteration!(model, iterations_dict, it_num = it, parallel=parallel)
    end
    return
end


function parse_variation(value)
    value_split = split(value, " ")[2:end]
    if value_split[1] == "linear"
        return range(
            parse(Float64, value_split[2]),
            parse(Float64, value_split[3]),
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
