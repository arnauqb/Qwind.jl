export Model, run!, run_parallel!, do_iteration!, create_models_folders
import Qwind: create_and_run_integrator

using YAML, SharedArrays, Printf, ProgressMeter

mutable struct Model
    config::Dict
    wind_grid::Grid
    bh::BlackHole
    rad::Radiation
    rt::RadiativeTransfer
    ic::InitialConditions
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


function do_iteration!(model::Model, iterations_dict::Dict; it_num)
    if it_num ∉ keys(iterations_dict)
        iterations_dict[1] = Dict()
    end
    save_path = model.config[:integrator][:save_path]
    lines_range, lines_widths = get_initial_radii_and_linewidths(model.ic, model.bh.Rg)
    f(i) = create_and_run_integrator(
        model,
        linewidth = lines_widths[i],
        r0 = lines_range[i],
        atol = model.config[:integrator][:atol],
        rtol = model.config[:integrator][:rtol],
        line_id = i,
    )
    integrators = @showprogress pmap(f, 1:length(lines_range), batch_size=10)
    #integrators_future = Array{Future}(undef, length(lines_range))
    #for (i, (r0, lw)) in enumerate(zip(lines_range, lines_widths))
    #    integrators_future[i] =
    #        #@spawnat (i % nprocs() + 1) create_and_run_integrator(
    #        @spawn create_and_run_integrator(
    #            model,
    #            linewidth = lw,
    #            r0 = r0,
    #            atol = model.config[:integrator][:atol],
    #            rtol = model.config[:integrator][:rtol],
    #            line_id = i,
    #        )
    #end
    #println("done!")
    #integrators = Array{Sundials.IDAIntegrator}(undef, length(lines_range))
    #@sync for i in 1:length(integrators_future)
    #    @async integrators[i] = fetch(integrators_future[i])
    #end
    #integrators[:] .= fetch.(integrators_future)
    iterations_dict[it_num]["integrators"] = integrators
    @info "Integration of iteration $it_num ended!"
    @info "Saving results..."
    flush(stdout)
    save_wind(integrators, model, save_path, it_num)
    @info "Done"
    flush(stdout)
    radiative_transfer = update_radiative_transfer(model.rt, integrators)
    update_model!(model, radiative_transfer)
    iterations_dict[it_num + 1] = Dict()
    iterations_dict[it_num + 1]["radiative_transfer"] = model.rt
    return
end

function run!(model::Model, iterations_dict = nothing)
    if iterations_dict === nothing
        iterations_dict = Dict()
    end
    save_path = model.config[:integrator][:save_path]
    @info "Saving results to $save_path"
    flush(stdout)
    # iterations
    n_iterations = model.config[:integrator][:n_iterations]
    iterations_dict[1] = Dict()
    iterations_dict[1]["radiative_transfer"] = model.rt
    for it = 1:n_iterations
        do_iteration!(model, iterations_dict, it_num = it)
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
        return parse.(Float64, split(value_split[2], ","))
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

function create_running_script(path, n_cpus)
    text = """using Distributed, ClusterManagers
    addprocs_slurm($(n_cpus-1))
    @everywhere using DrWatson
    @everywhere @quickactivate \"Qwind\"
    @everywhere using Qwind, Printf
    args = parse_cl()
    model_num = args[\"model\"]
    model_name = @sprintf(\"model_%03d\", model_num)
    model_path = \"$path\" * \"/\$model_name/config.yaml\"
    model = Model(model_path)
    run!(model)
    """
    open(path * "/run_model.jl", "w") do io
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
    flush(stdout)
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
    create_running_script(save_folder, config[:system][:n_cpus])
    YAML.write_file(save_folder * "/all_configs.yaml", model_dict)
    make_cosma_scripts(length(configs), path = save_folder; configs[1][:system]...)
end

create_models_folders(config::String) =
    create_models_folders(YAML.load_file(config, dicttype = Dict{Symbol,Any}))
