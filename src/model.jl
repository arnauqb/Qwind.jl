export Model, run!, run_parallel!, run_iteration!
import Qwind: create_and_run_integrator

using YAML, Printf, ProgressMeter

mutable struct Model{T<:AbstractFloat}
    config::Dict
    wind_grid::Grid{T}
    bh::BlackHole{T}
    rad::Radiation{T}
    ic::InitialConditions{T}
end

function Model(config::Dict)
    bh = BlackHole(config)
    rad = Radiation(bh, config)
    wind_grid = Rectangular(config)
    ic = getfield(Qwind, Symbol(config[:initial_conditions][:mode]))(rad, config)
    save_path = config[:integrator][:save_path]
    return Model(config, wind_grid, bh, rad, ic)
end

Model(config_path::String) = Model(YAML.load_file(config_path, dicttype = Dict{Symbol,Any}))
update_model!(model::Model, radiation::Radiation) = (model.rad = radiation)

run_parallel!(config::String, iterations_dict) =
    run_parallel!(YAML.load_file(config, dicttype = Dict{Symbol,Any}), iterations_dict)
run!(config::String, iterations_dict = nothing) = run!(
    YAML.load_file(config, dicttype = Dict{Symbol,Any}),
    iterations_dict = iterations_dict,
)
run!(config::Dict, iterations_dict = nothing) =
    run!(Model(config), iterations_dict = iterations_dict)

initialize_integrators(model::Model) = initialize_integrators(
    model.rad,
    model.wind_grid,
    model.ic,
    atol = model.config[:integrator][:atol],
    rtol = model.config[:integrator][:rtol],
)

function run_integrators!(model::Model, iterations_dict::Dict; it_num, parallel=true)
    @info "Computing initial conditions for iteration $it_num..."
    lines_range, lines_widths = get_initial_radii_and_linewidths(model)
    f(i) = create_and_run_integrator(
        model,
        linewidth = lines_widths[i],
        r0 = lines_range[i],
        atol = model.config[:integrator][:atol],
        rtol = model.config[:integrator][:rtol],
        trajectory_id= i,
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
    @info "KL fraction $(wind_properties["kinetic_luminosity"] / wind_properties["bolometric_luminosity"])"
    flush()
    new_radiation = update_radiation(model.rad, integrators)
    update_model!(model, new_radiation)
    @info "Done"
    flush()
    iterations_dict[it_num + 1] = Dict()
    iterations_dict[it_num + 1]["rad"] = model.rad
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
    iterations_dict[1]["rad"] = model.rad
    for it = start_it:(start_it+n_iterations-1)
        run_iteration!(model, iterations_dict, it_num = it, parallel=parallel)
    end
    return
end
