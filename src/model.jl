export Model, run!, run_parallel!, run_iteration!
import Qwind: create_and_run_integrator

using YAML, Printf, ProgressMeter

mutable struct Model{T<:AbstractFloat,U<:Int,V<:Bool,S<:Flag,W<:String}
    parameters::Parameters{T,S,V,U,W}
    bh::BlackHole{T}
    wind::Wind{T,U,V}
    rad::Radiation{T}
    ic::InitialConditions{T}
end

function Model(config::Dict)
    parameters = Parameters(config)
    bh = BlackHole(parameters)
    rad = Radiation(bh, parameters)
    wind = Wind(parameters)
    ic = get_initial_conditions(rad, parameters)
    return Model(parameters, bh, wind, rad, ic)
end

Model(config_path::String) = Model(YAML.load_file(config_path, dicttype = Dict{Symbol,Any}))
update_wind!(model::Model, wind::Wind) = (model.wind = wind)

run_parallel!(config::String, iterations_dict) =
    run_parallel!(YAML.load_file(config, dicttype = Dict{Symbol,Any}), iterations_dict)
run!(config::String, iterations_dict = nothing; parallel = true) = run!(
    YAML.load_file(config, dicttype = Dict{Symbol,Any}),
    iterations_dict = iterations_dict,
    parallel = parallel,
)
run!(config::Dict, iterations_dict = nothing; parallel = true) =
    run!(Model(config), iterations_dict = iterations_dict, parallel = parallel)

initialize_integrators(model::Model) = initialize_integrators(
    model.rad,
    model.wind,
    model.ic,
    atol = model.parameters.integrator_atol,
    rtol = model.parameters.integrator_rtol,
)

function run_integrators!(model::Model, iterations_dict::Dict; it_num, parallel = true)
    @info "Computing initial conditions for iteration $it_num..."
    lines_range, lines_widths = get_initial_radii_and_linewidths(model)
    f(i) = create_and_run_integrator(
        model,
        linewidth = lines_widths[i],
        r0 = lines_range[i],
        trajectory_id = i,
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
    # Resolve intersections.
    max_times = get_intersection_times(integrators)
    streamlines = interpolate_integrators(
        integrators,
        max_times = max_times,
        n_timesteps = 1000,
        log = true,
    )
    iterations_dict[it_num]["integrators"] = integrators
    iterations_dict[it_num]["streamlines"] = streamlines
    return integrators, streamlines
end

function run_iteration!(model::Model, iterations_dict::Dict; it_num, parallel = true)
    if it_num ∉ keys(iterations_dict)
        iterations_dict[1] = Dict()
    end
    save_path = model.parameters.save_path
    integrators, streamlines =
        run_integrators!(model, iterations_dict, it_num = it_num, parallel = parallel)
    @info "Integration of iteration $it_num ended!"
    @info "Saving results..."
    flush()
    wind_properties = save_wind(integrators, streamlines, model, save_path, it_num)
    @info "Done"
    flush()
    @info "Wind properties"
    @info "Mass loss fraction $(wind_properties["mass_loss_fraction"])"
    @info "KL fraction $(wind_properties["kinetic_luminosity"] / wind_properties["bolometric_luminosity"])"
    flush()
    new_wind = update_wind(model.wind, streamlines)
    update_wind!(model, new_wind)
    @info "Done"
    flush()
    iterations_dict[it_num + 1] = Dict()
    iterations_dict[it_num + 1]["rad"] = model.rad
    return
end

function run!(
    model::Model,
    iterations_dict = nothing;
    start_it = 1,
    n_iterations = nothing,
    parallel = true,
)
    if iterations_dict === nothing
        iterations_dict = Dict()
    end
    save_path = model.parameters.save_path
    @info "Saving results to $save_path"
    flush()
    # iterations
    if n_iterations === nothing
        n_iterations = model.parameters.n_iterations
    end
    iterations_dict[1] = Dict()
    iterations_dict[1]["rad"] = model.rad
    for it = start_it:(start_it + n_iterations - 1)
        run_iteration!(model, iterations_dict, it_num = it, parallel = parallel)
    end
    return
end
