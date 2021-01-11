using DrWatson
@quickactivate "Qwind"
using RegionTrees, DataFrames, CSV, YAML, Printf
using Qwind

function parse_data(integrators)
    fields = ["r", "z", "vr", "vz"]
    ret = DataFrame()
    for (i, integrator) in enumerate(integrators)
        df = DataFrame()
        for field in fields
            df[!, :line] = i * ones(Int64, length(integrator.p.data.r))
            df[!, Symbol(field)] = getfield(integrator.p.data, Symbol(field))
        end
        append!(ret, df)
    end
    @info "Results saved!"
    return ret
end

function run_iterations!(
    iterations_dict,
    radiative_transfer,
    grid,
    initial_conditions,
    config,
)
    save_path = config["integrator"]["save_path"]
    mkpath(save_path)
    # iterations
    n_iterations = config["integrator"]["n_iterations"]
    for it = 1:n_iterations
        @info "Starting iteration $it of $n_iterations"
        iterations_dict[it] = Dict()
        integrators = initialize_integrators(
            radiative_transfer,
            grid,
            initial_conditions,
            atol = config["integrator"]["atol"],
            rtol = config["integrator"]["rtol"],
        )
        iterations_dict[it]["integrators"] = integrators
        iterations_dict[it]["radiative_transfer"] = radiative_transfer
        run_integrators!(integrators)
        radiative_transfer = update_radiative_transfer(radiative_transfer, integrators)
        df = parse_data(integrators)
        output_path = save_path * "/iteration_$(@sprintf "%03d" it).csv"
        CSV.write(output_path, df)
    end
end


config = YAML.load_file("scripts/config.yaml")

black_hole = BlackHole(config)

radiation = @eval $(Symbol(config["radiation"]["mode"]))(black_hole, config)

radiative_transfer =
    @eval $(Symbol(config["radiative_transfer"]["mode"]))(radiation, config)

grid = Grid(config)

initial_conditions =
    @eval $(Symbol(config["initial_conditions"]["mode"]))(radiation, black_hole, config)


iterations_dict = Dict()
iterations_dict =
    run_iterations!(iterations_dict, radiative_transfer, grid, initial_conditions, config)
