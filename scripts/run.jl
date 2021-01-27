using DrWatson
@quickactivate "Qwind"
using RegionTrees, DataFrames, CSV, YAML, Printf
using Qwind
using Profile
using PProf
using PyPlot
LogNorm = matplotlib.colors.LogNorm

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

function plot_smoothed_grid(smoothed_grid; vmin=1e2, vmax=1e9)
    fig, ax = plt.subplots()
    cm = ax.pcolormesh(
        smoothed_grid.r_range,
        smoothed_grid.z_range,
        smoothed_grid.grid',
        norm = LogNorm(vmin, vmax),
    )
    plt.colorbar(cm, ax = ax)
end

function plot_tau_x_grid(radiative_transfer; vmin=1e2, vmax=1e9)
    fig, ax = plt.subplots()
    r_range = range(0, 1000, length=250)
    z_range = range(0, 1000, length=250)
    tau_x_grid = 1e2 * ones((length(r_range), length(z_range)))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            tau_x_grid[i,j] = compute_xray_tau(radiative_transfer, r, z)
        end
    end
    cm = ax.pcolormesh(
        r_range,
        z_range,
        tau_x_grid',
        norm = LogNorm(vmin, vmax),
    )
    plt.colorbar(cm, ax = ax)
end


config = YAML.load_file("scripts/config.yaml")

black_hole = BlackHole(config)

radiation = @eval $(Symbol(config["radiation"]["mode"]))(black_hole, config)

radiative_transfer =
    @eval $(Symbol(config["radiative_transfer"]["mode"]))(radiation, config)

grid = Rectangular(config)

initial_conditions =
    @eval $(Symbol(config["initial_conditions"]["mode"]))(radiation, black_hole, config)


iterations_dict = Dict()

save_path = config["integrator"]["save_path"]
mkpath(save_path)
# iterations
n_iterations = config["integrator"]["n_iterations"]
for it in 1:n_iterations
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
    @info "Integration of iteration $it ended!"
    radiative_transfer = update_radiative_transfer(radiative_transfer, integrators)
    df = parse_data(integrators);
    output_path = save_path * "/iteration_$(@sprintf "%03d" it).csv";
    CSV.write(output_path, df);
end
