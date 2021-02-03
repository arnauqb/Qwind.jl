using DrWatson
@quickactivate "Qwind"
using RegionTrees, DataFrames, CSV, YAML, Printf
using Qwind
using Profile
using PProf
using PyPlot
LogNorm = matplotlib.colors.LogNorm

function plot_smoothed_grid(smoothed_grid; vmin = 1e2, vmax = 1e9)
    fig, ax = plt.subplots()
    cm = ax.pcolormesh(
        smoothed_grid.r_range,
        smoothed_grid.z_range,
        smoothed_grid.grid',
        norm = LogNorm(vmin, vmax),
    )
    plt.colorbar(cm, ax = ax)
end

function plot_tau_x_grid(radiative_transfer; vmin = 1e2, vmax = 1e9)
    fig, ax = plt.subplots()
    r_range = range(0, 1000, length = 250)
    z_range = range(0, 1000, length = 250)
    tau_x_grid = 1e2 * ones((length(r_range), length(z_range)))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            tau_x_grid[i, j] = compute_xray_tau(radiative_transfer, r, z)
        end
    end
    cm = ax.pcolormesh(r_range, z_range, tau_x_grid', norm = LogNorm(vmin, vmax))
    plt.colorbar(cm, ax = ax)
end


config = YAML.load_file("scripts/config.yaml", dicttype = Dict{Symbol,Any})

black_hole = BlackHole(config)

radiation = getfield(Qwind, Symbol(config[:radiation][:mode]))(black_hole, config)


wind_grid = Rectangular(config)

initial_conditions = getfield(Qwind, Symbol(config[:initial_conditions][:mode]))(
    radiation,
    black_hole,
    config,
)


for kernel_size in [0, 1, 2, 5, 10, 50]
    config = YAML.load_file("scripts/config.yaml", dicttype = Dict{Symbol,Any})
    config[:radiative_transfer][:kernel_size] = kernel_size
    radiative_transfer =
        getfield(Qwind, Symbol(config[:radiative_transfer][:mode]))(radiation, config)
    iterations_dict = Dict()
    save_path = "results_kernel_size_$kernel_size" #config[:integrator][:save_path]
    mkpath(save_path)
    # iterations
    n_iterations = config[:integrator][:n_iterations]
    for it = 1:n_iterations
        @info "Starting iteration $it of $n_iterations"
        iterations_dict[it] = Dict()
        output_path = save_path * "/iteration_$(@sprintf "%03d" it).csv"
        integrators = initialize_integrators(
            radiative_transfer,
            wind_grid,
            initial_conditions,
            atol = config[:integrator][:atol],
            rtol = config[:integrator][:rtol],
            save_path = output_path,
        )
        iterations_dict[it]["integrators"] = integrators
        iterations_dict[it]["radiative_transfer"] = radiative_transfer
        run_integrators!(integrators)
        @info "Integration of iteration $it ended!"
        radiative_transfer = update_radiative_transfer(radiative_transfer, integrators)
    end
end




