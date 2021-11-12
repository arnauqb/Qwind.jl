using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

function get_model(config)
    model = Model(config)
    try
        mv(model.config[:integrator][:save_path], "backup", force = true)
    catch
    end
    model = Model(config)
    iterations_dict = Dict()
    return model, iterations_dict
end
#model, iterations_dict = get_model("./configs/debug.yaml");
model, iterations_dict = get_model("./todebug.yaml");
run!(model, iterations_dict);

iterations_dict[1] = Dict()
integrators, streamlines =
    Qwind.run_integrators!(model, iterations_dict, it_num = 1, parallel = true);
wi = model.rad.wi;

#new_radiation = update_radiation(model.rad, streamlines)

hull = Hull(streamlines, rtol=1e-2);

density_grid = update_density_grid(wi.density_grid, wi.update_grid_flag, streamlines, hull);
velocity_grid =
    update_velocity_grid(wi.velocity_grid, wi.update_grid_flag, streamlines, hull);


#Profile.clear()
@time ionization_grid = IonizationGrid(
    density_grid,
    Rg = model.rad.bh.Rg,
    xray_luminosity = model.rad.xray_luminosity,
    z_xray = model.rad.z_xray,
    parallel = false,
);
#pprof()

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    ionization_grid.r_range,
    ionization_grid.z_range,
    ionization_grid.grid',
    norm = LogNorm(),
)
plt.colorbar(cm, ax = ax)
#
#fig, ax = plt.subplots()
#cm = ax.pcolormesh(
#    density_grid.r_range,
#    density_grid.z_range,
#    density_grid.grid',
#    norm = LogNorm(),
#)
#plt.colorbar(cm, ax = ax)

ion0 = Qwind.compute_initial_ionization_grid(
    density_grid,
    model.rad.xray_luminosity,
    model.bh.Rg,
    model.rad.z_xray,
);
absorbed_from_center = compute_luminosity_absorbed_grid(
    density_grid,
    ion0,
    Rg = model.bh.Rg,
    source_position = [0, 0, model.rad.z_xray],
    source_luminosity = model.rad.xray_luminosity,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
);


nodes,weights = gausslegendre(10);
i = 25;
j = 26;
#Profile.clear()
@time Qwind.compute_scattered_flux_in_cell(
    density_grid,
    ion0;
    scattered_luminosity_per_cell = absorbed_from_center,
    cell = Rectangle(
        density_grid.r_range[i],
        density_grid.r_range[i + 1],
        density_grid.z_range[j],
        density_grid.z_range[j + 1],
    ),
    cell_density = 1e8,
    absorption_opacity = BoostOpacity(),
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    Rg = model.bh.Rg,
    nodes=nodes,
    weights=weights,
)
#pprof()

times = zeros(length(density_grid.r_range)-1, length(density_grid.z_range)-1);
for i = 1:(length(density_grid.r_range) - 1)
    for j = 1:(length(density_grid.z_range) - 1)
        r_source = (density_grid.r_range[i + 1] + density_grid.r_range[i]) / 2
        z_source = (density_grid.z_range[j + 1] + density_grid.z_range[j]) / 2
        time = @elapsed compute_optical_depth(
            density_grid.iterator,
            density_grid,
            ion0,
            BoostOpacity(),
            ri = 0,
            phii = π/2,
            zi = 0,
            rf = density_grid.r_range[i],
            phif = 0,
            zf = density_grid.z_range[j],
            Rg = model.bh.Rg,
            mu_electron = 1.17,
            mu_nucleon = 0.61,
            max_tau = 30,
        )
        times[i, j] = time
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(density_grid.r_range[1:end-1], density_grid.z_range[1:end-1], times', norm=LogNorm())
plt.colorbar(cm, ax = ax)
