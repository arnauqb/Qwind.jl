using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools, ProgressMeter
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

function get_model(config)
    model = Model(config)
    try
        mv(model.parameters.save_path, "backup", force = true)
    catch
    end
    model = Model(config)
    iterations_dict = Dict()
    return model, iterations_dict
end
model, iterations_dict = get_model("./configs/debug.yaml");
run!(model, iterations_dict, parallel=true)
#integs, sls = Qwind.run_integrators(model, iterations_dict; it_num=1, parallel = true);

#run_iteration!(model, iterations_dict, parallel = true, it_num=1);


Profile.clear()
@profile integ = Qwind.create_and_run_integrator(
        model,
        linewidth = 1,
        r0 = 100,
        trajectory_id = 1,
    );
pprof()

@code_warntype Qwind.compute_radiation_acceleration(model.rad, [0.1, 0.1, 0.1, 0.1], [1000, 1, 0.1, 0.1], integ.p)

@code_warntype Qwind.compute_force_multiplier(1e-3, 1e2, FMInterp())

@code_warntype Qwind.compute_force_multiplier_k(1e2)


streamlines = iterations_dict[1]["streamlines"];
fig, ax = plt.subplots()
for sl in streamlines
    ax.plot(sl.r, sl.z)
end


r, z = Qwind.reduce_streamlines(streamlines, rtol=1e-2, atol=1e-2)
plt.scatter(r, z)

hull = Hull(r, z)

fig, ax = plt.subplots()
QwindPlotting.plot_wind_hull(hull, nr=500, nz=500, zmax=50, ax=ax, rmax=1500)
ax.scatter(r, z, s=1)

QwindPlotting.plot_density_grid(model.rad.wi.density_grid, zmax=50)

QwindPlotting.plot_wind_hull(model.rad.wi.wind_hull, zmax=50, nr=250, nz=250)

QwindPlotting.plot_streamlines(iterations_dict[1]["integrators"])

@everywhere function compute_xi(rad, r, z; include_scattering = true)
    n = get_density(rad.wi.density_grid, r, z)
    tau_x = compute_tau_xray(rad, ri = 0, phii = 0, zi = 0, rf = r, zf = z, phif = 0)
    ret = compute_ionization_parameter(
        r = r,
        z = z,
        vr = 0,
        vz = 0,
        number_density = n,
        tau_x = tau_x,
        xray_luminosity = rad.xray_luminosity,
        Rg = rad.bh.Rg,
        include_scattering = include_scattering,
        density_grid = rad.wi.density_grid,
        absorption_opacity = rad.xray_opacity,
        zh = rad.z_xray,
        mu_electron = rad.mu_electron,
        mu_nucleon = rad.mu_nucleon,
        scattered_luminosity_grid = rad.wi.scattered_lumin_grid,
    )
    return ret
end

@everywhere rr = collect(10 .^ range(log10(6), log10(1000), length = 50))
@everywhere zz = 10 .^ range(log10(1e-3), log10(1000), length = 51)
@everywhere rets = zeros(length(rr), length(zz))
@showprogress for (i, r) in enumerate(rr)
    rets[i, :] = pmap(
        z -> compute_xi(model.rad, r, z, include_scattering = true),
        zz,
        batch_size = 10,
    )
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, rets', norm = LogNorm(), shading="auto")
cbar = plt.colorbar(cm, ax = ax)
ax.set_xlabel("R [Rg]")
ax.set_ylabel("R [Rg]")
cbar.set_label("Ionization parameter", rotation=270)




r = 1000
z = 1
phi = 0
n = 1e9
tau_x = compute_tau_xray(model.rad, ri = 0.0, phii = 0, zi = 0, rf = r, phif = phi, zf = z)
compute_ionization_parameter(model.rad, 
    r = r,
    z = z,
    vr = 0,
    vz = 0,
    n = n,
    tau_x = tau_x,
)


rr = range(6, 1500, length = 500);
zz = 10 .^ range(-6, log10(1500), length = 500);
ret = zeros(length(rr), length(zz));
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        ret[i, j] = compute_tau_xray(
            model.rad,
            ri = 0.0,
            phii = 0.0,
            zi = model.rad.z_xray,
            rf = r,
            phif = 0.0,
            zf = z,
        )
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, ret', norm = LogNorm(vmin = 1e-2, vmax = 1e2))
plt.colorbar(cm, ax = ax)



scattered_grid = ScatteredLuminosityGrid()

rr = range(6, 1500, length = 500);
zz = 10 .^ range(-6, log10(1500), length = 500);
real_grid = zeros(length(rr), length(zz));
int_grid = zeros(length(rr), length(zz));
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        taux = compute_tau_xray(
            model.rad,
            ri = 0.0,
            phii = 0.0,
            zi = model.rad.z_xray,
            rf = r,
            phif = 0.0,
            zf = z,
        )
        density = interpolate_density(model.rad.wi.density_grid, r, z)
        xi = compute_ionization_parameter(model.rad, r, z, 1e-4, 1e-4, density, taux)
        int_grid[i, j] =
            interpolate_ionization_parameter(model.rad.wi.ionization_grid, r, z)
        real_grid[i, j] = xi
    end
end


fig, ax = plt.subplots(1, 2)
ax[1].pcolormesh(rr, zz, real_grid', norm = LogNorm())
ax[2].pcolormesh(rr, zz, int_grid', norm = LogNorm())

fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, (real_grid ./ int_grid)', norm = LogNorm())
plt.colorbar(cm, ax = ax)

iterations_dict[1] = Dict()
integrators, streamlines =
    Qwind.run_integrators!(model, iterations_dict, it_num = 1, parallel = true);
wi = model.rad.wi;

#new_radiation = update_radiation(model.rad, streamlines)

hull = Hull(streamlines, rtol = 1e-2);

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


nodes, weights = gausslegendre(10);
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
    nodes = nodes,
    weights = weights,
)
#pprof()

times = zeros(length(density_grid.r_range) - 1, length(density_grid.z_range) - 1);
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
            phii = π / 2,
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
cm = ax.pcolormesh(
    density_grid.r_range[1:(end - 1)],
    density_grid.z_range[1:(end - 1)],
    times',
    norm = LogNorm(),
)
plt.colorbar(cm, ax = ax)
