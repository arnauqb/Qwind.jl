using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter
using Qwind
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();
run!(model, iterations_dict, parallel = true)


# Test Relativistic vs Non-Relativistic
rad = model.rad;
rt = model.rt;
norel_radiation = QsosedRadiation(
    rad.disk_grid,
    rad.fuv_grid,
    rad.mdot_grid,
    rad.xray_luminosity,
    rad.efficiency,
    rad.spin,
    rad.isco,
    rad.z_xray,
    rad.Rg,
    NoRelativistic(),
);
norel_rt = RegularGrid(norel_radiation, rt.interpolator);

density_grid = model.rt.interpolator.density_grid;
velocity_grid = model.rt.interpolator.velocity_grid;

r_range = range(0, 1000, length = 50);
z_range = range(0, 1000, length = 51);
rel_grid_r = zeros((length(r_range), length(z_range)));
rel_grid_z = zeros((length(r_range), length(z_range)));
norel_grid_r = zeros((length(r_range), length(z_range)));
norel_grid_z = zeros((length(r_range), length(z_range)));
vel_grid_r = zeros((length(r_range), length(z_range)));
vel_grid_z = zeros((length(r_range), length(z_range)));
@showprogress for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        vr, vz = interpolate_velocity(velocity_grid, r, z)
        vel_grid_r[i,j] = vr
        vel_grid_z[i,j] = vz
        force_rel = compute_disc_radiation_field(rt, r, z, vr, vz)
        force_norel = compute_disc_radiation_field(norel_rt, r, z, vr, vz)
        rel_grid_r[i, j] = force_rel[1]
        rel_grid_z[i, j] = force_rel[2]
        norel_grid_r[i, j] = force_norel[1]
        norel_grid_z[i, j] = force_norel[2]
    end
end

using PyPlot
LogNorm = matplotlib.colors.LogNorm

vmin = 1e-8
vmax = 1e-4
fig, ax = plt.subplots(1, 2)
cm1 = ax[1].pcolormesh(
    r_range,
    z_range,
    rel_grid_r',
    norm = LogNorm(vmin = vmin, vmax = vmax),
    shading = "auto",
)
cm2 = ax[2].pcolormesh(
    r_range,
    z_range,
    norel_grid_r',
    norm = LogNorm(vmin = vmin, vmax = vmax),
    shading = "auto",
)
plt.colorbar(cm2, ax = ax[2])
for axis in ax
    QwindPlotting.plot_streamlines(iterations_dict[2]["integrators"], fig=fig, ax=axis, color="white", linewidth=0.1, markersize=0.1)
    axis.set_xlim(r_range[1], r_range[end])
    axis.set_ylim(z_range[1], z_range[end])
end

fig, ax = plt.subplots(1, 2)
cm = ax[1].pcolormesh(r_range, z_range, vel_grid_r', norm=LogNorm(vmin=1e-3, vmax=1e-1))
cm = ax[2].pcolormesh(r_range, z_range, vel_grid_z', norm=LogNorm(vmin=1e-3, vmax=1e-1))
for axis in ax
    QwindPlotting.plot_streamlines(iterations_dict[2]["integrators"], fig=fig, ax=axis, color="black", linewidth=0.1, markersize=0.1)
    axis.set_xlim(r_range[1], r_range[end])
    axis.set_ylim(z_range[1], z_range[end])
end

