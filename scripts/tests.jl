using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter, PyPlot, JLD2
using Qwind
LogNorm = matplotlib.colors.LogNorm
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config);
#iterations_dict = Dict();
#run!(model, iterations_dict)

#@load "iterations_dict.jld2" iterations_dict

Rg = model.bh.Rg;
xl = model.rad.xray_luminosity;

velocity_grid = iterations_dict[2]["radiative_transfer"].interpolator.velocity_grid;

dgrid = iterations_dict[2]["radiative_transfer"].interpolator.density_grid
fig, ax = plt.subplots()
cm = ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm())
plt.colorbar(cm, ax=ax)
ax.set_xlim(10, 12)

QwindPlotting.plot_xray_grid(dgrid, xl, Rg, rmax=1000, zmax=30, nr=250, nz=250)

integrators = iterations_dict[1]["integrators"];
for integ in integrators
    for (rp, zp, np) in zip(integ.p.data[:r], integ.p.data[:z], integ.p.data[:n])
        if rp > 200
            continue
        end
        density = get_density(dgrid, rp, zp)
        if abs(density - np) / np > 20
            println("r $rp z $zp")
            println("Should be \t $np")
            println("Is \t $density")
            println("---------------")
        end
    end
end

using ProgressMeter

rt = iterations_dict[2]["radiative_transfer"];
nr = 250
nz = 250
r_range = range(0, 1000, length=nr);
z_range = range(0, 10, length=nz);
force_grid_r = zeros((nr, nz));
force_grid_z = zeros((nr, nz));
xi_grid = 1e-20 .* ones((nr, nz));
tau_x_grid = zeros((nr, nz));
tau_x0_grid = zeros((nr, nz));
@showprogress for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        tau_x = compute_xray_tau(rt.interpolator.density_grid, 0.0, 0.0, r, z, xl, Rg)
        density = get_density(rt.interpolator.density_grid, r, z)
        xi = compute_ionization_parameter(r, z, density, tau_x, xl, Rg)
        xi_grid[i, j] = xi
        tau_x_grid[i, j] = tau_x
        tau_x0_grid[i, j] = log(xl / (1e5 * density * sqrt(r^2 + z^2) * Rg))
        #vr, vz = get_velocity(rt.interpolator.velocity_grid, r, z)
        #force = compute_disc_radiation_field(rt, r, z, vr, vz)
        #force_grid_r[i,j] = force[1]
        #force_grid_z[i,j] = force[2]
    end
end

fig, ax = plt.subplots(1, 1, figsize=(5,2))
cm = ax.pcolormesh(r_range, z_range, (tau_x_grid ./ tau_x0_grid)', norm=LogNorm(1e-2, 1e2))
plt.colorbar(cm, ax=ax)
#cm = ax[1].pcolormesh(r_range, z_range, tau_x_grid', norm=LogNorm())
#plt.colorbar(cm, ax=ax[1])
#cm = ax[2].pcolormesh(r_range, z_range, tau_x0_grid', norm=LogNorm())
#plt.colorbar(cm, ax=ax[2])
