using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using Qwind
using YAML, Profile, PProf, PyCall
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model1 = Model(config_path);
iterations_dict1 = Dict();
#run_iteration!(model1, iterations_dict1, it_num=1, parallel=false);
run!(model1, iterations_dict1, parallel=true)

fig, ax = QwindPlotting.plot_streamlines(iterations_dict1[2]["integrators"], linestyle="-")


integrators_raw = iterations_dict1[2]["integrators"]
max_times = get_intersection_times(integrators_raw)

integrators = Qwind.filter_close_trajectories(integrators_raw, 5e-2);
r0 = [integ.p.r0 for integ in integrators];
integrators_interpolated_linear = interpolate_integrators(
    integrators,
    max_times = max_times,
    n_timesteps = 50,
    log = true,
);
r, z, _, _, _ = reduce_integrators(integrators_interpolated_linear);

hull = Hull(r, z, r0, sigdigits = 6)

QwindPlotting.plot_wind_hull(hull, nr=500, nz=500)


fig, ax = plt.subplots()
for integ in integrators_interpolated_linear
    ax.plot(integ.r, integ.z)
end


#QwindPlotting.plot_wind_hull(model1.rt.interpolator.wind_hull, zmax=1)

integs = iterations_dict[1]["integrators"];


xray_luminosity = model1.rad.xray_luminosity
Rg = model1.bh.Rg

QwindPlotting.plot_density_grid(iterations_dict1[3]["radiative_transfer"].interpolator.density_grid, zmax=1, vmin=1e6, vmax=1e10)


integrators = iterations_dict1[1]["integrators"];
n0s = [integ.p.n0 for integ in integrators];

integ = integrators[10];
fig, ax = plt.subplots()
ax.semilogy(integ.p.data[:z], integ.p.data[:n])


# XRAY
it_num = 2
di = iterations_dict1[it_num]["radiative_transfer"].interpolator;
#di2 = iterations_dict2[it_num]["radiative_transfer"].density_interpolator;
fig, ax = QwindPlotting.plot_xray_grid(di.density_grid, xray_luminosity, Rg, rmin=0, rmax=1000, zmax=50, nr=250, nz=250, vmin=1e-1, vmax=1e4)
#ax.set_yscale("log")
#ax.set_xscale("log")
ax.set_title("New")

fig, ax = QwindPlotting.plot_xray_grid(di2.grid, xray_luminosity, Rg, rmin=0, rmax=200, nr=250, nz=250, vmin=1e-1, vmax=1e4)
ax.set_title("Old")
ax.set_yscale("log")
#ax.set_xscale("log")




# DENSITY
it_num = 3
di = iterations_dict1[it_num]["radiative_transfer"].density_interpolator;
#di2 = iterations_dict2[it_num]["radiative_transfer"].density_interpolator;
fig, ax = QwindPlotting.plot_density_grid(di.grid)
#ax.set_yscale("log")
#ax.set_xscale("log")
ax.set_title("New")

fig, ax = QwindPlotting.plot_density_grid(di2.grid, xlim=(0,1e3), ylim=(1e-6,1e3))
ax.set_title("Old")
ax.set_yscale("log")
#ax.set_xscale("log")


QwindPlotting.plot_streamlines(iterations_dict2[4]["integrators"])


integ2 = iterations_dict2[3]["integrators"];
old_rt = iterations_dict2[4]["radiative_transfer"];
old_di = old_rt.density_interpolator;
new_rt = update_radiative_transfer(model1.rt, integ2);
new_di = new_rt.density_interpolator;

fig, ax = QwindPlotting.plot_xray_grid(model1.rt.interpolator.density_grid, model1.rad.xray_luminosity, model1.bh.Rg, rmin=45, nr=500, nz=500, vmin=1e-2, vmax=1e2)
ax.set_yscale("log")
ax.set_xscale("log")

fig, ax = QwindPlotting.plot_xray_grid(old_di.grid, xray_luminosity, Rg, rmin=45, nr=500, nz=500, vmin=1e-2, vmax=1e2)
ax.set_title("Old")
ax.set_yscale("log")
ax.set_xscale("log")

it_num = 2
integ1 = iterations_dict1[it_num]["integrators"];
integ2 = iterations_dict2[it_num]["integrators"];
fig, ax = QwindPlotting.plot_streamlines(integ1)
ax.set_title("new")
fig, ax = QwindPlotting.plot_streamlines(integ2)
ax.set_title("old")
