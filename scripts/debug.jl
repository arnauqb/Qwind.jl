using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")

model = Model("configs/to_debug.yaml");
mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("configs/to_debug.yaml");
iterations_dict = Dict();
run!(model, iterations_dict)


integrators = iterations_dict[3]["integrators"];

max_times = get_intersection_times(integrators);

hull = Hull(integrators, max_times)


#hull = Hull(integrators, max_times);
integrators_filtered = Qwind.filter_close_trajectories(integrators, 5e-2);
r0 = [integ.p.r0 for integ in integrators_filtered];
integrators_interpolated_linear = interpolate_integrators(
    integrators_filtered,
    max_times = max_times,
    n_timesteps = 50,
    log = true,
);
r, z, _, _, _ = reduce_integrators(integrators_interpolated_linear);


hull = Hull(r, z, r0, sigdigits=5)

fig, ax = plt.subplots()
for integ in integrators_interpolated_linear
    ax.plot(integ.r, integ.z)
end
QwindPlotting.plot_wind_hull(hull, zmax=2,ax=ax)
