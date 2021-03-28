using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter, PyPlot, JLD2
using Qwind
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config);
iterations_dict = Dict();
run!(model, iterations_dict)

integrators = iterations_dict[1]["integrators"];
r0 = [integ.p.r0 for integ in integrators];
max_times = get_intersection_times(integrators);

integrators2 = Qwind.filter_close_trajectories(integrators, 5e-2);
max_times2 = [max_times[integrator.p.id] for integrator in integrators];
length(integrators2)

fig, ax = plt.subplots()
for integ in integrators2
    ax.plot(integ.p.data[:r], integ.p.data[:z])
end

integrators_interpolated_linear = interpolate_integrators(
    integrators2,
    max_times = max_times2,
    n_timesteps = 100,
    log = true,
);
r, z, _, _, _ = reduce_integrators(integrators_interpolated_linear);

fig, ax = plt.subplots()
for integ in integrators_interpolated_linear
    ax.scatter(integ.r, integ.z)
end

points = Hull(r, z, r0, sigdigits=6);

hull = ConcaveHull.concave_hull(points);

fig, ax = plt.subplots()
ps = reduce(hcat, points);
ax.scatter(ps[1,:], ps[2,:])
QwindPlotting.plot_wind_hull(hull, zmax=20, rmax=1500);

