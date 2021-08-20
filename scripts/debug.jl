using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools

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

model, iterations_dict = get_model("./configs/debug.yaml");
run!(model, iterations_dict)

lr, lw = Qwind.compute_lines_range(model, rin=6.1, rfi=1500, delta_mdot=0.01)

fig, ax = plt.subplots()
for l in lr
    ax.axvline(lr)
end


iterations_dict[1] = Dict()
integrators = run_integrators!(model, iterations_dict, it_num=1, parallel=true);



integrators = iterations_dict[1]["integrators"];
fig, ax = plt.subplots()
QwindPlotting.plot_streamlines(integrators, ax=ax, alpha=0.25, color="black")
ax.set_xlim(0,2000)
ax.set_ylim(0,2000)

max_times = get_intersection_times(integrators);

r0 = [integ.p.r0 for integ in integrators];
integrators_interpolated_linear = interpolate_integrators(
    integrators,
    max_times = max_times,
    n_timesteps = 100,
    log = true,
);
#QwindPlotting.plot_streamlines(integrators_interpolated_linear)
r, z, _, _, _ = reduce_integrators(integrators_interpolated_linear);
plt.scatter(r, z)


hull = Hull(r, z, r0, sigdigits = 6);

hull = iterations_dict[1]["rad"].wi.wind_hull;

integrators = iterations_dict[1]["integrators"];
whull = iterations_dict[2]["rad"].wi.wind_hull;

fig, ax = QwindPlotting.plot_wind_hull(whull, rmin=1, rmax=2000, zmax=100, nr =500, nz=500);
QwindPlotting.plot_streamlines(integrators, ax=ax, alpha=0.25)
