using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter, PyPlot, JLD2
using Qwind
LogNorm = matplotlib.colors.LogNorm
include("scripts/plotting.jl")

config_path = "configs/rin_scan.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config);
#iterations_dict = Dict();
#run!(model, iterations_dict)
#lr, lw = Qwind.compute_lines_range(model)
#fig, ax = plt.subplots()
#for l in lr
#    ax.axvline(l)
#end



integrators = iterations_dict[1]["integrators"];
max_times = get_intersection_times(integrators);

fig, ax =plt.subplots()
for integrator in integrators 
    ax.plot(integrator.p.data[:r], integrator.p.data[:z])
end


dense_integrators = DenseIntegrator.(integrators);

max_indices = Dict()
for integrator in dense_integrators
    max_indices[integrator.id] = searchsorted_nearest(integrator.t, max_times[integrator.id])
end


sliced_integs = []
for integ in dense_integrators
    push!(sliced_integs, Qwind.slice_integrator(integ, fi = max_indices[integ.id]))
end

fig, ax =plt.subplots()
for integ in sliced_integs
    ax.plot(integ.r, integ.z)
end

integrators_interpolated_linear = interpolate_integrators(
    integrators,
    max_times = max_times,
    n_timesteps = 100,
    log = true,
)


fig, ax = plt.subplots()
#for integ in integrators
#    ax.plot(integ.p.data[:r], integ.p.data[:z])
#end
for integ in integrators_interpolated_linear
    ax.plot(integ.r, integ.z)
end

max_times = get_intersection_times(integrators);

integs2 = Qwind.filter_close_trajectories(integrators, 0.5);
length(integs2)

points = Hull(integrators, max_times);

QwindPlotting.plot_wind_hull(hull, zmax=25);

fig, ax = plt.subplots()
ps = reduce(hcat, points)
ax.scatter(ps[1,:], ps[2,:])
