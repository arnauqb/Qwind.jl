using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter, PyPlot
using Qwind
include("scripts/plotting.jl")

rin = 46.100482541805334
#rin = 20.0
config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
config[:initial_conditions][:r_in] = rin
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config);
#iterations_dict = Dict();
#iterations_dict[1] = Dict();
#run!(model, iterations_dict, parallel = true)

lr, lw = Qwind.compute_lines_range(model);
fig, ax = plt.subplots()
for ll in lr
    ax.axvline(ll)
end


using JLD2

@load "./scratch/trajectories.jld2" integrators
#@load "./scratch/max_times.jld2" max_times

max_times = get_intersection_times(integrators)


integrators_interp = interpolate_integrators(integrators, max_times=max_times, n_timesteps=100, log=false);

using ConcaveHull

fig, ax = plt.subplots()
for integ in integrators_interp
    ax.plot(integ.r, integ.z)
end

fig,ax = plt.subplots()
for integ in integrators 
    ax.plot(integ.p.data[:r], integ.p.data[:z])
end

points = Hull(integrators, max_times, hull_sigdigits=6);

fig,ax = plt.subplots()
ps = reduce(hcat, points)
ax.scatter(ps[1,:], ps[2,:])
hull = ConcaveHull.concave_hull(points)
#ax.scatter(hr, hz)
ax.set_xlim(0,1000)
ax.set_ylim(0,1000)

#hull = Hull(integrators, max_times)

fig, ax = plt.subplots()
QwindPlotting.plot_wind_hull(hull, fig=fig, ax=ax, rmin=40, rmax=5000, zmax=5000, nr=500, nz=500)
#ps = reduce(hcat, points)
#ax.scatter(ps[1,:], ps[2,:], s=10, alpha=0.5)
#vs = reduce(hcat, hull.vertices)
#ax.scatter(vs[1,:], vs[2,:], color="red", alpha=0.5)
for integ in integrators_interp
    ax.plot(integ.r, integ.z)
end
ax.set_xlim(0, 5000)
ax.set_ylim(0, 5000)



