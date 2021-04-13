using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter, JLD2, PyPlot
using Qwind
LogNorm = matplotlib.colors.LogNorm
include("scripts/plotting.jl")

config_path = "configs/debug.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config);
iterations_dict = Dict();
run!(model, iterations_dict)

integrator = initialize_integrator(model.rt, model.wind_grid, model.ic, 25.0, 1);
solve!(integrator)

fig, ax = plt.subplots()
#ax.loglog(integrator.p.data[:z], integrator.p.data[:n])
ax.plot(integrator.p.data[:r], integrator.p.data[:z])



#update_radiative_transfer(model.rt, iterations_dict[1]["integrators"])
integrators = iterations_dict[6]["integrators"];
fig, ax = plt.subplots()
for integ in integrators
    #ax.loglog(integ.p.data[:z], integ.p.data[:n])
    ax.plot(integ.p.data[:r], integ.p.data[:z])
end


r0 = [integ.p.r0 for integ in integrators];
max_times = get_intersection_times(integrators);
hull = Hull(integrators, max_times);

density_grid_old = DensityGrid(integrators, max_times, hull, nr=1000, nz=500);

fig, ax = plt.subplots()
ax.pcolormesh(density_grid.r_range, density_grid.z_range, density_grid.grid', norm=LogNorm())

r_grid = density_grid.r_range' .* ones(length(density_grid.z_range));
z_grid = density_grid.z_range .* ones(length(density_grid.r_range))';

#@load "iterations_dict.jld2" iterations_dict

integrators = load_trajectories("./results.hdf5", 2);

fig, ax = plt.subplots()
for integ in integrators[1400:1500]
    ax.plot(integ.r, integ.z)
end


