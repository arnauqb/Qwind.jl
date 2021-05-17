using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using Qwind
using YAML, Profile, PProf, PyCall
include("scripts/plotting.jl")

config_path = "configs/debug.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();
#run_iteration!(model1, iterations_dict1, it_num=1, parallel=false);
run!(model, iterations_dict, parallel=true)

fig, ax = QwindPlotting.plot_streamlines(iterations_dict[2]["integrators"], linestyle="-")


using PyPlot
LogNorm = matplotlib.colors.LogNorm;

interp = iterations_dict[2]["radiative_transfer"].interpolator;
dgrid = interp.density_grid
fig, ax = plt.subplots()
cm = ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm(1e4, 1e10))
plt.colorbar(cm, ax=ax)
