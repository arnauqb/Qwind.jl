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
run_iteration!(model, iterations_dict, it_num=1, parallel=true);
#run!(model, iterations_dict, parallel=true)


fig, ax = QwindPlotting.plot_streamlines(iterations_dict[1]["integrators"], linestyle="-")


using PyPlot
LogNorm = matplotlib.colors.LogNorm;


interp = iterations_dict[2]["radiative_transfer"].interpolator;
vgrid = interp.velocity_grid
fig, ax = plt.subplots()
cm = ax.pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vr_grid')
plt.colorbar(cm, ax=ax)
