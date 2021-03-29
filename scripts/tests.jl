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

LogNorm = matplotlib.colors.LogNorm

integrators = iterations_dict[1]["integrators"];
Profile.clear()
@profile rt = update_radiative_transfer(model.rt, integrators);
pprof()

QwindPlotting.plot_wind_hull(rt.interpolator.wind_hull, zmax=25)

dgrid = rt.interpolator.velocity_grid; #rt.interpolator.density_grid;
fig, ax = plt.subplots()
ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.vz_grid', norm=LogNorm())
