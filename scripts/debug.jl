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

integrators = iterations_dict[10]["integrators"];

QwindPlotting.plot_streamlines(integrators[end:end])




fig,ax = plt.subplots()
for integ in integrators[76:end]
    ax.loglog(integ.p.data[:z], integ.p.data[:fm])
end



fig, ax = plt.subplots()
rr = range(6, 1500, length=5000)
ax.loglog(rr, getn0.(Ref(model), rr), "o-")
#ax.set_xlim(0,50)

lr, lw = Qwind.compute_lines_range(model);
length(lr)

fig, ax = plt.subplots()
for l in lr
    ax.axvline(l)
end
#ax.set_xlim(6, 10)

rad = iterations_dict[2]["rad"];
QwindPlotting.plot_xray_grid(rad.wi.density_grid, rad.xray_luminosity, rad.bh.Rg, vmin=1e-2, vmax=1e2)

integrators = iterations_dict[1]["integrators"];
fig, ax = plt.subplots()
QwindPlotting.plot_streamlines(integrators, ax=ax, alpha=0.25, color="black")
#ax.set_xlim(0,5000)
#ax.set_ylim(0,5000)
#

rad = iterations_dict[6]["rad"];
QwindPlotting.plot_density_grid(rad.wi.density_grid);

r0s = [integ.p.r0 for integ in integrators];
n0s = [integ.p.n0 for integ in integrators];
plt.loglog(r0s, n0s)

