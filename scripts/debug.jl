using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools
LogNorm = matplotlib.colors.LogNorm

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

fig, ax = plt.subplots()
dgrid = iterations_dict[2]["rad"].wi.density_grid
ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm())

QwindPlotting.plot_density_grid(iterations_dict[2]["rad"].wi.density_grid, rmax=200, zmax=50)

QwindPlotting.plot_streamlines(iterations_dict[2]["integrators"])

radiation = model.rad;

fig, ax = QwindPlotting.plot_xray_grid(radiation.wi.density_grid, radiation.xray_luminosity, radiation.bh.Rg, rmax=100, zmax=100, nr=500, nz=500, vmin=1e-2, vmax=1e2)


radiation = iterations_dict[2]["rad"];

model.rad = radiation;

integ1 = Qwind.create_and_run_integrator(model, r0=20.01, linewidth=10, trajectory_id=1, atol=1e-8, rtol=1e-3);
#integ2 = Qwind.create_and_run_integrator(model, r0=377, linewidth=10, trajectory_id=1, atol=1e-8, rtol=1e-3);

fig, ax = plt.subplots()
ax.plot(integ1.p.data[:r], integ1.p.data[:z])
#ax.plot(integ2.p.data[:z], integ2.p.data[:xi])
#ax.set_yscale("log")
#ax.set_xscale("log")
