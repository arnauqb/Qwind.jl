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

streamlines = Streamlines("./debugging/results.hdf5", 10);
fig, ax = plt.subplots()
for sl in streamlines
    ax.loglog(sl.z, sl.vz)
end

QwindPlotting.plot_streamlines(iterations_dict[5]["integrators"])

fig, ax = plt.subplots()
dgrid = iterations_dict[2]["rad"].wi.density_grid
ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm())

QwindPlotting.plot_density_grid(iterations_dict[10]["rad"].wi.density_grid, rmax=200, zmax=50, nr=200, nz=200)



fig, ax = plt.subplots()
for line in iterations_dict[5]["integrators"]
    if escaped(line)
        ax.plot(line.p.data[:r], line.p.data[:z])
    end
end

radiation = model.rad;
fig, ax = QwindPlotting.plot_xray_grid(radiation.wi.density_grid, radiation.xray_luminosity, radiation.bh.Rg, rmax=100, zmax=1, nr=250, nz=250, vmin=1e-2, vmax=1e2)


radiation = iterations_dict[5]["rad"];

model.rad = radiation;

integ1 = Qwind.create_and_run_integrator(model, r0=500, linewidth=10, trajectory_id=1, atol=1e-8, rtol=1e-3);
integ2 = Qwind.create_and_run_integrator(model, r0=700, linewidth=10, trajectory_id=1, atol=1e-8, rtol=1e-3);
fig, ax = plt.subplots()
ax.plot(integ1.p.data[:z], integ1.p.data[:])
ax.plot(integ2.p.data[:z], integ2.p.data[:])
ax.set_yscale("log")
ax.set_xscale("log")

tauuvs1 = []
tauuvs2 = []
for (r, z) in zip(integ1.p.data[:r], integ1.p.data[:z])
    tauuv = compute_tau_uv(radiation, rd=0.0, phid=0.0, r=r, z=z)
    push!(tauuvs1, tauuv)
end
for (r, z) in zip(integ2.p.data[:r], integ2.p.data[:z])
    tauuv = compute_tau_uv(radiation, rd=0.0, phid=0.0, r=r, z=z)
    push!(tauuvs2, tauuv)
end

plt.plot(tauuvs1)
plt.plot(tauuvs2)
