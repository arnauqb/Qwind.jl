using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

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

model, iterations_dict = get_model("./configs/config_test.yaml");
run!(model, iterations_dict)


sls = Streamlines("./example/results.hdf5");

fig, ax = plt.subplots()
for sl in sls
    ax.plot(sl.r, sl.z, linewidth=1)
end
ax.set_xlim(0,2500)
ax.set_ylim(0,2500)
ax.set_xlabel(L"$R$ [ $R_g$ ]")
ax.set_ylabel(L"$z$ [ $R_g$ ]")
fig.savefig("./docs/src/figures/streamlines.png", dpi=300)

dgrid = DensityGrid("./example/results.hdf5");

fig, ax = plt.subplots()
cm = ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm(1e4, 1e8))
plt.colorbar(cm, ax=ax)
ax.set_xlim(0,2500)
ax.set_ylim(0,2500)
ax.set_xlabel(L"$R$ [ $R_g$ ]")
ax.set_ylabel(L"$z$ [ $R_g$ ]")
fig.savefig("./docs/src/figures/density_grid.png", dpi=300)

vgrid = VelocityGrid("./example/results.hdf5");

fig, ax = plt.subplots(1, 3, figsize=(6,2), sharex=true, sharey=true)
ax[1].pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vr_grid', norm=Normalize(0, 0.5))
ax[2].pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vphi_grid', norm=Normalize(0, 0.5))
cm = ax[3].pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vz_grid', norm=Normalize(0, 0.5))
plt.colorbar(cm, ax=ax[3])
ax[1].set_ylim(0,2500)
ax[1].set_xlim(0,2500)
for i in 1:3
    ax[i].set_xlabel(L"$R$ [ $R_g$ ]")
end
ax[1].set_ylabel(L"$z$ [ $R_g$ ]")
ax[1].set_title("Vr")
ax[2].set_title("Vphi")
ax[3].set_title("Vz")
plt.subplots_adjust(wspace=0.05, hspace=0.05)
fig.savefig("./docs/src/figures/velocity_grid.png", dpi=300)






rr = range(10, 1500, length=100);
Ts = disk_temperature.(Ref(model.bh), rr) .^ 4 .* rr .^2
fuvs = uv_fractions(model.bh, rr)
plt.loglog(rr, Ts)
plt.loglog(rr, Ts .* fuvs)

plt.loglog(rr, fuvs)


streamlines = Streamlines("./1e8_0.5_rin75/results.hdf5", 1);

r,z = Qwind.reduce_streamlines(streamlines, rtol=1e-1, atol=1e-4);

hull = Hull(streamlines, atol=1e-4, rtol=1e-1);

QwindPlotting.plot_wind_hull(hull, zmax=10)

vertices = hull.vertices;
rv = [v[1] for v in vertices];
zv = [v[2] for v in vertices];

fig, ax =plt.subplots()
ax.scatter(r,z)

ax.scatter(rv, zv)




fig, ax = plt.subplots()
for sl in streamlines
    ax.plot(sl.r, sl.z, alpha=0.5, linewidth=1)
end
#QwindPlotting.plot_wind_hull(hull, ax=ax, zmax=20, rmin=20, rmax=25, nr=500, nz=500)
#ax.set_xlim(20,25)
#ax.set_ylim(0,20)


hull = Hull(streamlines)

QwindPlotting.plot_streamlines(iterations_dict[4]["integrators"])

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
