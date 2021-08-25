using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
LogNorm = matplotlib.colors.LogNorm
SymLogNorm = matplotlib.colors.SymLogNorm
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

model2, iterations_dict2 = get_model("./configs/debug.yaml");
run!(model2, iterations_dict2)


fig, ax = plt.subplots()
it_num = 5
integs1 = iterations_dict[it_num]["integrators"];
integs2 = iterations_dict2[it_num]["integrators"];
for integ in integs1
    ax.plot(integ.p.data[:r], integ.p.data[:z], color = "C0")
end
for integ in integs2
    ax.plot(integ.p.data[:r], integ.p.data[:z], color = "C1")
end




rr = 10 .^ range(log10(6.1), log10(50), length=100);
zz = 10 .^ range(log10(1e-2), log10(10), length=100);
r_force = zeros(length(rr), length(zz));
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        r_force[i, j] = compute_disc_radiation_field(model.rad, r=r, z=z, vr=0, vz=0)[1]
    end
end

fig, ax = plt.subplots(figsize=(1,1))
cm = ax.contourf(rr, zz, r_force', norm=SymLogNorm(linthresh=1e-5),  cmap="RdBu_r")
cbar = plt.colorbar(cm, ax=ax)
cbar.set_label("Radial radiation force", rotation=270, labelpad=20, fontsize=15)
ax.set_ylabel("z [ Rg ]", fontsize=15)
ax.set_xlabel("R [ Rg ]", fontsize=15)


integrators = iterations_dict[2]["integrators"];
QwindPlotting.plot_streamlines(integrators)




fig,ax = plt.subplots()
for integ in integrators[76:end]
    ax.loglog(integ.p.data[:z], integ.p.data[:fm])
end



fig, ax = plt.subplots()
rr = range(6, 1500, length=5000)
ax.loglog(rr, getn0.(Ref(model), rr), "o-")
#ax.set_xlim(0,50)

lr, lw = Qwind.compute_lines_range(model);
println(length(lr))

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

bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
r_range = 10 .^ range(log10(6), 3, length=500);
nts = disk_nt_rel_factors.(Ref(bh), r_range);
uvf = uv_fractions(bh, r_range) .* disk_flux.(Ref(bh), r_range) .* r_range;
#uvf = disk_flux.(Ref(bh), r_range);
plt.loglog(r_range, uvf)
