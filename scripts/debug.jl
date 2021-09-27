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

bh = BlackHole(1e8 * M_SUN, 0.5, 0.0);
r_range = 10 .^ range(log10(6), log10(1500), length=1000);
T_range = disk_temperature.(Ref(bh), r_range);

h5open("./sed.hdf5", "w") do file
    g = create_group(file, "temperature_profile")
    g["r_range"] = r_range
    g["T_range"] = T_range
    attributes(g)["Description"] = "Temperature profile for M=10^8, mdot = 0.5"
end

plt.loglog(r_range, T_range);

model, iterations_dict = get_model("./configs/debug.yaml");
run!(model, iterations_dict)

it_num = 20
itdict = iterations_dict[it_num];
integs=  itdict["integrators"];
radiation =  iterations_dict[it_num-1]["rad"];

fig, ax = plt.subplots()
for integ in integs
    if integ.p.r0 < 100
        continue
    end
    ax.loglog(integ.p.data[:z], integ.p.data[:taux])
end

integ = integs[1];
#plt.loglog(integ.p.data[:n])

d = integ.p.data
ns1 = Qwind.compute_density.(d[:r], d[:z], d[:vr], d[:vz], d[:r][1], d[:z][1], d[:vz][1], d[:n][1])

ns2 = Qwind.compute_density.(d[:r], d[:z], d[:vr], d[:vz], d[:r][1], d[:z][1], d[:vz][1], d[:n][1])

plt.loglog(ns1)
plt.loglog(ns2)

#ax.set_xlim(1e-1, 100)

QwindPlotting.plot_streamlines(integs);

QwindPlotting.plot_density_grid(radiation.wi.density_grid, rmax=500, zmax=50, nr=500, nz=500, vmin=1e6, vmax=1e13)

fig, ax = QwindPlotting.plot_xray_grid(radiation.wi.density_grid, radiation.xray_luminosity, radiation.bh.Rg, rmax=100, zmax=100, nr=500, nz=500, vmin=1e-2, vmax=1e2)
#QwindPlotting.plot_streamlines(integs, ax=ax, color = "black", alpha=0.5);
#ax.set_xlim(0,500)
#ax.set_ylim(0,500)



r_range = 10 .^ range(log10(6), log10(1500), length=50);
fluxes = disk_flux.(Ref(model.bh), r_range);
maxv, maxi = findmax(fluxes)

fig, ax = plt.subplots()
ax.loglog(r_range, fluxes)
ax.axvline(r_range[maxi])

n0s = [integ.p.n0 for integ in iterations_dict[20]["integrators"]];
r0s = [integ.p.r0 for integ in iterations_dict[20]["integrators"]];
plt.loglog(r0s, n0s, "o-")



integs = iterations_dict[20]["integrators"];
max_times = Qwind.get_intersection_times(integs);

integs_interp = Qwind.interpolate_integrators(integs, n_timesteps=1000, max_times=max_times, log=false);
fig, ax = plt.subplots()
for integ in integs_interp
    ax.plot(integ.r, integ.z)
end
ax.set_xlim(0,2000)
ax.set_ylim(0,2000)

QwindPlotting.plot_streamlines()

model2, iterations_dict2 = get_model("./configs/debug.yaml");
run!(model2, iterations_dict2)

for integ in reverse(iterations_dict[20]["integrators"])
    if escaped(integ)
        println("last escape at $(integ.p.data[:r][1])")
    end
end

for integ in reverse(iterations_dict2[10]["integrators"])
    if escaped(integ)
        println("last escape at $(integ.p.data[:r][1])")
    end
end


integs = iterations_dict2[10]["integrators"];
max_times = Qwind.get_intersection_times(integs);

integs_interp = Qwind.interpolate_integrators(integs, max_times=max_times, log=false);

fig, ax = plt.subplots()
for integ in integs_interp
    ax.plot(integ.r, integ.z)
end
for integ in integs
    ax.plot(integ.p.data[:r], integ.p.data[:z], color = "C1")
end



fig, ax = plt.subplots()
it_num = 9
integs1 = iterations_dict[it_num]["integrators"];
for integ in integs1
    #ax.plot(integ.p.data[:r], integ.p.data[:z], color = "C0")
end
integs2 = iterations_dict2[it_num]["integrators"];
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
rr = range(6, 1500, length=500)
ax.loglog(rr, getn0.(Ref(model), rr), "o-")
#ax.set_xlim(0,50)

lr, lw = Qwind.compute_lines_range(model);
println(length(lr))

fig, ax = plt.subplots()
for l in lr
    ax.axvline(l)
end
#ax.set_xlim(6, 10)




it_num = 2
rad1 = iterations_dict[it_num]["rad"];
rad2 = iterations_dict2[it_num]["rad"];

QwindPlotting.plot_density_grid(rad1.wi.density_grid, nr=100, nz=100, zmax=100, rmax=1000)

QwindPlotting.plot_density_grid(rad2.wi.density_grid, nr=100, nz=100, zmax=100, rmax=1000)

QwindPlotting.plot_xray_grid(rad1.wi.density_grid, rad1.xray_luminosity, rad1.bh.Rg, nr=100, nz=100, vmin=1e-2, vmax=1e2, rmax=100, zmax=10)

QwindPlotting.plot_xray_grid(rad2.wi.density_grid, rad2.xray_luminosity, rad2.bh.Rg, nr=100, nz=100, vmin=1e-2, vmax=1e2, rmax=100, zmax=10)

#QwindPlotting.plot_xray_grid(rad2.wi.density_grid, rad1.xray_luminosity, rad1.bh.Rg, ax=ax, nr=100, nz=100, vmin=1e-2, vmax=1e2)
#

QwindPlotting.plot_streamlines(iterations_dict[9]["integrators"]);

integs = iterations_dict[9]["integrators"];
max_times = Qwind.get_intersection_times(integs);

integs_interp = Qwind.interpolate_integrators(
    integs,
    max_times = max_times,
    n_timesteps = 100,
    log = true,
);
r, z, _, _, _ = reduce_integrators(integs_interp);
@info "Constructing wind hull"

plt.scatter(r, z)

hull = Hull(r, z)
