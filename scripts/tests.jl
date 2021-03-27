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
iterations_dict[1] = Dict();
run_integrators!(model, iterations_dict, it_num=1, parallel=true);
#run!(model, iterations_dict, parallel = false)

QwindPlotting.plot_wind_hull(
    iterations_dict[2]["radiative_transfer"].interpolator.wind_hull,
    nr = 500,
    nz = 500,
    rmax=20,
    zmax=20
)

integrators = iterations_dict[1]["integrators"];
dense_integrators = DenseIntegrator.(integrators);

fig, ax = plt.subplots()
for integ in dense_integrators
    ax.plot(integ.r, integ.z)
end

intersections_dict = Qwind.get_intersections(dense_integrators);
intersection_times = Dict(integrator.id => integrator.t[end] for integrator in dense_integrators);

#Profile.clear()
Qwind.resolve_intersections!(intersection_times, intersections_dict);
#pprof()

intersection_indices = [searchsortedfirst(integrator.t, intersection_times[integrator.id]) for integrator in dense_integrators];
sliced_integs = [Qwind.slice_integrator(integrator, fi=index) for (integrator, index) in zip(dense_integrators, intersection_indices)];
fig, ax = plt.subplots()
for integ in sliced_integs
    ax.plot(integ.r, integ.z)
end


#intersection_times = get_intersection_times(integrators);

integs_interpol = interpolate_integrators(
    integrators,
    max_times = intersection_times,
    n_timesteps = 200,
    log = true,
);
r, z, _, _, _ = reduce_integrators(integs_interpol);
r0 = [integ.p.r0 for integ in integrators];
fig, ax = plt.subplots()
for integ in integs_interpol
    ax.plot(integ.r, integ.z)
end
ax.set_xlim(0,100)

windh = Hull(r, z, r0, sigdigits = 2);


fig, ax = plt.subplots()
#QwindPlotting.plot_streamlines(integrators, fig=fig, ax=ax)
#for integ in integs_interpol
#    ax.plot(integ.r, integ.z, "o-", color = "black")
#end
#ps = reduce(hcat, points);
#ax.scatter(ps[1, :], ps[2, :])
QwindPlotting.plot_wind_hull(
    windh,
    zmax = 50,
    nr = 500,
    nz = 500,
    fig = fig,
    ax = ax,
)
ax.set_xlim(0, 500)
ax.set_ylim(0, 50)


fig, ax = plt.subplots()
r_range = 10 .^ range(log10(6), log10(1000), length=50);
ic = CAKIC(model.rad, model.bh, 50.0, 1500.0, 100, 0.0, 0.03, 0.6, 0.5, true) 
ax.loglog(r_range, getn0.(Ref(ic), r_range), label = "current")
icp = CAKIC(model.rad, model.bh, 50.0, 1500.0, 100, 0.0, 0.002, 0.6, 0.5, true) 
ax.loglog(r_range, getn0.(Ref(icp), r_range), label = "Pereyra")
#ic = CAKIC(model.rad, model.bh, 50.0, 1500.0, 100, 0.0, 0.007, 0.75, 0.5, true) 
#ax.loglog(r_range, getn0.(Ref(ic), r_range), label = "K = 0.03")
#ic2 = CAKIC(model.rad, model.bh, 50.0, 1500.0, 100, 0.0, 0.0026, 0.737, 0.5, true) 
#ax.loglog(r_range, getn0.(Ref(ic2), r_range), label = "K = 0.2")
#ic3 = CAKIC(model.rad, model.bh, 50.0, 1500.0, 100, 0.0, 0.0021, 0.8, 0.5, true) 
#ax.loglog(r_range, getn0.(Ref(ic3), r_range), label = "K = 0.002")
icauto = CAKIC(model.rad, model.bh, 50.0, 1500.0, 100, 0.0, "auto", 0.6, 0.5, true) 
ax.loglog(r_range, getn0.(Ref(icauto), r_range), label = "K = auto")
ax.legend()




integ1 = integrators[11];
integ2 = integrators[12];
integ3 = integrators[13];
integ11 = integs_interpol[11];
integ22 = integs_interpol[12];
integ33 = integs_interpol[13];
fig, ax = plt.subplots()
ax.plot(integ1.p.data[:r], integ1.p.data[:z], linewidth=5)
ax.plot(integ2.p.data[:r], integ2.p.data[:z], linewidth=5)
ax.plot(integ3.p.data[:r], integ3.p.data[:z], linewidth=5)
ax.plot(integ11.r, integ11.z, color="black", linewidth=2)
ax.plot(integ22.r, integ22.z, color="black", linewidth=2)
ax.plot(integ33.r, integ33.z, color="black", linewidth=2)

fig, ax = plt.subplots()
vt(integ) = sqrt.(integ.p.data[:vr] .^2 + integ.p.data[:vz] .^ 2)
ax.semilogy(integ2.p.data[:r], vt(integ2))
ax.semilogy(integ3.p.data[:r], vt(integ3))

momentum1 = integ2.p.data[:n]# .* vt(integ2) .^ 2
momentum2 = integ3.p.data[:n]# .* vt(integ3) .^ 2
fig, ax = plt.subplots()
ax.semilogy(momentum1)
ax.semilogy(momentum2)

