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
run!(model, iterations_dict, parallel = true)

QwindPlotting.plot_wind_hull(
    iterations_dict[2]["radiative_transfer"].interpolator.wind_hull,
    nr = 500,
    nz = 500,
    zmax=20
)

integrators = iterations_dict[1]["integrators"];
intersection_times = get_intersection_times(integrators);
integs_interpol = interpolate_integrators(
    integrators,
    max_times = intersection_times,
    n_timesteps = 100,
    log = true,
);
r, z, _, _, _ = reduce_integrators(integs_interpol);
r0 = [integ.p.r0 for integ in integrators];
windh = Hull(r, z, r0, sigdigits = 6);


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
