using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind

using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools, ProgressMeter
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

function get_model(config)
    model = Model(config)
    try
        mv(model.parameters.save_path, "backup", force = true)
    catch
    end
    model = Model(config)
    iterations_dict = Dict()
    return model, iterations_dict
end
model, iterations_dict = get_model("./configs/proga.yaml");
run!(model, iterations_dict, parallel = true)

integ = Qwind.create_and_run_integrator(model; r0 = 1000, linewidth = 1, trajectory_id = 1);

fig, ax = plt.subplots()
d = sqrt.(integ.p.data[:r] .^ 2 + integ.p.data[:z] .^ 2)
ax.plot(integ.p.data[:z], integ.p.data[:taux])


taux = compute_tau_xray(
    model.wind.density_grid,
    model.wind.grid_iterator,
    model.rad,
    model.parameters,
    ri = 0.0,
    zi = 1.5,
    phii = 0,
    rf = 1000,
    phif = 0,
    zf = 300,
)

fig, ax = plt.subplots()
dgrid = model.wind.density_grid;
cm = ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm())
plt.colorbar(cm, ax=ax)

