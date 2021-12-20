using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind

using PyPlot
#using YAML, HDF5, CSV, DataFrames
include("scripts/plotting.jl")
#using Profile, PProf, TimerOutputs, BenchmarkTools, ProgressMeter
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

Qwind.CAK_Σ(model.wind.density_grid, model.wind.grid_iterator, model.parameters, model.rad, 100, K=0.03)

fig, ax = plt.subplots()
r_range = 10 .^ range(log10(6.1), log10(1500), length=50);
n0 = getn0.(Ref(model), r_range);
ax.loglog(r_range, n0, "o-")

run!(model, iterations_dict, parallel = true)
