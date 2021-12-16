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
