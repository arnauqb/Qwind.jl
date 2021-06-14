using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")

#mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("configs/debug.yaml");
iterations_dict = Dict();
#run!(model, iterations_dict)
