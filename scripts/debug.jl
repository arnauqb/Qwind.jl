using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames
include("scripts/plotting.jl")
matplotlib.rcParams["figure.dpi"] = 300

mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("results/resolution_tests/config1.yaml");
iterations_dict = Dict();
run!(model, iterations_dict)
