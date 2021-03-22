using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall, ProgressMeter
using Qwind
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();
iterations_dict[1] = Dict();
run!(model, iterations_dict, parallel = true)


using HDF5

c = h5open("./runs/tests/results.hdf5", "r") do file
    read(file, "iteration_001/trajectories")
end

ts = load_trajectories("./runs/tests/results.hdf5")
