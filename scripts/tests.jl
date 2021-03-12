using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using Qwind
using YAML, Profile, PProf, PyCall
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();
run!(model, iterations_dict)

# Velocity interpolation
