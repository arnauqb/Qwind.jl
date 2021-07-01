using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using Qwind
using YAML, Profile, PProf, PyCall
#include("scripts/plotting.jl")

configs = YAML.load_file("./bh_scan/all_configs.yaml", dicttype = Dict{Symbol,Any})
for (key, config) in zip(keys(configs), values(configs))
    println(key)
    if key in [:model_007, :model_008, :model_009, :model_003, :model_006, :model_002]
        continue
    end
    model = Model(config)
    iterations_dict = Dict();
    run!(model, iterations_dict, parallel=true)
end
