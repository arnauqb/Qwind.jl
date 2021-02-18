using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML
include("scripts/plotting.jl")
matplotlib.rcParams["figure.dpi"] = 300

try
    mv(model1.config[:integrator][:save_path], "backup", force = true)
catch
end
model1 = Model("results/resolution_tests/config1.yaml");
iterations_dict1 = Dict();
run!(model1, iterations_dict1)
#do_iteration!(model1, iterations_dict1, it_num=1);

