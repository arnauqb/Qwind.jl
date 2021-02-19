using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML
include("scripts/plotting.jl")
matplotlib.rcParams["figure.dpi"] = 300

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();

using JLD2, Profile, PProf
@load "big_rt.jld2" rt

iterations_dict[1] = Dict();
model.rt = rt;
iterations_dict[1]["radiative_transfer"] = rt;


do_iteration!(model, iterations_dict, it_num=1);

lines_range, lines_widths = get_initial_radii_and_linewidths(model.ic, model.bh.Rg);
i = 100
integrator = Qwind.create_and_run_integrator(
    model,
    linewidth = lines_widths[i],
    r0 = lines_range[i],
    atol = model.config[:integrator][:atol],
    rtol = model.config[:integrator][:rtol],
    line_id = i,
)


#run!(model, iterations_dict1)
#do_iteration!(model, iterations_dict, it_num=1);
#do_iteration!(model, iterations_dict, it_num=2);

Profile.clear()
@profile rt = update_radiative_transfer(model.rt, iterations_dict[1]["integrators"]);
pprof()

rt = iterations_dict[2]["radiative_transfer"];

#using JLD2
#@save "big_rt.jld2" rt
