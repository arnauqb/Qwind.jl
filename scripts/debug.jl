using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf

function get_model(config)
    model = Model(config)
    try
        mv(model.config[:integrator][:save_path], "backup", force = true)
    catch
    end
    model = Model(config)
    iterations_dict = Dict()
    return model, iterations_dict
end

model, iterations_dict = get_model("./configs/debug.yaml");

run_iteration!(model, iterations_dict, it_num = 1, parallel = true);


@profile
function profile_integ()
    rr = range(20.0, 1000.0, length=50)
    zz = 10 .^ range(-6, 3.0, length=50)
    for r in rr
        for z in zz
            compute_disk_radiation_field(radiation, r, z, 0.0, 0.0)
        end
    end
end
