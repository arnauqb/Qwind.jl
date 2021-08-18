using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools

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


function profile_integ(model)
    rr = range(20.0, 1000.0, length=10)
    zz = 10 .^ range(-6, 3.0, length=10)
    for r in rr
        for z in zz
            compute_disc_radiation_field(model, r, z, 0.0, 0.0)
        end
    end
end

disable_timer!(timer)
            
enable_timer!(timer)

reset_timer!(timer); compute_disc_radiation_field(model.rad, r=100.0, z=1.0, vr=0.0, vz=0.0); timer


@btime compute_disc_radiation_field(model.rad, r=100.0, z=1.0, vr=0.0, vz=0.0)

profile_integ(model)

Profile.clear()
@profile profile_integ(model.rad)
pprof()
