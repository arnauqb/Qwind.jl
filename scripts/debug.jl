using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf

function get_model(config)
    model = Model(config);
    try
        mv(model.config[:integrator][:save_path], "backup",  force=true)
    catch
    end
    model = Model(config);
    iterations_dict = Dict();
    return model, iterations_dict
end
model, iterations_dict = get_model("./configs/debug.yaml");

run_iteration!(model, iterations_dict, it_num=1, parallel=true)

run_iteration!(model, iterations_dict, it_num=2, parallel=true)

#using JLD2
#@save "./it_dict.jld2" iterations_dict

density_grid = iterations_dict[2]["radiative_transfer"].interpolator.density_grid;
rt = iterations_dict[2]["radiative_transfer"];



fig, ax = plt.subplots()
integs1 = iterations_dict[1]["integrators"];
integs2 = iterations_dict2[1]["integrators"];
for integ in integs1
    #ax.plot(integ.p.data[:r], integ.p.data[:z], color = "C0")
end
for integ in integs2
    ax.plot(integ.p.data[:r], integ.p.data[:z], color = "C1")
end

function profile_tauuv(grid, iterator, Rg)
    #r_range = 10 .^ range(-3, 2, length=25)
    #z_range = 10 .^ range(-3, 2, length=25)
    #for rp in r_range
    #    for zp in z_range
    #        #compute_disc_radiation_field(rt, rp, zp, 0.0, 0.0)
    #        compute_disc_radiation_field(rt, rp, zp, 0.0, 0.0)
    #    end
    #end
    rp_range_test = range(6.0, 1500.0, length = 20)
    r_range_test = range(6.0, 1500.0, length = 20)
    z_range_test = 10 .^ range(-6, 3.0, length = 20)
    phid_range_test = range(0.0, π, length = 20)
    for rdp in rp_range_test
        for rp in r_range_test
            for zp in z_range_test
                for pd in phid_range_test
                    qwsol = compute_uv_tau(grid, iterator, rdp, pd, 0.0, rp, 0.0, zp, Rg)
                end
            end
        end
    end
end

Profile.clear()
@profile profile_tauuv(density_grid, density_grid.iterator, model.bh.Rg)
pprof()

#run!(model, iterations_dict, parallel=true)
Profile.clear()
@profile run_iteration!(model, iterations_dict, it_num=2, parallel=false)
pprof()

lr, lw = Qwind.compute_lines_range(model);

fig, ax = plt.subplots()
for r in lr
    ax.axvline(r)
end
