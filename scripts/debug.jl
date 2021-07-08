using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf

model = Model("configs/debug.yaml");
mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("configs/debug.yaml");

iterations_dict = Dict();
run_iteration!(model, iterations_dict, it_num=1, parallel=true)

run_iteration!(model, iterations_dict, it_num=2, parallel=true)

integs = iterations_dict[1]["integrators"];

times = Qwind.get_intersection_times(integs);


QwindPlotting.plot_streamlines(integs)

density_grid = iterations_dict[2]["radiative_transfer"].interpolator.density_grid;
rt = iterations_dict[2]["radiative_transfer"];

function profile_tauuv(rt)
    r_range = 10 .^ range(-3, 2, length=25)
    z_range = 10 .^ range(-3, 2, length=25)
    for rp in r_range
        for zp in z_range
            compute_disc_radiation_field(rt, rp, zp, 0.0, 0.0)
        end
    end
end

Profile.clear()
@profile profile_tauuv(rt)
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
