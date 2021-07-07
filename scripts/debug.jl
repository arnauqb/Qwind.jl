using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")

model = Model("configs/debug.yaml");
mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("configs/debug.yaml");
iterations_dict = Dict();
run!(model, iterations_dict, parallel=true)

lr, lw = Qwind.compute_lines_range(model);

fig, ax = plt.subplots()
for r in lr
    ax.axvline(r)
end
