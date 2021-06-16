using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")

model = Model("configs/to_debug.yaml");
mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("configs/to_debug.yaml");
iterations_dict = Dict();
run!(model, iterations_dict)


integrators = iterations_dict[3]["integrators"];
QwindPlotting.plot_streamlines(integrators)

n0s = [integ.p.n0 for integ in integrators];
r0s = [integ.p.r0 for integ in integrators];
loglog(r0s, n0s)

rt = iterations_dict[2]["radiative_transfer"];
xl = rt.radiation.xray_luminosity;
Rg = rt.radiation.Rg;
QwindPlotting.plot_xray_grid(rt.interpolator.density_grid, xl, Rg, zmax=10)
