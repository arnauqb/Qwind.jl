using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

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
iterations_dict[1] = Dict()
integrators, streamlines = Qwind.run_integrators!(model, iterations_dict, it_num = 1, parallel = false);
wi = model.rad.wi;

#new_radiation = update_radiation(model.rad, streamlines)


hull = Hull(streamlines);
density_grid = update_density_grid(wi.density_grid, wi.update_grid_flag, streamlines, hull);
velocity_grid = update_velocity_grid(wi.velocity_grid, wi.update_grid_flag, streamlines, hull);


ionization_grid = IonizationGrid(
    density_grid,
    Rg = model.rad.bh.Rg,
    xray_luminosity = model.rad.xray_luminosity,
    z_xray = model.rad.z_xray,
)
