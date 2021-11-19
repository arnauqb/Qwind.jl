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
run!(model, iterations_dict)


compute_tau_xray(
    model.rad.wi.density_grid,
    Boost(),
    ri = 0,
    zi = 0,
    rf = 1000,
    zf = 1,
    xray_luminosity = model.rad.xray_luminosity,
    Rg = model.bh.Rg,
    mu_nucleon = model.rad.mu_nucleon,
    mu_electron = model.rad.mu_electron,
)

compute_tau_uv(
    model.rad.wi.density_grid,
    ri = 0,
    phii=0,
    zi = 0,
    rf = 1000,
    phif=0,
    zf = 1,
    Rg = model.bh.Rg,
    mu_electron = model.rad.mu_electron,
    max_tau=Inf,
)

