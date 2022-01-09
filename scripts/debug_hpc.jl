ENV["JULIA_WORKER_TIMEOUT"] = 250
using Distributed, ClusterManagers
pids = addprocs_slurm(
    100,
    topology = :master_worker,
    p = "cosma7-shm",
    A = "dp004",
    t = "04:00:00",
    job_file_loc = "cpu_logs",
)

@everywhere using Pkg
@everywhere Pkg.activate("/cosma/home/dp004/dc-quer1/Qwind.jl")
#@everywhere pushfirst!(Base.DEPOT_PATH, "/tmp/julia.cache")
println("Running on $(nprocs()) cores.")
@everywhere using LinearAlgebra, YAML, CSV, DataFrames
@everywhere BLAS.set_num_threads(1)
@info "Compiling Qwind..."
flush(stdout)
flush(stderr)
@everywhere using Qwind, YAML
@info "Done"
flush(stdout)
flush(stderr)

using ProgressMeter, PyPlot
LogNorm = matplotlib.colors.LogNorm

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
run!(model, iterations_dict, parallel = true);

@everywhere function compute_xi(rad, r, z; include_scattering = true)
    n = get_density(rad.wi.density_grid, r, z)
    tau_x = compute_tau_xray(rad, ri = 0, phii = 0, zi = 0, rf = r, zf = z, phif = 0)
    ret = compute_ionization_parameter(
        r = r,
        z = z,
        vr = 0,
        vz = 0,
        number_density = n,
        tau_x = tau_x,
        xray_luminosity = rad.xray_luminosity,
        Rg = rad.bh.Rg,
        include_scattering = include_scattering,
        density_grid = rad.wi.density_grid,
        absorption_opacity = rad.xray_opacity,
        zh = rad.z_xray,
        mu_electron = rad.mu_electron,
        mu_nucleon = rad.mu_nucleon,
        scattered_luminosity_grid = rad.wi.scattered_lumin_grid,
    )
    return ret
end

@everywhere rr = collect(10 .^ range(log10(6), log10(1000), length = 200))
@everywhere zz = 10 .^ range(log10(1e-3), log10(1000), length = 201)
@everywhere rets = zeros(length(rr), length(zz))
@showprogress for (i, r) in enumerate(rr)
    rets[i, :] = pmap(
        z -> compute_xi(model.rad, r, z, include_scattering = true),
        zz,
        batch_size = 10,
    )
end


fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, rets', norm = LogNorm(vmin=1e-2, vmax=1e7), shading="auto")
cbar = plt.colorbar(cm, ax = ax)
ax.set_xlabel("R [Rg]")
ax.set_ylabel("R [Rg]")
cbar.set_label("Ionization parameter", rotation=270)
#fig.savefig("./tests_noscatt.pdf")
fig.savefig("./tests.pdf")

fig, ax = plt.subplots()
cm = ax.pcolormesh(model.rad.wi.density_grid.r_range, model.rad.wi.density_grid.z_range, model.rad.wi.density_grid.grid', norm=LogNorm())
fig.savefig("./density_tests.pdf")

@everywhere function f(i, model)
    ret = zeros(length(zz))
    r = rr[i]
    for (j, z) in enumerate(zz)
    end
    return ret
end
xi_scatt = @showprogress pmap(i -> f(i, model), 1:length(rr))
xi_scatt = hcat(xi_scatt)';

#iterations_dict[1] = Dict()
#integrators, streamlines = Qwind.run_integrators!(model, iterations_dict, it_num = 1, parallel = true);
