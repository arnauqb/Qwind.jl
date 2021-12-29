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
@everywhere using Qwind, YAML, Printf
@info "Done"
flush(stdout)
flush(stderr)

@everywhere function generate_model(M, mdot, spin)
    config_path = "/cosma/home/dp004/dc-quer1/Qwind.jl/configs/config_base.yaml"
    config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
    config[:black_hole][:M] = M
    config[:black_hole][:mdot] = mdot
    config[:black_hole][:spin] = spin
    model = Model(config)
    return model
end

@everywhere function generate_and_save_df(M, mdot)
    model = generate_model(M, mdot, 0.0)
    rr, mdots, zcs = calculate_wind_mdots(
        model.wind.density_grid,
        model.wind.grid_iterator,
        model.parameters,
        model.rad,
        rmin = 6.1,
        rmax = 1500,
        nr = 100,
        nz = 250,
    )
    df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    mdot_string = @sprintf "%.4f" mdot
    output_file = "/cosma6/data/dp004/dc-quer1/critical_points/M_$(M)_mdot_$(mdot_string).csv"
    CSV.write(output_file, df)
end

M_range = [1e6, 1e7, 1e8, 1e9, 1e10]
mdot_range = 10 .^ range(log10(0.05), log10(0.5), length = 5)
mdot_range = vcat(mdot_range, [0.1, 0.25])
mdot_range = vcat(mdot_range, [0.075, 0.2, 0.3, 0.4])

for M in M_range
    for mdot in mdot_range
        println("M $M mdot $mdot")
        flush(stdout)
        flush(stderr)
        generate_and_save_df(M, mdot)
    end
end
