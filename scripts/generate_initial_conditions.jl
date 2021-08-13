using Distributed, ClusterManagers
pids = addprocs_slurm(32,
                      topology=:master_worker,
                      p="cosma6",
                      A="dp004",
                      t="04:00:00",
                      job_file_loc="cpu_logs")
@everywhere pushfirst!(Base.DEPOT_PATH, "/tmp/julia.cache")
println("Running on $(nprocs()) cores.")
@everywhere using LinearAlgebra, YAML, CSV, DataFrames
@everywhere BLAS.set_num_threads(1)
println("Qwind single node")
using Qwind, Printf
println("Done")
@info "Compiling Qwind..."
flush(stdout)
flush(stderr)

@info "Done"
flush(stdout)
flush(stderr)

@everywhere function generate_model(M, mdot, spin)
    #config_path = "configs/config_base.yaml"
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
    rr, mdots, zcs = calculate_wind_mdots(model.rt, molde.bh)
    df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    output_file = "/cosma6/data/dp004/dc-quer1/critical_points_mu/M_$(M)_mdot_$(mdot).csv"
    CSV.write(output_file, df)
end

M_range = [1e6, 1e7, 1e8, 1e9, 1e10]
mdot_range = [0.01] #10 .^ range(log10(0.025), log10(0.5), length=5) 

for M in M_range
    #f(mdot) = generate_and_save_df(M, mdot)
    #pmap(f, mdot_range)
    for mdot in mdot_range
        println("M $M mdot $mdot")
        flush(stdout)
        flush(stderr)
        generate_and_save_df(M, mdot)
    end
end
