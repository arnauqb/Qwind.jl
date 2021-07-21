using Distributed, ClusterManagers
pids = addprocs_slurm(32,
                      topology=:master_worker,
                      p="cosma6",
                      A="dp004",
                      t="24:00:00",
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

@everywhere function calculate_wind_mdots(model)
    rr = 10 .^ range(log10(15.0), log10(1500), length = 50)
    mdots = []
    zcs = []
    f(r) = Qwind.find_nozzle_function_minimum(
            model.rt,
            model.bh,
            r,
            alpha = 0.6,
            zmax = 1e-2,
        )
    results = pmap(f, rr)
    zcs = [res[1] for res in results]
    mdots = [res[2] for res in results]
    #for r in rr
    #    zc, mdc =         println(r)
    #    println(mdc)
    #    push!(zcs, zc)
    #    push!(mdots, mdc)
    #end
    rr = rr[.!isnan.(mdots)]
    mdots = mdots[.!isnan.(mdots)]
    zcs = zcs[zcs .!= Inf]
    return rr, mdots, zcs
end

@everywhere function generate_and_save_df(M, mdot)
    model = generate_model(M, mdot, 0.0)
    rr, mdots, zcs = calculate_wind_mdots(model)
    df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    #output_file = "src/initial_conditions/critical_points_data/M_$(M)_mdot_$(mdot).csv"
    #output_file = "/cosma/home/dp004/dc-quer1/critical_points_mu/M_$(M)_mdot_$(mdot).csv"
    output_file = "/cosma6/data/dp004/dc-quer1/critical_points_mu/M_$(M)_mdot_$(mdot).csv"
    CSV.write(output_file, df)
end

M_range = [1e6, 1e7, 1e8, 1e9, 1e10]
mdot_range = 10 .^ range(log10(0.025), log10(0.5), length=5) 

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
