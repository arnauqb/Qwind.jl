using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind, YAML, DataFrames, CSV

@everywhere function generate_model(M, mdot, spin)
    config_path = "configs/config_base.yaml"
    config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
    config[:black_hole][:M] = M
    config[:black_hole][:mdot] = mdot
    config[:black_hole][:spin] = spin
    model = Model(config)
    return model
end

@everywhere function calculate_wind_mdots(model)
    rr = 10 .^ range(log10(6.0), log10(1500), length = 50)
    mdots = []
    zcs = []
    for r in rr
        zc, mdc = Qwind.find_nozzle_function_minimum(
            model.rt,
            model.bh,
            r,
            alpha = 0.6,
            zmax = 1e-2,
        )
        println(r)
        println(mdc)
        push!(zcs, zc)
        push!(mdots, mdc)
    end
    rr = rr[.!isnan.(mdots)]
    mdots = mdots[.!isnan.(mdots)]
    zcs = zcs[zcs .!= Inf]
    return rr, mdots, zcs
end

@everywhere function generate_and_save_df(M, mdot)
    model = generate_model(M, mdot, 0.0)
    rr, mdots, zcs = calculate_wind_mdots(model)
    df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    output_file = "src/initial_conditions/critical_points_data/M_$(M)_mdot_$(mdot).csv"
    CSV.write(output_file, df)
end

M_range = [1e6, 1e7, 1e8, 1e9, 1e10]
mdot_range = [0.025, 0.05, 0.075, 0.1, 0.25, 0.5]

for M in M_range
    f(mdot) = generate_and_save_df(M, mdot)
    pmap(f, mdot_range)
    #for mdot in mdot_range
    #    generate_and_save_df(M, mdot)
    #end
end
