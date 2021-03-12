using CSV, DataFrames, YAML, HDF5
export save_integrator, save_integrators, save_wind

function create_integrators_df(integrators, Rg)
    df = DataFrame()
    for integrator in integrators
        df2 = DataFrame()
        data = integrator.p.data
        line_id = pop!(data, :line_id)
        df2[!, :line_id] = line_id * ones(Int64, length(data[:r]))
        for key in keys(data)
            df2[!, key] = data[key]
        end
        # add integrator momentum
        integrator_momentum = compute_integrator_momentum(integrator, Rg)
        df2[!, "momentum"] = integrator_momentum
        append!(df, df2)
    end
    return df
end

function save_integrators(integrators, save_path, Rg)
    if save_path === nothing
        return
    end
    df = create_integrators_df(integrators, Rg)
    CSV.write(save_path, df)
end

function compute_integrator_mdot(integrator, Rg)
    if escaped(integrator)
        n0 = integrator.p.n0
        v0 = integrator.p.v0
        lw = integrator.p.lwnorm * integrator.p.r0
        mw = 2π * integrator.p.r0 * lw * Rg^2 * n0 * v0 * C * M_P
        return mw
    else
        return 0.
    end
end

function compute_integrator_momentum(integrator, Rg)
    n0 = integrator.p.n0
    v0 = integrator.p.v0
    lw = integrator.p.lwnorm * integrator.p.r0
    mw = 2π * integrator.p.r0 * lw * Rg^2 * n0 * v0 * C * M_P
    vt = sqrt.(integrator.p.data[:vr] .^ 2 + integrator.p.data[:vz] .^ 2) * C
    return mw * vt
end


function compute_kinetic_luminosity(integrator, Rg)
    if !escaped(integrator)
        return 0.0
    end
    vrf = integrator.p.data[:vr][end]
    vzf = integrator.p.data[:vz][end]
    vf = sqrt(vrf^2 + vzf^2) * C
    kin_lumin = 0.5 * compute_integrator_mdot(integrator, Rg) * vf^2
    return kin_lumin
end

function compute_wind_mdot(integrators::Vector, Rg)
    mdot_wind = 0.0
    for i in 1:length(integrators)
        isassigned(integrators, i) || continue
        integrator = integrators[i]
        mdot_wind += compute_integrator_mdot(integrator, Rg)
    end
    return mdot_wind
end

function compute_kinetic_luminosity(integrators::Vector, Rg)
    kin_lumin = 0.0
    for i in 1:length(integrators)
        isassigned(integrators, i) || continue
        integrator = integrators[i]
        kin_lumin += compute_kinetic_luminosity(integrator, Rg)
    end
    return kin_lumin
end

function compute_maximum_velocity(integrators)
    maxv = 0.
    for i in 1:length(integrators)
        isassigned(integrators, i) || continue
        integrator = integrators[i]
        escaped(integrator) || continue
        vrf = integrator.p.data[:vr][end]
        vzf = integrator.p.data[:vz][end]
        vf = sqrt(vrf^2 + vzf^2)
        if vf > maxv
            maxv = vf
        end
    end
    return maxv
end

function create_wind_properties(integrators, bh::BlackHole)
    ret = Dict()
    ret["eddington_luminosity"] = compute_eddington_luminosity(bh)
    ret["bolometric_luminosity"] = compute_bolometric_luminosity(bh)
    ret["mass_accretion_rate"] = compute_mass_accretion_rate(bh)
    ret["kinetic_luminosity"] = compute_kinetic_luminosity(integrators, bh.Rg)
    ret["mass_loss"] = compute_wind_mdot(integrators, bh.Rg)
    ret["max_velocity"] = compute_maximum_velocity(integrators)
    ret["mass_loss_fraction"] = ret["mass_loss"] / ret["mass_accretion_rate"]
    return ret
end

function save_wind_properties(integrators, save_path, bh::BlackHole)
    ret = create_wind_properties(integrators, bh)
    YAML.write_file(save_path, ret)
    return ret
end

function save_wind(integrators, model, save_path, it_num)
    mkpath(save_path)
    iteration_save_path = save_path * "/iteration_$(@sprintf "%03d" it_num)"
    mkpath(iteration_save_path)
    hdf5_save_path = save_path * "/results.hdf5"
    save_hdf5(integrators, model, hdf5_save_path, it_num)
    lines_save_path = iteration_save_path * "/streamlines.csv"
    properties_save_path = iteration_save_path * "/wind_properties.yaml"
    save_integrators(integrators, lines_save_path, model.bh.Rg)
    properties = save_wind_properties(integrators, properties_save_path, model.bh)
    return properties
end


function save_hdf5(integrators, model, hdf5_save_path, it_num)
    iteration = @sprintf "iteration_%03d" it_num
    density_grid = model.rt.interpolator.density_grid
    bh = model.bh
    h5open(hdf5_save_path,isfile(hdf5_save_path) ? "r+" : "w") do file
        g = create_group(file, iteration)
        dg = create_group(g, "density_grid")
        dg["r"] = density_grid.r_range
        dg["z"] = density_grid.z_range
        if density_grid.grid === nothing
            grid_to_save = zeros((length(density_grid.r_range), length(density_grid.z_range)))
        else
            grid_to_save = density_grid.grid
        end
        dg["grid"] = grid_to_save
        g["eddington_luminosity"] = compute_eddington_luminosity(bh)
        g["bolometric_luminosity"] = compute_bolometric_luminosity(bh)
        macc = compute_mass_accretion_rate(bh)
        g["mass_accretion_rate"] = macc
        g["kinetic_luminosity"] = compute_kinetic_luminosity(integrators, bh.Rg)
        mloss = compute_wind_mdot(integrators, bh.Rg)
        g["mass_loss"] = mloss
        g["max_velocity"] = compute_maximum_velocity(integrators)
        g["mass_loss_fraction"] = mloss / macc
    end
end
