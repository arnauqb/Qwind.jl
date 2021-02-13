using CSV, DataFrames, YAML, JLD2
export save_integrator, save_integrators, save_wind

function save_integrator(data::Dict, save_path)
    if save_path === nothing
        return
    end
    append = false
    if isfile(save_path)
        ret = CSV.read(save_path, DataFrame)
        append = true
    end
    df = DataFrame()
    line_id = pop!(data, :line_id)
    df[!, :line_id] = line_id * ones(Int64, length(data[:r]))
    for key in keys(data)
        df[!, key] = data[key]
    end
    if append
        append!(ret, df)
    else
        ret = df
    end
    CSV.write(save_path, ret)
end

function save_integrators(integrators, save_path)
    if save_path === nothing
        return
    end
    df = DataFrame()
    for integrator in integrators
        df2 = DataFrame()
        data = integrator.p.data
        line_id = pop!(data, :line_id)
        df2[!, :line_id] = line_id * ones(Int64, length(data[:r]))
        for key in keys(data)
            df2[!, key] = data[key]
        end
        append!(df, df2)
    end
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

function save_wind_properties(integrators, save_path, bh::BlackHole)
    ret = Dict()
    ret["eddington_luminosity"] = compute_eddington_luminosity(bh)
    ret["bolometric_luminosity"] = compute_bolometric_luminosity(bh)
    ret["mass_accretion_rate"] = compute_mass_accretion_rate(bh)
    ret["kinetic_luminosity"] = compute_kinetic_luminosity(integrators, bh.Rg)
    ret["mass_loss"] = compute_wind_mdot(integrators, bh.Rg)
    ret["max_velocity"] = compute_maximum_velocity(integrators)
    ret["mass_loss_fraction"] = ret["mass_loss"] / ret["mass_accretion_rate"]
    YAML.write_file(save_path, ret)
end

function save_radiative_transfer(rt::RadiativeTransfer, save_path)
    @save save_path * "/radiative_transfer.jld2" rt
end

function save_wind(integrators, model, save_path)
    mkpath(save_path)
    lines_save_path = save_path * "/streamlines.csv"
    properties_save_path = save_path * "/wind_properties.yaml"
    save_integrators(integrators, lines_save_path)
    save_wind_properties(integrators, properties_save_path, model.bh)
    save_radiative_transfer(model.rt, save_path)
end


