using CSV, DataFrames, YAML, HDF5, Printf, JLD2
export save_wind


function save_integrators(integrators, save_path)
    @save save_path integrators
end

function compute_streamline_mdot(streamline::Streamline, Rg)
    if escaped(streamline)
        n0 = streamline.n[1]
        v0 = sqrt(streamline.vr[1]^2 + streamline.vz[1]^2)
        lw = streamline.width
        mw = 2π * streamline.r[1] * lw * Rg^2 * n0 * v0 * C * M_P * 0.61
        return mw
    else
        return 0.
    end
end

function compute_streamline_momentum(streamline::Streamline, Rg)
    n0 = streamline.n[1]
    v0 = sqrt(streamline.vr[1]^2 + streamline.vz[1]^2)
    lw = streamline.width
    mw = 2π * streamline.r[1] * lw * Rg^2 * n0 * v0 * C * M_P * 0.61
    vt = sqrt.(streamline.vr .^ 2 + streamline.vz .^ 2) * C
    return mw * vt
end


function compute_kinetic_luminosity(streamline::Streamline, Rg)
    if !escaped(streamline)
        return 0.0
    end
    vrf = streamline.vr[end]
    vzf = streamline.vz[end]
    vf = sqrt(vrf^2 + vzf^2) * C
    kin_lumin = 0.5 * compute_integrator_mdot(streamline, Rg) * vf^2
    return kin_lumin
end

function compute_wind_mdot(streamlines::Streamlines, Rg)
    mdot_wind = 0.0
    for streamline in streamlines
        mdot_wind += compute_integrator_mdot(streamline, Rg)
    end
    return mdot_wind
end

function compute_kinetic_luminosity(streamlines::Streamlines, Rg)
    kin_lumin = 0.0
    for streamline in streamlines
        kin_lumin += compute_kinetic_luminosity(streamline, Rg)
    end
    return kin_lumin
end

function compute_momentum_density(streamline::Streamline, Rg)
    if !escaped(streamline)
        return 0.0
    end
    vrf = streamline.vr[end]
    vzf = streamline.vz[end]
    vf = sqrt(vrf^2 + vzf^2) * C
    ret = compute_streamline_mdot(streamline, Rg) * vf
    return ret
end

function compute_momentum_density(streamlines::Streamlines, Rg)
    ret = 0.0
    for streamline in streamlines
        ret += compute_momentum_density(streamline, Rg)
    end
    return ret 
end

function compute_maximum_velocity(streamlines::Streamlines)
    maxv = 0.
    for streamline in streamlines
        escaped(streamline) || continue
        vrf = streamline.vr[end]
        vzf = streamline.vz[end]
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
    ret["momentum"] = compute_momentum_density(integrators, bh.Rg)
    ret["mass_loss_fraction"] = ret["mass_loss"] / ret["mass_accretion_rate"]
    return ret
end

function save_wind_properties(integrators, save_path, bh::BlackHole)
    ret = create_wind_properties(integrators, bh)
    YAML.write_file(save_path, ret)
    return ret
end

function save_wind(integrators, streamlines, model, save_path, it_num)
    mkpath(save_path)
    iteration_save_path = save_path * "/iteration_$(@sprintf "%03d" it_num)"
    mkpath(iteration_save_path)
    hdf5_save_path = save_path * "/results.hdf5"
    save_hdf5(integrators, model, hdf5_save_path, it_num)
    properties_save_path = iteration_save_path * "/wind_properties.yaml"
    integrators_save_path = iteration_save_path * "/trajectories.jld2"
    save_integrators(integrators, integrators_save_path)
    properties = save_wind_properties(integrators, properties_save_path, model.bh)
    return properties
end

function save_density_grid!(density_grid::DensityGrid, group)
    dg = create_group(group, "density_grid")
    dg["r"] = density_grid.r_range
    dg["z"] = density_grid.z_range
    if density_grid.grid === nothing
        grid_to_save = zeros((length(density_grid.r_range), length(density_grid.z_range)))
    else
        grid_to_save = density_grid.grid
    end
    dg["grid"] = grid_to_save
    return
end

function save_velocity_grid!(velocity_grid::VelocityGrid, group)
    dg = create_group(group, "velocity_grid")
    dg["r"] = velocity_grid.r_range
    dg["z"] = velocity_grid.z_range
    dg["vr_grid"] = velocity_grid.vr_grid
    dg["vphi_grid"] = velocity_grid.vphi_grid
    dg["vz_grid"] = velocity_grid.vz_grid
    return
end

function save_streamlines!(integrators::Vector{<:Sundials.IDAIntegrator}, group)
    g = create_group(group, "trajectories")
    for (i, integrator) in enumerate(integrators)
        integ = Streamline(integrator)
        tgroup = create_group(g, "$i")
        tgroup["id"] = integ.id
        tgroup["t"] = integ.t
        tgroup["r"] = integ.r
        tgroup["z"] = integ.z
        tgroup["vr"] = integ.vr
        tgroup["vphi"] = integ.vphi
        tgroup["vz"] = integ.vz
        tgroup["n"] = integ.n
    end
    return
end

function save_wind_hull!(hull::ConcaveHull.Hull, group)
    g = create_group(group, "wind_hull")
    g["k"] = hull.k
    g["vertices_r"] = [v[1] for v in hull.vertices]
    g["vertices_z"] = [v[2] for v in hull.vertices]
    g["converged"] = hull.converged
end
save_wind_hull!(hull::Nothing, group) = nothing


function save_hdf5(streamlines, model, hdf5_save_path, it_num)
    iteration = @sprintf "iteration_%03d" it_num
    density_grid = model.rad.wi.density_grid
    velocity_grid = model.rad.wi.velocity_grid
    wind_hull = model.rad.wi.wind_hull
    bh = model.bh
    h5open(hdf5_save_path,isfile(hdf5_save_path) ? "r+" : "w") do file
        g = create_group(file, iteration)
        save_density_grid!(density_grid, g)
        save_velocity_grid!(velocity_grid, g)
        save_streamlines!(streamlines, g)
        save_wind_hull!(wind_hull, g)
        g["eddington_luminosity"] = compute_eddington_luminosity(bh)
        g["bolometric_luminosity"] = compute_bolometric_luminosity(bh)
        macc = compute_mass_accretion_rate(bh)
        g["mass_accretion_rate"] = macc
        g["kinetic_luminosity"] = compute_kinetic_luminosity(streamlines, bh.Rg)
        g["momentum"] = compute_momentum_density(streamlines, bh.Rg)
        mloss = compute_wind_mdot(streamlines, bh.Rg)
        g["mass_loss"] = mloss
        g["max_velocity"] = compute_maximum_velocity(streamlines)
        g["mass_loss_fraction"] = mloss / macc
    end
end
