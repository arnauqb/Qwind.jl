using CSV, DataFrames, YAML, HDF5, Printf, JLD2
export save_wind


function save_integrators(integrators, save_path)
    @save save_path integrators
end

function compute_streamline_mdot(streamline::Streamline, Rg; mu_nucleon = 0.61)
    n0 = streamline.n[1]
    v0 = sqrt(streamline.vr[1]^2 + streamline.vz[1]^2)
    lw = streamline.width[1]
    mw = 2π * streamline.r[1] * lw * Rg^2 * n0 * v0 * C * M_P * mu_nucleon
    return mw
end

function compute_wind_mdot(streamlines::Streamlines, Rg; mu_nucleon = 0.61)
    mdot_wind = 0.0
    for streamline in streamlines
        if escaped(streamline)
            mdot_wind += compute_streamline_mdot(streamline, Rg, mu_nucleon = mu_nucleon)
        end
    end
    return mdot_wind
end


function compute_kinetic_luminosity(streamline::Streamline, Rg; mu_nucleon = 0.61)
    vrf = streamline.vr[end]
    vzf = streamline.vz[end]
    vf = sqrt(vrf^2 + vzf^2) * C
    kin_lumin =
        0.5 * compute_streamline_mdot(streamline, Rg, mu_nucleon = mu_nucleon) * vf^2
    return kin_lumin
end

function compute_kinetic_luminosity(streamlines::Streamlines, Rg; mu_nucleon = 0.61)
    kin_lumin = 0.0
    for streamline in streamlines
        if escaped(streamline)
            kin_lumin += compute_kinetic_luminosity(streamline, Rg, mu_nucleon = mu_nucleon)
        end
    end
    return kin_lumin
end

function compute_momentum_rate(streamline::Streamline, Rg; mu_nucleon = 0.61)
    vrf = streamline.vr[end]
    vzf = streamline.vz[end]
    vf = sqrt(vrf^2 + vzf^2) * C
    ret = compute_streamline_mdot(streamline, Rg, mu_nucleon = mu_nucleon) * vf
    return ret
end

function compute_momentum_rate(streamlines::Streamlines, Rg; mu_nucleon = 0.61)
    ret = 0.0
    for streamline in streamlines
        if escaped(streamline)
            ret += compute_momentum_rate(streamline, Rg, mu_nucleon = mu_nucleon)
        end
    end
    return ret
end

function compute_maximum_velocity(streamlines::Streamlines)
    maxv = 0.0
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

function create_wind_properties(
    streamlines::Streamlines,
    bh::BlackHole;
    mu_nucleon = 0.61,
    mu_electron = 1.17,
)
    ret = Dict()
    ret["eddington_luminosity"] =
        compute_eddington_luminosity(bh, mu_electron = mu_electron)
    ret["bolometric_luminosity"] =
        compute_bolometric_luminosity(bh, mu_electron = mu_electron)
    ret["mass_accretion_rate"] = compute_mass_accretion_rate(bh, mu_electron = mu_electron)
    ret["kinetic_luminosity"] =
        compute_kinetic_luminosity(streamlines, bh.Rg, mu_nucleon = mu_nucleon)
    ret["mass_loss"] = compute_wind_mdot(streamlines, bh.Rg, mu_nucleon = mu_nucleon)
    ret["max_velocity"] = compute_maximum_velocity(streamlines)
    ret["momentum"] = compute_momentum_rate(streamlines, bh.Rg, mu_nucleon = mu_nucleon)
    ret["mass_loss_fraction"] = ret["mass_loss"] / ret["mass_accretion_rate"]
    return ret
end

function save_wind_properties(
    streamlines,
    save_path,
    bh::BlackHole;
    mu_nucleon = 0.61,
    mu_electron = 1.17,
)
    ret = create_wind_properties(
        streamlines,
        bh,
        mu_nucleon = mu_nucleon,
        mu_electron = mu_electron,
    )
    YAML.write_file(save_path, ret)
    return ret
end

function save_wind(integrators, streamlines, model, save_path, it_num)
    mkpath(save_path)
    iteration_save_path = save_path * "/iteration_$(@sprintf "%03d" it_num)"
    mkpath(iteration_save_path)
    hdf5_save_path = save_path * "/results.hdf5"
    save_hdf5(integrators, streamlines, model, hdf5_save_path, it_num)
    properties_save_path = iteration_save_path * "/wind_properties.yaml"
    integrators_save_path = iteration_save_path * "/trajectories.jld2"
    save_integrators(integrators, integrators_save_path)
    properties = save_wind_properties(
        streamlines,
        properties_save_path,
        model.bh,
        mu_nucleon = model.rad.mu_nucleon,
        mu_electron = model.rad.mu_electron,
    )
    return properties
end

function save_density_grid!(density_grid::DensityGrid, group)
    dg = create_group(group, "density_grid")
    dg["r"] = density_grid.r_range
    dg["z"] = density_grid.z_range
    dg["nr"] = density_grid.nr
    dg["nz"] = density_grid.nz
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
    dg["nr"] = velocity_grid.nr
    dg["nz"] = velocity_grid.nz
    dg["vr_grid"] = velocity_grid.vr_grid
    dg["vphi_grid"] = velocity_grid.vphi_grid
    dg["vz_grid"] = velocity_grid.vz_grid
    return
end

function save_streamlines_and_trajectories!(integrators, streamlines, group)
    g = create_group(group, "trajectories")
    gg = create_group(group, "streamlines")
    for (i, integrator) in enumerate(integrators)
        trajectory = Trajectory(integrator)
        tgroup = create_group(g, "$i")
        tgroup["id"] = trajectory.id
        tgroup["t"] = trajectory.t
        tgroup["r"] = trajectory.r
        tgroup["z"] = trajectory.z
        tgroup["vr"] = trajectory.vr
        tgroup["vphi"] = trajectory.vphi
        tgroup["vz"] = trajectory.vz
        tgroup["n"] = trajectory.n
    end

    for (i, streamline) in enumerate(streamlines)
        sgroup = create_group(gg, "$i")
        sgroup["id"] = streamline.id
        sgroup["t"] = streamline.t
        sgroup["r"] = streamline.r
        sgroup["z"] = streamline.z
        sgroup["vr"] = streamline.vr
        sgroup["vphi"] = streamline.vphi
        sgroup["vz"] = streamline.vz
        sgroup["n"] = streamline.n
        sgroup["line_width"] = streamline.width
    end
    return
end

function save_wind_hull!(hull::ConcaveHull.Hull, group)
    g = create_group(group, "wind_hull")
    g["k"] = hull.k
    g["vertices_r"] = [v[1] for v in hull.vertices]
    g["vertices_z"] = [v[2] for v in hull.vertices]
    #g["converged"] = hull.converged
end
save_wind_hull!(hull::Nothing, group) = nothing


function save_hdf5(integrators, streamlines, model, hdf5_save_path, it_num)
    iteration = @sprintf "iteration_%03d" it_num
    density_grid = model.rad.wi.density_grid
    velocity_grid = model.rad.wi.velocity_grid
    wind_hull = model.rad.wi.wind_hull
    bh = model.bh
    h5open(hdf5_save_path, isfile(hdf5_save_path) ? "r+" : "w") do file
        g = create_group(file, iteration)
        save_density_grid!(density_grid, g)
        save_velocity_grid!(velocity_grid, g)
        save_streamlines_and_trajectories!(integrators, streamlines, g)
        save_wind_hull!(wind_hull, g)
        g["eddington_luminosity"] = compute_eddington_luminosity(bh)
        g["bolometric_luminosity"] = compute_bolometric_luminosity(bh)
        macc = compute_mass_accretion_rate(bh)
        g["mass_accretion_rate"] = macc
        g["kinetic_luminosity"] = compute_kinetic_luminosity(streamlines, bh.Rg)
        g["momentum"] = compute_momentum_rate(streamlines, bh.Rg)
        mloss = compute_wind_mdot(streamlines, bh.Rg)
        g["mass_loss"] = mloss
        g["max_velocity"] = compute_maximum_velocity(streamlines)
        g["mass_loss_fraction"] = mloss / macc
    end
end
