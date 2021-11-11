export Radiation, update_radiation

struct Radiation{T<:AbstractFloat}
    bh::BlackHole{T}
    wi::WindInterpolator{T}
    disk_grid::Vector{T}
    fuv_grid::Vector{T}
    mdot_grid::Vector{T}
    xray_luminosity::T
    disk_r_in::T
    z_xray::T
    z_disk::T
    mu_nucleon::T
    mu_electron::T
    relativistic::RelativisticFlag
    uv_opacity::OpacityFlag
    xray_opacity::OpacityFlag
    xray_scattering::Bool
    tau_uv_calculation::TauUVCalculationFlag
    disk_integral_rtol::T
end

function Radiation(
    bh::BlackHole,
    wi::WindInterpolator;
    nr,
    fx,
    fuv,
    disk_r_in,
    z_xray,
    disk_height,
    mu_nucleon,
    mu_electron,
    relativistic,
    uv_opacity,
    xray_opacity,
    xray_scattering,
    tau_uv_calculation,
    disk_integral_rtol = 1e-3,
)
    rmin = bh.isco
    rmax = 1400.0
    disk_grid = 10 .^ range(log10(rmin), log10(rmax), length = nr)
    if fuv == "auto"
        uvf = uv_fractions(bh, disk_grid)
    else
        uvf = fuv .* ones(length(disk_grid))
    end
    if any(isnan.(uvf))
        error("UV fractions contain NaN, check radiation and boundaries")
    end
    mdot_grid = bh.mdot .* ones(length(disk_grid))
    xray_luminosity = fx * compute_bolometric_luminosity(bh)
    return Radiation(
        bh,
        wi,
        disk_grid,
        uvf,
        mdot_grid,
        xray_luminosity,
        disk_r_in,
        z_xray,
        disk_height,
        mu_nucleon,
        mu_electron,
        relativistic,
        uv_opacity,
        xray_opacity,
        xray_scattering,
        tau_uv_calculation,
        disk_integral_rtol,
    )
end

# read from config
function Radiation(bh::BlackHole, config::Dict)
    rc = config[:radiation]
    if rc[:relativistic]
        rel = Relativistic()
    else
        rel = NoRelativistic()
    end
    disk_r_in = rc[:disk_r_in]
    if disk_r_in == "isco"
        disk_r_in = bh.isco
    elseif disk_r_in == "r_in"
        disk_r_in = config[:initial_conditions][:r_in]
    end
    tau_uv_calc = rc[:tau_uv_calculation]
    if tau_uv_calc == "center"
        tau_uv_calc = TauUVCenter()
    elseif tau_uv_calc == "disk"
        tau_uv_calc = TauUVDisk()
    elseif tau_uv_calc == "no_tau_uv"
        tau_uv_calc = NoTauUV()
    else
        error("tau uv calc not recognised")
    end
    if rc[:xray_opacity] == "thomson"
        xray_opacity = ThomsonOpacity()
    else
        xray_opacity = BoostOpacity()
    end
    if rc[:uv_opacity] == "thomson"
        uv_opacity = ThomsonOpacity()
    else
        uv_opacity = BoostOpacity()
    end
    disk_rtol = get(rc, :disk_integral_rtol, 1e-3)
    wi = WindInterpolator(rc[:wind_interpolator])
    return Radiation(
        bh,
        wi,
        nr = rc[:n_r],
        fx = rc[:f_x],
        fuv = rc[:f_uv],
        disk_r_in = disk_r_in,
        z_xray = rc[:z_xray],
        disk_height = rc[:disk_height],
        mu_nucleon = get(rc, :mu_nucleon, 0.61),
        mu_electron = get(rc, :mu_electron, 1.17),
        relativistic = rel,
        uv_opacity = uv_opacity,
        xray_opacity = xray_opacity,
        xray_scattering = rc[:xray_scattering],
        tau_uv_calculation = tau_uv_calc,
        disk_integral_rtol = disk_rtol,
    )
end

# quick access to BH functions
disk_nt_rel_factors(radiation::Radiation, r) = disk_nt_rel_factors(radiation.bh, r)
compute_eddington_luminosity(radiation::Radiation) =
    compute_eddington_luminosity(radiation.bh)
compute_bolometric_luminosity(radiation::Radiation) =
    compute_bolometric_luminosity(radiation.bh)
get_Rg(radiation::Radiation) = radiation.bh.Rg
get_efficiency(radiation::Radiation) = radiation.bh.efficiency
get_spin(radiation::Radiation) = radiation.bh.spin
get_isco(radiation::Radiation) = radiation.bh.isco



# Disk functions

"""
Returns the fraction of UV emitted at disk anuli at position r,
and the local accretion rate.
"""
function get_fuv_mdot(radiation::Radiation, r)
    r_index = searchsorted_nearest(radiation.disk_grid, r)
    f_uv = radiation.fuv_grid[r_index]
    mdot = radiation.mdot_grid[r_index]
    return f_uv, mdot
end

# Quick access to wind interpolator functions

get_density(radiation::Radiation, r, z) = get_density(radiation.wi, r, z)

function update_radiation(radiation::Radiation, streamlines::Streamlines)
    @info "Updating radiation... "
    flush()
    #new_interp = update_wind_interpolator(radiation.wi, streamlines)
    wi = radiation.wi
    if maximum(wi.density_grid.grid) == wi.vacuum_density
        # first iteration, do not average
        new_interp = WindInterpolator(
            streamlines,
            nr = wi.density_grid.nr,
            nz = wi.density_grid.nz,
            vacuum_density = wi.vacuum_density,
            update_grid_flag = wi.update_grid_flag,
            xray_luminosity = radiation.xray_luminosity,
            z_xray = radiation.z_xray,
            mu_nucleon = radiation.mu_nucleon,
            mu_electron = radiation.mu_electron,
            Rg = radiation.bh.Rg,
        )
        return update_radiation(radiation, new_interp)
    end
    hull = Hull(streamlines)
    density_grid =
        update_density_grid(wi.density_grid, wi.update_grid_flag, streamlines, hull)
    velocity_grid =
        update_velocity_grid(wi.velocity_grid, wi.update_grid_flag, streamlines, hull)
    ionization_grid = IonizationGrid(
        density_grid,
        Rg = radiation.bh.Rg,
        xray_luminosity = radiation.xray_luminosity,
        z_xray = radiation.z_xray,
    )
    new_interp = WindInterpolator(
        hull,
        density_grid,
        velocity_grid,
        ionization_grid,
        wi.vacuum_density,
        wi.update_grid_flag,
    )
    return update_radiation(radiation, new_interp)
end

function update_radiation(radiation::Radiation, dgrid::DensityGrid)
    new_interp = update_wind_interpolator(radiation.wi, dgrid)
    return update_radiation(radiation, new_interp)
end

function update_radiation(radiation::Radiation, wind_interpolator::WindInterpolator)
    return Radiation(
        radiation.bh,
        wind_interpolator,
        radiation.disk_grid,
        radiation.fuv_grid,
        radiation.mdot_grid,
        radiation.xray_luminosity,
        radiation.disk_r_in,
        radiation.z_xray,
        radiation.z_disk,
        radiation.mu_nucleon,
        radiation.mu_electron,
        radiation.relativistic,
        radiation.uv_opacity,
        radiation.xray_opacity,
        radiation.xray_scattering,
        radiation.tau_uv_calculation,
        radiation.disk_integral_rtol,
    )
end

function set_tau_uv_calculation(radiation, tau_uv_calculation)
    return Radiation(
        radiation.bh,
        radiation.wi,
        radiation.disk_grid,
        radiation.fuv_grid,
        radiation.mdot_grid,
        radiation.xray_luminosity,
        radiation.disk_r_in,
        radiation.z_xray,
        radiation.z_disk,
        radiation.mu_nucleon,
        radiation.mu_electron,
        radiation.relativistic,
        radiation.uv_opacity,
        radiation.xray_opacity,
        radiation.xray_scattering,
        tau_uv_calculation,
        radiation.disk_integral_rtol,
    )
end

# Optical depths

function compute_tau_uv(radiation::Radiation; ri, phii, zi, rf, phif, zf, max_tau = 50)
    return compute_optical_depth(
        radiation.wi.density_grid.iterator,
        radiation.wi.density_grid,
        radiation.wi.ionization_grid,
        radiation.uv_opacity,
        ri = ri,
        phii = phii,
        zi = zi,
        rf = rf,
        phif = phif,
        zf = zf,
        max_tau = max_tau,
        mu_nucleon = radiation.mu_nucleon,
        mu_electron = radiation.mu_electron,
        Rg = radiation.bh.Rg,
    )
end

function compute_tau_xray(radiation::Radiation; ri, phii, zi, rf, phif, zf, max_tau = 50)
    return compute_optical_depth(
        radiation.wi.density_grid.iterator,
        radiation.wi.density_grid,
        radiation.wi.ionization_grid,
        radiation.xray_opacity,
        ri = ri,
        phii = phii,
        zi = zi,
        rf = rf,
        phif = phif,
        zf = zf,
        max_tau = max_tau,
        mu_nucleon = radiation.mu_nucleon,
        mu_electron = radiation.mu_electron,
        Rg = radiation.bh.Rg,
    )
end
