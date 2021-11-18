export Radiation, update_radiation, compute_tau_uv, compute_tau_xray

struct Radiation{T<:AbstractFloat}
    bh::BlackHole{T}
    disk_grid::Vector{T}
    fuv_grid::Vector{T}
    mdot_grid::Vector{T}
    xray_luminosity::T
end

function Radiation(bh::BlackHole; nr, fx, fuv, disk_r_in, disk_r_out)
    disk_grid = 10 .^ range(log10(disk_r_in), log10(disk_r_out), length = nr)
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
    return Radiation(bh, disk_grid, uvf, mdot_grid, xray_luminosity)
end

# read from config
function Radiation(bh::BlackHole, parameters::Parameters)
    return Radiation(
        bh,
        nr = parameters.disk_nr,
        fx = parameters.f_x,
        fuv = parameters.f_uv,
        disk_r_in = parameters.disk_r_in,
        disk_r_out = parameters.disk_r_out,
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

#function set_tau_uv_calculation(radiation, tau_uv_calculation)
#    return Radiation(
#        radiation.bh,
#        radiation.wi,
#        radiation.disk_grid,
#        radiation.fuv_grid,
#        radiation.mdot_grid,
#        radiation.xray_luminosity,
#        radiation.disk_r_in,
#        radiation.z_xray,
#        radiation.z_disk,
#        radiation.mu_nucleon,
#        radiation.mu_electron,
#        radiation.relativistic,
#        radiation.uv_opacity,
#        radiation.xray_opacity,
#        radiation.xray_scattering,
#        tau_uv_calculation,
#        radiation.disk_integral_rtol,
#    )
#end

# Optical depths

function compute_tau_uv(
    density_grid::DensityGrid,
    iterator::GridIterator,
    radiation::Radiation,
    parameters::Parameters;
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    max_tau = 50,
)
    return compute_optical_depth(
        density_grid,
        iterator,
        parameters.uv_opacity_flag,
        ri = ri,
        phii = phii,
        zi = zi,
        rf = rf,
        phif = phif,
        zf = zf,
        max_tau = max_tau,
        mu_nucleon = parameters.mu_nucleon,
        mu_electron = parameters.mu_electron,
        Rg = radiation.bh.Rg,
        source_luminosity = radiation.xray_luminosity,
    )
end


function compute_tau_xray(
    density_grid::DensityGrid,
    iterator::GridIterator,
    radiation::Radiation,
    parameters::Parameters;
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    max_tau = 50,
)
    return compute_optical_depth(
        density_grid,
        iterator,
        parameters.xray_opacity_flag,
        ri = ri,
        phii = phii,
        zi = zi,
        rf = rf,
        phif = phif,
        zf = zf,
        max_tau = max_tau,
        mu_nucleon = parameters.mu_nucleon,
        mu_electron = parameters.mu_electron,
        Rg = radiation.bh.Rg,
        source_luminosity = radiation.xray_luminosity,
    )
end
