import Qsosed
export Radiation, update_radiation, compute_tau_uv, compute_tau_xray

struct Radiation{T<:AbstractFloat}
    bh::BlackHole{T}
    disk_grid::Vector{T}
    fuv_grid::Vector{T}
    mdot_grid::Vector{T}
    xray_luminosity::T
    scattered_lumin_grid::ScatteredLuminosityGrid
    qsosed_model::Qsosed.QsosedModel
end

function Radiation(
    bh::BlackHole;
    disk_nr,
    fx,
    fuv,
    disk_r_in,
    disk_r_out,
    nr,
    nz,
    mu_nucleon,
)
    qsosed_parameters = Qsosed.Parameters(
        M = bh.M / M_SUN,
        mdot = bh.mdot,
        spin = bh.spin,
        mu_nucleon = mu_nucleon,
    )
    qsosed_model = Qsosed.QsosedModel(qsosed_parameters)
    if disk_r_in == "corona_radius"
        disk_r_in = qsosed_model.corona.radius
    elseif disk_r_in == "warm_radius"
        disk_r_in = qsosed_model.warm.radius
    end
    disk_grid = 10 .^ range(log10(disk_r_in), log10(disk_r_out), length = disk_nr)
    if fuv == "qsosed"
        _, uvf =
            Qsosed.radial_uv_fraction.(
                Ref(qsosed_model),
                r_min = disk_r_in,
                r_max = disk_r_out,
                n_r = disk_nr,
            )
    elseif fuv == "disk"
        uvf = Qsosed.disk_uv_fraction.(Ref(bh), disk_grid)
    else
        uvf = fuv .* ones(length(disk_grid))
    end
    if fx == "qsosed"
        fx = Qsosed.total_xray_fraction(qsosed_model)
    end
    if any(isnan.(uvf))
        error("UV fractions contain NaN, check radiation and boundaries")
    end
    mdot_grid = bh.mdot .* ones(length(disk_grid))
    xray_luminosity = fx * compute_bolometric_luminosity(bh)
    scatt_lumin = ScatteredLuminosityGrid(nr, nz, 0.0)
    return Radiation(
        bh,
        disk_grid,
        uvf,
        mdot_grid,
        xray_luminosity,
        scatt_lumin,
        qsosed_model,
    )
end

# read from config
function Radiation(bh::BlackHole, parameters::Parameters)
    return Radiation(
        bh,
        disk_nr = parameters.disk_nr,
        fx = parameters.f_x,
        fuv = parameters.f_uv,
        disk_r_in = parameters.disk_r_in,
        disk_r_out = parameters.disk_r_out,
        nr = parameters.radiation_grid_nr,
        nz = parameters.radiation_grid_nz,
        mu_nucleon = parameters.mu_nucleon,
    )
end
function Radiation(bh::BlackHole, config::Dict)
    parameters = Parameters(config)
    return Radiation(bh, parameters)
end

function update_radiation(radiation::Radiation, new_wind::Wind, parameters)
    new_lumin_grid = ScatteredLuminosityGrid(
        new_wind.density_grid,
        new_wind.grid_iterator,
        Rg = radiation.bh.Rg,
        source_luminosity = radiation.xray_luminosity,
        source_position = [0, 0, parameters.z_xray],
    )
    return Radiation(
        radiation.bh,
        radiation.disk_grid,
        radiation.fuv_grid,
        radiation.mdot_grid,
        radiation.xray_luminosity,
        new_lumin_grid,
        radiation.qsosed_model,
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
        Rg = radiation.bh.Rg,
        source_luminosity = radiation.xray_luminosity,
    )
end
