export QsosedRadiation, from_quadtree, get_fuv_mdot, relativistic_correction

struct QsosedRadiation{T} <: Radiation{T}
    disk_grid::Vector{T}
    fuv_grid::Vector{T}
    mdot_grid::Vector{T}
    xray_luminosity::T
    efficiency::T
    spin::T
    isco::T
    disk_r_in::T
    z_xray::T
    Rg::T
    flux_correction::FluxCorrection
end

QsosedRadiation() =
    QsosedRadiation([0.0], [0.0], [0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Relativistic())

function compute_ionization_parameter(
    radiation::QsosedRadiation,
    r,
    z,
    number_density,
    tau_x,
)
    d2 = (r^2 + (z - radiation.z_xray)^2) * radiation.Rg^2
    return max(radiation.xray_luminosity * exp(-tau_x) / (number_density * d2), 1e-20)
end


flux_correction(::Relativistic, beta) = ((1 - beta) / (1 + beta))^2
flux_correction(::NoRelativistic, beta) = 1.0

function compute_mass_accretion_rate(radiation::Radiation, r)
    r_idx = searchsorted_nearest(radiation.disk_grid, r)
    mdot = radiation.mdot_grid[r_idx]
    return mdot * 4 * π * radiation.Rg * M_P * C / (SIGMA_T * radiation.efficiency)
end


function QsosedRadiation(
    bh::BlackHole,
    nr::Int,
    fx::Float64,
    disk_r_in::Float64,
    z_xray::Float64,
    relativistic::FluxCorrection,
)
    rmin = bh.isco
    rmax = 1400.0
    disk_grid = 10 .^ range(log10(rmin), log10(rmax), length = nr)
    uvf = uv_fractions(bh, disk_grid)
    if any(isnan.(uvf))
        error("UV fractions contain NaN, check radiation and boundaries")
    end
    mdot_grid = bh.mdot .* ones(length(disk_grid))
    xray_luminosity = fx * compute_bolometric_luminosity(bh)
    return QsosedRadiation(
        disk_grid,
        uvf,
        mdot_grid,
        xray_luminosity,
        bh.efficiency,
        bh.spin,
        bh.isco,
        disk_r_in,
        z_xray,
        bh.Rg,
        relativistic,
    )
end
function QsosedRadiation(bh::BlackHole, config::Dict)
    radiation_config = config[:radiation]
    relativistic = radiation_config[:relativistic]
    if typeof(relativistic) == String
        relativistic = parse(Bool, relativistic)
    end
    if relativistic
        mode = Relativistic()
    else
        mode = NoRelativistic()
    end
    disk_r_in = radiation_config[:disk_r_in]
    if disk_r_in == "isco"
        disk_r_in = bh.isco
    elseif disk_r_in == "r_in"
        disk_r_in = config[:initial_conditions][:r_in]
    end
    return QsosedRadiation(
        bh,
        radiation_config[:n_r],
        radiation_config[:f_x],
        disk_r_in,
        radiation_config[:z_xray],
        mode,
    )
end

function get_fuv_mdot(radiation::QsosedRadiation, r)
    r_index = searchsorted_nearest(radiation.disk_grid, r)
    f_uv = radiation.fuv_grid[r_index]
    mdot = radiation.mdot_grid[r_index]
    return f_uv, mdot
end

compute_radiation_constant(radiation::QsosedRadiation) = 3 / (π * radiation.efficiency)
