export QsosedRadiation, from_quadtree, get_fuv_mdot, relativistic_correction

struct QsosedRadiation <: Radiation
    disk_grid::Vector{Float64}
    fuv_grid::Vector{Float64}
    mdot_grid::Vector{Float64}
    xray_luminosity::Float64
    efficiency::Float64
    spin::Float64
    isco::Float64
    Rg::Float64
    flux_correction::FluxCorrection
end

flux_correction(::Relativistic, beta) = ((1 - beta) / (1 + beta))^2
flux_correction(::NoRelativistic, beta) = 1.0

function compute_mass_accretion_rate(radiation::Radiation, r)
    r_idx = searchsorted_nearest(radiation.disk_grid, r)
    mdot = radiation.mdot_grid[r_idx]
    return mdot * 4 * π * radiation.Rg * M_P * C / (SIGMA_T * radiation.efficiency)
end


function QsosedRadiation(bh::BlackHole, nr::Int, fx::Float64, relativistic::FluxCorrection)
    rmin = bh.isco
    rmax = 1400
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
        bh.Rg,
        relativistic,
    )
end
function QsosedRadiation(bh::BlackHole, config::Dict)
    radiation_config = config[:radiation]
    if radiation_config[:relativistic]
        mode = Relativistic()
    else
        mode = NoRelativistic()
    end
    return QsosedRadiation(bh, radiation_config[:n_r], radiation_config[:f_x], mode)
end

function get_fuv_mdot(radiation::QsosedRadiation, r)
    r_index = searchsorted_nearest(radiation.disk_grid, r)
    f_uv = radiation.fuv_grid[r_index]
    mdot = radiation.mdot_grid[r_index]
    return f_uv, mdot
end

compute_radiation_constant(radiation::QsosedRadiation) = 3 / (π * radiation.efficiency)
