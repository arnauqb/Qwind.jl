using RegionTrees, PyCall
    QsosedRadiation, from_quadtree, get_fuv_mdot

struct QsosedRadiation <: Radiation
    disk_grid::Vector{Float64}
    fuv_grid::Vector{Float64}
    mdot_grid::Vector{Float64}
    xray_luminosity::Float64
    efficiency::Float64
    spin::Float64
    isco::Float64
    Rg::Float64
end

function QsosedRadiation(bh::BlackHole, nr::Int, fx::Float64)
    rmin = bh.isco
    rmax = 1400
    disk_grid = 10 .^ range(log10(rmin), log10(rmax), length=nr)
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
    )
end
function QsosedRadiation(bh::BlackHole, config::Dict)
    radiation_config = config["radiation"]
    return QsosedRadiation(bh, radiation_config["n_r"], radiation_config["f_x"])
end

function get_fuv_mdot(radiation::QsosedRadiation, r)
    r_index = searchsortednearest(radiation.disk_grid, r)
    f_uv = radiation.fuv_grid[r_index]
    mdot = radiation.mdot_grid[r_index]
    return f_uv, mdot
end
compute_radiation_constant(radiation::QsosedRadiation) = 6 / (8 * π * radiation.efficiency)
