using RegionTrees, PyCall
export ConstantFUV,
    QsosedRadiation, from_quadtree, compute_disk_temperature, get_fuv_mdot
struct ConstantFUV <: Flag end
struct ConstantMdot <: Flag end

struct QsosedRadiation <: Radiation
    disk_grid::Vector{Float64}
    fuv_grid::Vector{Float64}
    mdot_grid::Vector{Float64}
    xray_luminosity::Float64
    sed::PyObject
    efficiency::Float64
    spin::Float64
    isco::Float64
    Rg::Float64
end

function QsosedRadiation(bh::BlackHole, nr::Int, fx::Float64)
    qsosed_python = pyimport("qsosed")
    sed = qsosed_python.SED(M = bh.M / M_SUN, mdot = bh.mdot, number_bins_fractions = nr)
    rmin = bh.isco
    rmax = 1400
    disk_grid, uv_fractions = sed.compute_uv_fractions(
        inner_radius = rmin,
        outer_radius = rmax,
        return_all = false,
        log_spaced = true,
    )
    mdot_grid = bh.mdot .* ones(length(disk_grid))
    xray_luminosity = fx * compute_bolometric_luminosity(bh)
    return QsosedRadiation(
        disk_grid,
        uv_fractions,
        mdot_grid,
        xray_luminosity,
        sed,
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

compute_disk_temperature(radiation::QsosedRadiation, r) =
    radiation.sed.disk_nt_temperature4(r)^(1.0 / 4.0)
