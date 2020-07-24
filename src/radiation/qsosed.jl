using RegionTrees
export ConstantFUV, QsosedRadiation, from_quadtree, get_f_uv
struct ConstantFUV <: Flag end
struct ConstantMdot <: Flag end

struct QsosedRadiation <: Radiation
    disk_grid::Vector{Float64}
    fuv_grid::Vector{Float64}
    mdot_grid::Vector{Float64}
    xray_luminosity::Float64
    function QsosedRadiation(rmin, rmax, nr, xray_luminosity)
        nothing
    end
end

from_quadtree(radiation::QsosedRadiation, quadtree::Cell) = QsosedRadiation(
    radiation.fuv,
    radiation.xray_luminosity,
    radiation.radiative_efficiency,
    quadtree,
    radiation.Rg,
)

getridx(radiation::QsosedRadiation, r) = searchsortednearest(radiation.disk_grid, r)
getfuv(radiation::QsosedRadiation, flag::ConstantFUV) = radiation.fuv_grid[1]
getmdot(radiation::QsosedRadiation, flag::ConstantMdot) = radiation.mdot_grid[1]
