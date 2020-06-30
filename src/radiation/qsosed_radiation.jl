using RegionTrees
export ConstantFUV, QsosedRadiation, from_quadtree, get_f_uv
struct ContantFUV <: Flag end

struct QsosedRadiation <: Radiation
    f_uv::Vector{Float64}
    disk_grid::Vector{Float64}
    xray_luminosity::Float64
    radiative_efficiency::Float64
    quadtree::Cell
    Rg::Float64
    function QsosedRadiation(bh::BlackHole, quadtree::Cell, f_uv, f_x, mu = 0.5)
        bolometric_luminosity = compute_bolometric_luminosity(bh)
        xray_luminosity = f_x * bolometric_luminosity
        new(f_uv, xray_luminosity, bh.Rg)
    end
end

from_quadtree(radiation::QsosedRadiation, quadtree::Cell) = QsosedRadiation(
    radiation.f_uv,
    radiation.xray_luminosity,
    radiation.radiative_efficiency,
    quadtree,
    radiation.Rg,
)

get_f_uv(radiation::QsosedRadiation, r) = searchsortednearest(radiation.f_uv, r)
get_f_uv(radiation::QsosedRadiation, flag::ContantFUV) = radiation.f_uv[1]
