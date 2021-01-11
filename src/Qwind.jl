module Qwind
export Flag, FirstIter, Radiation, RadiativeTransfer, InitialConditions, Model
abstract type Flag end
struct FirstIter <: Flag end

abstract type Radiation end
abstract type RadiativeTransfer end
abstract type InitialConditions end

include("utils.jl")
include("black_hole.jl")
include("grid.jl")

struct Model
    grid::Grid
    bh::BlackHole
    rad::Radiation
    rt::RadiativeTransfer
    ic::InitialConditions
end

include("constants.jl")
include("tables.jl")
include("thermodynamics.jl")
# radiation
include("radiation/base.jl")
include("radiation/disk.jl")
include("radiation/qsosed.jl")
include("radiation/risaliti_elvis.jl")
# radiation transfer
include("radiation_transfer/risaliti_elvis.jl")
include("radiation_transfer/nntree.jl")
include("radiation_transfer/quadtree.jl")
include("radiation_transfer/adaptive_mesh.jl")
include("initial_conditions.jl")

# integrator
include("integrator.jl")
include("plotting.jl")

end # module
