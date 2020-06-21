module Qwind
export Flag, FirstIter, Radiation, RadiationTransfer, InitialConditions, Model
abstract type Flag end
struct FirstIter <: Flag end

abstract type Radiation end
abstract type RadiationTransfer end
abstract type InitialConditions end

include("utils.jl")
include("black_hole.jl")
include("grid.jl")

struct Model
    grid::Grid
    bh::BlackHole
    rad::Radiation
    rt::RadiationTransfer
    ic::InitialConditions
end

include("constants.jl")
include("tables.jl")
include("thermodynamics.jl")
include("radiation/base.jl")
include("radiation/re_radiation.jl")
include("radiation/qsosed_radiation.jl")

#include("Streamline.jl")
##include("QuadTree.jl")
#include("Radiation.jl")
#include("Qsosed.jl")
#include("Integration.jl")
#include("RadiativeTransfer.jl")
#include("Solver.jl")
#include("Plotting.jl")
#

end # module
