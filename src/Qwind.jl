module Qwind
export Radiation, RadiativeTransfer, InitialConditions, Model

include("utils.jl")
include("black_hole.jl")
include("grid.jl")

include("constants.jl")
include("tables.jl")
include("thermodynamics.jl")
# radiation
include("radiation/base.jl")
include("radiation/disk.jl")
include("radiation/qsosed.jl")
include("radiation/risaliti_elvis.jl")
# radiation transfer
include("radiative_transfer/risaliti_elvis.jl")
include("radiative_transfer/nntree.jl")
include("radiative_transfer/quadtree.jl")
include("radiative_transfer/adaptive_mesh.jl")
include("initial_conditions.jl")

# integrator
include("integrator.jl")
include("plotting.jl")

end # module
