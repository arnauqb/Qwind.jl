module Qwind
export Flag, FirstIter
abstract type Flag end
struct FirstIter <: Flag end


include("Tables.jl")
include("Grid.jl")
include("Thermodynamics.jl")
include("BlackHole.jl")
include("Constants.jl")
include("Streamline.jl")
#include("QuadTree.jl")
include("Radiation.jl")
include("Qsosed.jl")
include("Integration.jl")
include("RadiativeTransfer.jl")
include("Solver.jl")
include("Plotting.jl")


end # module
