module Qwind
export Flag
abstract type Flag end


include("Tables.jl")
include("Grid.jl")
include("Thermodynamics.jl")
include("BlackHole.jl")
include("Constants.jl")
include("Streamline.jl")
include("Radiation.jl")
include("Qsosed.jl")
include("Integration.jl")
include("Solver.jl")


end # module
