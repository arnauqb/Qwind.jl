module Qwind
export RE, QSOSED, RadiationMode, TestMode, NoRad, ConstantRad
abstract type Flag end
abstract type RadiationMode <: Flag end
struct RE <: RadiationMode end
struct QSOSED <: RadiationMode end
struct NoRad <: RadiationMode end
struct TestMode <: RadiationMode end
struct ConstantRad <: RadiationMode end

include("Tables.jl")
include("BlackHole.jl")
include("Constants.jl")
include("Qsosed.jl")
include("Streamline.jl")
include("Solver.jl")
include("Radiation.jl")
include("Integration.jl")
include("Grid.jl")


end # module
