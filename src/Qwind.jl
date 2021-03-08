module Qwind
export Radiation, RadiativeTransfer, InitialConditions, Model

include("black_hole.jl")
include("types.jl")
include("utils.jl")
include("cluster_utils.jl")
include("grid.jl")

include("constants.jl")
include("tables.jl")
include("thermodynamics.jl")
# radiation
include("radiation/base.jl")
include("radiation/disk.jl")
include("radiation/qsosed.jl")
include("radiation/risaliti_elvis.jl")
# radiative transfer
include("radiative_transfer/density_interpolators/grids.jl")
include("radiative_transfer/base.jl")
include("radiative_transfer/risaliti_elvis.jl")
include("radiative_transfer/regular_grid.jl")
include("radiative_transfer/density_interpolators/base.jl")
include("radiative_transfer/ray_tracing.jl")
include("radiative_transfer/density_interpolators/interpolation.jl")

# initial conditions
include("initial_conditions.jl")

# integrator
include("integrator.jl")
#include("plotting.jl")

# saver
include("saver.jl")

# model
include("model.jl")

end # module
