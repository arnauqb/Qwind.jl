module Qwind
export Radiation, RadiativeTransfer, InitialConditions, Model

include("black_hole.jl")
include("types.jl")
include("utils.jl")
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
include("radiative_transfer/base.jl")
include("radiative_transfer/risaliti_elvis.jl")
include("radiative_transfer/ray_tracing.jl")
include("radiative_transfer/regular_grid.jl")
include("radiative_transfer/density_interpolators/base.jl")
include("radiative_transfer/density_interpolators/nntree.jl")
#include("radiative_transfer/density_interpolators/quadtree.jl")
include("radiative_transfer/density_interpolators/grids.jl")
#include("radiative_transfer/density_interpolators/smoothed_grid.jl")
#include("radiative_transfer/density_interpolators/log_grid.jl")
include("radiative_transfer/density_interpolators/z_interpolation.jl")
#include("radiative_transfer/adaptive_mesh.jl")
include("initial_conditions.jl")

# integrator
include("integrator.jl")
#include("plotting.jl")

# saver
include("saver.jl")

# model
include("model.jl")

end # module
