module Qwind
export Radiation, RadiativeTransfer, InitialConditions, Model, scipy_interpolate


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
include("radiative_transfer/base.jl")
include("radiative_transfer/regular_grid.jl")
# risaliti elvis
include("radiative_transfer/risaliti_elvis.jl")
# wind interpolation
include("radiative_transfer/wind_interpolator/base.jl")
include("radiative_transfer/wind_interpolator/grids.jl")
include("radiative_transfer/wind_interpolator/integrators_interpolator.jl")
include("radiative_transfer/wind_interpolator/wind_hull.jl")
include("radiative_transfer/wind_interpolator/density_grid.jl")
include("radiative_transfer/wind_interpolator/velocity_grid.jl")
include("radiative_transfer/wind_interpolator/ray_tracing.jl")

# initial conditions
include("initial_conditions/initial_conditions.jl")

# integrator
include("integrator.jl")
include("intersections.jl")
#include("plotting.jl")
include("initial_conditions/cak_formalism.jl")

# saver
include("saver.jl")

# model
include("model.jl")


# initialisation
const scipy_interpolate = PyNULL()
function __init__()
    copy!(scipy_interpolate, pyimport_conda("scipy.interpolate", "scipy"))
    setup_logging()
end

end # module
