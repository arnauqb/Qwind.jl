module Qwind

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
include("radiation/distances.jl")
include("radiation/disk.jl")
include("radiation/force_multiplier.jl")
include("radiation/xrays.jl")
include("radiation/relativistic_correction.jl")
include("radiation/disk_integration.jl")

# wind interpolator
include("radiation/wind_interpolator/base.jl")
include("radiation/wind_interpolator/grids.jl")
include("radiation/wind_interpolator/integrators_interpolator.jl")
include("radiation/wind_interpolator/wind_hull.jl")
include("radiation/wind_interpolator/density_grid.jl")
include("radiation/wind_interpolator/velocity_grid.jl")
include("radiation/wind_interpolator/ray_tracing.jl")

# initial conditions
include("initial_conditions/initial_conditions.jl")

# integrator
include("integrator.jl")
include("intersections.jl")
include("initial_conditions/cak_formalism.jl")

# saver
include("saver.jl")

# model
include("model.jl")


# PyCall initialisations
const scipy_interpolate = PyNULL()
function __init__()
    copy!(scipy_interpolate, pyimport_conda("scipy.interpolate", "scipy"))
    setup_logging()
end

end # module
