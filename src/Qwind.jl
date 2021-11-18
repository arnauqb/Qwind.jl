module Qwind
using TimerOutputs

const timer = TimerOutput()
disable_timer!(timer)
export timer

include("types.jl")
include("constants.jl")
include("parameters.jl")
include("black_hole.jl")
include("distances.jl")
include("utils.jl")
include("cluster_utils.jl")

include("tables.jl")
include("thermodynamics.jl")
#
#
# Streamline struct
include("streamlines.jl")

# wind
include("wind/grids.jl")
include("wind/wind_hull.jl")
include("wind/density_grid.jl")
include("wind/velocity_grid.jl")
include("wind/iterator.jl")
include("wind/base.jl")

# radiation
include("radiation/disk.jl")
include("radiation/force_multiplier.jl")
include("radiation/scattering.jl")
include("radiation/base.jl")
include("radiation/xrays.jl")
include("radiation/relativistic_correction.jl")
include("radiation/disk_integration.jl")
include("radiation/optical_depths.jl")

# radiation

# initial conditions
include("initial_conditions/initial_conditions.jl")

# integrator
include("integrator.jl")
include("trajectory_intersections.jl")
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
