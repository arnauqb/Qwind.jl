using Qwind, Test
@time begin
include("utils.jl")
include("black_hole_tests.jl")

include("grid_tests.jl")

include("initial_conditions/initial_conditions_tests.jl")
include("initial_conditions/cak_formalism_tests.jl")

include("integrator_tests.jl")

include("radiation/base_tests.jl")
include("radiation/xrays_tests.jl")
include("radiation/force_multiplier_tests.jl")
include("radiation/disk_tests.jl")
include("radiation/wind_interpolator/ray_tracing_tests.jl")
include("radiation/wind_interpolator/optical_depth_tests.jl")

include("radiation/wind_interpolator/density_interpolators/interpolation_tests.jl")

include("thermodynamics_tests.jl")

include("utils_tests.jl")
end
