using Qwind, Test
@time begin
include("utils.jl")
include("black_hole_tests.jl")
include("parameters_tests.jl")

include("coordinate_interpolation_tests.jl")

include("initial_conditions/initial_conditions_tests.jl")
include("initial_conditions/cak_formalism_tests.jl")

include("integrator_tests.jl")

# Radiation
include("radiation/base_tests.jl")
include("radiation/xrays_tests.jl")
include("radiation/force_multiplier_tests.jl")
include("radiation/disk_tests.jl")
include("radiation/relativistic_corrections_tests.jl")
include("radiation/scattering_tests.jl")

# Wind 
include("wind/ray_tracing_tests.jl")
include("wind/optical_depth_tests.jl")
include("wind/wind_hull_tests.jl")
include("wind/density_grid_tests.jl")
include("wind/velocity_grid_tests.jl")
include("wind/base_tests.jl")

include("thermodynamics_tests.jl")
include("trajectory_intersection_tests.jl")

include("utils_tests.jl")

# saver tests
include("saver_tests.jl")
end
