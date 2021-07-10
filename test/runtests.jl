using Qwind, Test
@time begin
#    @time @safetestset "Quadtree tests" begin include("radiation/quadtree_tests.jl") end
include("utils.jl")
include("black_hole_tests.jl")

include("grid_tests.jl")

include("initial_conditions/initial_conditions_tests.jl")
include("initial_conditions/cak_formalism_tests.jl")

include("integrator_tests.jl")

include("radiation/base_tests.jl")
include("radiation/disk_tests.jl")
#include("radiation/nn_tree_tests.jl")
include("radiative_transfer/wind_interpolator/ray_tracing_tests.jl")
include("radiative_transfer/wind_interpolator/optical_depth_tests.jl")
include("radiative_transfer/re_radiation_tests.jl")

#include("radiative_transfer/density_interpolators/z_interpolation_tests.jl")
include("radiative_transfer/density_interpolators/interpolation_tests.jl")

include("thermodynamics_tests.jl")

include("utils_tests.jl")
end
