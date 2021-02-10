using Qwind, Test, SafeTestsets
@time begin
#    @time @safetestset "Quadtree tests" begin include("radiation/quadtree_tests.jl") end
include("black_hole_tests.jl")

include("grid_tests.jl")

include("initial_conditions_tests.jl")

include("integrator_tests.jl")

include("radiation/base_tests.jl")

include("radiation/re_radiation_tests.jl")

include("thermodynamics_tests.jl")

include("utils_tests.jl")

include("utils.jl")
end
