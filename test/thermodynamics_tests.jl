using Qwind
using Test

@testset "thermal velocity" begin
    @test compute_thermal_velocity(2e6) ≈ 12848657.328083131 # cm / s
    @test compute_thermal_velocity(25e3) ≈ 1436523.560259735# cm / s
end
