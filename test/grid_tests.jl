using Test
using Qwind

@testset "Test grid properties" begin
    grid = Grid(0.0, 2000.0, 0.0, 1000.0)
    @test grid.r_min == 0.0
    @test grid.z_min == 0.0
    @test grid.r_max == 2000.0
    @test grid.z_max == 1000.0
end
