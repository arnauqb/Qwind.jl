using Test
using Qwind

@testset "Test grid properties" begin
    grid = Grid(0.0, 2000.0, 0.0, 1000.0, 1e2)
    @test grid.r_min == 0.0
    @test grid.z_min == 0.0
    @test grid.r_max == 2000.0
    @test grid.z_max == 1000.0
    @test grid.vacuum_number_density == 1e2
    line = Streamline(0, 3000, 1, 1, 1, 1, 1)
    @test out_of_grid(grid, line) == true
end
