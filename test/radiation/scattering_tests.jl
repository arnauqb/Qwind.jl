using Qwind, Test

@testset "test rectangle" begin
    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
    @test cell.rmin == 5
    @test cell.rmax == 10
    @test cell.zmin == 4
    @test cell.zmax == 11
end

@testset "compute intersection with cell" begin

    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)

    @testset "angles" begin
        @test Qwind.get_lowest_intercepting_angle(cell) ≈ 0.427 rtol = 1e-2
        @test Qwind.get_highest_intercepting_angle(cell) ≈ 1.190 rtol = 1e-2
    end

    @testset "intersections" begin
        theta = 0.4638
        @test Qwind.get_first_intersection(cell, theta) ≈ [5, 10] rtol = 1e-2
        @test Qwind.get_second_intersection(cell, theta) ≈ [5.5, 11] rtol = 1e-2
        @test Qwind.get_intersection_size(cell, theta) ≈ 1.118 rtol = 1e-2
        @test Qwind.get_distance_to_first_intersection(cell, theta) ≈ 11.18 rtol = 1e-2
    end
end

@testset "scattering optical depths" begin
    rr = range(0, 100, length = 100)
    zz = range(0, 100, length = 100)
    grid_values = 1e8 .* ones(length(rr), length(zz))
    dgrid = DensityGrid(rr, zz, grid_values)
end
