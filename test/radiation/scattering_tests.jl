using Qwind, Test

@testset "compute intersection with cell" begin
    cell = Rectangle(rmin=5, rmax=10, zmin=4, zmax=11)
    @test get_lowest_intercepting_angle(cell) ≈ 0.38 rtol = 1e-2
    @test get_highest_intercepting_angle(cell) ≈ 1.144 rtol = 1e-2

    theta = 1.107
    @test get_first_intersection(cell, theta) ≈ [5,10]
    @test get_second_intersection(cell, theta) ≈ [5.5, 11]
    @test get_intersection_size(cell, theta) ≈ 1.118 rtol = 1e-2
    @test get_distance_to_first_intersection(cell, theta) ≈ 11.18 rtol 1e-2
end

@testset "scattering optical depths" begin
    rr = range(0, 100, length=100);
    zz = range(0, 100, length=100);
    grid_values = 1e8 .* ones(length(rr), length(zz));
    dgrid = DensityGrid(rr, zz, grid_values);
end
