using Qwind
using StaticArrays
using Test
using RegionTrees

@testset "Quadtree cell basic functions" begin
    cell = Cell(SVector(1, 2), SVector(3, 4), 10)
    @test getwidth(cell) == 3
    @test getheight(cell) == 4
    @test getcellsize(cell) == 5
    @test getrmin(cell) == 1
    @test getzmin(cell) == 2
    @test getrmax(cell) == 4
    @test getzmax(cell) == 6
    @test getlimits(cell) == [1, 4, 2, 6]
    cell = Cell(SVector(1, 2), SVector(3, 4), 10)
    @test compute_cell_optical_thickness(cell, 100) ≈ 5 * 10 * 100 * SIGMA_T
end

@testset "Ray intersection with cell" begin
    # all answers derived with geogebra
    cell = Cell(SVector(0, 0), SVector(2, 2), 10)
    @test Qwind.compute_cell_intersection(cell, (0, 0), (0, 0), (4, 4)) == [2, 2]
    @test Qwind.compute_cell_intersection(cell, (0, 0), (0, 0), (6, 1)) == [2, 1 / 3]
    @test Qwind.compute_cell_intersection(cell, (0, 0), (0, 0), (1, 6)) == [1 / 3, 2]

    # extreme cases (completely flat or vertical lines)
    @test Qwind.compute_cell_intersection(cell, (0, 1), (0, 1), (4, 1)) == [2,1]
    @test Qwind.compute_cell_intersection(cell, (1, 0), (1, -1), (1, 4)) == [1,2]

    cell = Cell(SVector(3, 4), SVector(2, 1), 10)
    @test Qwind.compute_cell_intersection(cell, (13 / 3, 4), (4, 0), (4.5, 6)) ≈ [53 / 12, 5]
    @test Qwind.compute_cell_intersection(cell, (3.43, 4), (8, 0), (0, 7)) ≈ [3, 4.37625]
end

@testset "Quadtree functions" begin
    @test Qwind.create_quadtree(1, 2, 3, 4, vacuum_density=1e3) == Cell(SVector(1, 3), SVector(1, 1), 1e3)
end