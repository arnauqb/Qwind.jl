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
    @test Qwind.compute_cell_intersection(cell, (0, 1), (0, 1), (4, 1)) == [2, 1]
    @test Qwind.compute_cell_intersection(cell, (1, 0), (1, -1), (1, 4)) == [1, 2]
    cell = Cell(SVector(3, 4), SVector(2, 1), 10)
    @test Qwind.compute_cell_intersection(cell, (13 / 3, 4), (4, 0), (4.5, 6)) ≈
          [53 / 12, 5]
    @test Qwind.compute_cell_intersection(cell, (3.43, 4), (8, 0), (0, 7)) ≈ [3, 4.37625]
end

@testset "Quadtree functions" begin
    @test Qwind.create_quadtree(1, 2, 3, 4, vacuum_density = 1e3) ==
          Cell(SVector(1, 3), SVector(1, 1), 1e3)
end

constant_density(r, z) = 1e8 * ones(length(r))
linear_density(r, z) = 1e8 .* (r + 2 * z)
powerlaw_density_1(r, z) = 1e8 ./ (r + z)
powerlaw_density_2(r, z) = 1e8 .* (r + z)
exponential_density(r, z) = 1e8 .* exp.(-(r + z) / 20)
density_functions = [
    constant_density,
    linear_density,
    powerlaw_density_1,
    powerlaw_density_2,
    exponential_density,
]

function create_test_bh()
    bh = BlackHole(1e8 * M_SUN, 0.5, 0)
    bh
end

function create_test_kdtree(density_function; r_range, z_range, width_range)
    zmax = 100. * ones(length(r_range))
    if width_range === nothing
        width_range = 20 .* ones(length(zmax))
    end
    density = density_function(r_range, z_range)
    kdtree = create_wind_kdtree(r_range, z_range, zmax, width_range, density)
    kdtree
end

function create_test_quadtree(
    density_function;
    r_range,
    z_range,
    atol = 1e-3,
    rtol = 1e-3,
    cell_min_size = 1e-6,
    width_range = nothing,
)
    kdtree = create_test_kdtree(
        density_function,
        r_range = r_range,
        z_range = z_range,
        width_range = width_range,
    )
    bh = create_test_bh()
    quadtree = create_and_refine_quadtree(
        kdtree,
        Rg = bh.Rg,
        atol = atol,
        rtol = rtol,
        cell_min_size = cell_min_size,
    )
    return quadtree
end

@testset "Test Quadtree refinement" begin
    function test_quadtree_get_density(density_func)
        r_range_quadtree = 10 .* rand(1000000)
        z_range_quadtree = 10 .* rand(1000000)
        quadtree = create_test_quadtree(
            density_func,
            r_range = r_range_quadtree,
            z_range = z_range_quadtree,
        )
        r_range = range(0, 10, length=100)
        z_range = range(0, 10, length=100)
        failed = 0
        failed_points = []
        failed_expected = []
        failed_values = []
        for r in r_range
            for z in z_range
                density = get_density(quadtree, r, z)
                expected = density_func(r, z)[1]
                if expected == Inf
                    continue
                end
                if !isapprox(density, expected, rtol = 1)
                    failed += 1
                    push!(failed_points, [r, z])
                    push!(failed_expected, expected)
                    push!(failed_values, density)
                end
            end
        end
        if failed > 0
            println("Failed at points")
            for (point, expected, value) in
                zip(failed_points, failed_expected, failed_values)
                println("point $point \t expected $expected \t $value")
            end
        end
        @test failed == 0
    end
    @testset "Density function $i" for i = 1:length(density_functions)
        density_func = density_functions[i]
        println(density_func)
        test_quadtree_get_density(density_func)
    end
end
