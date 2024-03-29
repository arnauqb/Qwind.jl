using Qwind
using Test

@testset "Test cell intersections" begin

    r_range = range(0, 10, step = 0.5)
    z_range = range(0, 10, step = 0.5)

    function check_equal(point, intersection; rtol = 1e-3)
        @test isapprox(point[1], intersection[1], rtol = rtol) &&
              isapprox(point[2], intersection[3], rtol = rtol)
    end

    @testset "Normal path" begin
        gi = GridIterator(r_range, z_range, 1.5, 0.0, 2.5, 1.5)
        y(x) = -2.25 + 1.5 * x
        x(y) = (2.25 + y) / 1.5
        expected_points = [[1.5, 0.0, 0.0], [x(0.5), 0.5], [2, y(2)], [x(1), 1], [2.5, 1.5]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Normal path backwards" begin
        gi = GridIterator(r_range, z_range, 2.5, 1.5, 1.5, 0.0)
        y(x) = -2.25 + 1.5 * x
        x(y) = (2.25 + y) / 1.5
        expected_points =
            [[2.5, 1.5], [x(1.0), 1.0], [2.0, y(2.0)], [x(0.5), 0.5], [1.5, 0.0]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        next_intersection!(gi)
        @test gi.finished
    end

    @testset "Purely vertical path" begin
        gi = GridIterator(r_range, z_range, 1.5, 0.0, 1.5, 2.0)
        expected_points = [[1.5, 0.0], [1.5, 0.5], [1.5, 1.0], [1.5, 1.5], [1.5, 2.0]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Purely vertical path backwards" begin
        gi = GridIterator(r_range, z_range, 1.5, 2.0, 1.5, 0.0)
        expected_points =
            reverse([[1.5, 0.0], [1.5, 0.5], [1.5, 1.0], [1.5, 1.5], [1.5, 2.0]])
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
    end

    @testset "Purely horizontal path" begin
        gi = GridIterator(r_range, z_range, 0.1, 2.0, 1.7, 2.0)
        expected_points = [[0.1, 2.0], [0.5, 2.0], [1.0, 2.0], [1.5, 2.0], [1.7, 2.0]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Purely horizontal path backwards" begin
        gi = GridIterator(r_range, z_range, 1.7, 2.0, 0.1, 2.0)
        expected_points =
            reverse([[0.1, 2.0], [0.5, 2.0], [1.0, 2.0], [1.5, 2.0], [1.7, 2.0]])
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Completely outside path" begin
        gi = GridIterator(r_range, z_range, -10, -20, -50, -60)
        @test gi.finished == true
    end

    @testset "Outside path falls short" begin
        gi = GridIterator(r_range, z_range, -10, -20, -5, -6)
        expected_points = [[-10, -20], [-5, -6]]
        d = 0
        while !gi.finished
            int = copy(gi.intersection)
            next_intersection!(gi)
            int2 = copy(gi.intersection)
            d += d_euclidean(int[1], int2[1], int[2], int2[2])
        end
        @test d == 0
        @test gi.finished == true
    end

    @testset "Paths starts outside grid bottom left" begin
        gi = GridIterator(r_range, z_range, -1.0, -0.5, 0.5, 0.5)
        y(x) = (0.25 + x) / 1.5
        x(y) = 1.5y - 0.25
        expected_points = [[0, y(0)], [0.5, 0.5]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths starts outside grid bottom right" begin
        gi = GridIterator(r_range, z_range, 10.69, -1.19, 8, 1)
        y(x) = (20.21 - 2.19x) / 2.69
        x(y) = (20.21 - 2.69y) / 2.19
        expected_points = [[x(0), 0.0], [9, y(9)], [x(0.5), 0.5], [8.5, y(8.5)], [8, 1]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection, rtol=5e-2)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths starts outside grid top right" begin
        gi = GridIterator(r_range, z_range, 11.99, 11.36, 9.46, 9.32)
        x(y) = (-4.26 + 2.53y) / 2.04
        y(x) = (-4.26 - 2.04x) / (-2.53)
        expected_points = [[10, y(10)], [x(9.5), 9.5], [9.5, y(9.5)], [9.46, 9.32]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection, rtol=5e-2)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths starts outside grid top left" begin
        gi = GridIterator(r_range, z_range, -1.63, 10.96, 3.2, 8.84)
        x(y) = (49.49 - 4.83y) / 2.12
        y(x) = (49.49 - 2.12x) / 4.83
        expected_points = [
            [0.5, y(0.5)],
            [1, y(1)],
            [1.5, y(1.5)],
            [x(9.5), 9.5],
            [2, y(2)],
            [2.5, y(2.5)],
            [x(9), 9],
            [3, y(3)],
            [3.2, 8.84],
        ]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection, rtol=2e-1)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths ends outside grid bottom left" begin
        gi = GridIterator(r_range, z_range, 0.5, 0.5, -1.0, -0.5)
        y(x) = (0.25 + x) / 1.5
        x(y) = 1.5y - 0.25
        expected_points = [[0.5, 0.5], [0, y(0)]]
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths end outside grid bottom right" begin
        gi = GridIterator(r_range, z_range, 8, 1, 10.69, -1.19)
        y(x) = (20.19 - 2.19x) / 2.69
        x(y) = (20.19 - 2.69y) / 2.19
        expected_points =
            reverse([[x(0), 0.0], [9, y(9)], [x(0.5), 0.5], [8.5, y(8.5)], [8, 1]])
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection, rtol=5e-2)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths ends outside grid top right" begin
        gi = GridIterator(r_range, z_range, 9.46, 9.32, 11.99, 11.36)
        x(y) = (-4.26 + 2.53y) / 2.04
        y(x) = (-4.26 - 2.04x) / (-2.53)
        expected_points = reverse([[10, y(10)], [x(9.5), 9.5], [9.5, y(9.5)], [9.46, 9.32]])
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection, rtol=1e-2)
            next_intersection!(gi)
        end
        @test gi.finished
    end

    @testset "Paths ends outside grid top left" begin
        gi = GridIterator(r_range, z_range, 3.2, 8.84, -1.63, 10.96)
        x(y) = (49.49 - 4.83y) / 2.12
        y(x) = (49.49 - 2.12x) / 4.83
        expected_points = reverse([
            [0.5, y(0.5)],
            [1, y(1)],
            [1.5, y(1.5)],
            [x(9.5), 9.5],
            [2, y(2)],
            [2.5, y(2.5)],
            [x(9), 9],
            [3, y(3)],
            [3.2, 8.84],
        ])
        for (i, point) in enumerate(expected_points)
            check_equal(point, gi.intersection, rtol=5e-1)
            next_intersection!(gi)
        end
        @test gi.finished
    end

end
