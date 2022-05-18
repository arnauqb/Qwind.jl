using Qwind, Test

@testset "Test line widths calculations" begin
    r1 = [1., 1., 1.]
    z1 = [1., 2., 3.]
    t = zero(r1)
    vr = zero(r1)
    vz = zero(r1)
    vphi = zero(r1)
    n = zero(r1)
    width = zero(r1)
    escaped = false
    sl1 = Streamline(0, t, r1, z1, vr, vphi, vz, n, width, escaped)
    sl_kdtree_1 = KDTree(sl1)

    r2 = [4., 2., 2.]
    z2 = [1., 2., 3.]
    sl2 = Streamline(0, t, r2, z2, vr, vphi, vz, n, width, escaped)
    sl_kdtree_2 = KDTree(sl2)

    r3 = [3., 3., 3.]
    z3 = [1., 2., 3.]
    sl3 = Streamline(0, t, r3, z3, vr, vphi, vz, n, width, escaped)
    sl_kdtree_3 = KDTree(sl3)

    @testset "Test KD-Tree" begin
        @test Qwind.get_nearest_point(sl_kdtree_1, [0.0, 0.0]) == [1., 1.]
        @test Qwind.get_nearest_point(sl_kdtree_1, [1.1, 2.1]) == [1., 2.]
        @test Qwind.get_nearest_point(sl_kdtree_1, [3.0, 3.0]) == [1., 3.]
        @test Qwind.get_nearest_point(sl_kdtree_2, [0.0, 0.0]) == [2., 2.]
        @test Qwind.get_nearest_point(sl_kdtree_2, [2.1, 2.1]) == [2., 2.]
        @test Qwind.get_nearest_point(sl_kdtree_2, [3.0, 3.0]) == [2., 3.]
    end

    @testset "Test distance to line" begin
        @test Qwind.get_distance_to_line(sl_kdtree_1, [0.0, 0.0]) == sqrt(2)
        @test Qwind.get_distance_to_line(sl_kdtree_2, [0.0, 0.0]) == sqrt(8)
        @test Qwind.get_distance_to_line(sl_kdtree_3, [0.0, 0.0]) == sqrt(10)
    end

    @testset "Test distance between lines" begin
        line_widths = Qwind.get_distances_between_lines(sl2, sl1)
        lw1 = 3
        lw2 = 1
        lw3 = 1
        @test line_widths == [lw1, lw2, lw3]
    end

    @testset "Test get all widths" begin
        streamlines = Streamlines([sl1, sl2, sl3])
        widths = get_streamlines_widths(streamlines)
        @test length(widths) == 2
        @test widths[1] == [sqrt(2), 1, 1]
        @test widths[2] == [1, 1, 1]
    end

end
