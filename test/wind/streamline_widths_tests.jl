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

    r2 = [2., 2., 2.]
    z2 = [1., 2., 3.]
    sl2 = Streamline(0, t, r2, z2, vr, vphi, vz, n, width, escaped)

    r3 = [3., 3., 3.]
    z3 = [1., 2., 3.]
    sl3 = Streamline(0, t, r3, z3, vr, vphi, vz, n, width, escaped)

    @testset "Test KD-Tree" begin
        sl_kdtree_1 = KDTree(sl1)
        @test get_nearest_point(sl_kdtree_1, [0.0, 0.0]) == [1., 1.]
        @test get_nearest_point(sl_kdtree_1, [1.1, 2.1]) == [1., 2.]
        @test get_nearest_point(sl_kdtree_1, [3.0, 3.0]) == [1., 3.]
        sl_kdtree_2 = KDTree(sl2)
        @test get_nearest_point(sl_kdtree_2, [0.0, 0.0]) == [2., 1.]
        @test get_nearest_point(sl_kdtree_1, [2.1, 2.1]) == [2., 2.]
        @test get_nearest_point(sl_kdtree_2, [3.0, 3.0]) == [2., 3.]
    end

    @testset "Test distance to closest line" begin
    end

end
