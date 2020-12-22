using Qwind
using Test

#@testset "Dense solution from integrators" begin
# I don't know how to test this in an easy way
#
#
#end
@testset "Wind KDTree" begin
    r = [1, 5, 10]
    z = [15, 20, 25]
    zmax = [25, 25, 25]
    density = [10, 20, 30]
    width = [6, 2, 1]
    kdtree = create_wind_kdtree(r, z, zmax, width, density)
    @test kdtree.r == r
    @test kdtree.z == z
    @test kdtree.zmax == zmax
    @test kdtree.width == width
    @test kdtree.n == density
    @testset "KDTree -- get closest point" begin
        rc, zc, zmaxc, widthc, nc = get_closest_points(kdtree, 8, 18, 2) 
        @test Set(rc) == Set([5, 10])
        @test Set(zc) == Set([20, 25])
        @test Set(zmaxc) == Set([25])
        @test Set(widthc) == Set([2,1])
        @test Set(nc) == Set([20, 30])
    end
    @testset "KDTree -- get density of closest point" begin
        @test get_density(kdtree, 2, 16) == 15
        @test get_density(kdtree, 15, 22) == 1e2
        @test get_density(kdtree, 5, 26) == 1e2
    end
end
