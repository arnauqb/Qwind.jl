using Qwind
using RegionTrees
using StaticArrays
using Test
using QuadGK

function tau_integrand(s, density_function, r0, z0, r1, z1, Rg)
    r = r0 + s * (r1 - r0)
    z = z0 + s * (z1 - z0)
    return density_function(r, z) * SIGMA_T * Rg
end

function compute_analytic_tau(density_function, r0, z0, r1, z1, Rg)
    arclength = sqrt((r1 - r0)^2 + (z1 - z0)^2)
    res, err = quadgk(
        x -> tau_integrand(x, density_function, r0, z0, r1, z1, Rg),
        0.0,
        1.0,
        rtol = 1e-12,
    )
    return res[1] * arclength
end

@testset "Test X-Ray tau in tree leaf" begin
    # Used GeoGebra for this
    bh = BlackHole(1e8 * M_SUN, 0.5, 0)
    quadtree = Cell(SVector(1, 2), SVector(3, 4), 1e2)
    point = SVector(1.9, 2)
    intersection = SVector(2.65, 6)
    taux_0 = 0
    xray_luminosity = 1e60 # high luminosity σ_x = σ_t
    Rg = bh.Rg
    tau_x = compute_analytic_tau((r, z) -> 1e2, 1.9, 2, 2.65, 6, Rg)
    @test tau_x ≈ Qwind.compute_xray_tau_leaf(
        quadtree,
        quadtree,
        point,
        intersection,
        taux_0,
        xray_luminosity,
        Rg,
    )
    xray_luminosity = 10 # low luminosity σ_x = 100 σ_t
    @test tau_x * 1e2 ≈ Qwind.compute_xray_tau_leaf(
        quadtree,
        quadtree,
        point,
        intersection,
        taux_0,
        xray_luminosity,
        Rg,
    )
end

@testset "Test Quadtree refinement" begin
    function test_quadtree_get_density(density_func, rmax = 100, zmax = 100)
        quadtree, bh = create_test_quadtree(density_func)
        r_range = range(0, rmax, length = 50)
        z_range = range(0, zmax, length = 50)
        failed = 0
        failed_points = []
        failed_expected = []
        failed_values = []
        for r in r_range
            for z in z_range
                density = get_density(quadtree, r, z)
                expected = density_func(r, z)[1]
                if !isapprox(density, expected, rtol=0.1)
                    failed += 1
                    push!(failed_points, [r,z])
                    push!(failed_expected, expected)
                    push!(failed_values, density)
                end
            end
        end
        if failed > 0
            println("Failed at points")
            for (point, expected, value) in zip(failed_points, failed_expected, failed_values)
                println("point $point \t expected $expected \t $value")
            end
        end
        @test failed == 0
    end
    @testset "Density function $i" for i = 1:2 #1:length(density_functions)
        density_func = density_functions[i]
        println(density_func)
        test_quadtree_get_density(density_func)
    end
    #@testset "constant density" begin
    #    test_quadtree_get_density(constant_density)
    #end
    #@testset "linear density" begin
    #    test_quadtree_get_density(linear_density)
    #end
    #@testset "constant density" begin
    #    test_quadtree_get_density(constant_density)
    #end
end


#@testset "Test X-Ray tau" begin
#    r0 = 1.52
#    z0 = 0
#    r = 2.65
#    z = 6
#    xray_luminosity_high = 1e50
#    xray_luminosity_low = 1e1
#    function test_tau_x_equality(density_func)
#        quadtree, bh = create_test_quadtree(density_func)
#        tau_x = compute_xray_tau(quadtree, r, z, r0, z0, xray_luminosity_high, bh.Rg)
#        expected = compute_analytic_tau(density_func, 1.52, 0, 2.65, 6, bh.Rg)
#        @test tau_x ≈ expected rtol = 1e-2
#        tau_x = compute_xray_tau(quadtree, r, z, r0, z0, xray_luminosity_low, bh.Rg)
#        expected = compute_analytic_tau(density_func, 1.52, 0, 2.65, 6, bh.Rg) * 1e2
#        @test tau_x ≈ expected rtol = 1e-2
#    end
#    @testset "Constant density" begin
#        test_tau_x_equality((r,z) -> 1e8 * ones(length(r)))
#    end
#    @testset "Linear density" begin
#        test_tau_x_equality((r,z) -> 1e8 ./ (r+z))
#    end
##    @testset "Powerlaw density" begin
##        test_tau_x_equality((r,z) -> 1e8 ./ (r.^2 + z.^2))
##    end
##    @testset "Exponential density" begin
##        test_tau_x_equality((r,z) -> 1e8 .* exp.(-(r+z)))
##    end
#end

#@testset "Test UV tau" begin
#    r0_range = 10 .^ range(-3, 3, length = 10)
#    r_range = 10 .^ range(-3, 3, length = 10)
#    z_range = 10 .^ range(-3, 3, length = 10)
#    function test_tau_uv_equality(density_func, rtol = 1e-2)
#        quadtree, bh = create_test_quadtree(density_func)
#        error = false
#        for r0 in r0_range
#            for r in r_range
#                for z in z_range
#                    tau_uv = compute_uv_tau(quadtree, r0, r, z, bh.Rg)
#                    expected = compute_analytic_tau(density_func, r0, 0, r, z, bh.Rg)
#                    #@test tau_uv ≈ expected rtol = rtol
#                    if !(1 - abs(tau_uv / expected) < rtol)
#                        println("r0 $r0, r $r, z $z")
#                        println("tau_uv $tau_uv")
#                        println("expected $expected")
#                        println("Error!")
#                        error = true
#                        break
#                    end
#                end
#                error && break
#            end
#            error && break
#        end
#    end
#    #@testset "Constant density" begin
#    #    test_tau_uv_equality((r, z) -> 1e8 * ones(length(r)))
#    #end
#    @testset "Linear density" begin
#        test_tau_uv_equality((r, z) -> 1e8 * (r + z))
#    end
#    #@testset "Powerlaw density" begin
#    #    test_tau_uv_equality((r, z) -> 1e8 ./ (r .^ 2 + z .^ 2))
#    #end
#    #@testset "Exponential density" begin
#    #    test_tau_uv_equality((r, z) -> 1e8 .* exp.(-(r + z)), 2e-2)
#    #end
#end
