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
end

