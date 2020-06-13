using Qwind
using Test


@testset "Ionization Parameter" begin
    @test ionization_parameter(1, 1, 10, 0, 1) ≈ 1 / 20
    @test ionization_parameter(1, 1, 10, 100, 1) ≈ 0 atol = 1e-20
    @test ionization_parameter(1, 1, 10, 0, 2) ≈ 1 / 10
end

@testset "X-Ray Opacity" begin
    @test xray_opacity(1) == 100 * SIGMA_T
    @test xray_opacity(9e4) == 100 * SIGMA_T
    @test xray_opacity(1e5) == SIGMA_T
    @test xray_opacity(1e6) == SIGMA_T
end

@testset "X-Ray optical depth" begin
    @test xray_tau(3, 4, 2, 1e6) == 10 * SIGMA_T
    @test xray_tau(3, 4, 2, 1e-2) ≈ 1000 * SIGMA_T
    @test xray_tau(3, 4, 4, 1e6) == 20 * SIGMA_T
    @test xray_tau(12, 16, 2, 1e-2) ≈ 4000 * SIGMA_T
    @test xray_tau(4, 5, 2, 1e6, 1, 1) == 10 * SIGMA_T
    @test xray_tau([1, 3], [4 / 3, 4], [2, 2], [1e6, 1e6]) == 10 * SIGMA_T
end

@testset "Ionization radius" begin
    xray_lumin = 1e6
    density = 10.0
    r_in = 1.0
    xi_0 = 1e5
    r_x = Qwind.ionization_radius(xray_lumin, density, r_in, atol = 1e-2)
    @test Qwind.ionization_radius_kernel(xray_lumin, density, r_in, r_x, xi_0) ≈
          0 atol = 1e-2
    tau_x = xray_tau(r_x, 0, density, xi_0)
    @test ionization_parameter(r_x, 0, density, tau_x, xray_lumin) ≈ xi_0
    r_x = Qwind.ionization_radius(xray_lumin, density, r_in, atol = 1e-4)
    @test Qwind.ionization_radius_kernel(xray_lumin, density, r_in, r_x, xi_0) ≈
          0 atol = 1e-4
    tau_x = xray_tau(r_x, 0, density, xi_0)
    @test ionization_parameter(r_x, 0, density, tau_x, xray_lumin) ≈ xi_0
end

@testset "UV optical depth" begin
    @test uv_tau(3, 4, 2) == 10 * SIGMA_T
    @test uv_tau(3, 4, 4) == 20 * SIGMA_T
    @test uv_tau(12, 16, 2) ≈ 40 * SIGMA_T
    @test uv_tau(4, 5, 2, 1, 1) == 10 * SIGMA_T
end

@testset "Force multiplier" begin
    @test Qwind.force_multiplier_k(1e-3) == 0.411
    @test Qwind.force_multiplier_k(10^3.51) == 0.045
    @test Qwind.force_multiplier_k(1e-5) == 0.411
    @test Qwind.force_multiplier_k(1e10) == 0.013
    @test Qwind.force_multiplier_eta(1e-3) == 10^6.95
    @test Qwind.force_multiplier_eta(10^3.50) == 10^1.58
    @test Qwind.force_multiplier_eta(1e-5) == 10^6.95
    @test Qwind.force_multiplier_eta(1e10) == 10^0.78
    @test force_multiplier(1e5, 1e5) < 1e-3
    @test force_multiplier(1e5, 1e-5) < 1e-3
    @test force_multiplier(1e-8, 1e-5) > 2000
    @test force_multiplier(1e-8, 1e-5) < 3000
    @test force_multiplier(1e-10, 1e-5) ≈ 2300 atol = 0 rtol = 0.1
    @test force_multiplier(1e-2, 1e-5) ≈ 6 atol = 0 rtol = 0.1
end

#@testset "Force radiation" begin
#    @test all(isapprox.(
#        integrate(100.0, 50.0, wind, include_tau_uv = false, maxevals = 10000),
#        [2.4e-6, 2e-6],
#        atol = 0,
#        rtol = 0.1,
#    ))
#    @test all(isapprox.(
#        integrate(5.0, 10.0, wind, include_tau_uv = false, maxevals = 10000),
#        [-3.3e-6, 3e-5],
#        atol = 0,
#        rtol = 0.1,
#    ))
#    # @test all(isapprox.(force_radiation(1000., 0.1, 0., wind, include_tau_uv=false), [1.9e-11, 9e-11], atol=0, rtol=0.2))
#end
