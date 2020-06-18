using Qwind
using Test


@testset "Ionization Parameter" begin
    @test compute_ionization_parameter(1, 1, 10, 0, 1, 1) ≈ 1 / 20
    @test compute_ionization_parameter(1, 1, 10, 100, 1, 1) ≈ 0 atol = 1e-20
    @test compute_ionization_parameter(1, 1, 10, 0, 2, 2) ≈ 1 / 10 / 4
end

@testset "X-Ray Opacity" begin
    @test compute_xray_opacity(1) == 100 * SIGMA_T
    @test compute_xray_opacity(9e4) == 100 * SIGMA_T
    @test compute_xray_opacity(1e5) == SIGMA_T
    @test compute_xray_opacity(1e6) == SIGMA_T
end

@testset "X-Ray optical depth" begin
    @test compute_xray_tau(3, 4, 2, 1e6, 1) == 10 * SIGMA_T
    @test compute_xray_tau(3, 4, 2, 1e-2, 2) ≈ 2000 * SIGMA_T
    @test compute_xray_tau(3, 4, 4, 1e6, 2) == 40 * SIGMA_T
    @test compute_xray_tau(12, 16, 2, 1e-2, 1) ≈ 4000 * SIGMA_T
    @test compute_xray_tau(4, 5, 2, 1e6, 1, 1, 1) == 10 * SIGMA_T
end

@testset "RE X-ray optical depth" begin
    # r<= r_in
    Rg = 2
    @test compute_xray_tau(4, 3, 10, 20, 5, 100, 2, Rg) / SIGMA_T ==
          0.0
    # r < r_0, r_x < r
    @test compute_xray_tau(4, 3, 10, 20, 1, 2, 6, Rg) / SIGMA_T ≈
          4010 / cos(4 / 5) atol = 0 rtol = 0.1
    # r < r_0, r_x > r
    @test compute_xray_tau(4, 3, 10, 20, 1, 5, 6, Rg) / SIGMA_T ≈
          60 / cos(4 / 5) atol = 0 rtol = 0.1
    # r > r_0, r_in < r_x < r_0
    @test compute_xray_tau(4, 3, 10, 20, 1, 2, 3, Rg) / SIGMA_T ≈
          6010 / cos(4 / 5) atol = 0 rtol = 0.1
    # r > r_0, r_in < r_x < r_0
    @test compute_xray_tau(4, 3, 10, 20, 1, 3.5, 3, Rg) / SIGMA_T ≈
          2030 / cos(4 / 5) atol = 0 rtol = 0.1
    # r > r_0, r_in < r_0 < r < r
    @test compute_xray_tau(4, 3, 10, 20, 1, 5, 2, Rg) / SIGMA_T ≈
          100 / cos(4 / 5) atol = 0 rtol = 0.1
    # r > r_0, r_in < r_0 < r_x < r
    @test compute_xray_tau(4, 3, 10, 20, 1, 3, 2, Rg) / SIGMA_T ≈
          4030 / cos(4 / 5) atol = 0 rtol = 0.1

end

@testset "Ionization radius" begin
    xray_lumin = 1e6
    density = 10.0
    r_in = 1.0
    xi_0 = 1e5
    r_x =
        Qwind.compute_ionization_radius(xray_lumin, density, r_in, atol = 1e-2)
    @test Qwind.ionization_radius_kernel(xray_lumin, density, r_in, r_x, xi_0) ≈
          0 atol = 1e-2
    tau_x = compute_xray_tau(r_x, 0, density, xi_0, 1)
    @test compute_ionization_parameter(r_x, 0, density, tau_x, xray_lumin, 1) ≈
          xi_0
    r_x =
        Qwind.compute_ionization_radius(xray_lumin, density, r_in, atol = 1e-4)
    @test Qwind.ionization_radius_kernel(xray_lumin, density, r_in, r_x, xi_0) ≈
          0 atol = 1e-4
    tau_x = compute_xray_tau(r_x, 0, density, xi_0, 1)
    @test compute_ionization_parameter(r_x, 0, density, tau_x, xray_lumin, 1) ≈
          xi_0
end

@testset "UV optical depth" begin
    @test compute_uv_tau(3, 4, 2) == 10 * SIGMA_T
    @test compute_uv_tau(3, 4, 4) == 20 * SIGMA_T
    @test compute_uv_tau(12, 16, 2) ≈ 40 * SIGMA_T
    @test compute_uv_tau(4, 5, 2, 1, 1) == 10 * SIGMA_T
end

@testset "Force multiplier" begin
    @test Qwind.compute_force_multiplier_k(1e-3) == 0.411
    @test Qwind.compute_force_multiplier_k(1e-3, NoInterp()) ≈ 0.411 rtol=0.1
    @test Qwind.compute_force_multiplier_k(10^3.51) == 0.045
    @test Qwind.compute_force_multiplier_k(1e-5) == 0.411
    @test Qwind.compute_force_multiplier_k(1e10) == 0.013
    @test Qwind.compute_force_multiplier_k(1e10, NoInterp()) == 0.03
    @test Qwind.compute_force_multiplier_eta(1e-3) == 10^6.95
    @test Qwind.compute_force_multiplier_eta(1e-3, NoInterp()) ≈ 9.33e6 rtol=0.1
    @test Qwind.compute_force_multiplier_eta(10^3.50) == 10^1.58
    @test Qwind.compute_force_multiplier_eta(1e-5) == 10^6.95
    @test Qwind.compute_force_multiplier_eta(1e10) == 10^0.78
    @test Qwind.compute_force_multiplier_eta(1e10, NoInterp()) == 1.0
    @test compute_force_multiplier(1e5, 1e5) < 1e-3
    @test compute_force_multiplier(1e5, 1e-5) < 1e-3
    @test compute_force_multiplier(1e-8, 1e-5) > 2000
    @test compute_force_multiplier(1e-8, 1e-5) < 3000
    @test compute_force_multiplier(1e-10, 1e-5) ≈ 2300 atol = 0 rtol = 0.1
    @test compute_force_multiplier(1e-2, 1e-5) ≈ 6 atol = 0 rtol = 0.1
    @test compute_tau_eff(1e8, 0, 1e6) == 1.0
end

@testset "Radiation structure" begin
      black_hole = BlackHole(1e8 * M_SUN, 0.5, 0.0)
      radiation = SimpleRadiation(black_hole, 0.70906799789733695, 0.14436810902317002, 1e8, 50.0)
      @test xray_luminosity(radiation) ≈ 9.074006146667613e+44
      @test eddington_luminosity(radiation) ≈ 1.2570651798467906e+46
      @test bolometric_luminosity(radiation) ≈ 6.285325899233953e+45
      @test compute_mass_accretion_rate(radiation) ≈ 1.2228101097715833e+26
      @test radiation.shielding_density ≈ 1e8
      @test radiation.r_in ≈ 50
      println(radiation.r_x)
end
