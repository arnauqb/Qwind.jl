using Qwind
using Test

@testset "Radiation structure" begin
      black_hole = BlackHole(1e8 * M_SUN, 0.5, 0.0)
      radiation = RERadiation(
            black_hole,
            0.70906799789733695,
            0.14436810902317002,
      )
      @test radiation.xray_luminosity ≈ 9.074006146667613e+44
end

@testset "RE X-ray optical depth" begin
      # r<= r_in
      Rg = 2.0
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 5.0, 100.0, 2.0, Rg) / SIGMA_T == 0.0
      # r < r_0, r_x < r
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 1.0, 2.0, 6.0, Rg) / SIGMA_T ≈
            4010 / cos(4 / 5) atol = 0.0 rtol = 0.1
      # r < r_0, r_x > r
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 1.0, 5.0, 6.0, Rg) / SIGMA_T ≈
            60 / cos(4 / 5) atol = 0.0 rtol = 0.1
      # r > r_0, r_in < r_x < .0r_0.0
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 1.0, 2.0, 3.0, Rg) / SIGMA_T ≈
            6010 / cos(4 / 5) atol = 0.0 rtol = 0.1
      # r > r_0, r_in < r_x < .0r_0.0
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 1.0, 3.5, 3.0, Rg) / SIGMA_T ≈
            2030 / cos(4 / 5) atol = 0.0 rtol = 0.1
      # r > r_0, r_in < r_0 < .0r <.0 r
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 1.0, 5.0, 2.0, Rg) / SIGMA_T ≈
            100 / cos(4 / 5) atol = 0.0 rtol = 0.1
      # r > r_0, r_in < r_0 < .0r_x.0 < r.0
      @test compute_xray_tau(4.0, 3.0, 10.0, 20.0, 1.0, 3.0, 2.0, Rg) / SIGMA_T ≈
            4030 / cos(4 / 5) atol = 0.0 rtol = 0.1
end

@testset "Ionization radius" begin
      xray_lumin = 1e46
      density = 1e8
      r_in = 10
      xi_0 = 1e5
      Rg = 1e14
      r_x = Qwind.compute_ionization_radius(
            xray_lumin,
            density,
            r_in,
            Rg,
            atol = 1e-2,
      )
      @test Qwind.ionization_radius_kernel(
            xray_lumin,
            density,
            r_in,
            r_x,
            xi_0,
            Rg,
      ) ≈ 0 atol = 1e-2
      tau_x = 1.0
      @test compute_ionization_parameter(
            r_x,
            0,
            density,
            tau_x,
            xray_lumin,
            Rg,
      ) ≈ xi_0 atol = 0 rtol = 0.2
end

@testset "UV optical depth" begin
      shielding_density = 1e8
      local_density = 5e7
      r_in = 100
      Rg = 1e14
      r = 300
      z = 400

      r_0 = 200
      # ((200 - 100) * 1e8 + (300-200) * 5e7) / cos(300/500)
      @test compute_uv_tau(
            r,
            z,
            shielding_density,
            local_density,
            r_in,
            r_0,
            Rg,
      ) ≈ 1.81744e10 * SIGMA_T * Rg atol=0 rtol=1e-2

      r_0 = 400
      @test compute_uv_tau(
            r,
            z,
            shielding_density,
            local_density,
            r_in,
            r_0,
            Rg,
      ) ≈ 2.42325e10 * SIGMA_T * Rg atol=0 rtol=1e-2
end
