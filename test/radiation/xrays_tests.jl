using Qwind
using Test

@testset "Ionization Parameter" begin
    @test compute_ionization_parameter(1, 1, 10, 0, 1, 1) ≈ 1 / 20
    @test compute_ionization_parameter(1, 1, 10, 100, 1, 1) ≈ 0 atol = 1e-20
    @test compute_ionization_parameter(1, 1, 10, 0, 2, 2) ≈ 1 / 10 / 4
    # test with rel. corrections
    @test compute_ionization_parameter(1, 1, 1/sqrt(2), 1/sqrt(2), 10, 0, 2, 2) ≈ 1e-20
end

@testset "X-Ray Opacity" begin
    @test compute_xray_opacity(1) == 100 * SIGMA_T
    @test compute_xray_opacity(9e4) == 100 * SIGMA_T
    @test compute_xray_opacity(1e5) == SIGMA_T
    @test compute_xray_opacity(1e6) == SIGMA_T
end
