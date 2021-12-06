using Qwind
using Test

@testset "Ionization Parameter" begin
    path = String(@__DIR__) * "/../../configs/config_test.yaml"
    model = Model(path)
    @test compute_ionization_parameter(
        model.rad,
        model.wind,
        model.parameters,
        r = 1,
        z = 1,
        number_density=10,
        tau_x=0,
    ) ≈ model.rad.xray_luminosity / 20 / model.bh.Rg^2

    @test compute_ionization_parameter(
        model.rad,
        model.wind,
        model.parameters,
        r = 1,
        z = 1,
        number_density=1e20,
        tau_x=100,
    ) ≈ 1e-20
end

@testset "X-Ray Opacity" begin
    @test compute_xray_opacity(1) == 100 * SIGMA_T
    @test compute_xray_opacity(9e4) == 100 * SIGMA_T
    @test compute_xray_opacity(1e5) == SIGMA_T
    @test compute_xray_opacity(1e6) == SIGMA_T
end
