using Qwind
using Test

@testset "Blackbody radiation" begin
    @testset "Test spectral radiance" begin
        # Tested with the Astropy BlackBody model.
        bb = BlackBody(1e5)
        @test spectral_radiance(bb, 0.1) ≈ 2.8706e23 / H_PLANCK rtol = 1e-3
        @test spectral_radiance(bb, 0.02) ≈ 2.74e25 / H_PLANCK rtol = 1e-2
        bb = BlackBody(5e3)
        @test spectral_radiance(bb, 0.001) ≈ 3.425e21 / H_PLANCK rtol = 1e-3
    end
    @testset "Test radiance" begin
    end
end

@testset "Disk radiation" begin
    bh = BlackHole(1e8, 0.5, 0)
    radiation = QsosedRadiation(bh, 1000, 0.15)
    @testset "Disk flux" begin
        @test compute_disk_temperature(radiation, 100)^4 ≈ 6.920770795914223e+17 rtol = 1e-5
        @test compute_disk_temperature(radiation, 10)^4 ≈ 1.2061136229423976e+20 rtol = 1e-5
    end

    @testset "BB spectral radiance" begin
        @test compute_bb_spectral_radiance()
    end

    @testset "Disk UV fractions" begin
    ckend
end

