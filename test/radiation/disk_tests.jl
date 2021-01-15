using Qwind
using Test

@testset "Disk radiation" begin
    @testset "Blackbody radiation" begin
        @testset "Test spectral radiance" begin
            bb = BlackBody(1e5)
            @test spectral_radiance_frequency(bb, 100) ≈ 3.07235837e-28 rtol = 1e-3
            @test spectral_radiance_frequency(bb, 1e15) ≈ 0.02393854 rtol = 1e-4
            @test spectral_radiance_frequency(bb, 1e18) ≈ 0. rtol = 1e-6 atol = 1e-15
            bb = BlackBody(5e3)
            @test spectral_radiance_frequency(bb, 100) ≈ 1.53617919e-29 rtol = 1e-4
            @test spectral_radiance_frequency(bb, 1e15) ≈ 1.00024067e-06
            @test spectral_radiance_frequency(bb, 1e18) ≈ 0
        end
        @testset "Test integrating energy band" begin
            bb = BlackBody(1e4)
            @test spectral_band_radiance_frequency(bb, 100e12, 1000e12) ≈ 1.30413e11 rtol = 1e-2
            @test spectral_band_radiance(bb, 0.00041357, 0.0041357) ≈ 1.30413e11 rtol = 1e-2
        end
        @testset "Integrating spectral radiance recovers SB law" begin
            bb = BlackBody(1e4)
            @test spectral_band_radiance(bb, 1e-5, 1e2) ≈ SIGMA_SB * bb.T^4 / π rtol = 1e-6
        end
        @testset "Energy band fraction" begin
            bb = BlackBody(1e4)
            @test spectral_band_fraction_frequency(bb, 1e12, 1e18) ≈ 1 rtol = 1e-8 
            @test spectral_band_fraction(bb, 1e-5, 1e2) ≈ 1 rtol = 1e-6
        end
    end
    @testset "Disk radiation" begin
        @testset "Test relativistic AD correction" begin
            bh = BlackHole(1e8, 0.5, 0.6)
            @test disk_nt_rel_factors(bh, 10) ≈ 0.2615157665830733
            @test disk_nt_rel_factors(bh, 50) ≈ 0.6086578218842061
            @test disk_nt_rel_factors(bh, 200) ≈ 0.7849369908504347
            bh = BlackHole(1e9, 0.5, 0.99)
            @test disk_nt_rel_factors(bh, 10) ≈ 0.49104425779640487
            @test disk_nt_rel_factors(bh, 50) ≈ 0.7072886586539658
            @test disk_nt_rel_factors(bh, 200) ≈ 0.8339124564014411
        end
        @testset "Test disc radiated power" begin
            bh = BlackHole(1e8 * M_SUN, 0.5, 0)
            @test disk_flux(bh, 10) ≈ 6839115834162428.0 rtol = 1e-6
            @test disk_flux(bh, 100) ≈ 3.924336168219069e13 rtol = 1e-6
            bh = BlackHole(1e9 * M_SUN, 0.5, 0.9)
            @test disk_flux(bh, 10) ≈ 889138044477300.2 rtol = 1e-6
            @test disk_flux(bh, 100) ≈ 1662600341587.0435 rtol = 1e-6
        end
        @testset "UV fraction" begin
            bh = BlackHole(1e8 * M_SUN, 0.5, 0)
            uvf1 = uv_fraction(bh, 6.01)
            uvf2 = uv_fraction(bh, 10)
            uvf3 = uv_fraction(bh, 100)
            @test uvf1 < uvf2
            @test uvf3 < uvf2
        end
    end
end
