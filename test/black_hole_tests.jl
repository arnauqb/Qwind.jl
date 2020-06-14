using Qwind
using Test

@testset "Test Black Hole properties" begin
    black_hole = BlackHole(1e8 * M_SUN, 0.5, 0.5)
    @test solar_mass(black_hole) ≈ 1e8
    @test eddington_luminosity(black_hole) ≈ 1.25750e46 rtol = 1e-3 #erg/s
    @test bolometric_luminosity(black_hole) ≈ 0.5 * eddington_luminosity(black_hole) rtol =
        1e-3
    #@test mass_accretion_rate(black_hole) ≈ 1.39916e26 rtol = 1e-4
    @test Rg(black_hole) ≈ 14766250380501.246 rtol = 1e-4
    @test Rs(black_hole) ≈ 2 * Rg(black_hole) rtol = 1e-4
    @test isco(black_hole) ≈ 62505575212347.3  rtol = 1e-4
    black_hole2 = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    @test isco(black_hole2) ≈ 88597502283007.47 rtol = 1e-4
    @test gravitational_acceleration(3, 4, black_hole) ≈ 1e8 * M_SUN * G * [3/5^3, 4/5^3]
end
