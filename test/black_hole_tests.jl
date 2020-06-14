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
    @test compute_gravitational_acceleration(3, 4, black_hole) ≈ - 1e8 * M_SUN * G * [3/5^3, 4/5^3]
    @test efficiency(black_hole2) ≈ 0.057190958417936644
    @test mass_accretion_rate(black_hole2) ≈ 1.2228101097715833e+26
    earth_radius = 6371e5 #cm
    earth_gravity = -981 # cm/s^2
    earth = BlackHole(5.972e27, 0.0, 0.0)
    @test compute_gravitational_acceleration(earth_radius, 0, earth) ≈ [earth_gravity, 0.0] atol=0.0 rtol=0.01
end
