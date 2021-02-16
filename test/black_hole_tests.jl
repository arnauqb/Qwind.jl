using Qwind
using Test

@testset "Test Black Hole properties" begin
    black_hole = BlackHole(1e8 * M_SUN, 0.5, 0.5)
    @test compute_eddington_luminosity(black_hole) ≈ 1.25750e46 rtol = 1e-3 #erg/s
    @test compute_bolometric_luminosity(black_hole) ≈
          0.5 * compute_eddington_luminosity(black_hole) rtol = 1e-3
    #@test mass_accretion_rate(black_hole) ≈ 1.39916e26 rtol = 1e-4
    @test black_hole.Rg ≈ 14766250380501.246 rtol = 1e-4
    @test compute_isco(black_hole) ≈ 4.2330025 rtol = 1e-4
    black_hole2 = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    @test compute_isco(black_hole2) ≈ 6.0 rtol = 1e-4
    @test compute_gravitational_acceleration(3, 4) ≈ -[3 / 5^3, 4 / 5^3]
    @test compute_efficiency(black_hole2) ≈ 0.057190958417936644
    @test compute_mass_accretion_rate(black_hole2) ≈ 1.2228101097715833e+26
    @test compute_escape_velocity(2) ≈ 1
    @test compute_angular_momentum(2) ≈ √(2)
    earth_radius = 6371e5 #cm
    earth_gravity = -981 # cm/s^2
    earth = BlackHole(5.972e27, 0.0, 0.0)
    earth_gravity_rgc2 = earth_gravity / C^2 * earth.Rg
    earth_radius_rg = earth_radius / earth.Rg
    @test compute_gravitational_acceleration(earth_radius_rg, 0) ≈ [earth_gravity_rgc2, 0.0] atol =
        0.0 rtol = 0.01
end

@testset "Initialise black holes" begin
    config = Dict(:black_hole => Dict(:M => 1e9, :mdot => 0.2, :spin => 0.6))
    bh = BlackHole(config)
    @test bh.M == 1e9 * M_SUN
    @test bh.mdot == 0.2
    @test bh.spin == 0.6
end
