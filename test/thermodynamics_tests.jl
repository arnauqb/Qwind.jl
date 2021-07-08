using Qwind
using Test

@testset "thermal velocity" begin
    @test compute_thermal_velocity(2e6, 0.5) ≈ sqrt(2) * 12848657.328083131 / C
    @test compute_thermal_velocity(25e3, 0.5) ≈ sqrt(2) * 1436523.560259735 / C
end

@testset "updating density" begin
    #line = Streamline(1, 2, 3, 4, 5, 6, 7)
    r = 10
    z = 20
    v_r = 3
    v_z = 4
    r_0 = 2
    v_0 = 1
    n_0 = 10
    z_0 = 0
    @test compute_density(r, z, v_r, v_z, r_0, z_0, v_0, n_0) ≈
          10 * (2 / sqrt(10^2 + 20^2))^2 * (1/5)
end
