using Test
using Qwind
import Qwind.radiation_force_integrand!

struct RadiationTest <: Radiation end

function radiation_force_integrand!(radiation::RadiationTest, v, r_d, phi_d, r, z)
    v[1] = 2 * r_d * cos(phi_d) * r_d * r
    v[2] = 2 * r_d^2 * sin(phi_d) * r_d * z
end
@testset "integrand integration" begin
    result, error = integrate_radiation_force_integrand(
        RadiationTest(),
        2.0,
        3.0,
        r_lims = (1.0, 2.0),
        phi_lims = (0.0, π / 4),
        rtol = 1e-4,
    )
    @test result[1] ≈ 14 / 3 * sqrt(2)
    @test error[1] < 5e-3 * result[1]
    @test result[2] ≈ 15 / 2 * (1 - sqrt(2) / 2) * 3
    @test error[2] < 5e-3 * result[2]
end


