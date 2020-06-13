using Test
using Qwind

function integrand_tests!(v, r, phi)
    v[1] = 2 * r * cos(phi) * r
    v[2] = 2 * r^2 * sin(phi) * r
end

@testset "2D integrator" begin
    result, error = integrate_2d(
        integrand_tests!,
        x_lims = (1, 2),
        y_lims = (0, π / 4),
        rtol = 1e-4,
    )
    @test result[1] ≈ 14 / 3 * sqrt(2) / 2
    @test error[1] < 5e-3 * result[1]
    @test result[2] ≈ 15 / 2 * (1 - sqrt(2) / 2)
    @test error[2] < 5e-3 * result[2]
end

@testset "Radiation force integrand" begin
    nt_func(r_d, phi_d, r, z) = 2 * r_d + z
    mdot_func(r_d, phi_d, r, z) = z + 2
    uv_func(r_d, phi_d, r, z) = r^2
    geometric_func(r_d, phi_d, r, z) = 1.0
    r_d_func(r_d, phi_d, r, z) = r_d * sin(phi_d)
    phi_d_func(r_d, phi_d, r, z) = r_d
    integrand = make_radiation_force_integrand(
        r_d_func,
        phi_d_func,
        geometric_func,
        nt_func,
        mdot_func,
        uv_func,
    )
    v = [0., 0.]
    integrand(v, 1, π/6, 3, 2)
    @test v[1] ≈ 72
    @test v[2] ≈ 72 * 2
    # now let's integrate this thing
    f(v, r, phi) = integrand(v, r, phi, 3, 2)
    res, err = integrate_2d(f, x_lims=(0,2), y_lims=(0, π/4))
    @test res[1] ≈ 336 * (1- sqrt(2)/2)
    @test res[2] ≈ 336 * π / 4
end

@testset "Default radiation force integrand" begin
    mdot_func(r_d, phi_d, r, z) = 1.0
    uv_func(r_d, phi_d, r, z) = 1.0
    nt_func(r_d, phi_d, r, z) = nt_rel_factors(r_d, 0.0, 6.0)
    geometric_func(r_d, phi_d, r, z) =
        1 / r_d^2 * (r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    r_d_func(r_d, phi_d, r, z) = r - r_d * sin(phi_d)
    phi_d_func(r_d, phi_d, r, z) = 1.0
    integrand = make_radiation_force_integrand(
        r_d_func,
        phi_d_func,
        geometric_func,
        nt_func,
        mdot_func,
        uv_func,
    )
    v = [0., 0.]
    f(v, r, phi) = integrand(v, r, phi, 3, 2)
    res, err = integrate_2d(f, x_lims=(0,2), y_lims=(0, π/4))
    integrand = make_radiation_force_integrand(
        mdot_func=mdot_func,
        uv_func=uv_func,
    )
    v = [0., 0.]
    f(v, r, phi) = integrand(v, r, phi, 3, 2)
    res2, err2 = integrate_2d(f, x_lims=(0,2), y_lims=(0, π/4))
    @test res2[1] == res2[1]
    @test res2[2] == res2[2]
end
