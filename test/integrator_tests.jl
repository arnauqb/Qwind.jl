using Qwind
using Test
import Qwind.compute_radiation_acceleration
import Qwind.compute_initial_acceleration
import Qwind.residual!
abstract type TestRadiation{T} <: RadiativeTransfer{T} end
struct TestRadiation1{T} <: TestRadiation{T}
    bh::BlackHole{T}
end
struct TestRadiation2{T} <: TestRadiation{T}
    bh::BlackHole{T}
end


@testset "Get initial radii linewidths" begin
    # linear
    ic = UniformIC(0.0, 100.0, 25, 2.0, 1e8, 1e6 / C, false)
    bh = BlackHole(1e8 *M_SUN, 0.5, 0.0)
    Rg = bh.Rg
    xray_lumin = 1e40
    lines_range, lines_widths = get_initial_radii_and_linewidths(ic, xray_lumin, Rg)
    @test length(lines_range) == 25
    @test length(lines_widths) == 25
    @test lines_widths == 4 .* ones(25)
    @test lines_range[1] == 2
    @test lines_range[2] == 6
    @test lines_range[end-1] == 94
    @test lines_range[end] == 98
    # log 
    ic = UniformIC(1.0, 1e4, 4, 2.0, 1e8, 1e6 / C, true)
    lines_range, lines_widths = get_initial_radii_and_linewidths(ic, xray_lumin, Rg)
    @test diff(log10.(lines_widths)) ≈ ones(3)
    @test lines_range[1] == 5.5
    @test lines_range[2] == 55.0
    @test lines_range[3] == 550.0
    @test lines_range[4] == 5500.0
    @test lines_widths[1] == 9
    @test lines_widths[2] == 90
    @test lines_widths[3] == 900
    @test lines_widths[4] == 9000
end

earth_mass = 5.972e27 # g
earth = BlackHole(earth_mass, 0, 0)
earth_radius = 6371e5 / earth.Rg #cm
earth_gravity = -981 / C^2 * earth.Rg# cm/s^2

function residual!(radiation::TestRadiation, out, du, u, p, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    gravitational_acceleration = compute_gravitational_acceleration(r, z)
    radiation_acceleration = compute_radiation_acceleration(radiation, du, u, p)
    a_r = gravitational_acceleration[1] + radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - v_r
    out[2] = z_dot - v_z
    out[3] = v_r_dot - a_r
    out[4] = v_z_dot - a_z
end

function compute_radiation_acceleration(radiation::TestRadiation1, du, u, p::Parameters)
    return [0.0, 0.0]
end
function compute_radiation_acceleration(radiation::TestRadiation2, du, u, p::Parameters)
    return [0.0, 400 / C^2 * earth.Rg]
end

function compute_initial_acceleration(
    radiation::TestRadiation,
    r,
    z,
    v_r,
    v_z,
    params::Parameters,
)
    u = [r, z, v_r, v_z]
    du = [v_r, v_z, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(r, z)
    radiation_acceleration = compute_radiation_acceleration(radiation, du, u, params)
    a_r = gravitational_acceleration[1] + radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    return [a_r, a_z]
end

@testset "Test Free fall" begin
    ic = UniformIC(0.0, 1000.0, 1, earth_radius + 1e4 / earth.Rg, 1e8, 0.0, false)
    grid = Rectangular(-5.0, 1e10, -5.0, 1e10)
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    radiation = TestRadiation1(earth)
    integrator = initialize_integrator(
        radiation,
        bh,
        grid,
        ic,
        0.0,
        0.0;
        atol = 1e-7,
        rtol = 1e-4,
        tmax = 1e2 * C / earth.Rg,
        save_results = false
    )
    run_integrator!(integrator)
    analytical_solution(t) = 0.5 * earth_gravity * t^2
    t_range = range(0, integrator.t, length = 50)
    for t in t_range
        z_solv = integrator.sol(t)[2] - getz0(ic, 0.0)
        z_analy = analytical_solution(t)
        @test z_solv ≈ z_analy rtol = 5e-2
    end
end

@testset "Test Free fall + constant radiation" begin
    ic = UniformIC(0.0, 1000.0, 1, earth_radius + 1e4 / earth.Rg, 1e8, 0.0, false)
    grid = Rectangular(-5.0, 1e10, -5.0, 1e10)
    radiation = TestRadiation2(earth)
    integrator = initialize_integrator(
        radiation,
        grid,
        ic,
        0.0,
        0.0;
        atol = 1e-7,
        rtol = 1e-4,
        tmax = 1e2 * C / earth.Rg,
        save_results = false
    )
    run_integrator!(integrator)
    analytical_solution(t) = 0.5 * ((400 / C^2 * earth.Rg) + earth_gravity) * t^2
    t_range = range(0, integrator.t, length = 50)
    for t in t_range
        z_solv = integrator.sol(t)[2] - getz0(ic, 0.0)
        z_analy = analytical_solution(t)
        @test z_solv ≈ z_analy rtol = 5e-2
    end
end
