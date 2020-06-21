using Qwind
using Test
import Qwind.compute_radiation_acceleration
import Qwind.compute_initial_acceleration
import Qwind.residual!
abstract type TestRadiation <: Radiation end
struct TestRadiation1 <: TestRadiation
    bh::BlackHole
end
struct TestRadiation2 <: TestRadiation
    bh::BlackHole
end

earth_mass = 5.972e27 # g
earth = BlackHole(earth_mass, 0, 0)
earth_radius = 6371e5 / earth.Rg #cm
earth_gravity = -981 / C^2 * earth.Rg# cm/s^2

function residual!(radiation::TestRadiation, out, du, u, p, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    gravitational_acceleration = compute_gravitational_acceleration(r, z, p.radiation.bh.M)
    radiation_acceleration = compute_radiation_acceleration(p.radiation, du, u, p)
    a_r = gravitational_acceleration[1] + radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - v_r
    out[2] = z_dot - v_z
    out[3] = v_r_dot - a_r
    out[4] = v_z_dot - a_z
end

function compute_radiation_acceleration(radiation::TestRadiation1, du, u, p)
    return [0.0, 0.0]
end
function compute_radiation_acceleration(radiation::TestRadiation2, du, u, p)
    return [0.0, 400.0]
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
    gravitational_acceleration =
        compute_gravitational_acceleration(r, z, params.radiation.bh.M)
    radiation_acceleration =
        compute_radiation_acceleration(radiation, du, u, params)
    a_r =
        gravitational_acceleration[1] +
        radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    return [a_r, a_z]
end

@testset "Test Free fall" begin
    line = Streamline(
        0.0,
        earth_radius + 1e5 / earth.Rg, # free fall from 1km
        0.0,
        0.0,
        1e8,
        1,
        earth_mass,
    )
    grid = Grid(0, 2, 0, 2)
    radiation = TestRadiation1(earth)
    params = Parameters(line, grid, radiation)
    solver = Qwind.initialize_solver(
        line::Streamline,
        params::Parameters;
        atol = 1e-7,
        rtol = 1e-3,
    )
    solve!(solver)
    analytical_solution(t) = z_0(line) - 0.5 * earth_gravity * t^2
    for (t, z) in zip(line.t, line.z)
        true_solution = analytical_solution(t)
        @test z ≈ true_solution rtol = 1e-2
    end
end

@testset "Test Free fall + constant radiation" begin
    line = Streamline(
        0.0,
        earth_radius + 1e5, # free fall from 1km
        0.0,
        0.0,
        1e8,
        1,
        earth_mass,
    )
    constant_radiation(du, u, p) = [0.0, 400]
    grid = Grid(0, 2, 0, 2)
    radiation = TestRadiation2(earth)
    params = Parameters(line, grid, radiation)
    solver = Qwind.initialize_solver(
        line::Streamline,
        params::Parameters;
        atol = 1e-7,
        rtol = 1e-3,
    )
    solve!(solver)
    analytical_solution(t) = z_0(line) - 0.5 * (400 + earth_gravity) * t^2
    for (t, z) in zip(line.t, line.z)
        true_solution = analytical_solution(t)
        @test z ≈ true_solution rtol = 1e-2
    end
end
