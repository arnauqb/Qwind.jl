using Qwind
using Test
using DifferentialEquations

earth_radius = 6371e5 #cm
earth_gravity = -981 # cm/s^2
earth_mass = 5.972e27 # g

@testset "Test Free fall" begin
    line = Streamline(
        0.0,
        earth_radius + 1e5, # free fall from 1km
        0.0,
        0.0,
        1e8,
        1,
        earth_mass)
    params = Parameters(
        line,
        earth_mass,
        0, 
        0.0,
        0.0,
        0.0,
        0,
        2 * earth_radius,
        earth_radius,
        2 * earth_radius,
        nothing,
        NoRad(),
    )
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
        @test z ≈ true_solution rtol=1e-3
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
        earth_mass)
    params = Parameters(
        line,
        earth_mass,
        0, 
        0.0,
        0.0,
        0.0,
        0,
        2 * earth_radius,
        earth_radius,
        2 * earth_radius,
        nothing,
        ConstantRad(),
    )
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
        @test z ≈ true_solution rtol=1e-3
    end
end
