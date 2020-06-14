using Test
using Qwind

@testset "Test structures" begin
    line = Streamline(1.0, 2.0, 3.0 ,4.0, 5.0, 6.0, 7.0)
    @test line.r == [1.0]
    @test line.z == [2.0]
    @test line.v_r == [3.0]
    @test line.v_z == [4.0]
    @test line.number_density == [5.0]
    @test line.tau_x == [0.0]
    @test line.tau_uv == [0.0]
    @test line.force_multiplier == [0.0]
    @test line.dv_dr == [0.0]
    @test line.ionization_parameter == [0.0]
    @test line.angular_momentum == sqrt(G * 7 * 1)
    @test line.line_width_norm == 6.0
    push!(line.r, 5)
    push!(line.z, 6)
    @test r(line) == 5.0
    @test r_0(line) == 1.0
    @test z(line) == 6.0
    @test z_0(line) == 2.0
    M = 1 / G
    @test escaped(line, M) == true
    push!(line.r, 2.5)
    push!(line.z, 2.6)
    @test max_r(line) == 5.0
    @test max_z(line) == 6.0
    push!(line.v_r, 20)
    push!(line.v_z, 30)
    @test v_r(line) == 20
    @test v_r_0(line) == 3.0
    @test v_z(line) == 30
    @test v_z_0(line) == 4.0
    @test line.angular_momentum == sqrt(7 * G / r_0(line)) * r_0(line)
    # Streamlines
    line2 = Streamline(8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0)
    lines = Streamlines([line, line2])
    push!(line.r, 2)
    push!(line.z, 1)
    @test initial_radii(lines) == [1.0, 8.0]
    @test r_in(lines) == 1.0
    @test max_r(lines) == 8.0
    @test max_z(lines) == 9.0
end
