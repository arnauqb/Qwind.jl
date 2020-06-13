using Test
using Qwind

@testset "Test structures" begin
    line = create_line()
    @test line.r == []
    @test line.z == []
    @test line.v_r == []
    @test line.v_z == []
    @test line.number_density == []
    @test line.tau_x == []
    @test line.tau_uv == []
    @test line.force_multiplier == []
    @test line.dv_dr == []
    @test line.ionization_parameter == []
    @test line.angular_momentum == 0.0
    @test line.line_width_norm == 0.0
    line2 = create_line(1.5, 2, 3, 4, 5)
    @test line2.r == [1.5]
    @test line2.z == [2.0]
    @test line2.v_z == [3.0]
    @test line2.number_density == [4.0]
    @test line2.line_width_norm == 5.0 / 1.5
    push!(line2.r, 5)
    push!(line2.z, 6)
    @test r(line2) == 5.0
    @test r_0(line2) == 1.5
    @test z(line2) == 6.0
    @test z_0(line2) == 2.0
    M = 1 / G
    @test escaped(line2, M) == true
    push!(line2.r, 2.5)
    push!(line2.z, 2.6)
    @test max_r(line2) == 5.0
    @test max_z(line2) == 6.0
    push!(line2.v_r, 20)
    push!(line2.v_z, 30)
    @test v_r(line2) == 20
    @test v_r_0(line2) == 0.0
    @test v_z(line2) == 30
    @test v_z_0(line2) == 3.0
    @test line2.angular_momentum == sqrt(1.0 / r_0(line2)) * r_0(line2)
    # Streamlines
    lines = Streamlines([line, line2])
    push!(line.r, 2)
    push!(line.z, 1)
    @test initial_radii(lines) == [2.0, 1.5]
    @test r_in(lines) == 2.0
    @test max_r(lines) == 5.0
    @test max_z(lines) == 6.0

end
