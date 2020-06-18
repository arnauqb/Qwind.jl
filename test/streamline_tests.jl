using Test
using Qwind

@testset "Test streamline structures" begin
    line = Streamline(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
    @test line.r == [1.0]
    @test line.z == [2.0]
    @test line.v_r == [3.0 / C]
    @test line.v_z == [4.0 / C]
    @test line.number_density == [5.0]
    @test line.tau_x == [0.0]
    @test line.tau_uv == [0.0]
    @test line.force_multiplier == [0.0]
    @test line.dv_dr == [0.0]
    @test line.ionization_parameter == [0.0]
    @test line.angular_momentum == 1.0
    @test line.line_width_norm == 6.0
    push!(line.r, 5)
    push!(line.z, 6)
    @test line_width(line, 5.0) == 30.0
    @test line_width(line) == 30.0
    @test r(line) == 5.0
    @test r_0(line) == 1.0
    @test z(line) == 6.0
    @test z_0(line) == 2.0
    M = 1e8 * M_SUN
    push!(line.v_r, C)
    push!(line.v_z, C)
    @test escaped(line, M) == true
    push!(line.r, 2.5)
    push!(line.z, 2.6)
    @test max_r(line) == 5.0
    @test max_z(line) == 6.0
    push!(line.v_r, 20)
    push!(line.v_z, 30)
    @test v_r(line) == 20
    @test v_r_0(line) == 3.0 / C
    @test v_z(line) == 30
    @test v_z_0(line) == 4.0 / C
    @test line.angular_momentum == sqrt(r_0(line))
    # Streamlines
    line2 = Streamline(8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0)
    lines = Streamlines([line, line2])
    push!(line.r, 2)
    push!(line.z, 1)
    @test initial_radii(lines) == [1.0, 8.0]
    @test get_r_in(lines) == 1.0
    @test max_r(lines) == 8.0
    @test max_z(lines) == 9.0
end

@testset "Create set of streamlines" begin
    r_in = 1
    r_fi = 1e4
    n_lines = 5
    velocity_generator(r0) = 2 * r0
    density_generator(r0) = r0 / 2
    streamlines = Streamlines(
        r_in,
        r_fi,
        n_lines,
        velocity_generator,
        density_generator,
        10,
        log_spaced = true,
    )
    line_beginnings = []
    for line in streamlines.lines
        line_beginning = r_0(line) - line_width(line) / 2.0
        @test v_r(line) == 0
        @test v_z(line) == r_0(line) * 2 / C
        @test number_density_0(line) == r_0(line) / 2
        @test line.angular_momentum == sqrt(r_0(line))
        push!(line_beginnings, line_beginning)
    end
    @test line_beginnings[1] == r_in
    last_line =streamlines.lines[end]
    @test r_0(last_line) + line_width(last_line) / 2 == r_fi
    streamlines = Streamlines(
        0,
        r_fi,
        n_lines,
        velocity_generator,
        density_generator,
        10,
        log_spaced = false,
    )
    line_beginnings = []
    for line in streamlines.lines
        line_beginning = r_0(line) - line_width(line) / 2.0
        @test v_r(line) == 0
        @test v_z(line) == r_0(line) * 2 / C
        @test number_density_0(line) == r_0(line) / 2
        @test line.angular_momentum == sqrt(r_0(line))
        push!(line_beginnings, line_beginning)
    end
    @test line_beginnings ≈ [0, 2000, 4000, 6000, 8000]
    @test initial_radii(streamlines) ≈ [1000, 3000, 5000, 7000, 9000]
end
