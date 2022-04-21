using Qwind, Test

@testset "Segment basic functions" begin
    s = Segment(1, 2, 3, 0)
    @test Qwind.min_x(s) == 1
    @test Qwind.max_x(s) == 3
    @test Qwind.min_y(s) == 0
    @test Qwind.max_y(s) == 2

    @test !Qwind.is_singular(s)
    s = Segment(1, 2, 1, 0)
    @test Qwind.is_singular(s)
    s = Segment(2, 1, 1, 1)
    @test Qwind.is_singular(s)

end

@testset "Intersection basic functions" begin
    i1 = Intersection("loser", 1, 2, 3.0, 4.0)
    i2 = Intersection("winer", 5, 6, 7.0, 8.0)
    @test i1 < i2
    i1 = Intersection("loser", 1, 2, 3.0, 7.0)
    i2 = Intersection("winer", 5, 6, 3.0, 4.0)
    @test i2 < i1
    i1 = Intersection("loser", 1, 2, 3.0, 4.0)
    i2 = Intersection("winer", 5, 6, 3.0, 4.0)
    @test i1 < i2
    i1 = Intersection("loser", 1, 2, 3.0, 4.0)
    i2 = Intersection("winer", 1, 2, 3.0, 4.0)
    @test !(i1 < i2)
    @test !(i2 < i1)
end

@testset "Test intersections" begin
    @testset "Trivial miss" begin
        s1 = Segment(1, 2, 3, 4)
        s2 = Segment(5, 2, 6, 4)
        @test Qwind.trivial_miss(s1, s2)
        s2 = Segment(1, 2, 3, 4)
        s1 = Segment(5, 2, 6, 4)
        @test Qwind.trivial_miss(s1, s2)
        s1 = Segment(1, 1, 3, 2)
        s2 = Segment(1, 5, 3, 6)
        @test Qwind.trivial_miss(s1, s2)
        s2 = Segment(1, 1, 3, 2)
        s1 = Segment(1, 5, 3, 6)
        @test Qwind.trivial_miss(s1, s2)
    end
    @testset "Simple segment intersection" begin
        A = zeros((2, 2))
        b = zeros(2)
        # do intersect
        segment1 = Segment(2.3, 7.99, 10.64, 3.93)
        segment2 = Segment(2.86, 3.45, 11.0, 7.0)
        do_intersect = intersect!(segment1, segment2, A, b)
        @test do_intersect == true
        # don't intersect
        segment1 = Segment(5, -2, 9.4, -2.79)
        segment2 = Segment(4.82, -5.83, 7.2, -3.41)
        do_intersect = intersect!(segment1, segment2, A, b)
        @test do_intersect == false
    end
    @testset "Test curves intersection" begin
        # do intersect
        r1 = [2.7, 4.08, 9.18, 11.72]
        z1 = [5.43, 12.17, 12.21, 6.65]
        r2 = [5.74, 7.6, 12.58, 15.8]
        z2 = [6.19, 14.51, 14.33, 6.09]
        index = intersect(r1, z1, r2, z2)
        @test index == 2
        # don't intersect
        r1 = [7.39, 11.88, 16.0, 11.64]
        z1 = [-2.7, 1.27, -2, -4.92]
        r2 = [8.85, 11.41, 14.52]
        z2 = [-2.28, -1.16, -1.97]
        index = intersect(r1, z1, r2, z2)
        @test index == 4
    end
end

@testset "Test get intersection times" begin
    trajectories = Qwind.Trajectory[]
    for i = 1:5
        ifloat = Float64(i)
        traj = Qwind.Trajectory(
            i,
            [0, ifloat, 2 * ifloat],
            [ifloat, ifloat, 0.5 * ifloat],
            [0, 1.0, 5.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1e8, 1e8, 1e8],
            false,
        )
        push!(trajectories, traj)
    end
    n = 10
    traj = Qwind.Trajectory(
        0,
        collect(range(0.0, 10.0, length = n)),
        collect(range(-1.0, 5.0, length = n)),
        0.8 .* collect(range(0, 5, length = n)),
        50.0 .* ones(n),
        50.0 .* ones(n),
        50.0 .* ones(n),
        1e8 .* ones(n),
        false,
    )
    push!(trajectories, traj)
    times, _ = Qwind.get_intersection_times(trajectories)
    for i = 1:5
        @test times[i] == i
    end
    @test times[0] == 10
end

@testset "Test determine line widths" begin
    A = zeros((2, 2))
    b = zeros(2)

    @testset "Test distance between segments" begin
        s1 = Segment(4, 0, 7, 7)
        s2 = Segment(15, 0, 12, 6)
        s3 = Segment(7, 7, 10, 10)
        distance = Qwind.get_distance_between_segments(s1, s2, A, b)
        @test distance == Inf
        distance = Qwind.get_distance_between_segments(s3, s2, A, b)
        @test distance ≈ sqrt(2) * 4.5
    end

    @testset "Test get distance pairs" begin
        sl1 = Streamline(
            1,
            zeros(5),
            [2.0, 2.35, 3.78, 6.12, 8.6],
            [0.0, 2.3, 4.58, 6.6, 7.99],
            zeros(5),
            zeros(5),
            zeros(5),
            zeros(5),
            zeros(5),
            true,
        )
        sl2 = Streamline(
            1,
            zeros(4),
            [5, 5.9, 7.48, 9.92],
            [0, 2.62, 4.48, 5.48],
            zeros(4),
            zeros(4),
            zeros(4),
            zeros(4),
            zeros(4),
            true,
        )
        sl3 = Streamline(
            1,
            zeros(4),
            [8, 9.67, 11.94, 14.54],
            [0, 2.16, 3.57, 4.08],
            zeros(4),
            zeros(4),
            zeros(4),
            zeros(4),
            zeros(4),
            true,
        )
        segments = Qwind.turn_streamline_to_segments(sl3)
        @test segments[1] == Segment(8, 0, 9.67, 2.16)
        @test segments[2] == Segment(9.67, 2.16, 11.94, 3.57)
        @test segments[3] == Segment(11.94, 3.57, 14.54, 4.08)
        linewidth_pairs = Qwind.get_width_pairs(sl2, sl1, sl3, A, b)
        @test linewidth_pairs ≈ [[3.24, 2.98], [2.72, 3.21], [2.69, 2.94]] rtol = 1e-2

    end

end
