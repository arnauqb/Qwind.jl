using Qwind, Test

@testset "Test segments intersection" begin
    A = zeros((2,2))
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
    index = intersect!(r1, z1, r2, z2)
    @test index == 2
    # don't intersect
    r1 = [7.39, 11.88, 16.0, 11.64]
    z1 = [-2.7, 1.27, -2, -4.92]
    r2 = [8.85, 11.41, 14.52]
    z2 = [-2.28, -1.16, -1.97]
    index = intersect!(r1, z1, r2, z2)
    @test index == 4
end
