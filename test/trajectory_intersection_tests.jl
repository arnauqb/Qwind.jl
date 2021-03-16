using Qwind, Test

@testset "Test segments intersection" begin
    # do intersect
    segment1 = Segment(2.3, 7.99, 10.64, 3.93)
    segment2 = Segment(2.86, 3.45, 11.0, 7.0)
    do_intersect, intersection = intersect(segment1, segment2)
    @test do_intersect == true
    @test intersection ≈  [7.48, 5.47] rtol = 0.01
    # don't intersect
    segment1 = Segment(5, -2, 9.4, -2.79)
    segment2 = Segment(4.82, -5.83, 7.2, -3.41)
    do_intersect, intersection = intersect(segment1, segment2)
    @test do_intersect == false
    @test intersection == [Inf, Inf]
end

@testset "Test curves intersection" begin
    # do intersect
    r1 = [2.7, 4.08, 9.18, 11.72]
    z1 = [5.43, 12.17, 12.21, 6.65]
    r2 = [5.74, 7.6, 12.58, 15.8]
    z2 = [6.19, 14.51, 14.33, 6.09]
    do_intersect, intersection = intersect(r1, z1, r2, z2)
    @test do_intersect == true
    @test intersection ≈ [7.08, 12.19] rtol = 0.1
    # don't intersect
    r1 = [7.39, 11.88, 16.0, 11.64]
    z1 = [-2.7, 1.27, -2, -4.92]
    r2 = [8.85, 11.41, 14.52]
    z2 = [-2.28, -1.16, -1.97]
    do_intersect, intersection = intersect(r1, z1, r2, z2)
    @test do_intersect == false
    @test intersection == [Inf, Inf]
end
