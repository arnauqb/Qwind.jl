using Test, Qwind

@testset "test cylindrical distances" begin
    @test compute_distance_cylindrical(10, 1, 2, 20, 2, 3) ≈ 16.884 rtol = 1e-2
end
