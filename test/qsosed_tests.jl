@testset "NT rel factors" begin
    @test nt_rel_factors.([6, 100, 1e8], 0, 6.0) ≈
          [0, 0.6522661452445981, 0.9960299145167933] rtol = 5e-3 atol = 0
    @test nt_rel_factors.(10, 0.99, 6.0) ≈ 0.22827160704457788 rtol = 5e-3 atol =
        0
    @test nt_rel_factors.(10, 0.99, 10) ≈ 0.0 rtol = 5e-3 atol = 0
    @test nt_rel_factors.(500, 0.5, 10) ≈ 0.8296494067756125 rtol = 5e-3 atol =
        0
end
