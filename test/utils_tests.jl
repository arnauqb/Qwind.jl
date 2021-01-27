using Qwind
using Test

@testset "test get_index" begin
    @test searchsorted_nearest([1, 2, 3], 2) == 2
    @test searchsorted_nearest([1, 5, 15], 7) == 2
    @test searchsorted_nearest([3, 8 ,10], 1) == 1
    @test searchsorted_nearest([3, 8 ,10], 25) == 3
    @test searchsorted_nearest([3, 7 ,10], 9) == 3
end
