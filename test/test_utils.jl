using Qwind
using Test

@testset "test get_index" begin
    @test searchsortednearest([1, 2, 3], 2) == 2
    @test searchsortednearest([1, 5, 15], 7) == 2
    @test searchsortednearest([3, 8 ,10], 1) == 1
    @test searchsortednearest([3, 8 ,10], 25) == 3
    @test searchsortednearest([3, 7 ,10], 9) == 3
end
