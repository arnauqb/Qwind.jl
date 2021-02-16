using Qwind
using Test

@testset "test euclidean distance" begin
    @test d_euclidean(0, 3, 0, 4) == 5
    @test d_euclidean(1, -2, 2, -4) == sqrt(9 + 36)
end

@testset "test get_index" begin
    # Nearest
    @test searchsorted_nearest([1, 2, 3], 2) == 2
    @test searchsorted_nearest([1, 5, 15], 7) == 2
    @test searchsorted_nearest([3, 8, 10], 1) == 1
    @test searchsorted_nearest([3, 8, 10], 25) == 3
    @test searchsorted_nearest([3, 7, 10], 9) == 3
    # First
    @test searchsorted_first([1, 2, 3], 1.8, 1) == 1
    @test searchsorted_first([1, 2, 3], 1.8, -1) == 2
    @test searchsorted_first([1, 2, 3], 5, -1) == 3
    @test searchsorted_first([1, 2, 3], 5, 1) == 3
    @test searchsorted_first([1, 2, 3], 0.5, -1) == 1
    @test searchsorted_first([1, 2, 3], 0.5, 1) == 1
    @test searchsorted_first([1, 2, 3], 2.1, 1) == 2
    @test searchsorted_first([1, 2, 3], 2.1, -1) == 3
end

@testset "count sign changes" begin
    array = [1.0, 1, 1, 1, -1, -1, 1, 1]
    @test countsignchanges(array) == 2
    array = Float64[]
    @test countsignchanges(array) == 0
    array = [1.0, 1, 1, 1, 1]
    @test countsignchanges(array) == 0
    array = [-1.0, -1, -1, -1, -1]
    @test countsignchanges(array) == 0
    array = [1.0, -1, 1, -1, 1]
    @test countsignchanges(array) == 4
end

@testset "Test dict utils" begin
    # iterate paths
    dict = Dict("a" => 1, "b" => 2, "c" => Dict("a" => [1, 2, 3], "b" => Dict("b" => 2)))
    iterator = Set(iter_paths(dict))
    @test length(iterator) == 6
    @test (["a"], 1) in iterator
    @test (["b"], 2) in iterator
    @test (["c", "a"], [1,2,3]) in iterator
    @test (["c", "b", "b"], 2) in iterator
    @test (["c", "b"], Dict("b" =>2)) in iterator
    @test (["c", "a"], [1,2,3]) in iterator
    # get value in path
    @test get_value_in_path(dict, (["a"])) == 1
    @test get_value_in_path(dict, (["b"])) == 2
    @test get_value_in_path(dict, (["c"])) == Dict("a" => [1,2,3], "b" => Dict("b" =>2))
    @test get_value_in_path(dict, (["c", "b", "b"])) == 2
    @test get_value_in_path(dict, (["c", "a"])) == [1,2,3]
    @test get_value_in_path(dict, (["c", "b"])) == Dict("b"=>2)
    # set value in path
    set_value_in_path!(dict, (["c", "b", "b"]), 10)
    @test dict["c"]["b"]["b"] == 10
    set_value_in_path!(dict, (["a"]), [1,2,3])
    @test dict["a"] == [1,2,3]
    # test does nothing without error.
    set_value_in_path!(dict, (["asdasda"]), [1,2,3])
end
