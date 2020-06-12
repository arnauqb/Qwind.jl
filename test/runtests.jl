using Qwind
using Test

@test Qwind.my_f(2,1)== 7

@testset "Qwind.jl" begin
    # Write your own tests here.
    Qwind.my_f(2,1)
end
