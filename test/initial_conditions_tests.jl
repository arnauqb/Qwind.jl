using Qwind
using Test

@testset "Uniform initial conditions" begin
        ic = UniformIC(1.0, 50.0, 30, 2.0, 1e8, 1e6 / C, false)
        @test getz0(ic, 5) == 2.0
        @test getrin(ic) == 1.0
        @test getrfi(ic) == 50
        @test getnlines(ic) == 30
        @test getn0(ic, 3) == 1e8
        @test getn0(ic, 30) == 1e8
        @test getv0(ic, 3) == 1e6 / C
        @test getv0(ic, 30) == 1e6 / C
end
