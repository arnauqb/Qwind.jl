using Qwind, Test


@testset "Test Nozzle functions" begin
    model = Model(String(@__DIR__) * "/config_test.yaml") 
    fuv = 0.5
    model.rad.fuv_grid .= fuv #simplifty fuv for testing.
    bh = model.bh
    rt = model.rt
    alpha = 0.6
    @testset "Test B0" begin
        @test Qwind.get_B0(2) == 1/4
        @test Qwind.get_B0(8) == 1/64
    end

    @testset "Test f" begin
        #fuv, _ = Qwind.get_fuv_mdot(model.rad)
        cc = 1 / (alpha^alpha * (1 - alpha)^(1 - alpha))
        # close to the disc f should be equal to f0
        # so the return value is just cc and constant.
        for r in range(6.1, 1000.0, length=50)
            @test Qwind.f(rt, bh, 1e-2, r=r, alpha=0.6) ≈ cc * fuv
        end
    end
    @testset "Test g" begin
        # radiation should be stronger at large z 
        for r in range(20.0, 1000.0, length=50)
            @test Qwind.g(rt, bh, 1000.0, r=r) < 0
        end
        # radiation should also be stronger at small z and smallish r
        for r in range(20.0, 100.0, length=50)
            @test Qwind.g(rt, bh, 1e-10, r=r) < 0
        end
    end
    @testset "Test nozzle function" begin
        # radiation should be stronger at large z 
        for r in range(20.0, 1000.0, length=50)
            @test Qwind.g(rt, bh, 1000.0, r=r) < 0
        end
        # radiation should also be stronger at small z and smallish r
        for r in range(20.0, 100.0, length=50)
            @test Qwind.g(rt, bh, 1e-10, r=r) < 0
        end
    end
end


