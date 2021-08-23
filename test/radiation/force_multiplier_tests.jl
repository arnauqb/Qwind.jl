using Qwind
using Test

@testset "Force multiplier" begin
    @test Qwind.compute_force_multiplier_k(1e-3) == 0.411
    @test Qwind.compute_force_multiplier_k(1e-3, FMNoInterp()) ≈ 0.411 rtol=0.1
    @test Qwind.compute_force_multiplier_k(10^3.51) == 0.045
    @test Qwind.compute_force_multiplier_k(1e-5) == 0.411
    @test Qwind.compute_force_multiplier_k(1e10) == 0.013
    @test Qwind.compute_force_multiplier_k(1e10, FMNoInterp()) == 0.03
    @test Qwind.compute_force_multiplier_eta(1e-3) == 10^6.95
    @test Qwind.compute_force_multiplier_eta(1e-3, FMNoInterp()) ≈ 9.33e6 rtol=0.1
    @test Qwind.compute_force_multiplier_eta(10^3.50) == 10^1.58
    @test Qwind.compute_force_multiplier_eta(1e-5) == 10^6.95
    @test Qwind.compute_force_multiplier_eta(1e10) == 10^0.78
    @test Qwind.compute_force_multiplier_eta(1e10, FMNoInterp()) == 1.0
    @test compute_force_multiplier(1e5, 1e5) < 1e-3
    @test compute_force_multiplier(1e5, 1e-5) < 1e-3
    @test compute_force_multiplier(1e-8, 1e-5) > 2000
    @test compute_force_multiplier(1e-8, 1e-5) < 3000
    @test compute_force_multiplier(1e-10, 1e-5) ≈ 2300 atol = 0 rtol = 0.1
    @test compute_force_multiplier(1e-2, 1e-5) ≈ 6 atol = 0 rtol = 0.1
    @test Qwind.compute_tau_eff(1e8, 0) == 1.0
end

