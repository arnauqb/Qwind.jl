using Qwind, Test

@testset "Relativistic corrections" begin

    @testset "beta" begin
        @test Qwind.compute_beta(0.3, 0.4) ≈ 0.5
        @test Qwind.compute_beta(10, 20) ≈ 1
        @test Qwind.compute_beta(0, 0) ≈ 1e-4
    end

    @testset "gamma" begin
        @test Qwind.compute_gamma(0.5) ≈ 1.154700538
    end

    @testset "rel. dispatch" begin
        r = 10
        z = 5
        z0 = 0
        vr = 0.3
        vz = 0.4
        beta = 0.5
        gamma = 1.15
        r_projection = 1
        delta = 1
        @test relativistic_correction(
            NoRelativistic(),
            r = r,
            z = z,
            z0 = z0,
            vr = vr,
            vz = vz,
            beta=beta,
            gamma = gamma,
            r_projection = r_projection,
            delta = delta,
        ) == 1.0
        beta = sqrt(0.3^2 + 0.3^2)
        @test relativistic_correction(0.0, 0.0, 0.0, 10, 10, 0.3, 0.3) ≈ ((1-beta)/(1+beta))^2
    end
end
