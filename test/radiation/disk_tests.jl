using Qwind, Test

@testset "Test relativistic AD correction" begin
    bh = BlackHole(1e8, 0.5, 0.6)
    @test disk_nt_rel_factors(bh, 1) ≈ 0.0
    @test disk_nt_rel_factors(bh, 10) ≈ 0.2615157665830733
    @test disk_nt_rel_factors(bh, 50) ≈ 0.6086578218842061
    @test disk_nt_rel_factors(bh, 200) ≈ 0.7849369908504347
    bh = BlackHole(1e9, 0.5, 0.99)
    @test disk_nt_rel_factors(bh, 10) ≈ 0.49104425779640487
    @test disk_nt_rel_factors(bh, 50) ≈ 0.7072886586539658
    @test disk_nt_rel_factors(bh, 200) ≈ 0.8339124564014411
end

@testset "gravity radius" begin
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    @test gravity_radius(bh) ≈ 1580 rtol=0.01
end
