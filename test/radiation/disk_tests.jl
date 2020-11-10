using Qwind
using Test

@testset "Disk temperature" begin
    bh = BlackHole(1e8, 0.5, 0)
    radiation = QsosedRadiation(bh, 1000, 0.15)
    @test compute_disk_temperature(radiation, 100)^4 ≈ 6.920770795914223e+17 rtol = 1e-5
    @test compute_disk_temperature(radiation, 10)^4 ≈ 1.2061136229423976e+20 rtol = 1e-5
end