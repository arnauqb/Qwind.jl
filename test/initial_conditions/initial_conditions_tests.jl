using Qwind
using Test
using YAML

@testset "Uniform initial conditions" begin
    ic = UniformIC(1.0, 50.0, 30, 2.0, 1e8, 1e6, "log")
    @test getz0(ic, 5) == 2.0
    @test getrin(ic) == 1.0
    @test getrfi(ic) == 50
    @test getnlines(ic) == 30
    @test getn0(ic, r=3) == 1e8
    @test getn0(ic, r=30) == 1e8
    @test getv0(nothing, ic, nothing, mu_nucleon = 3) == 1e6 / C
    @test getv0(nothing, ic, nothing, mu_nucleon = 30) == 1e6 / C
end

@testset "CAK initial conditons" begin
    config_test = YAML.load_file(
        String(@__DIR__) * "/cak_config_test.yaml",
        dicttype = Dict{Symbol,Any},
    )
    bh = BlackHole(config_test)
    radiation = Radiation(bh, config_test)
    ic = CAKIC(radiation, config_test)
    @test ic.rin == 20
    @test ic.rfi == 800
    @test ic.nlines == 10
    @test ic.K == 0.03
    @test ic.alpha == 0.6
    @test ic.z0 == 0
    @test getz0(ic, 50) == 0.0
    @test ic.trajs_spacing == "log"
    @test getv0(bh, ic, 50, mu_nucleon = 0.5) ≈
          compute_thermal_velocity(disk_temperature(bh, 50), 0.5)
end
