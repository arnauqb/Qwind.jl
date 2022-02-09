using Qwind, Test

@testset "Test read parameters" begin
    cbase_path = String(@__DIR__) * "/../configs/config_base.yaml"
    p = Parameters(cbase_path)
    @test p.M == 1e8
    @test p.mdot == 0.5
    @test p.spin == 0.0
    @test p.disk_r_in == 6.0
    @test p.disk_r_out == 1600
    @test p.disk_nr == 2000
    @test p.radiation_grid_nr == "auto"
    @test p.radiation_grid_nz == 250
    @test p.z_xray == 0.0
    @test p.z_disk == 0.0
    @test p.mu_nucleon == 0.61
    @test p.mu_electron == 1.17
    @test p.vacuum_density == 1e2
    @test p.integrator_r_min == 6.0
    @test p.integrator_r_max == 1e4
    @test p.integrator_z_min == 0
    @test p.integrator_z_max == 1e4
    @test p.wind_n_trajs == 50
    @test p.wind_trajs_spacing == "log"
    @test p.wind_r_in == 20
    @test p.wind_r_fi == 1500
    @test p.wind_z_0 == 0
    @test p.ic_K == 0.03
    @test p.ic_alpha == 0.6
    @test p.ic_use_precalculated == false
    @test p.n_iterations == 5
    @test p.save_path == "./example"
    @test p.relativistic_flag == Relativistic()
    @test p.uv_opacity_flag == ThomsonOpacity()
    @test p.xray_opacity_flag == BoostOpacity()
    @test p.xray_scattering_flag == NoScattering()
    @test p.tau_uv_calculation_flag == NoTauUV()
    @test p.update_grid_flag == AverageGrid()
    @test p.fm_interp_method_flag == FMInterp()
    @test p.initial_conditions_flag == CAKMode()
    @test p.integrator_atol == 1e-8
    @test p.integrator_rtol == 1e-3
    @test p.disk_integral_atol == 0
    @test p.disk_integral_rtol == 1e-3
    @test p.scattering_atol == 0
    @test p.scattering_rtol == 1e-3
end
