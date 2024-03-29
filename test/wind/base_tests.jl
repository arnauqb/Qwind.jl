using Qwind, Test

@testset "wind interpolator" begin
    config = Dict(
        :radiation => Dict(
            :grid_nr => 10,
            :grid_nz => 11,
            :vacuum_density => 1e3,
            :update_grid_method => "replace",
        ),
    )
    wi = Wind(config)
    @test wi.density_grid.nr == 10
    @test wi.density_grid.nz == 11
    @test wi.velocity_grid.nr == 10
    @test wi.velocity_grid.nz == 11
    @test wi.vacuum_density == 1e3
    @test wi.update_grid_flag == Qwind.ReplaceGrid()
end
