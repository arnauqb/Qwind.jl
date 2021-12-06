using Qwind, Test

@testset "xi grid construction tests" begin
    r_range = range(0,100,length=100)
    z_range = range(0,100,length=101)
    xi_grid = 1e5 .* ones(length(r_range), length(z_range))
    ion_grid = IonizationGrid(r_range, z_range, xi_grid)
    @test ion_grid.r_range == r_range
    @test ion_grid.z_range == z_range
    @test ion_grid.grid == xi_grid
    @test ion_grid.nr == 100
    @test ion_grid.nz == 101
end

@testset "xi grid interpolation tests" begin
    r_range = range(0,100,length=100)
    z_range = range(0,100,length=101)
    f(r, z) = 1e5 .* (r + z)
    xi_grid = zeros(length(r_range), length(z_range))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            xi_grid[i, j] = f(r, z)
        end
    end
    ion_grid = IonizationGrid(r_range, z_range, xi_grid)
    for r in 1:100
        for z in 1:100
            @test get_ionization_parameter(ion_grid, r, z) ≈ f(r, z)
            @test interpolate_ionization_parameter(ion_grid, r, z) ≈ f(r, z)
        end
    end
end
