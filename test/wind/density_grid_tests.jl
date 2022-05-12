using Qwind, Test

@testset "Get  spatial grid" begin
    r = [1.0, 2.0, 30.0, 4.0, 105.0]
    z = [1.0, 2.0, 3.0, 4.0, 3.0]
    r0s = [1.0, 2.0, 3.0]
    auto_r, auto_z = Qwind.get_spatial_grid(r, z, r0s, "auto", 20)
    @test auto_r[1:3] == r0s
    for i = 1:(length(auto_r) - 1)
        @test auto_r[i + 1] > auto_r[i]
    end
    @test auto_r[end] == 103.0
    @test auto_z[1] == 1.0
    @test auto_z[end] == 4
    @test length(auto_z) == 19 # because one is left for z = 0
    for i = 1:(length(auto_z) - 1)
        @test auto_z[i + 1] > auto_z[i]
    end
end

@testset "Test density grid" begin

    r_range = Float64[]
    z_range = Float64[]
    for r in range(10, 100, length = 50)
        for z in range(10, 100, length = 50)
            push!(r_range, r)
            push!(z_range, z)
        end
    end
    hull = Hull([0.01, 200, 0.01, 200], [0.01, 200, 200, 0.01])
    r0_range = collect(range(0.1, 100.0, step = 10))
    func(r, z) = 1e12 * exp(-z / 100) / r
    n_range = func.(r_range, z_range)
    grid = DensityGrid(r_range, z_range, n_range, r0_range, hull, nr = 50, nz = 50)

    @testset "Test gridpoints" begin
        # Grid points should be exact?
        for (i, r) in enumerate(grid.r_range)
            for (j, z) in enumerate(grid.z_range)
                @test grid.grid[i, j] ≈ func(r, z) rtol = 0.1
                @test get_density(grid, r, z) ≈ func(r, z) rtol = 0.1
            end
        end
    end

    @testset "Test interpolation" begin
        for r in range(grid.r_range[1], grid.r_range[end], length = 100)
            for z in range(grid.z_range[1], grid.z_range[end], length = 100)
                @test func(r, z) ≈ interpolate_density(grid, r, z) rtol = 0.1
            end
        end
    end

    @testset "Test update grid average method" begin
        new_hull = Hull([0.01, 50, 0.01, 50], [0.01, 200, 200, 0.01])
        new_density = 1e8 .* ones(length(r_range))
        new_grid = Qwind.update_density_grid(
            grid,
            Qwind.AverageGrid(),
            new_hull,
            r_range,
            z_range,
            r0_range,
            new_density,
        )
        for r in new_grid.r_range
            for z in new_grid.z_range
                if r <= 50
                    @test get_density(new_grid, r, z) ≈
                          10^((8 + log10(get_density(grid, r, z))) / 2)
                else
                    @test get_density(new_grid, r, z) ≈
                          10 .^ ((2 + log10(get_density(grid, r, z))) / 2)
                end
            end
        end
    end

    @testset "Test update grid replace method" begin
        new_hull = Hull([0.01, 50, 0.01, 50], [0.01, 200, 200, 0.01])
        new_density = 1e8 .* ones(length(r_range))
        new_grid = Qwind.update_density_grid(
            grid,
            Qwind.ReplaceGrid(),
            new_hull,
            r_range,
            z_range,
            r0_range,
            new_density,
        )
        for r in new_grid.r_range
            for z in new_grid.z_range
                if r <= 50
                    @test get_density(new_grid, r, z) ≈ 1e8
                else
                    @test get_density(new_grid, r, z) ≈ 1e2
                end
            end
        end
    end

end


