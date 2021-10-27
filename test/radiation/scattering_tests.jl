using Qwind, Test

@testset "test rectangle" begin
    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
    @test cell.rmin == 5
    @test cell.rmax == 10
    @test cell.zmin == 4
    @test cell.zmax == 11
end

@testset "compute intersection with cell" begin

    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)

    @testset "angles" begin
        @test Qwind.compute_lowest_intercepting_angle(cell) ≈ 0.427 rtol = 1e-2
        @test Qwind.compute_highest_intercepting_angle(cell) ≈ 1.190 rtol = 1e-2
    end

    @testset "intersections" begin
        theta = 0.4638
        @test Qwind.compute_first_intersection(cell, theta) ≈ [5, 10] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta) ≈ [5.5, 11] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta) ≈ 1.118 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta) ≈ 11.18 rtol = 1e-2
    end
end

@testset "scattering optical depths" begin
    rr = range(0, 100, length = 100)
    zz = range(0, 100, length = 100)
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    grid_values = 1e8 .* ones(length(rr), length(zz))
    dgrid = DensityGrid(rr, zz, grid_values)
    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
    theta = 0.4638
    @test Qwind.compute_cell_optical_depth(
        cell,
        theta,
        density = 1e8,
        Rg = bh.Rg,
        mu_electron = 1.17,
    ) ≈ 1.118 * SIGMA_T * 1.17 * 1e8 * bh.Rg rtol = 2e-1
    @test Qwind.compute_optical_depth_to_cell(
        dgrid,
        rectangle = cell,
        theta = theta,
        Rg = bh.Rg,
        mu_electron = 1.17,
    ) ≈ 11.18 * 1.17 * SIGMA_T * 1e8 * bh.Rg rtol = 1e-2

    @testset "cell absorption" begin
        lumin = compute_bolometric_luminosity(bh)
        cell = Rectangle(rmin = 50, rmax = 60, zmin = 0, zmax = 1e-4)
        theta = π / 2
        cell_density = 1e10
        tau_to_cell = Qwind.compute_optical_depth_to_cell(
            dgrid,
            rectangle = cell,
            theta = theta,
            Rg = bh.Rg,
            mu_electron = 1.17,
        )
        tau_cell = Qwind.compute_cell_optical_depth(
            cell,
            theta,
            density = cell_density,
            Rg = bh.Rg,
            mu_electron = 1.17,
        )
        delta_theta = atan(cell.zmax / cell.rmin)
        expected = lumin / (4π) * exp(-tau_to_cell) * (1 - exp(-tau_cell)) * delta_theta * 2π 
        @test Qwind.compute_luminosity_absorbed_by_cell(
            dgrid,
            cell = cell,
            Rg = bh.Rg,
            mu_electron = 1.17,
            cell_density = cell_density,
            source_luminosity = lumin,
        ) ≈ expected rtol = 1e-1
    end

end

