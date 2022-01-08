using Qwind, Test

@testset "rectangle" begin
    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
    @test cell.rmin == 5
    @test cell.rmax == 10
    @test cell.zmin == 4
    @test cell.zmax == 11
    new_rect = Qwind.change_origin(cell, 1, 2)
    @test new_rect.rmin == 4
    @test new_rect.rmax == 9
    @test new_rect.zmin == 2
    @test new_rect.zmax == 9
end


@testset "compute intersection with cell" begin


    @testset "angles" begin
        cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
        angles = Qwind.compute_intercepting_angles(cell)
        @test angles[1] ≈ 0.427 rtol = 1e-2
        @test angles[2] ≈ 1.19 rtol = 1e-2
        cell = Rectangle(rmin = 1, rmax = 4, zmin = 1, zmax = 6)
        angles = Qwind.compute_intercepting_angles(cell, [9.09, 7.27])
        @test angles[1] ≈ 3.823 rtol = 1e-2
        @test angles[2] ≈ 4.557 rtol = 1e-2
    end

    @testset "intersections" begin

        # first quadrant
        cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
        theta1 = 0.533
        @test Qwind.compute_first_intersection(cell, theta1) ≈ [5, 8.47] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta1) ≈ [6.49, 11] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta1) ≈ 2.94 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta1) ≈ 9.84 rtol = 1e-2
        theta2 = 1.061
        @test Qwind.compute_first_intersection(cell, theta2) ≈ [7.16, 4] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta2) ≈ [10, 5.59] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta2) ≈ 3.25 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta2) ≈ 8.2 rtol = 1e-2

        # second quadrant
        cell = Rectangle(rmin = 4, rmax = 10, zmin = -6, zmax = -4)
        theta1 = 2.018
        @test Qwind.compute_first_intersection(cell, theta1) ≈ [8.34, -4] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta1) ≈ [10, -4.8] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta1) ≈ 1.85 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta1) ≈ 9.25 rtol = 1e-2
        theta2 = 2.45
        @test Qwind.compute_first_intersection(cell, theta2) ≈ [4, -4.84] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta2) ≈ [4.96, -6] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta2) ≈ 1.51 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta2) ≈ 6.28 rtol = 1e-2

        # third quadrant
        cell = Rectangle(rmin = -5, rmax = -3, zmin = -10, zmax = -5)
        theta1 = 3.778
        @test Qwind.compute_first_intersection(cell, theta1) ≈ [-3.69, -5] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta1) ≈ [-5, -6.77] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta1) ≈ 2.2 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta1) ≈ 6.22 rtol = 1e-2
        theta2 = 3.504
        @test Qwind.compute_first_intersection(cell, theta2) ≈ [-3, -7.92] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta2) ≈ [-3.79, -10] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta2) ≈ 2.23 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta2) ≈ 8.47 rtol = 1e-2

        # fourth quadrant
        cell = Rectangle(rmin = -10, rmax = -7, zmin = 3, zmax = 6)
        theta1 = 5.378
        @test Qwind.compute_first_intersection(cell, theta1) ≈ [-7, 5.5] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta1) ≈ [-7.64, 6] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta1) ≈ 0.81 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta1) ≈ 8.9 rtol = 1e-2
        theta2 = 5.077
        @test Qwind.compute_first_intersection(cell, theta2) ≈ [-7.87, 3] rtol = 1e-2
        @test Qwind.compute_second_intersection(cell, theta2) ≈ [-10, 3.81] rtol = 1e-2
        @test Qwind.compute_intersection_size(cell, theta2) ≈ 2.28 rtol = 1e-2
        @test Qwind.compute_distance_to_first_intersection(cell, theta2) ≈ 8.42 rtol = 1e-2

    end

end

@testset "scattering optical depths" begin
    rr = range(0, 100, length = 100)
    zz = range(0, 100, length = 100)
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    grid_values = 1e8 .* ones(length(rr), length(zz))
    dgrid = DensityGrid(rr, zz, grid_values)
    iterator = GridIterator(rr, zz)
    cell = Rectangle(rmin = 5, rmax = 10, zmin = 4, zmax = 11)
    theta = 0.4638
    @test Qwind.compute_cell_optical_depth(cell, theta, density = 1e8, Rg = bh.Rg) ≈
          1.118 * SIGMA_T * 1e8 * bh.Rg rtol = 2e-1

    @test Qwind.compute_optical_depth_to_cell(
        dgrid,
        iterator,
        origin = [0, 0, 0],
        rectangle = cell,
        theta = theta,
        Rg = bh.Rg,
        source_luminosity = 1e40,
        absorption_opacity = ThomsonOpacity(),
    ) ≈ 11.18 * SIGMA_T * 1e8 * bh.Rg rtol = 1e-2

    @test Qwind.compute_optical_depth_to_cell(
        dgrid,
        iterator,
        origin = [13.89, 0, 11.25],
        rectangle = Rectangle(rmin = 3, rmax = 7, zmin = 4, zmax = 7),
        theta = 4.218,
        Rg = bh.Rg,
        source_luminosity = 1e40,
        absorption_opacity = ThomsonOpacity(),
    ) ≈ 8.96 * SIGMA_T * 1e8 * bh.Rg rtol = 1e-2

    @testset "cell absorption" begin
        lumin = compute_bolometric_luminosity(bh)
        cell = Rectangle(rmin = 50, rmax = 60, zmin = 0, zmax = 1e-4)
        theta = π / 2
        cell_density = 1e10
        tau_to_cell = Qwind.compute_optical_depth_to_cell(
            dgrid,
            iterator,
            origin = [0, 0, 0],
            rectangle = cell,
            theta = theta,
            Rg = bh.Rg,
            source_luminosity = 1e40,
            absorption_opacity = ThomsonOpacity(),
        )
        tau_cell = Qwind.compute_cell_optical_depth(
            cell,
            theta,
            density = cell_density,
            Rg = bh.Rg,
        )
        delta_theta = atan(cell.zmax / cell.rmin)
        expected = lumin / (4π) * exp(-tau_to_cell) * (1 - exp(-tau_cell)) * delta_theta
        @test Qwind.compute_luminosity_absorbed_by_cell(
            dgrid,
            iterator,
            source_position = [0, 0, 0],
            cell = cell,
            Rg = bh.Rg,
            cell_density = cell_density,
            source_luminosity = lumin,
            absorption_opacity = ThomsonOpacity(),
        ) ≈ expected rtol = 1e-1
    end
end

@testset "test ScatteredLuminosityGrid" begin
    dgrid_r = [1,5,10] #range(0, 10, length = 50)
    dgrid_z = [1,2,3] #range(0, 10, length = 50)
    dgrid = DensityGrid(dgrid_r, dgrid_z, [[1,4,7] [2,5,8] [3,6,9]])
    iterator = GridIterator(dgrid_r, dgrid_z)
    Rg = 1e14
    xl = 1e40
    sgrid = ScatteredLuminosityGrid(
        dgrid,
        iterator,
        Rg = Rg,
        source_luminosity = xl,
        source_position = (0, 0, 1),
    )
    @test sgrid.r_range == [3, 7.5]
    @test sgrid.z_range == [1.5, 2.5]
    @test sgrid.nr == 3
    @test sgrid.nz == 3
    @test Qwind.get_scattered_luminosity(sgrid, 10, 20) == 0
    @test Qwind.get_scattered_luminosity(sgrid, 5, 2) > 0
    @test Qwind.interpolate_scattered_luminosity(sgrid, 10, 20) == 0
    @test Qwind.interpolate_scattered_luminosity(sgrid, 5, 2) > 0
end
