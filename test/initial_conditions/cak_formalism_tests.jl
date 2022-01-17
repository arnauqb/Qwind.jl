using Qwind, Test


@testset "Test Nozzle functions" begin

    model = Model(String(@__DIR__) * "/config_test.yaml")
    parameters = model.parameters
    alpha = 0.6
    rr = range(6, 1000, length = 50)
    zz = range(0, 1000, length = 50)
    grid = 1e10 .* ones(length(rr), length(zz))
    dgrid = DensityGrid(rr, zz, grid)
    iterator = GridIterator(rr, zz)
    model.wind = Wind(
        model.wind.wind_hull,
        dgrid,
        model.wind.velocity_grid,
        model.wind.vacuum_density,
        model.wind.update_grid_flag,
    )
    r = 100
    fuv = 0.85
    rad = model.rad
    bh = model.bh

    @testset "Test B0" begin
        @test Qwind.get_B0(2) == 1 / 4
        @test Qwind.get_B0(8) == 1 / 64
    end

    @testset "Test gamma0" begin
        small_rad = Qwind.compute_disc_radiation_field_small_heights(
            dgrid,
            iterator,
            parameters,
            model.rad,
            r_wind = r,
            z_wind = 1.0,
            vr_wind = 0.0,
            vz_wind = 0.0,
        )[2]
        # check attenuation is removed
        gamma0 = Qwind.gamma0(dgrid, iterator, parameters, model.rad, r)
        @test small_rad ≈ small_rad
    end

    @testset "Test f" begin
        cc = 1 / (alpha^alpha * (1 - alpha)^(1 - alpha))
        # close to the disc f should be equal to f0
        # so the return value is just cc and constant.
        for r in range(6.1, 1000.0, length = 50)
            @test Qwind.f(dgrid, iterator, parameters, rad, 1e-2, r = r, alpha = 0.6) ≈
                  cc * fuv rtol = 5e-1
        end
    end

    @testset "Test g" begin
        # test grids are restored
        Qwind.g(dgrid, iterator, parameters, rad, 1, r = 100)
        @test unique(model.wind.density_grid.grid) == [1e10]
        @test rad.fuv_grid == fuv * ones(length(rad.fuv_grid))
        # at small heights is just disk gravity + constant rad, without attentuation
        model.wind.density_grid.grid .= 0.0
        z = 1
        disk_grav = compute_gravitational_acceleration(r, disk_height(bh, r) + z)[2]
        low_rad = Qwind.compute_disc_radiation_field_small_heights(
            dgrid,
            iterator,
            parameters,
            rad,
            r_wind = r,
            z_wind = z,
            vr_wind = 0.0,
            vz_wind = 0.0,
        )[2]
        @test Qwind.g(dgrid, iterator, parameters, rad, z, r = r) ≈
              -(disk_grav + low_rad / fuv) / Qwind.get_B0(100) rtol = 1e-1
        model.wind.density_grid.grid .= 1e10
    end

end

@testset "Test analytical result" begin
    model = Model(String(@__DIR__) * "/config_test.yaml")
    r_range = range(6.1, 1000, length = 50)
    for r in r_range
        calculated =
            disk_height(model.bh, r) + Qwind.find_nozzle_function_minimum(
                model.wind.density_grid,
                model.wind.grid_iterator,
                model.parameters,
                model.rad,
                r,
                zmax = Inf,
                include_a = false,
            )[1]
        expected = r / sqrt(2)
        @test calculated ≈ expected rtol = 1e-1
    end
end
