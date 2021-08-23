using Qwind, Test


@testset "Test Nozzle functions" begin
    model = Model(String(@__DIR__) * "/config_test.yaml")
    alpha = 0.6
    rr = range(6, 1000, length = 50)
    zz = range(0, 1000, length = 50)
    grid = 1e10 .* ones(length(rr), length(zz))
    dgrid = DensityGrid(rr, zz, grid)
    model.rad = update_radiation(model.rad, dgrid)
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
            model.rad,
            r = r,
            z = 1.0,
            vr = 0.0,
            vz = 0.0,
        )[2]
        # check attenuation is removed
        gamma0 = Qwind.gamma0(model.rad, r)
        @test small_rad ≈ small_rad
    end

    @testset "Test f" begin
        cc = 1 / (alpha^alpha * (1 - alpha)^(1 - alpha))
        # close to the disc f should be equal to f0
        # so the return value is just cc and constant.
        for r in range(6.1, 1000.0, length = 50)
            @test Qwind.f(rad, 1e-2, r = r, alpha = 0.6) ≈ cc * fuv rtol = 1e-1
        end
    end

    @testset "Test g" begin
        # test grids are restored
        Qwind.g(rad, 1, r = 100)
        @test unique(rad.wi.density_grid.grid) == [1e10]
        @test rad.fuv_grid == fuv * ones(length(rad.fuv_grid))
        # at small heights is just disk gravity + constant rad, without attentuation
        model.rad.wi.density_grid.grid .= 0.0
        z = 1e-3
        disk_grav = compute_gravitational_acceleration(r, disk_height(bh, r) + z)[2]
        low_rad = Qwind.compute_disc_radiation_field_small_heights(
            rad,
            r = r,
            z = z,
            vr = 0.0,
            vz = 0.0,
        )[2]
        @test Qwind.g(rad, z, r = r) ≈ -(disk_grav + low_rad / fuv) / Qwind.get_B0(100)
        model.rad.wi.density_grid.grid .= 1e10
    end
end


