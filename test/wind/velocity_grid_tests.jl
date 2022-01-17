using Qwind, Test


@testset "Test Velocity Grid" begin

    @testset "Empty grid" begin
        velocity_grid = VelocityGrid(100, 50, 10)
        @test velocity_grid.r_range == [-1, 0]
        @test velocity_grid.z_range == [-1, 0]
        @test velocity_grid.vr_grid == zeros(2, 2)
        @test velocity_grid.vphi_grid == zeros(2, 2)
        @test velocity_grid.vz_grid == zeros(2, 2)
        @test velocity_grid.nr == 100
        @test velocity_grid.nz == 50
    end
    f1(r, z) = 2 * r + 3 * z
    f2(r, z) = r + 2 * z
    rr = range(0, 100, length = 250)
    zz = range(0, 100, length = 250)
    vr_grid = zeros(length(rr), length(zz))
    vz_grid = zeros(length(rr), length(zz))
    vphi_grid = ones(length(rr), length(zz))
    for (i, r) in enumerate(rr)
        for (j, z) in enumerate(zz)
            vr_grid[i, j] = f1(r, z)
            vz_grid[i, j] = f2(r, z)
        end
    end
    vgrid = VelocityGrid(rr, zz, vr_grid, vphi_grid, vz_grid)

    @testset "Test interpolation" begin
        @test vgrid.r_range == rr
        @test vgrid.z_range == zz
        @test vgrid.vr_grid == vr_grid
        @test vgrid.vphi_grid == vphi_grid
        @test vgrid.vz_grid == vz_grid
        @test vgrid.nr == length(rr)
        @test vgrid.nz == length(zz)
        r_test = range(1, 90, length = 100)
        z_test = range(1, 90, length = 100)
        for r in r_test
            for z in z_test
                vel_interp = interpolate_velocity(vgrid, r, z)
                @test vel_interp[1] ≈ f1(r, z) rtol = 1e-2
                @test vel_interp[2] ≈ 1.0 rtol = 1e-2
                @test vel_interp[3] ≈ f2(r, z) rtol = 1e-2

                vel_get = get_velocity(vgrid, r, z)
                @test vel_get[1] ≈ f1(r, z) rtol = 5e-1
                @test vel_get[2] ≈ 1.0 rtol = 5e-1
                @test vel_get[3] ≈ f2(r, z) rtol = 5e-1
            end
        end
    end

    @testset "Test update grid" begin
        new_hull = Hull([0.0, 50, 0.0, 50], [0, 200, 200, 0.0])
        rr = []
        zz = []
        for r in vgrid.r_range
            for z in vgrid.z_range
                push!(rr, r)
                push!(zz, z)
            end
        end
        new_vr = 1 .* ones(length(rr))
        new_vphi = 2 .* ones(length(rr))
        new_vz = 3 .* ones(length(rr))
        new_grid = VelocityGrid(
            new_hull,
            rr,
            zz,
            new_vr,
            new_vphi,
            new_vz,
            rr,
            nr = 50,
            nz = 50,
            log = false,
        )
        for r in new_grid.r_range
            for z in new_grid.z_range[3:end-3] # dont include 0 as it is problematic for the interp.
                if r <= 50
                    @test get_velocity(new_grid, r, z) ≈ [1,2,3]
                else
                    @test get_velocity(new_grid, r, z) ≈ [0,0,0]
                end
            end
        end
    end

end
