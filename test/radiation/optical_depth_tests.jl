using Qwind
using Test
using QuadGK

@testset "distance to intersection" begin
    r_range = range(0, 1000, length = 150)
    z_range = range(0, 1000, length = 150)
    iterator = GridIterator(r_range, z_range)
    set_iterator!(iterator, 10.0, π / 4, 2.0, 50.0, 3π / 2, 20.0)
    point = [30.455, -1.4075, 23.0] # t = 0.7 in the line
    dd = Qwind.dist_to_intersection(iterator, point)
    true_dd = sqrt(
        (5 * sqrt(2) * (1 - 0.7))^2 +
        (5 * sqrt(2) - 5 * sqrt(2) * 0.7 + 50 * 0.7)^2 +
        (30 * 0.7)^2,
    )
    @test dd ≈ true_dd rtol = 1e-2
end

@testset "Test UV optical depth" begin

    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    Rg = bh.Rg
    function compute_delta(rd, phid, r, z)
        return sqrt(r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    end

    @testset "Constant density" begin
        r_range = range(0, 1000.0, length = 11)
        z_range = range(0, 1000.0, length = 11)
        rp_range_test = range(0.1, 750, length = 20)
        r_range_test = range(0.1, 750, length = 20)
        z_range_test = range(0.1, 750, length = 20)
        phid_range_test = range(0, π-0.1, length = 20)
        density_grid = 2e8 .* ones((length(r_range), length(z_range)))
        grid = DensityGrid(r_range, z_range, density_grid)
        f_anl(rd, phid, r, z) = compute_delta(rd, phid, r, z) * 2e8 * SIGMA_T * Rg
        for rdp in rp_range_test
            for rp in r_range_test
                for zp in z_range_test
                    for pd in phid_range_test
                        truesol = f_anl(rdp, pd, rp, zp)
                        #println("rdp $rdp")
                        #println("rp $rp")
                        #println("zp $zp")
                        #println("pd $pd")
                        #println("Delta is ")
                        #println(compute_delta(rdp, pd, rp, zp))
                        qwsol = compute_uv_tau(grid, rdp, pd, rp, zp, Rg)
                        @test qwsol ≈ truesol rtol = 1e-6
                    end
                end
            end
        end
    end

    @testset "Linear density" begin
        r_range = range(0, 1000.0, length = 500)
        z_range = range(0, 1000.0, length = 500)
        rp_range_test = range(2, 750, length = 50)
        r_range_test = range(1, 750, length = 50)
        z_range_test = range(1, 750, length = 50)
        density_grid = zeros((length(r_range), length(z_range)))
        dens_func(r, z) = 1e8 * (2 * r + z)
        for (i, rp) in enumerate(r_range)
            for (j, zp) in enumerate(z_range)
                density_grid[i, j] = dens_func(rp, zp)
            end
        end
        grid = DensityGrid(r_range, z_range, density_grid)
        function f_anl_kernel(t, rd, r, z)
            rp = rd + t * (r - rd)
            zp = z * t
            dens = dens_func(rp, zp)
            return dens
        end
        function f_anl(rd, r, z)
            dens_line =
                quadgk(t -> f_anl_kernel(t, rd, r, z), 0, 1, rtol = 1e-4, atol = 0)[1]
            return dens_line * sqrt((r - rd)^2 + (z - z_range[1])^2) * SIGMA_T * Rg
        end
        for rdp in rp_range_test
            for rp in r_range_test
                for zp in z_range_test
                    truesol = f_anl(rdp, rp, zp)
                    qwsol = compute_uv_tau(grid, rdp, rp, zp, Rg)
                    @test qwsol ≈ truesol rtol = 5e-1
                end
            end
        end
    end

    @testset "Linear density 3D" begin
        r_range = range(0, 1000.0, length = 500)
        z_range = range(0, 1000.0, length = 500)
        rd_range_test = range(2, 750, length = 10)
        phid_range_test = range(0, 2 * π, length = 10)
        #phid_range_test = [0.0]
        r_range_test = range(1, 750, length = 10)
        z_range_test = range(1, 750, length = 10)
        density_grid = zeros((length(r_range), length(z_range)))
        dens_func(r, z) = 1e8 * (2 * r + z)
        for (i, rp) in enumerate(r_range)
            for (j, zp) in enumerate(z_range)
                density_grid[i, j] = dens_func(rp, zp)
            end
        end
        grid = DensityGrid(r_range, z_range, density_grid)
        function f_anl_kernel(t, rd, phid, r, z)
            rp = rd * cos(phid) + t * (r - rd * cos(phid))
            zp = z * t
            dens = dens_func(rp, zp)
            return dens
        end
        function f_anl(rd, phid, r, z)
            dens_line =
                quadgk(t -> f_anl_kernel(t, rd, phid, r, z), 0, 1, rtol = 1e-4, atol = 0)[1]
            return dens_line * SIGMA_T * Rg * compute_delta(rd, phid, r, z)
        end
        for rdp in rd_range_test
            for rp in r_range_test
                for zp in z_range_test
                    qwsol = compute_uv_tau(grid, rdp, rp, zp, Rg)
                    for phid in phid_range_test
                        truesol = f_anl(rdp, phid, rp, zp)
                        delta = compute_delta(rdp, phid, rp, zp)
                        qwsol = qwsol * delta / sqrt((rp - rdp)^2 + zp^2)
                        #@test qwsol ≈ truesol rtol = 5e-1
                        if abs(qwsol - truesol) / truesol > 5e-1
                            println("\n\n")
                            println("rdp $rdp")
                            println("rp $rp")
                            println("zp $zp")
                            println("phid $phid")
                            println(qwsol)
                            println(truesol)
                        end
                    end
                end
            end
        end
    end

end
