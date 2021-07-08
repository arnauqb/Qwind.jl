using Qwind
using Test
using QuadGK

@testset "Test Tau approximation" begin
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    Rg = bh.Rg
    function compute_delta(rd, phid, r, z)
        return sqrt(r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    end
    r_range = range(0, 100.0, length = 500)
    z_range = range(0, 100.0, length = 500)
    rp_range_test = range(10, 90, length = 10)
    r_range_test = range(10, 90, length = 10)
    z_range_test = range(10, 90, length = 10)
    phid_range_test = range(0, π-0.1, length = 10)
    density_grid = zeros((length(r_range), length(z_range)))
    dens_func(r, z) = 1e8 * (2 * r + z)
    for (i, rp) in enumerate(r_range)
        for (j, zp) in enumerate(z_range)
            density_grid[i, j] = dens_func(rp, zp)
        end
    end
    grid = DensityGrid(r_range, z_range, density_grid);
    function f_anl_kernel(t, rd, phid, r, z)
        x = rd * cos(phid) * (t-1) + r * t 
        y = rd * sin(phid) * (t-1)
        zp = z * t
        rp = sqrt(x^2 + y^2)
        dens = dens_func(rp, zp)
        return dens
    end
    function f_anl(rd, phid, r, z)
        dens_line =
            quadgk(t -> f_anl_kernel(t, rd, phid, r, z), 0, 1, rtol = 1e-4, atol = 0)[1]
        return dens_line * compute_delta(rd, phid, r, z) * SIGMA_T * Rg
    end
    for rdp in rp_range_test
        for rp in r_range_test
            for zp in z_range_test
                for pd in phid_range_test
                    ansol = f_anl(rdp, pd, rp, zp)
                    qwsol = compute_uv_tau(grid, rdp, rp, zp, Rg)
                    ratio = compute_delta(rdp, pd, rp, zp) / sqrt((rp-rdp)^2 + zp^2)
                    qwsol = qwsol * ratio
                    @test qwsol ≈ ansol rtol = 5e-1
                end
            end
        end
    end
end
