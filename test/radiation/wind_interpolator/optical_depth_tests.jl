using Qwind
using Test
using QuadGK
using Roots

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

    function compute_delta(rd, phid, r, z)
        return sqrt(r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    end
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    Rg = bh.Rg

    @testset "Constant density" begin
        r_range = range(0, 1000.0, length = 11)
        z_range = range(0, 1000.0, length = 11)
        rp_range_test = range(0.0, 1000, length = 10)
        r_range_test = range(0.0, 1000, length = 10)
        z_range_test = range(0.0, 1000, length = 10)
        phid_range_test = range(0, π, length = 10)
        density_grid = 2e8 .* ones((length(r_range), length(z_range)))
        grid = DensityGrid(r_range, z_range, density_grid)
        mu_electron = 2
        f_anl(rd, phid, r, z) =
            compute_delta(rd, phid, r, z) * 2e8 * SIGMA_T * Rg * mu_electron
        for rdp in rp_range_test
            for rp in r_range_test
                for zp in z_range_test
                    for pd in phid_range_test
                        #println("rdp $rdp rp $rp zp $zp pd $pd")
                        truesol = f_anl(rdp, pd, rp, zp)
                        qwsol = compute_optical_depth(
                            grid,
                            ri = rdp,
                            phii = pd,
                            zi = 0.0,
                            rf = rp,
                            phif = 0.0,
                            zf = zp,
                            Rg = Rg,
                            max_tau = Inf,
                            mu_electron = mu_electron,
                        )
                        @test qwsol ≈ truesol rtol = 1e-6
                    end
                end
            end
        end
    end

    @testset "Linear density" begin
        r_range = range(0, 100.0, length = 500)
        z_range = range(0, 100.0, length = 500)
        rp_range_test = range(10, 90, length = 10)
        r_range_test = range(10, 90, length = 10)
        z_range_test = range(10, 90, length = 10)
        phid_range_test = range(0, π - 0.1, length = 10)
        density_grid = zeros((length(r_range), length(z_range)))
        mu_electron = 2
        dens_func(r, z) = 1e8 * (2 * r + z)
        for (i, rp) in enumerate(r_range)
            for (j, zp) in enumerate(z_range)
                density_grid[i, j] = dens_func(rp, zp)
            end
        end
        grid = DensityGrid(r_range, z_range, density_grid)
        function f_anl_kernel(t, rd, phid, r, z)
            x = rd * cos(phid) * (1 - t) + r * t
            y = rd * sin(phid) * (1 - t)
            zp = z * t
            rp = sqrt(x^2 + y^2)
            dens = dens_func(rp, zp)
            return dens
        end
        function f_anl(rd, phid, r, z)
            dens_line = quadgk(
                t -> f_anl_kernel(t, rd, phid, r, z),
                0,
                1,
                rtol = 1e-2,
                atol = 1e-8,
            )[1]
            return dens_line * mu_electron * compute_delta(rd, phid, r, z) * SIGMA_T * Rg
        end
        for rdp in rp_range_test
            for rp in r_range_test
                for zp in z_range_test
                    for pd in phid_range_test
                        #println("rdp $rdp rp $rp zp $zp pd $pd")
                        delta = compute_delta(rdp, pd, rp, zp)
                        #println("delta $delta")
                        truesol = f_anl(rdp, pd, rp, zp)
                        qwsol = compute_optical_depth(
                            grid,
                            ri = rdp,
                            phii = pd,
                            zi = 0.0,
                            rf = rp,
                            phif = 0.0,
                            zf = zp,
                            Rg = Rg,
                            max_tau = Inf,
                            mu_electron = mu_electron,
                        )
                        @test qwsol ≈ truesol rtol = 2e-2 atol = 1e-3
                    end
                end
            end
        end
    end

    @testset "Exponential density" begin
        r_range = range(0, 100.0, length = 500)
        z_range = range(0, 100.0, length = 500)
        rp_range_test = range(10, 90, length = 10)
        r_range_test = range(10, 90, length = 10)
        z_range_test = range(10, 90, length = 10)
        phid_range_test = range(0, π - 0.1, length = 10)
        density_grid = zeros((length(r_range), length(z_range)))
        mu_electron = 2
        dens_func(r, z) = 1e8 * exp(-r - z^2)
        for (i, rp) in enumerate(r_range)
            for (j, zp) in enumerate(z_range)
                density_grid[i, j] = dens_func(rp, zp)
            end
        end
        grid = DensityGrid(r_range, z_range, density_grid)
        function f_anl_kernel(t, rd, phid, r, z)
            x = rd * cos(phid) * (1 - t) + r * t
            y = rd * sin(phid) * (1 - t)
            zp = z * t
            rp = sqrt(x^2 + y^2)
            dens = dens_func(rp, zp)
            return dens
        end
        function f_anl(rd, phid, r, z)
            dens_line =
                quadgk(t -> f_anl_kernel(t, rd, phid, r, z), 0, 1, rtol = 1e-4, atol = 0)[1]
            return dens_line * mu_electron * compute_delta(rd, phid, r, z) * SIGMA_T * Rg
        end
        for rdp in rp_range_test
            for rp in r_range_test
                for zp in z_range_test
                    for pd in phid_range_test
                        truesol = f_anl(rdp, pd, rp, zp)
                        qwsol = compute_optical_depth(
                            grid,
                            ri = rdp,
                            phii = pd,
                            zi = 0,
                            rf = rp,
                            phif = 0,
                            zf = zp,
                            Rg = Rg,
                            max_tau = Inf,
                            mu_electron = mu_electron,
                        )
                        @test qwsol ≈ truesol rtol = 2e-2 atol = 1e-3
                    end
                end
            end
        end
    end

end

@testset "Test X-ray optical depth" begin

    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    Rg = bh.Rg
    xl = 1e45

    @testset "Constant density" begin

        r_range = range(0, 1000.0, length = 50)
        z_range = range(0, 1000.0, length = 50)
        r_range_test = range(0.0, 1000, length = 50)
        z_range_test = range(0.0, 1000, length = 50)
        density_grid = 2e8 .* ones((length(r_range), length(z_range)))
        grid = DensityGrid(r_range, z_range, density_grid)
        mu_electron = 2
        mu_nucleon = 5
        @testset "Thomson only" begin
            f_anl(r, z) = sqrt(r^2 + z^2) * 2e8 * SIGMA_T * Rg * mu_electron
            for rp in r_range_test
                for zp in z_range_test
                    truesol = f_anl(rp, zp)
                    qwsol = compute_tau_xray(
                        grid,
                        Thomson(),
                        ri = 0.0,
                        zi = 0.0,
                        rf = rp,
                        zf = zp,
                        xray_luminosity = xl,
                        Rg = Rg,
                        mu_nucleon = mu_nucleon,
                        mu_electron = mu_electron,
                    )
                    @test truesol ≈ qwsol rtol = 1e-3
                end
            end
        end

        @testset "Boost" begin
            function get_ion_radius(density)
                f(r) =
                    10^5 -
                    xl / (density / mu_nucleon * (r * Rg)^2) *
                    exp(-r * density * mu_electron * SIGMA_T * Rg)
                zero = find_zero(f, (10, 1000))
                return zero
            end
            function f_anl(r, z, rx)
                d = sqrt(r^2 + z^2)
                srg = SIGMA_T * Rg * 2e8 * mu_electron
                if d <= rx
                    return d * srg
                else
                    return rx * srg + 100 * (d - rx) * srg
                end
            end
            rx = get_ion_radius(2e8)
            for rp in r_range_test
                for zp in z_range_test
                    truesol = f_anl(rp, zp, rx)
                    qwsol = compute_tau_xray(
                        grid,
                        Boost(),
                        ri = 0.0,
                        zi = 0.0,
                        rf = rp,
                        zf = zp,
                        xray_luminosity = xl,
                        Rg = Rg,
                        mu_electron = mu_electron,
                        mu_nucleon = mu_nucleon,
                    )
                    @test truesol ≈ qwsol rtol = 1e-3
                end
            end
        end

    end
end

@testset "Test from center" begin
    function compute_delta(rd, phid, r, z)
        return sqrt(r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    end
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    Rg = bh.Rg
    r_range = range(0, 100.0, length = 500)
    z_range = range(0, 100.0, length = 500)
    r_range_test = range(10, 90, length = 100)
    z_range_test = range(10, 90, length = 100)
    density_grid = zeros((length(r_range), length(z_range)))
    mu_electron = 2
    dens_func(r, z) = 1e8 * exp(-r - z^2)
    for (i, rp) in enumerate(r_range)
        for (j, zp) in enumerate(z_range)
            density_grid[i, j] = dens_func(rp, zp)
        end
    end
    grid = DensityGrid(r_range, z_range, density_grid)
    function f_anl_kernel(t, rd, phid, r, z)
        x = rd * cos(phid) * (1 - t) + r * t
        y = rd * sin(phid) * (1 - t)
        zp = z * t
        rp = sqrt(x^2 + y^2)
        dens = dens_func(rp, zp)
        return dens
    end
    function f_anl(rd, phid, r, z)
        dens_line =
            quadgk(t -> f_anl_kernel(t, rd, phid, r, z), 0, 1, rtol = 1e-4, atol = 0)[1]
        return dens_line * mu_electron * compute_delta(rd, phid, r, z) * SIGMA_T * Rg
    end
    for rp in r_range_test
        for zp in z_range_test
            truesol = f_anl(0.0, 0.0, rp, zp)
            qwsol = compute_optical_depth(
                grid,
                ri = 0.0,
                phii = 0.0,
                zi = 0,
                rf = rp,
                phif = 0,
                zf = zp,
                Rg = Rg,
                max_tau = Inf,
                mu_electron = mu_electron,
            )
            @test qwsol ≈ truesol rtol = 2e-2 atol = 1e-3
        end
    end
end

