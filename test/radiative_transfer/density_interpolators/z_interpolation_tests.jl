import Distances
using QuadGK

function make_lines_kdtrees(
    density_func;
    r_norm = 1000,
    z_norm = 100,
    n_timesteps = 1000,
    n_lines = 10,
)
    ret = LineKDTree[]
    dr = 1000 / n_lines
    for i = 1:n_lines
        r_range = i * dr .* ones(n_timesteps)
        z_range = range(0, 100, length = n_timesteps)
        n_range = density_func.(r_range, z_range)
        width_range = dr .* ones(n_timesteps)
        zmax_range = 2000 .* ones(n_timesteps)
        z0_range = zeros(n_timesteps)
        lkt =
            LineKDTree(r_range, z_range, n_range, width_range, n_timesteps)
        push!(ret, lkt)
    end
    return ret
end

@testset "Test z interpolation" begin

    #@testset "Test reduce line" begin
    #    r = [1, 2, 3, 4, 1]
    #    z = [2, 4, 6, 3, 1]
    #    n = [1, 2, 3, 4, 5]
    #    widths = [5, 5, 5, 5, 5]
    #    rr, zr, nr, wr = Qwind.reduce_line(r, z, n, widths)
    #    @test rr == [1, 2, 3]
    #    @test zr == [2, 4, 6]
    #    @test nr == [1, 2, 3]
    #    @test wr == [5, 5, 5]
    #end

    #@testset "Test NEuclidean" begin
    #    norm = Qwind.NEuclidean(2, 4)
    #    @test Distances.evaluate(norm, [1, 2], [3, 4]) ≈ sqrt(1.25)
    #end

    #density_func = (r, z) -> 1e8 
    density_func = (r, z) -> 1e8 ./ (r^2 + z^2);
    lkdt = make_lines_kdtrees(density_func, n_lines = 1000, r_norm = 1000, z_norm = 100);

    @testset "Test lines kdtrees" begin
        for r in range(1, 1000, length = 50)
            for z in range(0, 100, length = 50)
                true_density = density_func(r, z)
                got_density = get_density(lkdt, r, z)
                @test true_density ≈ got_density rtol = 0.05
            end
        end
    end

    function tau_kernel(s, density_func, ri, zi, rf, zf)
        r = ri + s * (rf - ri)
        z = zi + s * (zf - zi)
        return density_func(r, z)
    end

    function compute_analytic_tau(density_func, ri, zi, rf, zf, Rg)
        integ, err = quadgk(s -> tau_kernel(s, density_func, ri, zi, rf, zf), 0, 1)
        arclen = sqrt((rf - ri)^2 + (zf - zi)^2)
        return integ * SIGMA_T * Rg * arclen
    end

    vi_grid = VIGrid(lkdt, 50, 50);

    @testset "Test density grid" begin
        bh = create_test_bh()
        Rg = bh.Rg
        zi = 0
        for ri in range(6, 900, length = 25)
            for rf in range(6, 900, length = 25)
                for zf in range(0, 100, length = 25)
                    tau_expected = compute_analytic_tau(density_func, ri, zi, rf, zf, Rg)
                    tau_grid = compute_uv_tau(vi_grid, ri, zi, rf, zf, Rg)
                    tau2_grid_highxr = compute_xray_tau(vi_grid, ri, zi, rf, zf, 1e80, Rg)
                    tau2_grid_lowxr = compute_xray_tau(vi_grid, ri, zi, rf, zf, 1, Rg)
                    #if ! (isapprox(tau2_grid, tau_expected, rtol=0.25))
                    #    println("ri $ri rf $rf zf $zf")
                    #    println("tau2 $tau2_grid exp $tau_expected")
                    #end
                    @test tau_grid ≈ tau_expected rtol=0.5
                    @test tau2_grid_highxr ≈ tau_expected rtol=0.5
                    @test tau2_grid_lowxr ≈ 100 * tau_expected rtol=0.5
                end
            end
        end
    end

end
