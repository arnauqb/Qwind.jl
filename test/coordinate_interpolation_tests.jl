using Qwind, Test

@testset "Spherical to cylindrical" begin

    rho_array = 10 .^ range(log10(6), log10(10), length=500);
    theta_array = range(0, π/2, length=500);
    f(r, theta) = r * cos(theta)
    values = zeros(length(rho_array), length(theta_array));
    for (i, rho) in enumerate(rho_array)
        for (j, theta) in enumerate(theta_array)
            values[i, j] = f(rho, theta)
        end
    end
    sph_grid = SphericalGrid(rho_array, theta_array, values)
    @test sph_grid.rho_range == rho_array
    @test sph_grid.theta_range == theta_array
    cgrid = CylindricalGrid(sph_grid)

    for (i, r) in enumerate(cgrid.r_range)
        for (j, z) in enumerate(cgrid.z_range)
            rho = sqrt(r^2 + z^2)
            theta = atan(r, z)
            if (rho < rho_array[1]) || (rho > rho_array[end])
                @test cgrid.grid[i,j] == 0
                continue
            end
            @test cgrid.grid[i, j] ≈ f(rho, theta) rtol = 1e-2 atol = 1e-2
        end
    end

end


#using PyPlot
#LogNorm = matplotlib.colors.LogNorm
#fig, ax = plt.subplots()
#cm = ax.pcolormesh(cgrid.r_range, cgrid.z_range, cgrid.grid')
#plt.colorbar(cm, ax=ax)
