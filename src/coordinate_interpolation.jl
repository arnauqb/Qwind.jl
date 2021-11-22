export SphericalGrid, CylindricalGrid

using Interpolations

struct SphericalGrid
    rho_range::Vector{Float64}
    theta_range::Vector{Float64}
    grid::Matrix{Float64}
end

struct CylindricalGrid
    r_range::Vector{Float64}
    z_range::Vector{Float64}
    grid::Matrix{Float64}
end

function CylindricalGrid(spherical::SphericalGrid)
    interpolator = Interpolations.interpolate(
        (spherical.rho_range, spherical.theta_range),
        spherical.grid,
        Gridded(Linear()),
    )
    interpolator = Interpolations.extrapolate(interpolator, 0)
    r_min = max(minimum(spherical.rho_range) * sin(minimum(spherical.theta_range)), 1e-1)
    r_max = maximum(spherical.rho_range) * sin(maximum(spherical.theta_range))
    z_min = max(minimum(spherical.rho_range) * cos(maximum(spherical.theta_range)), 1e-1)
    z_max = maximum(spherical.rho_range) * cos(minimum(spherical.theta_range))
    nr = length(spherical.rho_range)
    r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nr)
    grid = zeros(length(r_range), length(z_range))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            rho = sqrt(r^2 + z^2)
            theta = atan(r, z)
            grid[i, j] = interpolator(rho, theta)
        end
    end
    return CylindricalGrid(r_range, z_range, grid)
end

