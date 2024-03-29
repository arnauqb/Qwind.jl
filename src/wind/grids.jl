export point_outside_grid, get_spatial_grid

function point_outside_grid(grid::InterpolationGrid, r, z)
    if r < grid.r_range[1] || r > grid.r_range[end]
        return true
    end
    if z > grid.z_range[end]
        return true
    end
    return false
end

function get_spatial_grid(r, z, r0s, nr = "auto", nz = 50; log = true)
    r_min = max(6.0, minimum(r))
    z_min = max(minimum(z), 1e-6)
    r_max = min(1e4, maximum(r))
    z_max = min(1e4, maximum(z))
    if nr == "auto"
        r_range = r0s
        additional_r = collect(range(r_range[end], r_max, step = 100)[2:end])
        r_range = vcat(r_range, additional_r)
    else
        if log
            r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
        else
            r_range = range(r_min, r_max, length = nr)
        end
    end
    if log
        z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    else
        z_range = range(z_min, z_max, length = nz - 1)
    end
    r_range = collect(r_range)
    z_range = collect(z_range)
    #r_range = unique(round.(r_range, digits = 8))
    #z_range = unique(round.(z_range, digits = 8))
    return r_range, z_range
end
