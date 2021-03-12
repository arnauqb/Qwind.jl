export point_outside_grid

function point_outside_grid(grid::InterpolationGrid, r, z)
    if r < grid.r_range[1] || r > grid.r_range[end]
        return true
    end
    if z > grid.z_range[end]
        return true
    end
    return false
end
