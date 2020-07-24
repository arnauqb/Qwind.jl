export Grid, out_of_grid

struct Grid
    rmin::Float64
    rmax::Float64
    zmin::Float64
    zmax::Float64
end

function out_of_grid(grid::Grid, r, z)
    if grid.rmin <= r <= grid.rmax && grid.zmin <= z <= grid.zmax
        return false
    else
        return true
    end
end
