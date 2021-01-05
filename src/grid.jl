export Grid, out_of_grid

struct Grid
    r_min::Float64
    r_max::Float64
    z_min::Float64
    z_max::Float64
end

function Grid(config)
    gc = config["grid"]
    return Grid(gc["r_min"], gc["r_max"], gc["z_min"], gc["z_max"])
end

function out_of_grid(grid::Grid, r, z)
    if grid.r_min <= r <= grid.r_max && grid.z_min <= z <= grid.z_max
        return false
    else
        return true
    end
end
