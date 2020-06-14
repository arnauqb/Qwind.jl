export Grid, out_of_grid

struct Grid
    r_min::Float64
    r_max::Float64
    z_min::Float64
    z_max::Float64
    vacuum_number_density::Float64
end

out_of_grid(params::Parameters) =
    !(
        (params.r_min <= r(params.line) <= params.r_max) &&
        (params.z_min <= z(params.line) <= params.z_max)
    )

out_of_grid(grid::Grid, line::Streamline) = 
    !(
        (grid.r_min <= r(line) <= grid.r_max) &&
        (grid.z_min <= z(line) <= grid.z_max)
    )

