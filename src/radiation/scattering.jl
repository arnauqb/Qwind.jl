using QuadGK, ProgressMeter, LinearAlgebra, FastGaussQuadrature, Cubature

export Rectangle,
    compute_luminosity_absorbed_by_cell,
    compute_luminosity_absorbed_grid,
    compute_scattered_flux_in_cell,
    compute_total_flux_in_cell,
    compute_total_flux_grid,
    ScatteredLuminosityGrid,
    scattered_flux_in_point

struct Rectangle{T<:AbstractFloat}
    rmin::T
    rmax::T
    zmin::T
    zmax::T
end

function atan_positive(y, x)
    angle = atan(y, x)
    if angle < 0
        return angle + 2π
    else
        return angle
    end
end

function Rectangle(; rmin, rmax, zmin, zmax)
    return Rectangle(Float64(rmin), Float64(rmax), Float64(zmin), Float64(zmax))
end

function change_origin(rectangle, px, py)
    return Rectangle(
        rectangle.rmin - px,
        rectangle.rmax - px,
        rectangle.zmin - py,
        rectangle.zmax - py,
    )
end

"""
    compute_intercepting_angles(rectangle, point)

Considering all straight lines passing through the origin,
returns the lowest and highest angle of the lines that
intercept the rectangle.
"""
function compute_intercepting_angles(rectangle)
    angle1 = atan_positive(rectangle.rmin, rectangle.zmin)
    angle2 = atan_positive(rectangle.rmin, rectangle.zmax)
    angle3 = atan_positive(rectangle.rmax, rectangle.zmin)
    angle4 = atan_positive(rectangle.rmax, rectangle.zmax)
    amin = min(angle1, angle2, angle3, angle4)
    amax = max(angle1, angle2, angle3, angle4)
    return amin, amax
end
function compute_intercepting_angles(rectangle, point)
    new_rect = change_origin(rectangle, point[1], point[2])
    return compute_intercepting_angles(new_rect)
end

"""
    compute_first_intersection(rectangle, theta)

Computes the first point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function compute_first_intersection(rectangle, theta)
    cotheta = cot(theta)
    # first quadrant
    if 0 <= theta <= π / 2
        zi = cotheta * rectangle.rmin
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmin, zi]
        end
        # fourth quadrant
    elseif π / 2 <= theta <= π
        zi = cotheta * rectangle.rmin
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmin, zi]
        end
        # third quadrant
    elseif π <= theta <= 3π / 2
        zi = cotheta * rectangle.rmax
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmax, zi]
        end
        # second quadrant
    else
        zi = cotheta * rectangle.rmax
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmax, zi]
        end
    end
end
function compute_first_intersection(rectangle, theta, point)
    rectangle_centered_point = change_origin(rectangle, point[1], point[end])
    intersection = compute_first_intersection(rectangle_centered_point, theta)
    return intersection[1] + point[1], intersection[2] + point[end]
end

"""
    compute_second_intersection(rectangle, theta)

Computes the second point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function compute_second_intersection(rectangle, theta)
    cotheta = cot(theta)
    # first quadrant
    if 0 <= theta <= π / 2
        zi = cotheta * rectangle.rmax
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmax, zi]
        end
        # fourth quadrant
    elseif π / 2 <= theta <= π
        zi = cotheta * rectangle.rmax
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmax, zi]
        end
        # third quadrant
    elseif π <= theta <= 3π / 2
        zi = cotheta * rectangle.rmin
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmin, zi]
        end
        # second quadrant
    else
        zi = cotheta * rectangle.rmin
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmin, zi]
        end
    end

end
function compute_second_intersection(rectangle, theta, point)
    rectangle_centered_point = change_origin(rectangle, point[1], point[end])
    intersection = compute_second_intersection(rectangle_centered_point, theta)
    return intersection + point[1], intersection + point[end]
end

"""
    compute_intersection_size(rectangle, theta)

given a line anchored at the origin with angle theta (respect to z-axis),
computes the length of the line segment inside the given rectangle.
"""
function compute_intersection_size(rectangle, theta)
    first_intersection = compute_first_intersection(rectangle, theta)
    second_intersection = compute_second_intersection(rectangle, theta)
    return sqrt(
        (first_intersection[1] - second_intersection[1])^2 +
        (first_intersection[2] - second_intersection[2])^2,
    )
end

"""
    compute_distance_to_first_intersection(rectangle, theta)

given a line anchored at the origin with angle theta (respect to z-axis),
computes the distance from the centre to the first intersection.
"""
function compute_distance_to_first_intersection(rectangle, theta)
    first_intersection = compute_first_intersection(rectangle, theta)
    return norm(first_intersection)
end


"""
    compute_cell_optical_depth(rectangle, theta; density, Rg, mu_electron)

computes the optical depth through a cell
"""
function compute_cell_optical_depth(rectangle, theta; density, Rg, mu_electron = 1.17)
    intersection_size = compute_intersection_size(rectangle, theta)
    return intersection_size * density * mu_electron * SIGMA_T * Rg
end


"""
    compute_optical_depth_to_cell(density_grid; rectangle, theta, Rg, mu_electron)

computes the optical depth to the cell
"""
function compute_optical_depth_to_cell(
    density_grid;
    #ionization_grid;
    origin,
    rectangle,
    theta,
    Rg,
    source_luminosity,
    absorption_opacity = BoostOpacity(),
    mu_electron = 1.17,
    mu_nucleon = 0.61,
)
    first_intersection = compute_first_intersection(rectangle, theta, origin)
    return compute_optical_depth(
        density_grid.iterator,
        density_grid,
        #ionization_grid,
        absorption_opacity,
        ri = origin[1],
        phii = origin[2],
        zi = origin[3],
        rf = first_intersection[1],
        phif = 0.0,
        zf = first_intersection[2],
        Rg = Rg,
        mu_electron = mu_electron,
        mu_nucleon = mu_nucleon,
        source_luminosity = source_luminosity,
        max_tau = 30,
    )
end


"""
    compute_luminosity_absorbed_by_cell(density_grid; rectangle, Rg, mu_electron, source_luminosity)

Computes the amount of luminosity absorbed by a cell from a central source.
"""
function compute_luminosity_absorbed_by_cell(
    density_grid;
    #ionization_grid;
    source_position,
    cell,
    Rg,
    source_luminosity,
    cell_density,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    absorption_opacity = BoostOpacity(),
)
    cell_centered_source = change_origin(cell, source_position[1], source_position[2])
    theta_min, theta_max = compute_intercepting_angles(cell_centered_source)
    function f(theta)
        tau_to_cell = compute_optical_depth_to_cell(
            density_grid;
            #ionization_grid;
            origin = source_position,
            rectangle = cell,
            theta = theta,
            Rg = Rg,
            mu_electron = mu_electron,
            mu_nucleon = mu_nucleon,
            absorption_opacity = absorption_opacity,
            source_luminosity = source_luminosity,
        )
        tau_cell = compute_cell_optical_depth(
            cell,
            theta,
            density = cell_density,
            Rg = Rg,
            mu_electron = mu_electron,
        )
        return sin(theta) * exp(-tau_to_cell) * (1 - exp(-tau_cell))
    end
    integ, _ = quadgk(f, theta_min, theta_max, atol = 0, rtol = 1e-2)
    return source_luminosity / (4π) * integ
end

function compute_luminosity_absorbed_grid(
    density_grid;
    #ionization_grid;
    Rg,
    source_luminosity,
    source_position,
    mu_electron = 1.17,
    mu_nucleon = mu_nucleon,
)
    ret = zeros(length(density_grid.r_range) - 1, length(density_grid.z_range) - 1)
    function f(i)
        ret = zeros(length(density_grid.z_range) - 1)
        for j = 1:(length(density_grid.z_range) - 1)
            cell = Rectangle(
                rmin = density_grid.r_range[i],
                rmax = density_grid.r_range[i + 1],
                zmin = density_grid.z_range[j],
                zmax = density_grid.z_range[j + 1],
            )
            lumin_absorbed = compute_luminosity_absorbed_by_cell(
                density_grid;
                #ionization_grid;
                source_position = source_position,
                cell = cell,
                Rg = Rg,
                cell_density = density_grid.grid[i, j],
                source_luminosity = source_luminosity,
                mu_electron = mu_electron,
                mu_nucleon = mu_nucleon,
            )
            ret[j] = lumin_absorbed
        end
        return ret
    end
    rets = nothing
    if is_logging(stderr)
        rets = pmap(f, 1:(length(density_grid.r_range) - 1), batch_size = 10)
    else
        rets = @showprogress pmap(f, 1:(length(density_grid.r_range) - 1), batch_size = 10)
    end
    rets = reduce(hcat, rets)'
    return rets
end



function compute_scattered_flux_in_cell(
    density_grid;
    #ionization_grid;
    scattered_luminosity_per_cell,
    cell,
    cell_density,
    absorption_opacity,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    Rg,
    nodes,
    weights,
)
    density_grid = density_grid
    ret = 0.0
    r_target = (cell.rmin + cell.rmax) / 2
    z_target = (cell.zmin + cell.zmax) / 2
    for i = 1:(length(density_grid.r_range) - 1)
        for j = 1:(length(density_grid.z_range) - 1)
            r_source = (density_grid.r_range[i + 1] + density_grid.r_range[i]) / 2
            z_source = (density_grid.z_range[j + 1] + density_grid.z_range[j]) / 2
            source_luminosity = scattered_luminosity_per_cell[i, j]
            cell_density = density_grid.grid[i, j]
            function f(phi)
                phi = π * (phi + 1) / 2
                distance =
                    compute_distance_cylindrical(
                        r_source,
                        phi,
                        z_source,
                        r_target,
                        0,
                        z_target,
                    ) * Rg
                distance = max(distance, 1e-6)
                tau = compute_optical_depth(
                    density_grid.iterator,
                    density_grid,
                    #ionization_grid,
                    absorption_opacity,
                    ri = r_source,
                    phii = phi,
                    zi = z_source,
                    rf = r_target,
                    phif = 0,
                    zf = z_target,
                    Rg = Rg,
                    mu_electron = mu_electron,
                    mu_nucleon = mu_nucleon,
                    max_tau = 20,
                    source_luminosity = source_luminosity,
                )
                ret = source_luminosity * exp(-tau) / distance^2
                return ret
            end
            flux = π * weights' * f.(nodes)
            #flux, _ = 2 .* quadgk(f, 0, π, atol = 0, rtol = 1e-1)
            ret += flux
        end
    end
    return ret
end

function compute_ionizing_flux_from_center_in_cell(
    density_grid,
    ionization_grid;
    cell,
    absorption_opacity,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    source_luminosity,
    z_xray,
    Rg,
)
    rf = (cell.rmin + cell.rmax) / 2
    zf = (cell.zmin + cell.zmax) / 2
    distance_center_sq = (rf^2 + zf^2) * Rg^2
    tau_center = compute_optical_depth(
        density_grid.iterator,
        density_grid,
        ionization_grid,
        absorption_opacity,
        ri = 0,
        phii = 0,
        zi = z_xray,
        rf = rf,
        phif = 0,
        zf = zf,
        Rg = Rg,
        mu_electron = mu_electron,
        mu_nucleon = mu_nucleon,
        max_tau = 30,
    )
    from_center = source_luminosity * exp(-tau_center) / distance_center_sq
    return from_center
end

function compute_total_flux_in_cell(
    density_grid,
    ionization_grid;
    scattered_luminosity_per_cell,
    source_position,
    source_luminosity,
    cell,
    cell_density,
    absorption_opacity,
    z_xray,
    Rg,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    include_scattered = true,
    nodes,
    weights,
)
    ret = compute_ionizing_flux_from_center_in_cell(
        density_grid,
        ionization_grid,
        cell = cell,
        absorption_opacity = absorption_opacity,
        mu_electron = mu_electron,
        mu_nucleon = mu_nucleon,
        source_luminosity = source_luminosity,
        Rg = Rg,
        z_xray = z_xray,
    )
    if include_scattered
        ret += compute_scattered_flux_in_cell(
            density_grid,
            ionization_grid,
            scattered_luminosity_per_cell = scattered_luminosity_per_cell,
            cell = cell,
            cell_density = cell_density,
            absorption_opacity = absorption_opacity,
            Rg = Rg,
            nodes = nodes,
            weights = weights,
        )
    end
    return ret
end

function compute_total_flux_from_center_grid(
    density_grid,
    ionization_grid;
    Rg,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    source_luminosity,
)
    ret = zeros(length(density_grid.r_range) - 1, length(density_grid.z_range) - 1)
    function f(i)
        ret = zeros(length(density_grid.z_range) - 1)
        for j = 1:(length(density_grid.z_range) - 1)
            cell = Rectangle(
                density_grid.r_range[i],
                density_grid.r_range[i + 1],
                density_grid.z_range[j],
                density_grid.z_range[j + 1],
            )
            ret[j] = compute_ionizing_flux_from_center_in_cell(
                density_grid,
                ionization_grid,
                cell = cell,
                cell_density = density_grid.grid[i, j],
                Rg = Rg,
                mu_electron = mu_electron,
                mu_nucleon = mu_nucleon,
                source_luminosity = source_luminosity,
            )
        end
        return ret
    end
    rets = nothing
    if is_logging(stderr)
        rets = pmap(f, 1:(length(density_grid.r_range) - 1), batch_size = 10)
    else
        rets = @showprogress pmap(f, 1:(length(density_grid.r_range) - 1), batch_size = 10)
    end
    rets = reduce(hcat, rets)'
    return rets
end

#function compute_total_flux_grid(
#    density_grid,
#    ionization_grid;
#    Rg,
#    mu_nucleon = 0.61,
#    mu_electron = 1.17,
#    xray_luminosity,
#    parallel = true,
#    include_scattered = true,
#    z_xray,
#    absorption_opacity,
#    n_nodes = 10,
#)
#    ret = zeros(length(density_grid.r_range) - 1, length(density_grid.z_range) - 1)
#    absorbed_from_center = compute_luminosity_absorbed_grid(
#        density_grid,
#        ionization_grid,
#        Rg = Rg,
#        source_position = [0, 0, z_xray],
#        source_luminosity = xray_luminosity,
#        mu_electron = mu_electron,
#        mu_nucleon = mu_nucleon,
#    )
#    nodes, weights = gausslegendre(n_nodes)
#    function f(i)
#        ret = zeros(length(density_grid.z_range) - 1)
#        for j = 1:(length(density_grid.z_range) - 1)
#            cell = Rectangle(
#                density_grid.r_range[i],
#                density_grid.r_range[i + 1],
#                density_grid.z_range[j],
#                density_grid.z_range[j + 1],
#            )
#            ret[j] = compute_total_flux_in_cell(
#                density_grid,
#                ionization_grid,
#                scattered_luminosity_per_cell = absorbed_from_center,
#                cell = cell,
#                cell_density = density_grid.grid[i, j],
#                Rg = Rg,
#                mu_electron = mu_electron,
#                source_position = [0, 0, z_xray],
#                source_luminosity = xray_luminosity,
#                include_scattered = include_scattered,
#                absorption_opacity = absorption_opacity,
#                z_xray = z_xray,
#                nodes = nodes,
#                weights = weights,
#            )
#        end
#        return ret
#    end
#    rets = nothing
#    if parallel
#        if is_logging(stderr)
#            rets =
#                @showprogress pmap(f, 1:(length(density_grid.r_range) - 1), batch_size = 10)
#        else
#            rets = pmap(f, 1:(length(density_grid.r_range) - 1), batch_size = 10)
#        end
#    else
#        rets = f.(1:(length(density_grid.r_range) - 1))
#    end
#    rets = reduce(hcat, rets)'
#    return rets
#end

compute_total_flux_grid(radiation; parallel = true) = compute_total_flux_grid(
    radiation.wi.density_grid,
    radiation.wi.ionization_grid,
    Rg = radiation.bh.Rg,
    mu_nucleon = radiation.mu_nucleon,
    mu_electron = radiation.mu_electron,
    source_luminosity = radiation.xray_luminosity,
    include_scattered = radiation.xray_scattering,
    absorption_opacity = radiation.xray_opacity,
    parallel = parallel,
)

struct ScatteredLuminosityGrid{T,U,V} <: InterpolationGrid
    r_range::Vector{T}
    z_range::Vector{T}
    grid::Array{T,2}
    nr::Union{U,String}
    nz::U
    iterator::GridIterator{T,U,V}
    interpolator::Any
    function ScatteredLuminosityGrid(r_range, z_range, grid, nr = nothing, nz = nothing)
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        iterator = GridIterator(r_range, z_range)
        interpolator =
            Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
        interpolator = Interpolations.extrapolate(interpolator, 1e2)
        return new{typeof(r_range[1]),Int,Bool}(
            r_range,
            z_range,
            grid,
            nr,
            nz,
            iterator,
            interpolator,
        )
    end
end

function ScatteredLuminosityGrid(nr::Union{String,Int}, nz::Int, value::Number)
    r_range = [-1.0, 0.0]
    z_range = [-1.0, 0.0]
    density_grid = value .* [[1.0, 1.0] [1.0, 1.0]]
    return ScatteredLuminosityGrid(r_range, z_range, density_grid, nr, nz)
end

function ScatteredLuminosityGrid(
    density_grid::DensityGrid;
    Rg,
    source_luminosity,
    source_position,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
)
    grid = compute_luminosity_absorbed_grid(
        density_grid,
        Rg = Rg,
        source_luminosity = source_luminosity,
        source_position = source_position,
        mu_electron = mu_electron,
        mu_nucleon = mu_nucleon,
    )
    rr = density_grid.r_range[1:(end - 1)] + diff(density_grid.r_range) / 2
    zz = density_grid.z_range[1:(end - 1)] + diff(density_grid.z_range) / 2
    return ScatteredLuminosityGrid(rr, zz, grid, density_grid.nr, density_grid.nz)
end

ScatteredLuminosityGrid(density_grid, rad) = ScatteredLuminosityGrid(
    density_grid,
    Rg = rad.bh.Rg,
    source_luminosity = rad.xray_luminosity,
    source_position = (0, 0, rad.z_xray),
    mu_electron = rad.mu_electron,
    mu_nucleon = rad.mu_nucleon,
)

function interpolate_scattered_luminosity(grid::ScatteredLuminosityGrid, r, z)
    return grid.interpolator(r, z)
end

function get_scattered_luminosity(grid::ScatteredLuminosityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return 0.0
    end
    ridx = searchsorted_nearest(grid.r_range, r)
    zidx = searchsorted_nearest(grid.z_range, z)
    return grid.grid[ridx, zidx]
end

function scattered_flux_in_point_integrand!(
    v,
    scattered_luminosity_grid,
    density_grid;
    r_source,
    z_source,
    phi_source,
    r_target,
    z_target,
    Rg,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    absorption_opacity,
)
    source_luminosity =
        interpolate_scattered_luminosity(scattered_luminosity_grid, r_source, z_source)
    distance =
        compute_distance_cylindrical(
            r_source,
            phi_source,
            z_source,
            r_target,
            0,
            z_target,
        ) * Rg
    distance = max(distance, 1e-6)
    tau = compute_optical_depth(
        density_grid.iterator,
        density_grid,
        absorption_opacity,
        ri = r_source,
        phii = phi_source,
        zi = z_source,
        rf = r_target,
        phif = 0,
        zf = z_target,
        Rg = Rg,
        mu_electron = mu_electron,
        mu_nucleon = mu_nucleon,
        source_luminosity = source_luminosity,
        max_tau = 20,
    )
    v[1] = source_luminosity * exp(-tau) / distance^2
end

function scattered_flux_in_point(
    scattered_luminosity_grid,
    density_grid;
    r,
    z,
    absorption_opacity,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    Rg,
    #nodes,
    #weights,
    rtol = 1e-3,
    atol = 0,
    norm = Cubature.INDIVIDUAL,
    maxevals = 50000,
)
    f(x, v) = scattered_flux_in_point_integrand!(
        v,
        scattered_luminosity_grid,
        density_grid,
        r_source = x[1],
        phi_source = x[2],
        z_source = x[3],
        r_target = r,
        z_target = z,
        absorption_opacity = absorption_opacity,
        mu_electron = mu_electron,
        mu_nucleon = mu_nucleon,
        Rg = Rg,
    )
    integrated, error = hcubature(
        1,
        f,
        (density_grid.r_range[1], 0, density_grid.z_range[1]),
        (density_grid.r_range[end], π, density_grid.z_range[end]),
        abstol = atol,
        reltol = rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
    return integrated[1]
end
