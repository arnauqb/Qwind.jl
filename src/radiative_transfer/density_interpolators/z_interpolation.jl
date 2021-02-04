export VIGrid, get_density

struct LineDensity
    z::Vector{Float64}
    n::Vector{Float64}
    width::Vector{Float64}
end

struct VIGrid <: GridInterpolator
    grid::InterpolationGrid
    wind_kdtree::Union{KDTree,Nothing}
    line_densities::Union{Vector{LineDensity}, Nothing}
    #lines_kdtrees::Union{Array{LineKDTree,1},Nothing}
    vacuum_density::Float64
    n_timesteps::Int
    function VIGrid(
        r_range,
        z_range,
        grid;
        wind_kdtree,
        line_densities,
        vacuum_density = 1e2,
        n_timesteps = 10000,
    )
        return new(
            InterpolationGrid(r_range, z_range, grid),
            wind_kdtree,
            line_densities,
            vacuum_density,
            n_timesteps,
        )
    end
end

function VIGrid(integrators; nr = 250, nz = 250, vacuum_density = 1e2, n_timesteps = 10000)
    if integrators === nothing
        wind_kdtree = nothing
        lines_densities = nothing
    else
        wind_kdtree = create_wind_kdtree(integrators, n_timesteps, vacuum_density)
        lines_densities = create_lines_densities(wind_kdtree)
    end
    return VIGrid(
        wind_kdtree,
        lines_densities,
        nr,
        nz,
        vacuum_density = vacuum_density,
        n_timesteps = n_timesteps,
    )
end

function VIGrid(
    wind_kdtree::Union{KDTree,Nothing},
    lines_densities::Union{Array{LineDensity,1},Nothing},
    nr,
    nz;
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if wind_kdtree === nothing
        r_range = zeros(nr) #[0.0, 1.0]
        z_range = zeros(nz) #[0.0, 1.0]
        density_grid = nothing #vacuum_density .* ones((2, 2))
    else
        r_range, z_range = get_spatial_grid(wind_kdtree, nr, nz)
        density_grid = zeros(Float64, nr, nz)
        for (i, r) in enumerate(r_range)
            for (j, z) in enumerate(z_range)
                density_grid[i, j] = get_density(wind_kdtree, lines_densities, r, z)
            end
        end
    end
    return VIGrid(
        r_range,
        z_range,
        density_grid,
        wind_kdtree = wind_kdtree,
        lines_kdtrees = lines_kdtrees,
        vacuum_density = vacuum_density,
        n_timesteps = n_timesteps,
    )
end

function update_density_interpolator(interpolator::VIGrid, integrators)
    new_interpolator = VIGrid(
        integrators,
        nr = interpolator.grid.nr,
        nz = interpolator.grid.nz,
        vacuum_density = interpolator.vacuum_density,
        n_timesteps = interpolator.n_timesteps,
    )
    return new_interpolator
end

function get_density(kdtree::KDTree, lines_densities::Array{LineDensity,1}, r, z)
    rc, zc, zmax, z0, width, density, line_id, distance = get_closest_point(kdtree, r, z)
    if is_point_outside_wind(z, zmax, z0, distance, width)
        return kdtree.vacuum_density
    end
    distances = []
    densities = []
    point = [r, z]
    for (i, line_kdtree) in enumerate(lines_kdtrees)
        idx, distance = NearestNeighbors.nn(line_kdtree.line_tree, point)
        width = line_kdtree.width[idx]
        if distance > width
            continue
        end
        n = line_kdtree.n[idx]
        push!(densities, n)
        push!(distances, distance)
        #z_idx, dist2 = NearestNeighbors.nn(line_kdtree.line_tree_z, [z])
        #dens = line_kdtree.n[z_idx]
        #println("i $i")
        #println("z_idx $z_idx dens $dens")
        #println("r $(line_kdtree.r[z_idx]) z $(line_kdtree.z[z_idx])")
        #dist = Distances.evaluate(
        #    Euclidean(),
        #    point,
        #    [line_kdtree.r[z_idx], line_kdtree.z[z_idx]],
        #)
        #push!(densities, dens)
        #push!(distances, dist)
    end
    #println("densities $densities")
    #println("distances $distances")
    if length(distances) == 0
        return 1e2
    elseif length(distances) == 1
        return densities[1]
    else
        idcs = sortperm(distances)[1:2]
        d1, d2 = distances[idcs]
        n1, n2 = densities[idcs]
        println("r $r z $z")
        println("n1 $n1 n2 $n2 d1 $d1 d2 $d2")
        dens = (n1 * d2 + n2 * d1) / (d1 + d2)
        return dens
    end
end

function reduce_line(z, n, width)
    zs = [z[1]]
    ns = [[n[1]]]
    widths = [[n[1]]]
    for (zp, np, widthp) in zip(z,n, width)
        if zp > zs[end]
            push!(zs, zp)
            push!(ns, [np])
            push!(widths, [widthp])
        else
            idx = searchsorted_nearest(zs, zp)
            push!(ns[idx], np)
            push!(widths[idx], widthp)
        end
    end
    return zs, maximum.(ns), maximum.(widths)
end

function create_lines_densities(kdtree::KDTree)
    ret = []
    for line_id in unique(kdtree.line_ids)
        mask = kdtree.line_id == line_id
        z = kdtree.z[mask]
        n = kdtree.n[mask]
        width = kdtree.width[mask]
        zs, ns, widths = reduce_line(z,n, width)
        push!(ret, LineDensity(zs, ns, widths))
    end
    return ret
end
