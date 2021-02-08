export VIGrid, get_density
import Distances

struct NEuclidean <: Metric
    r_norm::Float64
    z_norm::Float64
end

#(norm::NEuclidean)(a, b) = Distances.evaluate(
#    Euclidean(),
#    [a[1] / norm.r_norm, a[2] / norm.z_norm],
#    [b[1] / norm.r_norm, b[2] / norm.z_norm],
#)
(norm::NEuclidean)(a, b) =
    sqrt(((a[1] - b[1]) / norm.r_norm)^2 + ((a[2] - b[2]) / norm.z_norm)^2)

struct LineKDTree <: NNInterpolator
    r::Vector{Float64}
    z::Vector{Float64}
    n::Vector{Float64}
    width::Vector{Float64}
    r_norm::Float64
    z_norm::Float64
    zmax::Float64
    z0::Float64
    line_tree::NearestNeighbors.BallTree
    n_timesteps::Int
    function LineKDTree(r, z, density, width, r_norm, z_norm, n_timesteps)
        points_line = hcat(r, z)'
        points_line = convert(Array{Float64,2}, points_line)
        line_tree = NearestNeighbors.BallTree(points_line, NEuclidean(r_norm, z_norm))
        return new(
            r,
            z,
            density,
            width,
            r_norm,
            z_norm,
            maximum(z),
            minimum(z),
            line_tree,
            n_timesteps,
        )
    end
end

function reduce_line(r, z, n, width)
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    widths = [width[1]]
    for (rp, zp, np, wp) in zip(r, z, n, width)
        if zp > zs[end]
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
            push!(widths, wp)
        end
    end
    return rs, zs, ns, widths
end

function create_lines_kdtrees(
    integrators::Vector{Sundials.IDAIntegrator};
    n_timesteps = 10000,
)
    ret = LineKDTree[]
    r_norm = maximum(maximum([integrator.p.data[:r] for integrator in integrators]))
    z_norm = maximum(maximum([integrator.p.data[:z] for integrator in integrators]))
    for integrator in integrators
        r, z, zmax, z0, width, density =
            get_dense_solution_from_integrator(integrator, n_timesteps)
        rs, zs, ns, widths = reduce_line(r, z, density, width)
        line_kdtree = LineKDTree(rs, zs, ns, widths, r_norm, z_norm, n_timesteps)
        push!(ret, line_kdtree)
    end
    return ret
end


function get_spatial_grid(lines_kdtrees::Vector{LineKDTree}, nr, nz)
    r_min = Inf
    r_max = 0
    z_min = Inf
    z_max = 0
    max_width = 0
    for lk in lines_kdtrees
        r_min = min(r_min, minimum(lk.r))
        r_max = max(r_max, maximum(lk.r))
        z_min = min(z_min, minimum(lk.z))
        z_max = max(z_max, maximum(lk.z))
        max_width = max(max_width, maximum(lk.width))
    end
    r_min = max(6, r_min)
    r_max += max_width
    z_min = max(1e-6, z_min)
    z_max += max_width
    r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    z_range = pushfirst!(z_range, 0.0)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    return r_range, z_range
end

struct VIGrid <: GridInterpolator
    grid::InterpolationGrid
    lines_kdtrees::Union{Array{LineKDTree,1},Nothing}
    vacuum_density::Float64
    n_timesteps::Int
    function VIGrid(
        r_range::Vector{Float64},
        z_range::Vector{Float64},
        grid::Union{Array{Float64,2},Nothing};
        lines_kdtrees,
        vacuum_density = 1e2,
        n_timesteps = 10000,
    )
        return new(
            InterpolationGrid(r_range, z_range, grid),
            lines_kdtrees,
            vacuum_density,
            n_timesteps,
        )
    end
end


function VIGrid(
    integrators;
    nr = 250,
    nz = 250,
    n = nothing,
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if n !== nothing
        nr = n
        nz = n
    end
    nr = Int(nr)
    nz = Int(nz)
    if integrators === nothing
        lines_kdtrees = nothing
    else
        lines_kdtrees = create_lines_kdtrees(integrators, n_timesteps = n_timesteps)
    end
    return VIGrid(
        lines_kdtrees,
        nr,
        nz,
        vacuum_density = vacuum_density,
        n_timesteps = n_timesteps,
    )
end


function VIGrid(
    lines_kdtrees::Union{Array{LineKDTree,1},Nothing},
    nr,
    nz;
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if lines_kdtrees === nothing
        r_range = zeros(nr) #[0.0, 1.0]
        z_range = zeros(nz) #[0.0, 1.0]
        density_grid = nothing #vacuum_density .* ones((2, 2))
    else
        r_range, z_range = get_spatial_grid(lines_kdtrees, nr, nz)
        density_grid = zeros(Float64, nr, nz)
        for (i, r) in enumerate(r_range)
            for (j, z) in enumerate(z_range)
                density_grid[i, j] = get_density(lines_kdtrees, r, z)
            end
        end
    end
    return VIGrid(
        r_range,
        z_range,
        density_grid,
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

function get_closest_point(line_kdtree::LineKDTree, point)
    idx, _ = NearestNeighbors.nn(line_kdtree.line_tree, point)
    n = line_kdtree.n[idx]
    width = line_kdtree.width[idx]
    r = line_kdtree.r[idx]
    z = line_kdtree.z[idx]
    distance = Distances.evaluate(Euclidean(), point, [r, z])
    return n, width, distance, idx
end

function get_density(lines_kdtrees::Array{LineKDTree,1}, r, z)
    distances = []
    densities = []
    point = [r, z]
    for (i, line_kdtree) in enumerate(lines_kdtrees)
        if z > line_kdtree.zmax || z < line_kdtree.z0
            continue
        end
        n, width, distance, idx = get_closest_point(line_kdtree, point)
        if distance > width
            continue
        end
        push!(densities, n)
        push!(distances, distance)
    end
    if length(distances) == 0
        return 1e2
    elseif length(distances) == 1
        return densities[1]
    else
        idcs = sortperm(distances)[1:2]
        d1, d2 = distances[idcs]
        n1, n2 = densities[idcs]
        dens = (n1 * d2 + n2 * d1) / (d1 + d2)
        return dens
    end
end

