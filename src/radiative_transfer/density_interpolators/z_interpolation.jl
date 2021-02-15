export VIGrid, get_density
import Distances, Interpolations

struct NEuclidean <: Metric
    r_norm::Float64
    z_norm::Float64
end

(norm::NEuclidean)(a, b) =
    sqrt(((a[1] - b[1]) / norm.r_norm)^2 + ((a[2] - b[2]) / norm.z_norm)^2)
#(norm::NEuclidean)(a, b) = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2)

struct LineKDTree <: NNInterpolator
    r::Vector{Float64}
    z::Vector{Float64}
    n::Vector{Float64}
    width::Vector{Float64}
    #r_norm::Float64
    #z_norm::Float64
    zmax::Float64
    r0::Float64
    z0::Float64
    line_tree::NearestNeighbors.NNTree
    n_timesteps::Int
    function LineKDTree(r, z, density, width, n_timesteps)
        r_pos, z_pos, width_pos = reduce_line_uniformely_spaced_distance(r, z, width)
        points_line = hcat(r_pos, z_pos)'
        points_line = convert(Array{Float64,2}, points_line)
        line_tree = NearestNeighbors.KDTree(points_line)
        r_den, z_den, n_den = reduce_line_uniformely_spaced_density(r, z, density)
        #line_tree = NearestNeighbors.BallTree(points_line, NEuclidean(r_norm, z_norm))
        #line_tree = NearestNeighbors.BallTree(points_line, NEuclidean(maximum(r), maximum(z)))
        return new(
            r_den,
            z_den,
            n_den,
            width_pos,
            maximum(z),
            r[1] - width[1]/2,
            minimum(z),
            line_tree,
            n_timesteps,
        )
    end
end

function reduce_line_uniformely_spaced_distance(r, z, width)
    rs = [r[1]]
    zs = [z[1]]
    #ns = [n[1]]
    widths = [width[1]]
    #for (rp, zp, np, wp) in zip(r, z, n, width)
    for (rp, zp, wp) in zip(r, z, width)
        #ns_eps = abs(np - ns[end]) / ns[end]
        d0 = d_euclidean(rp, rs[1], zp, zs[1])
        d = d_euclidean(rp, rs[end], zp, zs[end])
        d_lim = 0.01 * d0
        if zp > zs[end] && d > d_lim  # && np > 1e3 #ns_eps > 0.01
            push!(rs, rp)
            push!(zs, zp)
            #push!(ns, np)
            push!(widths, wp)
        end
    end
    return rs, zs, widths #, ns, widths
end

function reduce_line_uniformely_spaced_density(r, z, n)
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    #widths = [width[1]]
    #for (rp, zp, np, wp) in zip(r, z, n, width)
    for (rp, zp, np) in zip(r, z, n)
        ns_eps = abs(np - ns[end]) / ns[end]
        if zp > zs[end] && ns_eps > 0.01# && np > 1e3
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
        end
    end
    return rs, zs, ns
end

function create_lines_kdtrees(
    integrators::Vector{Sundials.IDAIntegrator};
    n_timesteps = 10000,
)
    ret = LineKDTree[]
    r_norm = maximum(maximum([integrator.p.data[:r] for integrator in integrators]))
    z_norm = maximum(maximum([integrator.p.data[:z] for integrator in integrators]))
    for integrator in integrators
        r, z, zmax, z0, width, n =
            get_dense_solution_from_integrator(integrator, n_timesteps)
        #rs, zs, ns, widths = reduce_line(r, z, density, width)
        #println(length(rs))
        #line_kdtree = LineKDTree(rs, zs, ns, widths, r_norm, z_norm, n_timesteps)
        line_kdtree = LineKDTree(r, z, n, width, n_timesteps)
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
        r_min = min(r_min, lk.r0)
        r_max = max(r_max, lk.r[1] + lk.width[1])
        z_min = min(z_min, lk.z0)
        z_max = max(z_max, lk.zmax)
        max_width = max(max_width, maximum(lk.width))
    end
    r_min = max(6, r_min)
    r_max += max_width
    z_min = max(1e-3, z_min)
    z_max += max_width
    r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    z_range = pushfirst!(z_range, 0.0)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    return r_range, z_range
end

struct VIGrid <: GridInterpolator
    grid::InterpolationGrid
    interpolator
    #lines_kdtrees::Union{Array{LineKDTree,1},Nothing}
    vacuum_density::Float64
    n_timesteps::Int
    function VIGrid(
        r_range::Vector{Float64},
        z_range::Vector{Float64},
        grid::Union{Array{Float64,2},Nothing};
        #lines_kdtrees,
        vacuum_density = 1e2,
        n_timesteps = 10000,
    )
        if grid === nothing
            itp = (r,z) -> vacuum_density
        else
            itp = Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
            itp = Interpolations.extrapolate(itp, 1e2)
        end
        return new(
            InterpolationGrid(r_range, z_range, grid),
            itp,
            #lines_kdtrees,
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
        r_range = zeros(nr)
        z_range = zeros(nz)
        density_grid = nothing
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
        #lines_kdtrees = lines_kdtrees,
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
    idx, distance = NearestNeighbors.nn(line_kdtree.line_tree, point)
    width = line_kdtree.width[idx]
    return distance, width, idx
end

function get_density(lines_kdtrees::Array{LineKDTree,1}, r, z)
    distances = []
    densities = []
    point = [r, z]
    for (i, line_kdtree) in enumerate(lines_kdtrees)
        if z > line_kdtree.zmax || z < line_kdtree.z0
            #println("i $i zmax $(line_kdtree.zmax)")
            continue
        end
        distance, width, idx = get_closest_point(line_kdtree, point)
        if distance > width
            #println("i $i distance $distance width $width idx $idx")
            continue
        end
        #println("i $i distance $distance n $n")
        #println("r $(line_kdtree.r[idx]) z $(line_kdtree.z[idx])")
        z_idx = searchsorted_nearest(line_kdtree.z, z)
        n = line_kdtree.n[z_idx]
        distance = d_euclidean(line_kdtree.r[z_idx], r, line_kdtree.z[z_idx], z)
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

