export VIGrid, get_density
import Distances, Interpolations
using ProgressMeter

struct NEuclidean <: Metric
    r_norm::Float64
    z_norm::Float64
end

(norm::NEuclidean)(a, b) =
    sqrt(((a[1] - b[1]) / norm.r_norm)^2 + ((a[2] - b[2]) / norm.z_norm)^2)

struct LineKDTree <: NNInterpolator
    r::Vector{Float64}
    z::Vector{Float64}
    n::Vector{Float64}
    r_pos::Vector{Float64}
    z_pos::Vector{Float64}
    width_pos::Vector{Float64}
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
        return new(
            r_den,
            z_den,
            n_den,
            r_pos,
            z_pos,
            width_pos,
            maximum(z),
            r[1] - width_pos[1] / 2,
            minimum(z),
            line_tree,
            n_timesteps,
        )
    end
end

function reduce_line_uniformely_spaced_distance(r, z, width)
    rs = [r[1]]
    zs = [z[1]]
    widths = [width[1]]
    for (rp, zp, wp) in zip(r, z, width)
        d0 = d_euclidean(rp, rs[1], zp, zs[1])
        d = d_euclidean(rp, rs[end], zp, zs[end])
        d_lim = 0.01 * d0
        if zp > zs[end] && d > d_lim
            push!(rs, rp)
            push!(zs, zp)
            push!(widths, wp)
        end
    end
    return rs, zs, widths
end

function reduce_line_uniformely_spaced_density(r, z, n)
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    for (rp, zp, np) in zip(r, z, n)
        ns_eps = abs(np - ns[end]) / ns[end]
        if zp > zs[end] && ns_eps > 0.01
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
        end
    end
    return rs, zs, ns
end

function create_lines_kdtrees(integrators; n_timesteps = 10000)
    ret = LineKDTree[]
    r_norm = maximum(maximum([integrator.p.data[:r] for integrator in integrators]))
    z_norm = maximum(maximum([integrator.p.data[:z] for integrator in integrators]))
    for (i, integrator) in enumerate(integrators)
        r, z, zmax, z0, width, n =
            get_dense_solution_from_integrator(integrator, n_timesteps)
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
        r_max = max(r_max, lk.r[1] + lk.width_pos[1])
        z_min = min(z_min, lk.z0)
        z_max = max(z_max, lk.zmax)
        max_width = max(max_width, maximum(lk.width_pos))
    end
    r_min = max(6, r_min)
    r_max += max_width
    z_min = max(1e-6, z_min)
    if nr == "auto"
        r_range = [line.r0 for line in lines_kdtrees]
        additional_r = collect(range(r_range[end], r_max, step = 100)[2:end])
        r_range = vcat(r_range, additional_r)
    else
        r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    end
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    z_range = pushfirst!(z_range, 0.0)
    r_range = round.(r_range, digits = 7)
    z_range = round.(z_range, digits = 7)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    return r_range, z_range
end

struct VIGrid <: GridInterpolator
    grid::InterpolationGrid
    interpolator::Any
    vacuum_density::Float64
    n_timesteps::Int
    function VIGrid(
        r_range::Union{Vector{Float64},Nothing},
        z_range::Union{Vector{Float64},Nothing},
        grid::Union{Array{Float64,2},Nothing};
        nr = nothing,
        nz = nothing,
        vacuum_density = 1e2,
        n_timesteps = 10000,
    )
        if grid === nothing
            itp = (r, z) -> vacuum_density
        else
            itp = Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
            itp = Interpolations.extrapolate(itp, 1e2)
        end
        return new(
            InterpolationGrid(r_range, z_range, grid, itp, nr, nz),
            itp,
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
    if nr != "auto"
        nr = Int(nr)
    end
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
    nr::Union{Number,String},
    nz::Number;
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if lines_kdtrees === nothing
        r_range = zeros(2)
        z_range = zeros(2)
        density_grid = nothing
    else
        r_range, z_range = get_spatial_grid(lines_kdtrees, nr, nz)
        density_grid = @showprogress pmap(
            z -> get_density.(Ref(lines_kdtrees), r_range, z),
            z_range,
            batch_size = 10,
        )
        density_grid = hcat(density_grid...)
    end
    return VIGrid(
        r_range,
        z_range,
        density_grid,
        nr = nr,
        nz = nz,
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
    #n_timesteps = interpolator.n_timesteps
    #nr = interpolator.grid.nr
    #nz = interpolator.grid.nz
    #lines_kdtrees = create_lines_kdtrees(integrators, n_timesteps = n_timesteps)
    #r_range, z_range = get_spatial_grid(lines_kdtrees, nr, nz)
    #f(r, z) =
    #    10 .^ (
    #        (
    #            log10.(get_density.(Ref(interpolator), r, z)) +
    #            log10.(get_density.(Ref(lines_kdtrees), r, z))
    #        ) / 2
    #    )
    #density_grid = @showprogress pmap(z -> f.(r_range, z), z_range, batch_size = 10)
    #density_grid = hcat(density_grid...)
    #return VIGrid(
    #    r_range,
    #    z_range,
    #    density_grid,
    #    nr = nr,
    #    nz = nz,
    #    vacuum_density = interpolator.vacuum_density,
    #    n_timesteps = n_timesteps,
    #)
end

function get_closest_point(line_kdtree::LineKDTree, point)
    idx, distance = NearestNeighbors.nn(line_kdtree.line_tree, point)
    width = line_kdtree.width_pos[idx]
    return distance, width, idx
end

function get_density(lines_kdtrees::Vector{LineKDTree}, r, z)
    distances = []
    densities = []
    point = [r, z]
    for (i, line_kdtree) in enumerate(lines_kdtrees)
        if z > line_kdtree.zmax || z < line_kdtree.z0
            continue
        end
        distance, width, idx = get_closest_point(line_kdtree, point)
        if distance > width #/ 2
            continue
        end
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

