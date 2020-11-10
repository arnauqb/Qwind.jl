using NearestNeighbors, RegionTrees, StaticArrays
using Statistics: mean, std
import Qwind.compute_xray_tau
import Qwind.compute_ionization_parameter
import RegionTrees: AbstractRefinery, needs_refinement, refine_data
import Base.length
import Base: copy, ==
export initialize_wind_tree,
    get_density_from_tree,
    compute_tau,
    compute_xray_tau,
    compute_ionization_parameter,
    compute_cell_intersection,
    refine_quadtree!,
    get_width,
    get_height,
    get_rmin,
    get_rmax,
    get_zmin,
    get_zmax,
    count_leaves,
    refine_data,
    needs_refinement

abstract type Refinery <: AbstractRefinery end

struct WindTree <: Refinery
    r::Vector{Float64}
    z::Vector{Float64}
    width::Vector{Float64}
    n::Vector{Float64}
    kdtree::KDTree
    quadtree::Cell
    tau_max_std::Float64
    cell_min_size::Float64
    Rg::Float64
end

"Method of the == function for the Cell type"
function ==(cell1::Cell, cell2::Cell)
    return cell1.boundary == cell2.boundary
end


function initialize_quadtree(r_min, r_max, z_min, z_max; vacuum_density = 1e2)
    quadtree =
        Cell(SVector(r_min, z_min), SVector(r_max, z_max), vacuum_density)
    return quadtree
end

function initialize_wind_tree(
    r::Vector{Float64},
    z::Vector{Float64},
    width::Vector{Float64},
    n::Vector{Float64},
    Rg::Float64;
    cell_min_size = 100.0,
    tau_max_std = 0.01,
)
    points = hcat(r, z)'
    kdtree = KDTree(points)
    quadtree =
        initialize_quadtree(minimum(r), maximum(r), minimum(z), maximum(z))
    return WindTree(
        r,
        z,
        width,
        n,
        kdtree,
        quadtree,
        tau_max_std,
        cell_min_size,
        Rg,
    )
end

function initialize_wind_tree(flag::FirstIter)
    kdtree = KDTree(rand(2, 1))
    quadtree = initialize_quadtree(0.0, 1.0, 0.0, 1.0)
    return WindTree(
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        kdtree,
        quadtree,
        NaN,
        NaN,
        NaN,
    )
end

function initialize_wind_tree(solvers; n_time = 10000, tau_max_std=0.01, cell_min_size=100.0)
    r_dense = Float64[]
    z_dense = Float64[]
    n_dense = Float64[]
    v_r_dense = Float64[]
    v_z_dense = Float64[]
    line_width_dense = Float64[]
    density_dense = Float64[]
    for solver in solvers
        t_range =
            10 .^ range(
                max(log10(solver.sol.t[1]), 1e-6),
                log10(solver.sol.t[end-2]),
                length = n_time,
            )
        dense_solution = solver.sol(t_range)
        r_dense = vcat(r_dense, dense_solution[1, :])
        z_dense = vcat(z_dense, dense_solution[2, :])
        v_r_dense = vcat(v_r_dense, dense_solution[3, :])
        v_z_dense = vcat(v_z_dense, dense_solution[4, :])
        line_width_dense = vcat(
            line_width_dense,
            dense_solution[1, :] .* solver.p.line.width_norm,
        )
        v_t = @. sqrt(dense_solution[3,:]^2 + dense_solution[4,:]^2)
        density_dense = vcat(
            density_dense,
            compute_density.(dense_solution[1,:], dense_solution[2,:], v_t, Ref(solver.p.line)),
        )
    end
    radiation = solvers[1].p.radiation
    return initialize_wind_tree(
        r_dense,
        z_dense,
        line_width_dense,
        density_dense,
        radiation.bh.Rg,
        cell_min_size=cell_min_size,
        tau_max_std=tau_max_std
    )
end

function get_density_from_tree(windtree::WindTree, r, z)
    idcs, distances = knn(windtree.kdtree, [r, z], 2)
    r1, r2 = windtree.r[idcs]
    z1, z2 = windtree.z[idcs]
    point = [r1, z1] + [r2-r1, z2-z1] / 2
    distance = evaluate(Euclidean(), point, [r,z])
    width = sum(windtree.width[idcs]) / 4 #2 of average and 2 of half width
    density = sum(windtree.n[idcs]) / 2
    if distance <= width #&& x < 0
        return density
    else
        return max(density * exp(-(distance - width)^2 / 2), 1e2)
    end
end

# optical depths

#function compute_tau_integrand(windtree::WindTree, tau, t, r0, r1, z0, z1)
#    r = 0.5 * ((r1 - r0) * t + (r1 + r0))
#    z = 0.5 * ((z1 - z0) * t + (z1 + z0))
#    return get_density_from_tree(windtree, r, z)
#end
#
#function compute_tau(windtree::WindTree, r, z, r0, z0, Rg)
#    func(t) = compute_tau_integrand(windtree, t, r0, r, z0, z)
#    (val, err) =
#        quadgk(func, -1, 1; rtol = 1e-3, atol = 0, maxevals = 10^4, order = 5)
#    return val * 0.5 * sqrt((r - r0)^2 + (z - z0)^2) * Rg * SIGMA_T
#end
#
#function compute_xray_tau(
#    windtree::WindTree,
#    r,
#    z,
#    xray_luminosity,
#    Rg,
#    r0 = 0.0,
#    z0 = 0.0;
#    n_points = 100,
#)
#    tau = 0.0
#    t_range = range(0, 1, length = n_points)
#    delta_t = 1 / n_points
#    xi_comput(r, z, n, tau) =
#        compute_ionization_parameter(r, z, n, tau, xray_luminosity, Rg)
#    for t in t_range
#        rp = r0 + t * (r - r0)
#        zp = z0 + t * (z - z0)
#        n = get_density_from_tree(windtree, rp, zp)
#        tau_x = tau * sqrt((rp - r0)^2 + (z - zp)^2)
#        ξ = xi_comput(rp, zp, n, tau_x)
#        sigma = SIGMA_T #compute_xray_opacity(ξ)
#        tau += n * sigma * Rg * delta_t
#    end
#    return tau * sqrt((r - r0)^2 + (z - z0^2))
#end
#
#function compute_ionization_parameter(
#    windtree::WindTree,
#    r,
#    z,
#    xray_luminosity,
#    Rg,
#    r0 = 0.0,
#    z0 = 0.0;
#    n_points = 100,
#)
#    tau = 0.0
#    ξ = 0.0
#    t_range = range(0, 1, length = n_points)
#    delta_t = 1 / n_points
#    xi_comput(r, z, n, tau) =
#        compute_ionization_parameter(r, z, n, tau, xray_luminosity, Rg)
#    for t in t_range
#        rp = r0 + t * (r - r0)
#        zp = z0 + t * (z - z0)
#        n = get_density_from_tree(windtree, rp, zp)
#        tau_x = tau * sqrt((rp - r0)^2 + (z - zp)^2)
#        ξ = xi_comput(rp, zp, n, tau_x)
#        if ξ < 1e-20
#            break
#        end
#        sigma = compute_xray_opacity(ξ)
#        tau += n * sigma * Rg * delta_t
#    end
#    return ξ
#end


function compute_cell_intersection(
    cell::Cell,
    currentpoint,
    initialpoint,
    finalpoint,
)
    r, z = currentpoint
    ri, zi = initialpoint
    rf, zf = finalpoint
    direction_r = sign(rf - ri)
    direction_z = sign(zf - zi)
    cell_vertices = vertices(cell)
    direction_r == 1 ? cell_r_boundary = cell_vertices[2, 2][1] :
    cell_r_boundary = cell_vertices[1, 1][1]
    direction_z == 1 ? cell_z_boundary = cell_vertices[2, 2][2] :
    cell_z_boundary = cell_vertices[1, 1][2]
    ri == rf ? lambda_r = Inf : lambda_r = (cell_r_boundary - r) / (rf - ri)
    zi == zf ? lambda_z = Inf : lambda_z = (cell_z_boundary - z) / (zf - zi)
    lambda = min(lambda_r, lambda_z)
    intersection = [r, z] + lambda .* [rf - ri, zf - zi]
    try
        @assert lambda >= 0.0
    catch
        println("lambda: $lambda")
        println("current point : $currentpoint")
        println("cell: $cell")
        println("initial point: $initialpoint")
        println("final point: $finalpoint")
        throw(DomainError)
    end
    return intersection
end

# Quadtree refinement

function getdata(cell, child_indices, Rg, windtree)
    r, z = quadtree.boundary.origin + quadtree.boundary.widths / 2
    density = get_density_from_tree(windtree, r, z)
    return density
end

getdata(cell, child_indices) =
    getdata(cell, child_indices, black_hole.Rg, windtree)


function needs_refinement(windtree::Refinery, cell)
    r, z = cell.boundary.origin + cell.boundary.widths / 2
    rs = range(
        cell.boundary.origin[1],
        cell.boundary.origin[1] + cell.boundary.widths[1],
        length = 10,
    )
    zs = range(
        cell.boundary.origin[2],
        cell.boundary.origin[2] + cell.boundary.widths[2],
        length = 10,
    )
    densities = []
    for rp in rs
        for zp in zs
            push!(densities, get_density_from_tree(windtree, rp, zp))
        end
    end
    delta_d = sqrt(sum(cell.boundary.widths .^ 2))
    taus = densities .* delta_d * windtree.Rg * SIGMA_T
    refine_condition =
        std(taus) > windtree.tau_max_std || delta_d > windtree.cell_min_size
    if !refine_condition
        cell.data = mean(densities)
    end
    refine_condition
end

function refine_data(windtree::Refinery, cell::Cell, indices)
    rs = range(
        cell.boundary.origin[1],
        cell.boundary.origin[1] + cell.boundary.widths[1],
        length = 10,
    )
    zs = range(
        cell.boundary.origin[2],
        cell.boundary.origin[2] + cell.boundary.widths[2],
        length = 10,
    )
    densities = []
    for rp in rs
        for zp in zs
            push!(densities, get_density_from_tree(windtree, rp, zp))
        end
    end
    mean(densities)
end

function refine_quadtree!(quadtree::Cell, windtree)
    adaptivesampling!(quadtree, windtree)
end

refine_quadtree!(windtree::WindTree) =
    refine_quadtree!(windtree.quadtree, windtree)
