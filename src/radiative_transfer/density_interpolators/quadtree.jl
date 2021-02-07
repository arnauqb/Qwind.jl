using RegionTrees, StaticArrays
using Statistics: mean, std
import RegionTrees: AbstractRefinery, needs_refinement, refine_data
import Base.length
import Base: copy, ==
export create_quadtree,
    refine_quadtree!,
    QuadtreeRefinery,
    create_and_refine_quadtree,
    getwidth,
    getheight,
    getrmin,
    getrmax,
    getzmin,
    getzmax,
    getcellsize,
    getlimits,
    countleaves,
    get_density,
    compute_cell_optical_thickness

# functions for quadtree cells

"Method of the copy function for Cell type"
function copy(cell::Cell)
    newcell = Cell(cell.boundary, cell.data)
    return newcell
end

"Method of the == function for the Cell type"
function ==(cell1::Cell, cell2::Cell)
    return cell1.boundary == cell2.boundary
end

getwidth(cell::Cell) = cell.boundary.widths[1]
getheight(cell::Cell) = cell.boundary.widths[2]
getcellsize(cell::Cell) = sqrt(sum(cell.boundary.widths .^ 2))
getrmin(cell::Cell) = cell.boundary.origin[1]
getrmax(cell::Cell) = cell.boundary.origin[1] + cell.boundary.widths[1]
getzmin(cell::Cell) = cell.boundary.origin[2]
getzmax(cell::Cell) = cell.boundary.origin[2] + cell.boundary.widths[2]
getcenter(cell::Cell) = cell.boundary.origin + cell.boundary.widths / 2
function is_out(cell::Cell, r, z)
    rmin, zmin = cell.boundary.origin
    rmax, zmax = cell.boundary.origin + cell.boundary.widths
    if r > rmax || z > zmax || r < rmin || z < zmin
        return true
    else
        return false
    end
end
function compute_cell_optical_thickness(cell::Cell, Rg)
    size = getcellsize(cell)
    return size * cell.data * SIGMA_T * Rg
end
getlimits(cell::Cell) = [getrmin(cell), getrmax(cell), getzmin(cell), getzmax(cell)]
countleaves(quadtree::Cell) = length([leaf for leaf in allleaves(quadtree)])
function get_density(cell::Cell, r, z)
    if is_out(cell, r, z)
        return 1e2
    else
        return findleaf(cell, [r, z]).data
    end
end

"""
Given a line defined by initialpoint and finalpoint, computes the intersection
with the next intersection with the cell boundary from the currentpoint
"""
function compute_cell_intersection(cell::Cell, currentpoint, initialpoint, finalpoint)
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
    #intersection = [r, z] + lambda .* [rf - ri, zf - zi]
    intersection = [r + lambda * (rf-ri), z + lambda * (zf - zi)]
    try
        @assert lambda >= 0.0
    catch
        throw(DomainError)
    end
    return intersection
end

function create_quadtree(r_min, r_max, z_min, z_max; vacuum_density = 1e2)
    quadtree =
        Cell(SVector(r_min, z_min), SVector(r_max - r_min, z_max - z_min), vacuum_density)
    return quadtree
end

# Quadtree refinement

struct QuadtreeRefinery <: AbstractRefinery
    density_interpolator::DensityInterpolator
    Rg::Float64
    atol::Float64
    rtol::Float64
    cell_min_size::Float64
end

function does_tau_need_refinement(taus; atol, rtol)
    taus = filter(!isinf, (filter(!isnan, taus)))
    mean_tau = mean(taus)
    std_tau = std(taus)
    normalized_std = std_tau / mean_tau
    if normalized_std > rtol && std_tau > atol
        return true
    else
        return false
    end
end

function sample_densities_in_cell(cell::Cell, density_interpolator::DensityInterpolator, n=250)
    rs = range(
        cell.boundary.origin[1],
        cell.boundary.origin[1] + cell.boundary.widths[1],
        length = n,
    )
    #zs = 10 .^ range(
    #    max(log10(cell.boundary.origin[2]), 1e-6),
    #    log10(cell.boundary.origin[2] + cell.boundary.widths[2]),
    #    length = n,
    #)
    zs = range(
        cell.boundary.origin[2],
        cell.boundary.origin[2] + cell.boundary.widths[2],
        length = n,
    )
    points = hcat(rs, zs)'
    densities = get_density_points(density_interpolator, points)
    return densities
end

function needs_refinement(refinery::QuadtreeRefinery, cell)
    densities = sample_densities_in_cell(cell, refinery.density_interpolator)
    delta_d = getcellsize(cell)
    taus = densities .* delta_d * refinery.Rg * SIGMA_T
    does_tau_vary =
        does_tau_need_refinement(taus, atol = refinery.atol, rtol = refinery.rtol) 
    is_cell_big_enough = delta_d > refinery.cell_min_size
    refine_condition = does_tau_vary && is_cell_big_enough
    if !refine_condition
        cell.data = mean(densities)
    end
    return refine_condition
end

function refine_data(refinery::QuadtreeRefinery, cell::Cell, indices)
    return mean(sample_densities_in_cell(cell, refinery.density_interpolator))
end

function refine_quadtree!(
    quadtree::Cell,
    density_interpolator::DensityInterpolator,
    Rg,
    atol,
    rtol,
    cell_min_size,
)
    ref = QuadtreeRefinery(density_interpolator, Rg, atol, rtol, cell_min_size)
    adaptivesampling!(quadtree, ref)
end

function create_and_refine_quadtree(
    density_interpolator::Union{DensityInterpolator, Nothing};
    Rg,
    atol,
    rtol,
    cell_min_size,
    vacuum_density = 1e2,
)
    if density_interpolator === nothing
        return Cell(SVector(0, 0), SVector(1, 1), vacuum_density)
    end
    rmax = r_max(density_interpolator)
    zmax = z_max(density_interpolator)
    quadtree = create_quadtree(0.0, rmax, 0.0, zmax, vacuum_density = vacuum_density)
    @info "Refining quadtree..."
    refine_quadtree!(quadtree, density_interpolator, Rg, atol, rtol, cell_min_size)
    @info "Done"
    return quadtree
end