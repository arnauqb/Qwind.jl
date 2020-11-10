using RegionTrees, StaticArrays
using Statistics: mean, std
import RegionTrees: AbstractRefinery, needs_refinement, refine_data
import Base.length
import Base: copy, ==
export getdensity,
    create_quadtree,
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
    getdensity,
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
function compute_cell_optical_thickness(cell::Cell, Rg)
    size = getcellsize(cell)
    return size * cell.data * SIGMA_T * Rg
end
getlimits(cell::Cell) =
    [getrmin(cell), getrmax(cell), getzmin(cell), getzmax(cell)]
countleaves(quadtree::Cell) = length([leaf for leaf in allleaves(quadtree)])
getdensity(cell::Cell, r, z) = findleaf(cell, [r, z]).data

"""
Given a line defined by initialpoint and finalpoint, computes the intersection
with the next intersection with the cell boundary from the currentpoint
"""
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
        throw(DomainError)
    end
    return intersection
end

function create_quadtree(r_min, r_max, z_min, z_max; vacuum_density = 1e2)
    quadtree = Cell(SVector(r_min, z_min), SVector(r_max - r_min, z_max - z_min), vacuum_density)
    return quadtree
end

# Quadtree refinement

struct QuadtreeRefinery <: AbstractRefinery
    windkdtree::WindKDTree
    Rg::Float64
    tau_max_std::Float64
    cell_min_size::Float64
end

function getdata(cell, child_indices, Rg, windtree)
    r, z = quadtree.boundary.origin + quadtree.boundary.widths / 2
    density = get_density_from_tree(windtree, r, z)
    return density
end

getdata(cell, child_indices) =
    getdata(cell, child_indices, black_hole.Rg, windtree)


function needs_refinement(refinery::QuadtreeRefinery, cell)
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
            push!(densities, getdensity(refinery.windkdtree, rp, zp))
        end
    end
    delta_d = getcellsize(cell)
    taus = densities .* delta_d * refinery.Rg * SIGMA_T
    refine_condition =
        std(taus) > refinery.tau_max_std && delta_d > refinery.cell_min_size
    if !refine_condition
        cell.data = mean(densities)
    end
    refine_condition
end

function refine_data(refinery::QuadtreeRefinery, cell::Cell, indices)
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
            push!(densities, getdensity(refinery.windkdtree, rp, zp))
        end
    end
    mean(densities)
end

function refine_quadtree!(
    quadtree::Cell,
    windkdtree::WindKDTree,
    Rg,
    tau_max_std,
    cell_min_size,
)
    ref = QuadtreeRefinery(windkdtree, Rg, tau_max_std, cell_min_size)
    adaptivesampling!(quadtree, ref)
end

function create_and_refine_quadtree(
    windkdtree::WindKDTree,
    Rg,
    tau_max_std,
    cell_min_size,
    vacuum_density = 1e2,
)
    rmax = maximum(windkdtree.r)
    zmax = maximum(windkdtree.z)
    quadtree =
        create_quadtree(0.0, rmax, 0.0, zmax, vacuum_density = vacuum_density)
    refine_quadtree!(quadtree, windkdtree, Rg, tau_max_std, cell_min_size)
    return quadtree
end
