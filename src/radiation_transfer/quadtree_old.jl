using RegionTrees
using Distances
using Printf
using StaticArrays: SVector
import Base: copy, ==
export CellData,
    initialize_leaf,
    initialize_quadtree,
    quadtree_fill_point,
    quadtree_fill_timestep,
    quadtree_fill_line,
    quadtree_fill_all_lines,
    quadtree_fill_horizontal,
    quadtree_erase_line,
    compute_tau_cell,
    quadtree_deltatau,
    quadtree_effective_density,
    fill_last_point,
    count_leaves,
    get_width,
    get_height,
    get_rmin,
    get_zmin,
    get_rmax,
    get_zmax,
    normalize_point,
    denormalize_point,
    compute_accumulative_taus,
    get_index


struct CellData
    line_id::Array{Int,1}
    z_positions::Array{Float64,1}
    densities::Array{Float64,1}
    taus::Array{Float64,1}
end

"Method of the copy function for Cell type"
function copy(cell::Cell)
    newcell = Cell(cell.boundary, cell.data)
    return newcell
end

"Method of the == function for the Cell type"
function ==(cell1::Cell, cell2::Cell)
    return cell1.boundary == cell2.boundary
end

get_width(cell::Cell) = cell.boundary.widths[1]
get_height(cell::Cell) = cell.boundary.widths[2]
get_rmin(cell::Cell) = cell.boundary.origin[1]
get_rmax(cell::Cell) = cell.boundary.origin[1] + cell.boundary.widths[1]
get_zmin(cell::Cell) = cell.boundary.origin[2]
get_zmax(cell::Cell) = cell.boundary.origin[2] + cell.boundary.widths[2]
count_leaves(quadtree::Cell) = length([leaf for leaf in allleaves(quadtree)])

"Following the line initialpoint -> finalpoint, computes the intersection after the current point
to the next tree leaf"
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

function initialize_leaf(cell, child_indices)
    data = CellData(Int[], Float64[], Float64[], Float64[])
    return data
end

function get_cell_coordinates(cell)
    vcs = vertices(cell)
    xmin = vcs[1][1]
    ymin = vcs[1][2]
    xmax = vcs[4][1]
    ymax = vcs[4][2]
    return xmin, ymin, xmax, ymax
end

function initialize_quadtree(r_min, r_max, z_min, z_max)
    quadtree = Cell(
        SVector(0.0, 0.0),
        SVector(r_max, z_max),
        CellData(Int[], Float64[], Float64[], Float64[]),
    )
    return quadtree
end

initialize_quadtree(grid::Grid) =
    initialize_quadtree(grid.r_min, grid.r_max, grid.z_min, grid.z_max)

function refine_leaf(quadtree::Cell, point, leaf, targetsize)
    while get_width(leaf) > targetsize
        split!(leaf, initialize_leaf)
        leaf = findleaf(quadtree, point)
    end
    return leaf
end

function compute_accumulative_taus(z_positions, densities, z0)
    taus = zero(densities)
    taus[1] = densities[1] * (z_positions[1] - z0)
    @assert(taus[1] >= 0)
    tautotal = taus[1]
    for i = 2:length(taus)
        tautotal += densities[i-1] * (z_positions[i] - z_positions[i-1])
        taus[i] = tautotal
        @assert(taus[i] >= 0)
    end
    return taus
end

function quadtree_fill_point(
    quadtree::Cell,
    line::Streamline,
    point,
    zstep,
    density,
    needs_refinement = true,
)
    leaf = findleaf(quadtree, point)
    z0 = leaf.boundary.origin[2]
    try
        @assert z0 <= point[2]
    catch
        println("z0 lower than point???")
        println("point: $point")
        println("leaf: $leaf")
        throw(DomainError)
    end
    #if line_id == 0 #erasing line
    #    leaf.data.line_id = Int[]
    #    leaf.data.z_positions = Float64[]
    #    leaf.data.densities = Float64[]
    #    return leaf
    #end
    if length(leaf.data.line_id) == 0 && needs_refinement # empty cell
        leaf = refine_leaf(quadtree, point, leaf, zstep)
        z0 = leaf.boundary.origin[2]
        push!(leaf.data.line_id, line.id)
        push!(leaf.data.z_positions, point[2])
        push!(leaf.data.densities, density)
        push!(leaf.data.taus, (point[2] - z0) * density)
        return leaf
    end
    push!(leaf.data.line_id, line.id)
    push!(leaf.data.z_positions, point[2])
    push!(leaf.data.densities, density)
    sorting_idx = sortperm(leaf.data.z_positions)
    leaf.data.line_id .= leaf.data.line_id[sorting_idx]
    leaf.data.z_positions .= leaf.data.z_positions[sorting_idx]
    leaf.data.densities .= leaf.data.densities[sorting_idx]
    try
        @assert z0 <= leaf.data.z_positions[1]
    catch
        println("leaf: $leaf")
        println("z0 : $z0")
        println("pos: $(leaf.data.z_positions[1])")
        throw(DomainError)
    end
    try
        push!(leaf.data.taus, 0.0)
        leaf.data.taus .= compute_accumulative_taus(
            leaf.data.z_positions,
            leaf.data.densities,
            z0,
        )
    catch
        println("acc tau failed")
        println(leaf.data.z_positions)
        println(point)
        throw(DomainError)
    end
    return leaf
end


function quadtree_fill_timestep(
    quadtree::Cell,
    line::Streamline,
    point1,
    point2,
    density;
    minimum_cell_size = 0.1,
)
    r1, z1 = point1
    r2, z2 = point2
    if point1 == point2
        return nothing
    end
    if r2 > get_width(quadtree) || z2 > get_height(quadtree)
        return nothing
    end
    if z2 < z1
        point1[2], point2[2] = point2[2], point1[2]#z1, z2 = z2, z1
    end
    if r2 < r1
        point1[1], point2[1] = point2[1], point1[1]#z1, z2 = z2, z1
    end
    #r1 > r2 ? backwards = true : zs = false
    currentleaf = findleaf(quadtree, point1)
    previousleaf = copy(currentleaf)
    currentpoint = copy(point1)
    deltatau = 0.1 / (density * line.Rg * SIGMA_T)
    height = abs(point2[2] - point1[2])#point1[2] / 100
    linewidth_minsize = line.width_norm * point1[1] / 2.0
    zstep = max(min(deltatau, height), minimum_cell_size)
    while true
        rleft = currentpoint[1] * (1.0 - line.width_norm ./ 2.0)
        rright = currentpoint[1] * (1.0 + line.width_norm / 2.0)
        quadtree_fill_horizontal(
            quadtree,
            line,
            rleft,
            rright,
            currentpoint[2],
            currentleaf,
            zstep,
            density,
            true,
        )
        currentleaf = findleaf(quadtree, currentpoint)
        currentpoint =
            compute_cell_intersection(currentleaf, currentpoint, point1, point2)
        #backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(quadtree, currentpoint)
        if currentpoint[2] > point2[2] || previousleaf == currentleaf
            break
        end
        previousleaf = currentleaf
    end
    rleft = point2[1] * (1.0 - line.width_norm / 2)
    rright = point2[1] * (1.0 + line.width_norm / 2)
    quadtree_fill_horizontal(
        quadtree,
        line,
        rleft,
        rright,
        point2[2],
        currentleaf,
        zstep,
        density,
        true,
    )
end

quadtree_fill_timestep(quadtree::Cell, line::Streamline) =
    quadtree_fill_timestep(
        quadtree,
        line,
        get_last_point(line),
        get_current_point(line),
        line.number_density[end],
    )

function quadtree_fill_horizontal(
    quadtree::Cell,
    line::Streamline,
    rleft,
    rright,
    z,
    leaf,
    zstep,
    density,
    needs_refinement = true,
)
    rl = [rleft, z]
    rr = [rright, z]
    currentpoint = rl
    previousleaf = copy(leaf)
    while true
        leaf = quadtree_fill_point(
            quadtree,
            line,
            currentpoint,
            zstep,
            density,
            needs_refinement,
        )
        currentpoint = compute_cell_intersection(leaf, currentpoint, rl, rr)
        if currentpoint[1] > min(rright, get_width(quadtree)) ||
           previousleaf == leaf
            break
        end
        previousleaf = copy(leaf)
    end
end

function quadtree_fill_line(quadtree::Cell, line::Streamline)
    for i = 1:length(line.r)-1
        point1 = [line.r[i], line.z[i]]
        point2 = [line.r[i+1], line.z[i+1]]
        density = line.number_density[i]
        quadtree_fill_timestep(quadtree, line, point1, point2, density)
    end
end

"Fills and refines data from all streamlines"
function quadtree_fill_all_lines(lines)
    for line in lines
        println(@sprintf("Line %02d of %02d", line.id, length(lines)))
        quadtree_fill_line(quadtree, line)
    end
end

function quadtree_erase_line(quadtree::Cell, lines, line_id)
    line = lines[line_id]
    #line.p.n_hist .= wind.grids.n_vacuum
    quadtree_fill_line(quadtree, line)
end

#function quadtree_deltatau(quadtree::Cell, lines, point1, point2)
#    leaf = findleaf(quadtree, point1)
#    line = lines[leaf.data.line_id]
#    compute_tau_cell(quadtree, line, point1, point2, leaf)
#end

function quadtree_effective_density(quadtree::Cell, point1, point2, n_vaccum=1e2)
    leaf = findleaf(quadtree, point1)
    tau = compute_tau_cell(quadtree, leaf, point1, point2, n_vaccum)
    deltad = evaluate(Euclidean(), point1, point2)
    return tau / deltad
end

function compute_tau_cell(quadtree::Cell, leaf, point1, point2, n_vacuum=1e2)
    if isapprox(point1[2], point2[2], rtol = 0, atol = 1e4 * eps(Float64))
        return 0.0
    end
    try
        @assert point2[2] > point1[2]
    catch
        println(point1)
        println(point2)
        println(leaf)
        throw(DomainError)
    end
    deltad = evaluate(Euclidean(), point1, point2)
    if length(leaf.data.line_id) == 0
        return n_vacuum * deltad
    end
    if point1[2] < leaf.data.z_positions[1]
        z_1 = leaf.boundary.origin[2]
        den_1 = leaf.data.densities[1]
        tau_1 = den_1 * (point1[2] - z_1)
    else
        z_idx_1 = get_index(leaf.data.z_positions, point1[2])
        z_1 = leaf.data.z_positions[z_idx_1]
        den_1 = leaf.data.densities[z_idx_1]
        tau_1 = leaf.data.taus[z_idx_1] + den_1 * (point1[2] - z_1)
    end
    if point2[2] < leaf.data.z_positions[1]
        z_2 = leaf.boundary.origin[2]
        den_2 = leaf.data.densities[1]
        tau_2 = den_2 * (point2[2] - z_2)
    else
        z_idx_2 = get_index(leaf.data.z_positions, point2[2])
        z_2 = leaf.data.z_positions[z_idx_2]
        den_2 = leaf.data.densities[z_idx_2]
        tau_2 = leaf.data.taus[z_idx_2] + den_2 * (point2[2] - z_2)
    end
    try
        @assert point2[2] >= z_2
    catch
        println("point2: $point2")
        println("z_2: $z_2")
        throw(DomainError)
    end
    deltaz = point2[2] - point1[2]
    deltatau = tau_2 - tau_1
    @assert deltatau >= 0
    tau = deltatau / deltaz * deltad
    return tau
end

function fill_last_point(quadtree::Cell, line)
    currentpoint = [line.r[end], line.z[end]]
    previouspoint = [line.r[end-1], line.z[end-1]]
    quadtree_fill_timestep(
        quadtree,
        line,
        currentpoint,
        previouspoint,
        line.number_density[end-1],
    )
end
