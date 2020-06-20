using RegionTrees
using StaticArrays: SVector
using RegionTrees
using VoronoiDelaunay
using Distances
using Printf
export initialize_leaf,
       get_cell_coordinates,
       quadtree_initialize,
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
       countleaves,
       normalize_point,
       denormalize_point,
       compute_accumulative_taus



function initialize_leaf(cell, child_indices)
    #data = [0, Array{Float64}(undef, 2, 0), Float64[]] # line_id 0 -> vaccum
    #data = [0, Float64[], Float64[], 0] # line_id 0 -> vaccum
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

function quadtree_initialize(wind)
    height = wind.config["grids"]["z_max"]
    width = wind.config["grids"]["r_max"]
    wind.quadtree = Cell(SVector(0., 0.),
                    SVector(width, height),
                    CellData(Int[], Float64[], Float64[], Float64[])
                    #[0, Array{Float64}(undef, 2, 0), Float64[], 0] # line_id position density
                    )
end

function refine_leaf(point, leaf, targetsize, wind)
    while cell_width(leaf) > targetsize
        split!(leaf, initialize_leaf)
        leaf = findleaf(wind.quadtree, point)
    end
    return leaf
end

function compute_accumulative_taus(z_positions, densities, z0)
    taus = zero(densities)
    taus[1] = densities[1] * (z_positions[1] - z0)
    @assert(taus[1]>=0)
    tautotal = taus[1]
    for i in 2:length(taus)
        tautotal += densities[i-1] * (z_positions[i] - z_positions[i-1])
        taus[i] = tautotal
        @assert(taus[i]>=0)
    end
    return taus
end

function quadtree_fill_point(point, zstep, density, line_id, line, wind; needs_refinement=true)
    leaf = findleaf(wind.quadtree, point)
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
        leaf = refine_leaf(point, leaf, zstep, wind)
        z0 = leaf.boundary.origin[2]
        push!(leaf.data.line_id, line_id)
        push!(leaf.data.z_positions, point[2])
        push!(leaf.data.densities, density)
        push!(leaf.data.taus, (point[2] - z0) * density)
        return leaf
    end
    push!(leaf.data.line_id, line_id)
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
        leaf.data.taus .= compute_accumulative_taus(leaf.data.z_positions, leaf.data.densities, z0)
    catch
        println("acc tau failed")
        println(leaf.data.z_positions)
        println(point)
        throw(DomainError)
    end
    return leaf
end

function quadtree_fill_timestep(point1, point2, density, linewidthnorm, line_id, line, wind::WindStruct)
    r1, z1 = point1
    r2, z2 = point2
    if point1 == point2
        return nothing
    end
    if r2 > wind.r_max || z2 > wind.z_max
        return nothing
    end
    if z2 < z1
        point1[2], point2[2] = point2[2], point1[2]#z1, z2 = z2, z1
    end
    if r2 < r1
        point1[1], point2[1] = point2[1], point1[1]#z1, z2 = z2, z1
    end
    #r1 > r2 ? backwards = true : backwards = false
    currentleaf = findleaf(wind.quadtree, point1)
    previousleaf = copy(currentleaf)
    currentpoint = copy(point1)
    deltatau = 0.1 / (density * wind.bh.R_g * SIGMA_T)
    height = abs(point2[2] - point1[2])#point1[2] / 100
    linewidth_minsize = linewidthnorm * point1[1] / 2.0
    zstep = min(max(deltatau, wind.config["grids"]["minimum_cell_size"]), linewidth_minsize)
    while true
        rleft = currentpoint[1] * (1.0 - linewidthnorm/2.0)
        rright = currentpoint[1] * (1.0 + linewidthnorm/2.0)
        quadtree_fill_horizontal(rleft, rright, currentpoint[2], currentleaf, zstep, density, line_id, line, wind, true)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        #backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentpoint[2] > point2[2] || previousleaf == currentleaf
            break
        end
        previousleaf = currentleaf
    end
    rleft = point2[1] * (1.0 - linewidthnorm/2)
    rright = point2[1] * (1.0 + linewidthnorm/2)
    quadtree_fill_horizontal(rleft, rright, point2[2], currentleaf, zstep, density, line_id, line, wind, true)
end

function quadtree_fill_horizontal(rleft, rright, z, leaf, zstep, density, line_id, line, wind, needs_refinement=true)
    rl = [rleft, z]
    rr = [rright, z]
    currentpoint = rl
    previousleaf = copy(leaf)
    while true
        leaf = quadtree_fill_point(currentpoint, zstep, density, line_id, line, wind, needs_refinement=needs_refinement)
        currentpoint = compute_cell_intersection(currentpoint, leaf, rl, rr)
        if currentpoint[1] > min(rright, wind.r_max) || previousleaf == leaf
            break
        end
        previousleaf = copy(leaf)
    end
end

function quadtree_fill_line(line, line_id, wind)
    for i in 1:length(line.p.u_hist[:,1])-1
        point1 = [line.p.u_hist[i,1], line.p.u_hist[i,2]]
        point2 = [line.p.u_hist[i+1,1], line.p.u_hist[i+1,2]]
        lw_norm = line.p.line_width / line.p.r_0
        density = line.p.n_hist[i]
        quadtree_fill_timestep(point1, point2, density, lw_norm, line_id, line, wind)
    end
end

"Fills and refines data from all streamlines"
function quadtree_fill_all_lines(wind::WindStruct)
    for i in 1:length(wind.lines)
        println(@sprintf("Line %02d of %02d", i, length(wind.lines)))
        isassigned(wind.lines, i) || continue
        wind.lines[i] === nothing && continue
        line = wind.lines[i]
        quadtree_fill_line(line, i, wind)
    end
end

function quadtree_erase_line(line_id, wind)
    line = wind.lines[line_id]
    #line.p.n_hist .= wind.grids.n_vacuum
    quadtree_fill_line(line, 0, wind)
end

function quadtree_deltatau(point1, point2, wind)
    leaf = findleaf(wind.quadtree, point1)
    line = wind.lines[leaf.data.line_id]
    compute_tau_cell(line, point1, point2, leaf, wind)
end

function quadtree_effective_density(point1, point2, wind)
    leaf = findleaf(wind.quadtree, point1)
    tau = compute_tau_cell(point1, point2, leaf, wind)
    deltad = evaluate(Euclidean(), point1, point2)
    n = tau / deltad
    return n
end

function compute_tau_cell(point1, point2, leaf, wind)
    if isapprox(point1[2], point2[2], rtol=0, atol=1e4*eps(Float64))
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
        return wind.grids.n_vacuum * deltad
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

function fill_last_point(line)
    currentpoint = line.p.u_hist[end,1:2]
    previouspoint = line.p.u_hist[end-1,1:2]
    lwnorm = line.p.line_width / line.p.r_0
    density = line.p.n_hist[end-1]
    line_id = line.p.line_id
    wind = line.p.wind
    quadtree_fill_timestep(currentpoint, previouspoint, density, lwnorm, line_id, line, wind)
end

function countleaves(wind)
    counter = 0
    for leaf in allleaves(wind.quadtree)
        counter += 1
    end
    return counter
end
