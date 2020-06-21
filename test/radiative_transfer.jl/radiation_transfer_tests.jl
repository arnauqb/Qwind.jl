using Qwind
using Test
using PyPlot
using RegionTrees
using StaticArrays
using Roots
import RegionTrees: AbstractRefinery, needs_refinement, refine_data
LogNorm = matplotlib.colors.LogNorm

# create black hole
M = 1e8 * M_SUN
mdot = 0.5
spin = 0.0
black_hole = BlackHole(M, mdot, spin)

# create radiation
f_uv = 0.85
f_x = 0.15
shielding_density = 2e8
r_in = 20.0
radiation = SimpleRadiation(black_hole, f_uv, f_x, shielding_density, r_in)

# create simulation grid
r_min = 0.0
r_max = 10000
z_min = 0.0
z_max = 10000
vacuum_density = 1e2
grid = Grid(r_min, r_max, z_min, z_max)

# initialize streamlines
r_fi = 1600
n_lines = 80
velocity_generator(r_0) = 5e7
density_generator(r_0) = 2e8
streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 1e-2,
    log_spaced = true,
)

# initialize quadtree

quadtree = initialize_quadtree(grid)

#run simulation

solvers = []
for line in streamlines.lines
    parameters = Parameters(line, grid, radiation, quadtree)
    solver = initialize_solver(line, parameters)
    push!(solvers, solver)
    run_solver!(solver)
end

windtree = initialize_wind_tree(solvers)

MIN_R = 0
MIN_Z = 0
MAX_R = 1000
MAX_Z = 1000#600#0.21
r_range = range(MIN_R, MAX_R, length = convert(Int, 500))
z_range = range(MIN_Z, MAX_Z, length = convert(Int, 501))
density_grid = 1e2 .* ones((length(r_range), length(z_range)))
xi_grid = 1e2 .* ones((length(r_range), length(z_range)))
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        #xi_grid[i, j] = compute_ionization_parameter(
        #    windtree,
        #    r,
        #    z,
        #    xray_luminosity(radiation),
        #    black_hole.Rg,
        #    n_points=500,
        #)
        density_grid[i, j] = get_density_from_tree(windtree, r, z)
    end
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid', norm = LogNorm())
#cm = ax.pcolormesh(r_range, z_range, xi_grid', norm = LogNorm())
plt.colorbar(cm, ax = ax)
plot_streamlines(streamlines.lines, fig, ax)
ax.set_xlim(MIN_R, MAX_R)
ax.set_ylim(MIN_Z, MAX_Z)
#ax.set_yscale("log")
gcf()




function getdata(cell, child_indices, Rg, windtree)
    r, z = quadtree.boundary.origin + quadtree.boundary.widths / 2
    density = get_density_from_tree(windtree, r, z)
    return density
end

getdata(cell, child_indices) =
    getdata(cell, child_indices, black_hole.Rg, windtree)

struct MyRefinery <: AbstractRefinery
    delta_tau_max::Float64
    windtree::Any
    Rg::Any
end

function needs_refinement(refinery::MyRefinery, cell)
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
            push!(densities, get_density_from_tree(refinery.windtree, rp, zp))
        end
    end
    delta_d = sqrt(sum(cell.boundary.widths .^ 2))
    taus = densities .* delta_d * refinery.Rg * SIGMA_T
    #densities = get_density_from_tree.(Ref(refinery.windtree), rs, zs)
    println("r : $r z : $z")
    println("std: $(std(taus))")
    refine_condition = std(taus) > 0.01 || delta_d > 100
    if !refine_condition
        cell.data = mean(densities)
    end
    refine_condition

    #density = maximum(densities)

    #density = get_density_from_tree(refinery.windtree, r, z)
    #tau = delta_d * refinery.Rg * SIGMA_T * density
    #tau > refinery.delta_tau_max# || delta_d > 50
end

function refine_data(refinery::MyRefinery, cell::Cell, indices)
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
            push!(densities, get_density_from_tree(refinery.windtree, rp, zp))
        end
    end
    mean(densities)
end

quadtree = Cell(SVector(0.0, 0.0), SVector(1000.0, 1000.0), 1e2)

refinery = MyRefinery(1, windtree, black_hole.Rg)

adaptivesampling!(quadtree, refinery)


function compute_xray_tau_leaf2(
    quadtree::Cell,
    leaf::Cell,
    point,
    intersection,
    taux0,
    xray_luminosity,
    Rg,
)
    deltad = evaluate(Euclidean(), point, intersection) * Rg # cell size
    d = evaluate(Euclidean(), [0.0, 0.0], intersection) * Rg # distance from the center
    #density = interpolate_density(leaf.data.line_id, point, leaf, wind) #leaf.data[1] / leaf.data[2]
    #deltatau = compute_tau_cell(quadtree, leaf, point, intersection, 1e2) * Rg
    #n_eff = quadtree_effective_density(quadtree, point, intersection)
    density = leaf.data
    deltatau = density * deltad
    xi0 = xray_luminosity / (density * d^2)
    f(t) =
        t - log(xi0) + taux0 + min(40, deltatau * compute_xray_opacity(exp(t)))
    if f(20) < 0
        xi = xi0
    elseif f(-20) > 0
        xi = 1e-20
    else
        t = find_zero(f, (-20, 20), Bisection(), atol = 0, rtol = 0.1)
        xi = exp(t)
    end
    #xi = 1e6
    taux = compute_xray_opacity(xi) * deltatau
    return taux
end

function compute_xray_tau2(quadtree::Cell, r, z, z_0, xray_luminosity, Rg)
    #return 40
    z = max(z, get_zmin(quadtree))
    point1 = [0.0, z_0]
    #coords_list = [point1]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r, z]
    point2leaf = findleaf(quadtree, point2)
    taux = 0.0
    if point1leaf == point2leaf
        #push!(coords_list, point2)
        taux = compute_xray_tau_leaf2(
            quadtree,
            point1leaf,
            point1,
            point2,
            0.0,
            xray_luminosity,
            Rg,
        )
        return taux
    end
    intersection = compute_cell_intersection2(point1leaf, point1, point1, point2)
    taux = compute_xray_tau_leaf2(
        quadtree,
        point1leaf,
        point1,
        intersection,
        taux,
        xray_luminosity,
        Rg,
    )
    currentpoint = intersection
    #push!(coords_list, intersection)
    currentleaf = findleaf(quadtree, currentpoint)
    previousleaf = copy(currentleaf)
    while currentleaf != point2leaf
        intersection =
            compute_cell_intersection2(currentleaf, currentpoint, point1, point2)
        #push!(coords_list, intersection)
        taux += compute_xray_tau_leaf2(
            quadtree,
            currentleaf,
            currentpoint,
            intersection,
            taux,
            xray_luminosity,
            Rg,
        )
        #if taux > 40
        #    return 40.0
        #end
        currentpoint = intersection
        currentleaf = findleaf(quadtree, currentpoint)
        if currentleaf == previousleaf
            break
        end
        previousleaf = copy(currentleaf)
    end
    if currentpoint[2] < point2[2]
        taux += compute_xray_tau_leaf2(
            quadtree,
            currentleaf,
            currentpoint,
            point2,
            taux,
            xray_luminosity,
            Rg,
        )
    end
    #push!(coords_list, point2)
    return taux
    #if return_coords
    #    return taux, coords_list
    #else
    #    return taux
    #end
end

function compute_cell_intersection2(
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


r_range = range(0, 1000, length = 100)
z_range = range(0, 1000, length = 100)

tau_x_grid = zeros(length(r_range), length(z_range))

for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        tau_x_grid[i,j] = compute_xray_tau2(quadtree, r, z, 0, xlumin, rg)
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, tau_x_grid', norm=LogNorm(1e2, 5e3))
plt.colorbar(cm ,ax=ax)
gcf()
