export AdaptiveMesh, compute_xray_tau, compute_uv_tau, update_radiative_transfer

struct AdaptiveMesh <: RadiativeTransfer
    radiation::Radiation
    windkdtree::Union{KDTree,Nothing}
    quadtree::Union{Cell,Nothing}
    vacuum_density::Float64
    atol::Float64
    rtol::Float64
    cell_min_size::Float64
    n_timesteps::Int64
end

function AdaptiveMesh(
    radiation::Radiation,
    integrators;
    vacuum_density = 1e2,
    atol = 1e-4,
    rtol = 1e-3,
    cell_min_size = 0.0001,
    n_timesteps = 10000,
)
    if integrators === nothing
        windkdtree = nothing
    else
        windkdtree = create_wind_kdtree(integrators, n_timesteps)
    end
    quadtree = create_and_refine_quadtree(
        windkdtree,
        Rg = radiation.Rg,
        atol = atol,
        rtol = rtol,
        cell_min_size = cell_min_size,
        vacuum_density = vacuum_density,
    )
    return AdaptiveMesh(
        radiation,
        windkdtree,
        quadtree,
        vacuum_density,
        atol,
        rtol,
        cell_min_size,
        n_timesteps,
    )
end

function AdaptiveMesh(radiation::Radiation, config::Dict)
    rtc = config["radiative_transfer"]
    return AdaptiveMesh(
        radiation,
        nothing,
        vacuum_density = rtc["vacuum_density"],
        atol = rtc["atol"],
        rtol = rtc["rtol"],
        cell_min_size = rtc["cell_min_size"],
        n_timesteps = rtc["n_timesteps"],
    )
end

function update_radiative_transfer(rt::AdaptiveMesh, integrators)
    @info "Updating radiative transfer... "
    return AdaptiveMesh(
        rt.radiation,
        integrators,
        vacuum_density = rt.vacuum_density,
        atol = rt.atol,
        rtol = rt.rtol,
        cell_min_size = rt.cell_min_size,
        n_timesteps = rt.n_timesteps,
    )
end


function compute_xray_tau_leaf(
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
    density = leaf.data
    deltatau = density * deltad
    xi0 = xray_luminosity / (density * d^2)
    f(t) = t - log(xi0) + taux0 + min(40, deltatau * compute_xray_opacity(exp(t)))
    if f(20) < 0
        xi = xi0
    elseif f(-20) > 0
        xi = 1e-20
    else
        t = find_zero(f, (-20, 20), Bisection(), atol = 0, rtol = 0.1)
        xi = exp(t)
    end
    taux = compute_xray_opacity(xi) * deltatau
    return taux
end

compute_xray_tau_leaf(adaptive_mesh::AdaptiveMesh, leaf::Cell, point, intersection, taux0) =
    compute_xray_tau_leaf(
        adaptive_mesh.quadtree,
        leaf,
        point,
        intersection,
        taux0,
        adaptive_mesh.radiation.xray_luminosity,
        adaptive_mesh.Rg,
    )

function compute_xray_tau(quadtree::Cell, r, z, r0, z0, xray_luminosity, Rg)
    z = max(z, getzmin(quadtree))
    point1 = [r0, z0]
    #coords_list = [point1]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r, z]
    point2leaf = findleaf(quadtree, point2)
    taux = 0.0
    if point1leaf == point2leaf
        #push!(coords_list, point2)
        taux = compute_xray_tau_leaf(
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
    intersection = compute_cell_intersection(point1leaf, point1, point1, point2)
    taux = compute_xray_tau_leaf(
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
        intersection = compute_cell_intersection(currentleaf, currentpoint, point1, point2)
        #push!(coords_list, intersection)
        taux += compute_xray_tau_leaf(
            quadtree,
            currentleaf,
            currentpoint,
            intersection,
            taux,
            xray_luminosity,
            Rg,
        )
        if taux > 40
            return 40.0
        end
        currentpoint = intersection
        currentleaf = findleaf(quadtree, currentpoint)
        if currentleaf == previousleaf
            break
        end
        previousleaf = copy(currentleaf)
    end
    if currentpoint[2] < point2[2]
        taux += compute_xray_tau_leaf(
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
    return taux#, coords_list
end

compute_xray_tau(adaptive_mesh::AdaptiveMesh, r, z, r0 = 0, z0 = 0) = compute_xray_tau(
    adaptive_mesh.quadtree,
    r,
    z,
    r0,
    z0,
    adaptive_mesh.radiation.xray_luminosity,
    adaptive_mesh.radiation.Rg,
)

function compute_uv_tau_leaf(pointleaf::Cell, point, intersection, Rg)
    deltad = evaluate(Euclidean(), point, intersection) * Rg
    return deltad * SIGMA_T * pointleaf.data
end

"""
Compute the UV optical depth from a disc patch located at (r_d, phi_d),
until a gas element at (r,z).
"""
function compute_uv_tau(quadtree::Cell, rd, r, z, Rg, maxtau = 50, z0 = 0.0)
    rd > r ? backwards = true : backwards = false
    point1 = [rd, z0]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r, z]
    point2leaf = findleaf(quadtree, point2)
    if point1leaf == point2leaf
        tauuv = compute_uv_tau_leaf(point1leaf, point1, point2, Rg)
        return tauuv
    end
    currentleaf = point1leaf
    currentpoint = point1
    previousleaf = copy(currentleaf)
    tauuv = 0.0
    #intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    #tauuv = compute_tauuv_leaf(point1, intersection, point1leaf, wind)
    #currentpoint = intersection
    #backwards && (currentpoint[1] -= 1e-8)
    #currentleaf = findleaf(wind.quadtree, currentpoint)
    #previousleaf = copy(currentleaf)
    while (currentleaf != point2leaf) && (currentpoint[2] <= point2[2])
        intersection = compute_cell_intersection(currentleaf, currentpoint, point1, point2)
        tauuv += compute_uv_tau_leaf(currentleaf, currentpoint, intersection, Rg)
        tauuv >= maxtau && return tauuv
        currentpoint = intersection
        if currentpoint == point2
            return tauuv
        end
        backwards && (currentpoint[1] -= 1e-10)
        currentleaf = findleaf(quadtree, currentpoint)
        if currentleaf == previousleaf
            break
        end
        previousleaf = copy(currentleaf)
    end
    if currentpoint[2] < point2[2]
        tauuv += compute_uv_tau_leaf(currentleaf, currentpoint, point2, Rg)
    end
    return tauuv
end

compute_uv_tau(rt::AdaptiveMesh, rd, r, z, maxtau = 50, z0 = 0.0) =
    compute_uv_tau(rt.quadtree, rd, r, z, rt.radiation.Rg, maxtau, z0)
