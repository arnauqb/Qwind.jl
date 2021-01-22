export AdaptiveMesh,
    compute_xray_tau,
    compute_uv_tau,
    compute_disc_radiation_field,
    update_radiative_transfer
abstract type IntegrationType end
struct IntegrationFromStreamline <: IntegrationType end
struct IntegrationFromCenter <: IntegrationType end

struct AdaptiveMesh <: RadiativeTransfer
    radiation::Radiation
    windkdtree::Union{WindKDTree,Nothing}
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
    coords_list = [point1]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r, z]
    point2leaf = findleaf(quadtree, point2)
    taux = 0.0
    if point1leaf == point2leaf
        push!(coords_list, point2)
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
    push!(coords_list, intersection)
    currentleaf = findleaf(quadtree, currentpoint)
    previousleaf = copy(currentleaf)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentleaf, currentpoint, point1, point2)
        push!(coords_list, intersection)
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
    push!(coords_list, point2)
    return taux, coords_list
end

compute_xray_tau(adaptive_mesh::RadiativeTransfer, r, z, r0 = 0, z0 = 0) = compute_xray_tau(
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

"""
Radiation force integrand for the RE model.
"""
function radiation_force_integrand!(
    adaptive_mesh::AdaptiveMesh,
    integration_type::IntegrationFromCenter,
    v,
    rd,
    phid,
    r,
    z,
)
    nt = disk_nt_rel_factors(adaptive_mesh.radiation, rd)
    tauuv = compute_uv_tau(adaptive_mesh.quadtree, rd, r, z, adaptive_mesh.radiation.Rg)
    fuv, mdot = get_fuv_mdot(adaptive_mesh.radiation, rd)
    r_projection = (r - rd * cos(phid))
    delta_sq = (r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    common_projection = 1.0 / (rd^2 * delta_sq^2)
    v[:] = exp(-tauuv) * fuv * mdot * nt * common_projection * [r_projection, z]
end

function radiation_force_integrand!(
    adaptive_mesh::AdaptiveMesh,
    integration_type::IntegrationFromStreamline,
    v,
    p,
    psi,
    r,
    z,
)
    cosψ = cos(psi)
    rd = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    if rd <= 6.0 || rd > 1600.0
        v[1] = 0.0
        v[2] = 0.0
        return
    end
    nt = disk_nt_rel_factors(adaptive_mesh.radiation, rd)
    fuv, mdot = get_fuv_mdot(adaptive_mesh.radiation, rd)
    delta2 = p^2 + z^2
    tauuv = compute_uv_tau(
        adaptive_mesh.quadtree,
        rd,
        r,
        z,
        adaptive_mesh.radiation.Rg,
    )
    tauuv = tauuv / sqrt((r - rd)^2 + z^2) * sqrt(delta2)
    common_part = nt * mdot * fuv * p / delta2^2 / rd^3 * exp(-tauuv)
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

#function integrate_fromstreamline(
#    r,
#    z,
#    wind::WindStruct;
#    include_tau_uv = true,
#    maxevals = 0,
#)
#    rgsigma = wind.bh.R_g * SIGMA_T
#    maxtau = 20 / rgsigma
#    r_max = wind.config["disk"]["outer_radius"]
#    xmin = (0.0, 0.0)
#    #xmax = (r+r_max, pi)
#    #xmax = (r+r_max, pi)#
#    xmax = (500, pi)#(min(max(10, z * 100), 1000), pi)
#    if include_tau_uv
#        (val, err) = hcubature(
#            2,
#            (x, v) -> integrate_fromstreamline_kernel(
#                v,
#                x[1],
#                x[2],
#                r,
#                z,
#                wind,
#                r_max,
#                rgsigma,
#                maxtau,
#            ),
#            xmin,
#            xmax,
#            reltol = wind.config["radiation"]["integral_rtol"],
#            abstol = 0, #abs(maximum(grav)/100),
#            maxevals = maxevals,
#            error_norm = Cubature.L1,
#        )
#    else
#        (val, err) = hcubature(
#            2,
#            (x, v) -> integrate_fromstreamline_notau_kernel(
#                v,
#                x[1],
#                x[2],
#                r,
#                z,
#                wind,
#                r_max,
#            ),
#            xmin,
#            xmax,
#            reltol = wind.config["radiation"]["integral_rtol"],
#            abstol = 0.0,
#            maxevals = maxevals,
#            error_norm = Cubature.L1,
#        )
#    end
#    return [2 * z, 2 * z^2] .* val
#end



"""
Integrates the radiation acceleration integrand on the disc.

# Parameters
- r : radial coordinate [Rg]
- z : height coordinate [Rg]
-
"""
function integrate_radiation_force_integrand(
    radiative_transfer::AdaptiveMesh,
    r,
    z,
    rmin,
    rmax = 1600;
    phi_min = 0.0,
    phi_max = π,
    atol = 0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 50000,
    zmax_fromstreamline = 1e-4,
)
    if z < zmax_fromstreamline
        return [0.0, 0.0], [0.0, 0.0]
        #integration_type = IntegrationFromStreamline()
    else
        integration_type = IntegrationFromCenter()
    end
    f(x, v) = radiation_force_integrand!(
        radiative_transfer,
        integration_type,
        v,
        x[1],
        x[2],
        r,
        z,
    )
    return hcubature(
        2,
        f,
        (rmin, phi_min),
        (rmax, phi_max),
        abstol = atol,
        reltol = rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
end

"""
Computes the disc radiation field at the point (r,z) by performing
an integral over the disc.

# Parameters
- r: radial coordinate [Rg]
- z: height coordinate [Rg]
- radiative_efficiency: accretion radiative efficiency
"""
function compute_disc_radiation_field(
    radiative_transfer::AdaptiveMesh,
    r,
    z;
    rmax = 1600,
    atol = 0,
    rtol = 5e-3,
    norm = Cubature.INDIVIDUAL,
    maxevals = 10000,
)
    #println("r : $r,\t z : $z")
    res, err = integrate_radiation_force_integrand(
        radiative_transfer,
        r,
        z,
        radiative_transfer.radiation.isco,
        rmax,
        atol = atol,
        rtol = rtol,
        norm = norm,
        maxevals = maxevals,
    )
    radiation_constant = compute_radiation_constant(radiative_transfer.radiation)
    force = z * radiation_constant .* res
    #println("force $force")
    return force
end
