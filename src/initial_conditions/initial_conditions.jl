using Roots, CSV, DataFrames
using Optim: optimize, Brent
export InitialConditions,
    UniformIC, CAKIC, getz0, getrin, getrfi, getn0, getv0, getnlines, getl0

getl0(ic::InitialConditions, r) = sqrt(r)
getrin(ic::InitialConditions) = ic.rin
getrfi(ic::InitialConditions) = ic.rfi
getnlines(ic::InitialConditions) = ic.nlines

struct UniformIC{T} <: InitialConditions{T}
    rin::T
    rfi::T
    nlines::Int
    z0::T
    n0::T
    v0::T
    logspaced::Bool
end

function UniformIC(radiation, black_hole, config)
    icc = config[:initial_conditions]
    if :launch_range in keys(icc)
        rin, rfi = icc[:launch_range]
    else
        rin = icc[:r_in]
        rfi = icc[:r_fi]
    end
    nlines = icc[:n_lines]
    if nlines != "auto"
        nlines = Int(nlines)
    end
    return UniformIC(
        rin,
        rfi,
        nlines,
        icc[:z_0],
        icc[:n_0],
        icc[:v_0] / C,
        icc[:log_spaced],
    )
end


getz0(ic::UniformIC, r) = ic.z0
getn0(ic::UniformIC, r) = ic.n0
getv0(ic::UniformIC, r) = ic.v0
## CAK

struct CAKIC{T} <: InitialConditions{T}
    radiation::Radiation
    bh::BlackHole
    rin::T
    rfi::T
    nlines::Union{Int,String}
    z0::T
    K::Union{T,String}
    alpha::Union{T,String}
    mu::T
    logspaced::Bool
    critical_points_df::DataFrame
end

function CAKIC(radiation, black_hole, config)
    icc = config[:initial_conditions]
    M = config[:black_hole][:M]
    mdot = config[:black_hole][:mdot]
    filename = "critical_points_data/M_$(M)_mdot_$(mdot).csv"
    critical_points_df =
        CSV.read(joinpath(@__DIR__, filename), DataFrame)
    if :launch_range in keys(icc)
        rin, rfi = icc[:launch_range]
    else
        rin = icc[:r_in]
        rfi = icc[:r_fi]
    end
    nlines = icc[:n_lines]
    if nlines != "auto"
        nlines = Int(nlines)
    end
    return CAKIC(
        radiation,
        black_hole,
        rin,
        rfi,
        nlines,
        icc[:z_0],
        icc[:K],
        icc[:alpha],
        icc[:mu],
        icc[:log_spaced],
        critical_points_df,
    )
end

getz0(ic::CAKIC, r0) = ic.z0
function getn0(ic::CAKIC, rt::RadiativeTransfer, bh::BlackHole, r0; K="auto")
    rv, ridx = findmin(abs.(ic.critical_points_df.r .- r0))
    zc = ic.critical_points_df.zc[ridx]
    mdot = ic.critical_points_df.mdot[ridx]
    if K == "auto"
        taux = compute_xray_tau(rt, rt.radiation.z_xray, r0, zc)
        density = get_density(rt.interpolator.density_grid, r0, zc)
        ξ = compute_ionization_parameter(rt.radiation, r0, zc, density, taux)
        K = compute_force_multiplier_k(1e10)
    end
    n = get_initial_density(rt, bh, r0, mdot; K = K, alpha = 0.6, mu = 0.5)
    return n
end
getn0(model, r0; K="auto") = getn0(model.ic, model.rt, model.bh, r0; K=K)
#getn0(ic::CAKIC, r0) = @. 10^(
#    -1.8110675134268326 * log10(r0)^2 + 3.788166202329894 * log10(r0) + 11.537466980651164
#) #cak_density(ic.radiation, ic.bh, r0, ic.K, ic.alpha)
getv0(ic::CAKIC, r0) = compute_thermal_velocity(disk_temperature(ic.bh, r0))

"""
Computes k(T) and alpha(T) from CAK 1975 paper.
"""
function compute_cak_K(T)
    if T < 30000
        return 0.0076
    elseif T < 40000
        return 0.0026
    else
        return 0.0021
    end
end
compute_cak_K(bh::BlackHole, r0) = compute_cak_K(disk_temperature(bh, r0))
function compute_cak_alpha(T)
    if T < 30000
        return 0.742
    elseif T < 40000
        return 0.737
    else
        return 0.811
    end
end
compute_cak_alpha(bh::BlackHole, r0) = compute_cak_alpha(disk_temperature(bh, r0))

export cak_surface_mloss,
    cak_density, cak_normalized_mdot, cak_characteristic_mloss, cak_nozzle_function

"Nozzle function defined in Pereyra et al. (2004) (Paper I)"
function cak_nozzle_function(radiation::QsosedRadiation, bh::BlackHole, z, r_0, alpha)
    x = z[1] / r_0
    T = disk_temperature(bh, r_0)
    b = compute_thermal_velocity(25e3) * C
    M = bh.M
    R0 = r_0 * bh.Rg
    numerator = (1 + x^2)^(1 / alpha)
    denom_1 = x / (1 + x^2)^(3 / 2)
    denom_2 = -SIGMA_E * SIGMA_SB * T^4 / (G * M * C) * R0^2
    denom_3 = -4 * b^2 * R0 * x / (G * M)
    if denom_1 + denom_2 + denom_3 <= 0
        return Inf
    end
    denom = (denom_1 + denom_2 + denom_3)^((1 - alpha) / alpha)
    return numerator / denom
end

function cak_characteristic_mloss(radiation::QsosedRadiation, bh::BlackHole, r_0, K, alpha)
    T = disk_temperature(bh, r_0)
    b = compute_thermal_velocity(25e3) * C
    M = bh.M
    R0 = r_0 * bh.Rg
    constant = alpha * (1 - alpha)^((1 - alpha) / alpha) / (b * SIGMA_E)
    term_1 = G * M / R0^2
    term_2 = (SIGMA_E * SIGMA_SB * T^4 * K * R0^2 / (C * G * M))^(1 / alpha)
    return constant * term_1 * term_2
end

function cak_normalized_mdot(radiation::QsosedRadiation, bh::BlackHole, r_0, alpha)
    f(z) = cak_nozzle_function(radiation, bh, z, r_0, alpha)
    mdot = optimize(f, 0, 2 * r_0, Brent()).minimum
    return mdot
end

function cak_surface_mloss(radiation::QsosedRadiation, bh::BlackHole, r_0, K, alpha)
    f(z) = cak_nozzle_function(radiation, bh, z, r_0, alpha)
    mdot = optimize(f, 0, 2 * r_0, Brent()).minimum
    Sigma = cak_characteristic_mloss(radiation, bh, r_0, K, alpha) * mdot
    return Sigma
end

function cak_density(
    radiation::QsosedRadiation,
    bh::BlackHole,
    r_0,
    K = "auto",
    alpha = "auto",
)
    (K == "auto") && (K = compute_cak_K(bh, r_0))
    (alpha == "auto") && (alpha = compute_cak_alpha(bh, r_0))
    Sigma = cak_surface_mloss(radiation, bh, r_0, K, alpha)
    T = disk_temperature(bh, r_0)
    b = compute_thermal_velocity(T) * C
    return Sigma / b / M_P
end

