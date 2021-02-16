using Roots
using Optim: optimize, Brent
export InitialConditions,
    UniformIC, CAKIC, getz0, getrin, getrfi, getn0, getv0, getnlines, getl0

getl0(ic::InitialConditions, r) = sqrt(r)
getrin(ic::InitialConditions) = ic.rin
getrfi(ic::InitialConditions) = ic.rfi
getnlines(ic::InitialConditions) = ic.nlines

struct UniformIC <: InitialConditions
    rin::Float64
    rfi::Float64
    nlines::Int64
    z0::Float64
    n0::Float64
    v0::Float64
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
    return UniformIC(
        rin,
        rfi,
        Int(icc[:n_lines]),
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

struct CAKIC <: InitialConditions
    radiation::Radiation
    bh::BlackHole
    rin::Float64
    rfi::Float64
    nlines::Union{Int64, String}
    z0::Float64
    K::Float64
    alpha::Float64
    mu::Float64
    logspaced::Bool
end

function CAKIC(radiation, black_hole, config)
    icc = config[:initial_conditions]
    if :launch_range in keys(icc)
        rin, rfi = icc[:launch_range]
    else
        rin = icc[:r_in]
        rfi = icc[:r_fi]
    end
    return CAKIC(
        radiation,
        black_hole,
        rin,
        rfi,
        Int(icc[:n_lines]),
        icc[:z_0],
        icc[:K],
        icc[:alpha],
        icc[:mu],
        icc[:log_spaced],
    )
end

getz0(ic::CAKIC, r0) = ic.z0
getn0(ic::CAKIC, r0) = cak_density(ic.radiation, ic.bh, r0, ic.K)
getv0(ic::CAKIC, r0) = compute_thermal_velocity(disk_temperature(ic.bh, r0))

const ALPHA = 0.6
export cak_surface_mloss,
    cak_density, cak_normalized_mdot, cak_characteristic_mloss, cak_nozzle_function

"Nozzle function defined in Pereyra et al. (2004) (Paper I)"
function cak_nozzle_function(radiation::QsosedRadiation, bh::BlackHole, z, r_0)
    x = z[1] / r_0
    T = disk_temperature(bh, r_0)
    b = compute_thermal_velocity(25e3) * C
    M = bh.M
    R0 = r_0 * bh.Rg
    numerator = (1 + x^2)^(1 / ALPHA)
    denom_1 = x / (1 + x^2)^(3 / 2)
    denom_2 = -SIGMA_E * SIGMA_SB * T^4 / (G * M * C) * R0^2
    denom_3 = -4 * b^2 * R0 * x / (G * M)
    denom = (denom_1 + denom_2 + denom_3)^((1 - ALPHA) / ALPHA)
    return numerator / denom
end

function cak_characteristic_mloss(radiation::QsosedRadiation, bh::BlackHole, r_0, K = 0.03)
    T = disk_temperature(bh, r_0)
    b = compute_thermal_velocity(25e3) * C
    M = bh.M
    R0 = r_0 * bh.Rg
    constant = ALPHA * (1 - ALPHA)^((1 - ALPHA) / ALPHA) / (b * SIGMA_E)
    term_1 = G * M / R0^2
    term_2 = (SIGMA_E * SIGMA_SB * T^4 * K * R0^2 / (C * G * M))^(1 / ALPHA)
    return constant * term_1 * term_2
end

function cak_normalized_mdot(radiation::QsosedRadiation, bh::BlackHole, r_0)
    f(z) = cak_nozzle_function(radiation, bh, z, r_0)
    mdot = optimize(f, 0, 2 * r_0, Brent()).minimum
    return mdot
end

function cak_surface_mloss(radiation::QsosedRadiation, bh::BlackHole, r_0, K = 0.03)
    f(z) = cak_nozzle_function(radiation, bh, z, r_0)
    mdot = optimize(f, 0, 2 * r_0, Brent()).minimum
    Sigma = cak_characteristic_mloss(radiation, bh, r_0) * mdot
    return Sigma
end

function cak_density(radiation::QsosedRadiation, bh::BlackHole, r_0, K = nothing)
    if K === nothing
        K = 0.03
    end
    Sigma = cak_surface_mloss(radiation, bh, r_0, K)
    T = disk_temperature(bh, r_0)
    b = compute_thermal_velocity(T) * C
    return Sigma / b / M_P
end

