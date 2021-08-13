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

function UniformIC(radiation, radiative_transfer, black_hole, config)
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
getn0(ic::UniformIC, rt::RadiativeTransfer, bh::BlackHole, r0) = ic.n0
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
    logspaced::Bool
    critical_points_df::DataFrame
end

function CAKIC(radiation, radiative_transfer, black_hole, config)
    icc = config[:initial_conditions]
    M = config[:black_hole][:M]
    mdot = config[:black_hole][:mdot]
    if icc[:use_precalculated]
        filename = "critical_points_data/M_$(M)_mdot_$(mdot).csv"
        critical_points_df = CSV.read(joinpath(@__DIR__, filename), DataFrame)
    else
        rr, mdots, zcs = calculate_wind_mdots(radiative_transfer, black_hole)
        critical_points_df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    end
    if :launch_range in keys(icc)
        rin, rfi = icc[:launch_range]
    else
        rin = icc[:r_in]
        rfi = icc[:r_fi]
    end
    rin = max(rin, minimum(critical_points_df[!, :r]))
    rfi = min(rfi, maximum(critical_points_df[!, :r]))
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
        icc[:log_spaced],
        critical_points_df,
    )
end

getz0(ic::CAKIC, r0) = ic.z0
function getn0(ic::CAKIC, rt::RadiativeTransfer, bh::BlackHole, r0)
    rv, ridx = findmin(abs.(ic.critical_points_df.r .- r0))
    zc = ic.critical_points_df.zc[ridx]
    mdot = ic.critical_points_df.mdot[ridx]
    if ic.K == "auto"
        taux = compute_xray_tau(rt, rt.radiation.z_xray, r0, zc)
        density = get_density(rt.interpolator.density_grid, r0, zc)
        ξ = compute_ionization_parameter(rt.radiation, r0, zc, density, taux)
        K = compute_force_multiplier_k(ξ)
    end
    n = get_initial_density(rt, bh, r0, mdot; K = ic.K, alpha = ic.alpha)
    return n
end
getn0(model, r0) = getn0(model.ic, model.rt, model.bh, r0)
getv0(ic::CAKIC, r0) = compute_thermal_velocity(disk_temperature(ic.bh, r0))
