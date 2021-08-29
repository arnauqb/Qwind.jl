using Roots, CSV, DataFrames, Interpolations
using LsqFit: curve_fit
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

function UniformIC(radiation, config)
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
getn0(ic::UniformIC, rad, r0) = ic.n0
getv0(ic::UniformIC, r) = ic.v0
## CAK

struct CAKIC{T} <: InitialConditions{T}
    radiation::Radiation
    rin::T
    rfi::T
    nlines::Union{Int,String}
    z0::T
    K::Union{T,String}
    alpha::Union{T,String}
    logspaced::Bool
    critical_points_df::DataFrame
    zc_interpolator::Any
    log_density_interpolator::Any
end

function CAKIC(radiation, config)
    icc = config[:initial_conditions]
    M = config[:black_hole][:M]
    mdot = config[:black_hole][:mdot]
    if icc[:use_precalculated]
        filename = "critical_points_data/M_$(M)_mdot_$(mdot).csv"
        critical_points_df = CSV.read(joinpath(@__DIR__, filename), DataFrame)
    else
        rr, mdots, zcs = calculate_wind_mdots(radiation)
        critical_points_df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    end
    if :launch_range in keys(icc)
        rin, rfi = icc[:launch_range]
    else
        rin = icc[:r_in]
        rfi = icc[:r_fi]
    end
    rs = critical_points_df[!, :r]
    mdots = critical_points_df[!, :mdot]
    zc_interpolator = LinearInterpolation(critical_points_df[!, :r], critical_points_df[!, :zc], extrapolation_bc=Line())
    inter_mode = get(icc, :interapolation_mode, "interpolation")
    n0s = [get_initial_density(radiation, r=r, mdot=mdot, K=0.03, alpha=icc[:alpha]) for (r, mdot) in zip(rs, mdots)]
    if inter_mode == "fit"
        @. tofit(x,p) = p[1] / x^2 + p[2] / x + p[3] + p[4] * x + p[5] * x^2
        p0 = ones(5)
        fit = curve_fit(tofit, log10.(rs), log10.(n0s), p0)
        fitted(x) = tofit(x, fit.param)
        log_density_interpolator = fitted
    elseif inter_mode == "interpolation"
        log_density_interpolator = LinearInterpolation(log10.(rs), log10.(n0s), extrapolation_bc=Line())
    else
        error("inter mode not recognized")
    end
    rin = max(rin, minimum(critical_points_df[!, :r]))
    rfi = min(rfi, maximum(critical_points_df[!, :r]))
    nlines = icc[:n_lines]
    if nlines != "auto"
        nlines = Int(nlines)
    end
    return CAKIC(
        radiation,
        rin,
        rfi,
        nlines,
        icc[:z_0],
        icc[:K],
        icc[:alpha],
        icc[:log_spaced],
        critical_points_df,
        zc_interpolator,
        log_density_interpolator,
    )
end

getz0(ic::CAKIC, r0) = ic.z0
function getn0(ic::CAKIC, radiation::Radiation, r0)
    n_k003 = 10 ^ ic.log_density_interpolator(log10(r0))
    zc = ic.zc_interpolator(r0)
    mdot = nozzle_function(radiation, zc, r=r0)
    return get_initial_density(radiation, r=r0, mdot=mdot, K=ic.K)
    #K = ic.K
    #if K == "auto"
    #    taux = compute_tau_xray(radiation, r=r0, z=zc)
    #    density = get_density(radiation.wi.density_grid, r0, zc)
    #    ξ = compute_ionization_parameter(radiation, r0, zc, 0.0, 0.0, density, taux)
    #    K = compute_force_multiplier_k(ξ, FMNoInterp())
    #end
    #return n_k003 * (K / 0.03)^(1 / ic.alpha)
end
getn0(model, r0) = getn0(model.ic, model.rad, r0)
getv0(ic::CAKIC, r0) = compute_thermal_velocity(disk_temperature(ic.radiation.bh, r0))
getv0(model, r0) = getv0(model.ic, r0)
