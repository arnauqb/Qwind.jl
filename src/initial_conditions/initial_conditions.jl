using Roots, CSV, DataFrames, Interpolations
using Optim: optimize, Brent
export InitialConditions,
    UniformIC,
    CAKIC,
    getz0,
    getrin,
    getrfi,
    getn0,
    getv0,
    getnlines,
    getl0,
    get_initial_conditions

getl0(ic::InitialConditions, r) = sqrt(r)
getrin(ic::InitialConditions) = ic.rin
getrfi(ic::InitialConditions) = ic.rfi
getnlines(ic::InitialConditions) = ic.nlines

struct UniformIC{T} <: InitialConditions{T}
    rin::T
    rfi::T
    nlines::Union{Int,String}
    z0::T
    n0::T
    v0::T
    trajs_spacing::String
end

function UniformIC(radiation, parameters)
    return UniformIC(
        parameters.wind_r_in,
        parameters.wind_r_fi,
        parameters.wind_n_trajs,
        parameters.wind_z_0,
        parameters.wind_n_0,
        parameters.wind_v_0,
        parameters.wind_trajs_spacing,
    )
end


getz0(ic::UniformIC, r) = ic.z0
getn0(ic::UniformIC, r) = ic.n0
getn0(ic::UniformIC, rad, r0) = ic.n0
getv0(ic::UniformIC, parameters::Parameters, r) = ic.v0
## CAK

struct CAKIC{T} <: InitialConditions{T}
    radiation::Radiation
    rin::T
    rfi::T
    nlines::Union{Int,String}
    z0::T
    K::Union{T,String}
    alpha::Union{T,String}
    trajs_spacing::String
    critical_points_df::DataFrame
    zc_interpolator::Any
    mdot_interpolator::Any
end

function CAKIC(radiation, parameters)
    M = parameters.M
    mdot = parameters.mdot
    if parameters.ic_use_precalculated
        filename = "critical_points_data/M_$(M)_mdot_$(mdot).csv"
        critical_points_df = CSV.read(joinpath(@__DIR__, filename), DataFrame)
    else
        rr, mdots, zcs = calculate_wind_mdots(radiation, parameters)
        critical_points_df = DataFrame(:r => rr, :mdot => mdots, :zc => zcs)
    end
    zc_interpolator = LinearInterpolation(
        critical_points_df[!, :r],
        critical_points_df[!, :zc],
        extrapolation_bc = Line(),
    )
    mdot_interpolator = LinearInterpolation(
        critical_points_df[!, :r],
        critical_points_df[!, :mdot],
        extrapolation_bc = Line(),
    )
    rin = max(parameters.wind_r_in, minimum(critical_points_df[!, :r]))
    rfi = min(parameters.wind_r_fi, maximum(critical_points_df[!, :r]))
    return CAKIC(
        radiation,
        rin,
        rfi,
        parameters.wind_n_trajs,
        parameters.wind_z_0,
        parameters.ic_K,
        parameters.ic_alpha,
        parameters.wind_trajs_spacing,
        critical_points_df,
        zc_interpolator,
        mdot_interpolator,
    )
end

getz0(ic::CAKIC, r0) = ic.z0
function getn0(ic::CAKIC, radiation::Radiation, parameters, r0)
    mdot = ic.mdot_interpolator(r0)
    return get_initial_density(radiation, parameters, r = r0, mdot = mdot)
end
getn0(model, r0) = getn0(model.ic, model.rad, model.parameters, r0)
getv0(ic::CAKIC, r0; mu_nucleon = 0.61) =
    compute_thermal_velocity(disk_temperature(ic.radiation.bh, r0), mu_nucleon)
getv0(model, r0) = getv0(model.ic, r0, mu_nucleon = model.parameters.mu_nucleon)

function get_initial_conditions(radiation, parameters)
    if parameters.initial_conditions_flag == CAKMode()
        return CAKIC(radiation, parameters)
    elseif parameters.initial_conditions_flag == UniformMode()
        return UniformMode(radiation, parameters)
    end
end
