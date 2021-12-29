using Roots, CSV, DataFrames, Interpolations
using Optim: optimize, Brent
export InitialConditions,
    UniformIC,
    CAKIC,
    getrin,
    getrfi,
    getn0,
    getz0,
    getv0,
    getnlines,
    getl0,
    get_initial_conditions

getl0(ic::InitialConditions, r) = sqrt(r)
getz0(ic::InitialConditions, r) = ic.z0
getv0(bh, ic::InitialConditions, r; mu_nucleon = 0.61) = ic.v0
getn0(ic::InitialConditions; rad, wind, parameters, r) = ic.n0
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

struct SSIC{T} <: InitialConditions{T}
    rin::T
    rfi::T
    nlines::Union{Int,String}
    z0::T
    n0::Union{T,String}
    v0::Union{T,String}
    trajs_spacing::String
end
function SSIC(radiation, parameters)
    return SSIC(
        parameters.wind_r_in,
        parameters.wind_r_fi,
        parameters.wind_n_trajs,
        parameters.wind_z_0,
        parameters.wind_n_0,
        parameters.wind_v_0,
        parameters.wind_trajs_spacing,
    )
end
function getn0(ic::SSIC; radiation, wind, parameters, r)
    return compute_ss_number_density(radiation.bh, r; mu_nucleon = parameters.mu_nucleon)
end

## CAK

struct CAKIC{T} <: InitialConditions{T}
    radiation::Radiation
    rin::T
    rfi::T
    nlines::Union{Int,String}
    z0::T
    v0::Union{T,String}
    K::Union{T,String}
    alpha::Union{T,String}
    trajs_spacing::String
    critical_points_df::DataFrame
    zc_interpolator::Any
    mdot_interpolator::Any
end

function CAKIC(radiation::Radiation, parameters::Parameters, wind)
    M = parameters.M
    mdot = parameters.mdot
    if parameters.ic_use_precalculated
        mdot_string = @sprintf "%.4f" mdot
        filename = "critical_points_data/M_$(M)_mdot_$(mdot_string).csv"
        critical_points_df = CSV.read(joinpath(@__DIR__, filename), DataFrame)
    else
        rr, mdots, zcs = calculate_wind_mdots(
            wind.density_grid,
            wind.grid_iterator,
            parameters,
            radiation,
            rmax=parameters.wind_r_fi
        )
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
    if parameters.wind_r_in == "warm_radius"
        wind_r_in = radiation.qsosed_model.warm.radius
    else
        wind_r_in = parameters.wind_r_in
    end
    rin = max(wind_r_in, minimum(critical_points_df[!, :r]))
    rfi = min(parameters.wind_r_fi, maximum(critical_points_df[!, :r]))
    return CAKIC(
        radiation,
        rin,
        rfi,
        parameters.wind_n_trajs,
        parameters.wind_z_0,
        parameters.wind_v_0,
        parameters.ic_K,
        parameters.ic_alpha,
        parameters.wind_trajs_spacing,
        critical_points_df,
        zc_interpolator,
        mdot_interpolator,
    )
end
function CAKIC(radiation::Radiation, config::Dict)
    parameters = Parameters(config)
    wind = Wind(parameters)
    return CAKIC(radiation, parameters, wind)
end

function getn0(ic::CAKIC; radiation, wind, parameters, r)
    mdot = ic.mdot_interpolator(r)
    zc = ic.zc_interpolator(r)
    return get_initial_density(radiation, wind, parameters, r = r, mdot = mdot, zc = zc)
end

function getv0(bh, ic::Union{CAKIC,SSIC}, r0; mu_nucleon = 0.61)
    if ic.v0 == "thermal"
        compute_thermal_velocity(disk_temperature(bh, r0), mu_nucleon)
    else
        return ic.v0
    end
end

# model calls
getv0(model, r0) = getv0(model.bh, model.ic, r0, mu_nucleon = model.parameters.mu_nucleon)
getn0(model, r0) = getn0(
    model.ic,
    radiation = model.rad,
    wind = model.wind,
    parameters = model.parameters,
    r = r0,
)

function get_initial_conditions(radiation, parameters, wind)
    if parameters.initial_conditions_flag == CAKMode()
        return CAKIC(radiation, parameters, wind)
    elseif parameters.initial_conditions_flag == UniformMode()
        return UniformIC(radiation, parameters)
    elseif parameters.initial_conditions_flag == SSMode()
        return SSIC(radiation, parameters)
    else
        throw()
    end
end
