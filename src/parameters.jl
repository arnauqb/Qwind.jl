using YAML
export Parameters, change_parameter

struct Parameters{T<:AbstractFloat,S<:Flag,U<:Bool,V<:Int,W<:String}
    # T float, S flag, U bool, V int, W string
    M::T
    mdot::T
    spin::T
    disk_r_in::Union{T,W}
    disk_r_out::T
    disk_nr::V
    f_uv::Union{T,W}
    f_x::Union{T,W}
    radiation_grid_nr::Union{V,W}
    radiation_grid_nz::V
    z_xray::T
    z_disk::T
    mu_nucleon::T
    mu_electron::T
    vacuum_density::T
    wind_n_trajs::Union{V,W}
    wind_trajs_spacing::W
    wind_r_in::Union{T,W}
    wind_r_fi::T
    wind_z_0::T
    wind_n_0::Union{T,W}
    wind_v_0::Union{T,W}
    ic_K::Union{T,W}
    ic_alpha::T
    ic_use_precalculated::U
    integrator_r_min::T
    integrator_r_max::T
    integrator_z_min::T
    integrator_z_max::T
    n_iterations::V
    save_path::W
    # flags
    relativistic_flag::RelativisticFlag
    uv_opacity_flag::OpacityFlag
    xray_opacity_flag::OpacityFlag
    xray_scattering_flag::XRayScatteringFlag
    tau_uv_calculation_flag::TauUVCalculationFlag
    update_grid_flag::UpdateGridFlag
    fm_interp_method_flag::FMInterpolationFlag
    initial_conditions_flag::ICFlag
    # tolerances
    integrator_atol::T
    integrator_rtol::T
    disk_integral_atol::T
    disk_integral_rtol::T
    scattering_atol::T
    scattering_rtol::T
    function Parameters(;
        M,
        mdot,
        spin,
        disk_r_in,
        disk_r_out,
        disk_nr,
        f_uv,
        f_x,
        radiation_grid_nr,
        radiation_grid_nz,
        z_xray,
        z_disk,
        mu_nucleon,
        mu_electron,
        vacuum_density,
        wind_n_trajs,
        wind_trajs_spacing,
        wind_r_in,
        wind_r_fi,
        wind_z_0,
        wind_n_0,
        wind_v_0,
        ic_K,
        ic_alpha,
        ic_use_precalculated,
        integrator_r_min,
        integrator_r_max,
        integrator_z_min,
        integrator_z_max,
        n_iterations,
        save_path,
        # flag
        relativistic_flag,
        uv_opacity_flag,
        xray_opacity_flag,
        xray_scattering_flag,
        tau_uv_calculation_flag,
        update_grid_flag,
        fm_interp_method_flag,
        initial_conditions_flag,
        # tolerances
        integrator_atol,
        integrator_rtol,
        disk_integral_atol,
        disk_integral_rtol,
        scattering_atol,
        scattering_rtol,
    )
        return new{Float64,Flag,Bool,Int,String}(
            # floats
            M,
            mdot,
            spin,
            disk_r_in,
            disk_r_out,
            disk_nr,
            f_uv,
            f_x,
            radiation_grid_nr,
            radiation_grid_nz,
            z_xray,
            z_disk,
            mu_nucleon,
            mu_electron,
            vacuum_density,
            wind_n_trajs,
            wind_trajs_spacing,
            wind_r_in,
            wind_r_fi,
            wind_z_0,
            wind_n_0,
            wind_v_0,
            ic_K,
            ic_alpha,
            ic_use_precalculated,
            integrator_r_min,
            integrator_r_max,
            integrator_z_min,
            integrator_z_max,
            n_iterations,
            save_path,
            # flags
            relativistic_flag,
            uv_opacity_flag,
            xray_opacity_flag,
            xray_scattering_flag,
            tau_uv_calculation_flag,
            update_grid_flag,
            fm_interp_method_flag,
            initial_conditions_flag,
            # tolerances
            integrator_atol,
            integrator_rtol,
            disk_integral_atol,
            disk_integral_rtol,
            scattering_atol,
            scattering_rtol,
        )
    end
end

function Parameters(config::Dict)
    bhp = get(config, :black_hole, Dict())
    rp = get(config, :radiation, Dict())
    if get(rp, :relativistic, true)
        relativistic_flag = Relativistic()
    else
        relativistic_flag = NoRelativistic()
    end
    opacity_dict = Dict("thomson" => ThomsonOpacity(), "boost" => BoostOpacity())
    uv_opacity_flag = opacity_dict[get(rp, :uv_opacity, "thomson")]
    xray_opacity_flag = opacity_dict[get(rp, :xray_opacity, "thomson")]
    if get(rp, :xray_scattering, true)
        xray_scattering_flag = Scattering()
    else
        xray_scattering_flag = NoScattering()
    end
    tau_uv_dict =
        Dict("center" => TauUVCenter(), "disk" => TauUVDisk(), "no_tau_uv" => NoTauUV())
    tau_uv_calculation_flag = tau_uv_dict[get(rp, :tau_uv_calculation, "no_tau_uv")]
    update_grid_dict = Dict("average" => AverageGrid(), "replace" => ReplaceGrid())
    update_grid_flag = update_grid_dict[get(rp, :update_grid_method, "average")]
    fm_interp_method_flag =
        Dict("interpolate" => FMInterp(), "analytic" => FMNoInterp())[get(
            rp,
            :fm_interp_method_flag,
            "interpolate",
        )]

    np = get(config, :numerical_tolerances, Dict())
    ip = get(config, :initial_conditions, Dict())
    initial_conditions_flag =
        Dict("cak" => CAKMode(), "uniform" => UniformMode(), "ss" => SSMode())[get(
            ip,
            :mode,
            "cak",
        ),]
    intc = get(config, :integrator, Dict())

    return Parameters(
        # floats
        M = Float64(get(bhp, :M, 1e8)),
        mdot = Float64(get(bhp, :mdot, 0.5)),
        spin = Float64(get(bhp, :spin, 0.0)),
        disk_r_in = get(rp, :disk_r_in, 6.0),
        disk_r_out = get(rp, :disk_r_out, 1600.0),
        disk_nr = get(rp, :n_r, 1000),
        f_uv = get(rp, :f_uv, "auto"),
        f_x = get(rp, :f_x, 0.15),
        radiation_grid_nr = get(rp, :grid_nr, "auto"),
        radiation_grid_nz = get(rp, :grid_nz, 500),
        z_xray = get(rp, :z_xray, 0.0),
        z_disk = get(rp, :z_disk, 0.0),
        mu_nucleon = get(rp, :mu_nucleon, 0.61),
        mu_electron = get(rp, :mu_electron, 1.17),
        vacuum_density = get(rp, :vacuum_density, 1e2),

        # initial conditions
        wind_n_trajs = get(ip, :n_trajs, "auto"),
        wind_trajs_spacing = get(ip, :trajs_spacing, "log"),
        wind_r_in = get(ip, :r_in, 20.0),
        wind_r_fi = get(ip, :r_fi, 1500.0),
        wind_z_0 = get(ip, :z_0, 0.0),
        wind_n_0 = get(ip, :n_0, 1e8),
        wind_v_0 = get(ip, :v_0, 1e7),
        ic_K = get(ip, :K, 0.03),
        ic_alpha = get(ip, :alpha, 0.6),
        ic_use_precalculated = get(ip, :use_precalculated, true),

        # integrators
        integrator_r_min = get(intc, :integrator_r_min, 6.0),
        integrator_r_max = get(intc, :integrator_r_max, 5e4),
        integrator_z_min = get(intc, :integrator_z_min, 0.0),
        integrator_z_max = get(intc, :integrator_z_max, 5e4),
        n_iterations = get(intc, :n_iterations, 50),
        save_path = get(intc, :save_path, "./debug"),

        # flags
        relativistic_flag = relativistic_flag,
        uv_opacity_flag = uv_opacity_flag,
        xray_opacity_flag = xray_opacity_flag,
        xray_scattering_flag = xray_scattering_flag,
        tau_uv_calculation_flag = tau_uv_calculation_flag,
        update_grid_flag = update_grid_flag,
        fm_interp_method_flag = fm_interp_method_flag,
        initial_conditions_flag = initial_conditions_flag,

        # tolerances
        disk_integral_atol = Float64(get(np, :disk_integral_atol, 0)),
        disk_integral_rtol = Float64(get(np, :disk_integral_rtol, 1e-3)),
        integrator_atol = Float64(get(np, :integrator_atol, 1e-8)),
        integrator_rtol = Float64(get(np, :integrator_rtol, 1e-3)),
        scattering_atol = Float64(get(np, :scattering_atol, 0)),
        scattering_rtol = Float64(get(np, :scattering_rtol, 1e-3)),
    )
end

Parameters(config_file::String) =
    Parameters(YAML.load_file(config_file, dicttype = Dict{Symbol,Any}))

function change_parameter(parameters::Parameters, key, value)
    dd = Dict()
    for field in fieldnames(Parameters)
        if key == field
            dd[field] = value
        else
            dd[field] = getfield(parameters, field)
        end
    end
    return Parameters(; dd...)
end
