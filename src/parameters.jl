using YAML
export Parameters

struct Parameters{T<:AbstractFloat,S<:Flag,U<:Bool,V<:Int,W<:String}
    # T float, S flag, U bool, V int, W string
    M::T
    mdot::T
    spin::T
    disk_r_in::T
    disk_r_out::T
    disk_nr::V
    radiation_grid_nr::Union{V,W}
    radiation_grid_nz::V
    z_xray::T
    z_disk::T
    mu_nucleon::T
    mu_electron::T
    vacuum_density::T
    wind_r_min::T
    wind_r_max::T
    wind_z_min::T
    wind_z_max::T
    wind_n_trajs::Union{V,W}
    wind_trajs_spacing::W
    wind_r_in::T
    wind_r_fi::T
    wind_z_0::T
    ic_K::T
    ic_alpha::T
    ic_use_precalculated::U
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
        radiation_grid_nr,
        radiation_grid_nz,
        z_xray,
        z_disk,
        mu_nucleon,
        mu_electron,
        vacuum_density,
        wind_r_min,
        wind_r_max,
        wind_z_min,
        wind_z_max,
        wind_n_trajs,
        wind_trajs_spacing,
        wind_r_in,
        wind_r_fi,
        wind_z_0,
        ic_K,
        ic_alpha,
        ic_use_precalculated,
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
        return new{typeof(M), Flag, Bool, Int, String}(
            # floats
            M,
            mdot,
            spin,
            disk_r_in,
            disk_r_out,
            disk_nr,
            radiation_grid_nr,
            radiation_grid_nz,
            z_xray,
            z_disk,
            mu_nucleon,
            mu_electron,
            vacuum_density,
            wind_r_min,
            wind_r_max,
            wind_z_min,
            wind_z_max,
            wind_n_trajs,
            wind_trajs_spacing,
            wind_r_in,
            wind_r_fi,
            wind_z_0,
            ic_K,
            ic_alpha,
            ic_use_precalculated,
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
    rp = config[:radiation]
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

    np = config[:numerical_tolerances]
    ip = config[:initial_conditions]
    initial_conditions_flag =
        Dict("cak" => CAKMode(), "uniform" => UniformMode())[get(ip, :mode, "cak")]

    boundaries = get(config, :wind_boundaries, Dict())

    return Parameters(
        # floats
        M = Float64(config[:black_hole][:M]),
        mdot = Float64(config[:black_hole][:mdot]),
        spin = Float64(config[:black_hole][:spin]),
        disk_r_in = get(rp, :disk_r_in, 6.0),
        disk_r_out = get(rp, :disk_r_in, 1600.0),
        disk_nr = get(rp, :n_r, 1000),
        radiation_grid_nr = get(rp, :nr, "auto"),
        radiation_grid_nz = get(rp, :nz, 500),
        z_xray = get(rp, :z_xray, 0.0),
        z_disk = get(rp, :z_disk, 0.0),
        mu_nucleon = get(rp, :mu_nucleon, 0.61),
        mu_electron = get(rp, :mu_electron, 1.17),
        vacuum_density = get(rp, :vacuum_density, 1e2),
        wind_r_min = get(boundaries, :wind_r_min, 6.0),
        wind_r_max = get(boundaries, :wind_r_max, 5e4),
        wind_z_min = get(boundaries, :wind_z_min, 0.0),
        wind_z_max = get(boundaries, :wind_z_max, 5e4),

        # initial conditions
        wind_n_trajs = get(ip, :n_lines, "auto"),
        wind_trajs_spacing = get(ip, :trajs_spacing, "log"),
        wind_r_in = get(ip, :r_in, 20.0),
        wind_r_fi = get(ip, :r_fi, 1500.0),
        wind_z_0 = get(ip, :z0, 0.0),
        ic_K = get(ip, :K, 0.03),
        ic_alpha = get(ip, :alpha, 0.6),
        ic_use_precalculated = get(ip, :use_precalculated, true),

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
