using DifferentialEquations
using Sundials
export Parameters

struct Parameters
    line::Streamline
    bh_mass::Float64
    r_in::Float64
    r_x::Float64
    xray_luminosity::Float64
    v_th::Float64
    r_min::Float64
    r_max::Float64
    z_min::Float64
    z_max::Float64
    integrand::Any
    radiation_mode::RadiationMode
end

out_of_domain(params::Parameters) =
    !(
        (params.r_min <= r(params.line) <= params.r_max) &&
        (params.z_min <= z(params.line) <= params.z_max)
    )


"Updates the density of the streamline giving its current position and velocity,
using mass conservation."
function compute_density(r, z, v_t, line::Streamline)
    @assert r >= 0
    @assert z >= 0
    d = sqrt(r^2 + z^2)
    radial = (r_0(line) / d)^2
    v_ratio = v_0(line) / v_t
    n = number_density_0(line) * radial * v_ratio
    return n
end

function compute_initial_acceleration(r, z, v_r, v_z, params::Parameters)
    u = [r, z, v_r, v_z]
    du = [v_r, v_z, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(r, z, params.bh_mass)
    radiation_acceleration = radiation_calculator(du, u, params, params.radiation_mode)
    centrifugal_term = params.line.angular_momentum^2 / r^3
    isnan(centrifugal_term) ? centrifugal_term = 0 : nothing
    a_r = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    return [a_r, a_z]
end

function radiation_calculator(du, u, p, mode::RE)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = du
    v_t = sqrt(v_r^2 + v_z^2)
    a_t = sqrt(a_r^2 + a_z^2)
    dv_dr = a_t / v_t
    density = compute_density(r, z, v_t, p.line)
    tau_x = compute_xray_tau(r, z, density, p, RE())
    tau_uv = compute_uv_tau(r, z, density, p, RE())
    ξ = compute_ionization_parameter(r, z, density, tau_x, p.xray_luminosity)
    tau_eff = compute_tau_eff(density, dv_dr, p.v_th)
    force_multiplier = compute_force_multiplier(tau_eff, ξ)
    force_radiation =
        compute_radiation_force(p.integrand, r, z, force_multiplier, p.bh_mass)
    return force_radiation
end

radiation_calculator(du, u, p, mode::NoRad) = [0.0, 0.0]
radiation_calculator(du, u, p, mode::ConstantRad) = [0.0, 400]

function termination_condition(u, t, integrator)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = integrator.du
    line = integrator.p.line
    in_grid =
        (integrator.p.r_min <= r <= integrator.p.r_max) &&
        (integrator.p.z_min <= z <= integrator.p.z_max)
    out_of_grid_condition = !in_grid
    return out_of_grid_condition
end

function affect!(integrator)
    if escaped(integrator.p.line, integrator.p.bh_mass)
        print(" \U1F4A8")
    elseif out_of_domain(integrator.p)
        print(" \U2753")
    else
        print(" \U1F4A5")
    end
    terminate!(integrator)
end

function save(u, t, integrator)
    r, z, v_r, v_z = u
    save_position!(integrator.p.line, u, t)
end

function residual!(out, du, u, p, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    gravitational_acceleration = compute_gravitational_acceleration(r, z, p.bh_mass)
    radiation_acceleration = radiation_calculator(du, u, p, p.radiation_mode)
    centrifugal_term = p.line.angular_momentum^2 / r^3
    isnan(centrifugal_term) ? centrifugal_term = 0 : nothing
    a_r = gravitational_acceleration[1] + radiation_acceleration[1] + centrifugal_term
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - v_r
    out[2] = z_dot - v_z
    out[3] = v_r_dot - a_r
    out[4] = v_z_dot - a_z
end

function initialize_solver(line::Streamline, params::Parameters; atol = 1e-7, rtol = 1e-3)
    termination_callback =
        DiscreteCallback(termination_condition, affect!, save_positions = (false, false))
    saving_callback = SavingCallback(save, SavedValues(Float64, Array{Float64,1}))
    callback_set = CallbackSet(termination_callback, saving_callback)
    a₀ = compute_initial_acceleration(
        r(line),
        z(line),
        v_r(line),
        v_z(line),
        params,
    )
    du₀ = [v_r_0(line), v_z_0(line), a₀[1], a₀[2]]
    u₀ = [r_0(line), z_0(line), v_r_0(line), v_z_0(line)]
    println(a₀)
    tspan = (0.0, 1e8)
    dae_problem = DAEProblem(
        residual!,
        du₀,
        u₀,
        tspan,
        params,
        differential_vars = [true, true, true, true],
    )
    integrator = init(dae_problem, IDA(init_all=true), callback = callback_set)
    integrator.opts.abstol = atol
    integrator.opts.reltol = rtol
    return integrator
end

#function initialize_solver(line::Streamline, bh::BlackHole, grid::Grid, mode::RE)
#    params = Parameters(
#        line.line_width_norm,
#        mass(bh),
#        r_in(grid),
#        number_density_0(line),
#        r_min,
#        r_max,
#        z_min,
#        z_max,
#    )
#end
#
