using Cubature
export integrate_2d, make_radiation_force_integrand

function integrate_2d(
    func;
    x_lims,
    y_lims,
    norm = Cubature.INDIVIDUAL,
    rtol = 1e-3,
    atol = eps(Float64),
    maxevals = typemax(Int),
)
    f(x, v) = func(v, x[1], x[2])
    xmin = (x_lims[1], y_lims[1])
    xmax = (x_lims[2], y_lims[2])
    return hcubature(
        2,
        f,
        xmin,
        xmax;
        error_norm = norm,
        reltol = rtol,
        abstol = atol,
        maxevals = maxevals,
    )
end

"""
Composes different functions to create a radiation force integrand.
An integrand is a function that modifies a vector v (2d) in place,
by passing the arguments (r_d, phi_d, r, z) to a collection of functions.
All functions should expect the arguments (r_d, phi_d, r, z)

# Parameters
- r_d_func:  function that it is only saved in v[1]
- phi_d_func: function that it is only saved in v[2]
- geometric_func : general geometric terms
- uv_func : uv radiation profile (constant/qsosed profile)
- nt_func : Relativistic corrections to the disc profile (NT/Newton/SS)
- mdot_func : mass accretion rate profile.

# Example
See radiation force integrand test in ``integration_tests.jl``.
"""
function make_radiation_force_integrand(
    r_d_func,
    phi_d_func,
    geometric_func,
    nt_func,
    mdot_func,
    uv_func,
)
    function integrand(v, r_d, phi_d, r, z)
        common =
            geometric_func(r_d, phi_d, r, z) *
            nt_func(r_d, phi_d, r, z) *
            mdot_func(r_d, phi_d, r, z) *
            uv_func(r_d, phi_d, r, z)
        v[1] = common * r_d_func(r_d, phi_d, r, z)
        v[2] = common * phi_d_func(r_d, phi_d, r, z)
    end
    return integrand
end

function make_radiation_force_integrand(;
    mdot_func=mdot_func,
    uv_func=uv_func,
    spin::Float64 = 0.0,
    isco::Float64 = 6.0,
)
    r_d_func(r_d, phi_d, r, z) = (r - r_d * cos(phi_d))
    phi_d_func(r_d, phi_d, r, z) = 1.0
    geometric_func(r_d, phi_d, r, z) =
        1 / r_d^2 * (r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    nt_func(r_d, phi_d, r, z) = nt_rel_factors(r_d, spin, isco)
    return make_radiation_force_integrand(
        r_d_func,
        phi_d_func,
        geometric_func,
        nt_func,
        mdot_func,
        uv_func,
    )
end

function make_qwind1_integrand(mass_accretion_rate, spin, isco)
    mdot_func(r_d, phi_d, r, z) = mass_accretion_rate
    uv_func(r_d, phi_d, r, z) = 1.0
    return make_radiation_force_integrand(mdot_func, uv_func)
end
