using Cubature
export integrate_radiation_force_integrand, radiation_force_integrand!

function radiation_force_integrand!(
    radiation::SimpleRadiation,
    v,
    r_d,
    phi_d,
    r,
    z,
)
    #println("-----")
    #println("r_d : $r_d, phi_d : $phi_d, r : $r, z $z")
    #println("\n")
    nt = nt_rel_factors(radiation, r_d)
    r_projection = (r - r_d * cos(phi_d))
    delta_sq = (r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    common_projection = 1.0 / (r_d^2 * delta_sq^2)
    v[:] = nt * common_projection * [r_projection, z]
    #println("nt $nt")
    #println("r_proj $r_projection")
    #println("delta sq : $delta_sq")
    #println("common $common_projection")
    #println("-----")
    #println("v : $v")
end

"""
The signature of the input integrand has to be (radiation::Radiation, v::Array{Float64, 2}, r_d::Float64, phi_d::Float64, r::Float64, z::Float64)
"""


function integrate_radiation_force_integrand(
    radiation::Radiation,
    r,
    z;
    r_lims,
    phi_lims,
    atol = 0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = convert(Int, 5e4),
)
    f(x, v) = radiation_force_integrand!(radiation, v, x[1], x[2], r, z)
    return hcubature(
        2,
        f,
        (r_lims[1], phi_lims[1]),
        (r_lims[2], phi_lims[2]),
        abstol = atol,
        reltol = rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
end
