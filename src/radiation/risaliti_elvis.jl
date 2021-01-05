"""
This radiation module implements the radiation quantities as computed
in the first version of Qwind (Risaliti & Elvis 2010, Quera-Bofarull et al. 2020)
"""

using Roots
export RERadiation,
    radiation_force_integrand!, compute_disc_radiation_field, compute_radiation_constant

struct RERadiation <: Radiation
    mdot::Float64
    efficiency::Float64
    isco::Float64
    spin::Float64
    fuv::Float64
    xray_luminosity::Float64
    Rg::Float64
    function RERadiation(bh::BlackHole, fuv, fx)
        bolometric_luminosity = compute_bolometric_luminosity(bh)
        xray_luminosity = fx * bolometric_luminosity
        new(
            bh.mdot,
            compute_efficiency(bh),
            compute_isco(bh),
            bh.spin,
            fuv,
            xray_luminosity,
            bh.Rg,
        )
    end
end

function RERadiation(bh::BlackHole, config::Dict)
    radiation_config = config["radiation"]
    return RERadiation(bh, radiation_config["f_uv"], radiation_config["f_x"])
end

compute_radiation_constant(radiation::RERadiation) =
    6 * radiation.mdot * radiation.fuv / (8 * π * radiation.efficiency)

