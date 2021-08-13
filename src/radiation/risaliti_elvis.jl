"""
This radiation module implements the radiation quantities as computed
in the first version of Qwind (Risaliti & Elvis 2010, Quera-Bofarull et al. 2020)
"""

using Roots
export RERadiation,
    radiation_force_integrand!,
    compute_disc_radiation_field,
    compute_radiation_constant,
    get_fuv_mdot

struct RERadiation{T} <: Radiation{T}
    mdot::T
    efficiency::T
    isco::T
    spin::T
    fuv::T
    z_xray::T
    disk_height::T
    xray_opacity::XRayOpacity
    xray_luminosity::T
    Rg::T
    function RERadiation(
        bh::BlackHole;
        fuv::T,
        fx::T,
        z_xray::T,
        disk_height::T,
        xray_opacity::XRayOpacity,
    ) where {T<:AbstractFloat}
        bolometric_luminosity = compute_bolometric_luminosity(bh)
        xray_luminosity = fx * bolometric_luminosity
        new{typeof(fx)}(
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
    radiation_config = config[:radiation]
    if radiation_config[:xray_opacity] == "thomson"
        xray_opacity = Thomson()
    else
        xray_opacity = Boost()
    end
    return RERadiation(
        bh,
        fuv = radiation_config[:f_uv],
        fx = radiation_config[:f_x],
        z_xray = radiation_config[:z_xray],
        disk_height = radiation_config[:disk_height],
        xray_opacity = xray_opacity,
    )
end

compute_radiation_constant(radiation::RERadiation) =
    3 * radiation.mdot * radiation.fuv / (π * radiation.efficiency)

get_fuv_mdot(radiation::RERadiation, r) = radiation.fuv, radiation.mdot

