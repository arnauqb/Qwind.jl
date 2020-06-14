export BlackHole,
    eddington_luminosity,
    bolometric_luminosity,
    isco,
    efficiency,
    mass_accretion_rate,
    Rg,
    Rs,
    solar_mass,
    gravitational_acceleration

struct BlackHole
    M::Float64
    mdot::Float64
    spin::Float64
end

function eddington_luminosity(bh::BlackHole)
    return 4π * G * bh.M * M_P * C / SIGMA_T
end

function bolometric_luminosity(bh::BlackHole)
    return bh.mdot * eddington_luminosity(bh)
end

Rg(bh::BlackHole) = G * bh.M / C^2
Rs(bh::BlackHole) = 2 * Rg(bh)

function isco(bh::BlackHole)
    z1 =
        1 +
        (1 - bh.spin^2)^(1 / 3) *
        ((1 + abs(bh.spin))^(1 / 3) + (1 - abs(bh.spin))^(1 / 3))
    z2 = sqrt(3 * bh.spin^2 + z1^2)
    rms = 3 + z2 - sign(bh.spin) * sqrt((3 - z1) * (3 + z1 + 2 * z2))
    return rms * Rg(bh)
end

solar_mass(bh::BlackHole) = bh.M / M_SUN

function efficiency(bh::BlackHole)
    return 1 - sqrt(1 - 2 / (3 * isco(bh)))
end

function mass_accretion_rate(bh::BlackHole)
    return bolometric_luminosity(bh) / (efficiency(bh) * C^2)
end

function gravitational_acceleration(r, z, bh::BlackHole)
    d = sqrt(r^2 + z^2)
    return G * bh.M / d^3 * [r, z]
end

