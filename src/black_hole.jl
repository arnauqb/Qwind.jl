export BlackHole,
    eddington_luminosity,
    bolometric_luminosity,
    isco,
    efficiency,
    compute_mass_accretion_rate,
    Rg,
    Rs,
    solar_mass,
    mass,
    compute_gravitational_acceleration,
    compute_escape_velocity


struct BlackHole
    M::Float64
    mdot::Float64
    spin::Float64
    Rg::Float64
    isco::Float64
    efficiency::Float64
    mass_accretion_rate::Float64
    function BlackHole(M, mdot, spin)
        Rg = G * M / C^2
        isco = compute_isco(spin)
        efficiency = compute_blackhole_efficiency(isco, Rg)
        mass_accretion_rate = compute_mass_accretion_rate(M, mdot, spin, efficiency)
        new(M, mdot, spin, Rg, isco, efficiency, mass_accretion_rate)
    end
end

function eddington_luminosity(bh_mass)
    return 4π * G * bh_mass * M_P * C / SIGMA_T
end

function eddington_luminosity(bh::BlackHole)
    return eddington_luminosity(bh.M)
end

function bolometric_luminosity(bh_mass, eddinton_rate)
    return eddinton_rate * eddington_luminosity(bh_mass)
end
function bolometric_luminosity(bh::BlackHole)
    return bolometric_luminosity(bh.M, bh.mdot)
end

Rg(bh::BlackHole) = bh.Rg
Rg(M::Float64) = G * M / C^2
Rs(bh::BlackHole) = 2 * bh.Rg
compute_escape_velocity(distance) = sqrt(2 / distance)

function compute_isco(spin)
    z1 =
        1 +
        (1 - spin^2)^(1 / 3) * ((1 + abs(spin))^(1 / 3) + (1 - abs(spin))^(1 / 3))
    z2 = sqrt(3 * spin^2 + z1^2)
    rms = 3 + z2 - sign(spin) * sqrt((3 - z1) * (3 + z1 + 2 * z2))
    return rms
end

mass(bh::BlackHole) = bh.M
solar_mass(bh::BlackHole) = bh.M / M_SUN

function compute_blackhole_efficiency(isco, Rg)
    return 1 - sqrt(1 - 2 / (3 * isco))
end

function compute_mass_accretion_rate(M, mdot, spin, efficiency)
    return bolometric_luminosity(M, mdot) / (efficiency * C^2)
end

function compute_mass_accretion_rate(bh::BlackHole)
    return bolometric_luminosity(bh) / (bh.efficiency * C^2)
end

function compute_gravitational_acceleration(r, z, M)
    d = sqrt(r^2 + z^2)
    force = - 1.0 / d^3 * [r, z]
    return force
end

compute_gravitational_acceleration(r, z, bh::BlackHole) =
    compute_gravitational_acceleration(r, z, bh.M)
