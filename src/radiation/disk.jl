export disk_nt_rel_factors, gravity_radius, disk_height, disk_temperature

"""
Novikov-Thorne relativistic factors for the AD spectrum.
If the passed radius is smaller than ISCO, it returns 0.

Parameters
----------
radius
    disc radius measured in R_g
spin
    normalized black hole spin (-1, 1)
isco
    innermost stable circular orbit measured in Rg
"""
function disk_nt_rel_factors(radius, spin, isco)
    if radius <= isco
        return 0.0
    end
    yms = sqrt(isco)
    y1 = 2 * cos((acos(spin) - pi) / 3)
    y2 = 2 * cos((acos(spin) + pi) / 3)
    y3 = -2 * cos(acos(spin) / 3)
    y = sqrt(radius)
    C = 1 - 3 / radius + 2 * spin / radius^1.5
    B = 3 * (y1 - spin)^2 * log((y - y1) / (yms - y1)) / (y * y1 * (y1 - y2) * (y1 - y3))
    B += 3 * (y2 - spin)^2 * log((y - y2) / (yms - y2)) / (y * y2 * (y2 - y1) * (y2 - y3))
    B += 3 * (y3 - spin)^2 * log((y - y3) / (yms - y3)) / (y * y3 * (y3 - y1) * (y3 - y2))
    A = 1 - yms / y - 3 * spin * log(y / yms) / (2 * y)
    factor = (A - B) / C
    return factor
end
disk_nt_rel_factors(bh::BlackHole, r) = disk_nt_rel_factors(r, bh.spin, bh.isco)

"""
Computes the disk flux of the accretion disc of black hole ``bh``
at a given radius ``r`` (in Rg). See equation 5 of https://arxiv.org/pdf/2001.04720.pdf
"""
function disk_flux(bh::BlackHole, r)
    Mdot = compute_mass_accretion_rate(bh)
    NT_factors = disk_nt_rel_factors(bh, r)
    return 3 * G * bh.M * Mdot / (8 * π * (r * bh.Rg)^3) * NT_factors
end

"""
Computes the temperature of the disk at a radius r (in units of Rg).
"""
function disk_temperature(bh::BlackHole, r)
    flux = disk_flux(bh, r)
    return (flux / SIGMA_SB)^(1 / 4)
end

"""
Self-gravity radius as described by Laor & Netzer (1989).
"""
function gravity_radius(bh::BlackHole)
    mass = (bh.M / (1e9 * M_SUN))
    alpha = 0.1 # assumption
    r_sg = 2150 * mass^(-2 / 9) * bh.mdot^(4 / 9) * alpha^(2 / 9)# / bh.Rg
    return r_sg
end

"""
Disk height
"""
function disk_height(bh::BlackHole, r; mu_e = 1.17)
    return SIGMA_E * mu_e * SIGMA_SB * disk_temperature(bh, r)^4 * (r * bh.Rg)^3 /
           (G * bh.M * C) / bh.Rg
end

function compute_proga_density(bh, r; mu_nucleon = 0.61)
    temp = disk_temperature(bh, r)
    cs = compute_thermal_velocity(temp) * C
    z0 = 3 * bh.mdot * bh.Rg
    r = r * bh.Rg
    ttheta = tan(r/z0)
    rho0 = 1e-9 # g / cm^3
    rho = rho0 * exp(-(G * bh.M) / (2 * cs^2 * r * ttheta^2))
    return rho / (M_P * mu_nucleon)
end
