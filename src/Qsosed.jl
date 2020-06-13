export nt_rel_factors

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
function nt_rel_factors(radius, spin, isco)
    if radius <= isco
        return 0.
    end
    yms = sqrt(isco)
    y1 = 2 * cos((acos(spin) - pi) / 3)
    y2 = 2 * cos((acos(spin) + pi) / 3)
    y3 = -2 * cos(acos(spin) / 3)
    y = sqrt(radius);
    C = 1 - 3 / radius + 2 * spin / radius^1.5
    B = 3 * (y1 - spin)^2 * log((y - y1) / (yms - y1)) / (y * y1 * (y1 - y2) * (y1 - y3))
    B += 3 * (y2 - spin)^2 * log((y - y2) / (yms - y2)) / (y * y2 * (y2 - y1) * (y2 - y3))
    B += 3 * (y3 - spin)^2 * log((y - y3) / (yms - y3)) / (y * y3 * (y3 - y1) * (y3 - y2))
    A = 1 - yms / y - 3 * spin * log(y / yms) / (2 * y)
    factor = (A - B) / C
    return factor
end
