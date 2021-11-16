export compute_distance_cylindrical


function compute_distance_cylindrical(rd, phid, zd, r, phi, z)
    try
        return sqrt(r^2 + rd^2 + (z - zd)^2 - 2 * r * rd * cos(phi - phid))
    catch DomainError # numerical overflow to negative
        return 0.0
    end
end
