export compute_distance_cylindrical

compute_distance_cylindrical(rd, phid, zd, r, phi, z) = sqrt(r^2 + rd^2 + (z - zd)^2 - 2 * r * rd * cos(phi - phid))
