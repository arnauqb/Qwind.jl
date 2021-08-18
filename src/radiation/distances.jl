export distance_from_disk

distance_from_disk(rd, phid, zd, r, phi, z)= sqrt(r^2 + rd^2 + (z - zd)^2 - 2 * r * rd * cos(phi-phid))
