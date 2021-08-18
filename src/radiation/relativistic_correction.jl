export relativistic_correction, beta, gamma

compute_beta(vr, vz) = max(min(1.0, sqrt(vr^2 + vz^2)), 1e-4)
compute_gamma(beta) = 1.0 / sqrt(1 - beta^2)

relativistic_correction(::NoRelativistic; r, z, vr, vz, beta, gamma, r_projection, delta) =
    1.0
function relativistic_correction(
    ::Relativistic;
    r,
    z,
    z0,
    vr,
    vz,
    beta,
    gamma,
    r_projection,
    delta,
)
    cosθ = (r_projection * vr + (z - z0) * vz) / (delta * beta)
    return 1.0 / (gamma * (1 + beta * cosθ))^4
end

# for testing
function relativistic_correction(rd, phid, zd, r, z, vr, vz)
    beta = compute_beta(vr, vz)
    gamma = compute_gamma(beta) 
    r_projection = (r - rd * cos(phid))
    delta = distance_from_disk(rd, phid, zd, r, 0.0, z)
    return relativistic_correction(
        Relativistic(),
        r = r,
        z = z,
        vr = vr,
        vz = vz,
        beta = beta,
        gamma = gamma,
        r_projection = r_projection,
        delta = delta,
    )
end
