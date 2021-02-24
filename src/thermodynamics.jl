export compute_thermal_velocity, compute_density

"""
Computes gas thermal velocity.

# Parameters
- temperature: gas temperature in K
- mu : gas mean molecular weight. Default = 0.5 (ionized hydrogen gas)
"""
function compute_thermal_velocity(temperature, mu)
    return sqrt(K_B * temperature / (mu * M_P)) / C
end
compute_thermal_velocity(temperature) = compute_thermal_velocity(temperature, 0.5)

"""
Updates the density of the streamline giving its current position and velocity,
using mass conservation.

# Parameters
- r : radial coordinate [Rg]
- z : height coordinate [Rg]
- v_t :
"""
function compute_density(r, z, v_r, v_z, r_0, v_0, n_0)
    d = sqrt(r^2 + z^2)
    v_t = sqrt(v_r^2 + v_z^2)
    radial = (r_0 / d)^2
    v_ratio = v_0 / v_t
    return n_0 * radial * v_ratio
end
