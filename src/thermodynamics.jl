export compute_thermal_velocity, compute_density

"""
Computes gas thermal velocity.

# Parameters
- temperature: gas temperature in K
- mu : gas mean molecular weight. Default = 0.61 (ionized solar abundance gas)
"""
function compute_thermal_velocity(temperature, mu_nucleon)
    return sqrt(K_B * temperature / (mu_nucleon * M_P)) / C
end
compute_thermal_velocity(temperature) = compute_thermal_velocity(temperature, 0.61)

"""
Updates the density of the streamline giving its current position and velocity,
using mass conservation.

# Parameters
- r : radial coordinate [Rg]
- z : height coordinate [Rg]
- v_t :
"""
function compute_density(r, z, v_r, v_z, r_0, z_0, v_0, n_0)
    d = sqrt(r^2 + z^2)
    d_0 = sqrt(r_0^2 + z_0^2)
    v_t = sqrt(v_r^2 + v_z^2)
    radial = (d_0 / d)^2
    if v_0 == v_t
        v_ratio = 1.0
    else
        v_ratio = v_0 / v_t
    end
    return n_0 * v_ratio * radial 
end
