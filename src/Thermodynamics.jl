export compute_thermal_velocity

function compute_thermal_velocity(temperature, mu)
    return sqrt(K_B * temperature / (mu * M_P)) / C
end

compute_thermal_velocity(temperature) = compute_thermal_velocity(temperature, 1.0)
