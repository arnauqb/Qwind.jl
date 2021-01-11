export BlackBody, spectral_radiance

struct BlackBody
    T::Float64
end

function spectral_radiance(bb::BlackBody, energy)
    a = 2 / (H_PLANCK^2 * C^2)
    # energy to ergs
    energy_erg = energy / ERG_TO_KEV
    b = energy_erg^3 / (exp(energy_erg / (K_B * bb.T)) - 1)
    ret = a * b
    # convert back to keV
    ret = ret * ERG_TO_KEV / HZ_TO_KEV
end
