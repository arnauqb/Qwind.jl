using DrWatson
using Cubature
using PyPlot
@quickactivate "Qwind"
using Qwind

function compute_fr(integrator)
    data = integrator.p.data
    fr_r = data[:disc_radiation_field_r] .* exp.(-data[:tauuv]) .* (1 .+ data[:fm])
    fr_z = data[:disc_radiation_field_z] .* exp.(-data[:tauuv]) .* (1 .+ data[:fm])
    return hcat(fr_r, fr_z)
end

model = Model("scripts/config_re.yaml")
#iterations_dict = Dict()
#run!(model, iterations_dict)


function compute_vertical_flux(r)
    disk_flux_norel(model.bh, r) 
end

function estimate_from_force(r, z=1)
    acc_units = model.bh.M *G /model.bh.Rg/model.bh.Rg; #Acceleration scaling;
    res, err = integrate_radiation_force_integrand(
        model.rt,
        r,
        z,
        6,
        1400,
        atol = 0,
        rtol = 1e-3,
        norm = Cubature.INDIVIDUAL,
        maxevals = 0,
    )
    radiation_constant = compute_radiation_constant(model.rt.radiation) / model.rt.radiation.fuv;
    force = z * radiation_constant .* res[2];
    a_z = force * acc_units;
    F = a_z * C / (SIGMA_E) 
    return F
end

function find_closest_fraction(number)
    distance = Inf
    num = 0
    den = 0
    for i in 1:9
        for j in 1:9
            dist = abs(number - i/j) 
            if dist < distance
                distance = dist
                num = i
                den = j
            end
        end
    end
    return "$num/$den"
end

r_range = range(200, 1000, length=50);
flux = compute_vertical_flux.(r_range);
estimate = estimate_from_force.(r_range, 1);
fig, ax = plt.subplots()
ax.plot(r_range, flux, label="SS Flux", linewidth=2)
ax.plot(r_range, estimate, label="Qwind", linewidth=4, alpha=0.5)
#ratios = estimate ./ flux
#ax.plot(r_range, ratios, label="Qwind")
#ax.set_yscale("log")
ax.legend()
#ax.set_ylim(0.0, 2)

# lets check we recover total luminosity.

function S(r)
    return 3 * 0.5 / (8 * π * model.bh.efficiency) / r^2 # * ( 1- sqrt(6/r)) / r^2
end
using QuadGK
integral, err = quadgk(S, 6, 2000, rtol=1e-8)
println(integral * 2 * π * 0.85)


println(compute_bolometric_luminosity(model.bh))

