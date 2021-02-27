using Roots
ic = model1.ic
rin = 50
rfi = 1200
Rg = model1.bh.Rg
xray_luminosity = model1.rad.xray_luminosity

lines_range = [];
lines_widths = [];
r_range = 10 .^ range(log10(rin), log10(rfi), length = 100);
z_range = [0.0, 1.0];
density_grid = zeros((length(r_range), length(z_range)));
density_grid[:, 1] .= getn0.(Ref(ic), r_range);
interp_grid = InterpolationGrid(r_range, z_range, density_grid);
tau_x = 0;
tau_uv = 0;
fx(delta_r, rc, delta_tau, tau_x) =
    delta_tau -
    (compute_xray_tau(interp_grid, 0, 0, rc + delta_r, 0, xray_luminosity, Rg) - tau_x);
fuv(delta_r, rc, delta_tau, tau_uv) =
    delta_tau -
    (compute_uv_tau(interp_grid, 0, 0, rc + delta_r, 0, Rg) - tau_uv);
rc = rin;
while rc < rfi
    if tau_x < 50
        if tau_x < 10
            delta_tau = 0.05;
        else
            delta_tau = 0.5;
        end
        tau_x = compute_xray_tau(interp_grid, 0, 0, rc, 0, xray_luminosity, Rg);
        delta_r = find_zero(delta_r ->fx(delta_r, rc, delta_tau, tau_x), 0.1);
    elseif tau_uv < 50
        if tau_uv < 10
            delta_tau = 0.1;
        else
            delta_tau = 1;
        end
        tau_uv = compute_uv_tau(interp_grid, 0, 0, rc, 0, Rg);
        delta_r = find_zero(delta_r ->fuv(delta_r, rc, delta_tau, tau_uv) , 1);
    else
        break
    end
    push!(lines_range, rc + delta_r / 2);
    push!(lines_widths, delta_r);
    rc += delta_r;
end;
# distribute remaining ones logarithmically;
additional_range = 10 .^ range(log10(rc), log10(rfi), length=50);
additional_widths = diff(additional_range);
pushfirst!(additional_widths, additional_range[1] - lines_range[end]);
lines_range = vcat(lines_range, additional_range);
lines_widths = vcat(lines_widths, additional_widths);




compute_xray_tau(interp_grid, 0, 0, 51, 0, xray_luminosity, Rg)

r_range = range(0, 1000, length=5000);
taux = compute_xray_tau.(Ref(interp_grid), 0, 0, r_range, 0, xray_luminosity, Rg);
tauuv = compute_uv_tau.(Ref(interp_grid), 0, 0, r_range, 0, Rg);

loglog(r_range, taux)
loglog(r_range, tauuv)

rr, ww = Qwind.compute_lines_range(model1);
length(rr)
