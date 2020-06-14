export Streamline,
    Streamlines,
    r,
    z,
    v_r,
    v_r_0,
    v_z,
    v_z_0,
    v_0,
    line_width,
    number_density,
    number_density_0,
    max_r,
    max_z,
    r_0,
    z_0,
    r_in,
    create_line,
    initial_radii,
    escaped,
    save_position!

struct Streamline
    t::Vector{Float64}
    r::Vector{Float64}
    z::Vector{Float64}
    v_r::Vector{Float64}
    v_z::Vector{Float64}
    number_density::Vector{Float64}
    tau_x::Vector{Float64}
    tau_uv::Vector{Float64}
    force_multiplier::Vector{Float64}
    dv_dr::Vector{Float64}
    ionization_parameter::Vector{Float64}
    angular_momentum::Float64
    line_width_norm::Float64
    function Streamline(r, z, v_r, v_z, number_density, line_width, bh_mass)
        tau_x = [0.0]
        tau_uv = [0.0]
        force_multiplier = [0.0]
        dv_dr = [0.0]
        ionization_parameter = [0.0]
        angular_momentum = sqrt(G * bh_mass * r)
        line_width_norm = line_width / r
        new(
            [0.0],
            [r],
            [z],
            [v_r],
            [v_z],
            [number_density],
            tau_x,
            tau_uv,
            force_multiplier,
            dv_dr,
            ionization_parameter,
            angular_momentum,
            line_width_norm,
        )
    end
end

r(line::Streamline) = line.r[end]
z(line::Streamline) = line.z[end]
line_width(line::Streamline, r) = line.line_width_norm * r
line_width(line::Streamline) = line.line_width_norm * r(line)
v_r(line::Streamline) = line.v_r[end]
v_r_0(line::Streamline) = line.v_r[1]
v_z(line::Streamline) = line.v_z[end]
number_density_0(line::Streamline) = line.number_density[1]
number_density(line::Streamline) = line.number_density[end]
v_z_0(line::Streamline) = line.v_z[1]
v_0(line::Streamline) = sqrt(line.v_z[1]^2 + line.v_r[1]^2)
max_r(line::Streamline) = maximum(line.r)
max_z(line::Streamline) = maximum(line.z)
r_0(line::Streamline) = line.r[1]
z_0(line::Streamline) = line.z[1]

function save_position!(line::Streamline, u, t)
    r, z, v_r, v_z = u
    push!(line.t, t)
    push!(line.r, r)
    push!(line.z, z)
    push!(line.v_r, v_r)
    push!(line.v_z, v_z)
end

function escaped(line::Streamline, bh_mass)
    d = sqrt.(line.r .^ 2 + line.z .^ 2)
    v_T = sqrt.(line.v_r .^ 2 + line.v_z .^ 2)
    v_esc = sqrt.((2 * G * bh_mass) ./ d)
    return any(v_T > v_esc)
end

struct Streamlines
    lines::Vector{Streamline}
end
Base.length(iter::Streamlines) = length(iter.lines)
initial_radii(lines::Streamlines) = [line.r[1] for line in lines.lines]
r_in(lines::Streamlines) = lines.lines[1].r[1]
max_r(lines::Streamlines) = maximum([max_r(line) for line in lines.lines])
max_z(lines::Streamlines) = maximum([max_z(line) for line in lines.lines])


function Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    bh_mass;
    z_0 = 0.0,
    log_spaced = true,
)
    if log_spaced
        line_delimiters = 10 .^ range(log10(r_in), log10(r_fi), length = n_lines + 1)
        lines_range = []
        for i = 1:n_lines
            r0 = line_delimiters[i] + (line_delimiters[i+1] - line_delimiters[i]) / 2.0
            push!(lines_range, r0)
        end
        lines_widths = diff(line_delimiters)
    else
        dr = (r_fi - r_in) / n_lines
        lines_range = [r_in + (i + 0.5) * dr for i = 0:n_lines-1]
        lines_widths = diff([lines_range; r_fi + dr / 2])
    end
    lines = []
    for (r0, line_width) in zip(lines_range, lines_widths)
        v0 = velocity_generator(r0)
        density = density_generator(r0)
        line = Streamline(r0, 0.0, 0.0, v0, density, line_width, bh_mass)
        push!(lines, line)
    end
    return Streamlines(lines)
end
