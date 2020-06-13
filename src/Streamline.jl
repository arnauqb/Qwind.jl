export Streamline,
    Streamlines,
    r,
    z,
    v_r,
    v_r_0,
    v_z,
    v_z_0,
    max_r,
    max_z,
    r_0,
    z_0,
    r_in,
    create_line,
    initial_radii,
    escaped

struct Streamline
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
end

create_line() = Streamline(
    Float64[], #r
    Float64[], #z
    Float64[], #v_r
    Float64[], #v_z
    Float64[], #density
    Float64[], #tau_x
    Float64[], #tau_uv
    Float64[], #force_multiplier
    Float64[], #dv_dr
    Float64[], #xi
    0.0, #angular_momentum
    0.0, #line_width_norm
)
function create_line(
    r_0,
    z_0,
    v_z_0,
    number_density_0,
    line_width;
    tau_x_0 = 0,
    tau_uv_0 = 0,
    force_multiplier_0 = 0,
    ionization_parameter_0 = 0,
    dv_dr_0 = 0,
    v_r_0 = 0,
)
    angular_momentum = sqrt(1.0 / r_0) * r_0
    line_width_norm = line_width / r_0
    return Streamline(
        Float64[r_0],
        Float64[z_0],
        Float64[v_r_0],
        Float64[v_z_0],
        Float64[number_density_0],
        Float64[tau_x_0],
        Float64[tau_uv_0],
        Float64[force_multiplier_0],
        Float64[dv_dr_0],
        Float64[ionization_parameter_0],
        angular_momentum,
        line_width_norm,
    )
end
r(line::Streamline) = line.r[end]
z(line::Streamline) = line.z[end]
v_r(line::Streamline) = line.v_r[end]
v_r_0(line::Streamline) = line.v_r[1]
v_z(line::Streamline) = line.v_z[end]
v_z_0(line::Streamline) = line.v_z[1]
max_r(line::Streamline) = maximum(line.r)
max_z(line::Streamline) = maximum(line.z)
r_0(line::Streamline) = line.r[1]
z_0(line::Streamline) = line.z[1]

function escaped(line::Streamline, bh_mass)
    d = sqrt.(line.r .^ 2 + line.z .^ 2)
    v_T = sqrt.(line.v_r .^ 2 + line.v_z .^ 2)
    v_esc = sqrt.((2 * G * bh_mass) ./ d)
    return any(v_T > v_esc)
end

struct Streamlines
    lines::Vector{Streamline}
end
initial_radii(lines::Streamlines) = [line.r[1] for line in lines.lines]
r_in(lines::Streamlines) = lines.lines[1].r[1]
max_r(lines::Streamlines) = maximum([max_r(line) for line in lines.lines])
max_z(lines::Streamlines) = maximum([max_z(line) for line in lines.lines])
