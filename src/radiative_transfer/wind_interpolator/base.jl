import ConcaveHull, Interpolations, Sundials
export WindInterpolator, get_density, update_interpolator


struct WindInterpolator{T} <: Interpolator{T}
    wind_hull::Union{ConcaveHull.Hull,Nothing}
    density_grid::InterpolationGrid{T}
    velocity_grid::Union{InterpolationGrid{T}, Nothing}
    vacuum_density::T
    n_timesteps::Int
end

function WindInterpolator(
    integrators;
    nr = "auto",
    nz = 50,
    vacuum_density = 1e2,
    n_timesteps = 1000,
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if integrators === nothing
        hull = nothing
        density_grid = construct_density_grid(nr, nz)
        velocity_grid = nothing 
    else
        r0 = [integ.p.r0 for integ in integrators]
        @info "Constructing wind hull..."
        hull = construct_wind_hull(integrators)
        flush()
        @info "Constructing interpolation grid..."
        flush()
        r, z, n = reduce_integrators(integrators, n_timesteps = 1000)
        density_grid = construct_density_grid(r, z, n, r0, hull, nr = nr, nz = nz)
        velocity_grid = nothing
    end
    return WindInterpolator(hull, density_grid, velocity_grid, vacuum_density, n_timesteps)
end

function reduce_line(r::Vector{<:Number}, z::Vector{<:Number}, n::Vector{<:Number})
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    for (rp, zp, np) in zip(r, z, n)
        if (zp > zs[end])
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
        end
    end
    return rs, zs, ns
end

function reduce_line(integrator::Sundials.IDAIntegrator; n_timesteps = 10000)
    dense_integrator = DenseIntegrator(integrator, n_timesteps=n_timesteps, log=true)
    return reduce_line(dense_integrator.r, dense_integrator.z, dense_integrator.n)
end

function reduce_integrators(integrators; n_timesteps = 10000)
    rs = Float64[]
    zs = Float64[]
    ns = Float64[]
    for integ in integrators
        r, z, n = reduce_line(integ, n_timesteps = n_timesteps)
        rs = vcat(rs, r)
        zs = vcat(zs, z)
        ns = vcat(ns, n)
    end
    return rs, zs, ns
end


function update_interpolator(wi::WindInterpolator, integrators)
    if maximum(wi.density_grid.grid) == wi.vacuum_density
        # first iteration, do not average
        return WindInterpolator(
            integrators,
            nr = wi.density_grid.nr,
            nz = wi.density_grid.nz,
            vacuum_density = wi.vacuum_density,
            n_timesteps = wi.n_timesteps,
        )
    end
    @info "Constructing wind hull..."
    flush()
    new_hull = construct_wind_hull(integrators)
    @info "Constructing interpolation grid..."
    density_grid = update_density_grid(wi.density_grid, new_hull, integrators)
    velocity_grid = nothing
    #old_density_grid = wi.density_grid
    flush()
    return WindInterpolator(new_hull, density_grid, velocity_grid, wi.vacuum_density, wi.n_timesteps)
end

function get_density(wi::WindInterpolator, r, z)
    if is_point_in_wind(wi, [r, z])
        return wi.grid.interpolator(r, z)
    else
        return wi.vacuum_density
    end
end
get_density(wi::WindInterpolator, point) = get_density(wi, point[1], point[2])
