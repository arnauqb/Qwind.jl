export InitialConditions,
    UniformIC, CAKIC, getz0, getrin, getrfi, getn0, getv0, getnlines, getl0
abstract type InitialConditions end

struct UniformIC <: InitialConditions
    rin::Float64
    rfi::Float64
    nlines::Int64
    z0::Float64
    n0::Float64
    v0::Float64
    logspaced::Bool
end

struct CAKIC <: InitialConditions
    rin::Float64
    rfi::Float64
    nlines::Int64
    z0::Float64
    K::Float64
    alpha::Float64
    mu::Float64
    logspaced::Bool
end

getrin(ic::UniformIC) = ic.rin
getrfi(ic::UniformIC) = ic.rfi
getnlines(ic::UniformIC) = ic.nlines
getz0(ic::UniformIC, r) = ic.z0
getn0(ic::UniformIC, r) = ic.n0
getv0(ic::UniformIC, r) = ic.v0
getl0(ic::InitialConditions, r) = sqrt(r)
