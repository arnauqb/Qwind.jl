# general pieces of a model
abstract type Radiation end
abstract type RadiativeTransfer end
abstract type InitialConditions end

struct Model
    grid::Grid
    bh::BlackHole
    rad::Radiation
    rt::RadiativeTransfer
    ic::InitialConditions
end

# Disk integration modes
abstract type IntegrationType end
struct IntegrationFromStreamline <: IntegrationType end
struct IntegrationFromCenter <: IntegrationType end

# Force multiplier interpolation
abstract type FMInterpolationType <: Flag end
struct FMNoInterp <: InterpolationType end
struct FMInterp <: InterpolationType end
