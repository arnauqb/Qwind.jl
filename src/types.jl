# general pieces of a model

export DensityInterpolator 

abstract type Radiation end
abstract type RadiativeTransfer end
abstract type DensityInterpolator end
abstract type InitialConditions end
abstract type Grid end

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
abstract type Flag end
abstract type FMInterpolationType <: Flag end
struct FMNoInterp <: FMInterpolationType end
struct FMInterp <: FMInterpolationType end

# Radiative transfer
abstract type CellIterator end
#abstract type LogGridIterator <: CellIterator end
