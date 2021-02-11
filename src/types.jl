# general pieces of a model

export DensityInterpolator, GridInterpolator, FMNoInterp, FMInterp

abstract type Radiation end
abstract type RadiativeTransfer end
abstract type DensityInterpolator end
abstract type GridInterpolator <: DensityInterpolator end
abstract type NNInterpolator <: DensityInterpolator end
abstract type InitialConditions end
abstract type Grid end

# Disk integration modes
abstract type IntegrationType end
struct IntegrationFromStreamline <: IntegrationType end
struct IntegrationFromCenter <: IntegrationType end

# Radiation
# Force multiplier interpolation
abstract type Flag end
abstract type FMInterpolationType <: Flag end
struct FMNoInterp <: FMInterpolationType end
struct FMInterp <: FMInterpolationType end
# Relativistic correction to luminosity
abstract type FluxCorrection <: Flag end
struct Relativistic <: FluxCorrection end
struct NoRelativistic <: FluxCorrection end

# Radiative transfer
abstract type CellIterator end
#abstract type LogGridIterator <: CellIterator end
