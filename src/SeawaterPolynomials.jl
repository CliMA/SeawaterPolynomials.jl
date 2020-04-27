module SeawaterPolynomials

"""
    AbstractSeawaterPolynomial

Abstract type for Boussinesq equations of state expressed as polynomial functions
of conservative temperature, absolute salinity, and geopotential depth.
"""
abstract type AbstractSeawaterPolynomial end

"""
    ρ(Θ, Sᴬ, Z, equation_of_state)

Returns the total density of a seawater parcel with conservative tempertuare `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, using the Boussinesq `equation_of_state`.

The geopotential height Z(x, y) is defined such that `Z(x, y) = 0` at sea level and 
*decreases* downwards to negative values, towards the bottom of the ocean.

This function aliases `total_density`.
"""
@inline ρ(Θ, Sᴬ, Z, eos) = eos.reference_density + ρ′(Θ, Sᴬ, Z, eos)

"""
    ρ′(Θ, Sᴬ, Z, equation_of_state)

Returns the total density of a seawater parcel with conservative tempertuare `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, using the Boussinesq `equation_of_state`.

The geopotential height Z(x, y) is defined such that `Z(x, y) = 0` at sea level and 
*decreases* downwards to negative values, towards the bottom of the ocean.

The function aliases `density_anomaly`.
"""
function ρ′ end

"""
    thermal_expansion(Θ, Sᴬ, Z, equation_of_state)

Returns the thermal expansion coefficient,

``α = - ∂ρ / ∂Θ``

for a seawater parcel with conservative tempertuare `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, and using the 
Boussinesq `equation_of_state`.
The thermal expansion coefficient measures how much seawater density changes
when conservative temperature is changed.
'Thermal expansion' is so named because, due to sign convention, positive values
reflect decreasing seawater density with increasing conservative temperature, and thus
an 'expansion' of oceanic fluid parcels.
In many, but not all conditions in Earth's ocean (at temperatures greater than 4ᵒC in freshwater),
the thermal expansion coefficient is positive.

The geopotential height Z(x, y) is defined such that `Z(x, y) = 0` at sea level and 
*decreases* downwards to negative values, towards the bottom of the ocean.
"""
function thermal_expansion end

"""
    haline_contraction(Θ, Sᴬ, Z, equation_of_state)

Returns the haline contraction coefficient,

``β = ∂ρ / ∂Sᴬ``

for a seawater parcel with absolute salinity `Sᴬ` and at fixed 
conservative tempertuare `Θ` and geopotential height `Z`, using the 
Boussinesq `equation_of_state`.
The haline contraction coefficient measures how much seawater density changes
when absolute salinity is changed.
'Haline contraction' is so named because, due to sign convention, positive values
reflect increasing seawater density with increasing absolute salinity, and thus
a slight `contraction' of oceanic fluid parcels.

The geopotential height Z(x, y) is defined such that `Z(x, y) = 0` at sea level and 
*decreases* downwards to negative values, towards the bottom of the ocean.
"""
function haline_contraction end

const total_density = ρ
const density_anomaly = ρ′

struct BoussinesqEquationOfState{P, FT}
    seawater_polynomial :: P
      reference_density :: FT
end

@inline reference_density(eos) = eos.reference_density

include("SecondOrderSeawaterPolynomials.jl")
include("TEOS10.jl")

using .SecondOrderSeawaterPolynomials
using .TEOS10

end # module
