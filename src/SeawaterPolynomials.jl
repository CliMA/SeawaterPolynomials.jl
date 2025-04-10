module SeawaterPolynomials

import Base

"""
    thermal_sensitivity(Θ, Sᴬ, Z, equation_of_state)

Return the "Boussinesq thermal sensitivity coefficient" for a seawater parcel with
conservative temperature `Θ`, at fixed absolute salinity `Sᴬ`, and geopotential height `Z`
using the Boussinesq `equation_of_state`. The thermal sensitivity coefficient is

```math
a(Θ, Sᴬ, Z) = - \\left.\\frac{∂ρ}{∂Θ}\\right|_{Sᴬ, Z} ,
```

and measures how much seawater density changes when conservative temperature is changed.
It has units of [kg/m³/(g/kg)], and differs from `thermal_expansion` (``α``) by a factor of the reference density ``ρᵣ``.
In many, but not all conditions in Earth's ocean (at temperatures greater than
4ᵒC in freshwater), the thermal sensitivity coefficient is positive.

The geopotential height is defined such that ``Z(x, y) = 0`` at sea level and *decreases*
downwards to negative values, towards the bottom of the ocean.
"""
function thermal_sensitivity end

"""
    haline_sensitivity(Θ, Sᴬ, Z, equation_of_state)

Return the "Boussinesq haline sensitivity coefficient" for a seawater parcel with absolute
salinity `Sᴬ`, at fixed conservative temperature `Θ`, and geopotential height `Z`, using
the Boussinesq `equation_of_state`. The haline sensitivity coefficient is

```math
b(Θ, Sᴬ, Z) = \\left.\\frac{∂ρ}{∂Sᴬ}\\right|_{Θ, Z} ,
```

and measures how much seawater density changes when absolute salinity is changed.
It has units of [kg/m³/(g/kg)], and differs from `haline_contraction` (``β``) by a factor of the reference density ``ρᵣ``.
The sign of the haline sensitivity coefficient is defined to be consistent with the haline contraction coefficient.

The geopotential height is defined such that ``Z(x, y) = 0`` at sea level and *decreases*
downwards to negative values, towards the bottom of the ocean.
"""
function haline_sensitivity end

struct BoussinesqEquationOfState{P, FT}
    seawater_polynomial :: P
      reference_density :: FT

    function BoussinesqEquationOfState(seawater_polynomial, reference_density)
        FT = eltype(seawater_polynomial)
        P = typeof(seawater_polynomial)
        reference_density = convert(FT, reference_density)
        return new{P, FT}(seawater_polynomial, reference_density)
    end
end

const BEOS = BoussinesqEquationOfState

Base.eltype(::BEOS{<:Any, FT}) where FT = FT
with_float_type(FT, eos::BEOS) = BoussinesqEquationOfState(with_float_type(FT, eos.seawater_polynomial),
                                                           convert(FT, eos.reference_density))

Base.summary(eos::BoussinesqEquationOfState{P, FT}) where {P, FT} =
    string("BoussinesqEquationOfState{$FT}")

function Base.show(io::IO, eos::BoussinesqEquationOfState)
    print(io, summary(eos), ":", '\n')
    print(io, "├── seawater_polynomial: ", summary(eos.seawater_polynomial), '\n')
    print(io, "└── reference_density: ", eos.reference_density)
end

@inline reference_density(eos) = eos.reference_density

"""
    ρ(Θ, Sᴬ, Z, equation_of_state)

Return the total density of a seawater parcel with conservative temperature `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, using the Boussinesq `equation_of_state`.

The geopotential height is defined such that ``Z(x, y) = 0`` at sea level and *decreases*
downwards to negative values, towards the bottom of the ocean.

This function aliases `total_density`.
"""
@inline ρ(Θ, Sᴬ, Z, eos) = eos.reference_density + ρ′(Θ, Sᴬ, Z, eos)

"""
    ρ′(Θ, Sᴬ, Z, equation_of_state)

Return the density anomaly of a seawater parcel with conservative temperature `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, using the Boussinesq `equation_of_state`.

The geopotential height is defined such that ``Z(x, y) = 0`` at sea level and *decreases*
downwards to negative values, towards the bottom of the ocean.

The function aliases `density_anomaly`.
"""
function ρ′ end

const total_density = ρ
const density_anomaly = ρ′

"""
    thermal_expansion(Θ, Sᴬ, Z, equation_of_state)

Return the "Boussinesq thermal expansion coefficient" for a seawater parcel with reference
density `ρᵣ` and conservative temperture `Θ`, at fixed absolute salinity `Sᴬ` and
geopotential height `Z`, using the Boussinesq `equation_of_state`. The thermal expansion
coefficient is

```math
α(Θ, Sᴬ, Z) = - \\left.\\frac{1}{ρᵣ} \\frac{∂ρ}{∂Θ}\\right|_{Sᴬ, Z} ,
```

and describes seawater density changes due to changes in conservative temperature.
'Thermal expansion' is so named because, due to sign convention, positive values reflect decreasing
seawater density with increasing conservative temperature, and thus an 'expansion' of oceanic
fluid parcels. In many, but not all conditions in Earth's ocean (at temperatures greater than
4ᵒC in freshwater), the thermal expansion coefficient is positive.

The geopotential height is defined such that ``Z(x, y) = 0`` at sea level and *decreases*
downwards to negative values, towards the bottom of the ocean.
"""
@inline thermal_expansion(Θ, Sᴬ, Z, eos::BEOS{<:Any, FT}) where FT =
    thermal_sensitivity(FT(Θ), FT(Sᴬ), FT(Z), eos) / eos.reference_density

"""
    haline_contraction(Θ, Sᴬ, Z, equation_of_state)

Return the "Boussinesq haline contraction coefficient" for a seawater parcel with reference
density `ρᵣ` and absolute salinity `Sᴬ`, at fixed conservative temperture `Θ` and
geopotential height `Z`, using the Boussinesq `equation_of_state`. The haline contraction
coefficient is

```math
β(Θ, Sᴬ, Z) = \\left.\\frac{1}{ρᵣ} \\frac{∂ρ}{∂Sᴬ}\\right|_{Θ, Z},
```

and describes changes in seawater density due to changes in absolute salinity.
'Haline contraction' is so named because, due to sign convention, positive values reflect increasing
seawater density with increasing absolute salinity, and thus a slight 'contraction' of oceanic
fluid parcels.

The geopotential height is defined such that ``Z(x, y) = 0`` at sea level and *decreases*
downwards to negative values, towards the bottom of the ocean.
"""
@inline haline_contraction(Θ, Sᴬ, Z, eos::BEOS) = haline_sensitivity(Θ, Sᴬ, Z, eos) / eos.reference_density

"""
    AbstractSeawaterPolynomial

Abstract type for Boussinesq equations of state expressed as polynomial functions
of conservative temperature, absolute salinity, and geopotential depth.
"""
abstract type AbstractSeawaterPolynomial end

include("SecondOrderSeawaterPolynomials.jl")
include("TEOS10.jl")

using .SecondOrderSeawaterPolynomials
using .TEOS10

end # module
