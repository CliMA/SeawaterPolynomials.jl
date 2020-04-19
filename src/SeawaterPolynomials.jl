module SeawaterPolynomials

"""
    AbstractSeawaterPolynomial

Abstract type for Boussinesq equations of state expressed as polynomial functions
of conservative temperature, absolute salinity, and geopotential depth.
"""
abstract type AbstractSeawaterPolynomial end

"""
    ρ(Θ, Sᴬ, Z, equation_of_state)

Returns the total density of a seawater parcel that has conservative tempertuare `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, using the Boussinesq `equation_of_state`.

The geopotential height Z(x, y) is defined such that `Z(x, y) = 0` at sea level and 
*decreases* downwards to negative values, towards the bottom of the ocean.

This function aliases `total_density`.
"""
function ρ end

"""
    ρ′(Θ, Sᴬ, Z, equation_of_state)

Returns the total density of a seawater parcel that has conservative tempertuare `Θ`, absolute
salinity `Sᴬ`, at the geopotential height `Z`, using the Boussinesq `equation_of_state`.

The geopotential height Z(x, y) is defined such that `Z(x, y) = 0` at sea level and 
*decreases* downwards to negative values, towards the bottom of the ocean.

The function aliases `density_anomaly`.
"""
function ρ′ end

const total_density = ρ
const density_anomaly = ρ′

struct BoussinesqEquationOfState{P, FT}
    seawater_polynomial :: P
      reference_density :: FT
end

include("SecondOrderSeawaterPolynomials.jl")

end # module
