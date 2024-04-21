module SecondOrderSeawaterPolynomials

export
    SecondOrderSeawaterPolynomial,
    RoquetSeawaterPolynomial,
    RoquetEquationOfState

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ′, thermal_sensitivity, haline_sensitivity

"""
    struct SecondOrderSeawaterPolynomial{FT} <: AbstractSeawaterPolynomial

Container of coefficients for a second-order polynomial function of
absolute salinity `Sᴬ`, conservative temperature `Θ`, and geopotential
depth `Z` for seawater density.

The coefficients have the form

```math
Rᵦᵪᵩ ,
```

where ``β, χ, φ`` denote the order of the term to which the coefficent corresponds:
``β`` is the polynomial order of absolute salinity, ``Sᴬ``, ``χ`` is the polynomial order
of conservative temperature, ``Θ``, and ``φ`` is the order of geopotential height, ``Z``.

For a `SecondOrderSeawaterPolynomial`, ``β + χ + φ < 3``.

The coefficient ``R₁₁₀`` arises in the seawater polynomial as

```
seawater_polynomial(Θ, Sᴬ, Z) = ⋯ + R₁₁₀ * Sᴬ * Θ + ⋯
```
"""
@Base.kwdef struct SecondOrderSeawaterPolynomial{FT} <: AbstractSeawaterPolynomial
    R₁₀₀ :: FT = 0
    R₀₁₀ :: FT = 0
    R₁₀₁ :: FT = 0
    R₀₁₁ :: FT = 0
    R₁₁₀ :: FT = 0
    R₀₂₀ :: FT = 0
    R₂₀₀ :: FT = 0
end

const EOS₂ = BoussinesqEquationOfState{<:SecondOrderSeawaterPolynomial}

Base.eltype(::SecondOrderSeawaterPolynomial{FT}) where FT = FT
Base.summary(::SecondOrderSeawaterPolynomial{FT}) where FT = "SecondOrderSeawaterPolynomial{$FT}"

signstr(x) = sign(x) < 0 ? " - " : " + "

function Base.show(io::IO, eos::SecondOrderSeawaterPolynomial)
    print(io, "ρ' = ")
    print(io, eos.R₁₀₀, " Sᴬ")
    print(io, signstr(eos.R₀₁₀), abs(eos.R₀₁₀), " Θ")
    print(io, signstr(eos.R₀₂₀), abs(eos.R₀₂₀), " Θ²")
    print(io, signstr(eos.R₀₁₁), abs(eos.R₀₁₁), " Θ Z")
    print(io, signstr(eos.R₂₀₀), abs(eos.R₂₀₀), " Sᴬ²")
    print(io, signstr(eos.R₁₀₁), abs(eos.R₁₀₁), " Sᴬ Z")
    print(io, signstr(eos.R₁₁₀), abs(eos.R₁₁₀), " Sᴬ Θ")
end

@inline ρ′(Θ, Sᴬ, Z, eos::EOS₂) = (  eos.seawater_polynomial.R₁₀₀ * Sᴬ
                                   + eos.seawater_polynomial.R₀₁₀ * Θ
                                   + eos.seawater_polynomial.R₀₂₀ * Θ^2
                                   - eos.seawater_polynomial.R₀₁₁ * Θ * Z
                                   + eos.seawater_polynomial.R₂₀₀ * Sᴬ^2
                                   - eos.seawater_polynomial.R₁₀₁ * Sᴬ * Z
                                   + eos.seawater_polynomial.R₁₁₀ * Sᴬ * Θ )

@inline thermal_sensitivity(Θ, Sᴬ, Z, eos::EOS₂) = (      eos.seawater_polynomial.R₀₁₀
                                                    + 2 * eos.seawater_polynomial.R₀₂₀ * Θ
                                                    -     eos.seawater_polynomial.R₀₁₁ * Z
                                                    +     eos.seawater_polynomial.R₁₁₀ * Sᴬ )

@inline haline_sensitivity(Θ, Sᴬ, Z, eos::EOS₂) = (      eos.seawater_polynomial.R₁₀₀
                                                   + 2 * eos.seawater_polynomial.R₂₀₀ * Sᴬ
                                                   -     eos.seawater_polynomial.R₁₀₁ * Z
                                                   +     eos.seawater_polynomial.R₁₁₀ * Θ )

"""
    RoquetSeawaterPolynomial([FT=Float64,] coefficient_set=:SecondOrder)

Return a `SecondOrderSeawaterPolynomial` with coefficients optimized by

> Roquet et al., "Defining a Simplified yet 'Realistic' Equation of State for Seawater", Journal of Physical Oceanography (2015).

The `coefficient_set` is a symbol or string that selects one of the "sets" of
optimized second order coefficients.

Coefficient sets
================

- `:Linear`: a linear equation of state, ``ρ = ρᵣ + R₁₀₀ Θ + R₀₁₀ Sᴬ``.

- `:Cabbeling`: includes quadratic temperature term, ``ρ = ρᵣ + R₀₁₀ Θ + R₁₀₀ Sᴬ + R₀₂₀ Θ²``.

- `:CabbelingThermobaricity`: includes 'thermobaricity' term,
                              ``ρ = ρᵣ + R₀₁₀ Θ + R₁₀₀ Sᴬ + R₀₂₀ Θ² - R₀₁₁ Θ Z``.

- `:Freezing`: same as `:cabbeling_thermobaricity` with modified constants to increase
               accuracy near freezing.

- `:SecondOrder`: includes quadratic salinity, halibaricity, and thermohaline term,
                  ``ρ = ρᵣ + R₁₀₀ Sᴬ + R₀₁₀ Θ + R₀₂₀ Θ² - R₀₁₁ Θ Z
                           + R₂₀₀ (Sᴬ)² - R₁₀₁ Sᴬ Z + R₁₁₀ Sᴬ Θ``.

- `:SimplestRealistic`: the proposed simplest though "realistic" equation of state for
                        seawater from Rouquet et al. (2015),
                        ``ρ = ρᵣ + R₁₀₀ Sᴬ  + R₀₁₀ Θ - R₀₂₀ Θ² - R₀₁₁ Θ Z``

The optimized coefficients are reported in Table 3 of Roquet et al., "Defining a Simplified
yet 'Realistic' Equation of State for Seawater", Journal of Physical Oceanography (2015), and
further discussed around equations (12)--(15). The optimization minimizes errors in estimated
horizontal density gradient estiamted from climatological temperature and salinity distributions
between the 5 simplified forms chosen by Roquet et. al and the full-fledged
[TEOS-10](http://www.teos-10.org) equation of state.
The `:SimplestRealistic` equation of state is equation (17) in Rouquet et al. (2015) which
they propose is the simplest yet "realistic" form for the equation of state for the density
of seawater.
**Note:** In equation (17) from Rouquet et al. (2015), the `:SimplestRealistic` equation
of state has an extra term ``R₀₀₀ = -C_{b}Θ₀²/2`` which is not included in the coefficient set
above. This is because this term has no effect on ocean dynamics.
"""
RoquetSeawaterPolynomial(FT::DataType, coefficient_set=:SecondOrder) =
    eval(Symbol(coefficient_set, :RoquetSeawaterPolynomial))(FT)

# For Float64 default
RoquetSeawaterPolynomial(coefficient_set=:SecondOrder) =
    eval(Symbol(coefficient_set, :RoquetSeawaterPolynomial))(Float64)

"""
    RoquetEquationOfState([FT=Float64,] coefficient_set=:SecondOrder; reference_density=1024.6)

Return an `BoussinesqEquationOfState` with a `RoquetSeawaterPolynomial` corresponding to
`coefficient_set` and with `reference density = 1024.6 kg m⁻³`, the average density of seawater
at the surface of the world ocean.

See [`RoquetSeawaterPolynomial`](@ref) for options for the `coefficient_set`.
The optimzed coefficient sets for the `RoquetSeawaterPolynomial` are listed in Table 3 in

> Roquet et al., "Defining a Simplified yet 'Realistic' Equation of State for Seawater", Journal of Physical Oceanography (2015).
"""
RoquetEquationOfState(FT::DataType, coefficient_set=:SecondOrder; reference_density=1024.6) =
    BoussinesqEquationOfState(RoquetSeawaterPolynomial(FT, coefficient_set), reference_density)

RoquetEquationOfState(coefficient_set=:SecondOrder; reference_density=1024.6) =
    RoquetEquationOfState(Float64, coefficient_set; reference_density)

"""
    LinearRoquetSeawaterPolynomial([FT=Float64])

Parameters for a linear equation of state optimized for the 'current' oceanic
temperature and salinity distribution.

For more information see [`RoquetSeawaterPolynomial`](@ref).
"""
LinearRoquetSeawaterPolynomial(FT=Float64) =
    SecondOrderSeawaterPolynomial{FT}(R₀₁₀ = - 1.775e-1,
                                      R₁₀₀ =   7.718e-1)

"""
    CabbelingRoquetSeawaterPolynomial([FT=Float64])

Parameters for a minimal equation of state that describes cabbeling,
optimized for the 'current' oceanic temperature and salinity distribution.

For more information see [`RoquetSeawaterPolynomial`](@ref).
"""
CabbelingRoquetSeawaterPolynomial(FT=Float64) =
    SecondOrderSeawaterPolynomial{FT}(R₀₁₀ = - 0.844e-1,
                                      R₁₀₀ =   7.718e-1,
                                      R₀₂₀ = - 4.561e-3)

"""
    CabbelingThermobaricityRoquetSeawaterPolynomial([FT=Float64])

Parameters for a minimal equation of state that describes cabbeling and thermobaric
effects on sewater density, optimized for the 'current' oceanic temperature and salinity
distribution.

For more information see [`RoquetSeawaterPolynomial`](@ref).
"""
CabbelingThermobaricityRoquetSeawaterPolynomial(FT=Float64) =
    SecondOrderSeawaterPolynomial{FT}(R₀₁₀ = - 0.651e-1,
                                      R₁₀₀ =   7.718e-1,
                                      R₀₂₀ = - 5.027e-3,
                                      R₀₁₁ = - 2.5681e-5)

"""
    FreezingRoquetSeawaterPolynomial(FT=Float64)

Parameters for a minimal equation of state that describes seawater density near its
freezing point, optimized for the 'current' oceanic temperature and salinity
distribution.

For more information see [`RoquetSeawaterPolynomial`](@ref).
"""
FreezingRoquetSeawaterPolynomial(FT=Float64) =
    SecondOrderSeawaterPolynomial{FT}(R₀₁₀ = - 0.491e-1,
                                      R₁₀₀ =   7.718e-1,
                                      R₀₂₀ = - 5.027e-3,
                                      R₀₁₁ = - 2.5681e-5)

"""
    SecondOrderRoquetSeawaterPolynomial([FT=Float64])

Parameters for a fully second-order equation of state for seawater,
optimized for the 'current' oceanic temperature and salinity  distribution.

For more information see [`RoquetSeawaterPolynomial`](@ref).
"""
SecondOrderRoquetSeawaterPolynomial(FT=Float64) =
    SecondOrderSeawaterPolynomial{FT}(R₀₁₀ =   0.182e-1,
                                      R₁₀₀ =   8.078e-1,
                                      R₀₂₀ = - 4.937e-3,
                                      R₀₁₁ = - 2.4677e-5,
                                      R₂₀₀ = - 1.115e-4,
                                      R₁₀₁ = - 8.241e-6,
                                      R₁₁₀ = - 2.446e-3)

"""
    SimplestRealisticRoquetSeawaterPolynomial([FT=Float64])

Parameters for the simplest yet "realistic" equation of state for seawater
from Roquet et al. (2015) (see equation (17)), optimized for the 'current' oceanic
temperature and salinity distribution.

For more information see [`RoquetSeawaterPolynomial`](@ref).
"""
function SimplestRealisticRoquetSeawaterPolynomial(FT=Float64)

    Cb = 0.011   # kg m⁻³ K⁻²
    Tₕ = 2.5e-5  # kg m⁻⁴ K⁻¹
    b₀ = 0.77    # kg m⁻³ (g kg⁻¹)⁻¹
    Θ₀ = -4.5    # °C

    return SecondOrderSeawaterPolynomial{FT}(R₁₀₀ = b₀,
                                             R₀₁₀ = Cb * Θ₀,
                                             R₀₂₀ = -Cb / 2,
                                             R₀₁₁ = -Tₕ)
end

end # module
