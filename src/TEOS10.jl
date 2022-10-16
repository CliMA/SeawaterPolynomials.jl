module TEOS10

export 
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ′, thermal_sensitivity, haline_sensitivity

#####
##### The TEOS-10 polynomial approximation implemented in this file has been translated
##### into Julia from https://github.com/fabien-roquet/polyTEOS/blob/master/polyTEOS10.py
#####

"""
    struct TEOS10SeawaterPolynomial{FT} <: AbstractSeawaterPolynomial end

A 55-term polynomial approximation to the TEOS-10 standard equation of state for seawater.
"""
struct TEOS10SeawaterPolynomial{FT} <: AbstractSeawaterPolynomial end

Base.eltype(::TEOS10SeawaterPolynomial{FT}) where FT = FT
Base.summary(::TEOS10SeawaterPolynomial{FT}) where FT = "TEOS10SeawaterPolynomial{$FT}"

"""
    TEOS10SeawaterPolynomial(FT=Float64)

Returns an object representing a 55-term polynomial approximation to the TEOS-10 standard equation
of state for seawater. See

> Roquet et al., "Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard", Ocean Modelling (2015).
"""
TEOS10SeawaterPolynomial(FT=Float64) = TEOS10SeawaterPolynomial{FT}()

const EOS₁₀ = BoussinesqEquationOfState{<:TEOS10SeawaterPolynomial}

"""
    TEOS10EquationOfState(FT=Float64; reference_density=1020)

Returns an `BoussinesqEquationOfState` with a `TEOS10SeawaterPolynomial` of float type `FT`
with `reference density = 1020 kg m⁻³`, the value used by 

> Roquet et al., "Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard", Ocean Modelling (2015).

when fitting polynomial coefficients to the full TEOS-10 standard equation of state.
See the discussion prior to equation 8 in Roquet et al. (2015).

Note that according to Roquet et al. (2015):

> "In a Boussinesq model, the choice of the ρ₀ value is important, yet it varies significantly among OGCMs, as it is a matter of personal preference."
"""
TEOS10EquationOfState(FT=Float64; reference_density=FT(1020)) =
    BoussinesqEquationOfState(TEOS10SeawaterPolynomial{FT}(), reference_density)

#####
##### Reference values chosen using TEOS-10 recommendation
#####

const Sᵤ_64 = Float64(40 * 35.16504 / 35)
const Θᵤ_64 = Float64(40.0)
const Zᵤ_64 = Float64(1e4)
const ΔS_64 = Float64(32.0)

const Sᵤ_32 = Float32(40 * 35.16504 / 35)
const Θᵤ_32 = Float32(40.0)
const Zᵤ_32 = Float32(1e4)
const ΔS_32 = Float32(32.0)

#####
##### Coordinate transformations from (Θ, Sᴬ, p) to (τ, s, ζ)
#####

@inline τ(Θ::Float64) = Θ / Θᵤ_64
@inline s(Sᴬ::Float64) = √((Sᴬ + ΔS_64) / Sᵤ_64)
@inline ζ(Z::Float64) = - Z / Zᵤ_64

@inline τ(Θ::Float32) = Θ / Θᵤ_32
@inline s(Sᴬ::Float32) = √((Sᴬ + ΔS_32) / Sᵤ_32)
@inline ζ(Z::Float32) = - Z / Zᵤ_32

#####
##### Vertical reference profile of density
#####

const R₀₀_64 = Float64(+4.6494977072e+01) 
const R₀₁_64 = Float64(-5.2099962525e+00)
const R₀₂_64 = Float64(+2.2601900708e-01)
const R₀₃_64 = Float64(+6.4326772569e-02)
const R₀₄_64 = Float64(+1.5616995503e-02)
const R₀₅_64 = Float64(-1.7243708991e-03)

const R₀₀_32 = Float32(+4.3294977072e+01) 
const R₀₁_32 = Float32(-5.2099962525e+00)
const R₀₂_32 = Float32(+2.2601900708e-01)
const R₀₃_32 = Float32(+6.4326772569e-02)
const R₀₄_32 = Float32(+1.5616995503e-02)
const R₀₅_32 = Float32(-1.7243708991e-03)

@inline r₀(ζ::Float64) = (((((R₀₅_64 * ζ + R₀₄_64) * ζ + R₀₃_64) * ζ + R₀₂_64) * ζ + R₀₁_64) * ζ + R₀₀_64) * ζ
@inline r₀(ζ::Float32) = (((((R₀₅_32 * ζ + R₀₄_32) * ζ + R₀₃_32) * ζ + R₀₂_32) * ζ + R₀₁_32) * ζ + R₀₀_32) * ζ

#####
##### Density anomaly fit
#####

const R₀₀₀_64 = Float64(+8.0189615746e+02) 
const R₁₀₀_64 = Float64(+8.6672408165e+02)
const R₂₀₀_64 = Float64(-1.7864682637e+03)
const R₃₀₀_64 = Float64(+2.0375295546e+03)
const R₄₀₀_64 = Float64(-1.2849161071e+03)
const R₅₀₀_64 = Float64(+4.3227585684e+02)
const R₆₀₀_64 = Float64(-6.0579916612e+01)
const R₀₁₀_64 = Float64(+2.6010145068e+01)
const R₁₁₀_64 = Float64(-6.5281885265e+01)
const R₂₁₀_64 = Float64(+8.1770425108e+01)
const R₃₁₀_64 = Float64(-5.6888046321e+01)
const R₄₁₀_64 = Float64(+1.7681814114e+01)
const R₅₁₀_64 = Float64(-1.9193502195e+00)
const R₀₂₀_64 = Float64(-3.7074170417e+01)
const R₁₂₀_64 = Float64(+6.1548258127e+01)
const R₂₂₀_64 = Float64(-6.0362551501e+01)
const R₃₂₀_64 = Float64(+2.9130021253e+01)
const R₄₂₀_64 = Float64(-5.4723692739e+00)
const R₀₃₀_64 = Float64(+2.1661789529e+01)
const R₁₃₀_64 = Float64(-3.3449108469e+01)
const R₂₃₀_64 = Float64(+1.9717078466e+01)
const R₃₃₀_64 = Float64(-3.1742946532e+00)
const R₀₄₀_64 = Float64(-8.3627885467e+00)
const R₁₄₀_64 = Float64(+1.1311538584e+01)
const R₂₄₀_64 = Float64(-5.3563304045e+00)
const R₀₅₀_64 = Float64(+5.4048723791e-01)
const R₁₅₀_64 = Float64(+4.8169980163e-01)
const R₀₆₀_64 = Float64(-1.9083568888e-01)
const R₀₀₁_64 = Float64(+1.9681925209e+01)
const R₁₀₁_64 = Float64(-4.2549998214e+01)
const R₂₀₁_64 = Float64(+5.0774768218e+01)
const R₃₀₁_64 = Float64(-3.0938076334e+01)
const R₄₀₁_64 = Float64(+6.6051753097e+00)
const R₀₁₁_64 = Float64(-1.3336301113e+01)
const R₁₁₁_64 = Float64(-4.4870114575e+00)
const R₂₁₁_64 = Float64(+5.0042598061e+00)
const R₃₁₁_64 = Float64(-6.5399043664e-01)
const R₀₂₁_64 = Float64(+6.7080479603e+00)
const R₁₂₁_64 = Float64(+3.5063081279e+00)
const R₂₂₁_64 = Float64(-1.8795372996e+00)
const R₀₃₁_64 = Float64(-2.4649669534e+00)
const R₁₃₁_64 = Float64(-5.5077101279e-01)
const R₀₄₁_64 = Float64(+5.5927935970e-01)
const R₀₀₂_64 = Float64(+2.0660924175e+00)
const R₁₀₂_64 = Float64(-4.9527603989e+00)
const R₂₀₂_64 = Float64(+2.5019633244e+00)
const R₀₁₂_64 = Float64(+2.0564311499e+00)
const R₁₁₂_64 = Float64(-2.1311365518e-01)
const R₀₂₂_64 = Float64(-1.2419983026e+00)
const R₀₀₃_64 = Float64(-2.3342758797e-02)
const R₁₀₃_64 = Float64(-1.8507636718e-02)
const R₀₁₃_64 = Float64(+3.7969820455e-01)

const R₀₀₀_32 = Float32(+8.0189615746e+02) 
const R₁₀₀_32 = Float32(+8.6672408165e+02)
const R₂₀₀_32 = Float32(-1.7832682637e+03)
const R₃₀₀_32 = Float32(+2.0375295546e+03)
const R₄₀₀_32 = Float32(-1.2849161071e+03)
const R₅₀₀_32 = Float32(+4.3227585684e+02)
const R₆₀₀_32 = Float32(-6.0579916612e+01)
const R₀₁₀_32 = Float32(+2.6010145068e+01)
const R₁₁₀_32 = Float32(-6.5281885265e+01)
const R₂₁₀_32 = Float32(+8.1770425108e+01)
const R₃₁₀_32 = Float32(-5.6888046321e+01)
const R₄₁₀_32 = Float32(+1.7681814114e+01)
const R₅₁₀_32 = Float32(-1.9193502195e+00)
const R₀₂₀_32 = Float32(-3.7074170417e+01)
const R₁₂₀_32 = Float32(+6.1548258127e+01)
const R₂₂₀_32 = Float32(-6.0362551501e+01)
const R₃₂₀_32 = Float32(+2.9130021253e+01)
const R₄₂₀_32 = Float32(-5.4723692739e+00)
const R₀₃₀_32 = Float32(+2.1661789529e+01)
const R₁₃₀_32 = Float32(-3.3449108469e+01)
const R₂₃₀_32 = Float32(+1.9717078466e+01)
const R₃₃₀_32 = Float32(-3.1742946532e+00)
const R₀₄₀_32 = Float32(-8.3627885467e+00)
const R₁₄₀_32 = Float32(+1.1311538584e+01)
const R₂₄₀_32 = Float32(-5.3563304045e+00)
const R₀₅₀_32 = Float32(+5.4048723791e-01)
const R₁₅₀_32 = Float32(+4.8169980163e-01)
const R₀₆₀_32 = Float32(-1.9083568888e-01)
const R₀₀₁_32 = Float32(+1.9681925209e+01)
const R₁₀₁_32 = Float32(-4.2549998214e+01)
const R₂₀₁_32 = Float32(+5.0774768218e+01)
const R₃₀₁_32 = Float32(-3.0938076334e+01)
const R₄₀₁_32 = Float32(+6.6051753097e+00)
const R₀₁₁_32 = Float32(-1.3336301113e+01)
const R₁₁₁_32 = Float32(-4.4870114575e+00)
const R₂₁₁_32 = Float32(+5.0042598061e+00)
const R₃₁₁_32 = Float32(-6.5399043632e-01)
const R₀₂₁_32 = Float32(+6.7080479603e+00)
const R₁₂₁_32 = Float32(+3.5063081279e+00)
const R₂₂₁_32 = Float32(-1.8795372996e+00)
const R₀₃₁_32 = Float32(-2.4329669534e+00)
const R₁₃₁_32 = Float32(-5.5077101279e-01)
const R₀₄₁_32 = Float32(+5.5927935970e-01)
const R₀₀₂_32 = Float32(+2.0660924175e+00)
const R₁₀₂_32 = Float32(-4.9527603989e+00)
const R₂₀₂_32 = Float32(+2.5019633244e+00)
const R₀₁₂_32 = Float32(+2.0532311499e+00)
const R₁₁₂_32 = Float32(-2.1311365518e-01)
const R₀₂₂_32 = Float32(-1.2419983026e+00)
const R₀₀₃_32 = Float32(-2.3342758797e-02)
const R₁₀₃_32 = Float32(-1.8507636718e-02)
const R₀₁₃_32 = Float32(+3.7969820455e-01)

@inline r′₃(τ::Float64, s::Float64) = R₀₁₃_64 * τ + R₁₀₃_64 * s + R₀₀₃_64
@inline r′₂(τ::Float64, s::Float64) = (R₀₂₂_64 * τ + R₁₁₂_64 * s + R₀₁₂_64) * τ + (R₂₀₂_64 * s + R₁₀₂_64) * s + R₀₀₂_64

@inline r′₃(τ::Float32, s::Float32) = R₀₁₃_32 * τ + R₁₀₃_32 * s + R₀₀₃_32
@inline r′₂(τ::Float32, s::Float32) = (R₀₂₂_32 * τ + R₁₁₂_32 * s + R₀₁₂_32) * τ + (R₂₀₂_32 * s + R₁₀₂_32) * s + R₀₀₂_32

@inline r′₁(τ::Float64, s::Float64) =
    (((R₀₄₁_64 * τ + R₁₃₁_64  * s + R₀₃₁_64) * τ +
      (R₂₂₁_64 * s + R₁₂₁_64) * s + R₀₂₁_64) * τ +
     ((R₃₁₁_64 * s + R₂₁₁_64) * s + R₁₁₁_64) * s + R₀₁₁_64) * τ +
    (((R₄₀₁_64 * s + R₃₀₁_64) * s + R₂₀₁_64) * s + R₁₀₁_64) * s + R₀₀₁_64

@inline r′₁(τ::Float32, s::Float32) =
    (((R₀₄₁_32 * τ + R₁₃₁_32  * s + R₀₃₁_32) * τ +
      (R₂₂₁_32 * s + R₁₂₁_32) * s + R₀₂₁_32) * τ +
     ((R₃₁₁_32 * s + R₂₁₁_32) * s + R₁₁₁_32) * s + R₀₁₁_32) * τ +
    (((R₄₀₁_32 * s + R₃₀₁_32) * s + R₂₀₁_32) * s + R₁₀₁_32) * s + R₀₀₁_32

@inline r′₀(τ::Float64, s::Float64) =
    (((((R₀₆₀_64 * τ + R₁₅₀_64  * s + R₀₅₀_64) * τ +
        (R₂₄₀_64 * s + R₁₄₀_64) * s + R₀₄₀_64) * τ +
       ((R₃₃₀_64 * s + R₂₃₀_64) * s + R₁₃₀_64) * s + R₀₃₀_64) * τ +
      (((R₄₂₀_64 * s + R₃₂₀_64) * s + R₂₂₀_64) * s + R₁₂₀_64) * s + R₀₂₀_64) * τ +
     ((((R₅₁₀_64 * s + R₄₁₀_64) * s + R₃₁₀_64) * s + R₂₁₀_64) * s + R₁₁₀_64) * s + R₀₁₀_64) * τ +
    (((((R₆₀₀_64 * s + R₅₀₀_64) * s + R₄₀₀_64) * s + R₃₀₀_64) * s + R₂₀₀_64) * s + R₁₀₀_64) * s + R₀₀₀_64

@inline r′₀(τ::Float32, s::Float32) =
    (((((R₀₆₀_32 * τ + R₁₅₀_32  * s + R₀₅₀_32) * τ +
        (R₂₄₀_32 * s + R₁₄₀_32) * s + R₀₄₀_32) * τ +
       ((R₃₃₀_32 * s + R₂₃₀_32) * s + R₁₃₀_32) * s + R₀₃₀_32) * τ +
      (((R₄₂₀_32 * s + R₃₂₀_32) * s + R₂₂₀_32) * s + R₁₂₀_32) * s + R₀₂₀_32) * τ +
     ((((R₅₁₀_32 * s + R₄₁₀_32) * s + R₃₁₀_32) * s + R₂₁₀_32) * s + R₁₁₀_32) * s + R₀₁₀_32) * τ +
    (((((R₆₀₀_32 * s + R₅₀₀_32) * s + R₄₀₀_32) * s + R₃₀₀_32) * s + R₂₀₀_32) * s + R₁₀₀_32) * s + R₀₀₀_32

@inline r′(τ, s, ζ) = ((r′₃(τ, s) * ζ + r′₂(τ, s)) * ζ + r′₁(τ, s)) * ζ + r′₀(τ, s)

#####
##### Density perturbation
#####

"""
    ρ(Θ, Sᴬ, Z, ::BoussinesqEquationOfState{<:TEOS10EquationOfState})

Returns the in-situ density of seawater with state `(Θ, Sᴬ, Z)` using the 55-term polynomial
approximation to TEOS-10 described in Roquet et al. (§3.1, 2014).

# Inputs
    Θ: conservative temperature (ITS-90) [°C]
    Sᴬ: absolute salinity [g/kg]
    Z: geopotential depth [m]

# Output
    ρ: in-situ density [kg/m³]

# References
- Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate polynomial expressions
  for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling.
"""
@inline ρ(Θ, Sᴬ, Z, ::EOS₁₀) = _ρ(τ(Θ), s(Sᴬ), ζ(Z))

@inline ρ′(Θ, Sᴬ, Z, eos::EOS₁₀) = ρ(Θ, Sᴬ, Z, eos) - eos.reference_density 

@inline _ρ(τ, s, ζ) = r₀(ζ) + r′(τ, s, ζ)

#####
##### Thermal expansion fit
#####

const α₀₀₀_64 = Float64(-6.5025362670e-01) 
const α₁₀₀_64 = Float64(+1.6320471316e+00)
const α₂₀₀_64 = Float64(-2.0442606277e+00)
const α₃₀₀_64 = Float64(+1.4222011580e+00)
const α₄₀₀_64 = Float64(-4.4204535284e-01)
const α₅₀₀_64 = Float64(+4.7983755487e-02)
const α₀₁₀_64 = Float64(+1.8537085209e+00)
const α₁₁₀_64 = Float64(-3.0774129064e+00)
const α₂₁₀_64 = Float64(+3.0181275751e+00)
const α₃₁₀_64 = Float64(-1.4565010626e+00)
const α₄₁₀_64 = Float64(+2.7361846370e-01)
const α₀₂₀_64 = Float64(-1.6246342147e+00)
const α₁₂₀_64 = Float64(+2.5086831352e+00)
const α₂₂₀_64 = Float64(-1.4787808849e+00)
const α₃₂₀_64 = Float64(+2.3807209899e-01)
const α₀₃₀_64 = Float64(+8.3627885467e-01)
const α₁₃₀_64 = Float64(-1.1311538584e+00)
const α₂₃₀_64 = Float64(+5.3563304045e-01)
const α₀₄₀_64 = Float64(-6.7560904739e-02)
const α₁₄₀_64 = Float64(-6.0212475204e-02)
const α₀₅₀_64 = Float64(+2.8625353333e-02)
const α₀₀₁_64 = Float64(+3.3340752782e-01)
const α₁₀₁_64 = Float64(+1.1217528644e-01)
const α₂₀₁_64 = Float64(-1.2510649515e-01)
const α₃₀₁_64 = Float64(+1.6349760916e-02)
const α₀₁₁_64 = Float64(-3.3540239802e-01)
const α₁₁₁_64 = Float64(-1.7531540640e-01)
const α₂₁₁_64 = Float64(+9.3976864981e-02)
const α₀₂₁_64 = Float64(+1.8487252150e-01)
const α₁₂₁_64 = Float64(+4.1307825959e-02)
const α₀₃₁_64 = Float64(-5.5927935970e-02)
const α₀₀₂_64 = Float64(-5.1410778748e-02)
const α₁₀₂_64 = Float64(+5.3278413794e-03)
const α₀₁₂_64 = Float64(+6.2099915132e-02)
const α₀₀₃_64 = Float64(-9.4924551138e-03)

const α₀₀₀_32 = Float32(-6.5025362670e-01) 
const α₁₀₀_32 = Float32(+1.6320471316e+00)
const α₂₀₀_32 = Float32(-2.0442606277e+00)
const α₃₀₀_32 = Float32(+1.4222011580e+00)
const α₄₀₀_32 = Float32(-4.4204535284e-01)
const α₅₀₀_32 = Float32(+4.7983755487e-02)
const α₀₁₀_32 = Float32(+1.8537085209e+00)
const α₁₁₀_32 = Float32(-3.0774129032e+00)
const α₂₁₀_32 = Float32(+3.0181275751e+00)
const α₃₁₀_32 = Float32(-1.4565010626e+00)
const α₄₁₀_32 = Float32(+2.7361846370e-01)
const α₀₂₀_32 = Float32(-1.6246342147e+00)
const α₁₂₀_32 = Float32(+2.5086831352e+00)
const α₂₂₀_32 = Float32(-1.4787808849e+00)
const α₃₂₀_32 = Float32(+2.3807209899e-01)
const α₀₃₀_32 = Float32(+8.3627885467e-01)
const α₁₃₀_32 = Float32(-1.1311538584e+00)
const α₂₃₀_32 = Float32(+5.3563304045e-01)
const α₀₄₀_32 = Float32(-6.7560904739e-02)
const α₁₄₀_32 = Float32(-6.0212475204e-02)
const α₀₅₀_32 = Float32(+2.8625353333e-02)
const α₀₀₁_32 = Float32(+3.3340752782e-01)
const α₁₀₁_32 = Float32(+1.1217528324e-01)
const α₂₀₁_32 = Float32(-1.2510329515e-01)
const α₃₀₁_32 = Float32(+1.6349760916e-02)
const α₀₁₁_32 = Float32(-3.3540239802e-01)
const α₁₁₁_32 = Float32(-1.7531540320e-01)
const α₂₁₁_32 = Float32(+9.3976832981e-02)
const α₀₂₁_32 = Float32(+1.8487252150e-01)
const α₁₂₁_32 = Float32(+4.1307825959e-02)
const α₀₃₁_32 = Float32(-5.5927935970e-02)
const α₀₀₂_32 = Float32(-5.1410778748e-02)
const α₁₀₂_32 = Float32(+5.3278413794e-03)
const α₀₁₂_32 = Float32(+6.2099915132e-02)
const α₀₀₃_32 = Float32(-9.4924551138e-03)

"""
    thermal_sensitivity(Θ, Sᴬ, Z, ::TEOS10)

Returns the Boussinesq thermal expansion coefficient ``-∂ρ/∂Θ`` [kg/m³/K] computed using
the 55-term polynomial approximation to TEOS-10 described in Roquet et al. (§3.1, 2014).

# Inputs
    Θ: conservative temperature (ITS-90) [°C]
    Sᴬ: absolute salinity [g/kg]
    Z: geopotential depth [m]

# Output
    a: Boussinesq thermal expansion coefficient -∂ρ/∂Θ [kg/m³/K]

# References
- Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate polynomial expressions
  for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling.
"""
@inline thermal_sensitivity(Θ, Sᴬ, Z, ::EOS₁₀) = _a(τ(Θ), s(Sᴬ), ζ(Z))

@inline _a(τ::Float64, s::Float64, ζ::Float64) =
      ((α₀₀₃_64 * ζ + α₀₁₂_64  * τ + α₁₀₂_64  * s + α₀₀₂_64) * ζ +
      ((α₀₃₁_64 * τ + α₁₂₁_64  * s + α₀₂₁_64) * τ +
       (α₂₁₁_64 * s + α₁₁₁_64) * s + α₀₁₁_64) * τ +
      ((α₃₀₁_64 * s + α₂₀₁_64) * s + α₁₀₁_64) * s + α₀₀₁_64) * ζ +
    ((((α₀₅₀_64 * τ + α₁₄₀_64  * s + α₀₄₀_64) * τ +
       (α₂₃₀_64 * s + α₁₃₀_64) * s + α₀₃₀_64) * τ +
      ((α₃₂₀_64 * s + α₂₂₀_64) * s + α₁₂₀_64) * s + α₀₂₀_64) * τ +
     (((α₄₁₀_64 * s + α₃₁₀_64) * s + α₂₁₀_64) * s + α₁₁₀_64) * s + α₀₁₀_64) * τ +
    ((((α₅₀₀_64 * s + α₄₀₀_64) * s + α₃₀₀_64) * s + α₂₀₀_64) * s + α₁₀₀_64) * s + α₀₀₀_64

@inline _a(τ::Float32, s::Float32, ζ::Float32) =
      ((α₀₀₃_32 * ζ + α₀₁₂_32  * τ + α₁₀₂_32  * s + α₀₀₂_32) * ζ +
      ((α₀₃₁_32 * τ + α₁₂₁_32  * s + α₀₂₁_32) * τ +
       (α₂₁₁_32 * s + α₁₁₁_32) * s + α₀₁₁_32) * τ +
      ((α₃₀₁_32 * s + α₂₀₁_32) * s + α₁₀₁_32) * s + α₀₀₁_32) * ζ +
    ((((α₀₅₀_32 * τ + α₁₄₀_32  * s + α₀₄₀_32) * τ +
       (α₂₃₀_32 * s + α₁₃₀_32) * s + α₀₃₀_32) * τ +
      ((α₃₂₀_32 * s + α₂₂₀_32) * s + α₁₂₀_32) * s + α₀₂₀_32) * τ +
     (((α₄₁₀_32 * s + α₃₁₀_32) * s + α₂₁₀_32) * s + α₁₁₀_32) * s + α₀₁₀_32) * τ +
    ((((α₅₀₀_32 * s + α₄₀₀_32) * s + α₃₀₀_32) * s + α₂₀₀_32) * s + α₁₀₀_32) * s + α₀₀₀_32

#####
##### Haline contraction
#####

const β₀₀₀_64 = Float64(+1.0783203594e+01) 
const β₁₀₀_64 = Float64(-4.4452095908e+01)
const β₂₀₀_64 = Float64(+7.6048755820e+01)
const β₃₀₀_64 = Float64(-6.3944280668e+01)
const β₄₀₀_64 = Float64(+2.6890441098e+01)
const β₅₀₀_64 = Float64(-4.5221697773e+00)
const β₀₁₀_64 = Float64(-8.1219372432e-01)
const β₁₁₀_64 = Float64(+2.0346663041e+00)
const β₂₁₀_64 = Float64(-2.1232895170e+00)
const β₃₁₀_64 = Float64(+8.7994140485e-01)
const β₄₁₀_64 = Float64(-1.1939638360e-01)
const β₀₂₀_64 = Float64(+7.6574242289e-01)
const β₁₂₀_64 = Float64(-1.5019813020e+00)
const β₂₂₀_64 = Float64(+1.0872489522e+00)
const β₃₂₀_64 = Float64(-2.7233429080e-01)
const β₀₃₀_64 = Float64(-4.1615152308e-01)
const β₁₃₀_64 = Float64(+4.9061350869e-01)
const β₂₃₀_64 = Float64(-1.1847737788e-01)
const β₀₄₀_64 = Float64(+1.4073062708e-01)
const β₁₄₀_64 = Float64(-1.3327978879e-01)
const β₀₅₀_64 = Float64(+5.9929880134e-03)
const β₀₀₁_64 = Float64(-5.2937873009e-01)
const β₁₀₁_64 = Float64(+1.2634116779e+00)
const β₂₀₁_64 = Float64(-1.1547328025e+00)
const β₃₀₁_64 = Float64(+3.2870876279e-01)
const β₀₁₁_64 = Float64(-5.5824407214e-02)
const β₁₁₁_64 = Float64(+1.2451933313e-01)
const β₂₁₁_64 = Float64(-2.4409539932e-02)
const β₀₂₁_64 = Float64(+4.3623149752e-02)
const β₁₂₁_64 = Float64(-4.6767901790e-02)
const β₀₃₁_64 = Float64(-6.8523260060e-03)
const β₀₀₂_64 = Float64(-6.1618945251e-02)
const β₁₀₂_64 = Float64(+6.2255521644e-02)
const β₀₁₂_64 = Float64(-2.6514181169e-03)
const β₀₀₃_64 = Float64(-2.3025968587e-04)

const β₀₀₀_32 = Float32(+1.0783203594e+01) 
const β₁₀₀_32 = Float32(-4.4452095908e+01)
const β₂₀₀_32 = Float32(+7.6048755820e+01)
const β₃₀₀_32 = Float32(-6.3944280668e+01)
const β₄₀₀_32 = Float32(+2.6890441098e+01)
const β₅₀₀_32 = Float32(-4.5221697773e+00)
const β₀₁₀_32 = Float32(-8.1219372432e-01)
const β₁₁₀_32 = Float32(+2.0346663041e+00)
const β₂₁₀_32 = Float32(-2.1232895170e+00)
const β₃₁₀_32 = Float32(+8.7994140485e-01)
const β₄₁₀_32 = Float32(-1.1939638360e-01)
const β₀₂₀_32 = Float32(+7.6574242289e-01)
const β₁₂₀_32 = Float32(-1.5019813020e+00)
const β₂₂₀_32 = Float32(+1.0872489522e+00)
const β₃₂₀_32 = Float32(-2.7233429080e-01)
const β₀₃₀_32 = Float32(-4.1615152308e-01)
const β₁₃₀_32 = Float32(+4.9061350869e-01)
const β₂₃₀_32 = Float32(-1.1847737788e-01)
const β₀₄₀_32 = Float32(+1.4073062708e-01)
const β₁₄₀_32 = Float32(-1.3327978879e-01)
const β₀₅₀_32 = Float32(+5.9929880134e-03)
const β₀₀₁_32 = Float32(-5.2937873009e-01)
const β₁₀₁_32 = Float32(+1.2634116779e+00)
const β₂₀₁_32 = Float32(-1.1547328025e+00)
const β₃₀₁_32 = Float32(+3.2870876279e-01)
const β₀₁₁_32 = Float32(-5.5824407214e-02)
const β₁₁₁_32 = Float32(+1.2451933313e-01)
const β₂₁₁_32 = Float32(-2.4409539932e-02)
const β₀₂₁_32 = Float32(+4.3623149752e-02)
const β₁₂₁_32 = Float32(-4.6767901790e-02)
const β₀₃₁_32 = Float32(-6.8523260060e-03)
const β₀₀₂_32 = Float32(-6.1618945251e-02)
const β₁₀₂_32 = Float32(+6.2255521324e-02)
const β₀₁₂_32 = Float32(-2.6514181169e-03)
const β₀₀₃_32 = Float32(-2.3025968587e-04)

"""
    haline_sensitivity(Θ, Sᴬ, Z, ::TEOS10)

Returns the Boussinesq haline contraction coefficient ``∂ρ/∂Sᴬ`` [kg/m³/(g/kg)] computed using
the 55-term polynomial approximation to TEOS-10 described in Roquet et al. (§3.1, 2014).

# Inputs
    Θ: conservative temperature (ITS-90) [°C]
    Sᴬ: absolute salinity [g/kg]
    Z: geopotential depth [m]

# Output
    β: Boussinesq haline contraction coefficient ∂ρ/∂Sᴬ [kg/m³/(g/kg)]

# References
- Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate polynomial expressions
  for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling.
"""
@inline haline_sensitivity(Θ, Sᴬ, Z, ::EOS₁₀) = _b(τ(Θ), s(Sᴬ), ζ(Z)) / s(Sᴬ)

@inline _b(τ::Float64, s::Float64, ζ::Float64) =
      ((β₀₀₃_64 * ζ + β₀₁₂_64  * τ + β₁₀₂_64  * s + β₀₀₂_64) * ζ +
      ((β₀₃₁_64 * τ + β₁₂₁_64  * s + β₀₂₁_64) * τ +
       (β₂₁₁_64 * s + β₁₁₁_64) * s + β₀₁₁_64) * τ +
      ((β₃₀₁_64 * s + β₂₀₁_64) * s + β₁₀₁_64) * s + β₀₀₁_64) * ζ +
    ((((β₀₅₀_64 * τ + β₁₄₀_64  * s + β₀₄₀_64) * τ +
       (β₂₃₀_64 * s + β₁₃₀_64) * s + β₀₃₀_64) * τ +
      ((β₃₂₀_64 * s + β₂₂₀_64) * s + β₁₂₀_64) * s + β₀₂₀_64) * τ +
     (((β₄₁₀_64 * s + β₃₁₀_64) * s + β₂₁₀_64) * s + β₁₁₀_64) * s + β₀₁₀_64) * τ +
    ((((β₅₀₀_64 * s + β₄₀₀_64) * s + β₃₀₀_64) * s + β₂₀₀_64) * s + β₁₀₀_64) * s + β₀₀₀_64

@inline _b(τ::Float32, s::Float32, ζ::Float32) =
      ((β₀₀₃_32 * ζ + β₀₁₂_32  * τ + β₁₀₂_32  * s + β₀₀₂_32) * ζ +
      ((β₀₃₁_32 * τ + β₁₂₁_32  * s + β₀₂₁_32) * τ +
       (β₂₁₁_32 * s + β₁₁₁_32) * s + β₀₁₁_32) * τ +
      ((β₃₀₁_32 * s + β₂₀₁_32) * s + β₁₀₁_32) * s + β₀₀₁_32) * ζ +
    ((((β₀₅₀_32 * τ + β₁₄₀_32  * s + β₀₄₀_32) * τ +
       (β₂₃₀_32 * s + β₁₃₀_32) * s + β₀₃₀_32) * τ +
      ((β₃₂₀_32 * s + β₂₂₀_32) * s + β₁₂₀_32) * s + β₀₂₀_32) * τ +
     (((β₄₁₀_32 * s + β₃₁₀_32) * s + β₂₁₀_32) * s + β₁₁₀_32) * s + β₀₁₀_32) * τ +
    ((((β₅₀₀_32 * s + β₄₀₀_32) * s + β₃₀₀_32) * s + β₂₀₀_32) * s + β₁₀₀_32) * s + β₀₀₀_32

end # module
