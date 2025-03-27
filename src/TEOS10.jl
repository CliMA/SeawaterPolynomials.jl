module TEOS10

export 
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ, ρ′, thermal_sensitivity, haline_sensitivity

#####
##### The TEOS-10 polynomial approximation implemented in this file has been translated
##### into Julia from https://github.com/fabien-roquet/polyTEOS/blob/master/polyTEOS10.py
#####
##### Note: this implementation codes the polynomial weights as `const`, thus capturing the
##### polynomial weights implicitly in function closures. As a result of this design, the
##### polynomial weights cannot be modified without changing this file. To implement a
##### polynomial approximation with settable weights, the struct `TEOS10SeawaterPolynomial`
##### should be modified to carry the weights explicitly.
#####

"""
    struct TEOS10SeawaterPolynomial{FT} <: AbstractSeawaterPolynomial end

A 55-term polynomial approximation to the TEOS-10 standard equation of state for seawater.
"""
struct TEOS10SeawaterPolynomial{FT} <: AbstractSeawaterPolynomial end

# The constant, reference heat capacity that acts as a conversion factor between the TEOS10 
# conservative temperature and potential enthalpy. See equation 3.3.3 (section 3.3, page 27)
# in the TEOS10 manual: http://www.teos-10.org/pubs/TEOS-10_Manual.pdf
const teos10_reference_heat_capacity = 3991.86795711963 # J kg⁻¹ K⁻¹

Base.eltype(::TEOS10SeawaterPolynomial{FT}) where FT = FT
Base.summary(::TEOS10SeawaterPolynomial{FT}) where FT = "TEOS10SeawaterPolynomial{$FT}"
with_float_type(FT::DataType, ::TEOS10SeawaterPolynomial) = TEOS10SeawaterPolynomial{FT}()

"""
    TEOS10SeawaterPolynomial(FT=Float64)

Return an object representing a 55-term polynomial approximation to the TEOS-10 standard equation
of state for seawater. See

> Roquet et al., "Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard", Ocean Modelling (2015).
"""
TEOS10SeawaterPolynomial(FT=Float64) = TEOS10SeawaterPolynomial{FT}()

const EOS₁₀ = BoussinesqEquationOfState{<:TEOS10SeawaterPolynomial}
const TEOS10EquationOfState{FT} = BoussinesqEquationOfState{<:TEOS10SeawaterPolynomial, FT} where FT

"""
    TEOS10EquationOfState(FT=Float64; reference_density=1020)

Return an `BoussinesqEquationOfState` with a `TEOS10SeawaterPolynomial` of float type `FT`
with `reference density = 1020 kg m⁻³`, the value used by 

> Roquet et al., "Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard", Ocean Modelling (2015).

when fitting polynomial coefficients to the full TEOS-10 standard equation of state.
See the discussion prior to equation 8 in Roquet et al. (2015).

Note that according to Roquet et al. (2015):

> "In a Boussinesq model, the choice of the ``ρ₀`` value is important, yet it varies significantly among OGCMs, as it is a matter of personal preference."
"""
TEOS10EquationOfState(FT=Float64; reference_density=1020) =
    BoussinesqEquationOfState(TEOS10SeawaterPolynomial{FT}(),
                              convert(FT, reference_density))

#####
##### Reference values chosen using TEOS-10 recommendation
#####

const Sₐᵤ = 40 * 35.16504 / 35
const Θᵤ  = 40.0
const Zᵤ  = 1e4
const ΔS  = 32.0

#####
##### Coordinate transformations from (Θ, Sᴬ, p) to (τ, s, ζ)
#####

@inline τ(Θ::FT) where FT = Θ / FT(Θᵤ)
@inline s(Sᴬ::FT) where FT = √((Sᴬ + FT(ΔS)) / FT(Sₐᵤ))
@inline ζ(Z::FT) where FT = - Z / FT(Zᵤ)

#####
##### Vertical reference profile of density
#####

const R₀₀ =  4.6494977072e+01
const R₀₁ = -5.2099962525e+00
const R₀₂ =  2.2601900708e-01
const R₀₃ =  6.4326772569e-02
const R₀₄ =  1.5616995503e-02
const R₀₅ = -1.7243708991e-03

@inline r₀(ζ) = (((((R₀₅ * ζ + R₀₄) * ζ + R₀₃) * ζ + R₀₂) * ζ + R₀₁) * ζ + R₀₀) * ζ

#####
##### Density anomaly fit
#####

const R₀₀₀ =  8.0189615746e+02
const R₁₀₀ =  8.6672408165e+02
const R₂₀₀ = -1.7864682637e+03
const R₃₀₀ =  2.0375295546e+03
const R₄₀₀ = -1.2849161071e+03
const R₅₀₀ =  4.3227585684e+02
const R₆₀₀ = -6.0579916612e+01
const R₀₁₀ =  2.6010145068e+01
const R₁₁₀ = -6.5281885265e+01
const R₂₁₀ =  8.1770425108e+01
const R₃₁₀ = -5.6888046321e+01
const R₄₁₀ =  1.7681814114e+01
const R₅₁₀ = -1.9193502195e+00
const R₀₂₀ = -3.7074170417e+01
const R₁₂₀ =  6.1548258127e+01
const R₂₂₀ = -6.0362551501e+01
const R₃₂₀ =  2.9130021253e+01
const R₄₂₀ = -5.4723692739e+00
const R₀₃₀ =  2.1661789529e+01
const R₁₃₀ = -3.3449108469e+01
const R₂₃₀ =  1.9717078466e+01
const R₃₃₀ = -3.1742946532e+00
const R₀₄₀ = -8.3627885467e+00
const R₁₄₀ =  1.1311538584e+01
const R₂₄₀ = -5.3563304045e+00
const R₀₅₀ =  5.4048723791e-01
const R₁₅₀ =  4.8169980163e-01
const R₀₆₀ = -1.9083568888e-01
const R₀₀₁ =  1.9681925209e+01
const R₁₀₁ = -4.2549998214e+01
const R₂₀₁ =  5.0774768218e+01
const R₃₀₁ = -3.0938076334e+01
const R₄₀₁ =  6.6051753097e+00
const R₀₁₁ = -1.3336301113e+01
const R₁₁₁ = -4.4870114575e+00
const R₂₁₁ =  5.0042598061e+00
const R₃₁₁ = -6.5399043664e-01
const R₀₂₁ =  6.7080479603e+00
const R₁₂₁ =  3.5063081279e+00
const R₂₂₁ = -1.8795372996e+00
const R₀₃₁ = -2.4649669534e+00
const R₁₃₁ = -5.5077101279e-01
const R₀₄₁ =  5.5927935970e-01
const R₀₀₂ =  2.0660924175e+00
const R₁₀₂ = -4.9527603989e+00
const R₂₀₂ =  2.5019633244e+00
const R₀₁₂ =  2.0564311499e+00
const R₁₁₂ = -2.1311365518e-01
const R₀₂₂ = -1.2419983026e+00
const R₀₀₃ = -2.3342758797e-02
const R₁₀₃ = -1.8507636718e-02
const R₀₁₃ =  3.7969820455e-01

@inline r′₃(τ::FT, s::FT) where FT = FT(R₀₁₃) * τ + FT(R₁₀₃) * s + FT(R₀₀₃)

@inline r′₂(τ::FT, s::FT) where FT = (FT(R₀₂₂) * τ + FT(R₁₁₂) * s + FT(R₀₁₂)) * τ + (FT(R₂₀₂) * s + FT(R₁₀₂)) * s + FT(R₀₀₂)

@inline r′₁(τ::FT, s::FT) where FT =
(((FT(R₀₄₁) * τ + FT(R₁₃₁) * s + FT(R₀₃₁)) * τ +
  (FT(R₂₂₁) * s + FT(R₁₂₁)) * s + FT(R₀₂₁)) * τ +
 ((FT(R₃₁₁) * s + FT(R₂₁₁)) * s + FT(R₁₁₁)) * s + FT(R₀₁₁)) * τ +
(((FT(R₄₀₁) * s + FT(R₃₀₁)) * s + FT(R₂₀₁)) * s + FT(R₁₀₁)) * s + FT(R₀₀₁)

@inline r′₀(τ::FT, s::FT) where FT =
(((((FT(R₀₆₀) * τ + FT(R₁₅₀) * s + FT(R₀₅₀)) * τ +
    (FT(R₂₄₀) * s + FT(R₁₄₀)) * s + FT(R₀₄₀)) * τ +
   ((FT(R₃₃₀) * s + FT(R₂₃₀)) * s + FT(R₁₃₀)) * s + FT(R₀₃₀)) * τ +
  (((FT(R₄₂₀) * s + FT(R₃₂₀)) * s + FT(R₂₂₀)) * s + FT(R₁₂₀)) * s + FT(R₀₂₀)) * τ +
 ((((FT(R₅₁₀) * s + FT(R₄₁₀)) * s + FT(R₃₁₀)) * s + FT(R₂₁₀)) * s + FT(R₁₁₀)) * s + FT(R₀₁₀)) * τ +
(((((FT(R₆₀₀) * s + FT(R₅₀₀)) * s + FT(R₄₀₀)) * s + FT(R₃₀₀)) * s + FT(R₂₀₀)) * s + FT(R₁₀₀)) * s + FT(R₀₀₀)

@inline r′(τ, s, ζ) = ((r′₃(τ, s) * ζ + r′₂(τ, s)) * ζ + r′₁(τ, s)) * ζ + r′₀(τ, s)

#####
##### Density perturbation
#####

"""
    ρ(Θ, Sᴬ, Z, ::BoussinesqEquationOfState{<:TEOS10EquationOfState})

Return the in-situ density of seawater with state `(Θ, Sᴬ, Z)` using the 55-term polynomial
approximation to TEOS-10 described in Roquet et al. (§3.1, 2014).

# Inputs
- `Θ`: conservative temperature (ITS-90) [°C]
- `Sᴬ`: absolute salinity [g/kg]
- `Z`: geopotential depth [m]

# Output
- `ρ`: in-situ density [kg/m³]

# References
- Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate polynomial expressions
  for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling.
"""
@inline ρ(Θ, Sᴬ, Z, ::EOS₁₀) = _ρ(τ(Θ), s(Sᴬ), ζ(Z))
@inline ρ′(Θ, Sᴬ, Z, eos::EOS₁₀) = r′(τ(Θ), s(Sᴬ), ζ(Z)) - eos.reference_density 

@inline _ρ(τ, s, ζ) = r₀(ζ) + r′(τ, s, ζ)

#####
##### Thermal sensitivity fit
#####

const α₀₀₀ = -6.5025362670e-01
const α₁₀₀ =  1.6320471316e+00
const α₂₀₀ = -2.0442606277e+00
const α₃₀₀ =  1.4222011580e+00
const α₄₀₀ = -4.4204535284e-01
const α₅₀₀ =  4.7983755487e-02
const α₀₁₀ =  1.8537085209e+00
const α₁₁₀ = -3.0774129064e+00
const α₂₁₀ =  3.0181275751e+00
const α₃₁₀ = -1.4565010626e+00
const α₄₁₀ =  2.7361846370e-01
const α₀₂₀ = -1.6246342147e+00
const α₁₂₀ =  2.5086831352e+00
const α₂₂₀ = -1.4787808849e+00
const α₃₂₀ =  2.3807209899e-01
const α₀₃₀ =  8.3627885467e-01
const α₁₃₀ = -1.1311538584e+00
const α₂₃₀ =  5.3563304045e-01
const α₀₄₀ = -6.7560904739e-02
const α₁₄₀ = -6.0212475204e-02
const α₀₅₀ =  2.8625353333e-02
const α₀₀₁ =  3.3340752782e-01
const α₁₀₁ =  1.1217528644e-01
const α₂₀₁ = -1.2510649515e-01
const α₃₀₁ =  1.6349760916e-02
const α₀₁₁ = -3.3540239802e-01
const α₁₁₁ = -1.7531540640e-01
const α₂₁₁ =  9.3976864981e-02
const α₀₂₁ =  1.8487252150e-01
const α₁₂₁ =  4.1307825959e-02
const α₀₃₁ = -5.5927935970e-02
const α₀₀₂ = -5.1410778748e-02
const α₁₀₂ =  5.3278413794e-03
const α₀₁₂ =  6.2099915132e-02
const α₀₀₃ = -9.4924551138e-03

"""
    thermal_sensitivity(Θ, Sᴬ, Z, ::TEOS10)

Return the Boussinesq thermal sensitivity coefficient ``-∂ρ/∂Θ`` [kg/m³/K] computed using
the 55-term polynomial approximation to TEOS-10 described in Roquet et al. (§3.1, 2014).

# Inputs
- `Θ`: conservative temperature (ITS-90) [°C]
- `Sᴬ`: absolute salinity [g/kg]
- `Z`: geopotential depth [m]

# Output
- `a`: Boussinesq thermal sensitivity coefficient ``-∂ρ/∂Θ`` [kg/m³/K]

# References
- Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate polynomial expressions
  for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling.
"""
@inline thermal_sensitivity(Θ, Sᴬ, Z, ::EOS₁₀) = _a(τ(Θ), s(Sᴬ), ζ(Z))

@inline _a(τ::FT, s::FT, ζ::FT) where FT =
     ((FT(α₀₀₃) * ζ  + FT(α₀₁₂)  * τ + FT(α₁₀₂)  * s + FT(α₀₀₂)) * ζ +
     ((FT(α₀₃₁) * τ  + FT(α₁₂₁)  * s + FT(α₀₂₁)) * τ +
      (FT(α₂₁₁) * s  + FT(α₁₁₁)) * s + FT(α₀₁₁)) * τ +
     ((FT(α₃₀₁) * s  + FT(α₂₀₁)) * s + FT(α₁₀₁)) * s + FT(α₀₀₁)) * ζ +
    ((((FT(α₀₅₀) * τ + FT(α₁₄₀)  * s + FT(α₀₄₀)) * τ +
       (FT(α₂₃₀) * s + FT(α₁₃₀)) * s + FT(α₀₃₀)) * τ +
      ((FT(α₃₂₀) * s + FT(α₂₂₀)) * s + FT(α₁₂₀)) * s + FT(α₀₂₀)) * τ +
     (((FT(α₄₁₀) * s + FT(α₃₁₀)) * s + FT(α₂₁₀)) * s + FT(α₁₁₀)) * s + FT(α₀₁₀)) * τ +
    ((((FT(α₅₀₀) * s + FT(α₄₀₀)) * s + FT(α₃₀₀)) * s + FT(α₂₀₀)) * s + FT(α₁₀₀)) * s + FT(α₀₀₀)

#####
##### Saline sensitivity
#####

const β₀₀₀ =  1.0783203594e+01
const β₁₀₀ = -4.4452095908e+01
const β₂₀₀ =  7.6048755820e+01
const β₃₀₀ = -6.3944280668e+01
const β₄₀₀ =  2.6890441098e+01
const β₅₀₀ = -4.5221697773e+00
const β₀₁₀ = -8.1219372432e-01
const β₁₁₀ =  2.0346663041e+00
const β₂₁₀ = -2.1232895170e+00
const β₃₁₀ =  8.7994140485e-01
const β₄₁₀ = -1.1939638360e-01
const β₀₂₀ =  7.6574242289e-01
const β₁₂₀ = -1.5019813020e+00
const β₂₂₀ =  1.0872489522e+00
const β₃₂₀ = -2.7233429080e-01
const β₀₃₀ = -4.1615152308e-01
const β₁₃₀ =  4.9061350869e-01
const β₂₃₀ = -1.1847737788e-01
const β₀₄₀ =  1.4073062708e-01
const β₁₄₀ = -1.3327978879e-01
const β₀₅₀ =  5.9929880134e-03
const β₀₀₁ = -5.2937873009e-01
const β₁₀₁ =  1.2634116779e+00
const β₂₀₁ = -1.1547328025e+00
const β₃₀₁ =  3.2870876279e-01
const β₀₁₁ = -5.5824407214e-02
const β₁₁₁ =  1.2451933313e-01
const β₂₁₁ = -2.4409539932e-02
const β₀₂₁ =  4.3623149752e-02
const β₁₂₁ = -4.6767901790e-02
const β₀₃₁ = -6.8523260060e-03
const β₀₀₂ = -6.1618945251e-02
const β₁₀₂ =  6.2255521644e-02
const β₀₁₂ = -2.6514181169e-03
const β₀₀₃ = -2.3025968587e-04

"""
    haline_sensitivity(Θ, Sᴬ, Z, ::TEOS10)

Return the Boussinesq haline sensitivity coefficient ``∂ρ/∂Sᴬ`` [kg/m³/(g/kg)] computed using
the 55-term polynomial approximation to TEOS-10 described in Roquet et al. (§3.1, 2014).

# Inputs
- `Θ`: conservative temperature (ITS-90) [°C]
- `Sᴬ`: absolute salinity [g/kg]
- `Z`: geopotential depth [m]

# Output
- `b`: Boussinesq haline sensitivity coefficient ``∂ρ/∂Sᴬ`` [kg/m³/(g/kg)]

# References
- Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate polynomial expressions
  for the density and specific volume of seawater using the TEOS-10 standard. Ocean Modelling.
"""
@inline haline_sensitivity(Θ, Sᴬ, Z, ::EOS₁₀) = _b(τ(Θ), s(Sᴬ), ζ(Z)) / s(Sᴬ)

@inline _b(τ::FT, s::FT, ζ::FT) where FT =
      ((FT(β₀₀₃) * ζ + FT(β₀₁₂)  * τ + FT(β₁₀₂) * s + FT(β₀₀₂)) * ζ +
      ((FT(β₀₃₁) * τ + FT(β₁₂₁)  * s + FT(β₀₂₁)) * τ +
       (FT(β₂₁₁) * s + FT(β₁₁₁)) * s + FT(β₀₁₁)) * τ +
      ((FT(β₃₀₁) * s + FT(β₂₀₁)) * s + FT(β₁₀₁)) * s + FT(β₀₀₁)) * ζ +
    ((((FT(β₀₅₀) * τ + FT(β₁₄₀)  * s + FT(β₀₄₀)) * τ +
       (FT(β₂₃₀) * s + FT(β₁₃₀)) * s + FT(β₀₃₀)) * τ +
      ((FT(β₃₂₀) * s + FT(β₂₂₀)) * s + FT(β₁₂₀)) * s + FT(β₀₂₀)) * τ +
     (((FT(β₄₁₀) * s + FT(β₃₁₀)) * s + FT(β₂₁₀)) * s + FT(β₁₁₀)) * s + FT(β₀₁₀)) * τ +
    ((((FT(β₅₀₀) * s + FT(β₄₀₀)) * s + FT(β₃₀₀)) * s + FT(β₂₀₀)) * s + FT(β₁₀₀)) * s + FT(β₀₀₀)

end # module
