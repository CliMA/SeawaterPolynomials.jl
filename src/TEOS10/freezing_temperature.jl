#####
##### TEOS-10 freezing temperature of seawater, translated from `gsw_ct_freezing_poly` and
##### `gsw_ct_freezing_first_derivatives_poly` of `gsw_oceanographic_toolbox.c` of https://github.com/TEOS-10/GSW-C.
#####
##### The fit is written in the reduced variables  Sᵣ = Sᴬ / 100,  x = √Sᵣ,  pᵣ = p / 10⁴  (McDougall et al., 2014).
#####

#####
##### Coefficients (gsw_internal_const.h, GSW_FREEZING_POLY_COEFFICIENTS)
#####

const cᶠ₀  =  0.017947064327968736
const cᶠ₁  = -6.076099099929818
const cᶠ₂  =  4.883198653547851
const cᶠ₃  = -11.88081601230542
const cᶠ₄  =  13.34658511480257
const cᶠ₅  = -8.722761043208607
const cᶠ₆  =  2.082038908808201
const cᶠ₇  = -7.389420998107497
const cᶠ₈  = -2.110913185058476
const cᶠ₉  =  0.2295491578006229
const cᶠ₁₀ = -0.9891538123307282
const cᶠ₁₁ = -0.08987150128406496
const cᶠ₁₂ =  0.3831132432071728
const cᶠ₁₃ =  1.054318231187074
const cᶠ₁₄ =  1.065556599652796
const cᶠ₁₅ = -0.7997496801694032
const cᶠ₁₆ =  0.3850133554097069
const cᶠ₁₇ = -2.078616693017569
const cᶠ₁₈ =  0.8756340772729538
const cᶠ₁₉ = -2.079022768390933
const cᶠ₂₀ =  1.596435439942262
const cᶠ₂₁ =  0.1338002171109174
const cᶠ₂₂ =  1.242891021876471

# Dissolved-air correction, (a₀ - aᶠ Sᴬ) (1 + bᶠ (1 - Sᴬ / Sₒ)) K per unit saturation fraction
const a₀ = 2.4e-3
const aᶠ = 1.4289763856964e-5
const bᶠ = 0.057000649899720

"""
    Θᶠ(Sᴬ, p, saturation_fraction = 1)

Return the TEOS-10 conservative temperature ``Θᶠ`` at which seawater of absolute salinity ``Sᴬ`` freezes at sea
pressure ``p``, from the polynomial fit of McDougall et al. (2014), which is accurate to within
``-5 × 10⁻⁴`` K and ``6 × 10⁻⁴`` K of the exact freezing temperature. Dissolved air lowers the freezing
temperature by up to about ``2.4 × 10⁻³`` K, scaled by `saturation_fraction`.
Direct translation of `gsw_ct_freezing_poly` of https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴬ` : absolute salinity                                           [g/kg]
- `p`  : sea pressure (absolute pressure - 10.1325 dbar)             [dbar]
- `saturation_fraction` : saturation fraction of dissolved air, 0 (air-free) to 1 (air-saturated)

# Output
- `Θᶠ` : conservative temperature at freezing                        [°C]

# References
- McDougall, T. J., P. M. Barker, R. Feistel and B. K. Galton-Fenzi, 2014: Melting of ice and sea ice into
  seawater and frazil ice formation. Journal of Physical Oceanography, 44, 1751–1775.
"""
@inline function Θᶠ(Sᴬ, p, saturation_fraction = 1)
    Sᴬ, p, saturation_fraction = map(float, promote(Sᴬ, p, saturation_fraction))
    FT = typeof(Sᴬ)
    Sᵣ = Sᴬ / 100
    x  = sqrt(Sᵣ)
    pᵣ = p * FT(rp)

    Θᵖ = FT(cᶠ₀) + Sᵣ * (FT(cᶠ₁) + x * (FT(cᶠ₂) + x * (FT(cᶠ₃) + x * (FT(cᶠ₄) + x * (FT(cᶠ₅) + FT(cᶠ₆) * x))))) +
         pᵣ * (FT(cᶠ₇) + pᵣ * (FT(cᶠ₈) + FT(cᶠ₉) * pᵣ)) +
         Sᵣ * pᵣ * (FT(cᶠ₁₀) + pᵣ * (FT(cᶠ₁₂) + pᵣ * (FT(cᶠ₁₅) + FT(cᶠ₂₁) * Sᵣ)) +
                    Sᵣ * (FT(cᶠ₁₃) + FT(cᶠ₁₇) * pᵣ + FT(cᶠ₁₉) * Sᵣ) +
                    x * (FT(cᶠ₁₁) + pᵣ * (FT(cᶠ₁₄) + FT(cᶠ₁₈) * pᵣ) + Sᵣ * (FT(cᶠ₁₆) + FT(cᶠ₂₀) * pᵣ + FT(cᶠ₂₂) * Sᵣ)))

    air_correction = saturation_fraction * (FT(a₀) - FT(aᶠ) * Sᴬ) * (1 + FT(bᶠ) * (1 - Sᴬ / FT(Sₒ)))

    return Θᵖ - air_correction
end

"""
    ∂Θᶠ∂Sᴬ(Sᴬ, p, saturation_fraction = 1)

Return ``∂Θᶠ/∂Sᴬ``, the derivative of [`Θᶠ`](@ref) with respect to absolute
salinity at fixed sea pressure, in K/(g/kg). This is the exact derivative of the polynomial fit. It matches
`gsw_ct_freezing_first_derivatives_poly` of https://github.com/TEOS-10/GSW-C for air-free seawater, and
differs by the sign of one term of the dissolved-air correction (about ``3 × 10⁻⁶`` K kg/g) otherwise.

# Inputs
- `Sᴬ` : absolute salinity                                           [g/kg]
- `p`  : sea pressure (absolute pressure - 10.1325 dbar)             [dbar]
- `saturation_fraction` : saturation fraction of dissolved air, 0 (air-free) to 1 (air-saturated)
"""
@inline function ∂Θᶠ∂Sᴬ(Sᴬ, p, saturation_fraction = 1)
    Sᴬ, p, saturation_fraction = map(float, promote(Sᴬ, p, saturation_fraction))
    FT = typeof(Sᴬ)
    x  = sqrt(Sᴬ / 100)
    pᵣ = p * FT(rp)

    ∂Θᶠ∂Sᵣ = FT(cᶠ₁) + x * (3 * FT(cᶠ₂) / 2 + x * (2 * FT(cᶠ₃) + x * (5 * FT(cᶠ₄) / 2 + x * (3 * FT(cᶠ₅) + 7 * FT(cᶠ₆) / 2 * x)))) +
             pᵣ * (FT(cᶠ₁₀) + x * (3 * FT(cᶠ₁₁) / 2 + x * (2 * FT(cᶠ₁₃) + x * (5 * FT(cᶠ₁₆) / 2 + x * (3 * FT(cᶠ₁₉) + 7 * FT(cᶠ₂₂) / 2 * x)))) +
             pᵣ * (FT(cᶠ₁₂) + x * (3 * FT(cᶠ₁₄) / 2 + x * (2 * FT(cᶠ₁₇) + 5 * FT(cᶠ₂₀) / 2 * x)) +
             pᵣ * (FT(cᶠ₁₅) + x * (3 * FT(cᶠ₁₈) / 2 + 2 * FT(cᶠ₂₁) * x))))

    ∂air∂Sᴬ = - saturation_fraction * (FT(aᶠ) * (1 + FT(bᶠ) * (1 - Sᴬ / FT(Sₒ))) +
                                       FT(bᶠ) * (FT(a₀) - FT(aᶠ) * Sᴬ) / FT(Sₒ))

    return ∂Θᶠ∂Sᵣ / 100 - ∂air∂Sᴬ
end

"""
    ∂Θᶠ∂p(Sᴬ, p)

Return ``∂Θᶠ/∂p``, the derivative of [`Θᶠ`](@ref) with respect to sea pressure
at fixed absolute salinity, in K/dbar. It does not depend on the saturation fraction.
Translation of `gsw_ct_freezing_first_derivatives_poly` of https://github.com/TEOS-10/GSW-C,
which returns the same derivative in K/Pa.

# Inputs
- `Sᴬ` : absolute salinity                                           [g/kg]
- `p`  : sea pressure (absolute pressure - 10.1325 dbar)             [dbar]
"""
@inline function ∂Θᶠ∂p(Sᴬ, p)
    Sᴬ, p = map(float, promote(Sᴬ, p))
    FT = typeof(Sᴬ)
    Sᵣ = Sᴬ / 100
    x  = sqrt(Sᵣ)
    pᵣ = p * FT(rp)

    ∂Θᶠ∂pᵣ = FT(cᶠ₇) + Sᵣ * (FT(cᶠ₁₀) + x * (FT(cᶠ₁₁) + x * (FT(cᶠ₁₃) + x * (FT(cᶠ₁₆) + x * (FT(cᶠ₁₉) + FT(cᶠ₂₂) * x))))) +
             pᵣ * (2 * FT(cᶠ₈) + Sᵣ * (2 * FT(cᶠ₁₂) + x * (2 * FT(cᶠ₁₄) + x * (2 * FT(cᶠ₁₇) + 2 * FT(cᶠ₂₀) * x))) +
             pᵣ * (3 * FT(cᶠ₉) + Sᵣ * (3 * FT(cᶠ₁₅) + x * (3 * FT(cᶠ₁₈) + 3 * FT(cᶠ₂₁) * x))))

    return ∂Θᶠ∂pᵣ * FT(rp)
end
