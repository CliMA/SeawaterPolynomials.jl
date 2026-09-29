#####
##### TEOS-10 temperature conversions (Θ ↔ θᴾ ↔ T), translated from `gsw_oceanographic_toolbox.c` of
##### https://github.com/TEOS-10/GSW-C.
#####
##### Notation (TEOS-10 manual, §A.1):
#####   Sᴬ : absolute salinity                      [g/kg]
#####   T  : in-situ temperature, ITS-90            [°C]
#####   θᴾ : potential temperature, p_ref = 0 dbar  [°C]
#####   Θ  : conservative temperature               [°C]
#####   p  : sea pressure (gauge)                   [dbar]
#####
##### The Gibbs-function polynomials are written in the non-dimensional TEOS-10 variables
#####   x = √(rS Sᴬ),  y = rT T,  z = rp p,
##### and their coefficients are named after the monomial they multiply: ηᵢⱼₖ multiplies xⁱ yʲ zᵏ.
#####

#####
##### Constants (gsw_internal_const.h, GSW_TEOS10_CONSTANTS macro)
#####

const T₀  = 273.15          # gsw_t0   [K]      Celsius zero point
const Sₒ  = 35.16504        # gsw_sso  [g/kg]   standard ocean salinity
const uₚₛ = Sₒ / 35         # gsw_ups  [g/kg]   PSS-78 to absolute-salinity scale factor
const rS  = 1 / (40 * uₚₛ)  # gsw_sfac [kg/g]   x² = rS Sᴬ
const rT  = 0.025           #          [1/°C]   y  = rT T
const rT² = 0.000625        #          [1/°C²]  ∂²/∂T² = rT² ∂²/∂y²
const rp  = 1e-4            #          [1/dbar] z  = rp p

# cᵖ ≈ cᵖ⁰ / (1 - heat_capacity_salinity_coefficient (1 - Sᴬ / Sₒ)), used for the first ∂θᴾ/∂η estimate
const heat_capacity_salinity_coefficient = 0.05

#####
##### Specific entropy (gsw_entropy_part, gsw_entropy_part_zerop):
#####   η(Sᴬ, T, p) = - rT Σ ηᵢⱼₖ xⁱ yʲ zᵏ, excluding the terms that do not depend on T
#####

# Water part (x⁰)
const η₀₀₁ = -270.983805184062
const η₀₀₂ =  776.153611613101
const η₀₀₃ = -196.51255088122
const η₀₀₄ =  28.9796526294175
const η₀₀₅ = -2.13290083518327
const η₀₁₀ = -24715.571866078
const η₀₁₁ =  2910.0729080936
const η₀₁₂ = -1513.116771538718
const η₀₁₃ =  546.959324647056
const η₀₁₄ = -111.1208127634436
const η₀₁₅ =  8.68841343834394
const η₀₂₀ =  2210.2236124548363
const η₀₂₁ = -2017.52334943521
const η₀₂₂ =  1498.081172457456
const η₀₂₃ = -718.6359919632359
const η₀₂₄ =  146.4037555781616
const η₀₂₅ = -4.9892131862671505
const η₀₃₀ = -592.743745734632
const η₀₃₁ =  1591.873781627888
const η₀₃₂ = -1207.261522487504
const η₀₃₃ =  608.785486935364
const η₀₃₄ = -105.4993508931208
const η₀₄₀ =  290.12956292128547
const η₀₄₁ = -973.091553087975
const η₀₄₂ =  602.603274510125
const η₀₄₃ = -276.361526170076
const η₀₄₄ =  32.40953340386105
const η₀₅₀ = -113.90630790850321
const η₀₅₁ =  381.06836198507096
const η₀₅₂ = -133.7383902842754
const η₀₅₃ =  49.023632509086724
const η₀₆₀ =  21.35571525415769
const η₀₆₁ = -67.41756835751434

# Saline part (x², x³, x⁴)
const η₂₀₁ =  729.116529735046
const η₂₀₂ = -343.956902961561
const η₂₀₃ =  124.687671116248
const η₂₀₄ = -31.656964386073
const η₂₀₅ =  7.04658803315449
const η₂₁₀ =  1760.062705994408
const η₂₁₁ = -1721.528607567954
const η₂₁₂ =  674.819060538734
const η₂₁₃ = -356.629112415276
const η₂₁₄ =  88.4080716616
const η₂₁₅ = -15.84003094423364
const η₂₂₀ = -675.802947790203
const η₂₂₁ =  2082.7344423998043
const η₂₂₂ = -614.668925894709
const η₂₂₃ =  340.685093521782
const η₂₂₄ = -33.3848202979239
const η₂₃₀ =  365.7041791005036
const η₂₃₁ = -1190.914967948748
const η₂₃₂ =  298.904564555024
const η₂₃₃ = -145.9491676006352
const η₂₄₀ = -108.30162043765552
const η₂₅₀ =  12.78101825083098
const η₃₀₁ = -175.292041186547
const η₃₀₂ =  83.1923927801819
const η₃₀₃ = -29.483064349429
const η₃₁₀ = -86.1329351956084
const η₃₁₁ =  766.116132004952
const η₃₁₂ = -108.3834525034224
const η₃₁₃ =  51.2796974779828
const η₃₂₀ = -30.0682112585625
const η₃₂₁ = -1380.9597954037708
const η₃₃₀ =  3.50240264723578
const η₃₃₁ =  938.26075044542
const η₄₀₁ = -22.6683558512829
const η₄₁₀ = -137.1145018408982
const η₄₂₀ =  148.10030845687618
const η₄₃₀ = -68.5590309679152
const η₄₄₀ =  12.4848504784754

@inline function entropy_part(Sᴬ::FT, T::FT, p::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = T * FT(rT)
    z  = p * FT(rp)

    # Water part as a polynomial in y whose yʲ-coefficient ηᵂⱼ is a polynomial in z
    ηᵂ₀ =                  z*(FT(η₀₀₁) + z*(FT(η₀₀₂) + z*(FT(η₀₀₃) + z*(FT(η₀₀₄) + z*FT(η₀₀₅)))))
    ηᵂ₁ = FT(η₀₁₀) + z*(FT(η₀₁₁) + z*(FT(η₀₁₂) + z*(FT(η₀₁₃) + z*(FT(η₀₁₄) + z*FT(η₀₁₅)))))
    ηᵂ₂ = FT(η₀₂₀) + z*(FT(η₀₂₁) + z*(FT(η₀₂₂) + z*(FT(η₀₂₃) + z*(FT(η₀₂₄) + z*FT(η₀₂₅)))))
    ηᵂ₃ = FT(η₀₃₀) + z*(FT(η₀₃₁) + z*(FT(η₀₃₂) + z*(FT(η₀₃₃) + z*FT(η₀₃₄))))
    ηᵂ₄ = FT(η₀₄₀) + z*(FT(η₀₄₁) + z*(FT(η₀₄₂) + z*(FT(η₀₄₃) + z*FT(η₀₄₄))))
    ηᵂ₅ = FT(η₀₅₀) + z*(FT(η₀₅₁) + z*(FT(η₀₅₂) + z*FT(η₀₅₃)))
    ηᵂ₆ = FT(η₀₆₀) + z*FT(η₀₆₁)

    ηᵂ = ηᵂ₀ + y*(ηᵂ₁ + y*(ηᵂ₂ + y*(ηᵂ₃ + y*(ηᵂ₄ + y*(ηᵂ₅ + y*ηᵂ₆)))))

    # Saline part, factored as x² (ηᶻ + x ηˣ + y ηʸ)
    ηᶻ = z*(FT(η₂₀₁) + z*(FT(η₂₀₂) + z*(FT(η₂₀₃) + z*(FT(η₂₀₄) + z*FT(η₂₀₅)))))

    ηˣ = x*(y*(FT(η₄₁₀) + y*(FT(η₄₂₀) + y*(FT(η₄₃₀) + y*FT(η₄₄₀)))) + z*FT(η₄₀₁)) +
         z*(FT(η₃₀₁) + z*(FT(η₃₀₂) + z*FT(η₃₀₃))) +
         y*(FT(η₃₁₀) + z*(FT(η₃₁₁) + z*(FT(η₃₁₂) + z*FT(η₃₁₃))) +
         y*(FT(η₃₂₀) + z*FT(η₃₂₁) + y*(FT(η₃₃₀) + z*FT(η₃₃₁))))

    ηʸ = FT(η₂₁₀) + y*(FT(η₂₂₀) +
                    y*(FT(η₂₃₀) + y*(FT(η₂₄₀) + y*FT(η₂₅₀)) +
                       z*(FT(η₂₃₁) + z*(FT(η₂₃₂) + z*FT(η₂₃₃)))) +
                    z*(FT(η₂₂₁) + z*(FT(η₂₂₂) + z*(FT(η₂₂₃) + z*FT(η₂₂₄))))) +
                    z*(FT(η₂₁₁) + z*(FT(η₂₁₂) + z*(FT(η₂₁₃) + z*(FT(η₂₁₄) + z*FT(η₂₁₅)))))

    ηˢ = x² * (ηᶻ + x*ηˣ + y*ηʸ)

    return -(ηᵂ + ηˢ) * FT(rT)
end

@inline function zero_pressure_entropy_part(Sᴬ::FT, θᴾ::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = θᴾ * FT(rT)

    ηᵂ = y*(FT(η₀₁₀) + y*(FT(η₀₂₀) + y*(FT(η₀₃₀) + y*(FT(η₀₄₀) + y*(FT(η₀₅₀) + y*FT(η₀₆₀))))))

    ηˣ = x*(y*(FT(η₄₁₀) + y*(FT(η₄₂₀) + y*(FT(η₄₃₀) + FT(η₄₄₀)*y)))) +
         y*(FT(η₃₁₀) + y*(FT(η₃₂₀) + y*FT(η₃₃₀)))

    ηʸ = FT(η₂₁₀) + y*(FT(η₂₂₀) + y*(FT(η₂₃₀) + y*(FT(η₂₄₀) + y*FT(η₂₅₀))))

    ηˢ = x² * (x*ηˣ + y*ηʸ)

    return -(ηᵂ + ηˢ) * FT(rT)
end

#####
##### Gibbs-function curvature at the sea surface (gsw_gibbs_pt0_pt0):
#####   ∂²g/∂T²(Sᴬ, θᴾ, 0) = rT² Σ gᵀᵀᵢⱼ xⁱ yʲ
#####

const gᵀᵀ₀₀ = -24715.571866078
const gᵀᵀ₀₁ =  4420.4472249096725
const gᵀᵀ₀₂ = -1778.231237203896
const gᵀᵀ₀₃ =  1160.5182516851419
const gᵀᵀ₀₄ = -569.531539542516
const gᵀᵀ₀₅ =  128.13429152494615
const gᵀᵀ₂₀ =  1760.062705994408
const gᵀᵀ₂₁ = -1351.605895580406
const gᵀᵀ₂₂ =  1097.1125373015109
const gᵀᵀ₂₃ = -433.20648175062206
const gᵀᵀ₂₄ =  63.905091254154904
const gᵀᵀ₃₀ = -86.1329351956084
const gᵀᵀ₃₁ = -60.136422517125
const gᵀᵀ₃₂ =  10.50720794170734
const gᵀᵀ₄₀ = -137.1145018408982
const gᵀᵀ₄₁ =  296.20061691375236
const gᵀᵀ₄₂ = -205.67709290374563
const gᵀᵀ₄₃ =  49.9394019139016

@inline function zero_pressure_gibbs_curvature(Sᴬ::FT, θᴾ::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = θᴾ * FT(rT)

    gᵀᵀᵂ = FT(gᵀᵀ₀₀) + y*(FT(gᵀᵀ₀₁) + y*(FT(gᵀᵀ₀₂) + y*(FT(gᵀᵀ₀₃) + y*(FT(gᵀᵀ₀₄) + y*FT(gᵀᵀ₀₅)))))

    gᵀᵀˣ = FT(gᵀᵀ₃₀) + x*(FT(gᵀᵀ₄₀) + y*(FT(gᵀᵀ₄₁) + y*(FT(gᵀᵀ₄₂) + y*FT(gᵀᵀ₄₃)))) + y*(FT(gᵀᵀ₃₁) + y*FT(gᵀᵀ₃₂))
    gᵀᵀʸ = FT(gᵀᵀ₂₁) + y*(FT(gᵀᵀ₂₂) + y*(FT(gᵀᵀ₂₃) + y*FT(gᵀᵀ₂₄)))
    gᵀᵀˢ = x² * (FT(gᵀᵀ₂₀) + x*gᵀᵀˣ + y*gᵀᵀʸ)

    return (gᵀᵀᵂ + gᵀᵀˢ) * FT(rT²)
end

#####
##### Potential enthalpy (gsw_ct_from_pt):
#####   h⁰(Sᴬ, θᴾ) = Σ hᵢⱼ xⁱ yʲ
#####

const h₀₀ =  61.01362420681071
const h₀₁ =  168776.46138048015
const h₀₂ = -2735.2785605119625
const h₀₃ =  2574.2164453821433
const h₀₄ = -1536.6644434977543
const h₀₅ =  545.7340497931629
const h₀₆ = -50.91091728474331
const h₀₇ = -18.30489878927802
const h₂₀ =  268.5520265845071
const h₂₁ = -12019.028203559312
const h₂₂ =  3734.858026725145
const h₂₃ = -2046.7671145057618
const h₂₄ =  465.28655623826234
const h₂₅ = -0.6370820302376359
const h₂₆ = -10.650848542359153
const h₃₀ =  937.2099110620707
const h₃₁ =  588.1802812170108
const h₃₂ =  248.39476522971285
const h₃₃ = -3.871557904936333
const h₃₄ = -2.6268019854268356
const h₄₀ = -1687.914374187449
const h₄₁ =  936.3206544460336
const h₄₂ = -942.7827304544439
const h₄₃ =  369.4389437509002
const h₄₄ = -33.83664947895248
const h₄₅ = -9.987880382780322
const h₅₀ =  246.9598888781377
const h₆₀ =  123.59576582457964
const h₇₀ = -48.5891069025409

"""
    Θ_from_θᴾ(Sᴬ, θᴾ)

Return the TEOS-10 conservative temperature ``Θ`` from absolute salinity ``Sᴬ`` and potential temperature
``θᴾ`` referenced to ``p = 0`` dbar, computed as ``Θ = h⁰(Sᴬ, θᴾ) / cᵖ⁰`` where ``h⁰`` is the potential-enthalpy
polynomial and ``cᵖ⁰ = 3991.867957...`` J/kg/K is the TEOS-10 reference heat capacity.
Direct translation of `gsw_ct_from_pt` of https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴬ`: absolute salinity                              [g/kg]
- `θᴾ`: potential temperature, ITS-90, p_ref = 0 dbar  [°C]

# Output
- `Θ` : conservative temperature                       [°C]

# References
- IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater – 2010: Calculation and
  use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
  UNESCO. http://www.teos-10.org/pubs/TEOS-10_Manual.pdf
"""
@inline function Θ_from_θᴾ(Sᴬ, θᴾ)
    Sᴬ, θᴾ = map(float, promote(Sᴬ, θᴾ))
    FT = typeof(Sᴬ)
    cᵖ⁰ = teos10_reference_heat_capacity
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = θᴾ * FT(rT)

    hᵂ = FT(h₀₀) + y*(FT(h₀₁) + y*(FT(h₀₂) + y*(FT(h₀₃) + y*(FT(h₀₄) + y*(FT(h₀₅) + y*(FT(h₀₆) + y*FT(h₀₇)))))))

    # Saline part as a polynomial in x whose xⁱ-coefficient hᵢ is a polynomial in y
    h₂ = FT(h₂₀) + y*(FT(h₂₁) + y*(FT(h₂₂) + y*(FT(h₂₃) + y*(FT(h₂₄) + y*(FT(h₂₅) + y*FT(h₂₆))))))
    h₃ = FT(h₃₀) + y*(FT(h₃₁) + y*(FT(h₃₂) + y*(FT(h₃₃) + y*FT(h₃₄))))
    h₄ = FT(h₄₀) + y*(FT(h₄₁) + y*(FT(h₄₂) + y*(FT(h₄₃) + y*(FT(h₄₄) + y*FT(h₄₅)))))

    hˢ = x² * (h₂ + x*(h₃ + x*(h₄ + x*(FT(h₅₀) + x*(FT(h₆₀) + x*FT(h₇₀))))))

    return (hᵂ + hˢ) / FT(cᵖ⁰)
end

#####
##### Potential temperature from in-situ temperature (gsw_pt0_from_t), with a first guess
#####   θᴾ ≈ T + Σ θᴾᵢⱼₖ s₁ⁱ Tʲ pᵏ,  s₁ = Sᴬ / uₚₛ
#####

const θᴾ₀₀₁ =  8.65483913395442e-6
const θᴾ₁₀₁ = -1.41636299744881e-6
const θᴾ₀₀₂ = -7.38286467135737e-9
const θᴾ₀₁₁ = -8.38241357039698e-6
const θᴾ₁₁₁ =  2.83933368585534e-8
const θᴾ₀₂₁ =  1.77803965218656e-8
const θᴾ₀₁₂ =  1.71155619208233e-10

"""
    θᴾ_from_T(Sᴬ, T, p)

Return the TEOS-10 potential temperature ``θᴾ`` referenced to ``p = 0`` dbar from absolute salinity,
in-situ temperature and sea pressure. ``θᴾ`` is the root of ``η(Sᴬ, θᴾ, 0) = η(Sᴬ, T, p)``, found with a
polynomial first guess followed by two modified Newton–Raphson iterations (McDougall and Wotherspoon, 2014).
Direct translation of `gsw_pt0_from_t` of https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴬ`: absolute salinity                                       [g/kg]
- `T` : in-situ temperature, ITS-90                             [°C]
- `p` : sea pressure (gauge: absolute pressure − 10.1325 dbar)  [dbar]

# Output
- `θᴾ`: potential temperature, p_ref = 0 dbar                   [°C]

# References
- IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater – 2010.
  http://www.teos-10.org/pubs/TEOS-10_Manual.pdf
- McDougall, T. J. and S. J. Wotherspoon, 2014: A simple modification of Newton's method to achieve
  convergence of order 1 + √2. Applied Mathematics Letters, 29, 20–25.
"""
@inline function θᴾ_from_T(Sᴬ, T, p)
    Sᴬ, T, p = map(float, promote(Sᴬ, T, p))
    FT = typeof(Sᴬ)
    cᵖ⁰ = teos10_reference_heat_capacity
    s₁ = Sᴬ / FT(uₚₛ)

    θᴾ = T + p*(FT(θᴾ₀₀₁) + s₁*FT(θᴾ₁₀₁) + p*FT(θᴾ₀₀₂) +
             T*(FT(θᴾ₀₁₁) + s₁*FT(θᴾ₁₁₁) + T*FT(θᴾ₀₂₁) + p*FT(θᴾ₀₁₂)))

    dθᴾdη = (FT(T₀) + θᴾ) * (1 - FT(heat_capacity_salinity_coefficient) * (1 - Sᴬ / FT(Sₒ))) / FT(cᵖ⁰)

    η = entropy_part(Sᴬ, T, p)

    # Modified Newton–Raphson on η(Sᴬ, θᴾ, 0) = η(Sᴬ, T, p), with ∂θᴾ/∂η re-evaluated at the midpoint of each step
    for _ in 1:2
        θᴾⁿ   = θᴾ
        Δη    = zero_pressure_entropy_part(Sᴬ, θᴾⁿ) - η
        θᴾ    = θᴾⁿ - Δη * dθᴾdη
        θᴾᵐ   = (θᴾ + θᴾⁿ) / 2
        dθᴾdη = - 1 / zero_pressure_gibbs_curvature(Sᴬ, θᴾᵐ)
        θᴾ    = θᴾⁿ - Δη * dθᴾdη
    end

    return θᴾ
end

#####
##### Conservative temperature from in-situ temperature (gsw_ct_from_t)
#####

"""
    Θ_from_T(Sᴬ, T, p)

Return the TEOS-10 conservative temperature ``Θ`` from absolute salinity, in-situ temperature and sea
pressure, computed as `Θ_from_θᴾ(Sᴬ, θᴾ_from_T(Sᴬ, T, p))`. Direct translation of `gsw_ct_from_t` of
https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴬ`: absolute salinity            [g/kg]
- `T` : in-situ temperature, ITS-90  [°C]
- `p` : sea pressure                 [dbar]

# Output
- `Θ` : conservative temperature     [°C]

# References
- IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater – 2010.
  http://www.teos-10.org/pubs/TEOS-10_Manual.pdf
"""
@inline Θ_from_T(Sᴬ, T, p) = Θ_from_θᴾ(Sᴬ, θᴾ_from_T(Sᴬ, T, p))

#####
##### Potential temperature from conservative temperature (gsw_pt_from_ct), with a rational first guess
#####   θᴾ ≈ Σ νᵢⱼ s₁ⁱ Θʲ / Σ δᵢⱼ s₁ⁱ Θʲ,  s₁ = Sᴬ / uₚₛ
#####

const ν₀₀ = -1.446013646344788e-2
const ν₁₀ = -3.305308995852924e-3
const ν₂₀ =  1.062415929128982e-4
const ν₀₁ =  9.477566673794488e-1
const ν₁₁ =  2.166591947736613e-3
const ν₀₂ =  3.828842955039902e-3
const δ₀₀ =  1.000000000000000e0
const δ₁₀ =  6.506097115635800e-4
const δ₀₁ =  3.830289486850898e-3
const δ₀₂ =  1.247811760368034e-6

"""
    θᴾ_from_Θ(Sᴬ, Θ)

Return the TEOS-10 potential temperature ``θᴾ`` referenced to ``p = 0`` dbar from absolute salinity and
conservative temperature. ``θᴾ`` is the root of ``Θ_from_θᴾ(Sᴬ, θᴾ) = Θ``, found with a rational-polynomial
first guess followed by one and a half modified Newton–Raphson iterations (McDougall and Wotherspoon, 2014).
Direct translation of `gsw_pt_from_ct` of https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴬ`: absolute salinity                          [g/kg]
- `Θ` : conservative temperature                   [°C]

# Output
- `θᴾ`: potential temperature, p_ref = 0 dbar      [°C]

# References
- IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater – 2010.
  http://www.teos-10.org/pubs/TEOS-10_Manual.pdf
- McDougall, T. J. and S. J. Wotherspoon, 2014: A simple modification of Newton's method to achieve
  convergence of order 1 + √2. Applied Mathematics Letters, 29, 20–25.
"""
@inline function θᴾ_from_Θ(Sᴬ, Θ)
    Sᴬ, Θ = map(float, promote(Sᴬ, Θ))
    FT = typeof(Sᴬ)
    cᵖ⁰ = teos10_reference_heat_capacity
    s₁ = Sᴬ / FT(uₚₛ)

    ν₀₂Θ = FT(ν₀₂) * Θ
    δ₀₂Θ = FT(δ₀₂) * Θ

    χ  = FT(ν₀₁) + FT(ν₁₁) * s₁ + ν₀₂Θ
    ν  = FT(ν₀₀) + s₁ * (FT(ν₁₀) + FT(ν₂₀) * s₁) + Θ * χ
    δ  = FT(δ₀₀) + FT(δ₁₀) * s₁ + Θ * (FT(δ₀₁) + δ₀₂Θ)
    θᴾ = ν / δ

    dνdΘ  = χ + ν₀₂Θ
    dδdΘ  = FT(δ₀₁) + δ₀₂Θ + δ₀₂Θ
    dθᴾdΘ = (dνdΘ - dδdΘ * θᴾ) / δ

    # Modified Newton–Raphson on Θ_from_θᴾ(Sᴬ, θᴾ) = Θ, with ∂θᴾ/∂Θ re-evaluated at the midpoint of the first step
    ΔΘ  =Θ_from_θᴾ(Sᴬ, θᴾ) - Θ
    θᴾⁿ = θᴾ
    θᴾ  = θᴾⁿ - ΔΘ * dθᴾdΘ

    θᴾᵐ   = (θᴾ + θᴾⁿ) / 2
    dθᴾdΘ = - FT(cᵖ⁰) / ((θᴾᵐ + FT(T₀)) * zero_pressure_gibbs_curvature(Sᴬ, θᴾᵐ))

    θᴾ  = θᴾⁿ - ΔΘ * dθᴾdΘ
    ΔΘ  = Θ_from_θᴾ(Sᴬ, θᴾ) - Θ
    θᴾⁿ = θᴾ
    θᴾ  = θᴾⁿ - ΔΘ * dθᴾdΘ

    return θᴾ
end
