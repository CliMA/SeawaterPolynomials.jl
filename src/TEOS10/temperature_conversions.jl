#####
##### TEOS-10 temperature conversions (Θ ↔ θᴾ ↔ T).
#####
##### Direct Julia translations of routines from `gsw_oceanographic_toolbox.c` of
##### https://github.com/TEOS-10/GSW-C, the canonical C implementation maintained by the TEOS-10 team.
##### The salinity conversions (`Sᴬ_from_Sᴾ`) and the SAAR atlas live in `salinity_conversions.jl`.
#####
##### Notation (TEOS-10 manual, §A.1):
#####   Sᴬ : absolute salinity                      [g/kg]
#####   Sᴾ : practical salinity, PSS-78             [unitless]
#####   T  : in-situ temperature, ITS-90            [°C]
#####   θᴾ : potential temperature, p_ref = 0 dbar  [°C]
#####   Θ  : conservative temperature               [°C]
#####   p  : sea pressure (gauge)                   [dbar]
#####   λ  : longitude                              [°E]
#####   φ  : latitude                               [°N]
#####
##### Public conversions follow the pattern `<output>_from_<input>`, e.g. `Θ_from_θᴾ`, `θᴾ_from_T`, `Sᴬ_from_Sᴾ`. 
#####

#####
##### Constants (gsw_internal_const.h, GSW_TEOS10_CONSTANTS macro)
#####

const T₀ = 273.15        # gsw_t0   [K]: celsius zero point
const Sₒ = 35.16504      # gsw_sso  [g/kg]: standard ocean salinity
const uₚₛ = Sₒ / 35       # gsw_ups: PSS-78 → absolute-salinity scale factor
const rS = 1 / (40 * uₚₛ) # gsw_rS: x² = rS · Sᴬ

#####
##### Polynomial kernels for the temperature conversions.
##### `_η`, `_η⁰` and `_∂θᴾ²_g` are direct translations of the GSW-C routines `gsw_entropy_part`,
##### `gsw_entropy_part_zerop` and `gsw_gibbs_pt0_pt0`. The last is `∂²g/∂(θᴾ)²` evaluated at `p = 0`,
##### and `_η⁰` denotes specific entropy evaluated at `p = 0` (i.e. ``η(Sᴬ, θᴾ, 0)``).
#####

@inline function _η(Sᴬ::FT, T::FT, p::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = T * FT(0.025)
    z  = p * FT(1e-4)

    # g₀₃: salinity-independent terms, organised as a Horner in y whose yʲ-coefficient is itself a Horner in z.
    g₃₀ =                           z*(FT(-270.983805184062)  + z*(FT(776.153611613101)   + z*(FT(-196.51255088122)   + z*(FT(28.9796526294175)   - z*FT(2.13290083518327) ))))
    g₃₁ = FT(-24715.571866078)    + z*(FT(2910.0729080936)    + z*(FT(-1513.116771538718) + z*(FT(546.959324647056)   + z*(FT(-111.1208127634436) + z*FT(8.68841343834394) ))))
    g₃₂ = FT(2210.2236124548363)  + z*(FT(-2017.52334943521)  + z*(FT(1498.081172457456)  + z*(FT(-718.6359919632359) + z*(FT(146.4037555781616)  - z*FT(4.9892131862671505)))))
    g₃₃ = FT(-592.743745734632)   + z*(FT(1591.873781627888)  + z*(FT(-1207.261522487504) + z*(FT(608.785486935364)   - z* FT(105.4993508931208))))
    g₃₄ = FT(290.12956292128547)  + z*(FT(-973.091553087975)  + z*(FT(602.603274510125)   + z*(FT(-276.361526170076)  + z* FT(32.40953340386105))))
    g₃₅ = FT(-113.90630790850321) + z*(FT(381.06836198507096) + z*(FT(-133.7383902842754) + z* FT(49.023632509086724)))
    g₃₆ = FT(21.35571525415769)   - z* FT(67.41756835751434)

    g₀₃ = g₃₀ + y*(g₃₁ + y*(g₃₂ + y*(g₃₃ + y*(g₃₄ + y*(g₃₅ + y*g₃₆)))))

    # g₀₈: salinity-dependent terms, factored as x² · (g₈ᶻ + x · g₈ˣ + y · g₈ʸ).
    g₈ᶻ = z*(FT(729.116529735046) + z*(FT(-343.956902961561) + z*(FT(124.687671116248) + z*(FT(-31.656964386073) + z*FT(7.04658803315449)))))

    g₈ˣ = x*(y*(FT(-137.1145018408982) + y*(FT(148.10030845687618) + y*(FT(-68.5590309679152)  + y*FT(12.4848504784754)))) - z*FT(22.6683558512829)) +
          z*(FT(-175.292041186547)     + z*(FT(83.1923927801819)   - z*FT(29.483064349429)))   +
          y*(FT(-86.1329351956084)     + z*(FT(766.116132004952)   + z*(FT(-108.3834525034224) + z*FT(51.2796974779828))) +
          y*(FT(-30.0682112585625)     - z*FT(1380.9597954037708)  + y*(FT(3.50240264723578)   + z*FT(938.26075044542))))

    g₈ʸ =    FT(1760.062705994408)   + y*(FT(-675.802947790203) +
          y*(FT(365.7041791005036)   +
          y*(FT(-108.30162043765552) + y*FT(12.78101825083098)) +
          z*(FT(-1190.914967948748)  + z*(FT(298.904564555024)  - z*FT(145.9491676006352)))) +
          z*(FT(2082.7344423998043)  + z*(FT(-614.668925894709) + z*(FT(340.685093521782)    - z*FT(33.3848202979239))))) +
          z*(FT(-1721.528607567954)  + z*(FT(674.819060538734)  + z*(FT(-356.629112415276)   + z*(FT(88.4080716616)       - z*FT(15.84003094423364)))))

    g₀₈ = x² * (g₈ᶻ + x*g₈ˣ + y*g₈ʸ)

    return -(g₀₃ + g₀₈) * FT(0.025)
end

@inline function _η⁰(Sᴬ::FT, θ₀::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = θ₀ * FT(0.025)

    # g₀₃ at z = 0: pure Horner in y, starting at y¹.
    g₀₃ = y*(FT(-24715.571866078)    + y*(FT(2210.2236124548363) +
          y*(FT(-592.743745734632)   + y*(FT(290.12956292128547) +
          y*(FT(-113.90630790850321) + y* FT(21.35571525415769))))))

    # g₀₈ at z = 0: x² · (g₈ᶻ + x · g₈ˣ + y · g₈ʸ) with z dropped.
    g₈ˣ = x*(y*(FT(-137.1145018408982) + y*(FT(148.10030845687618) +
          y*(FT(-68.5590309679152)     +    FT(12.4848504784754)*y)))) +
          y*(FT(-86.1329351956084)     + y*(FT(-30.0682112585625) +
                                         y* FT(3.50240264723578)))

    g₈ʸ = FT(1760.062705994408)    + y*(FT(-675.802947790203) +
          y*(FT(365.7041791005036) + y*(FT(-108.30162043765552) +
                                     y* FT(12.78101825083098))))

    g₀₈ = x²*(x*g₈ˣ + y*g₈ʸ)

    return -(g₀₃ + g₀₈) * FT(0.025)
end

@inline function _∂θᴾ²_g(Sᴬ::FT, θ₀::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = θ₀ * FT(0.025)

    # g₀₃ second derivative w.r.t. θ at z = 0: pure Horner in y.
    g₀₃ =    FT(-24715.571866078)   + y*(FT(4420.4472249096725) +
          y*(FT(-1778.231237203896) + y*(FT(1160.5182516851419) +
          y*(FT(-569.531539542516)  + y*FT(128.13429152494615)))))

    # g₀₈ second derivative w.r.t. θ at z = 0: x² · (g₈ᵃ + x · g₈ˣ + y · g₈ʸ).
    g₈ᵃ = FT(1760.062705994408)
    g₈ˣ = FT(-86.1329351956084)     + x*(FT(-137.1145018408982)  +
          y*(FT(296.20061691375236) + y*(FT(-205.67709290374563) +
                                      y* FT(49.9394019139016)))) +
          y*(FT(-60.136422517125)   + y* FT(10.50720794170734))
    
    g₈ʸ = FT(-1351.605895580406)     + y*(FT(1097.1125373015109) +
          y*(FT(-433.20648175062206) + y* FT(63.905091254154904)))

    g₀₈ = x²*(g₈ᵃ + x*g₈ˣ + y*g₈ʸ)

    return (g₀₃ + g₀₈) * FT(0.000625)
end

#####
##### Conservative temperature from potential temperature (gsw_ct_from_pt)
#####

"""
    Θ_from_θᴾ(Sᴬ, θᴾ)

Return the TEOS-10 conservative temperature ``Θ`` from absolute salinity ``Sᴬ`` and potential temperature
``θᴾ`` referenced to ``p = 0`` dbar, computed as ``Θ = h⁰(Sᴬ, θᴾ) / cₚ⁰`` where ``h⁰`` is the 27-term
potential-enthalpy polynomial and ``cₚ⁰ = 3991.867957...`` J/kg/K is the TEOS-10 reference heat capacity.
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
    FT = promote_type(typeof(Sᴬ), typeof(θᴾ))
    return _Θ_from_θᴾ(convert(FT, Sᴬ), convert(FT, θᴾ))
end

@inline function _Θ_from_θᴾ(Sᴬ::FT, θᴾ::FT) where FT
    x² = FT(rS) * Sᴬ
    x  = sqrt(x²)
    y  = θᴾ * FT(0.025)

    # Salinity-dependent (x²) branch, organised by x-power (x², x³, x⁴, x⁵, x⁶).
    h₂ =      FT(268.5520265845071)   + y * (FT(-12019.028203559312) +
         y * (FT(3734.858026725145)   + y * (FT(-2046.7671145057618) +
         y * (FT(465.28655623826234)  + y * (FT(-0.6370820302376359) -
         y *  FT(10.650848542359153))))))

    h₃ =      FT(937.2099110620707)   + y * (FT(588.1802812170108)  +
         y * (FT(248.39476522971285)  + y * (FT(-3.871557904936333) -
         y *  FT(2.6268019854268356))))

    h₄ =      FT(-1687.914374187449)  + y * (FT(936.3206544460336) +
         y * (FT(-942.7827304544439)  + y * (FT(369.4389437509002) +
         y * (FT(-33.83664947895248)  - y *  FT(9.987880382780322)))))

    h₅ = FT(246.9598888781377)
    h₆ = FT(123.59576582457964) - x*FT(48.5891069025409)

    # Potential enthalpy: salinity-independent y-Horner plus the x²-leading branch in x.
    h⁰ =      FT(61.01362420681071)   + y * (FT(168776.46138048015) +
         y * (FT(-2735.2785605119625) + y * (FT(2574.2164453821433) +
         y * (FT(-1536.6644434977543) + y * (FT(545.7340497931629)  +
         y * (FT(-50.91091728474331)  - y *  FT(18.30489878927802)))))))   +
         x² * (h₂ + x * (h₃ + x * (h₄ + x * (h₅ + x * h₆))))

    return h⁰ / FT(cₚ⁰)
end

#####
##### Potential temperature from in-situ temperature (gsw_pt0_from_t)
#####

"""
    θᴾ_from_T(Sᴬ, T, p)

Return the TEOS-10 potential temperature ``θᴾ`` referenced to ``p = 0`` dbar from absolute salinity,
in-situ temperature and sea pressure. A polynomial first guess is followed by two iterations of a
modified Newton–Raphson update on the specific entropy, yielding agreement with the GSW-C reference to
within machine precision throughout the oceanographic state space. Direct translation of `gsw_pt0_from_t`
of https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴬ`: absolute salinity                                       [g/kg]
- `T` : in-situ temperature, ITS-90                             [°C]
- `p` : sea pressure (gauge: absolute pressure − 10.1325 dbar)  [dbar]

# Output
- `θᴾ`: potential temperature, p_ref = 0 dbar                   [°C]

# References
- IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater – 2010.
  http://www.teos-10.org/pubs/TEOS-10_Manual.pdf
"""
@inline function θᴾ_from_T(Sᴬ, T, p)
    FT = promote_type(typeof(Sᴬ), typeof(T), typeof(p))
    return _θᴾ_from_T(convert(FT, Sᴬ), convert(FT, T), convert(FT, p))
end

@inline function _θᴾ_from_T(Sᴬ::FT, T::FT, p::FT) where FT
    s₁ = Sᴬ / FT(uₚₛ)

    θ₀ = T + p*(FT(8.65483913395442e-6)  -
             s₁*FT(1.41636299744881e-6)  -
             p *FT(7.38286467135737e-9)  +
             T *(FT(-8.38241357039698e-6) +
             s₁*FT(2.83933368585534e-8)  +
             T *FT(1.77803965218656e-8)  +
             p *FT(1.71155619208233e-10)))

    dθ_η = ((FT(T₀) + θ₀) * (one(FT) - FT(0.05) * (one(FT) - Sᴬ/FT(Sₒ)))) / FT(cₚ⁰)

    η = _η(Sᴬ, T, p)

    for _ in 1:2
        θ₀ⁿ  = θ₀
        Δη   = _η⁰(Sᴬ, θ₀ⁿ) - η
        θ₀   = θ₀ⁿ - Δη * dθ_η
        θ₀ᵐ  = FT(0.5) * (θ₀ + θ₀ⁿ)
        dθ_η = - one(FT) / _∂θᴾ²_g(Sᴬ, θ₀ᵐ)
        θ₀   = θ₀ⁿ - Δη * dθ_η
    end

    return θ₀
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
@inline function Θ_from_T(Sᴬ, T, p)
    θᴾ = θᴾ_from_T(Sᴬ, T, p)
    return Θ_from_θᴾ(Sᴬ, θᴾ)
end
