#####
##### Absolute salinity from practical salinity, translated from `gsw_sa_from_sp` of https://github.com/TEOS-10/GSW-C:
##### an Absolute Salinity Anomaly Ratio (SAAR) atlas interpolated trilinearly, plus an analytic Baltic Sea correction.
#####
##### Notation:
#####   Sᴾ : practical salinity, PSS-78  [unitless]
#####   λ  : longitude                   [°E]
#####   φ  : latitude                    [°N]
#####

const gsw_error_limit         = 1.0e10
const gsw_neighbour_threshold = 100.0   # cut-off between valid SAAR ratios (O(0.01)) and the ±9e90 land flags
const panama_inset            = 0.001   # keeps samples strictly inside the Panama barrier domain

# Index of the interval of the sorted `xs` that brackets `x`, clamped to the first and last intervals (`gsw_util_indx`)
@inline interval_index(xs, x) = clamp(count(≤(x), xs), 1, length(xs) - 1)

@inline function linear_interpolate(xs, ys, x)
    k = interval_index(xs, x)
    r = (x - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] + r * (ys[k+1] - ys[k])
end

#####
##### SAAR reference atlas
#####
##### The atlas is a headerless little-endian Float64 blob holding, back to back, `p [Np]`, `φ [Nφ]`, `λ [Nλ]`,
##### `saar [Np × Nφ × Nλ]` and `ndepth [Nφ × Nλ]`. `saar[k, j, i]` is GSW-C's `saar_ref[idz0 + Np*(idy0 + Nφ*idx0)]`
##### with `(k, j, i) = (idz0, idy0, idx0) .+ 1`. Land points carry the flag 9e90.
#####

const Nλ = 91
const Nφ = 45
const Np = 45

"""
    SAARAtlas{V, A, M}

TEOS-10 Absolute Salinity Anomaly Ratio reference atlas.

# Fields
- `p`     : pressure levels  `[dbar]`,  length `Np`
- `φ`     : latitudes        `[°N]`,    length `Nφ`
- `λ`     : longitudes       `[°E]`,    length `Nλ`
- `saar`  : absolute salinity anomaly ratio,                size `Np × Nφ × Nλ`
- `ndepth`: maximum valid depth-index per `(φ, λ)` column,  size `Nφ × Nλ`
"""
struct SAARAtlas{V, A, M}
    p      :: V
    φ      :: V
    λ      :: V
    saar   :: A
    ndepth :: M
end

Adapt.@adapt_structure SAARAtlas

function SAARAtlas(path::AbstractString)
    p      = Array{Float64, 1}(undef, Np)
    φ      = Array{Float64, 1}(undef, Nφ)
    λ      = Array{Float64, 1}(undef, Nλ)
    saar   = Array{Float64, 3}(undef, Np, Nφ, Nλ)
    ndepth = Array{Float64, 2}(undef, Nφ, Nλ)
    open(path, "r") do io
        read!(io, p)
        read!(io, φ)
        read!(io, λ)
        read!(io, saar)
        read!(io, ndepth)
    end
    return SAARAtlas(p, φ, λ, saar, ndepth)
end

const SAAR_ATLAS = Ref{SAARAtlas{Vector{Float64}, Array{Float64, 3}, Matrix{Float64}}}()

#####
##### Panama isthmus barrier and atlas cell corners
#####

const panama_longitudes = (260.00, 272.59, 276.50, 278.65, 280.73, 292.00)
const panama_latitudes  = ( 19.55,  13.97,   9.60,   8.10,   9.33,   3.40)

const longitude_corner_offsets = (0, 1, 1, 0)
const latitude_corner_offsets  = (0, 0, 1, 1)

# Replace the flagged corners of `d` by the mean of the valid ones (`gsw_add_mean`)
@inline function mean_of_valid_neighbours(d::NTuple{4, FT}) where FT
    N = 0
    S = zero(FT)
    for k in 1:4
        N += ifelse(abs(d[k]) <= gsw_neighbour_threshold, 1, 0)
        S += ifelse(abs(d[k]) <= gsw_neighbour_threshold, d[k], zero(FT))
    end
    m = ifelse(N == 0, zero(FT), S / N)
    return ntuple(k -> ifelse(abs(d[k]) >= gsw_neighbour_threshold, m, d[k]), Val(4))
end

# Replace the corners of `d` that are flagged or lie across the Panama isthmus from `(λ, φ)` by the mean of the rest
# (`gsw_add_barrier`). `(λr, φr)` is the south-west corner of the atlas cell and `(Δλ, Δφ)` its extent.
@inline function apply_panama_barrier(d::NTuple{4, FT}, λ, φ, λr, φr, Δλ, Δφ) where FT
    λP = panama_longitudes
    φP = panama_latitudes

    sample_side = linear_interpolate(λP, φP, λ) ≤ φ
    φᵂ = linear_interpolate(λP, φP, λr)
    φᴱ = linear_interpolate(λP, φP, λr + Δλ)
    side = (φᵂ ≤ φr, φᴱ ≤ φr, φᴱ ≤ φr + Δφ, φᵂ ≤ φr + Δφ)

    N = 0
    S = zero(FT)
    for k in 1:4
        valid = abs(d[k]) <= gsw_neighbour_threshold && sample_side == side[k]
        N += ifelse(valid, 1, 0)
        S += ifelse(valid, d[k], zero(FT))
    end
    m = ifelse(N == 0, zero(FT), S / N)

    return ntuple(k -> ifelse(abs(d[k]) >= gsw_error_limit || sample_side != side[k], m, d[k]), Val(4))
end

#####
##### Baltic correction (gsw_sa_from_sp_baltic)
#####

const baltic_west_longitudes = (12.6, 7.0, 26.0)
const baltic_west_latitudes  = (50.0, 59.0, 69.0)
const baltic_east_longitudes = (45.0, 26.0)
const baltic_east_latitudes  = (50.0, 69.0)

"""
    baltic_absolute_salinity(Sᴾ, λ, φ)

Return the absolute salinity [g/kg] inside the Baltic Sea polygon, or `NaN` outside it. Mirrors `gsw_sa_from_sp_baltic`.

# References
- Feistel, R., Marion, G. M., Pawlowicz, R., Wright, D. G., 2010: Thermophysical property anomalies of
  Baltic seawater. Ocean Science 6, 949–981.
"""
function baltic_absolute_salinity(Sᴾ, λ, φ)
    λᵂ = baltic_west_longitudes
    φᵂ = baltic_west_latitudes
    λᴱ = baltic_east_longitudes
    φᴱ = baltic_east_latitudes

    inside_polygon = λᵂ[2] < λ < λᴱ[1] && φᵂ[1] < φ < φᵂ[3]
    inside_polygon || return NaN

    λᴸ = linear_interpolate(φᵂ, λᵂ, φ)
    λᴿ = linear_interpolate(φᴱ, λᴱ, φ)

    return ifelse(λᴸ ≤ λ ≤ λᴿ, ((Sₒ - 0.087) / 35) * Sᴾ + 0.087, NaN)
end

#####
##### Absolute Salinity Anomaly Ratio lookup (gsw_saar)
#####

"""
    saar(atlas, p, λ, φ)

Return the Absolute Salinity Anomaly Ratio at sea pressure `p` [dbar], longitude `λ ∈ [0, 360)` [°E] and latitude
`φ` [°N], interpolated trilinearly from `atlas`. Returns `NaN` where the atlas cannot resolve the point and `0` far
from any ocean column. Mirrors `gsw_saar`.

# References
- McDougall, T. J., Jackett, D. R., Millero, F. J., Pawlowicz, R., Barker, P. M., 2012: A global
  algorithm for estimating absolute salinity. Ocean Science 8, 1123–1134.
"""
function saar(atlas::SAARAtlas, p, λ, φ)
    (isnan(p) || isnan(λ) || isnan(φ) || φ < -86 || φ > 90) && return NaN

    λP = panama_longitudes
    φP = panama_latitudes
    δi = longitude_corner_offsets
    δj = latitude_corner_offsets

    i = min(floor(Int, (Nλ - 1) * (λ - atlas.λ[1]) / (atlas.λ[end] - atlas.λ[1])) + 1, Nλ - 1)
    j = min(floor(Int, (Nφ - 1) * (φ - atlas.φ[1]) / (atlas.φ[end] - atlas.φ[1])) + 1, Nφ - 1)

    Nᵐᵃˣ = -1.0
    for c in 1:4
        n = atlas.ndepth[j + δj[c], i + δi[c]]
        if 0 < n < 1e90
            Nᵐᵃˣ = max(Nᵐᵃˣ, n)
        end
    end

    Nᵐᵃˣ == -1 && return 0.0

    p = min(p, atlas.p[Int(Nᵐᵃˣ)])
    k = interval_index(atlas.p, p)

    ξ = (λ - atlas.λ[i]) / (atlas.λ[i+1] - atlas.λ[i])
    η = (φ - atlas.φ[j]) / (atlas.φ[j+1] - atlas.φ[j])
    ζ = (p - atlas.p[k]) / (atlas.p[k+1] - atlas.p[k])

    inside_panama = (λP[1] ≤ λ ≤ λP[end] - panama_inset) && (φP[end] ≤ φ ≤ φP[1])

    rS⁺ = level_interpolation(atlas, i, j, k,   λ, φ, ξ, η, inside_panama)
    rS⁻ = level_interpolation(atlas, i, j, k+1, λ, φ, ξ, η, inside_panama)
    rS⁻ = ifelse(abs(rS⁻) >= gsw_error_limit, rS⁺, rS⁻)

    out = rS⁺ + ζ * (rS⁻ - rS⁺)
    return ifelse(abs(out) >= gsw_error_limit, NaN, out)
end

# Bilinear interpolation on pressure level `k`, after filling flagged or cross-isthmus corners
@inline function level_interpolation(atlas, i, j, k, λ, φ, ξ, η, inside_panama)
    λr = atlas.λ[i]
    φr = atlas.φ[j]
    Δλ = atlas.λ[2] - atlas.λ[1]
    Δφ = atlas.φ[2] - atlas.φ[1]
    δi = longitude_corner_offsets
    δj = latitude_corner_offsets

    corners = ntuple(c -> atlas.saar[k, j + δj[c], i + δi[c]], Val(4))

    corners = if inside_panama
        apply_panama_barrier(corners, λ, φ, λr, φr, Δλ, Δφ)
    elseif abs(sum(corners)) >= gsw_error_limit
        mean_of_valid_neighbours(corners)
    else
        corners
    end

    return (1 - η) * (corners[1] + ξ * (corners[2] - corners[1])) +
                η  * (corners[4] + ξ * (corners[3] - corners[4]))
end

#####
##### Absolute salinity from practical salinity (gsw_sa_from_sp)
#####

"""
    Sᴬ_from_Sᴾ(Sᴾ, p, λ, φ, atlas = SAAR_ATLAS[])

Return the TEOS-10 absolute salinity ``Sᴬ`` from practical salinity, sea pressure, longitude and latitude,

```math
Sᴬ = (Sₒ/35)\\, Sᴾ\\, (1 + r_{\\mathrm{SAAR}}) ,
```

where ``Sₒ = 35.16504`` g/kg is the standard ocean reference salinity and ``r_{\\mathrm{SAAR}}`` is interpolated
from `atlas`. Inside the Baltic Sea an analytic correction replaces the atlas. Inside a GPU kernel, pass the atlas
moved to the device, e.g. `Adapt.adapt(CuArray, SAAR_ATLAS[])`. Translation of `gsw_sa_from_sp` of
https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴾ`   : practical salinity, PSS-78  [unitless]
- `p`    : sea pressure                [dbar]
- `λ`    : longitude                   [°E]
- `φ`    : latitude                    [°N]
- `atlas`: `SAARAtlas`

# Output
- `Sᴬ`: absolute salinity [g/kg], `NaN` where the atlas cannot resolve the point

# References
- McDougall, T. J., Jackett, D. R., Millero, F. J., Pawlowicz, R., Barker, P. M., 2012: A global
  algorithm for estimating absolute salinity. Ocean Science 8, 1123–1134.
"""
function Sᴬ_from_Sᴾ(Sᴾ, p, λ, φ, atlas = SAAR_ATLAS[])
    FT = typeof(float(Sᴾ))
    λ  = mod(float(λ), 360)

    baltic_salinity = baltic_absolute_salinity(Sᴾ, λ, φ)
    isnan(baltic_salinity) || return convert(FT, baltic_salinity)

    r = saar(atlas, float(p), λ, float(φ))
    return convert(FT, uₚₛ * Sᴾ * (1 + r))
end
