#####
##### Absolute salinity from practical salinity (`Sᴬ_from_Sᴾ`).
#####
##### Provides the TEOS-10 mapping ``Sᴬ(Sᴾ, p, λ, φ)`` via
#####   1. an Absolute Salinity Anomaly Ratio (SAAR) atlas, trilinearly interpolated on a global grid
#####      held in a `SAARAtlas` struct,
#####   2. an analytic Baltic Sea correction.
##### Direct port of `gsw_sa_from_sp` of https://github.com/TEOS-10/GSW-C, with the atlas wrapped in a
##### struct rather than as global mutable arrays.
#####

#####
##### Sentinel values returned by GSW-C atlas routines for points outside the ocean (originally Matlab
##### NaNs replaced with ±9e90 for `ndepth_ref` / `saar_ref`, and 9e15 for invalid scalar returns).
#####

const gsw_invalid_value = 9.0e15
const gsw_error_limit   = 1.0e10

@inline invalid(x) = abs(x) >= gsw_error_limit

#####
##### SAAR reference atlas
#####
##### The bundled binary blob (little-endian Float64, no header) packs five arrays back-to-back:
#####
#####     p_ref      [Np]                      pressure levels                                [dbar]
#####     φ_ref      [Nφ]                      latitudes                                      [°N]
#####     λ_ref      [Nλ]                      longitudes                                     [°E]
#####     saar_ref   [Np × Nφ × Nλ]            absolute salinity anomaly ratio                [—]
#####     ndepth_ref [Nφ × Nλ]                 maximum valid depth-index per (φ, λ) column    [—]
#####
##### The Julia column-major layout `arr[k, j, i]` mirrors the GSW-C linear index
##### `idz0 + Np*(idy0 + Nφ*idx0)` once 1-based.
#####

const Nλ = 91
const Nφ = 45
const Np = 45

"""
    SAARAtlas{V, A, M}

Container for the TEOS-10 Absolute Salinity Anomaly Ratio reference atlas.

# Fields
- `p`     : pressure levels                                [dbar]   (length `Np`)
- `φ`     : latitudes                                      [°N]     (length `Nφ`)
- `λ`     : longitudes                                     [°E]     (length `Nλ`)
- `saar`  : absolute salinity anomaly ratio                [—]      (`Np × Nφ × Nλ`)
- `ndepth`: maximum valid depth-index per `(φ, λ)` column  [—]      (`Nφ × Nλ`)
"""
struct SAARAtlas{V <: AbstractVector{Float64},
                 A <: AbstractArray{Float64,3},
                 M <: AbstractMatrix{Float64}}
    p      :: V
    φ      :: V
    λ      :: V
    saar   :: A
    ndepth :: M
end

"""
    load_saar_atlas(path)

Read a SAAR atlas from the binary blob at `path` and return a fully populated [`SAARAtlas`](@ref). The
expected blob layout is described in the file header.
"""
function load_saar_atlas(path)
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

const SAAR_ATLAS = Ref{SAARAtlas{Vector{Float64}, Array{Float64,3}, Matrix{Float64}}}()
const saar_data_path = joinpath(@__DIR__, "../../", "data", "gsw_saar_data.bin")

#####
##### Atlas geometry: Panama isthmus barrier and corner offsets for the bilinear cell.
#####

const panama_longitudes = [260.00, 272.59, 276.50, 278.65, 280.73, 292.00]
const panama_latitudes  = [ 19.55,  13.97,   9.60,   8.10,   9.33,   3.40]

# Corner offsets used to build the (λ, φ) interpolation cell.
const atlas_cell_offsets_i = [0, 1, 1, 0]
const atlas_cell_offsets_j = [0, 0, 1, 1]

#####
##### Atlas helpers
#####

"""
    mean_of_valid_neighbours(d)

Replace each ±9e90 invalid-flag entry of the four-element atlas-corner tuple `d` with the arithmetic
mean of the valid neighbours. Mirrors `gsw_add_mean`.

# Inputs
- `d`: four atlas-corner values (NTuple{4})

# Output
- `d′`: 4-tuple with invalid corners replaced by the mean of valid ones
"""
@inline function mean_of_valid_neighbours(d::NTuple{4, FT}) where FT
    N = 0
    S = zero(FT)
    for k in 1:4
        N += ifelse(abs(d[k]) <= 100, 1, 0)
        S += ifelse(abs(d[k]) <= 100, d[k], zero(FT))
    end
    m = ifelse(N == 0, zero(FT), S / N)
    return ntuple(k -> ifelse(abs(d[k]) >= 100, m, d[k]), Val(4))
end

"""
    apply_panama_barrier(d, lon, lat, λr, φr, Δλ, Δφ)

In the Central American region the four atlas neighbours of a single `(lon, lat)` point may straddle the
Panama isthmus. Corners that lie on the wrong side of the barrier — or that carry a ±9e90 invalid flag —
are replaced by the arithmetic mean of the valid same-side corners. Mirrors `gsw_add_barrier`.

# Inputs
- `d`      : four atlas-corner values (NTuple{4})
- `λ,  φ`  : sample location                       [°E, °N]
- `λr, φr` : south-west corner of the atlas cell   [°E, °N]
- `Δλ, Δφ` : cell extents                          [°]

# Output
- `d′`: 4-tuple with wrong-side or invalid corners replaced by the mean of the rest
"""
@inline function apply_panama_barrier(d::NTuple{4, T}, λ, φ, λr, φr, Δλ, Δφ) where T

    λP = panama_longitudes
    φP = panama_latitudes

    # Side of the isthmus the (λ, φ) sample lies on.
    k  = searchsortedlast(λP, λ)
    r  = (λ - λP[k]) / (λP[k+1] - λP[k])
    sample_side = (φP[k] + r * (φP[k+1] - φP[k])) ≤ φ

    # Side of the isthmus for the cell corners (1, 4) on the western edge.
    k  = searchsortedlast(λP, λr)
    r  = (λr - λP[k]) / (λP[k+1] - λP[k])
    φW = φP[k] + r * (φP[k+1] - φP[k])

    # Side of the isthmus for the cell corners (2, 3) on the eastern edge.
    k  = searchsortedlast(λP, λr + Δλ)
    r  = (λr + Δλ - λP[k]) / (λP[k+1] - λP[k])
    φE = φP[k] + r * (φP[k+1] - φP[k])

    side1 = φW ≤  φr
    side4 = φW ≤ (φr + Δφ)
    side2 = φE ≤  φr
    side3 = φE ≤ (φr + Δφ)

    side = (side1, side2, side3, side4)

    nmean = 0
    s = zero(T)

    for k in 1:4
        same_side = sample_side == side[k]
        valid     = abs(d[k]) <= 100 && same_side
        nmean    += ifelse(valid, 1, 0)
        s        += ifelse(valid, d[k], zero(T))
    end

    m = ifelse(nmean == 0, zero(T), s / nmean)

    return ntuple(k -> (abs(d[k]) >= T(1e10) || sample_side != side[k]) ? m : d[k], Val(4))
end

"""
    linear_interpolate(xs::NTuple, ys::NTuple, x0)

Piecewise-linear interpolation of the small static table `(xs, ys)` at `x0`. Used for the Baltic Sea
boundary polygons in [`baltic_absolute_salinity`](@ref). Mirrors `gsw_util_xinterp1` for tuple inputs.

# Inputs
- `xs`: monotonically increasing abscissae (NTuple{N})
- `ys`: ordinates at `xs` (NTuple{N})
- `x0`: query abscissa

# Output
- `y0`: linear interpolant of `(xs, ys)` at `x0`
"""
@inline function linear_interpolate(xs, ys, x0)
    k = searchsortedlast(xs, x0)
    return ys[k] + (x0 - xs[k]) * (ys[k+1] - ys[k]) / (xs[k+1] - xs[k])
end

#####
##### Baltic correction (gsw_sa_from_sp_baltic)
#####

# Baltic Sea boundary polygons, as in `gsw_sa_from_sp_baltic`.
const baltic_west_longitudes = [12.6, 7.0, 26.0]
const baltic_west_latitudes  = [50.0, 59.0, 69.0]
const baltic_east_longitudes = [45.0, 26.0]
const baltic_east_latitudes  = [50.0, 69.0]

"""
    baltic_absolute_salinity(Sᴾ, λ, φ)

Return the analytic absolute-salinity correction inside the Baltic Sea polygon, or `gsw_invalid_value`
if `(λ, φ)` lies outside the polygon. Mirrors `gsw_sa_from_sp_baltic`.

# Inputs
- `Sᴾ`: practical salinity, PSS-78  [unitless]
- `λ` : longitude                   [°E]
- `φ` : latitude                    [°N]

# Output
- `Sᴬ`: absolute salinity inside the Baltic polygon [g/kg], or `gsw_invalid_value` outside it

# References
- Feistel, R., Marion, G. M., Pawlowicz, R., Wright, D. G., 2010: Thermophysical property anomalies of
  Baltic seawater. Ocean Science 6, 949–981.
"""
function baltic_absolute_salinity(Sᴾ, λ, φ)
    λ  = mod(λ, 360.0)
    λᵂ = baltic_west_longitudes
    φᵂ = baltic_west_latitudes
    λᴱ = baltic_east_longitudes
    φᴱ = baltic_east_latitudes

    inside_polygon = λᵂ[2] < λ && λ < λᴱ[1] && φᵂ[1] < φ && φ < φᵂ[3]
    inside_polygon || return gsw_invalid_value

    λᴸ = linear_interpolate(φᵂ, λᵂ, φ)
    λᴿ = linear_interpolate(φᴱ, λᴱ, φ)

    return ifelse(λᴸ ≤ λ ≤ λᴿ, ((Sₒ - 0.087) / 35.0) * Sᴾ + 0.087, gsw_invalid_value)
end

#####
##### Absolute Salinity Anomaly Ratio lookup (gsw_saar)
#####

"""
    saar(p, λ, φ)

Return the Absolute Salinity Anomaly Ratio at sea pressure, longitude and latitude. The atlas is
trilinearly interpolated on the bundled `Nλ × Nφ × Np` grid, with NaN-flag handling via
[`mean_of_valid_neighbours`](@ref) and Panama-isthmus handling via [`apply_panama_barrier`](@ref).
Mirrors `gsw_saar`.

# Inputs
- `p`: sea pressure  [dbar]
- `λ`: longitude     [°E]
- `φ`: latitude      [°N]

# Output
- `r_SAAR`: absolute salinity anomaly ratio [unitless], `gsw_invalid_value` for points the atlas cannot
  resolve, or `0.0` far from any ocean column

# References
- McDougall, T. J., Jackett, D. R., Millero, F. J., Pawlowicz, R., Barker, P. M., 2012: A global
  algorithm for estimating absolute salinity. Ocean Science 8, 1123–1134.
"""
function saar(p, λ, φ)
    
    if isnan(φ) || isnan(λ) || isnan(p)
        return gsw_invalid_value
    end
    if φ < -86.0 || φ > 90.0
        return gsw_invalid_value
    end


    λP  = panama_longitudes
    φP  = panama_latitudes
    Atl = SAAR_ATLAS[]
    λ   = mod(λ, 360.0)
    Δλ  = Atl.λ[2] - Atl.λ[1]
    Δφ  = Atl.φ[2] - Atl.φ[1]

    # 1-based lower-corner indices for the (λ, φ) cell.
    i = floor(Int, (Nλ - 1) * (λ - Atl.λ[1]) / (Atl.λ[end] - Atl.λ[1])) + 1
    j = floor(Int, (Nφ - 1) * (φ - Atl.φ[1]) / (Atl.φ[end] - Atl.φ[1])) + 1
    
    # wrap around
    i = ifelse(i == Nλ, Nλ - 1, i)
    j = ifelse(j == Nφ, Nφ - 1, j)

    δi = atlas_cell_offsets_i
    δj = atlas_cell_offsets_j

    # Maximum depth-index of the four atlas corners
    Nmax = -1.0
    for c in 1:4
        nd = Atl.ndepth[j + δj[c], i + δi[c]]
        if 0.0 < nd < 1e90
            Nmax = max(Nmax, nd)
        end
    end

    Nmax == -1.0 && return 0.0   # well outside any ocean column → SAAR = 0

    p = ifelse(p > Atl.p[Int(Nmax)], Atl.p[Int(Nmax)], p)
    k = searchsortedlast(Atl.p, p)

    ξ = (λ - Atl.λ[i]) / (Atl.λ[i+1] - Atl.λ[i])
    η = (φ - Atl.φ[j]) / (Atl.φ[j+1] - Atl.φ[j])
    ζ = (p - Atl.p[k]) / (Atl.p[k+1] - Atl.p[k])

    inside_panama = (λP[1] ≤ λ ≤ λP[end] - 0.001) && (φP[end] ≤ φ ≤ φP[1])

    rS⁺ = level_interpolation(Atl, i, j, k,   λ, φ, ξ, η, inside_panama)
    rS⁻ = level_interpolation(Atl, i, j, k+1, λ, φ, ξ, η, inside_panama)

    rS⁻ = ifelse(abs(rS⁻) >= gsw_error_limit, rS⁺, rS⁻)

    out = rS⁺ + ζ * (rS⁻ - rS⁺)
    return ifelse(abs(out) >= gsw_error_limit, gsw_invalid_value, out)
end

# Bilinear interpolation of the four `(atlas_cell_offsets_i, atlas_cell_offsets_j)` corners 
# of one pressure level, with Panama-barrier or invalid-neighbour fill-in.
@inline function level_interpolation(Atl, i, j, k, λ, φ, ξ, η, inside_panama)

    λr = Atl.λ[i]
    φr = Atl.φ[j]
    Δλ = Atl.λ[2] - Atl.λ[1]
    Δφ = Atl.φ[2] - Atl.φ[1]
    δi = atlas_cell_offsets_i
    δj = atlas_cell_offsets_j

    corners = ntuple(c -> Atl.saar[k, j + δj[c], i + δi[c]], Val(4))

    corners = if inside_panama
        apply_panama_barrier(corners, λ, φ, λr, φr, Δλ, Δφ)
    elseif abs(corners[1] + corners[2] + corners[3] + corners[4]) >= gsw_error_limit
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
    Sᴬ_from_Sᴾ(Sᴾ, p, λ, φ)

Return the TEOS-10 absolute salinity ``Sᴬ`` from practical salinity, sea pressure, longitude and
latitude, computed as

```math
Sᴬ = (Sₒ/35)\\, Sᴾ\\, (1 + r_{\\mathrm{SAAR}}) ,
```

where ``Sₒ = 35.16504`` g/kg is the standard ocean reference salinity and ``r_{\\mathrm{SAAR}}`` is
interpolated from the bundled global atlas (see [`saar`](@ref)). An analytic correction is applied
inside the Baltic Sea polygon (see [`baltic_absolute_salinity`](@ref)). Translation of `gsw_sa_from_sp`
of https://github.com/TEOS-10/GSW-C.

# Inputs
- `Sᴾ`: practical salinity, PSS-78  [unitless]
- `p` : sea pressure                [dbar]
- `λ` : longitude                   [°E]
- `φ` : latitude                    [°N]

# Output
- `Sᴬ`: absolute salinity [g/kg], `NaN` when the atlas cannot resolve the lookup point

# References
- McDougall, T. J., Jackett, D. R., Millero, F. J., Pawlowicz, R., Barker, P. M., 2012: A global
  algorithm for estimating absolute salinity. Ocean Science 8, 1123–1134.
"""
function Sᴬ_from_Sᴾ(Sᴾ, p, λ, φ)
    Sᴬ = baltic_absolute_salinity(Sᴾ, λ, φ)
    invalid(Sᴬ) || return oftype(float(Sᴾ), Sᴬ)

    rsaar = saar(float(p), float(λ), float(φ))
    invalid(rsaar) && return oftype(float(Sᴾ), NaN)

    return oftype(float(Sᴾ), uₚₛ * Sᴾ * (1 + rsaar))
end
