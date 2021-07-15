# SeawaterPolynomials.jl

This package provides approximations to the Boussinesq equation of state for seawater expressed as polynomial functions of conservative temperature, absolute salinity, and geopotential height. 

Computationally efficient polynomial approximations to the [Boussinesq seawater](https://doi.org/10.1175/2009JPO4294.1) equation of state are crucial components of ocean modeling software.

## The Seawater Boussinesq approximation

In the seawater Boussinesq approximation, the density of seawater is expanded around a constant reference value, ``ρᵣ``,

```math
ρ = ρ_r + ρ'(Θ, Sᴬ, Z) ,
```

where the anomaly, ``ρ'``, is a function of conservative temperature ``Θ``, absolute salinity ``Sᴬ``, and geopotential height ``Z``.
One choice for ``ρ_r`` is the average density at the surface of the world ocean, ``ρ_r = 1024.6 \, \text{kg} \, \text{m}^{-3}``, according to [Roquet et al. (2015)](https://doi.org/10.1175/JPO-D-15-0080.1).

## The TEOS-10 standard

The [Thermodynamic Equation of SeaWater (TEOS-10)](http://www.teos-10.org) is a Gibbs function formulation of seawater themodynamics.

The error between the polynomials implemented in this package and the TEOS-10 is minimized for current 'climatological' ocean distributions of conservative temperature and absolute salinity. For more information, see

* [Roquet et al., "Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard", Ocean Modelling (2015)](https://doi.org/10.1016/j.ocemod.2015.04.002)

* [Roquet et al., "Defining a Simplified Yet “Realistic” Equation of State for Seawater", Journal of Physical Oceanography (2015)](https://doi.org/10.1175/JPO-D-15-0080.1)

## Related packages

* [GibbsSeaWater.jl](https://github.com/TEOS-10/GibbsSeaWater.jl)

# References

* Roquet et al., "[Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard](https://doi.org/10.1016/j.ocemod.2015.04.002)", _Ocean Modelling_ (2015).

* Roquet et al., "[Defining a Simplified Yet “Realistic” Equation of State for Seawater](https://doi.org/10.1175/JPO-D-15-0080.1)", _Journal of Physical Oceanography_ (2015).

* Young, W. R., "[Dynamic Enthalpy, Conservative Temperature, and the Seawater Boussinesq Approximation](https://doi.org/10.1175/2009JPO4294.1)", _Journal of Physical Oceanography_ (2010)

