"""
    EntityUnits{T <: Real}

A struct to hold the normalized units for a plasma simulation. 
The fields are:
- `skindepth0`: The skin depth in normalized units.
- `larmor0`: The Larmor radius in normalized units.
- `sigma0`: The magnetization parameter, calculated as (skindepth0/larmor0)^2.
- `B0`: The magnetic field strength in normalized units, calculated as 1/larmor0.
- `omegaB0`: The cyclotron frequency in normalized units, calculated as B0.
- `q0`: The charge normalization factor, calculated as 1/(n0 * skindepth0^2).
- `V0`: The volume normalization factor, calculated as dx0^size(extent, 1).
- `n0`: The number density normalization factor, calculated as ppc0 / V0.
"""
struct EntityUnits{T}

    skindepth0::T
    larmor0::T
    sigma0::T
    B0::T
    omegaB0::T
    q0::T
    V0::T
    n0::T

    """
        EntityUnits(skindepth0::T, larmor0::T, sigma0::T, B0::T, omegaB0::T, q0::T, V0::T, n0::T) where T <: Real

    Constructs an `EntityUnits` struct from the provided parameters.
    """
    EntityUnits(skindepth0::T, larmor0::T, sigma0::T, B0::T, omegaB0::T, q0::T, V0::T, n0::T) where T <: Real = new{T}(skindepth0, larmor0, sigma0, B0, omegaB0, q0, V0, n0)

    """
        EntityUnits(skindepth0::T, larmor0::T) where T <: Real

    Constructs an `EntityUnits` struct from the provided skin depth and Larmor radius. 
    The other fields are calculated based on these two parameters.
    Charges and densities are set to 1.
    """
    function EntityUnits(skindepth0::T, larmor0::T) where T <: Real
        sigma0 = (skindepth0/larmor0)^2
        B0 = 1/larmor0
        new{T}(skindepth0, larmor0, sigma0, B0, B0, T(1), T(1), T(1))
    end

    """
        EntityUnits(skindepth0::T, larmor0::T, ppc0::T, extent::Vector{T}, resolution::Integer, metric::String) where T <: Real

    Constructs an `EntityUnits` struct from the provided parameters.
    - `skindepth0`: The skin depth in normalized units.
    - `larmor0`: The Larmor radius in normalized units.
    - `ppc0`: The number of particles per cell.
    - `extent`: The spatial extent of the simulation domain.
    - `resolution`: The number of grid points in each dimension.
    - `metric`: The metric of the simulation (currently only "minkowski" is supported).

    """
    function EntityUnits(skindepth0::T, larmor0::T, ppc0::T, extent::Vector{Any}, resolution::Vector{<:Integer}, metric::String) where T <: Real
        
        if metric != "minkowski"
            error("Only minkowski metric is currently supported.")
        end

        dx0 = (extent[1][2] - extent[1][1]) / resolution[1]
        V0 = dx0^size(extent, 1)
        n0 = ppc0 / V0
        q0 = 1/ (n0 * skindepth0^2)
        sigma0 = (skindepth0/larmor0)^2
        B0 = 1/larmor0

        new{T}(skindepth0, larmor0, sigma0, B0, B0, q0, V0, n0)
    end

end

