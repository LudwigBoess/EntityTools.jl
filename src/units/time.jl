"""

omega_p = sqrt(4*pi*n*q^2/m)
"""
function omega_i(unit::EntityUnits, mi::T=T(100), me::T=T(1)) where T <: Real
    return unit.larmor0 * mi / me 
end