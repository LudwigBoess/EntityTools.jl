"""
    rLi(unit::EntityUnits, mi, me)

Larmor radius of a particle with mass mi in normalized units.
"""
function rLi(unit::EntityUnits, mi::T, me::T=T(1)) where T <: Real
    return sqrt(2/unit.sigma0)*sqrt(mi/me)*unit.skindepth0
end