"""
    vA(unit::EntityUnits)

Alfvén speed in normalized units.
"""
function vA(unit::EntityUnits)
    return 1 / sqrt(unit.sigma0)
end