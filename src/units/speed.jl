"""
    vA(unit::EntityUnits)

Fiducial Alfvén speed in units of c, v_A = Ω_B0/ω_p0 = skindepth0/larmor0 =
sqrt(`sigma0`) — the magnetisation σ0 is v_A²/c², so this is the speed for one
particle mass (m = m0 = 1) at the fiducial density n0.

For the Alfvén speed of the plasma that was actually injected, which carries the
ion mass and the injected density, use the `v_A` field of
[`EntityUnits`](@ref): `u.v_A = sqrt(sigma0/(n_frac*mu_rho))`. The two agree
when `mass_ratio = 1` and `n_frac = 1`.
"""
function vA(unit::EntityUnits)
    return sqrt(unit.sigma0)
end
