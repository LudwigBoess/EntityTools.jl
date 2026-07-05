# Plasma normalisations and heat-flux conversions.
#
# Entity normalises EM fields to B0 = 1/larmor0 and lengths to the fiducial cell
# size; the helpers below turn the stored `EntityUnits` (skindepth0, larmor0,
# sigma0, B0, ...) into the dimensional quantities used in the whistler /
# heat-conduction problem of Komarov et al. (2018, https://arxiv.org/abs/1711.11462).
#
# Conventions (matching the conduction_komarov pgen):
#   - temperatures T are in units of kB T / (m c^2);
#   - the non-relativistic thermal speed is v_th = sqrt(T/m)  (so v_th1 = sqrt(T_hot));
#   - the electron Larmor radius rho_e = v_th / Omega_e = sqrt(T) * larmor0;
#   - with larmor0 = 1/v_th1, the hot-side rho_e1 = 1 code length.

"""
    omega_pe(u::EntityUnits)

Electron plasma frequency in inverse code-time units, ω_pe = c / d_e = 1 / skindepth0.
"""
omega_pe(u::EntityUnits) = 1 / u.skindepth0

"""
    omega_ce(u::EntityUnits)

Electron cyclotron frequency Ω_e = |q| B0 / m = B0 = 1 / larmor0 (≡ `u.omegaB0`).
Multiply a code-unit time by this to get Ω_e t.
"""
omega_ce(u::EntityUnits) = u.omegaB0

"""
    v_th(T; mass = 1)

Non-relativistic thermal speed (in units of c) for temperature `T` [kB T / (m c^2)]:
v_th = sqrt(T / mass). Komarov et al. use v_th1 = sqrt(T_hot).
"""
v_th(T::Real; mass::Real = 1.0) = sqrt(T / mass)

"""
    rho_e(u::EntityUnits, T; mass = 1)

Electron Larmor radius at temperature `T` in code-length units,
ρ_e = v_th / Ω_e = sqrt(T / mass) * larmor0.
"""
rho_e(u::EntityUnits, T::Real; mass::Real = 1.0) = v_th(T; mass) * u.larmor0

"""
    betae(u::EntityUnits, T)

Electron plasma beta at temperature `T`,
β_e = 2 n T / B0² = 2 T / σ0  (n = n0 = 1, B0 = 1 / larmor0, σ0 = (skindepth0/larmor0)²).

For the conduction_komarov runs this is the hot-side β evaluated with T = T_hot.
"""
betae(u::EntityUnits, T::Real) = 2 * T / u.sigma0

"""
    q_free(T; n = 1, mass = 1, prefactor = 1)

Free-streaming reference heat flux q0 = prefactor · n · m · v_th³ = prefactor · n · T^{3/2}/√m.

This is the per-n0 normalisation that matches the Entity `Qx`/`Qy`/`Qz` custom output,
which already carries the 1/n0 factor and the ½ m of the third moment ½ m ⟨v_∥ v²⟩.
Komarov et al. (2018) normalise the measured parallel heat flux q_∥ by this q0 and find
q_∥/q0 ≈ 1.5 / β_e. An O(1) `prefactor` (e.g. √(2/π) ≈ 0.8 or ½) can be applied to match
a particular free-streaming convention.
"""
q_free(T::Real; n::Real = 1.0, mass::Real = 1.0, prefactor::Real = 1.0) =
    prefactor * n * mass * v_th(T; mass)^3
