# =============================================================================
#  Normalised-unit conversion factors for an Entity run
# =============================================================================
#
#  Entity's `[scales]` block is anchored to a *fiducial* population: number
#  density n0 = ppc0/V0, mass m0 = 1, charge q0 and field B0, with
#
#    omega_p0 = 1/skindepth0,   Omega_B0 = 1/larmor0 = B0,
#    sigma0   = (skindepth0/larmor0)^2
#
#  The plasma actually injected is generally *not* at n0: `ComputeNumInject`
#  fills a cell with ppc0*nd/2 particles of each species, so a two-species
#  injection with tot_number_density = 1 leaves every species at n_s = n0/2.
#  `n_frac` = n_s/n0 carries that factor and every density-dependent scale below
#  picks up the corresponding power of it. `n_frac = 1` (the default) reproduces
#  the fiducial scales Entity itself reports; pass `n_frac = 0.5` for the
#  standard two-species SRPIC injection.
#
#  Unit chain (c = 1, masses in m_e, temperatures in m_e c^2, mu = m_i/m_e,
#  fn = n_frac, T_i = temperature_ratio * T_e):
#
#    omega_pe  = sqrt(fn)/skindepth0        omega_pi = sqrt(fn/mu)/skindepth0
#    Omega_ce  = 1/larmor0                  Omega_ci = 1/(larmor0*mu)
#    d_e       = skindepth0/sqrt(fn)        d_i      = sqrt(mu/fn)*skindepth0
#    rho_e     = sqrt(T_e)*larmor0          rho_i    = sqrt(T_i*mu)*larmor0
#    lambda_De = sqrt(T_e/fn)*skindepth0
#    v_A       = sqrt(sigma0/(fn*mu_rho))   c_s      = sqrt(gamma*(T_e+T_i)/mu_rho)
#    beta_s    = 2*fn*T_s/sigma0
#
#  mu_rho is the mass density in units of n_i*m_e: mu for `mass_density = :ion`
#  (ions only, the shock-literature convention) or mu+1 for `:total` (the plasma
#  actually simulated). The two shift v_A and c_s together, so beta and c_s/v_A
#  are unaffected.
#
#  The chain assumes the SRPIC (full-PIC) normalisation with charges |q_s| = 1.
# =============================================================================

"""
    EntityUnits{T}

Conversion factors between Entity's normalised code units and the physical
plasma scales of a run.

Lengths, times and speeds are all in code units (c = 1, so a code time and a
code length are the same number). The `x_*`, `t_*` and `v_*` fields are
*multiplicative* conversion factors: multiply a code-unit value by them to get
the value in the unit named by the suffix, e.g.

    x_di   = positions in ion inertial lengths     x = x_code * u.x_di
    t_Ωci  = time in inverse ion cyclotron times   t = t_code * u.t_Omega_ci

Field output is already normalised by Entity itself (B and E in units of B0,
`N_s` in units of n0), so those need no conversion.

See the source header for the full unit chain, and
[`EntityUnits(::EntityData)`](@ref) / [`EntityUnits(::String)`](@ref) for the
constructors that pull the scales out of a run.

# Fields

## Entity `[scales]`
|Name | Meaning |
|:--- |:--- |
|`skindepth0::T` | fiducial skin depth c/ω_p0 [code length] |
|`larmor0::T` | fiducial Larmor radius c/Ω_B0 at v = c [code length] |
|`sigma0::T` | magnetisation (skindepth0/larmor0)² |
|`B0::T` | fiducial field strength 1/larmor0 |
|`omegaB0::T` | fiducial cyclotron frequency 1/larmor0 [1/code time] |
|`q0::T` | fiducial charge V0/(ppc0·skindepth0²) |
|`V0::T` | fiducial cell volume [code length^dim] |
|`n0::T` | fiducial number density ppc0/V0 |
|`ppc0::T` | fiducial particles per cell |

## Grid
|Name | Meaning |
|:--- |:--- |
|`dim::Int` | number of spatial dimensions |
|`dx::T` | fiducial cell size [code length] |

## Plasma composition
|Name | Meaning |
|:--- |:--- |
|`mass_ratio::T` | μ = m_i/m_e |
|`n_frac::T` | fn = n_s/n0 carried by one species |
|`mu_rho::T` | mass density in units of n_i·m_e (μ or μ+1) |
|`temperature::T` | T_e in units of m_e c² (`NaN` if unknown) |
|`temperature_ratio::T` | T_i/T_e |
|`gamma_ad::T` | adiabatic index behind c_s |

## Lengths [code length]
|Name | Meaning |
|:--- |:--- |
|`d_e::T` | electron skin depth c/ω_pe |
|`d_i::T` | ion inertial length c/ω_pi |
|`rho_e::T` | thermal electron Larmor radius v_the/Ω_ce |
|`rho_i::T` | thermal ion Larmor radius v_thi/Ω_ci |
|`lambda_De::T` | electron Debye length v_the/ω_pe |

## Frequencies [1/code time]
|Name | Meaning |
|:--- |:--- |
|`omega_pe::T` | electron plasma frequency |
|`omega_pi::T` | ion plasma frequency |
|`Omega_ce::T` | electron cyclotron frequency |
|`Omega_ci::T` | ion cyclotron frequency |

## Speeds [c]
|Name | Meaning |
|:--- |:--- |
|`v_A::T` | Alfvén speed |
|`c_s::T` | sound speed |
|`v_the::T` | electron thermal speed sqrt(T_e) |
|`v_thi::T` | ion thermal speed sqrt(T_i/μ) |

## Plasma beta
|Name | Meaning |
|:--- |:--- |
|`beta_e::T` | electron beta 2·fn·T_e/σ0 |
|`beta_i::T` | ion beta 2·fn·T_i/σ0 |
|`beta::T` | total beta β_e + β_i |

## Length conversion factors (multiply a code length)
|Name | Meaning |
|:--- |:--- |
|`x_cells::T` | position in cells |
|`x_de::T` | position in electron skin depths |
|`x_di::T` | position in ion inertial lengths |
|`x_rho_e::T` | position in electron Larmor radii |
|`x_rho_i::T` | position in ion Larmor radii |
|`x_lambda_De::T` | position in electron Debye lengths |

## Time conversion factors (multiply a code time)
|Name | Meaning |
|:--- |:--- |
|`t_omega_pe::T` | time in inverse electron plasma frequencies |
|`t_omega_pi::T` | time in inverse ion plasma frequencies |
|`t_Omega_ce::T` | time in inverse electron cyclotron times |
|`t_Omega_ci::T` | time in inverse ion cyclotron times |

## Speed conversion factors (multiply a code velocity, which is in c)
|Name | Meaning |
|:--- |:--- |
|`v_vA::T` | speed in Alfvén speeds (Alfvén Mach number) |
|`v_cs::T` | speed in sound speeds (sonic Mach number) |
"""
struct EntityUnits{T}

    # --- Entity [scales] ---
    skindepth0::T
    larmor0::T
    sigma0::T
    B0::T
    omegaB0::T
    q0::T
    V0::T
    n0::T
    ppc0::T

    # --- grid ---
    dim::Int
    dx::T

    # --- plasma composition ---
    mass_ratio::T
    n_frac::T
    mu_rho::T
    temperature::T
    temperature_ratio::T
    gamma_ad::T

    # --- lengths [code length] ---
    d_e::T
    d_i::T
    rho_e::T
    rho_i::T
    lambda_De::T

    # --- frequencies [1/code time] ---
    omega_pe::T
    omega_pi::T
    Omega_ce::T
    Omega_ci::T

    # --- speeds [c] ---
    v_A::T
    c_s::T
    v_the::T
    v_thi::T

    # --- plasma beta ---
    beta_e::T
    beta_i::T
    beta::T

    # --- length conversion factors ---
    x_cells::T
    x_de::T
    x_di::T
    x_rho_e::T
    x_rho_i::T
    x_lambda_De::T

    # --- time conversion factors ---
    t_omega_pe::T
    t_omega_pi::T
    t_Omega_ce::T
    t_Omega_ci::T

    # --- speed conversion factors ---
    v_vA::T
    v_cs::T
end

"""
    EntityUnits(skindepth0::Real, larmor0::Real; kwargs...)

Build the conversion factors from Entity's two fundamental scales.

Everything that does not depend on the density or the temperature (`sigma0`,
`B0`, `rho_e`, `rho_i`, `Omega_ce`, `Omega_ci`) follows from `skindepth0` and
`larmor0` alone. The remaining keywords describe the plasma that was actually
injected and the grid it lives on.

# Keyword arguments
- `dim = 1`               : number of spatial dimensions.
- `dx = 1`                : fiducial cell size [code length].
- `V0 = dx^dim`           : fiducial cell volume.
- `ppc0 = 1`              : fiducial particles per cell.
- `n0 = ppc0/V0`          : fiducial number density.
- `q0 = V0/(ppc0*skindepth0^2)` : fiducial charge.
- `mass_ratio = 1`        : μ = m_i/m_e.
- `n_frac = 1`            : fn = n_s/n0 for one species. Entity injects
  `ppc0*nd/2` particles per species, so a two-species run with
  `tot_number_density = 1` wants `n_frac = 0.5` — with `n_frac = 1` the derived
  scales are the fiducial ones Entity reports, and d_e, d_i, λ_De, v_A and β
  are then off by a factor √2 or 2.
- `mass_density = :ion`   : mass density behind v_A and c_s, `:ion` (ρ = n_i m_i)
  or `:total` (ρ = n_i m_i + n_e m_e).
- `temperature = NaN`     : T_e in units of m_e c². Left unset, every
  temperature-dependent field comes out `NaN` rather than silently zero.
- `temperature_ratio = 1` : T_i/T_e.
- `gamma_ad = 5/3`        : adiabatic index behind c_s.
"""
function EntityUnits(skindepth0::Real, larmor0::Real;
                     dim::Integer          = 1,
                     dx::Real              = one(skindepth0),
                     V0::Real              = dx^dim,
                     ppc0::Real            = one(skindepth0),
                     n0::Real              = ppc0 / V0,
                     q0::Real              = V0 / (ppc0 * skindepth0^2),
                     mass_ratio::Real      = 1,
                     n_frac::Real          = 1,
                     mass_density::Symbol  = :ion,
                     temperature::Real     = NaN,
                     temperature_ratio::Real = 1,
                     gamma_ad::Real        = 5 / 3)

    # ----------------------------- validation -----------------------------
    skindepth0 > 0 || error("skindepth0 must be positive (got $skindepth0)")
    larmor0 > 0    || error("larmor0 must be positive (got $larmor0)")
    dim in (1, 2, 3) || error("dim must be 1, 2 or 3 (got $dim)")
    mass_ratio > 0 || error("mass_ratio must be positive (got $mass_ratio)")
    n_frac > 0     || error("n_frac must be positive (got $n_frac)")
    temperature_ratio > 0 || error("temperature_ratio must be positive (got $temperature_ratio)")
    isnan(temperature) || temperature ≥ 0 ||
        error("temperature must be non-negative (got $temperature)")
    mass_density in (:ion, :total) ||
        error("mass_density must be :ion or :total (got $mass_density)")

    # analysis is always done in at least double precision, whatever the run wrote
    TT = float(promote_type(Float64, typeof(skindepth0), typeof(larmor0),
                            typeof(dx), typeof(V0), typeof(n0), typeof(q0)))
    d0 = TT(skindepth0)
    r0 = TT(larmor0)
    μ  = TT(mass_ratio)
    fn = TT(n_frac)
    τ  = TT(temperature_ratio)
    γ  = TT(gamma_ad)
    Te = TT(temperature)
    Ti = τ * Te
    # mass density in units of n_i*m_e: ions only, or the full plasma
    μρ = mass_density === :ion ? μ : μ + one(TT)

    # ------------------------------ scales --------------------------------
    σ0 = (d0 / r0)^2
    B0 = one(TT) / r0

    # frequencies [1/code time]
    ω_pe = sqrt(fn) / d0
    ω_pi = sqrt(fn / μ) / d0
    Ω_ce = B0
    Ω_ci = B0 / μ

    # thermal speeds [c]
    v_the = sqrt(Te)
    v_thi = sqrt(Ti / μ)

    # lengths [code length]
    d_e  = one(TT) / ω_pe
    d_i  = one(TT) / ω_pi
    ρ_e  = v_the / Ω_ce
    ρ_i  = v_thi / Ω_ci
    λ_De = v_the / ω_pe

    # bulk speeds [c]
    v_A = sqrt(σ0 / (fn * μρ))
    c_s = sqrt(γ * (Te + Ti) / μρ)

    # plasma beta
    β_e = 2 * fn * Te / σ0
    β_i = 2 * fn * Ti / σ0

    return EntityUnits{TT}(
        d0, r0, σ0, B0, B0, TT(q0), TT(V0), TT(n0), TT(ppc0),
        Int(dim), TT(dx),
        μ, fn, μρ, Te, τ, γ,
        d_e, d_i, ρ_e, ρ_i, λ_De,
        ω_pe, ω_pi, Ω_ce, Ω_ci,
        v_A, c_s, v_the, v_thi,
        β_e, β_i, β_e + β_i,
        one(TT) / TT(dx), one(TT) / d_e, one(TT) / d_i,
        one(TT) / ρ_e, one(TT) / ρ_i, one(TT) / λ_De,
        ω_pe, ω_pi, Ω_ce, Ω_ci,
        one(TT) / v_A, one(TT) / c_s)
end

function Base.show(io::IO, ::MIME"text/plain", u::EntityUnits{T}) where {T}
    @printf(io, "EntityUnits{%s}  dim = %d   μ = %.4g   n_frac = %.4g   T_e = %.4g   T_i/T_e = %.4g\n",
            string(T), u.dim, u.mass_ratio, u.n_frac, u.temperature, u.temperature_ratio)
    @printf(io, "  scales   skindepth0 = %-10.4g larmor0 = %-10.4g sigma0 = %-10.4g B0 = %.4g\n",
            u.skindepth0, u.larmor0, u.sigma0, u.B0)
    @printf(io, "           dx = %-18.4g V0 = %-13.4g n0 = %-13.4g q0 = %.4g\n",
            u.dx, u.V0, u.n0, u.q0)
    @printf(io, "  lengths  d_e = %.4g (%.3g dx)   d_i = %.4g (%.3g dx)   rho_e = %.4g   rho_i = %.4g   lambda_De = %.4g\n",
            u.d_e, u.d_e * u.x_cells, u.d_i, u.d_i * u.x_cells, u.rho_e, u.rho_i, u.lambda_De)
    @printf(io, "  times    1/omega_pe = %.4g   1/omega_pi = %.4g   1/Omega_ce = %.4g   1/Omega_ci = %.4g\n",
            1 / u.omega_pe, 1 / u.omega_pi, 1 / u.Omega_ce, 1 / u.Omega_ci)
    @printf(io, "  speeds   v_A = %.4g c   c_s = %.4g c   c/v_A = %.4g   c_s/v_A = %.4g\n",
            u.v_A, u.c_s, 1 / u.v_A, u.c_s / u.v_A)
    @printf(io, "  beta     beta_e = %.4g   beta_i = %.4g   beta = %.4g",
            u.beta_e, u.beta_i, u.beta)
    return nothing
end

Base.show(io::IO, u::EntityUnits{T}) where {T} =
    @printf(io, "EntityUnits{%s}(skindepth0 = %.4g, larmor0 = %.4g, mass_ratio = %.4g, n_frac = %.4g)",
            string(T), u.skindepth0, u.larmor0, u.mass_ratio, u.n_frac)

# =============================================================================
#  The same quantities in cgs-Gaussian units
# =============================================================================
#
#  A PIC run in normalised units is a family of physical systems, not one: the
#  code fixes c = 1, the mass unit to the lightest species' mass and the charge
#  unit to |e|, which leaves exactly *one* dimensional freedom — the length of a
#  code length in cm. Fixing it (equivalently: naming a density, or the field
#  strength of B0) fixes every other cgs value.
#
#  Writing L0 [cm] for the code length, T0 = L0/c [s], M0 [g] for the mass unit
#  and Q0 [esu] for the charge unit, the derived units are
#
#    field   B  = M0 c / (Q0 T0)   [G]     (E in statvolt/cm uses the same factor)
#    energy  E  = M0 c^2           [erg]
#    density n0 = omega_p0^2 M0 / (4 pi Q0^2), omega_p0 = (1/skindepth0)/T0
#
#  Dimensionless quantities (sigma0, the betas, the mass ratio) are the same
#  numbers in both systems, and every ratio of two lengths, times or speeds is
#  preserved — the cgs struct only ever rescales.
# =============================================================================

"""
    EntityUnitsCGS{T}

The quantities of [`EntityUnits`](@ref) in cgs-Gaussian units: lengths in cm,
times in s, speeds in cm/s, frequencies in rad/s, fields in gauss, charges in
esu, temperatures in kelvin.

Field names and order match [`EntityUnits`](@ref) up to the trailing block of
code-to-cgs conversion factors, so code written against one struct also reads
the other. The `x_*`, `t_*` and `v_*` fields stay *multiplicative* factors, now
applied to a cgs value: `x_di` [1/cm] turns a length in cm into a length in ion
inertial lengths, `t_Omega_ci` [rad/s] turns a time in seconds into Ω_ci·t, and
`v_vA` [s/cm] turns a speed in cm/s into an Alfvén Mach number.

To go the other way — from a code-unit array to cgs — multiply by the unit
fields: `x_cm = x_code * u.length_unit`, `B_gauss = B_code * u.Bfield_unit`.

See [`EntityUnitsCGS(::EntityUnits)`](@ref) for the anchor that sets the scale.

# Fields

## Entity `[scales]`
|Name | Meaning |
|:--- |:--- |
|`skindepth0::T` | fiducial skin depth c/ω_p0 [cm] |
|`larmor0::T` | fiducial Larmor radius c/Ω_B0 at v = c [cm] |
|`sigma0::T` | magnetisation (skindepth0/larmor0)², dimensionless |
|`B0::T` | fiducial field strength [G] |
|`omegaB0::T` | fiducial cyclotron frequency q·B0/(m·c) [rad/s] |
|`q0::T` | charge of one fiducial macroparticle [esu] |
|`V0::T` | fiducial cell volume [cm³] |
|`n0::T` | fiducial number density [cm⁻³] |
|`ppc0::T` | fiducial particles per cell |

`V0` is the *three-dimensional* cell volume: in `dim < 3` the code's unit
thickness along the missing directions maps to `length_unit`, so
`V0 = V0_code · length_unit³`. The macroparticle weight `w = n0·V0/ppc0` is
generally not 1, so unlike in code units `n0 ≠ ppc0/V0`; `q0 = w·charge_unit` is
the charge that weight carries.

## Grid
|Name | Meaning |
|:--- |:--- |
|`dim::Int` | number of spatial dimensions |
|`dx::T` | fiducial cell size [cm] |

## Plasma composition
|Name | Meaning |
|:--- |:--- |
|`mass_ratio::T` | μ = m_i/m_e |
|`n_frac::T` | fn = n_s/n0 carried by one species |
|`mu_rho::T` | mass density in units of n_i·m_e (μ or μ+1) |
|`temperature::T` | T_e [K] (`NaN` if unknown) |
|`temperature_ratio::T` | T_i/T_e |
|`gamma_ad::T` | adiabatic index behind c_s |

## Lengths [cm]
|Name | Meaning |
|:--- |:--- |
|`d_e::T` | electron skin depth c/ω_pe |
|`d_i::T` | ion inertial length c/ω_pi |
|`rho_e::T` | thermal electron Larmor radius v_the/Ω_ce |
|`rho_i::T` | thermal ion Larmor radius v_thi/Ω_ci |
|`lambda_De::T` | electron Debye length v_the/ω_pe |

## Frequencies [rad/s]
|Name | Meaning |
|:--- |:--- |
|`omega_pe::T` | electron plasma frequency sqrt(4π n_e q²/m_e) |
|`omega_pi::T` | ion plasma frequency |
|`Omega_ce::T` | electron cyclotron frequency q·B0/(m_e·c) |
|`Omega_ci::T` | ion cyclotron frequency |

## Speeds [cm/s]
|Name | Meaning |
|:--- |:--- |
|`v_A::T` | Alfvén speed B0/sqrt(4π ρ) |
|`c_s::T` | sound speed |
|`v_the::T` | electron thermal speed sqrt(kB T_e/m_e) |
|`v_thi::T` | ion thermal speed sqrt(kB T_i/m_i) |

## Plasma beta (dimensionless)
|Name | Meaning |
|:--- |:--- |
|`beta_e::T` | electron beta 8π n_e kB T_e/B0² |
|`beta_i::T` | ion beta |
|`beta::T` | total beta β_e + β_i |

## Length conversion factors [1/cm] (multiply a length in cm)
|Name | Meaning |
|:--- |:--- |
|`x_cells::T` | position in cells |
|`x_de::T` | position in electron skin depths |
|`x_di::T` | position in ion inertial lengths |
|`x_rho_e::T` | position in electron Larmor radii |
|`x_rho_i::T` | position in ion Larmor radii |
|`x_lambda_De::T` | position in electron Debye lengths |

## Time conversion factors [rad/s] (multiply a time in s)
|Name | Meaning |
|:--- |:--- |
|`t_omega_pe::T` | time in inverse electron plasma frequencies |
|`t_omega_pi::T` | time in inverse ion plasma frequencies |
|`t_Omega_ce::T` | time in inverse electron cyclotron times |
|`t_Omega_ci::T` | time in inverse ion cyclotron times |

## Speed conversion factors [s/cm] (multiply a speed in cm/s)
|Name | Meaning |
|:--- |:--- |
|`v_vA::T` | speed in Alfvén speeds (Alfvén Mach number) |
|`v_cs::T` | speed in sound speeds (sonic Mach number) |

## Code-to-cgs conversion factors (multiply a code-unit value)
|Name | Meaning |
|:--- |:--- |
|`length_unit::T` | cm per code length |
|`time_unit::T` | s per code time, `length_unit`/c |
|`speed_unit::T` | cm/s per code speed, i.e. c |
|`mass_unit::T` | g per code mass |
|`charge_unit::T` | esu per code charge |
|`Bfield_unit::T` | G per code field (statvolt/cm for E) |
|`energy_unit::T` | erg per code energy, `mass_unit`·c² |
"""
struct EntityUnitsCGS{T}

    # --- Entity [scales] ---
    skindepth0::T
    larmor0::T
    sigma0::T
    B0::T
    omegaB0::T
    q0::T
    V0::T
    n0::T
    ppc0::T

    # --- grid ---
    dim::Int
    dx::T

    # --- plasma composition ---
    mass_ratio::T
    n_frac::T
    mu_rho::T
    temperature::T
    temperature_ratio::T
    gamma_ad::T

    # --- lengths [cm] ---
    d_e::T
    d_i::T
    rho_e::T
    rho_i::T
    lambda_De::T

    # --- frequencies [rad/s] ---
    omega_pe::T
    omega_pi::T
    Omega_ce::T
    Omega_ci::T

    # --- speeds [cm/s] ---
    v_A::T
    c_s::T
    v_the::T
    v_thi::T

    # --- plasma beta ---
    beta_e::T
    beta_i::T
    beta::T

    # --- length conversion factors [1/cm] ---
    x_cells::T
    x_de::T
    x_di::T
    x_rho_e::T
    x_rho_i::T
    x_lambda_De::T

    # --- time conversion factors [rad/s] ---
    t_omega_pe::T
    t_omega_pi::T
    t_Omega_ce::T
    t_Omega_ci::T

    # --- speed conversion factors [s/cm] ---
    v_vA::T
    v_cs::T

    # --- code -> cgs conversion factors ---
    length_unit::T
    time_unit::T
    speed_unit::T
    mass_unit::T
    charge_unit::T
    Bfield_unit::T
    energy_unit::T
end

"""
    EntityUnitsCGS(u::EntityUnits; <one anchor>, kwargs...)

Attach physical units to the normalised scales of `u`.

Exactly one anchor keyword must be given; it sets the length of a code length in
cm, and everything else follows:

- `n_e`         : number density of the light (unit-mass) species [cm⁻³], i.e.
  `n_frac`·n0 — the density that sets ω_pe.
- `n0`          : fiducial number density [cm⁻³], the total density of a
  two-species injection with `n_frac = 0.5`.
- `B0`          : field strength of Entity's fiducial B0 [G].
- `length_unit` : cm per code length, if the mapping is already known.

# Keyword arguments
- `mass_unit = EntityTools.m_e_cgs`   : grams per code mass unit, i.e. the mass
  of the unit-mass species. Set it to `EntityTools.m_p_cgs` to read the light
  species as protons (a hybrid run), or to any value for a reduced-mass setup.
- `charge_unit = EntityTools.q_e_cgs` : esu per code charge unit, the charge of a
  `q = 1` species.

```jldoctest
julia> u  = EntityUnits(1.0, 10.0; mass_ratio = 100.0, n_frac = 0.5, temperature = 1e-3);

julia> uc = EntityUnitsCGS(u; n_e = 5.0);   # 5 electrons per cm^3

julia> isapprox(uc.omega_pe, 5.64e4 * sqrt(5.0); rtol = 1e-3)   # NRL formulary
true
```
"""
function EntityUnitsCGS(u::EntityUnits{T};
                        n_e                = nothing,
                        n0                 = nothing,
                        B0                 = nothing,
                        length_unit        = nothing,
                        mass_unit::Real    = m_e_cgs,
                        charge_unit::Real  = q_e_cgs) where {T}

    # ----------------------------- validation -----------------------------
    anchors = (; n_e, n0, B0, length_unit)
    given   = [k for (k, v) in pairs(anchors) if v !== nothing]
    length(given) == 1 ||
        error("give exactly one of n_e, n0, B0, length_unit (got " *
              (isempty(given) ? "none" : join(given, ", ")) * ")")
    val = anchors[only(given)]
    val > 0 || error("$(only(given)) must be positive (got $val)")
    mass_unit > 0   || error("mass_unit must be positive (got $mass_unit)")
    charge_unit > 0 || error("charge_unit must be positive (got $charge_unit)")

    TT = float(promote_type(Float64, T, typeof(mass_unit), typeof(charge_unit), typeof(val)))
    M0 = TT(mass_unit)
    Q0 = TT(charge_unit)
    cc = TT(c_cgs)

    # --------------------- the one dimensional freedom ---------------------
    # every anchor reduces to L0, the number of cm in a code length
    L0 = if length_unit !== nothing
        TT(length_unit)
    elseif n_e !== nothing || n0 !== nothing
        # omega_p0 at the fiducial density fixes the code time unit
        nfid  = n_e !== nothing ? TT(n_e) / u.n_frac : TT(n0)
        ω_p0  = sqrt(4 * TT(π) * nfid * Q0^2 / M0)
        cc * (one(TT) / u.skindepth0) / ω_p0
    else
        # Omega_ce = q B0 / (m c) does the same job
        Ω_ce = Q0 * TT(B0) / (M0 * cc)
        cc * u.omegaB0 / Ω_ce
    end

    # ---------------------------- derived units ---------------------------
    T0 = L0 / cc                    # s per code time
    Bu = M0 * cc / (Q0 * T0)        # G per code field
    Eu = M0 * cc^2                  # erg per code energy

    # fiducial density and the macroparticle weight it implies
    n0_cgs = (one(TT) / u.skindepth0 / T0)^2 * M0 / (4 * TT(π) * Q0^2)
    V0_cgs = TT(u.V0) * L0^3        # 3D cell volume, unit thickness in dim < 3
    q0_cgs = n0_cgs * V0_cgs / u.ppc0 * Q0

    return EntityUnitsCGS{TT}(
        u.skindepth0 * L0, u.larmor0 * L0, TT(u.sigma0), u.B0 * Bu, u.omegaB0 / T0,
        q0_cgs, V0_cgs, n0_cgs, TT(u.ppc0),
        u.dim, u.dx * L0,
        TT(u.mass_ratio), TT(u.n_frac), TT(u.mu_rho),
        u.temperature * Eu / TT(k_B_cgs), TT(u.temperature_ratio), TT(u.gamma_ad),
        u.d_e * L0, u.d_i * L0, u.rho_e * L0, u.rho_i * L0, u.lambda_De * L0,
        u.omega_pe / T0, u.omega_pi / T0, u.Omega_ce / T0, u.Omega_ci / T0,
        u.v_A * cc, u.c_s * cc, u.v_the * cc, u.v_thi * cc,
        TT(u.beta_e), TT(u.beta_i), TT(u.beta),
        u.x_cells / L0, u.x_de / L0, u.x_di / L0,
        u.x_rho_e / L0, u.x_rho_i / L0, u.x_lambda_De / L0,
        u.t_omega_pe / T0, u.t_omega_pi / T0, u.t_Omega_ce / T0, u.t_Omega_ci / T0,
        u.v_vA / cc, u.v_cs / cc,
        L0, T0, cc, M0, Q0, Bu, Eu)
end

function Base.show(io::IO, ::MIME"text/plain", u::EntityUnitsCGS{T}) where {T}
    @printf(io, "EntityUnitsCGS{%s}  dim = %d   μ = %.4g   n_frac = %.4g   T_e = %.4g K (%.4g eV)   T_i/T_e = %.4g\n",
            string(T), u.dim, u.mass_ratio, u.n_frac, u.temperature,
            u.temperature * k_B_cgs / eV_cgs, u.temperature_ratio)
    @printf(io, "  scales   skindepth0 = %-10.4g cm   larmor0 = %-10.4g cm   sigma0 = %-10.4g B0 = %.4g G\n",
            u.skindepth0, u.larmor0, u.sigma0, u.B0)
    @printf(io, "           dx = %-18.4g cm   V0 = %-13.4g cm^3   n0 = %-10.4g cm^-3   q0 = %.4g esu\n",
            u.dx, u.V0, u.n0, u.q0)
    @printf(io, "  units    1 code length = %.4g cm   1 code time = %.4g s   1 code B = %.4g G\n",
            u.length_unit, u.time_unit, u.Bfield_unit)
    @printf(io, "  lengths  d_e = %.4g cm (%.3g dx)   d_i = %.4g cm (%.3g dx)   rho_e = %.4g cm   rho_i = %.4g cm   lambda_De = %.4g cm\n",
            u.d_e, u.d_e * u.x_cells, u.d_i, u.d_i * u.x_cells, u.rho_e, u.rho_i, u.lambda_De)
    @printf(io, "  times    1/omega_pe = %.4g s   1/omega_pi = %.4g s   1/Omega_ce = %.4g s   1/Omega_ci = %.4g s\n",
            1 / u.omega_pe, 1 / u.omega_pi, 1 / u.Omega_ce, 1 / u.Omega_ci)
    @printf(io, "  speeds   v_A = %.4g cm/s   c_s = %.4g cm/s   c/v_A = %.4g   c_s/v_A = %.4g\n",
            u.v_A, u.c_s, u.speed_unit / u.v_A, u.c_s / u.v_A)
    @printf(io, "  beta     beta_e = %.4g   beta_i = %.4g   beta = %.4g",
            u.beta_e, u.beta_i, u.beta)
    return nothing
end

Base.show(io::IO, u::EntityUnitsCGS{T}) where {T} =
    @printf(io, "EntityUnitsCGS{%s}(n0 = %.4g cm^-3, B0 = %.4g G, length_unit = %.4g cm, mass_ratio = %.4g)",
            string(T), u.n0, u.B0, u.length_unit, u.mass_ratio)
