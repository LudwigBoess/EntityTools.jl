# =============================================================================
#  Collisionless-shock setup generator for Entity (SRPIC / full-PIC mode)
# =============================================================================
#
#  Writes a ready-to-run Entity parameter file for a magnetised collisionless
#  shock, choosing the normalised scales so that user-supplied *physical*
#  numbers (mass ratio, temperature, sonic Mach number, plasma beta, runtime in
#  ion cyclotron times) come out correctly.
#
#  Three problem generators are supported via the `pgen` keyword:
#
#    :downstream    ->  LMBgens/shock_downstream    (key `drift_ux`)
#                       simulation frame = downstream rest frame; upstream is
#                       injected from the right and reflects off the left wall,
#                       the shock propagates to the right in a *fixed* box.
#
#    :upstream      ->  LMBgens/shock_upstream      (key `piston_velocity`)
#                       simulation frame = upstream rest frame; a piston pushes
#                       into plasma at rest and a *moving window* follows it.
#
#    :upstream_fill ->  LMBgens/shock_upstream_fill (adds `fill_fraction`,
#                       `wall_velocity`) as :upstream, plus the option to fill
#                       only the left `fill_fraction` of the box at t = 0 and
#                       close it off with a conducting wall that runs ahead of
#                       the shock at `wall_velocity`, refilling behind itself.
#
#  In all cases the injection/piston speed equals the velocity jump across the
#  shock, Δv, and the shock recedes from the box reference frame at Δv/(r-1) —
#  for the moving-window pgens because the window itself tracks the piston at Δv.
#
#  --------------------------------------------------------------------------
#  Unit chain (Entity SRPIC, c = 1)
#
#  Entity's [scales] are anchored to a *fiducial* population: number density n0
#  (ppc0 particles per fiducial cell), mass m0 = 1, charge q0 = 1, field B0:
#
#    ω_p0 = 1/skindepth0,   ω_B0 = 1/larmor0 = B0,   σ0 = (skindepth0/larmor0)^2
#
#  The plasma actually injected is not at n0. The pgens call
#  InjectUniformMaxwellians with tot_number_density = 1, and ComputeNumInject
#  puts ppc0/2 particles of *each* species in a cell, so both the electrons and
#  the ions sit at n_s = fn·n0 with fn = 1/2. Every density-dependent scale
#  therefore differs from the fiducial one by a power of fn:
#
#    ω_ps^2 = ω_p0^2 · fn · q_s^2/m_s
#    d_e  = skindepth0/sqrt(fn)      = sqrt(2)·skindepth0   (electron skin depth)
#    d_i  = sqrt(μ/fn)·skindepth0    = sqrt(2μ)·skindepth0  (ion inertial length)
#    v_A  = sqrt(σ0/(fn·μ_ρ))   [c]  (μ_ρ = μ for :ion mass density, μ+1 for :total)
#    β_e  = 2·fn·T_e/σ0,   β_i = 2·fn·T_i/σ0,   β_total = β_e + β_i
#    λ_De = sqrt(T_e/fn)·skindepth0
#
#  Density-free (both P and ρ scale with n, so fn cancels):
#
#    c_s  = sqrt(γ(T_e+T_i)/μ_ρ)  [c]
#    ρ_i  = larmor0·sqrt(T_i·μ)   [code length]  (thermal ion Larmor radius)
#    Ω_ci^-1 = larmor0·μ          [code time]    (c = 1 → time = length)
#
#  Identities locked by test/shock_setup.jl:
#
#    v_A = Ω_ci·d_i (for :ion),   d_i/d_e = sqrt(μ),   λ_De/d_e = sqrt(T_e)
#    ρ_i/d_i = sqrt(β_i/2),       c_s/v_A = sqrt(γ·β_total/2)
#
#  `temperature` is kB T_e/(m_e c^2) — the value written to `setup.temperature`,
#  which the pgen passes to InjectUniformMaxwellians (it divides by the species
#  mass), so T_e = temperature·m_e c^2 and T_i = temperature_ratio·T_e.
#  (Equivalently kB T_e/(m_i c^2) = temperature/μ.)
#
#  `mach_sonic` is the *flow* sonic Mach number, M_s = Δv/c_s (so Δv = M_s·c_s).
#  This matches the reference shock_3D/beta50 toml exactly. The corresponding
#  upstream-frame shock-Mach numbers are reported in the summary; if a setup is
#  quoted with the shock-frame Mach number instead, scale `mach_sonic` by the
#  ratio printed there.
# =============================================================================

# Fraction of the fiducial density n0 carried by each species. The shock pgens
# inject tot_number_density = 1 and Entity's ComputeNumInject fills a cell with
# ppc0·nd/2 particles per species, so n_e = n_i = n0/2. Changing this requires
# changing the injection call in the pgen.
const FN_PER_SPECIES = 0.5

"""
    write_shock_setup(path; mass_ratio, temperature, mach_sonic, beta, runtime, kwargs...)

Generate an Entity (full-PIC / SRPIC) collisionless-shock parameter file at
`path` and return a `NamedTuple` of the derived physical quantities.

# Required (keyword) arguments
- `mass_ratio`  : ion/electron mass ratio μ = m_i/m_e.
- `temperature` : kB T_e / (m_e c²) (electron-rest-mass units); T_i = `temperature_ratio`·T_e.
- `mach_sonic`  : flow sonic Mach number M_s = Δv/c_s.
- `beta`        : plasma beta (see `beta_definition`).
- `runtime`     : total run time in ion cyclotron times Ω_ci⁻¹.

# Geometry / discretisation
- `dim = 2`                 : spatial dimensionality (1, 2 or 3).
- `max_ratio = 0.5`         : the shock reaches at most this fraction of L_x at
                              `runtime` (sizes the x-extent).
- `rL = 8.0`                : transverse box size in thermal ion Larmor radii
                              ρ_i (used for dim ≥ 2; both transverse dims).
- `cells_per_skindepth = 4` : cells per `skindepth0`, i.e. dx = skindepth0/this.
                              `skindepth0` is Entity's *fiducial* skin depth (the
                              one it prints as `[N dx0]`), not the electron skin
                              depth: with the ppc0/2-per-species split the true
                              d_e = √2·skindepth0, so 4 here is 5.66 cells per d_e.
                              The summary reports both.
- `pow2 = true`             : round each resolution up to the next power of two.

# Problem generator / conventions
- `pgen = :downstream`        : `:downstream` (drift_ux), `:upstream` (piston_velocity)
      or `:upstream_fill` (piston_velocity + fill_fraction/wall_velocity).
- `beta_definition = :total`  : `:total` (β_i+β_e) or `:ion` (β_i) for `beta`.
- `mass_density = :ion`       : mass density behind v_A and c_s. `:ion` uses ρ = n_i m_i
      (the shock-literature convention, where m_e is negligible at the real mass ratio);
      `:total` uses ρ = n_i m_i + n_e m_e, correct for the plasma actually simulated.
      The two differ by sqrt(1+1/μ) (0.5% at μ = 100) and shift v_A and c_s together,
      so β and the ratio c_s/v_A are identical either way.
- `mach_convention = :flow`   : how `mach_sonic` maps to the velocity jump Δv:
      `:flow`  → Δv = M_s·c_s              (M_s = Δv/c_s; reproduces the beta50 ref);
      `:shock` → Δv = M_s·c_s·(r-1)/r      (M_s = standard upstream-frame shock Mach).
  The summary always reports both the flow and shock-frame Mach numbers.
- `compression_ratio = 0.0`   : shock compression ratio r (sizes the box and camera
      via v_shock = Δv/(r-1)). `0` (default) computes r analytically from the
      hydrodynamic Rankine-Hugoniot conditions: from the standard jump for the
      shock-frame Mach number (`mach_convention = :shock`), or from the self-consistent
      piston/reflecting-wall solution for the flow Mach number (`:flow`). Pass any
      value > 1 to override with a fixed r.
- `gamma_ad = 5/3`            : adiabatic index for the sound speed.

See the source header for the full unit chain. Remaining keywords (`skindepth0`,
`Bmag`, `Btheta`, `Bphi`, `ppc0`, `CFL`, `current_filters`, `n_outputs`,
`decomposition`, `name`, downstream-only `injection_frequency` /
`thermal_bath_*`, `match_ds`, `maxnpart_safety`, `verbose`) carry sensible
defaults.

# Moving wall (`pgen = :upstream_fill` only)
- `fill_fraction = 1.0`  : fraction of the box filled with plasma at t = 0. `1.0`
      disables the wall and the setup is identical to `:upstream`.
- `wall_velocity = 0.0`  : lab-frame speed of the wall closing off the filled region
      [c]. `0.0` (default) picks 2Δv, comfortably above the upstream-frame shock speed
      v_sh = Δv·r/(r-1) that it has to outrun. The summary reports the wall-to-shock
      separation and the time at which the wall retires at x1max.

# In-situ render / panning camera
- `render = true`             : emit an `[output.render]` block (N, Bmag, T22 scenes).
- `n_render_outputs = 1000`   : number of rendered frames over the run.
- `render_width = 2048`       : frame width in px (height follows the window aspect).
- `render_window_di = 100.0`  : x-window width in ion inertial lengths d_i.
- `render_shock_frac = 1/3`   : fractional position the shock front is held at.
- `render_view_azimuth = -55`, `render_view_elevation = 11.6` : (3D only) eye
      direction [deg] for the fixed 3/4 orthographic view; the box orientation in the
      frame is then independent of box size. `render_fill = 0.95` sets how much of the
      frame the region fills. In 3D the frame height and `[output.render.camera]`
      (position / look_at / up / ortho_height) are computed from the box geometry;
      preview with `entity/render_preview.py <toml>`.
  `camera_velocity` is prefilled with v_sh = Δv/(compression_ratio-1) (= [c] since
  Entity c=1) and `camera_start_time` with frac·window/v_sh. Both assume the hydro
  compression ratio, so they are ESTIMATES — measure the true front speed with
  `shock_3D/scripts/shock_speed.jl` and update `camera_velocity` for a long run.
"""
function write_shock_setup(path::AbstractString;
        # --- required physics ---
        mass_ratio::Real,
        temperature::Real,
        mach_sonic::Real,
        beta::Real,
        runtime::Real,
        # --- geometry / discretisation ---
        dim::Integer              = 2,
        max_ratio::Real           = 0.5,
        rL::Real                  = 8.0,
        cells_per_skindepth::Real = 4,
        pow2::Bool                = true,
        # --- pgen / conventions ---
        pgen::Symbol              = :downstream,
        beta_definition::Symbol   = :total,
        mach_convention::Symbol   = :flow,
        mass_density::Symbol      = :ion,
        compression_ratio::Real   = 0.0,
        gamma_ad::Real            = 5 / 3,
        # --- fixed scales / knobs ---
        skindepth0::Real          = 1.0,
        temperature_ratio::Real   = 1.0,
        Bmag::Real                = 1.0,
        Btheta::Real              = 0.0,
        Bphi::Real                = 0.0,
        ppc0::Real                = 8.0,
        CFL0::Real                 = 0.5,
        current_filters::Integer  = 4,
        name::AbstractString      = "_shock",
        n_outputs::Integer        = 100,
        decomposition             = nothing,
        match_ds::Real            = -1.0,
        maxnpart::Real            = 3.0e8,
        # --- downstream-only knobs ---
        injection_frequency::Integer = 100,
        thermal_bath_offset::Real    = 20.0,
        thermal_bath_width::Real     = 10.0,
        # --- upstream_fill-only knobs (moving wall) ---
        fill_fraction::Real          = 1.0,
        wall_velocity::Real          = 0.0,
        # --- in-situ render / panning camera ---
        render::Bool              = true,
        n_render_outputs::Integer = 1000,
        render_width::Integer     = 2048,
        render_window_di::Real    = 100.0,
        render_shock_frac::Real   = 1 / 3,
        # --- 3D orthographic camera (dim == 3 only) ---
        render_view_azimuth::Real   = -55.0,   # eye azimuth in the xy-plane [deg] (0 → +x, 90 → +y)
        render_view_elevation::Real = 11.6,    # eye elevation above the xy-plane [deg]
        render_fill::Real           = 0.95,    # region fills this fraction of the frame
        verbose::Bool             = true,
    )

    # ----------------------------- validation -----------------------------
    dim in (1, 2, 3)              || error("dim must be 1, 2 or 3 (got $dim)")
    pgen in (:downstream, :upstream, :upstream_fill) ||
        error("pgen must be :downstream, :upstream or :upstream_fill (got $pgen)")
    beta_definition in (:total, :ion) ||
        error("beta_definition must be :total or :ion (got $beta_definition)")
    mach_convention in (:flow, :shock) ||
        error("mach_convention must be :flow or :shock (got $mach_convention)")
    mass_density in (:ion, :total) ||
        error("mass_density must be :ion or :total (got $mass_density)")
    (temperature > 0 && beta > 0 && mass_ratio > 0 && mach_sonic > 0 && runtime > 0) ||
        error("mass_ratio, temperature, mach_sonic, beta and runtime must be positive")
    temperature_ratio > 0 || error("temperature_ratio must be positive")
    compression_ratio == 0 || compression_ratio > 1 ||
        error("compression_ratio must be 0 (auto: analytic Rankine-Hugoniot) or > 1 (got $compression_ratio)")
    0 < fill_fraction ≤ 1 || error("fill_fraction must be in (0, 1] (got $fill_fraction)")
    fill_fraction == 1 || pgen === :upstream_fill ||
        error("fill_fraction < 1 needs pgen = :upstream_fill (got $pgen)")

    μ  = float(mass_ratio)
    T  = float(temperature)              # T_e
    τ  = float(temperature_ratio)        # T_i/T_e
    d0 = float(skindepth0)

    # ------------------------- density split ------------------------------
    # Entity's [scales] are defined at the fiducial density n0, but the pgens
    # inject tot_number_density = 1 as ppc0/2 particles of each species per cell,
    # so each species sits at n_s = FN·n0. Everything below that carries a
    # density picks up the corresponding power of FN.
    fn = FN_PER_SPECIES
    # mass density in units of n_i·m_e: ions only, or the full plasma
    μρ = mass_density === :ion ? μ : μ + 1

    # -------------------- field scale from β and T (set larmor0) ----------
    #   β_s = 2·fn·T_s/σ0 ,  σ0 = (d0/larmor0)^2
    #     :total → β = β_e + β_i = 2·fn·T_e·(1+τ)/σ0
    #     :ion   → β = β_i       = 2·fn·T_e·τ/σ0
    fβ      = 2 * fn * (beta_definition === :total ? (1 + τ) : τ)
    σ0      = fβ * T / beta
    larmor0 = d0 / sqrt(σ0)

    # ----------------------------- speeds (c = 1) -------------------------
    vA = sqrt(σ0 / (fn * μρ))              # Alfvén speed       [c]
    cs = sqrt(gamma_ad * T * (1 + τ) / μρ) # sound speed        [c]

    # --- compression ratio r ---
    # If `compression_ratio == 0` (the default) it is computed analytically from the
    # hydrodynamic Rankine-Hugoniot conditions; otherwise the supplied value is used.
    #   :shock convention → `mach_sonic` IS the upstream (shock-frame) sonic Mach number,
    #                       so r follows from the standard jump (`_rh_shock`).
    #   :flow  convention → `mach_sonic` is the flow Mach M0 = Δv/c_s of the piston /
    #                       reflecting wall; downstream is at rest, so r is the
    #                       self-consistent piston solution (`_rh_piston`).
    # These are UNMAGNETISED jumps — an oblique/kinetic PIC shock (finite β, Btheta≠0)
    # will differ, so the render camera it feeds is still refined with shock_speed.jl.
    r_auto = compression_ratio == 0
    r_comp = r_auto ? (mach_convention === :shock ? _rh_shock(mach_sonic, gamma_ad) :
                                                    _rh_piston(mach_sonic, gamma_ad)) :
                      float(compression_ratio)
    r_comp > 1 || error("computed compression ratio r = $r_comp ≤ 1; check mach_sonic/gamma_ad")

    # velocity jump Δv = drift/piston speed.  `mach_sonic` is interpreted as:
    #   :flow  → M_s = Δv/c_s                 (Δv = M_s·c_s; matches beta50 ref)
    #   :shock → M_s = upstream-frame shock   (Δv = M_s·c_s·(r-1)/r, Rankine-Hugoniot)
    Δv = mach_convention === :flow ? mach_sonic * cs :
         mach_sonic * cs * (r_comp - 1) / r_comp

    # ------------------------ times / lengths (code units) ----------------
    Ωci_inv = larmor0 * μ                  # ion cyclotron time  [code time]
    t_run   = runtime * Ωci_inv            # total run time      [code time]
    dx      = d0 / cells_per_skindepth     # cell size           [code length]
    CFL     = CFL0 * sqrt(dim)             # entity CFL convention
    dt      = CFL * dx                     # timestep            [code time]

    d_e  = d0 / sqrt(fn)                   # electron skin depth [code length]
    d_i  = sqrt(μ / fn) * d0               # ion inertial length [code length]
    ρ_i  = larmor0 * sqrt(τ * T * μ)       # thermal ion Larmor radius [code length]
    λ_De = sqrt(T / fn) * d0               # electron Debye length [code length]

    # ------------------------------ box sizing ----------------------------
    # The shock recedes from the box reference frame at v_sh = Δv/(r-1); we want
    # it to reach max_ratio·L_x at t = t_run.
    v_sh_box = Δv / (r_comp - 1)
    Lx_min   = v_sh_box * t_run / max_ratio

    _res(L) = (n = ceil(Int, L / dx); pow2 ? nextpow(2, n) : (iseven(n) ? n : n + 1))

    Nx = _res(Lx_min)
    Lx = Nx * dx
    resolution = [Nx]
    extent     = [[0.0, Lx]]

    Ltr = rL * ρ_i                         # transverse physical size
    if dim ≥ 2
        Ny = _res(Ltr); Ly = Ny * dx
        push!(resolution, Ny); push!(extent, [-Ly / 2, Ly / 2])
    end
    if dim == 3
        Nz = _res(Ltr); Lz = Nz * dx
        push!(resolution, Nz); push!(extent, [-Lz / 2, Lz / 2])
    end

    # ------------------------------ moving wall ---------------------------
    # (:upstream_fill only) The wall starts at fill_fraction·L_x and moves at
    # v_wall in the lab frame while the window carries the box right at Δv, so it
    # gains on the box at (v_wall - Δv) and on the front at (v_wall - v_sh_up).
    # It retires once it reaches x1max, after which the box is filled everywhere.
    v_sh_up  = Δv * r_comp / (r_comp - 1)   # shock speed in the upstream frame
    v_wall   = wall_velocity > 0 ? float(wall_velocity) : 2 * Δv
    has_wall = pgen === :upstream_fill && fill_fraction < 1
    if has_wall
        v_wall > v_sh_up ||
            error("wall_velocity = $v_wall must exceed the shock speed " *
                  "v_sh = Δv·r/(r-1) = $v_sh_up, or the shock overruns the wall")
        v_wall < 1 || error("wall_velocity = $v_wall must be < c")
    end
    t_retire = has_wall ? (1 - fill_fraction) * Lx / (v_wall - Δv) : 0.0
    wall_gap = has_wall ? fill_fraction * Lx + (v_wall - v_sh_up) * t_retire : 0.0

    # ---------------------------- particle budget -------------------------
    ncells   = prod(resolution)

    # -------------------------- boundary conditions -----------------------
    if match_ds == -1.0
        match_ds = pgen === :downstream ? thermal_bath_offset + thermal_bath_width : 1.0
    end

    # ------------------------ render window / panning camera --------------
    # The shock recedes at v_sh = Δv/(r-1), the same speed used to size the box.
    # Entity code units have c = 1, so v_sh in [c] IS the coordinate pan rate the
    # renderer adds to x1_lim each step. The window is a fixed physical size; the
    # camera starts when the front reaches `render_shock_frac` of it (t = frac·W/v_sh)
    # and then pans at v_sh to hold it there. v_sh is only an ESTIMATE — it assumes
    # `compression_ratio` and hydro jump conditions, so an oblique/kinetic PIC shock
    # will differ (measure with scripts/shock_speed.jl and refine camera_velocity).
    cam_vel         = v_sh_box
    Wx_render       = min(round(render_window_di * d_i / dx) * dx, extent[1][2])
    cam_t0          = render_shock_frac * Wx_render / cam_vel
    render_interval = t_run / n_render_outputs

    # Frame height + (3D only) orthographic 3/4 camera. In 3D the height is set to the
    # projected aspect of the region so it wraps the box; the camera direction/up are
    # fixed (independent of box size) so the box orientation in the frame is invariant.
    cam3d = nothing
    if dim == 3
        region3 = ((0.0, Wx_render), (extent[2][1], extent[2][2]), (extent[3][1], extent[3][2]))
        u_eye, right, upc, center, diag, Δr, Δu =
            _cam3d_basis(region3, render_view_azimuth, render_view_elevation, (0.0, 0.0, 1.0))
        render_height = clamp(round(Int, render_width * Δu / Δr), 100, 4 * render_width)
        aspect        = render_width / render_height
        ortho_height  = 2 * max(0.5 * Δu / render_fill, 0.5 * Δr / (render_fill * aspect))
        D             = 2 * diag                       # eye a safe distance out (ortho: distance is cosmetic)
        cam3d = (; position = ntuple(d -> center[d] + D * u_eye[d], 3),
                   look_at = center, up = (0.0, 0.0, 1.0), ortho_height)
    else
        Ly_span       = dim ≥ 2 ? (extent[2][2] - extent[2][1]) : Wx_render
        render_height = clamp(round(Int, render_width * Ly_span / Wx_render), 100, render_width)
    end

    # ----------------------------- diagnostics ----------------------------
    derived = (; pgen, dim,
        larmor0, skindepth0 = d0, sigma0 = σ0, B0 = 1 / larmor0,
        mass_ratio = μ, temperature = T, temperature_ion = τ * T,
        temperature_ratio = τ, n_frac = fn, mass_density,
        compression_ratio = r_comp, compression_auto = r_auto,
        vA, cs, drift = Δv,
        c_over_vA   = 1 / vA,
        M_A_flow    = Δv / vA,             # Alfvén Mach (flow, = Δv/v_A)
        M_s_flow    = Δv / cs,             # sonic  Mach (flow, = Δv/c_s)
        M_A_shock   = Δv / vA * r_comp / (r_comp - 1),
        M_s_shock   = Δv / cs * r_comp / (r_comp - 1),
        beta_e = 2 * fn * T / σ0, beta_i = 2 * fn * τ * T / σ0,
        beta_total = 2 * fn * T * (1 + τ) / σ0,
        cs_over_vA = cs / vA,
        Omega_ci_inv = Ωci_inv, t_run, dx, dt,
        nsteps = round(Int, t_run / dt),
        d_e, d_i, rho_i = ρ_i,
        lambda_De = λ_De,
        fill_fraction = float(fill_fraction), wall_velocity = has_wall ? v_wall : 0.0,
        has_wall, t_retire, wall_gap, v_sh_upstream = v_sh_up,
        resolution, extent, ncells, maxnpart,
        v_sh = v_sh_box, cam_vel, cam_t0, Wx_render,
        render_interval, render_width, render_height, n_render_outputs, cam3d,
    )

    # ------------------------------ write toml -----------------------------
    open(path, "w") do io
        write(io, _shock_toml(io_ctx(; name, pgen, dim, beta_definition, beta,
            resolution, extent, decomposition, match_ds,
            CFL, current_filters, ppc0, maxnpart, μ,
            T, temperature_ratio, Bmag, Btheta, Bphi, Δv,
            t_run, n_outputs, injection_frequency,
            thermal_bath_offset, thermal_bath_width, derived,
            render, render_shock_frac, render_fill,
            render_view_azimuth, render_view_elevation)))
    end

    verbose && _print_shock_summary(path, derived; runtime, max_ratio, rL,
                                    cells_per_skindepth, pow2,
                                    gamma_ad, beta, beta_definition)
    return derived
end

# ----------------------------------------------------------------------------
#  Helpers
# ----------------------------------------------------------------------------

# Hydrodynamic Rankine-Hugoniot compression ratio for an upstream (shock-frame)
# sonic Mach number M and adiabatic index γ:  r = (γ+1)M² / ((γ-1)M² + 2).
_rh_shock(M, γ) = (γ + 1) * M^2 / ((γ - 1) * M^2 + 2)

# Self-consistent compression ratio for a piston / reflecting-wall shock driven at
# *flow* Mach number M0 = Δv/c_s. The downstream is at rest against the wall, so mass
# conservation gives v_sh = Δv/(r-1) and hence a shock-frame Mach M1 = M0·r/(r-1).
# Substituting M1 into `_rh_shock` and clearing denominators gives a quadratic in r,
#     [(γ-1)M0² + 2] r² − [(γ+1)M0² + 4] r + 2 = 0,
# whose physical (larger) root is taken. r → 1 as M0 → 0 and r → (γ+1)/(γ-1) as M0 → ∞.
function _rh_piston(M0, γ)
    a =  (γ - 1) * M0^2 + 2
    b = -((γ + 1) * M0^2 + 4)
    (-b + sqrt(max(b^2 - 8a, 0.0))) / (2a)      # 4·a·c with c = 2  ⇒  8a
end

# --- 3D orthographic camera geometry (mirrors entity renderer.cpp / composite.h) ---
_cross(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])
_unit(a) = (n = sqrt(a[1]^2 + a[2]^2 + a[3]^2); n > 1e-30 ? (a[1] / n, a[2] / n, a[3] / n) : a)

# projected spread (max − min) of a box's 8 corners along a screen axis
function _corner_spread(region, axis)
    lo, hi = Inf, -Inf
    for m in 0:7
        p = (iszero(m & 1) ? region[1][1] : region[1][2]) * axis[1] +
            (iszero(m & 2) ? region[2][1] : region[2][2]) * axis[2] +
            (iszero(m & 4) ? region[3][1] : region[3][2]) * axis[3]
        lo = min(lo, p); hi = max(hi, p)
    end
    hi - lo
end

# Build the orthographic camera basis for a 3D box `region` ((lo,hi) per axis) viewed
# from eye azimuth/elevation (deg, in the xy-plane / above it) with world-up `up`.
# Returns the eye-offset unit vector (look_at → eye), the screen right/up axes, the
# region center, its diagonal, and the projected corner spreads Δright/Δup. Because
# the projection is orthographic these depend only on the view *direction*, not the
# eye distance — so a fixed (az, el) reproduces the same box orientation at any size.
function _cam3d_basis(region, az_deg, el_deg, up)
    az, el = deg2rad(az_deg), deg2rad(el_deg)
    u_eye  = (cos(el) * cos(az), cos(el) * sin(az), sin(el))   # look_at → eye (unit)
    fwd    = (-u_eye[1], -u_eye[2], -u_eye[3])
    right  = _unit(_cross(fwd, up))
    upc    = _cross(right, fwd)                                # camera up (unit)
    center = ntuple(d -> 0.5 * (region[d][1] + region[d][2]), 3)
    diag   = sqrt(sum(ntuple(d -> (region[d][2] - region[d][1])^2, 3)))
    return u_eye, right, upc, center, diag,
           _corner_spread(region, right), _corner_spread(region, upc)
end

# format a real as a TOML float (always carries a '.' or exponent)
function _ftoml(x::Real)
    s = @sprintf("%.10g", float(x))
    (occursin('.', s) || occursin('e', s) || occursin('E', s)) ? s : s * ".0"
end

_fvec(v) = "[" * join(_ftoml.(v), ", ") * "]"
_iextent(ext) = "[" * join(["[" * _ftoml(e[1]) * ", " * _ftoml(e[2]) * "]" for e in ext], ", ") * "]"
_ires(res)    = "[" * join(string.(res), ", ") * "]"

function _bc_strings(dim::Integer)
    transverse = repeat(", [\"PERIODIC\"]", dim - 1)
    f = "[[\"CONDUCTOR\", \"MATCH\"]$transverse]"
    p = "[[\"REFLECT\", \"REFLECT\"]$transverse]"
    return f, p
end

# bundle the toml fields (keeps _shock_toml readable)
io_ctx(; kwargs...) = (; kwargs...)

function _shock_toml(c)
    fbc, pbc = _bc_strings(c.dim)
    decomp   = c.decomposition === nothing ? vcat(-1, fill(1, c.dim - 1)) : collect(c.decomposition)
    decomp_s = "[" * join(string.(decomp), ", ") * "]"
    interval = c.t_run / c.n_outputs
    bcomment = c.beta_definition === :total ? "total β = β_i + β_e" : "ion β = β_i"

    # --- in-situ render block with a shock-tracking camera (optional) ---
    d = c.derived
    render = if !c.render
        ""
    else
        pct  = round(Int, 100 * c.render_shock_frac)
        fpct = round(Int, 100 * c.render_fill)
        Nmax = round(1.5 * d.compression_ratio, digits = 1)   # density jumps to ~r behind the front
        # 3D renders composite the field volumetrically, so each scene needs an opacity
        # ramp: transparent at the scene's colorbar min, fully opaque at its max. 2D slice
        # renders ignore alpha, so only emit it for the 3D orthographic camera.
        _alpha3d(mn, mx) = d.cam3d === nothing ? "" :
            "\n      alpha    = [[$(_ftoml(mn)), 0.0], [$(_ftoml(mx)), 1.0]]"
        hcom = d.cam3d === nothing ? "~matches window aspect (x1 window : transverse extent)" :
                                     "matches the projected 3D box aspect"
        # 3D gets a fixed 3/4 orthographic camera; 2D uses the default slice window.
        camera = d.cam3d === nothing ? "" : """

    [output.render.camera]
      orthographic = true
      # Fixed 3/4 view (eye az=$(_ftoml(c.render_view_azimuth))°, el=$(_ftoml(c.render_view_elevation))°, up=+z): +x reads across the frame and the
      # box orientation is identical for any box size. look_at = region centre; ortho_height
      # fits the x1_lim region to ~$(fpct)% of the frame (aspect-aware). Preview: entity/render_preview.py.
      position     = $(_fvec(d.cam3d.position))
      look_at      = $(_fvec(d.cam3d.look_at))
      up           = $(_fvec(d.cam3d.up))
      ortho_height = $(_ftoml(d.cam3d.ortho_height))
      fov          = 15.0         # vertical FOV [deg] (perspective only)
"""
        """

  [output.render]
    enable            = true
    interval_time     = $(_ftoml(d.render_interval))   # $(d.n_render_outputs) frames over the run
    width             = $(d.render_width)
    height            = $(d.render_height)              # $hcom
    x1_lim            = [0.0, $(_ftoml(d.Wx_render))]  # window at the wall (wall + forming shock)
    axes              = true
    time_label        = true
    # Camera pans with the front so it stays at ~$(pct)% of the frame.
    # camera_velocity = v_sh = Δv/(r-1) = $(@sprintf("%.4g", d.cam_vel)) [c]  (Entity c=1 ⇒ this is also the x1_lim pan rate/time).
    # This is an ESTIMATE from r=$(@sprintf("%.3g", d.compression_ratio))$(d.compression_auto ? " (analytic RH)" : ""); an oblique/kinetic PIC shock will differ.
    # Refine by running PIC/Paper/shock_3D/scripts/shock_speed.jl on the field output.
    camera_velocity   = [$(_ftoml(d.cam_vel)), 0.0]
    camera_start_time = $(_ftoml(d.cam_t0))            # front reaches ~$(pct)% of the window here (t = frac·W/v_sh)
$camera
    [[output.render.scenes]]
      field    = "N"
      min      = 0.0
      max      = $(_ftoml(Nmax))
      colormap = "viridis"$(_alpha3d(0.0, Nmax))

    [[output.render.scenes]]
      field    = "Bmag"
      prefix   = "Bmag_"
      min      = 0.0
      max      = $(_ftoml(Nmax))                       # oblique field compresses ~like density; tune per run
      log      = false
      colormap = "dusk"$(_alpha3d(0.0, Nmax))

    [[output.render.scenes]]
      field    = "T22"
      prefix   = "T22_"
      min      = 0.01
      max      = 5.0                                    # tune per run
      log      = true
      colormap = "sunburst"$(_alpha3d(0.01, 5.0))
"""
    end

    # --- the [setup] block differs between the pgens ---
    if c.pgen === :upstream || c.pgen === :upstream_fill
        # the moving wall is only emitted when it is actually in use
        wall = if !d.has_wall
            c.pgen === :upstream_fill ?
                "\n  fill_fraction     = 1.0   # plasma fills the whole box at t = 0, no moving wall" : ""
        else
            "\n  fill_fraction     = $(_ftoml(d.fill_fraction))   # filled with plasma at t = 0; a wall closes off the rest\n" *
            "  wall_velocity     = $(_ftoml(d.wall_velocity))   # wall speed [c]; exceeds v_sh = $(@sprintf("%.4g", d.v_sh_upstream)) so the shock never catches it\n" *
            "  # the wall reaches x1max and retires at t = $(@sprintf("%.4g", d.t_retire)) = $(@sprintf("%.3g", d.t_retire / d.Omega_ci_inv)) Omega_ci^-1,\n" *
            "  # by which point it leads the front by $(@sprintf("%.4g", d.wall_gap)) = $(@sprintf("%.3g", d.wall_gap / d.d_i)) d_i"
        end
        setup = """
        [setup]
          temperature       = $(_ftoml(c.T))   # kB T_e / (m_e c^2)
          temperature_ratio = $(_ftoml(c.temperature_ratio))   # T_i/T_e
          Bmag              = $(_ftoml(c.Bmag))
          Btheta            = $(_ftoml(c.Btheta))   # field angle in plane [deg]
          Bphi              = $(_ftoml(c.Bphi))   # field angle out of plane [deg]
          piston_velocity   = $(_ftoml(c.Δv))   # = mach_sonic * c_s  [c]$wall
        """
    else
        setup = """
        [setup]
          drift_ux            = $(_ftoml(c.Δv))   # inflow speed toward wall [c] = mach_sonic * c_s
          temperature         = $(_ftoml(c.T))   # kB T / (m_e c^2); T_e = T_i
          temperature_ratio   = $(_ftoml(c.temperature_ratio))
          Bmag                = $(_ftoml(c.Bmag))
          Btheta              = $(_ftoml(c.Btheta))   # field angle in plane [deg]
          Bphi                = $(_ftoml(c.Bphi))   # field angle out of plane [deg]
          filling_fraction    = 1.0
          injector_velocity   = 0.0
          injection_start     = 0.0
          injection_frequency = $(c.injection_frequency)
          thermal_bath_offset = $(_ftoml(c.thermal_bath_offset))
          thermal_bath_width  = $(_ftoml(c.thermal_bath_width))
        """
    end

    return """
# Entity SRPIC collisionless-shock setup  (pgen = $(c.pgen))
# Generated by EntityTools.write_shock_setup
#   c/v_A = $(@sprintf("%.4g", d.c_over_vA))   M_A(flow) = $(@sprintf("%.4g", d.M_A_flow))   M_s(flow) = $(@sprintf("%.4g", d.M_s_flow))   $bcomment = $(@sprintf("%.4g", c.beta_definition === :total ? d.beta_total : d.beta_i))
#   each species is at n0/2, so d_e = $(@sprintf("%.4g", d.d_e / d.dx)) dx (skindepth0 = $(@sprintf("%.4g", d.skindepth0 / d.dx)) dx)
#   d_i = $(@sprintf("%.4g", d.d_i / d.dx)) dx   rho_i = $(@sprintf("%.4g", d.rho_i / d.dx)) dx   lambda_De = $(@sprintf("%.4g", d.lambda_De / d.dx)) dx
[simulation]
  name    = "$(c.name)"
  engine  = "srpic"
  runtime = $(_ftoml(c.t_run))

  [simulation.domain]
    decomposition = $decomp_s

  [simulation.domain.load_balance]
    enable     = true
    interval   = 10
    dimensions = [1]
    tolerance  = 0.1
    max_shift  = 2

[grid]
  resolution = $(_ires(c.resolution))
  extent     = $(_iextent(c.extent))

  [grid.metric]
    metric = "minkowski"

  [grid.boundaries]
    fields    = $fbc
    particles = $pbc

  [grid.boundaries.match]
    ds = $(_ftoml(c.match_ds))

[scales]
  larmor0    = $(_ftoml(c.derived.larmor0))
  skindepth0 = $(_ftoml(c.derived.skindepth0))

[algorithms]
  current_filters = $(c.current_filters)

  [algorithms.timestep]
    CFL = $(_ftoml(c.CFL))

[particles]
  ppc0 = $(_ftoml(c.ppc0))

  [[particles.species]]
    label    = "e-"
    mass     = 1.0
    charge   = -1.0
    maxnpart = $(@sprintf("%.4e", c.maxnpart))

  [[particles.species]]
    label    = "p+"
    mass     = $(_ftoml(c.μ))
    charge   = 1.0
    maxnpart = $(@sprintf("%.4e", c.maxnpart))

$setup
[output]
  interval_time = $(_ftoml(interval))
  format        = "BPFile"

  [output.fields]
    quantities = ["N", "B", "E", "T22"]

  [output.particles]
    enable = false
    stride = 512

  [output.spectra]
    enable = false
$render
[diagnostics]
  log_level       = "WARNING"
  blocking_timers = true

[checkpoint]
  interval = 100000
  keep     = 1
"""
end

function _print_shock_summary(path, d; runtime, max_ratio, rL, cells_per_skindepth,
                              pow2, gamma_ad, beta, beta_definition)
    println("┌─ write_shock_setup → $path")
    @printf("│  pgen = %-10s  dim = %d   r = %.3g%s   γ = %.4g\n", d.pgen, d.dim,
            d.compression_ratio, d.compression_auto ? " (RH)" : "", gamma_ad)
    println("│  ── requested ──")
    @printf("│   μ = %.4g   T_e = %.4g [m_e c²]   T_i/T_e = %.4g   M_s(flow) = %.4g   β(%s) = %.4g   t = %.4g Ω_ci⁻¹\n",
            d.mass_ratio, d.temperature, d.temperature_ratio, d.M_s_flow,
            beta_definition, beta, runtime)
    println("│  ── derived scales ──")
    @printf("│   larmor0 = %.6g   skindepth0 = %.6g   σ0 = %.4g   B0 = %.4g\n",
            d.larmor0, d.skindepth0, d.sigma0, d.B0)
    @printf("│   n_e = n_i = %.3g n0 (ppc0/2 per species)   ρ from %s masses\n",
            d.n_frac, d.mass_density)
    @printf("│   c/v_A = %.4g   v_A = %.4g c   c_s = %.4g c   c_s/v_A = %.4g   Δv = %.4g c\n",
            d.c_over_vA, d.vA, d.cs, d.cs_over_vA, d.drift)
    @printf("│   M_A(flow) = %.4g  M_s(flow) = %.4g   |  shock-frame: M_A = %.4g  M_s = %.4g\n",
            d.M_A_flow, d.M_s_flow, d.M_A_shock, d.M_s_shock)
    @printf("│   β_i = %.4g   β_e = %.4g   β_total = %.4g\n", d.beta_i, d.beta_e, d.beta_total)
    # delivered-vs-requested: β has to come back out of the written scales, and the
    # two independent routes to c_s/v_A have to agree
    βdel = beta_definition === :total ? d.beta_total : d.beta_i
    @printf("│   check: β(%s) delivered = %.6g vs requested %.6g   [%s]\n",
            beta_definition, βdel, beta, isapprox(βdel, beta; rtol = 1e-12) ? "ok" : "MISMATCH")
    @printf("│   check: c_s/v_A = %.6g vs sqrt(γ β_total/2) = %.6g   [%s]\n",
            d.cs_over_vA, sqrt(gamma_ad * d.beta_total / 2),
            isapprox(d.cs_over_vA, sqrt(gamma_ad * d.beta_total / 2); rtol = 1e-12) ?
                "ok" : "MISMATCH")
    println("│  ── lengths / time ──")
    @printf("│   d_e = %.4g (= %.4g dx)   d_i = %.4g (= %.4g dx)   ρ_i = %.4g (= %.4g dx)   λ_De = %.4g (= %.4g dx)\n",
            d.d_e, d.d_e / d.dx, d.d_i, d.d_i / d.dx,
            d.rho_i, d.rho_i / d.dx, d.lambda_De, d.lambda_De / d.dx)
    @printf("│   Ω_ci⁻¹ = %.6g   t_run = %.6g   dx = %.4g   dt = %.4g   steps = %d\n",
            d.Omega_ci_inv, d.t_run, d.dx, d.dt, d.nsteps)
    println("│  ── grid ──")
    @printf("│   resolution = %s   extent = %s\n", string(d.resolution), string(d.extent))
    @printf("│   cells = %.4g   maxnpart/species ≈ %.4e   (cells_per_skindepth=%g → %.3g per d_e, pow2=%s)\n",
            d.ncells, d.maxnpart, cells_per_skindepth, d.d_e / d.dx, pow2)
    if d.has_wall
        @printf("│   moving wall: fill_fraction = %.3g   v_wall = %.4g c (v_sh = %.4g c)   retires at t = %.4g (%.3g Ω_ci⁻¹), %.3g d_i ahead of the front\n",
                d.fill_fraction, d.wall_velocity, d.v_sh_upstream, d.t_retire,
                d.t_retire / d.Omega_ci_inv, d.wall_gap / d.d_i)
    end
    if d.lambda_De < d.dx
        @printf("│  ⚠ λ_De = %.3g < dx = %.3g (%.2f cells): risk of finite-grid heating; raise cells_per_skindepth or T.\n",
                d.lambda_De, d.dx, d.lambda_De / d.dx)
    end
    @printf("│   shock front reaches ≈ %.2g·L_x at t_run (target max_ratio = %.2g, rL = %g ρ_i)\n",
            d.drift / (d.compression_ratio - 1) * d.t_run / d.extent[1][2], max_ratio, rL)
    println("│  ── render / camera (ESTIMATE — refine v_sh with scripts/shock_speed.jl) ──")
    @printf("│   v_sh ≈ %.4g c = camera_velocity   window x1_lim = [0, %.4g]   camera_start_time = %.4g (%.3g Ω_ci⁻¹)\n",
            d.cam_vel, d.Wx_render, d.cam_t0, d.cam_t0 / d.Omega_ci_inv)
    @printf("│   %d frames  %dx%d px\n", d.n_render_outputs, d.render_width, d.render_height)
    if d.cam3d !== nothing
        @printf("│   3D ortho camera: look_at=%s  ortho_height=%.4g  (fixed 3/4 view; preview: entity/render_preview.py)\n",
                string(round.(d.cam3d.look_at; digits = 1)), d.cam3d.ortho_height)
    end
    println("└─")
    return nothing
end
