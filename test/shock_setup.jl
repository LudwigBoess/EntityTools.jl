# =============================================================================
#  Tests for EntityTools.write_shock_setup
# =============================================================================
#
#  The setup generator turns physical inputs (μ, T, M_s, β) into Entity's
#  [scales]. The risk is that the density normalisation drifts out of step with
#  what the pgen actually injects, which is silent: the run still starts, it is
#  just at the wrong β and M_A.
#
#  So the checks here never reuse the generator's algebra. They re-derive every
#  quantity from the *written toml* using Entity's primitive relations only:
#
#    ω_p0 = 1/skindepth0        B0 = 1/larmor0 = ω_B0
#    n_s  = ppc0/2 per cell     ω_ps² = ω_p0²·(n_s/n0)·q_s²/m_s
#    ω_Bs = q_s·B0/m_s          σ0 = (skindepth0/larmor0)²
#
#  and compare against what was asked for.
# =============================================================================

using Test
using TOML
using EntityTools

const TMPDIR = mktempdir()

# `_ftoml` writes 10 significant digits, so anything re-derived from a written
# file carries ~1e-10 of relative error. Identities evaluated entirely within
# one source are held to 1e-12.
const FTOL = 1e-9

# Re-derive the plasma parameters of a written setup from first principles.
# `nd` is the total number density injected by the pgen in units of n0 (the shock
# pgens all use 1), split evenly over the two species by ComputeNumInject.
function plasma_from_toml(path; nd = 1.0, gamma_ad = 5 / 3, mass_density = :ion)
    t   = TOML.parsefile(path)
    d0  = t["scales"]["skindepth0"]
    rL0 = t["scales"]["larmor0"]
    sp  = t["particles"]["species"]
    me, mi = sp[1]["mass"], sp[2]["mass"]
    qe, qi = sp[1]["charge"], sp[2]["charge"]
    Te  = t["setup"]["temperature"]                 # kB T_e / (m_e c²)
    τ   = t["setup"]["temperature_ratio"]
    Ti  = τ * Te

    ωp0 = 1 / d0                                    # fiducial plasma frequency
    B0  = 1 / rL0                                   # ω_B0 for m = q = 1
    σ0  = (d0 / rL0)^2
    fs  = nd / 2                                    # n_s/n0 for each species

    ωpe = ωp0 * sqrt(fs * qe^2 / me)                # per-species plasma frequency
    ωpi = ωp0 * sqrt(fs * qi^2 / mi)
    ωci = abs(qi) * B0 / mi                         # ion cyclotron frequency

    ρ   = fs * (mass_density === :ion ? mi : mi + me)   # mass density [n0·m_e]
    vA  = B0 / sqrt(ρ)                              # B²/ρ in Entity units (4π absorbed)
    cs  = sqrt(gamma_ad * fs * (Te + Ti) / ρ)       # γP/ρ

    return (; d0, larmor0 = rL0, sigma0 = σ0, B0, mass_ratio = mi / me, Te, Ti, τ,
              d_e = 1 / ωpe, d_i = 1 / ωpi, omega_ci = ωci, vA, cs,
              beta_e = 2 * fs * Te / σ0, beta_i = 2 * fs * Ti / σ0,
              beta_total = 2 * fs * (Te + Ti) / σ0,
              lambda_De = sqrt(Te / me) / ωpe,
              rho_i = mi * sqrt(Ti / mi) / (abs(qi) * B0))
end

@testset "write_shock_setup" begin

    # Park, Caprioli & Spitkovsky 2015: μ = 100, T = 1.12e-3 m_e c², v_u = 0.1 c,
    # and the paper's own Mach numbers M_s = v_sh/sqrt(T_p/m_p) = 40, M_A = 20,
    # which together fix β_total = 4(M_A/M_s)² = 1.
    μ, T, γ, r = 100.0, 1.12e-3, 5 / 3, 4.0
    β_req      = 1.0
    v_u        = 0.1
    cs_req     = sqrt(2 * γ * T / μ)          # :ion mass density, T_e = T_i
    path       = joinpath(TMPDIR, "park.toml")

    d = write_shock_setup(path;
        mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req,
        beta = β_req, beta_definition = :total, mach_convention = :flow,
        compression_ratio = r, gamma_ad = γ, runtime = 310.0,
        dim = 1, cells_per_skindepth = 4, pow2 = false, max_ratio = 0.386,
        pgen = :upstream_fill, Btheta = 30.0, Bphi = 90.0,
        ppc0 = 400.0, CFL0 = 0.45, render = false, verbose = false)

    p = plasma_from_toml(path; gamma_ad = γ)

    @testset "written scales reproduce the request" begin
        # β must survive the round-trip through larmor0 — this is the check that
        # the n0/2 density split is handled consistently on both sides.
        @test p.beta_total ≈ β_req         rtol = FTOL
        @test p.beta_i     ≈ β_req / 2     rtol = FTOL
        @test p.beta_e     ≈ β_req / 2     rtol = FTOL
        @test p.cs         ≈ cs_req        rtol = FTOL
        @test p.mass_ratio == μ
        @test p.Te         == T
        # the piston/drift speed is the requested velocity jump
        @test TOML.parsefile(path)["setup"]["piston_velocity"] ≈ v_u rtol = FTOL
        # and the paper's Mach numbers come back out of the written file
        v_sh = v_u * r / (r - 1)
        @test v_sh / sqrt(p.Ti / p.mass_ratio) ≈ 40.0 rtol = 0.01
        @test v_sh / p.vA                      ≈ 20.0 rtol = 0.01
    end

    @testset "generator agrees with the first-principles rederivation" begin
        @test d.sigma0     ≈ p.sigma0     rtol = FTOL
        @test d.vA         ≈ p.vA         rtol = FTOL
        @test d.cs         ≈ p.cs         rtol = FTOL
        @test d.beta_i     ≈ p.beta_i     rtol = FTOL
        @test d.beta_e     ≈ p.beta_e     rtol = FTOL
        @test d.beta_total ≈ p.beta_total rtol = FTOL
        @test d.d_e        ≈ p.d_e        rtol = FTOL
        @test d.d_i        ≈ p.d_i        rtol = FTOL
        @test d.rho_i      ≈ p.rho_i      rtol = FTOL
        @test d.lambda_De  ≈ p.lambda_De  rtol = FTOL
        @test d.Omega_ci_inv ≈ 1 / p.omega_ci rtol = FTOL
    end

    @testset "plasma identities" begin
        # v_A = Ω_ci·d_i holds exactly for the ion-mass Alfvén speed
        @test p.vA ≈ p.omega_ci * p.d_i rtol = 1e-12
        # skin-depth ordering, and the split that makes d_e ≠ skindepth0
        @test p.d_i / p.d_e ≈ sqrt(μ)         rtol = 1e-12
        @test p.d_e ≈ sqrt(2) * p.d0          rtol = 1e-12
        @test p.lambda_De / p.d_e ≈ sqrt(T)   rtol = 1e-12
        # gyroradius / inertial length ratio, and the two routes to c_s/v_A
        @test p.rho_i / p.d_i ≈ sqrt(p.beta_i / 2)          rtol = 1e-12
        @test p.cs / p.vA ≈ sqrt(γ * p.beta_total / 2)      rtol = 1e-12
        # density-free quantities must not move when the split changes
        @test p.rho_i ≈ p.larmor0 * sqrt(T * μ)             rtol = 1e-12
        @test 1 / p.omega_ci ≈ p.larmor0 * μ                rtol = 1e-12
    end

    @testset "grid and timestep" begin
        t  = TOML.parsefile(path)
        dx = (t["grid"]["extent"][1][2] - t["grid"]["extent"][1][1]) /
             t["grid"]["resolution"][1]
        @test dx ≈ p.d0 / 4 rtol = 1e-12                     # cells_per_skindepth = 4
        @test d.dt ≈ t["algorithms"]["timestep"]["CFL"] * dx rtol = 1e-12   # Entity: dt = CFL·dx0
        @test d.nsteps == round(Int, d.t_run / d.dt)
        # the front is placed at max_ratio·L_x at the end of the run
        @test d.drift / (r - 1) * d.t_run / t["grid"]["extent"][1][2] ≈ 0.386 rtol = 0.02
    end

    @testset "temperature_ratio enters β and c_s" begin
        # T_i = 2 T_e at fixed β_total must give β_i = 2 β_e and the same c_s/v_A
        p2 = joinpath(TMPDIR, "tratio.toml")
        d2 = write_shock_setup(p2;
            mass_ratio = μ, temperature = T, temperature_ratio = 2.0,
            mach_sonic = 10.0, beta = β_req, beta_definition = :total,
            compression_ratio = r, gamma_ad = γ, runtime = 10.0, dim = 1,
            cells_per_skindepth = 4, pow2 = false, pgen = :upstream,
            ppc0 = 400.0, render = false, verbose = false)
        q = plasma_from_toml(p2; gamma_ad = γ)
        @test q.beta_total ≈ β_req      rtol = FTOL
        @test q.beta_i ≈ 2 * q.beta_e   rtol = FTOL
        @test d2.beta_i ≈ q.beta_i      rtol = FTOL
        @test d2.cs ≈ q.cs              rtol = FTOL
        @test d2.cs / d2.vA ≈ sqrt(γ * q.beta_total / 2) rtol = FTOL
        @test d2.rho_i ≈ d2.larmor0 * sqrt(2 * T * μ)    rtol = FTOL
    end

    @testset ":ion vs :total mass density" begin
        # switching the mass density shifts v_A and c_s together, so β and the
        # ratio c_s/v_A are untouched and only the Mach numbers move by √(1+1/μ)
        p3 = joinpath(TMPDIR, "mtotal.toml")
        d3 = write_shock_setup(p3;
            mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req,
            beta = β_req, mass_density = :total, compression_ratio = r,
            gamma_ad = γ, runtime = 310.0, dim = 1, cells_per_skindepth = 4,
            pow2 = false, pgen = :upstream, ppc0 = 400.0, render = false, verbose = false)
        @test d3.sigma0 ≈ d.sigma0                    rtol = 1e-12
        @test d3.beta_total ≈ d.beta_total            rtol = 1e-12
        @test d3.cs_over_vA ≈ d.cs_over_vA            rtol = 1e-12
        @test d3.vA / d.vA ≈ 1 / sqrt(1 + 1 / μ)      rtol = 1e-12
        @test d3.cs / d.cs ≈ 1 / sqrt(1 + 1 / μ)      rtol = 1e-12
    end

    @testset "beta_definition = :ion" begin
        p4 = joinpath(TMPDIR, "betaion.toml")
        write_shock_setup(p4;
            mass_ratio = μ, temperature = T, mach_sonic = 10.0, beta = 0.5,
            beta_definition = :ion, compression_ratio = r, gamma_ad = γ,
            runtime = 10.0, dim = 1, cells_per_skindepth = 4, pow2 = false,
            pgen = :upstream, ppc0 = 400.0, render = false, verbose = false)
        q = plasma_from_toml(p4; gamma_ad = γ)
        @test q.beta_i ≈ 0.5     rtol = FTOL
        @test q.beta_total ≈ 1.0 rtol = FTOL
    end

    @testset "moving wall" begin
        p5 = joinpath(TMPDIR, "wall.toml")
        d5 = write_shock_setup(p5;
            mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req, beta = β_req,
            compression_ratio = r, gamma_ad = γ, runtime = 310.0, dim = 1,
            cells_per_skindepth = 4, pow2 = false, max_ratio = 0.386,
            pgen = :upstream_fill, fill_fraction = 0.3, ppc0 = 400.0,
            render = false, verbose = false)
        s = TOML.parsefile(p5)["setup"]
        @test s["fill_fraction"] ≈ 0.3
        # default wall speed is 2Δv, which has to beat the upstream-frame shock speed
        @test s["wall_velocity"] ≈ 2 * v_u rtol = 1e-12
        @test s["wall_velocity"] > d5.v_sh_upstream
        # the wall must still be ahead of the front when it retires
        @test d5.t_retire > 0
        @test d5.wall_gap > 0.3 * d5.extent[1][2]
        # a wall slower than the shock is rejected, and fill_fraction < 1 needs the pgen
        @test_throws ErrorException write_shock_setup(joinpath(TMPDIR, "bad1.toml");
            mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req, beta = β_req,
            compression_ratio = r, runtime = 10.0, dim = 1, pgen = :upstream_fill,
            fill_fraction = 0.3, wall_velocity = 0.5 * v_u, render = false, verbose = false)
        @test_throws ErrorException write_shock_setup(joinpath(TMPDIR, "bad2.toml");
            mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req, beta = β_req,
            compression_ratio = r, runtime = 10.0, dim = 1, pgen = :upstream,
            fill_fraction = 0.3, render = false, verbose = false)
        # fill_fraction = 1 emits no wall keys at all
        @test !haskey(TOML.parsefile(path)["setup"], "wall_velocity")
    end

    @testset "pgen setup blocks" begin
        @test haskey(TOML.parsefile(path)["setup"], "piston_velocity")
        p6 = joinpath(TMPDIR, "down.toml")
        write_shock_setup(p6;
            mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req, beta = β_req,
            compression_ratio = r, gamma_ad = γ, runtime = 10.0, dim = 1,
            cells_per_skindepth = 4, pow2 = false, pgen = :downstream,
            ppc0 = 400.0, render = false, verbose = false)
        s = TOML.parsefile(p6)["setup"]
        @test haskey(s, "drift_ux")
        @test s["drift_ux"] ≈ v_u rtol = 1e-12
        @test haskey(s, "thermal_bath_width")
    end
end
