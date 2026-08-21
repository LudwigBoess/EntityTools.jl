# =============================================================================
#  Tests for EntityTools.EntityUnits
# =============================================================================
#
#  The conversion factors are only useful if they are the *same* unit chain the
#  setup generator uses to pick Entity's [scales] — a run configured for β = 1
#  and M_A = 20 has to read back as β = 1 and M_A = 20. So the checks here
#
#    1. compare `EntityUnits` against `write_shock_setup`, which derives the
#       identical quantities by an independent route (physics -> scales, versus
#       scales -> physics here);
#    2. re-derive every scale from Entity's primitive relations only
#         ω_p0 = 1/skindepth0   B0 = 1/larmor0 = ω_B0   σ0 = (skindepth0/larmor0)²
#         ω_ps² = ω_p0²·(n_s/n0)·q_s²/m_s        ω_Bs = q_s·B0/m_s
#       and hold the dimensionless identities that follow from them;
#    3. check that the multiplicative conversion factors are the reciprocals of
#       the scales they name, which is the part that is easy to get backwards.
# =============================================================================

using Test
using TOML
using EntityTools

const UTMP = mktempdir()

# `_ftoml` writes 10 significant digits, so anything routed through a written
# parameter file carries ~1e-10 of relative error.
const UTOL = 1e-9

@testset "EntityUnits" begin

    # Park, Caprioli & Spitkovsky 2015: μ = 100, T = 1.12e-3 m_e c², v_u = 0.1 c,
    # M_s = 40, M_A = 20 and hence β_total = 4(M_A/M_s)² = 1.
    μ, T, γ, r = 100.0, 1.12e-3, 5 / 3, 4.0
    β_req, v_u = 1.0, 0.1
    cs_req     = sqrt(2 * γ * T / μ)
    fn         = 0.5                    # ppc0/2 per species, as the shock pgens inject
    path       = joinpath(UTMP, "units_park.toml")

    d = write_shock_setup(path;
        mass_ratio = μ, temperature = T, mach_sonic = v_u / cs_req,
        beta = β_req, beta_definition = :total, mach_convention = :flow,
        compression_ratio = r, gamma_ad = γ, runtime = 310.0,
        dim = 1, cells_per_skindepth = 4, pow2 = false, max_ratio = 0.386,
        pgen = :upstream_fill, Btheta = 30.0, Bphi = 90.0,
        ppc0 = 400.0, CFL0 = 0.45, render = false, verbose = false)

    u = EntityUnits(path; n_frac = fn, gamma_ad = γ)

    @testset "parameter file round-trip" begin
        # the species masses and the temperature come out of the written file
        @test u.mass_ratio        == μ
        @test u.temperature       ≈ T   rtol = UTOL
        @test u.temperature_ratio ≈ 1.0 rtol = UTOL
        @test u.dim               == 1
        @test u.n_frac            == fn
        # the grid geometry: cells_per_skindepth = 4 and V0 = dx^dim
        @test u.dx ≈ u.skindepth0 / 4 rtol = UTOL
        @test u.V0 ≈ u.dx^u.dim       rtol = UTOL
        @test u.n0 ≈ u.ppc0 / u.V0    rtol = UTOL
    end

    @testset "agrees with write_shock_setup" begin
        @test u.skindepth0 ≈ d.skindepth0   rtol = UTOL
        @test u.larmor0    ≈ d.larmor0      rtol = UTOL
        @test u.sigma0     ≈ d.sigma0       rtol = UTOL
        @test u.B0         ≈ d.B0           rtol = UTOL
        @test u.d_e        ≈ d.d_e          rtol = UTOL
        @test u.d_i        ≈ d.d_i          rtol = UTOL
        @test u.rho_i      ≈ d.rho_i        rtol = UTOL
        @test u.lambda_De  ≈ d.lambda_De    rtol = UTOL
        @test u.v_A        ≈ d.vA           rtol = UTOL
        @test u.c_s        ≈ d.cs           rtol = UTOL
        @test u.beta_e     ≈ d.beta_e       rtol = UTOL
        @test u.beta_i     ≈ d.beta_i       rtol = UTOL
        @test u.beta       ≈ d.beta_total   rtol = UTOL
        @test 1 / u.Omega_ci ≈ d.Omega_ci_inv rtol = UTOL
        # and the physics that was asked for is what comes back out
        @test u.beta ≈ β_req rtol = UTOL
        v_sh = v_u * r / (r - 1)
        @test v_sh * u.v_vA ≈ 20.0 rtol = 0.01
        @test v_sh / sqrt(u.temperature * u.temperature_ratio / μ) ≈ 40.0 rtol = 0.01
    end

    @testset "primitive relations" begin
        d0, r0 = u.skindepth0, u.larmor0
        @test u.sigma0   ≈ (d0 / r0)^2      rtol = 1e-12
        @test u.B0       == u.omegaB0
        @test u.q0       ≈ u.V0 / (u.ppc0 * d0^2) rtol = 1e-12
        # per-species plasma / cyclotron frequencies
        @test u.omega_pe ≈ sqrt(fn) / d0       rtol = 1e-12
        @test u.omega_pi ≈ sqrt(fn / μ) / d0   rtol = 1e-12
        @test u.Omega_ce ≈ 1 / r0              rtol = 1e-12
        @test u.Omega_ci ≈ 1 / (r0 * μ)        rtol = 1e-12
        # inertial lengths are the inverse frequencies (c = 1)
        @test u.d_e ≈ 1 / u.omega_pe rtol = 1e-12
        @test u.d_i ≈ 1 / u.omega_pi rtol = 1e-12
        # thermal speeds and the gyroradii built from them
        @test u.v_the ≈ sqrt(u.temperature)                       rtol = 1e-12
        @test u.v_thi ≈ sqrt(u.temperature * u.temperature_ratio / μ) rtol = 1e-12
        @test u.rho_e ≈ u.v_the / u.Omega_ce rtol = 1e-12
        @test u.rho_i ≈ u.v_thi / u.Omega_ci rtol = 1e-12
        @test u.lambda_De ≈ u.v_the / u.omega_pe rtol = 1e-12
    end

    @testset "dimensionless identities" begin
        @test u.v_A ≈ u.Omega_ci * u.d_i             rtol = 1e-12
        @test u.d_i / u.d_e ≈ sqrt(μ)                rtol = 1e-12
        @test u.d_e ≈ sqrt(2) * u.skindepth0         rtol = 1e-12   # n_frac = 1/2
        @test u.lambda_De / u.d_e ≈ sqrt(T)          rtol = 1e-12
        @test u.rho_e / u.d_e ≈ sqrt(u.beta_e / 2)   rtol = 1e-12
        @test u.rho_i / u.d_i ≈ sqrt(u.beta_i / 2)   rtol = 1e-12
        @test u.c_s / u.v_A ≈ sqrt(γ * u.beta / 2)   rtol = 1e-12
        @test u.Omega_ce / u.omega_pe ≈ sqrt(u.sigma0 / fn) rtol = 1e-12
        # density-free quantities: unchanged by the n0/2 split
        @test u.rho_i ≈ u.larmor0 * sqrt(T * μ) rtol = 1e-12
        @test 1 / u.Omega_ci ≈ u.larmor0 * μ    rtol = 1e-12
        @test u.c_s ≈ sqrt(2 * γ * T / μ)       rtol = 1e-12
    end

    @testset "conversion factors are the reciprocal scales" begin
        # a code position of exactly one d_i must convert to 1.0
        @test u.d_i * u.x_di               ≈ 1.0 rtol = 1e-12
        @test u.d_e * u.x_de               ≈ 1.0 rtol = 1e-12
        @test u.rho_e * u.x_rho_e          ≈ 1.0 rtol = 1e-12
        @test u.rho_i * u.x_rho_i          ≈ 1.0 rtol = 1e-12
        @test u.lambda_De * u.x_lambda_De  ≈ 1.0 rtol = 1e-12
        @test u.dx * u.x_cells             ≈ 1.0 rtol = 1e-12
        # times: multiplying by the factor gives t in inverse-frequency units
        @test u.t_Omega_ci ≈ u.Omega_ci rtol = 1e-12
        @test u.t_Omega_ce ≈ u.Omega_ce rtol = 1e-12
        @test u.t_omega_pe ≈ u.omega_pe rtol = 1e-12
        @test u.t_omega_pi ≈ u.omega_pi rtol = 1e-12
        @test (1 / u.Omega_ci) * u.t_Omega_ci ≈ 1.0 rtol = 1e-12
        # speeds
        @test u.v_A * u.v_vA ≈ 1.0 rtol = 1e-12
        @test u.c_s * u.v_cs ≈ 1.0 rtol = 1e-12
        # the whole point: the box in d_i and the runtime in Ω_ci⁻¹
        t = TOML.parsefile(path)
        @test t["grid"]["extent"][1][2] * u.x_di ≈ d.extent[1][2] / d.d_i rtol = UTOL
        @test t["simulation"]["runtime"] * u.t_Omega_ci ≈ 310.0           rtol = UTOL
    end

    @testset "n_frac scaling" begin
        # only the density-dependent scales move with n_frac, by its power
        u1 = EntityUnits(path; n_frac = 1.0, gamma_ad = γ)
        @test u1.d_e ≈ u.d_e / sqrt(2)   rtol = 1e-12
        @test u1.d_i ≈ u.d_i / sqrt(2)   rtol = 1e-12
        @test u1.lambda_De ≈ u.lambda_De / sqrt(2) rtol = 1e-12
        @test u1.v_A ≈ u.v_A / sqrt(2)   rtol = 1e-12
        @test u1.beta ≈ 2 * u.beta       rtol = 1e-12
        @test u1.omega_pe ≈ u.omega_pe * sqrt(2) rtol = 1e-12
        # density-free scales must not budge
        @test u1.rho_e ≈ u.rho_e         rtol = 1e-12
        @test u1.rho_i ≈ u.rho_i         rtol = 1e-12
        @test u1.c_s ≈ u.c_s             rtol = 1e-12
        @test u1.Omega_ci ≈ u.Omega_ci   rtol = 1e-12
        @test u1.sigma0 ≈ u.sigma0       rtol = 1e-12
        # d_e = skindepth0 exactly at the fiducial density
        @test u1.d_e ≈ u1.skindepth0     rtol = 1e-12
    end

    @testset "mass_density convention" begin
        # :total shifts v_A and c_s together, leaving β and c_s/v_A untouched
        ut = EntityUnits(path; n_frac = fn, gamma_ad = γ, mass_density = :total)
        @test ut.mu_rho == μ + 1
        @test ut.v_A / u.v_A ≈ 1 / sqrt(1 + 1 / μ) rtol = 1e-12
        @test ut.c_s / u.c_s ≈ 1 / sqrt(1 + 1 / μ) rtol = 1e-12
        @test ut.beta ≈ u.beta                     rtol = 1e-12
        @test ut.c_s / ut.v_A ≈ u.c_s / u.v_A      rtol = 1e-12
    end

    @testset "legacy helpers stay consistent" begin
        # the T-parameterised helpers must agree with the stored fields
        @test rho_e(u, u.temperature) ≈ u.rho_e   rtol = 1e-12
        @test v_th(u.temperature)     ≈ u.v_the   rtol = 1e-12
        @test betae(u, u.temperature) ≈ u.beta_e  rtol = 1e-12
        @test omega_pe(u)             == u.omega_pe
        @test omega_ce(u)             == u.Omega_ce
        # v_A/c = Ω_B/ω_p from first principles (v_A² = B²/4πnm, ω_p² = 4πnq²/m),
        # so the magnetisation σ0 is v_A²/c² for m = m0 = 1 at n = n0. `vA` is
        # that fiducial speed; `u.v_A` carries the ion mass and n_frac.
        @test EntityTools.vA(u) ≈ u.skindepth0 / u.larmor0 rtol = 1e-12
        @test EntityTools.vA(u)^2 ≈ u.sigma0               rtol = 1e-12
        @test EntityTools.vA(u) ≈ u.v_A * sqrt(u.n_frac * u.mu_rho) rtol = 1e-12
        @test EntityTools.vA(u) < 1                        # it is a speed, in units of c
        # and the two definitions coincide for a unit-mass species at n0
        u4 = EntityUnits(u.skindepth0, u.larmor0; mass_ratio = 1, n_frac = 1)
        @test EntityTools.vA(u4) ≈ u4.v_A rtol = 1e-12
        # the two-argument constructor still works and is fiducial
        u2 = EntityUnits(u.skindepth0, u.larmor0)
        @test u2.sigma0 ≈ u.sigma0 rtol = 1e-12
        @test u2.d_e    ≈ u2.skindepth0 rtol = 1e-12
        @test isnan(u2.temperature)
        @test isnan(u2.beta)
    end

    @testset "temperature-free runs report NaN" begin
        # a setup without a temperature must not silently claim β = 0
        u3 = EntityUnits(1.0, 10.0; mass_ratio = 100.0, n_frac = 0.5)
        @test isnan(u3.beta) && isnan(u3.rho_i) && isnan(u3.c_s)
        @test isnan(u3.x_rho_i) && isnan(u3.v_cs)
        # while the temperature-free scales are all finite
        @test all(isfinite, (u3.d_e, u3.d_i, u3.v_A, u3.Omega_ci, u3.x_di, u3.t_Omega_ci))
    end

    @testset "species table from a parameter file" begin
        sp = [(; label = "e-", mass = 1.0, charge = -1.0),
              (; label = "p+", mass = 100.0, charge = 1.0)]
        @test mass_ratio_from_species(sp) == 100.0
        @test mass_ratio_from_species(sp[1:1]) == 1.0
        @test mass_ratio_from_species(NamedTuple[]) === nothing
        # an explicit mass_ratio overrides whatever the file says
        @test EntityUnits(path; mass_ratio = 1.0).mass_ratio == 1.0
    end

    @testset "info-file species parser" begin
        info = joinpath(UTMP, "run.info")
        # the per-domain memory footprint uses the same `Species #N` prefix and
        # must not be mistaken for a species declaration
        write(info, """
            Particles
            - Species #1:
              - Label: e-
              - Mass: 1.0
              - Charge: -1.0
              - Pusher: boris
            - Species #2:
              - Label: p+
              - Mass: 100.0
              - Charge: 1.0
              - Pusher: boris
            Domains
            - Domain #0:
              - Memory footprint:
                - Species #1 (e-): 9.50 GB
                - Species #2 (p+): 9.50 GB
            """)
        sp = read_species(info)
        @test length(sp) == 2
        @test sp[1].label == "e-" && sp[1].mass == 1.0 && sp[1].charge == -1.0
        @test sp[2].label == "p+" && sp[2].mass == 100.0 && sp[2].charge == 1.0
        @test mass_ratio_from_species(sp) == 100.0
        # a missing file is not an error, just an empty table
        @test isempty(read_species(joinpath(UTMP, "does_not_exist.info")))
    end

    @testset "invalid input is rejected" begin
        @test_throws ErrorException EntityUnits(-1.0, 10.0)
        @test_throws ErrorException EntityUnits(1.0, -10.0)
        @test_throws ErrorException EntityUnits(1.0, 10.0; n_frac = 0.0)
        @test_throws ErrorException EntityUnits(1.0, 10.0; mass_ratio = -100.0)
        @test_throws ErrorException EntityUnits(1.0, 10.0; temperature = -1.0)
        @test_throws ErrorException EntityUnits(1.0, 10.0; mass_density = :electron)
        @test_throws ErrorException EntityUnits(1.0, 10.0; dim = 4)
        @test_throws ErrorException EntityUnits(joinpath(UTMP, "no_such.toml"))
    end
end
