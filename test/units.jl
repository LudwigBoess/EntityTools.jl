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

    # -------------------------------------------------------------------------
    #  EntityUnitsCGS: the same run in cgs-Gaussian units.
    #
    #  `EntityUnitsCGS` is built by rescaling `EntityUnits`, so the checks below
    #  re-derive every field from the *cgs* plasma relations
    #    omega_p = sqrt(4 pi n q^2/m)   Omega_B = q B/(m c)   v_A = B/sqrt(4 pi rho)
    #  instead, plus the NRL-formulary numbers, so an error in the unit chain
    #  cannot cancel out of both sides.
    # -------------------------------------------------------------------------
    @testset "EntityUnitsCGS" begin

        cc, qe, me = EntityTools.c_cgs, EntityTools.q_e_cgs, EntityTools.m_e_cgs
        kB, eV     = EntityTools.k_B_cgs, EntityTools.eV_cgs

        n_e = 5.0                       # electrons per cm^3, solar-wind-like
        uc  = EntityUnitsCGS(u; n_e)

        @testset "field parity with EntityUnits" begin
            fu, fc = fieldnames(EntityUnits), fieldnames(EntityUnitsCGS)
            @test fc[1:length(fu)] == fu
            @test fc[length(fu)+1:end] == (:length_unit, :time_unit, :speed_unit,
                                           :mass_unit, :charge_unit, :Bfield_unit,
                                           :energy_unit)
            @test uc isa EntityUnitsCGS{Float64}
            # both show methods run
            @test !isempty(sprint(show, uc))
            @test !isempty(sprint(show, MIME"text/plain"(), uc))
        end

        @testset "dimensionless fields are unchanged" begin
            @test uc.sigma0            == u.sigma0
            @test uc.beta_e            == u.beta_e
            @test uc.beta_i            == u.beta_i
            @test uc.beta              == u.beta
            @test uc.mass_ratio        == u.mass_ratio
            @test uc.n_frac            == u.n_frac
            @test uc.mu_rho            == u.mu_rho
            @test uc.temperature_ratio == u.temperature_ratio
            @test uc.gamma_ad          == u.gamma_ad
            @test uc.ppc0              == u.ppc0
            @test uc.dim               == u.dim
        end

        @testset "cgs plasma relations" begin
            n_i = uc.n_frac * uc.n0             # ion density [cm^-3]
            m_i = uc.mass_ratio * me            # ion mass [g]
            ρ   = uc.mu_rho * n_i * me          # mass density [g/cm^3]
            Te  = kB * uc.temperature           # [erg]
            Ti  = uc.temperature_ratio * Te

            # the anchor comes back out
            @test uc.n_frac * uc.n0 ≈ n_e rtol = 1e-12

            # frequencies
            @test uc.omega_pe ≈ sqrt(4π * n_e * qe^2 / me)  rtol = 1e-12
            @test uc.omega_pi ≈ sqrt(4π * n_i * qe^2 / m_i) rtol = 1e-12
            @test uc.Omega_ce ≈ qe * uc.B0 / (me * cc)      rtol = 1e-12
            @test uc.Omega_ci ≈ qe * uc.B0 / (m_i * cc)     rtol = 1e-12
            @test uc.omegaB0  == uc.Omega_ce

            # lengths
            @test uc.d_e       ≈ cc / uc.omega_pe            rtol = 1e-12
            @test uc.d_i       ≈ cc / uc.omega_pi            rtol = 1e-12
            @test uc.lambda_De ≈ sqrt(Te / (4π * n_e * qe^2)) rtol = 1e-12
            @test uc.rho_e     ≈ uc.v_the / uc.Omega_ce      rtol = 1e-12
            @test uc.rho_i     ≈ uc.v_thi / uc.Omega_ci      rtol = 1e-12

            # speeds
            @test uc.v_the ≈ sqrt(Te / me)                    rtol = 1e-12
            @test uc.v_thi ≈ sqrt(Ti / m_i)                   rtol = 1e-12
            @test uc.v_A   ≈ uc.B0 / sqrt(4π * ρ)             rtol = 1e-12
            @test uc.c_s   ≈ sqrt(γ * (Te + Ti) / (uc.mu_rho * me)) rtol = 1e-12
            @test uc.v_A < cc && uc.c_s < cc

            # beta, from the cgs definition rather than 2 fn T/sigma0
            @test uc.beta_e ≈ 8π * n_e * Te / uc.B0^2      rtol = 1e-12
            @test uc.beta_i ≈ 8π * n_i * Ti / uc.B0^2      rtol = 1e-12

            # temperature: the code unit is m_e c^2
            @test kB * uc.temperature ≈ u.temperature * me * cc^2 rtol = 1e-12
            @test uc.temperature * kB / eV ≈ T * me * cc^2 / eV    rtol = 1e-12

            # the fiducial scales are the skin depth and Larmor radius at v = c
            @test uc.skindepth0 ≈ cc / sqrt(4π * uc.n0 * qe^2 / me) rtol = 1e-12
            @test uc.larmor0    ≈ cc / uc.Omega_ce                  rtol = 1e-12
            @test uc.sigma0     ≈ (uc.skindepth0 / uc.larmor0)^2    rtol = 1e-12

            # macroparticle bookkeeping: q0 is the weight the fiducial cell carries
            @test uc.V0 ≈ u.V0 * uc.length_unit^3     rtol = 1e-12
            @test uc.q0 ≈ uc.n0 * uc.V0 / uc.ppc0 * qe rtol = 1e-12
        end

        @testset "NRL formulary numbers" begin
            T_eV = uc.temperature * kB / eV
            @test uc.omega_pe  ≈ 5.64e4 * sqrt(n_e)            rtol = 1e-3
            @test uc.Omega_ce  ≈ 1.76e7 * uc.B0                rtol = 1e-3
            @test uc.d_e       ≈ 5.31e5 / sqrt(n_e)            rtol = 2e-3
            @test uc.lambda_De ≈ 7.43e2 * sqrt(T_eV / n_e)     rtol = 1e-3
        end

        @testset "all four anchors agree" begin
            for other in (EntityUnitsCGS(u; n0 = uc.n0),
                          EntityUnitsCGS(u; B0 = uc.B0),
                          EntityUnitsCGS(u; length_unit = uc.length_unit))
                for f in fieldnames(EntityUnitsCGS)
                    a, b = getfield(uc, f), getfield(other, f)
                    @test a isa Int ? a == b : isapprox(a, b; rtol = 1e-12, nans = true)
                end
            end
        end

        @testset "code -> cgs conversion factors" begin
            @test uc.length_unit  ≈ uc.d_i / u.d_i        rtol = 1e-12
            @test uc.time_unit    ≈ uc.length_unit / cc   rtol = 1e-12
            @test uc.speed_unit   == cc
            @test uc.mass_unit    == me
            @test uc.charge_unit  == qe
            @test uc.energy_unit  ≈ me * cc^2             rtol = 1e-12
            @test uc.Bfield_unit  ≈ uc.B0 / u.B0          rtol = 1e-12
            # a code length/time/field converts to the cgs scale it names
            @test u.rho_i * uc.length_unit ≈ uc.rho_i     rtol = 1e-12
            @test (1 / u.Omega_ci) * uc.time_unit ≈ 1 / uc.Omega_ci rtol = 1e-12
            @test u.v_A * uc.speed_unit ≈ uc.v_A          rtol = 1e-12
        end

        @testset "conversion factors are the reciprocal scales" begin
            @test uc.d_e * uc.x_de              ≈ 1.0 rtol = 1e-12
            @test uc.d_i * uc.x_di              ≈ 1.0 rtol = 1e-12
            @test uc.rho_e * uc.x_rho_e         ≈ 1.0 rtol = 1e-12
            @test uc.rho_i * uc.x_rho_i         ≈ 1.0 rtol = 1e-12
            @test uc.lambda_De * uc.x_lambda_De ≈ 1.0 rtol = 1e-12
            @test uc.dx * uc.x_cells            ≈ 1.0 rtol = 1e-12
            @test uc.t_omega_pe == uc.omega_pe
            @test uc.t_Omega_ci == uc.Omega_ci
            @test (1 / uc.Omega_ci) * uc.t_Omega_ci ≈ 1.0 rtol = 1e-12
            @test uc.v_A * uc.v_vA ≈ 1.0 rtol = 1e-12
            @test uc.c_s * uc.v_cs ≈ 1.0 rtol = 1e-12
        end

        @testset "ratios survive the rescaling" begin
            @test uc.d_i / uc.d_e   ≈ u.d_i / u.d_e   rtol = 1e-12
            @test uc.rho_i / uc.d_i ≈ u.rho_i / u.d_i rtol = 1e-12
            @test uc.c_s / uc.v_A   ≈ u.c_s / u.v_A   rtol = 1e-12
            @test uc.v_A / cc       ≈ u.v_A           rtol = 1e-12
            @test uc.Omega_ce / uc.omega_pe ≈ u.Omega_ce / u.omega_pe rtol = 1e-12
            # the box in d_i is a code-unit number and must not move
            t = TOML.parsefile(path)
            L_code = t["grid"]["extent"][1][2]
            @test (L_code * uc.length_unit) * uc.x_di ≈ L_code * u.x_di rtol = 1e-12
        end

        @testset "mass and charge unit overrides" begin
            # read the unit-mass species as protons: a hybrid-style mapping
            up = EntityUnitsCGS(u; n_e, mass_unit = EntityTools.m_p_cgs)
            mp = EntityTools.m_p_cgs
            @test up.omega_pe ≈ sqrt(4π * n_e * qe^2 / mp) rtol = 1e-12
            @test up.Omega_ce ≈ qe * up.B0 / (mp * cc)     rtol = 1e-12
            @test up.energy_unit ≈ mp * cc^2               rtol = 1e-12
            @test kB * up.temperature ≈ u.temperature * mp * cc^2 rtol = 1e-12
            # heavier unit mass at fixed density: lower omega_pe, longer d_e
            @test up.omega_pe < uc.omega_pe && up.d_e > uc.d_e
            # dimensionless quantities cannot notice
            @test up.beta ≈ uc.beta rtol = 1e-12

            # doubling the charge unit at fixed density doubles omega_pe
            u2 = EntityUnitsCGS(u; n_e, charge_unit = 2qe)
            @test u2.omega_pe ≈ 2 * uc.omega_pe rtol = 1e-12
            @test u2.beta ≈ uc.beta             rtol = 1e-12
        end

        @testset "run directory and parameter file" begin
            # the parameter-file constructor passes n_frac through to EntityUnits
            uf = EntityUnitsCGS(path; n_frac = fn, gamma_ad = γ, n_e)
            @test uf.n_frac == fn
            @test uf.d_i ≈ uc.d_i             rtol = 1e-12
            @test uf.B0  ≈ uc.B0              rtol = 1e-12
            @test uf.temperature ≈ uc.temperature rtol = UTOL
        end

        @testset "temperature-free runs report NaN" begin
            u3 = EntityUnitsCGS(EntityUnits(1.0, 10.0; mass_ratio = 100.0, n_frac = 0.5);
                                n_e = 1.0)
            @test isnan(u3.temperature) && isnan(u3.beta) && isnan(u3.rho_i)
            @test isnan(u3.x_rho_i) && isnan(u3.v_cs)
            @test all(isfinite, (u3.d_e, u3.d_i, u3.v_A, u3.Omega_ci, u3.B0,
                                 u3.n0, u3.q0, u3.length_unit))
        end

        @testset "invalid input is rejected" begin
            @test_throws ErrorException EntityUnitsCGS(u)
            @test_throws ErrorException EntityUnitsCGS(u; n_e, B0 = 1.0)
            @test_throws ErrorException EntityUnitsCGS(u; n_e = -1.0)
            @test_throws ErrorException EntityUnitsCGS(u; n_e, mass_unit = 0.0)
            @test_throws ErrorException EntityUnitsCGS(u; n_e, charge_unit = -1.0)
        end
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
