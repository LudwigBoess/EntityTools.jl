# =============================================================================
#  Physical constants in cgs-Gaussian units
# =============================================================================
#
#  CODATA 2022 values. Charge is in statcoulomb (esu) and magnetic field in
#  gauss, so the Gaussian forms of the plasma relations hold with no vacuum
#  permittivity:
#
#    omega_p = sqrt(4 pi n q^2 / m)     Omega_B = q B / (m c)
#    v_A     = B / sqrt(4 pi rho)       beta    = 8 pi n kB T / B^2
#
#  The names are not exported; use them as `EntityTools.c_cgs` etc.
# =============================================================================

const c_cgs   = 2.99792458e10       # speed of light [cm/s]
const q_e_cgs = 4.80320471257e-10   # elementary charge [esu]
const m_e_cgs = 9.1093837139e-28    # electron mass [g]
const m_p_cgs = 1.67262192595e-24   # proton mass [g]
const k_B_cgs = 1.380649e-16        # Boltzmann constant [erg/K]
const eV_cgs  = 1.602176634e-12     # electronvolt [erg]
