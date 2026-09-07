"""Device parameters and the pressure/temperature -> (lambda, gamma) map for
study 3.

Everything here is the pre-ionization state of the fill gas and the antenna
calibration; none of it comes from the RMF solver. The point is just to turn
"argon at 10-100 mTorr, room temperature, T_e = 0.5-10 eV" into the two
dimensionless numbers rmf_solver.py actually takes.

Device (as given, not from the papers)
---------------------------------------
    f_RMF = 250 kHz                  ->  w = 2 pi f
    R     = 10 cm                    plasma column radius
    B_w   = 6 G                      vacuum RMF amplitude (no plasma)
            cross-check: coils are 0.163 G/A, vacuum antenna current ~35 A
            -> 5.7 G, consistent with the measured 6 G to 5%.
    with plasma coupled: antenna current drops to ~10 A -> boundary field
            ~1.6-2 G, and the field measured at the chamber centre drops to
            ~1 G. Kept here only as MEASURED_* for comparison with the
            solver's own B_r(1), B_r(0) once penetration is known.

Resistivity
-----------
Spitzer resistivity for a Z=1 (singly ionised) plasma, electron-electron
collisions included (Chen, "Introduction to Plasma Physics and Controlled
Fusion", Eq. 5-70; NRL Plasma Formulary Lorentz-gas value is ~2x higher):

    eta(T_e) = 5.2e-5 * Z * lnLambda / T_e[eV]^1.5      Ohm m

lnLambda from the NRL Plasma Formulary, thermal e-i regime (T_i m_e/m_i <
T_e < 10 Z^2 eV, which covers all of T_e = 0.5-10 eV at Z=1):

    lnLambda = 23 - ln( sqrt(n_e[cm^-3]) * Z * T_e[eV]^-1.5 )

Both n and T_e only enter the solver through eta -> gamma and eta -> lambda;
no magnetic-field-dependent (perpendicular) transport is modelled -- see the
README caveat in this study's write-up.
"""
import numpy as np

# ------------------------------------------------------------------- device
F_RMF = 250e3                      # Hz, antenna drive frequency
W = 2*np.pi*F_RMF                  # rad/s
R = 0.10                           # m, plasma column radius
B_W_G = 6.0                        # G, vacuum RMF amplitude (measured)
B_W = B_W_G*1e-4                   # T

# antenna calibration cross-check (not used in the lambda/gamma map itself)
COIL_GAUSS_PER_AMP = 0.163
I_VACUUM_A = 35.0
I_PLASMA_A = 10.0
B_W_FROM_CALIBRATION_G = COIL_GAUSS_PER_AMP*I_VACUUM_A            # ~5.7 G

# measured, with plasma coupled -- for comparing against the solver's own
# B_r(R), B_r(0) once (lambda, gamma) are known, not inputs to the map below
MEASURED_B_WALL_G = COIL_GAUSS_PER_AMP*I_PLASMA_A                 # ~1.6 G
MEASURED_B_CENTER_G = 1.0

# --------------------------------------------------------------------- gas
T_ROOM_K = 298.15                  # 25 C
K_B = 1.380649e-23                 # J/K
TORR_TO_PA = 133.322
Z_ARGON = 1

# -------------------------------------------------------------------- E&M
MU_0 = 4*np.pi*1e-7
E_CHARGE = 1.602176634e-19


def density_from_pressure(p_mtorr, T_K=T_ROOM_K):
    """Neutral fill density (m^-3) from pressure (mTorr) via the ideal gas law.

    Used as the plasma density n_e, i.e. assuming full single ionisation of
    the pre-breakdown fill -- the standard order-of-magnitude estimate when
    the actual ionisation fraction is not otherwise known.
    """
    p_pa = p_mtorr*1e-3*TORR_TO_PA
    return p_pa/(K_B*T_K)


def coulomb_log_ei(n_m3, Te_eV, Z=Z_ARGON):
    """NRL Plasma Formulary thermal e-i Coulomb log (valid for Te < 10 Z^2 eV)."""
    n_cm3 = n_m3*1e-6
    return 23.0 - np.log(np.sqrt(n_cm3)*Z*Te_eV**-1.5)


def eta_spitzer(n_m3, Te_eV, Z=Z_ARGON):
    """Isotropic Spitzer resistivity (Ohm m), Chen Eq. (5-70)."""
    lnL = coulomb_log_ei(n_m3, Te_eV, Z)
    return 5.2e-5*Z*lnL/Te_eV**1.5


def lam_of(Te_eV, n_m3, R=R, w=W, Z=Z_ARGON):
    """lambda = R sqrt(mu0 w / (2 eta)); depends on n only through lnLambda."""
    eta = eta_spitzer(n_m3, Te_eV, Z)
    return R*np.sqrt(MU_0*w/(2.0*eta))


def gam_of(n_m3, Te_eV, B_w=B_W, Z=Z_ARGON):
    """gamma = B_w/(e n eta)  (= w_ce/nu_ei, README)."""
    eta = eta_spitzer(n_m3, Te_eV, Z)
    return B_w/(E_CHARGE*n_m3*eta)


def skin_time_ms(Te_eV, n_m3, R=R, Z=Z_ARGON):
    """Classical resistive diffusion time mu0 R^2/eta, in ms.

    This is the L/R time of the transverse-field skin effect: the same
    combination as 2 lambda^2/w with w cancelled out, so (unlike lambda and
    gamma individually) it does not depend on the drive frequency at all.
    """
    eta = eta_spitzer(n_m3, Te_eV, Z)
    return MU_0*R**2/eta*1e3


def period_ms(w=W):
    return 2*np.pi/w*1e3


# ---------------------------------------------------------- solver settings
def choose_Nr(lam, ppd=6.0, nr_min=32, nr_max=180):
    """Radial points: ppd cells across the skin depth delta = 1/lambda."""
    return int(min(nr_max, max(nr_min, np.ceil(ppd*lam) + 1)))


def choose_dt(lam, Nr, diff_c=0.3, cap=0.002):
    """Explicit diffusion stability limit dt <= diff_c*lambda^2*dr^2.

    No Hall-term constraint here (study2's other limit, ~5*dr^2): gamma is
    << gamma_c everywhere in this study (see the grid), so there is no stiff
    nonlinear term to worry about, only linear diffusion.
    """
    dr = 1.0/(Nr - 1)
    return min(cap, diff_c*lam**2*dr**2)


def choose_periods(tau_ms, factor=3.0, floor=30, ceil=500):
    """Periods to reach steady state: `factor` times the classical skin time."""
    return int(min(ceil, max(floor, np.ceil(factor*tau_ms/period_ms()))))
