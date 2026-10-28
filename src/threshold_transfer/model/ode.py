"""Constants, growth law and deterministic ODE shared by the threshold construction, the ODE arm and the CTMC.

Untreated growth is dN/dt = r_S phi(N) N - d N, with the crowding fraction
phi(N) = clip((c N^(-1/2) + phi_bulk) (1 - N / n_max), 0, 1), the fraction of cells with an empty neighbour.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

L = 208                    # lattice edge
N_MAX = float(L * L)       # lattice sites
N0 = 3001.0                # initial tumour size
N_R0 = 30                  # initial resistant cells (1% of N0)
PROG_MUL = 1.2             # progression at N >= 1.2 N0
GROWTH = (3.5340, 0.2712)  # (c, phi_bulk), fitted to drug-free growth of the lattice model
S, RA, ESC = 0, 1, 2       # species order
N_SPECIES = 3


@dataclass
class ODE2Params:
    """Per-cell rates (per day) of the agent-based model, used by the ODE and the CTMC."""
    r_S: float = 0.040
    r_Ra: float = 0.035
    r_escape: float = 0.025
    base_death: float = 0.005
    d_A_S: float = 0.15
    d_B_S: float = 0.10
    d_B_Ra: float = 0.25
    mu_a: float = 1e-4
    mu_escape: float = 1e-5
    C_eff: float = 0.10 / (0.10 + 0.05)   # well-mixed drug level at an occupied site, kappa / (kappa + kappa_cell)
    growth_params: tuple = GROWTH
    dt: float = 1.0

    def growth_rates(self):
        return np.array([self.r_S, self.r_Ra, self.r_escape], dtype=np.float64)

    def kill_matrix(self):
        """Kill coefficients [species, drug]; multiply by the drug concentration."""
        return np.array([[self.d_A_S, self.d_B_S],
                         [0.0, self.d_B_Ra],
                         [0.0, 0.0]], dtype=np.float64)


P = ODE2Params(mu_a=0.0, mu_escape=0.0)


def phi(N, c=GROWTH[0], pb=GROWTH[1], nmax=N_MAX):
    """Crowding fraction for a scalar or array N."""
    return np.clip((c * N ** -0.5 + pb) * np.clip(1.0 - N / nmax, 0.0, 1.0), 0.0, 1.0)


def phi_vec(N, nmax=N_MAX, c=GROWTH[0], pb=GROWTH[1]):
    """phi for an array of totals, 0 where N = 0."""
    out = np.zeros_like(N, dtype=float)
    m = N > 0
    if m.any():
        out[m] = np.clip(phi(N[m], c, pb, nmax), 0.0, 1.0)
    return out


def regrow_time(n_start, c=GROWTH[0], pb=GROWTH[1], nmax=N_MAX, dt=0.02, tmax=30000.0):
    """Days for untreated growth from n_start to PROG_MUL * N0 (RK4, crossing interpolated); nan if never."""
    prog = PROG_MUL * N0
    dN = lambda N: P.r_S * phi(N, c, pb, nmax) * N - P.base_death * N
    N = float(n_start); t = 0.0
    for _ in range(int(tmax / dt)):
        k1 = dN(N); k2 = dN(N + dt/2*k1); k3 = dN(N + dt/2*k2); k4 = dN(N + dt*k3)
        Np = N; N = N + dt/6*(k1 + 2*k2 + 2*k3 + k4); t += dt
        if N >= prog:
            return t - dt * (N - prog) / (N - Np)
    return np.nan


def nstar(tau, c=GROWTH[0], pb=GROWTH[1], nmax=N_MAX):
    """Threshold N*(tau): the size whose untreated regrowth reaches progression in exactly tau days."""
    return brentq(lambda x: regrow_time(x, c, pb, nmax) - tau, 0.5, PROG_MUL * N0 - 1e-6, xtol=1e-7)


def nstar_closed(tau):
    """Closed form of N*(tau) without the capacity factor (Appendix B); for reference only."""
    A, B = P.r_S * GROWTH[0], P.r_S * GROWTH[1] - P.base_death
    return (((A + B * np.sqrt(PROG_MUL * N0)) * np.exp(-(B / 2) * tau) - A) / B) ** 2


def ode_cell(tau, thr, c=GROWTH[0], pb=GROWTH[1], t_end=750):
    """Three-species ODE under the size-only policy: every tau days treat if N >= thr, else hold.

    RK4 with step 0.2 d. Returns prem (1 if progression happens during a hold) and t_prog (nan if none).
    """
    y = np.array([N0 - N_R0, float(N_R0), 0.0]); prog = PROG_MUL * N0
    dk = P.kill_matrix(); r = P.growth_rates(); on = True; below = False
    for k in range(1, t_end + 1):
        if (k - 1) % tau == 0:
            N = y.sum(); on = (N >= thr); below = (N < thr)
        kill = dk[:, 0] * P.C_eff if on else np.zeros(3)

        def dv(v):
            N = max(v.sum(), 1e-12)
            ph = min(max((c * N ** -0.5 + pb) * max(1.0 - N / N_MAX, 0.0), 0.0), 1.0)
            return r * ph * v - (P.base_death + kill) * v
        h = 0.2
        for _ in range(int(1 / h)):
            k1 = dv(y); k2 = dv(y + h/2*k1); k3 = dv(y + h/2*k2); k4 = dv(y + h*k3)
            y = np.maximum(y + h/6*(k1 + 2*k2 + 2*k3 + k4), 0.0)
        if y.sum() >= prog:
            return dict(prem=int(below), t_prog=float(k))
    return dict(prem=0, t_prog=float("nan"))


def ode_boundary(tau, c=GROWTH[0], pb=GROWTH[1], lo_m=0.70, hi_m=2.20):
    """The ODE's own safe threshold: largest threshold with no premature progression (bisection, 34 steps).

    Returns (threshold, t_prog at the upper end of the final bracket).
    """
    ref = nstar(tau, c, pb)
    lo, hi = lo_m * ref, hi_m * ref
    while ode_cell(tau, lo, c, pb)["prem"] == 1: lo *= 0.9
    while ode_cell(tau, hi, c, pb)["prem"] == 0: hi *= 1.05
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ode_cell(tau, mid, c, pb)["prem"] == 0 else (lo, mid)
    return lo, ode_cell(tau, hi, c, pb)["t_prog"]
