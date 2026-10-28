"""Design of the main threshold-transfer run: intervals, tolerance, seeds and the grid of multipliers m.

The threshold is m * N*(tau). Multipliers with m >= 1.2 N0 / N*(tau) are dropped, since the threshold would
then reach the progression limit and the policy would never treat.
"""
from threshold_transfer.model.ode import N0

TAUS = (7, 14, 28, 42, 56, 84)
Q = 0.05                 # tolerance on p_prem; arbitrary, no clinical basis
T_END = 750              # run length, days
R = 500                  # seeds per cell
BLOCK_ABM = 1_200_000    # first seed of the main run
CORE = (0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10)
EXT = {7: (), 14: (), 28: (1.15, 1.20),
       42: (1.15, 1.20, 1.25, 1.30, 1.35),
       56: (1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50),
       84: (1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90)}


def m_max(nstar):
    return 1.2 * N0 / nstar


def grid(nstar_by_tau):
    """Multipliers per tau: CORE plus EXT[tau], below m_max(tau)."""
    g = {}
    for t in TAUS:
        mm = m_max(nstar_by_tau[t])
        g[t] = tuple(m for m in sorted(set(CORE) | set(EXT[t])) if m < mm)
    return g
