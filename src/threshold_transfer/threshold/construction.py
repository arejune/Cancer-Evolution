"""Initial tumour and the size-only treatment policy in the agent-based model."""
import numpy as np

from threshold_transfer.model.ode import L, N_MAX, PROG_MUL
from threshold_transfer.model.simulator import TumorSimulator, SimulatorParams

R_IC = 31        # radius of the initial disk (3001 sites)
RFRAC = 0.01     # resistant fraction at day 0
DOSE = 1.0
SENSITIVE, RA, EMPTY = 1, 2, 0


def build(seed):
    """Simulator at day 0: a filled disk of radius R_IC with RFRAC of its cells made resistant, no mutation.

    The resistant cells are chosen with a separate generator (seed + 313,000,000) so the simulator's own
    random stream is untouched. Returns (simulator, cells, resistant cells).
    """
    sim = TumorSimulator(SimulatorParams(grid_size=(L, L), mu_a=0.0, mu_escape=0.0), seed=int(seed))
    sim.initialize_tumor(cluster_radius=R_IC, n_Ra_init=0, n_escape_init=0)
    g = sim.grid
    idx = np.argwhere(g == SENSITIVE)
    n_tot = len(idx)
    n_R = int(round(RFRAC * n_tot))
    pick = np.random.default_rng(int(seed) + 313_000_000).choice(n_tot, size=n_R, replace=False)
    for k in pick:
        r, cc = idx[k]; g[r, cc] = RA
    return sim, n_tot, n_R


def run(seed, policy, tau=None, thr=None, N0=3001.0, horizon=750):
    """One run from build(seed).

    "AT": every tau days from day 0, treat for the next tau days if N >= thr, otherwise hold.
    "CT": treat every day. The run stops at progression (N >= 1.2 N0), extinction or `horizon` days.
    prem = 1 if progression happens during a hold. Also returns burden (sum of daily N), t_prog (nan if none),
    N_end, N_max and the simulator's division counters `led`.
    """
    sim, _, _ = build(seed)
    prog = PROG_MUL * N0
    burden = 0.0; on = True; prem = 0; t_prog = np.nan
    Ntrace_max = 0.0; below_at_epoch = False
    for k in range(1, horizon + 1):
        if policy == "AT" and (k - 1) % tau == 0:
            V = float(sim.get_state()["true_total_volume"])
            on = (V >= thr)
            below_at_epoch = (V < thr)
        sim.step(u_A=DOSE if on else 0.0, u_B=0.0)
        N = float(sim.get_state()["true_total_volume"])
        burden += N
        if N > Ntrace_max: Ntrace_max = N
        if N >= prog:
            if np.isnan(t_prog): t_prog = float(k)
            if policy == "AT" and below_at_epoch: prem = 1
            break
        if N <= 0: break
    return dict(burden=burden, t_prog=t_prog, prem=prem, N_end=N, N_max=Ntrace_max, led=sim.led.copy())
