"""Well-mixed stochastic model (CTMC): the ODE's rates and crowding fraction with random daily events.

R replicates run in lockstep from one seed. Each day, events ~ Binomial(cells, 1 - exp(-(death + division) dt))
are split into deaths and division attempts by hazard share, and an attempt succeeds with probability phi(N).
pk=None is CTMC-nominal (drug at the constant C_eff); pk=(kappa_pk, k_el) is CTMC-pk, a one-compartment drug
level dC/dt = kappa_pk u - k_el C fitted to the agent-based model's mean concentration.
"""
import numpy as np

from threshold_transfer.model.ode import GROWTH, L, N_MAX, N0, N_R0, PROG_MUL, N_SPECIES, S, RA, ESC, phi_vec


def initial_counts(R):
    c = np.zeros((R, N_SPECIES), dtype=np.int64)
    c[:, S] = int(N0) - N_R0; c[:, RA] = N_R0
    return c


def batch(p, R, block_seed, policy="AT", tau=None, thr=None, nmax=N_MAX,
          horizon_days=750.0, y0=None, pk=None, growth=GROWTH):
    """R replicates from np.random.default_rng(block_seed) under policy "AT" (threshold thr, interval tau) or
    "CT" (always treat). growth is (c, phi_bulk). Returns per-replicate prem, t_prog (nan if none), burden,
    auc, peak and the final counts.
    """
    rng = np.random.default_rng(block_seed)
    dkill = p.kill_matrix(); r = p.growth_rates(); dt = p.dt
    n_steps = int(round(horizon_days / dt))
    prog = PROG_MUL * N0

    counts = initial_counts(R) if y0 is None else np.array(y0, dtype=np.int64)
    on = np.ones(R, dtype=bool)
    below = np.zeros(R, dtype=bool)
    prem = np.zeros(R, dtype=np.int64)
    t_prog = np.full(R, np.nan)
    burden = np.zeros(R); auc = np.zeros(R); peak = counts.sum(1).astype(float)
    C = np.zeros(R)
    if pk is not None:
        kap, kel = pk
        c_inf = kap / kel; decay = np.exp(-kel * dt)
        mean_w = (1.0 - decay) / kel   # weight for the day's mean concentration

    for k in range(1, n_steps + 1):
        N = counts.sum(1).astype(float)
        if policy == "AT":
            if (k - 1) % tau == 0:
                on = (N >= thr); below = (N < thr)
            u = on.astype(float)
            if pk is None:
                kill = np.where(on[:, None], dkill[:, 0] * p.C_eff, 0.0)
            else:
                Cbar = c_inf * u + (C - c_inf * u) * mean_w
                C = c_inf * u + (C - c_inf * u) * decay
                kill = dkill[None, :, 0] * Cbar[:, None]
        else:
            if pk is None:
                kill = np.broadcast_to(dkill[:, 0] * p.C_eff, (R, N_SPECIES))
            else:
                Cbar = c_inf + (C - c_inf) * mean_w
                C = c_inf + (C - c_inf) * decay
                kill = dkill[None, :, 0] * Cbar[:, None]

        lam_death = p.base_death + kill
        lam_total = lam_death + r
        p_event = 1.0 - np.exp(-lam_total * dt)
        p_dge = lam_death / lam_total

        phi = phi_vec(N, nmax, *growth)
        n_event = rng.binomial(counts, p_event)
        n_death = rng.binomial(n_event, p_dge)
        n_birth = rng.binomial(n_event - n_death, phi[:, None])

        new = counts - n_death
        to_ra = rng.binomial(n_birth[:, S], p.mu_a)
        new[:, S] += n_birth[:, S] - to_ra
        new[:, RA] += to_ra
        to_esc = rng.binomial(n_birth[:, RA], p.mu_escape)
        new[:, RA] += n_birth[:, RA] - to_esc
        new[:, ESC] += to_esc
        new[:, ESC] += n_birth[:, ESC]
        counts = np.maximum(new, 0)

        N_new = counts.sum(1).astype(float)
        burden += N_new; auc += 0.5*(N+N_new)*dt; peak = np.maximum(peak, N_new)
        hit = np.isnan(t_prog) & (N_new >= prog)
        if hit.any():
            t_prog[hit] = k
            if policy == "AT":
                prem[hit & below] = 1
    return dict(prem=prem, t_prog=t_prog, burden=burden, auc=auc, peak=peak, counts=counts)
