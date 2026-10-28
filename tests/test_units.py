"""Checks that need no stored data: the simulator, the initial tumour and the threshold N*."""
import numpy as np
import pytest

from threshold_transfer.model.simulator import SimulatorParams, TumorSimulator, warmup_jit
from threshold_transfer.threshold import construction as P, grid as G
from threshold_transfer.model.ode import nstar, nstar_closed


def test_simulator_selftest():
    params = SimulatorParams()
    params.validate()
    warmup_jit(params)
    sim = TumorSimulator(params, seed=42)
    sim.initialize_tumor(cluster_radius=6, n_Ra_init=5)
    for _ in range(30):
        sim.step(u_A=params.MTD, u_B=0.0)
    s = sim.get_state()
    assert 0.0 <= s["true_p_a"] <= 1.0 and sim.time == 30.0
    # with mu_escape = 1, resistant divisions produce ESCAPE cells
    esc = TumorSimulator(SimulatorParams(mu_escape=1.0), seed=1)
    esc.initialize_tumor(cluster_radius=6, n_Ra_init=10)
    for _ in range(10):
        esc.step(u_A=1.0, u_B=1.0)
    assert esc.get_state()["n_Escape"] > 0
    # the same (seed, stream_id) gives the same draws, a different stream_id does not
    a1 = TumorSimulator(params, seed=7, stream_id=0)._rng.random(4).tolist()
    a2 = TumorSimulator(params, seed=7, stream_id=0)._rng.random(4).tolist()
    b = TumorSimulator(params, seed=7, stream_id=1)._rng.random(4).tolist()
    assert a1 == a2 and a1 != b


def test_initial_tumour():
    sim, n_tot, n_R = P.build(1_200_000)
    assert (n_tot, n_R) == (3001, 30)
    assert int((sim.grid != 0).sum()) == 3001


def test_nstar_closed_form_matches_root_without_capacity():
    # with a huge lattice the capacity factor is 1, so the numerical root equals the closed form
    for tau in (7.0, 28.0, 84.0):
        assert nstar(tau, nmax=1e15) == pytest.approx(nstar_closed(tau), rel=1e-5)


def test_nstar_decreases_with_tau():
    ns = [nstar(float(t)) for t in G.TAUS]
    assert all(a > b for a, b in zip(ns, ns[1:]))
    assert np.all(np.array(ns) < 1.2 * G.N0)
