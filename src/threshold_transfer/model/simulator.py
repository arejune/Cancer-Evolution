"""Agent-based tumour model on a square lattice with two diffusing drugs.

Site codes: 0 empty, 1 sensitive, 2 resistant to drug A (RA), 3 resistant to both (ESCAPE).
Each drug field follows dC/dt = D lap(C) - (kappa + kappa_cell * occupied) C + kappa u(t), solved by
explicit finite differences. Each day every cell, in random order, either dies (rate base_death + kill),
divides into a random empty Moore neighbour, or does nothing; with no empty neighbour it cannot divide.
The paper uses drug A only and no mutation.

`led` counts division attempts (row 0 sensitive, row 1 resistant; columns: attempts, succeeded, blocked,
blocked slots at the lattice edge, held by sensitive cells, held by resistant cells). Counting draws no
random numbers, so it does not change any run.
"""
from dataclasses import dataclass

import numpy as np
from numba import njit, prange

EMPTY, SENSITIVE, RA, ESCAPE = np.int8(0), np.int8(1), np.int8(2), np.int8(3)
_DOMAIN_TAG = 0x02A5A1E0  # mixed into every seed


def _make_rng(seed, stream_id, domain_tag):
    """Independent generator for (seed, stream_id); a Generator or SeedSequence is used as given."""
    if isinstance(seed, np.random.Generator):
        return seed
    if isinstance(seed, np.random.SeedSequence):
        return np.random.default_rng(seed)
    return np.random.default_rng(np.random.SeedSequence([int(seed), domain_tag, int(stream_id)]))


@dataclass
class SimulatorParams:
    """Rates are per day; D is in sites^2 per day."""
    D_A: float = 0.10
    kappa_A: float = 0.10
    kappa_cell_A: float = 0.05
    D_B: float = 0.10
    kappa_B: float = 0.10
    kappa_cell_B: float = 0.05
    MTD: float = 1.0
    r_S: float = 0.040
    r_Ra: float = 0.035
    r_escape: float = 0.025
    base_death: float = 0.005
    d_A_S: float = 0.15     # kill of sensitive cells by drug A
    d_B_S: float = 0.10
    d_B_Ra: float = 0.25
    mu_a: float = 1e-4      # S -> RA per division
    mu_escape: float = 1e-5  # RA -> ESCAPE per division
    dt: float = 1.0
    grid_size: tuple = (100, 100)

    def validate(self):
        """Raise if dt breaks the stability bound of either drug field; return the bound."""
        dt_cfl = min(1.0 / (4.0 * self.D_A + self.kappa_A + self.kappa_cell_A),
                     1.0 / (4.0 * self.D_B + self.kappa_B + self.kappa_cell_B))
        if self.dt > dt_cfl + 1e-12:
            raise ValueError(f"dt={self.dt:.6f} exceeds the stability bound {dt_cfl:.6f}")
        return dt_cfl

    def to_array(self):
        # [0] dt [1] r_S [2] r_Ra [3] base_death [4] d_A_S [5] d_B_S [6] d_B_Ra [7] mu_a [8] r_escape [9] mu_escape
        return np.array([self.dt, self.r_S, self.r_Ra, self.base_death, self.d_A_S, self.d_B_S,
                         self.d_B_Ra, self.mu_a, self.r_escape, self.mu_escape], dtype=np.float64)


@njit(fastmath=True, parallel=True)
def _pde_step(C, C_new, grid, D, kappa, kappa_cell, u_input, dt_pde):
    """One explicit diffusion-uptake step with zero-flux boundaries, written into C_new."""
    Nx, Ny = C.shape
    source = kappa * u_input
    for i in prange(Nx):
        im = i - 1 if i > 0 else 1
        ip = i + 1 if i < Nx - 1 else Nx - 2
        for j in range(Ny):
            jm = j - 1 if j > 0 else 1
            jp = j + 1 if j < Ny - 1 else Ny - 2
            laplacian = (C[im, j] + C[ip, j] + C[i, jm] + C[i, jp] - 4.0 * C[i, j])
            kappa_eff = kappa + kappa_cell * (1.0 if grid[i, j] > 0 else 0.0)
            C_new[i, j] = C[i, j] + dt_pde * (D * laplacian - kappa_eff * C[i, j] + source)


@njit(fastmath=True)
def _abm_step(grid, C_A, C_B, p, occ_i, occ_j, rng, led):
    """One day of deaths and divisions, cells visited in random order. p is SimulatorParams.to_array()."""
    Nx, Ny = grid.shape
    dt = p[0]

    n_occ = 0
    for i in range(Nx):
        for j in range(Ny):
            if grid[i, j] != 0:
                occ_i[n_occ] = i
                occ_j[n_occ] = j
                n_occ += 1
    if n_occ == 0:
        return

    for k in range(n_occ - 1, 0, -1):  # Fisher-Yates shuffle
        r_idx = int(rng.random() * (k + 1))
        tmp_i = occ_i[k]; occ_i[k] = occ_i[r_idx]; occ_i[r_idx] = tmp_i
        tmp_j = occ_j[k]; occ_j[k] = occ_j[r_idx]; occ_j[r_idx] = tmp_j

    di = np.array([-1, -1, -1, 0, 0, 1, 1, 1], dtype=np.int32)
    dj = np.array([-1, 0, 1, -1, 1, -1, 0, 1], dtype=np.int32)
    en_i = np.empty(8, dtype=np.int32)
    en_j = np.empty(8, dtype=np.int32)

    for idx in range(n_occ):
        i = occ_i[idx]; j = occ_j[idx]
        ct = grid[i, j]
        if ct == 0:
            continue

        ca = C_A[i, j]; cb = C_B[i, j]
        lam_death = p[3]
        if ct == 1:
            lam_death += p[4] * ca + p[5] * cb
        elif ct == 2:
            lam_death += p[6] * cb

        if ct == 1:
            r_c = p[1]
        elif ct == 2:
            r_c = p[2]
        else:
            r_c = p[8]
        lam_total = lam_death + r_c

        p_event = 1.0 - np.exp(-lam_total * dt)
        r1 = rng.random()
        if r1 >= p_event:
            continue

        r2 = rng.random()
        if r2 < lam_death / lam_total:
            grid[i, j] = 0
            continue

        n_empty = 0
        nb_wall = 0.0; nb_S = 0.0; nb_R = 0.0
        for q in range(8):
            ni = i + di[q]; nj = j + dj[q]
            if ni < 0 or ni >= Nx or nj < 0 or nj >= Ny:
                nb_wall += 1.0
            elif grid[ni, nj] == 0:
                en_i[n_empty] = ni; en_j[n_empty] = nj
                n_empty += 1
            elif grid[ni, nj] == 1:
                nb_S += 1.0
            else:
                nb_R += 1.0
        lr = 0 if ct == 1 else 1
        led[lr, 0] += 1.0
        if n_empty > 0:
            led[lr, 1] += 1.0
        else:
            led[lr, 2] += 1.0
            led[lr, 3] += nb_wall
            led[lr, 4] += nb_S
            led[lr, 5] += nb_R
        if n_empty == 0:
            continue

        slot = int(rng.random() * n_empty)
        ti, tj = en_i[slot], en_j[slot]

        daughter = ct
        if ct == 1 and rng.random() < p[7]:
            daughter = 2
        elif ct == 2 and rng.random() < p[9]:
            daughter = 3
        grid[ti, tj] = daughter


class TumorSimulator:
    """Lattice tumour with two drug fields. `step(u_A, u_B)` advances one day and returns `get_state()`."""

    def __init__(self, params, seed=0, stream_id=0):
        params.validate()
        self.params = params
        self.Nx, self.Ny = params.grid_size
        self._p_arr = params.to_array()

        dt_cfl_A = 1.0 / (4.0 * params.D_A + params.kappa_A + params.kappa_cell_A)
        dt_cfl_B = 1.0 / (4.0 * params.D_B + params.kappa_B + params.kappa_cell_B)
        self._n_pde = max(1, int(np.ceil(params.dt / min(dt_cfl_A, dt_cfl_B))))
        self._dt_pde = params.dt / self._n_pde

        self.grid = np.zeros((self.Nx, self.Ny), dtype=np.int8)
        self.C_A = np.zeros((self.Nx, self.Ny), dtype=np.float64)
        self.C_B = np.zeros((self.Nx, self.Ny), dtype=np.float64)
        self.time = 0.0
        self.led = np.zeros((2, 6), dtype=np.float64)

        self._C_A_buf = np.zeros((self.Nx, self.Ny), dtype=np.float64)
        self._C_B_buf = np.zeros((self.Nx, self.Ny), dtype=np.float64)
        self._occ_i = np.empty(self.Nx * self.Ny, dtype=np.int32)
        self._occ_j = np.empty(self.Nx * self.Ny, dtype=np.int32)
        self._rng = _make_rng(seed, stream_id, _DOMAIN_TAG)

    def initialize_tumor(self, cluster_radius=6, n_Ra_init=0, n_escape_init=0):
        """Filled disk of sensitive cells at the centre, plus resistant cells on random empty sites."""
        self.grid[:] = 0
        self.C_A[:] = 0.0
        self.C_B[:] = 0.0
        self._C_A_buf[:] = 0.0
        self._C_B_buf[:] = 0.0
        self.time = 0.0

        cx, cy = self.Nx // 2, self.Ny // 2
        for i in range(self.Nx):
            for j in range(self.Ny):
                if (i - cx) ** 2 + (j - cy) ** 2 <= cluster_radius ** 2:
                    self.grid[i, j] = SENSITIVE

        for cell_type, n_seed in ((RA, n_Ra_init), (ESCAPE, n_escape_init)):
            if n_seed <= 0:
                continue
            empty_idx = np.argwhere(self.grid == EMPTY)
            if len(empty_idx) == 0:
                break
            n_seed = min(n_seed, len(empty_idx))
            chosen = self._rng.choice(len(empty_idx), size=n_seed, replace=False)
            for k in chosen:
                r, c = empty_idx[k]
                self.grid[r, c] = cell_type

    def step(self, u_A=0.0, u_B=0.0):
        u_A = float(np.clip(u_A, 0.0, self.params.MTD))
        u_B = float(np.clip(u_B, 0.0, self.params.MTD))
        p = self.params
        for _ in range(self._n_pde):
            _pde_step(self.C_A, self._C_A_buf, self.grid, p.D_A, p.kappa_A, p.kappa_cell_A, u_A, self._dt_pde)
            self.C_A, self._C_A_buf = self._C_A_buf, self.C_A
            _pde_step(self.C_B, self._C_B_buf, self.grid, p.D_B, p.kappa_B, p.kappa_cell_B, u_B, self._dt_pde)
            self.C_B, self._C_B_buf = self._C_B_buf, self.C_B
        _abm_step(self.grid, self.C_A, self.C_B, self._p_arr, self._occ_i, self._occ_j, self._rng, self.led)
        self.time += p.dt
        return self.get_state()

    def get_state(self):
        n_S = int(np.sum(self.grid == SENSITIVE))
        n_Ra = int(np.sum(self.grid == RA))
        n_Escape = int(np.sum(self.grid == ESCAPE))
        V = n_S + n_Ra + n_Escape
        return {
            "n_S": n_S, "n_Ra": n_Ra, "n_Escape": n_Escape,
            "true_total_volume": V,
            "true_p_a": n_Ra / V if V > 0 else 0.0,
            "true_resistant_fraction": (V - n_S) / V if V > 0 else 0.0,
            "true_biomarker": float(n_Ra),
            "time": self.time,
        }


def warmup_jit(params=None):
    """Compile both kernels on a small lattice so that later timings exclude compilation."""
    if params is None:
        params = SimulatorParams(grid_size=(8, 8), dt=0.5)
    grid = np.zeros((8, 8), dtype=np.int8)
    grid[3:5, 3:5] = SENSITIVE
    C = np.full((8, 8), 0.1, dtype=np.float64)
    occ = np.empty(64, dtype=np.int32)
    _pde_step(C, np.empty((8, 8)), grid, params.D_A, params.kappa_A, params.kappa_cell_A, 0.5, 0.1)
    _abm_step(grid, C, C, params.to_array(), occ, occ.copy(), np.random.default_rng(0), np.zeros((2, 6)))
