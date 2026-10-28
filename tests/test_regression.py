"""The checks of tests/regression.py under pytest. Uses the results in data/results.
Quick checks run by default; set TT_FULL=1 to add the slower ones."""
import os
import pytest
import regression as R

pytestmark = pytest.mark.skipif(not (R.RES / "p_t_main_raw.npz").exists(), reason="raw results missing from data/results")


@pytest.mark.parametrize("exp", R.SIM_EXPERIMENTS)
def test_simulation_subset(exp):
    r = R.worker(exp)
    n, bad = R.ARCHIVE[exp](r["out"])
    assert n >= 200 and bad == 0
    assert r["sha256"] == R.EXPECTED["subsets"][exp]


@pytest.mark.parametrize("name,fn", R.QUICK_ANALYSES, ids=[n for n, _ in R.QUICK_ANALYSES])
def test_analysis_quick(name, fn):
    assert fn()


@pytest.mark.skipif(os.environ.get("TT_FULL") != "1", reason="set TT_FULL=1")
@pytest.mark.parametrize("name,fn", R.FULL_ANALYSES, ids=[n for n, _ in R.FULL_ANALYSES])
def test_analysis_full(name, fn):
    assert fn()
