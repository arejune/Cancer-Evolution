"""Repository paths and the number of worker processes.

TT_RESULTS, TT_FIGURES and TT_TABLES redirect the output folders (the tests use this); TT_WORKERS sets the
number of worker processes (default: one fewer than the CPU count). Results do not depend on it.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "threshold_transfer"
CONFIG = ROOT / "config"
RESULTS = Path(os.environ.get("TT_RESULTS", ROOT / "data" / "results")).resolve()
FIGURES = Path(os.environ.get("TT_FIGURES", ROOT / "figures")).resolve()
TABLES = Path(os.environ.get("TT_TABLES", ROOT / "tables")).resolve()
WORKERS = int(os.environ.get("TT_WORKERS", max(1, (os.cpu_count() or 2) - 1)))
