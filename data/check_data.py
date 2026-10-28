"""Check the raw results in data/results against data/MD5SUMS.

  python data/check_data.py      prints the files that are missing or changed; exit status 1 if any
"""
import hashlib, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def check():
    """Number of files listed in MD5SUMS that are missing or whose MD5 differs."""
    rows = [r for r in (HERE / "MD5SUMS").read_text().split("\n") if r]
    bad = 0
    for row in rows:
        want, rel = row.split("  ", 1); p = ROOT / rel
        if not p.exists() or hashlib.md5(p.read_bytes()).hexdigest() != want:
            bad += 1; print("missing or changed:", rel)
    print(f"{len(rows) - bad} of {len(rows)} files match data/MD5SUMS")
    return bad


if __name__ == "__main__":
    sys.exit(1 if check() else 0)
