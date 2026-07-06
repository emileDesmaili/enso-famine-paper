"""
Master runner for the data-preprocessing pipeline.

Executes each preprocessing script in order. All outputs (cleaned CSVs)
land under ``processed data/``, consumed by the analysis pipeline in
``analysis/run_all.R`` and ``analysis/run_all.py``.

Run from the repo root with the famine-enso conda env:
    conda run -n famine-enso python "preprocessing codes/run_all.py"
    # or (Windows absolute path):
    C:/Users/emile/anaconda3/envs/famine-enso/python.exe \
        "preprocessing codes/run_all.py"

Order:
    1. process_Famines.py             – famine records & region panel
    2. process_lungqvist_prices.py    – Ljungqvist price data
    3. process_ljungqvist_yields.py   – Ljungqvist v2 yield panel
    4. process_federico_prices.py     – Federico wheat price series
    5. process_fishprices.py          – historical fish prices
    6. process_fishcatch.py           – historical fish catches
"""

import subprocess
import sys
from pathlib import Path

PREPROC_DIR = Path(__file__).parent

steps: list[Path] = [
    PREPROC_DIR / "process_Famines.py",
    PREPROC_DIR / "process_lungqvist_prices.py",
    PREPROC_DIR / "process_ljungqvist_yields.py",
    PREPROC_DIR / "process_federico_prices.py",
    PREPROC_DIR / "process_fishprices.py",
    PREPROC_DIR / "process_fishcatch.py",
]


def _banner(msg: str) -> None:
    print("=" * 60)
    print(msg)
    print("=" * 60)


def run_script(path: Path) -> None:
    _banner(f"Running script: {path.name}")
    result = subprocess.run(
        [sys.executable, str(path)],
        check=False,
    )
    if result.returncode != 0:
        print(f"WARNING: {path.name} exited with code {result.returncode}")
    else:
        print(f"Done: {path.name}\n")


if __name__ == "__main__":
    for path in steps:
        if not path.exists():
            print(f"SKIP: {path.name} not found")
            continue
        run_script(path)

    print("Preprocessing pipeline complete. Outputs in processed data/.")
