"""
Master runner for the data-preprocessing pipeline.

Executes each preprocessing script in order. All outputs (cleaned CSVs)
land under ``processed data/``, consumed by the analysis pipeline in
``analysis/run_all.R`` and ``analysis/run_all.py``.

Run with the famine-enso conda environment, from any directory:
    conda run -n famine-enso python "preprocessing codes/run_all.py"

Each script is run from this folder, because the scripts read and write
relative to it (../data, ../processed data).

Order:
    1. process_Famines.py             – famine chronology, region-year panel
    2. process_lungqvist_prices.py    – Ljungqvist & Seim grain prices
    3. process_ljungqvist_yields.py   – grain harvest panel
    4. process_federico_prices.py     – Federico et al. grain prices (SI)
    5. process_fishprices.py          – cod and herring prices
    6. process_fishcatch.py           – cod and herring catches (SI)
"""

import subprocess
import sys
from pathlib import Path

PREPROC_DIR = Path(__file__).resolve().parent

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
        cwd=PREPROC_DIR,
        check=False,
    )
    if result.returncode != 0:
        sys.exit(f"ERROR: {path.name} exited with code {result.returncode}")
    print(f"Done: {path.name}\n")


if __name__ == "__main__":
    for path in steps:
        run_script(path)

    print("Preprocessing pipeline complete. Outputs in processed data/.")
