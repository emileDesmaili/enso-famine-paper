"""
Master runner – executes assemble_figures.py, which produces all
5 combined Nature-style figures (Figs 1–5).

Run *after* the R pipeline:
    Rscript analysis/run_all.R        # models + CSV + LaTeX exports
    python  analysis/run_all.py       # reads CSVs, runs ML, saves PDFs

Outputs land in:
    analysis/output/figures/main/
"""

import subprocess
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent

scripts = [
    "assemble_figures.py",
    "08_spatial_maps.py",
]


def run_script(path: Path) -> None:
    sep = "=" * 60
    print(sep)
    print(f"Running: {path.name}")
    print(sep)
    result = subprocess.run(
        [sys.executable, str(path)],
        check=False,
    )
    if result.returncode != 0:
        print(f"WARNING: {path.name} exited with code {result.returncode}")
    else:
        print(f"Done: {path.name}\n")


if __name__ == "__main__":
    for script in scripts:
        run_script(ANALYSIS_DIR / script)

    print("All Python scripts complete.")
