"""
Master runner – produces all 5 main figures + appendix maps.

Pipeline:
  assemble_figures.py  – Figs 1–5 (internally imports 07_ml_onset_survival.py
                          for ML models, which also saves ML appendix PDFs)
  08_spatial_maps.py   – supplementary teleconnection maps

Run *after* the R pipeline:
    Rscript analysis/run_all.R        # models + CSV + LaTeX exports
    python  analysis/run_all.py       # reads CSVs, runs ML, saves PDFs

NOTE: run with the famine-enso conda env, not system Python:
    C:/Users/emile/anaconda3/envs/famine-enso/python.exe analysis/run_all.py

Outputs land in:
    analysis/output/figures/main/   – main figures
    paper/figures/plots/            – appendix ML + teleconnection maps
"""

import subprocess
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent

# 07_ml_onset_survival.py is NOT listed here – it is imported directly by
# assemble_figures.py (make_fig2 / make_fig4) so ML runs exactly once.
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
