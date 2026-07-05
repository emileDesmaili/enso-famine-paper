"""
Master runner – produces main figures, SI/appendix maps, and ED figures.

Pipeline (Python side):
  assemble_figures.py         – Figs 1–5 main (internally imports
                                  07_ml_onset_survival.py for ML models,
                                  which also saves ML appendix PDFs)
  08_spatial_maps.py          – SI teleconnection map (NINO3.4)
  analysis/mechanisms.ipynb   – ED mechanism maps (figED_mechanisms_maps.pdf)
  analysis/maps_data_coverage.ipynb – SI data-coverage maps

The following legacy scripts are kept in the code base and executed as part
of the pipeline, but their outputs now land in
analysis/output/figures/investigation/ rather than the ED folder:
  09_enso_precip_lags.py
  33_fig_v2_onset_panels.py

Run *after* the R pipeline:
    Rscript analysis/run_all.R        # models + CSV + LaTeX + ED Rmds
    python  analysis/run_all.py       # reads CSVs, runs ML, saves PDFs

NOTE: run with the famine-enso conda env, not system Python:
    C:/Users/emile/anaconda3/envs/famine-enso/python.exe analysis/run_all.py

Outputs land in:
    analysis/output/figures/main/           – main figures
    analysis/output/figures/appendix/       – SI figures + ML plots
    analysis/output/figures/extended data/  – ED figures (mechanism, volcanic,
                                              NAO/JSL/ENSO from R side)
    analysis/output/figures/investigation/  – demoted legacy figures
"""

import subprocess
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent
ROOT_DIR     = ANALYSIS_DIR.parent

# .py scripts run via subprocess with the current interpreter.
# 07_ml_onset_survival.py is NOT listed here – it is imported directly by
# assemble_figures.py (make_fig2 / make_fig4) so ML runs exactly once.
scripts = [
    "assemble_figures.py",         # Figs 1-5 (main)
    "08_spatial_maps.py",          # SI NINO3.4 teleconnection map
    "09_enso_precip_lags.py",      # demoted → investigation/
    "33_fig_v2_onset_panels.py",   # demoted → investigation/
]

# .ipynb notebooks executed in-place via jupyter nbconvert.
notebooks = [
    ANALYSIS_DIR / "mechanisms.ipynb",           # ED mechanism maps
    ANALYSIS_DIR / "maps_data_coverage.ipynb",   # SI coverage maps
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


def run_notebook(path: Path) -> None:
    sep = "=" * 60
    print(sep)
    print(f"Executing notebook: {path.name}")
    print(sep)
    result = subprocess.run(
        [sys.executable, "-m", "jupyter", "nbconvert",
         "--to", "notebook", "--execute", "--inplace",
         str(path)],
        check=False,
    )
    if result.returncode != 0:
        print(f"WARNING: {path.name} exited with code {result.returncode}")
    else:
        print(f"Done: {path.name}\n")


if __name__ == "__main__":
    for script in scripts:
        run_script(ANALYSIS_DIR / script)

    for nb in notebooks:
        run_notebook(nb)

    print("All Python scripts complete.")
