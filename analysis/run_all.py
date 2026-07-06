"""
Master runner (Python side) – produces main figures, SI/appendix maps,
ED figures, and rebuttal-R2 investigation outputs.

Paper pipeline
--------------
  assemble_figures.py     Figs 1-5 for the main text.
                          Internally imports ml_onset_survival.py to run
                          ML models and also save appendix ML PDFs.
  spatial_maps.py         SI teleconnection map (NINO3.4).
  mechanisms.py           ED mechanism maps (figED_mechanisms_maps.pdf).
  maps_data_coverage.py   SI data-coverage map + coverage table.

Investigation / rebuttal-R2 pipeline
------------------------------------
Outputs land in analysis/output/figures/investigation/ instead of the ED
folder, keeping the paper build self-contained.
  investigation/enso_precip_lags.py   figED_enso_precip_lags.pdf
  investigation/v2_onset_panels.py    figED_v2_onset_panels.pdf
  investigation/r2_sim_mediator.py    sim_mediator_coeffs.pdf

Usage
-----
Run *after* the R pipeline (which produces the CSV/TeX inputs that
assemble_figures.py consumes):

    Rscript analysis/run_all.R        # models + CSVs + LaTeX tables + ED Rs
    python  analysis/run_all.py       # figures 1-5 + ML + SI + ED + invest.

IMPORTANT: use the `famine-enso` conda env (it ships scikit-survival,
xarray, cartopy, etc. — the base env does *not*):

    conda run -n famine-enso python analysis/run_all.py
    # or, absolute path (Windows example):
    C:/Users/emile/anaconda3/envs/famine-enso/python.exe analysis/run_all.py

Output tree
-----------
    analysis/output/figures/main/            5 main figures (fig{1..5}_combined.pdf)
    analysis/output/figures/appendix/        SI / appendix figures (ML, LP robustness)
    analysis/output/figures/extended data/   ED figures + ED table
    analysis/output/figures/investigation/   rebuttal-only figures (not in paper)
    analysis/output/tables/                  LaTeX tables
    analysis/output/data/                    tidy CSVs consumed by figures
"""

import subprocess
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent
ROOT_DIR     = ANALYSIS_DIR.parent

# .py scripts run via subprocess with the current interpreter.
# ml_onset_survival.py is NOT listed here – it is imported directly by
# assemble_figures.py (make_fig2 / make_fig4) so ML runs exactly once.
paper_scripts = [
    "assemble_figures.py",       # Figs 1-5 main + ML appendix plots
    "spatial_maps.py",           # SI NINO3.4 teleconnection map
    "mechanisms.py",             # ED mechanism maps (figED_mechanisms_maps.pdf)
    "maps_data_coverage.py",     # SI data-coverage maps
]

# Investigation / rebuttal-R2 scripts – outputs go to
# analysis/output/figures/investigation/. Kept separate from the paper
# pipeline for a clean replication package.
investigation_scripts = [
    "investigation/enso_precip_lags.py",   # figED_enso_precip_lags.pdf
    "investigation/v2_onset_panels.py",    # figED_v2_onset_panels.pdf
    "investigation/r2_sim_mediator.py",    # sim_mediator_coeffs.pdf
]

# (notebooks are now .py scripts in paper_scripts — no separate notebook runner)


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
    for script in paper_scripts:
        run_script(ANALYSIS_DIR / script)

    for script in investigation_scripts:
        run_script(ANALYSIS_DIR / script)

    print("All Python scripts complete.")
