"""
Master runner for the data-preprocessing pipeline.

Executes each notebook in place via ``jupyter nbconvert --execute`` (so the
notebook cells stay authoritative and reviewable) and runs the standalone
Ljungqvist v2 yield preprocessor as a script. All notebook outputs (CSVs)
land under ``processed data/``, consumed by the analysis pipeline in
``analysis/run_all.R`` and ``analysis/run_all.py``.

Run from the repo root with the famine-enso conda env, for example:
    C:/Users/emile/anaconda3/envs/famine-enso/python.exe \
        "preprocessing codes/run_all.py"

Order:
    1. process_Famines.ipynb                 – famine records & region panel
    2. process_lungqvist_prices.ipynb        – Ljungqvist price data
    3. process_ljungqvist_yields_v2.py       – v2 Ljungqvist yield panel
    4. process_federico_prices.ipynb         – Federico wheat price series
    5. process_fishprices.ipynb              – historical fish prices
    6. process_fishcatch.ipynb               – historical fish catches
"""

import subprocess
import sys
from pathlib import Path

PREPROC_DIR = Path(__file__).parent

# (kind, path) with kind ∈ {"nb", "py"}
steps: list[tuple[str, Path]] = [
    ("nb", PREPROC_DIR / "process_Famines.ipynb"),
    ("nb", PREPROC_DIR / "process_lungqvist_prices.ipynb"),
    ("py", PREPROC_DIR / "process_ljungqvist_yields_v2.py"),
    ("nb", PREPROC_DIR / "process_federico_prices.ipynb"),
    ("nb", PREPROC_DIR / "process_fishprices.ipynb"),
    ("nb", PREPROC_DIR / "process_fishcatch.ipynb"),
]


def _banner(msg: str) -> None:
    print("=" * 60)
    print(msg)
    print("=" * 60)


def run_notebook(path: Path) -> None:
    _banner(f"Executing notebook: {path.name}")
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
    for kind, path in steps:
        if not path.exists():
            print(f"SKIP: {path.name} not found")
            continue
        if kind == "nb":
            run_notebook(path)
        else:
            run_script(path)

    print("Preprocessing pipeline complete. Outputs in processed data/.")
