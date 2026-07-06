# El Niño Amplified Food Insecurity in Early Modern Europe

Replication code and data for **"El Niño amplified food insecurity in early modern Europe"** (Esmaili et al., 2026).

The pipeline reproduces every figure and table in the main text, the Extended Data, and the Supplementary Information from the raw sources on disk.

---

## Repository layout

```
enso-famine-paper/
├── data/                    Raw external data (ENSO reconstruction, OWDA, ModE-RA,
│                            volcanic forcing, famine chronologies, prices, yields, ...)
├── processed data/          Cleaned analysis-ready CSVs produced by preprocessing
├── preprocessing codes/     One script per raw-data source
│   └── run_all.py             ← runs all preprocessing steps in order
├── analysis/                Statistical analysis + figure/table generation
│   ├── run_all.R              ← R pipeline (models, CSVs, LaTeX tables)
│   ├── run_all.py             ← Python pipeline (figs 1-5, ML, ED maps)
│   ├── fig{2..5}_*.R          ← per-figure R analyses (CSV exports)
│   ├── ed_*.R                 ← R scripts producing Extended Data figures
│   ├── si_nonlinearity.R      ← SI non-linearity check
│   ├── tables.R               ← SI LaTeX tables
│   ├── assemble_figures.py    ← builds fig1_combined.pdf … fig5_combined.pdf
│   ├── ml_onset_survival.py   ← ML pipeline for fig 2 (C-E) and fig 4 (C-D)
│   ├── spatial_maps.py        ← SI NINO3.4 teleconnection map
│   ├── mechanisms.py          ← ED mechanism maps
│   └── maps_data_coverage.py  ← SI data-coverage map + coverage table
├── emileRegs.R              Shared R helpers (LP panel + interaction utilities)
├── environment.yml          Conda spec for the `famine-enso` env
└── README.md
```

---

## How to replicate

The full pipeline runs in three sequential steps.

### 0. Environment

```bash
git clone https://github.com/emileDesmaili/enso-famine-paper.git
cd enso-famine-paper
conda env create -f environment.yml
conda activate famine-enso
```

You also need a working R installation with the CRAN packages listed at the top of the R scripts (`dplyr`, `fixest`, `ggplot2`, `cowplot`, `scales`, `patchwork`, `broom`, `modelsummary`, `haven`, `purrr`, `furrr`, `future`, `readxl`, `glue`, `stringr`, `sf`).

### 1. Preprocess raw data (optional if `processed data/*.csv` is already populated)

```bash
python "preprocessing codes/run_all.py"
```

Runs each preprocessing script in order (famines → prices → yields → wheat prices → fish prices → fish catches) and writes the cleaned CSVs into `processed data/`.

### 2. Statistical analysis (R) → CSVs + LaTeX tables

```bash
Rscript analysis/run_all.R
```

Produces:
* `analysis/output/data/*.csv` – tidy data consumed by the Python figure builder.
* `analysis/output/tables/*.tex` – all Supplementary tables.
* `analysis/output/figures/appendix/` – LP-based SI figures (LOO, bootstrap, robustness).
* `analysis/output/figures/extended data/` – `figED_volcanic.pdf`, `figED_NAO_JSL_ENSO.pdf`, `figED_yield_geopartition.pdf`, `famine_starts_reg.tex`.

### 3. Figures (Python) → main + ED + SI PDFs

```bash
# use the famine-enso env (scikit-survival, xarray, cartopy, ...)
conda run -n famine-enso python analysis/run_all.py
# on Windows, equivalently:
#   C:/Users/emile/anaconda3/envs/famine-enso/python.exe analysis/run_all.py
```

Produces:
* `analysis/output/figures/main/fig{1..5}_combined.pdf` – the five main figures.
* `analysis/output/figures/appendix/MLOnset_*.pdf`, `MLSurvival_*.pdf` – ML SI plots.
* `analysis/output/figures/extended data/figED_mechanisms_maps.pdf` – ED mechanism composite.
* `analysis/output/figures/appendix/map_data_coverage.pdf`, `teleconnections_nino34.pdf` – SI maps.

---

## Common gotchas

* **Correct Python env.** `assemble_figures.py` imports `ml_onset_survival.py` which requires **scikit-survival**. That lives in the `famine-enso` env; the base Anaconda env does not have it. Running `run_all.py` with the wrong interpreter silently skips fig 2–5 (only fig 1 gets produced) — always launch from `famine-enso`.
* **Path independence.** Every analysis script uses `SCRIPT_DIR` (R) or `Path(__file__)` (Python) to locate the repo root, so you can invoke them from any working directory.
* **External datasets.** The three big NetCDFs (OWDA, ModE-RA temperature, ModE-RA precipitation) are not in this repository — download them from their public sources (see below) into `data/`.

## Data sources

* **ENSO reconstruction (NINO3.4)**: Cook et al. 2024, tree-ring-based (`cook2024-*.txt` in `data/`).
* **OWDA (Old World Drought Atlas)**: <https://www.ldeo.columbia.edu/res/fac/trl/owda> – `owda.nc` in `data/`.
* **ModE-RA (500-year climate reanalysis)**: <https://mode-ra.unibe.ch/> – temp/precip NetCDFs in `data/`.
* **Volcanic forcing**: Sigl et al. 2015 (`volcanic_data.xlsx` in `data/`).
* **Famine chronology and price/yield series**: Ljungqvist et al. 2024 (yields), Ljungqvist & Seim 2024 (prices), Alfani (famines) — full list of source references in the manuscript.

---

## Citation

```bibtex
@article{esmaili2026enso,
  title        = {El Ni\~no amplified food insecurity in early modern Europe},
  author       = {Esmaili, Emile and Puma, Michael J. and Ludlow, Francis and
                  Jobbov{\'a}, Eva and Kumar, Janavi and Holm, Poul and
                  Ljungqvist, Fredrik Charpentier and Matthews, John Alphonsus and
                  Dahl, Johannes Rom and Seim, Andrea},
  year         = {2026},
  note         = {in review}
}
```

## Contact

For questions about this research or code, please open an issue or contact:
Emile Esmaili — ee4561@princeton.edu
