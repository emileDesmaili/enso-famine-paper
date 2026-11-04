# Replication package: "El Niño amplified food insecurity in early modern Europe from 1500 to 1800"

Esmaili, E., Puma, M. J., Ludlow, F., Jobbová, E., Kumar, J., Holm, P., Ljungqvist, F. C., Matthews, J. A., Dahl, J. R. & Seim, A. *Nature Food* (2026). https://doi.org/10.1038/s43016-026-01460-8

This repository contains the data and code that reproduce every figure, table and number reported in the main text, the Extended Data and the Supplementary Information of the article.

---

## Overview

The package runs in three stages:

| Stage | Command | What it does | Run time |
|---|---|---|---|
| 1. Preprocessing (optional) | `python "preprocessing codes/run_all.py"` | Builds the analysis datasets in `processed data/` from the raw sources in `data/` | ~1 min |
| 2. Statistical analysis (R) | `Rscript analysis/run_all.R` | Estimates all models; writes estimates, LaTeX tables and Supplementary/Extended Data figures to `analysis/output/` | ~5 min |
| 3. Figures and machine learning (Python) | `python analysis/run_all.py` | Runs the machine-learning models and climate maps; assembles Figs. 1–5 and the remaining Extended Data and Supplementary figures | ~4 min |

The processed datasets are included, so stage 1 can be skipped. Stage 1 regenerates them byte for byte from the raw data (it needs the three large climate files listed under [Data not included](#data-not-included-download-before-running)).

---

## Data availability

All data are from published sources. The raw files used are in `data/`, except three large climate fields that must be downloaded (below). The analysis-ready datasets are in `processed data/`.

### Data included in this repository

| Data | File(s) | Source | Used for |
|---|---|---|---|
| Famine chronology (famine years by region, 1500–1800) | `data/FaminesFoodDisruptions_v2.xlsx` (sheet `Ljungqvist2024_TableA1`) | Ljungqvist, Seim & Collet (2024), *WIREs Clim. Change* 15, e859, Table A1 | Famine onset and duration |
| ENSO reconstruction (NINO3.4, boreal winter) | `data/cook2024-R15-ENSO-Rec-1500-2000.txt`, `data/cook2024-ENSO-latlon.txt` | Cook & Cane (2024), *Geophys. Res. Lett.* 51, e2024GL109759 | All analyses |
| North Atlantic Oscillation reconstruction | `data/nao_reconst.xlsx` | Ortega et al. (2015), *Nature* 523, 71–74 | Controls |
| Jet stream latitude reconstruction | `data/reconstructed EU JSL.xlsx` | Xu et al. (2024), *Nature* 634, 600–608 | Controls |
| Armed conflicts | `data/Conflict-Catalog-18-vars.xlsx`, `data/Brecke-Pre-1400-European-Conflicts.xlsx` | Brecke (1999), Conflict Catalog | Controls |
| Volcanic forcing and ice-core sulfate | `data/volcanic_data.xlsx` | Sigl et al. (2015), *Nature* 523, 543–549 | Extended Data Fig. 4 |
| Grain prices (44 cities) | `data/Grain_prices_series_Ljungqvist_and_Seim.xlsx` | Ljungqvist & Seim (2024), *Eur. Econ. Rev.* 170, 104852 | Fig. 5a |
| Grain prices (alternative) | `data/FSV prices.xlsx` | Federico, Schulze & Volckart (2021), *J. Econ. Hist.* 81, 276–308 | Supplementary Fig. 14a |
| Grain harvests (tithes, yield ratios, yields) | `data/Data_Ljungqvist_et_al_2023.csv`, `data/fcl_yields_new/` | Ljungqvist et al. (2023), *Clim. Past* 19, 2463–2491, and the sources in Supplementary Table 2 | Fig. 3, Extended Data Fig. 2 |
| Fish prices (cod, herring; 8 cities) | `data/Big_DF.csv` | Matthews (2021), PhD thesis, Trinity College Dublin; Allaire & Holm (2022), *Food Hist.* 20, 41–66; Allen & Unger (2019), *Res. Data J. Humanit. Soc. Sci.* 4, 81–90 | Fig. 5b |
| Fish catches (cod, herring) | `data/fish_catch_europe.csv` | Holm et al. (2021), *Fish Fish.*; Norfish data collection, https://doi.org/10.6084/m9.figshare.c.5514351.v1 | Supplementary Fig. 12 |
| Wheat prices and yield ratios (alternative data) | `processed data/wheatprices_dataset.dta`, `processed data/yieldratio_dataset.dta` | Processed versions of Allen (2001) and Slicher van Bath (1963) from the replication files of Waldinger (2022), *J. Polit. Econ.* 130, 2275–2314 | Supplementary Fig. 14b,c |
| Coordinates of price cities and harvest locations; country of each harvest record | `data/coordinates_*.csv`, `data/countries_harvest_records.csv` | This study (geocoded once with OpenStreetMap Nominatim and stored, so the pipeline needs no network access) | Climate extraction |

### Data not included (download before running)

These files are too large for the repository. Put them in `data/` with the file names shown.

| Data | File name in `data/` | Download | Citation |
|---|---|---|---|
| Old World Drought Atlas (JJA scPDSI, 0.5°) | `owda.nc` | https://www1.ncdc.noaa.gov/pub/data/paleo/treering/reconstructions/europe/owda.nc (NOAA study 19419, https://doi.org/10.25921/rjm6-mq74) | Cook et al. (2015), *Sci. Adv.* 1, e1500561 |
| ModE-RA ensemble-mean 2 m temperature, monthly anomalies | `ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc` | https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/ModE/ModE-RA/ensstat (also WDCC, https://doi.org/10.26050/WDCC/ModE-RA_s14203-18501) | Valler et al. (2024), *Sci. Data* 11, 36 |
| ModE-RA ensemble-mean total precipitation, monthly anomalies | `ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc` | as above | as above |

They are needed for preprocessing (stage 1) and for the climate maps (Fig. 3a,b, Extended Data Fig. 3, Supplementary Fig. 7) in stage 3. The R analysis (stage 2) does not need them.

### Analysis datasets (`processed data/`)

| File | Content | Built by |
|---|---|---|
| `famine_region_data.csv` | Region × year panel (9 regions, 1500–1800): famine onset and duration, ENSO, regional climate, NAO, JSL, conflict | `process_Famines.py` |
| `famine_survival.csv` | One row per famine episode (160), with episode-average predictors | `analysis/fig4_survival.R` |
| `price_2023_enso.csv` | City × year grain prices with local climate and teleconnection flags | `process_lungqvist_prices.py` |
| `yield_ljungqvist_v2.csv` | Record × year grain harvests with local climate and teleconnection flags | `process_ljungqvist_yields.py` |
| `fishprice_enso.csv` | City × species × year fish prices | `process_fishprices.py` |
| `federico_prices_enso.csv` | City × year grain prices (Federico et al.) | `process_federico_prices.py` |
| `fishcatch_enso.csv` | Fish catch series × year | `process_fishcatch.py` |

Rights to the original data remain with their authors; please cite the original sources when using them.

---

## Computational requirements

### Software

- **Python 3.10** with the packages pinned in `environment.yml` (numpy 2.2.6, pandas 2.3.3, scipy 1.15.2, scikit-learn 1.7.1, scikit-survival 0.25.0, xarray 2025.6.1, netCDF4 1.7.2, cartopy 0.25.0, matplotlib 3.10.0, seaborn 0.13.2, openpyxl 3.1.5).
- **R 4.5.1** with the packages pinned in `install_R_packages.R` (main ones: fixest 0.14.0 for the linear probability, logit and local-projection models; survival 3.8-3 for the Cox models).

### Hardware and run time

The results were produced on a Windows 11 laptop (Intel Core i7-14700HX, 64 GB RAM). Total run time is about 10 minutes. Most of the time goes to the permutation tests (1,000 draws for famine onset and duration, 500 for prices and harvests) and the nested cross-validation of the survival models. The climate files need about 1.2 GB of disk space.

### Randomness

All random procedures are seeded: the permutation tests (`set.seed(42)` in R, `furrr_options(seed = TRUE)` for parallel draws) and the machine-learning models (`random_state=42`). Results are reproducible on re-runs.

---

## Instructions

```bash
git clone https://github.com/emileDesmaili/enso-famine-paper.git
cd enso-famine-paper

# Python environment
conda env create -f environment.yml
conda activate famine-enso

# R packages (add --exact to install the exact versions used)
Rscript install_R_packages.R

# Download owda.nc and the two ModE-RA files into data/ (see "Data not included")

# 1. Preprocessing (optional: processed data/ is already populated)
python "preprocessing codes/run_all.py"

# 2. Statistical analysis
Rscript analysis/run_all.R

# 3. Machine learning, maps and figure assembly
python analysis/run_all.py
```

Stage 3 reads the estimates written by stage 2, so run them in this order. All scripts locate the repository root themselves and can be launched from any directory. Use the `famine-enso` environment for stage 3: the figures need scikit-survival, which a base Python installation does not have.

The machine-learning results are cached in `analysis/output/data/_ml_*_cache.pkl` after the first run. Delete these files to re-run the models.

---

## Repository structure

```
├── data/                    raw data (see Data availability)
├── processed data/          analysis datasets
├── preprocessing codes/     one script per data source; run_all.py runs them in order
├── analysis/
│   ├── run_all.R            runs the R scripts below
│   ├── fig2_onset.R         famine onset: linear probability and logit models, permutation and leave-one-out tests
│   ├── fig3_harvest.R       grain harvests: local projections, robustness checks
│   ├── fig4_survival.R      famine duration: Cox models, survival curves, parametric and fixed-effects models
│   ├── fig5_prices.R        grain and fish prices, fish catches: local projections, robustness checks
│   ├── ed_nao_jsl_enso.R    Extended Data Fig. 1
│   ├── ed_volcanic.R        Extended Data Fig. 4
│   ├── si_nonlinearity.R    Supplementary Fig. 13
│   ├── tables.R             summary statistics and fixed-effects tables
│   ├── run_all.py           runs the Python scripts below
│   ├── assemble_figures.py  Figs. 1–5 (calls ml_onset_survival.py)
│   ├── ml_onset_survival.py gradient boosting attribution and survival machine learning
│   ├── mechanisms.py        Extended Data Fig. 3
│   ├── spatial_maps.py      Supplementary Fig. 7
│   ├── maps_data_coverage.py Supplementary Fig. 1 and Supplementary Table 8
│   └── output/              created by the pipeline (data/, tables/, figures/)
├── emileRegs.R              local-projection helper functions used by the R scripts
├── environment.yml          Python environment
├── install_R_packages.R     R packages
└── LICENSE                  MIT licence (code)
```

---

## List of exhibits

Output paths are relative to `analysis/output/`. Figures are written as PDF; tables as LaTeX.

### Main text

| Exhibit | Output | Produced by |
|---|---|---|
| Fig. 1 | `figures/main/fig1_combined.pdf` | `assemble_figures.py` |
| Fig. 2 | `figures/main/fig2_combined.pdf` | `fig2_onset.R` (a, b), `ml_onset_survival.py` (c–e), `assemble_figures.py` |
| Fig. 3 | `figures/main/fig3_combined.pdf` | `assemble_figures.py` (a, b), `fig3_harvest.R` (c, d) |
| Fig. 4 | `figures/main/fig4_combined.pdf` | `fig4_survival.R` (a, b), `ml_onset_survival.py` (c, d), `assemble_figures.py` |
| Fig. 5 | `figures/main/fig5_combined.pdf` | `fig5_prices.R`, `assemble_figures.py` |
| Table 1 (classifier grid) | grid in `ml_onset_survival.py` | — |

### Extended Data

| Exhibit | Output | Produced by |
|---|---|---|
| Extended Data Table 1 | `figures/extended data/famine_starts_reg.tex` | `fig2_onset.R` |
| Extended Data Fig. 1 | `figures/extended data/figED_NAO_JSL_ENSO.pdf` | `ed_nao_jsl_enso.R` |
| Extended Data Fig. 2 | `figures/extended data/figED_yield_geopartition.pdf` | `fig3_harvest.R` |
| Extended Data Fig. 3 | `figures/extended data/figED_mechanisms_maps.pdf` | `mechanisms.py` |
| Extended Data Fig. 4 | `figures/extended data/figED_volcanic.pdf` | `ed_volcanic.R` |

### Supplementary Information

| Exhibit | Output | Produced by |
|---|---|---|
| Suppl. Tables 1–2 (data sources) | — (source lists) | — |
| Suppl. Table 3 | `tables/summary_famine.tex` | `fig4_survival.R` |
| Suppl. Tables 4–7 | `tables/summary_price.tex`, `summary_yield_TI.tex`, `summary_yield_YR.tex`, `summary_fishprice.tex` | `tables.R` |
| Suppl. Table 8 | `tables/data_coverage_table.tex` | `maps_data_coverage.py` |
| Suppl. Table 9 | `tables/dyn_survival_models_reg.tex` | `fig4_survival.R` |
| Suppl. Table 10 (survival grids) | grids in `ml_onset_survival.py` | — |
| Suppl. Table 11 (robustness summary) | — (text) | — |
| Suppl. Table 12 | `tables/survival_models_parametric.tex` | `fig4_survival.R` |
| Suppl. Tables 13–14 | `tables/famine_duration_fe_avgnino.tex`, `famine_duration_fe_maxnino.tex` | `fig4_survival.R` |
| Suppl. Tables 15–17 | `tables/Ljungqvist_ENSO_price_main.tex`, `ENSO_fishprice_main.tex`, `FE_ENSO_yieldPDSI.tex` | `tables.R` |
| Suppl. Fig. 1 | `figures/appendix/map_data_coverage.pdf` | `maps_data_coverage.py` |
| Suppl. Figs. 2–5 | `figures/appendix/MLOnset_<predictor set>_*.pdf` | `ml_onset_survival.py` |
| Suppl. Fig. 6 | `figures/appendix/figA_ML_onset_confusion_roc_All_features.pdf` | `ml_onset_survival.py` |
| Suppl. Fig. 7 | `figures/appendix/teleconnections_nino34.pdf` | `spatial_maps.py` |
| Suppl. Fig. 8 | `figures/appendix/IRF_ENSO_yield_WR_teleco.pdf`, `figA_irf_yield_WR_AMJJ_P.pdf` | `fig3_harvest.R` |
| Suppl. Figs. 9–10 | `figures/appendix/figA_surv_strata.pdf`, `figA_surv_inter.pdf` | `fig4_survival.R` |
| Suppl. Fig. 11 | `figures/appendix/MLSurvival_*.pdf` | `ml_onset_survival.py` |
| Suppl. Fig. 12 | `figures/appendix/figA_irf_catches.pdf` | `fig5_prices.R` |
| Suppl. Fig. 13 | `figures/appendix/FigA_nonlin.pdf` | `si_nonlinearity.R` |
| Suppl. Fig. 14 | `figures/appendix/figA_irf_FSV.pdf`, `figA_irf_allen.pdf`, `figA_irf_yield63.pdf` | `fig5_prices.R` (a, b), `fig3_harvest.R` (c) |
| Suppl. Fig. 15 | `figures/appendix/figA_permutation_cox.pdf`, `figA_permutation_onset.pdf` | `fig4_survival.R`, `fig2_onset.R` |
| Suppl. Fig. 16 | `figures/appendix/figA_irf_bootstrap_{price,WR,fishprice}.pdf` | `fig5_prices.R`, `fig3_harvest.R` |
| Suppl. Fig. 17 | `figures/appendix/figA_loo_onset.pdf`, `figA_loo_duration.pdf` | `fig2_onset.R`, `fig4_survival.R` |
| Suppl. Fig. 18 | `figures/appendix/figA_irf_loo_{loc,20y}_{price,WR,fishprice}.pdf` | `fig5_prices.R`, `fig3_harvest.R` |
| Suppl. Fig. 19 | `figures/appendix/figA_irf_serobust_{price,fishprice,WR}.pdf` | `fig5_prices.R`, `fig3_harvest.R` |

The numbers quoted in the text come from the estimate files in `analysis/output/data/` (for example `fig5A_irf_grain_price.csv` for the grain price responses and `fig4A_cox_coefs.csv` for the hazard ratios).

---

## Methodological notes

- Standard errors are Driscoll–Kraay (`vcov = DK ~ Year` in fixest) unless stated otherwise. The spatial HAC (Conley) robustness check in Supplementary Fig. 19 uses a 300 km cutoff for prices and 500 km for harvests.
- Local projections use 3 lags of NINO3.4 as controls, decade fixed effects and record- or city-specific fixed effects with linear trends. Horizons are 0–3 years for harvests and 0–10 years for prices.
- Teleconnected locations are those whose JJA scPDSI is correlated with NINO3.4 at the 10% level (computed in preprocessing, column `teleco_PDSI_10`).

---

## Licence and citation

The code is released under the MIT licence (`LICENSE`). The data remain subject to the terms of their original sources (see Data availability).

```bibtex
@article{esmaili2026enso,
  title   = {El Ni{\~n}o amplified food insecurity in early modern Europe from 1500 to 1800},
  author  = {Esmaili, Emile and Puma, Michael J. and Ludlow, Francis and
             Jobbov{\'a}, Eva and Kumar, Janavi and Holm, Poul and
             Ljungqvist, Fredrik Charpentier and Matthews, John Alphonsus and
             Dahl, Johannes Rom and Seim, Andrea},
  journal = {Nature Food},
  year    = {2026},
  doi     = {10.1038/s43016-026-01460-8}
}
```

Questions: open an issue or contact Emile Esmaili (ee4561@princeton.edu).
