# Replication package: "El Niño amplified food insecurity in early modern Europe from 1500 to 1800"

Esmaili, E., Puma, M. J., Ludlow, F., Jobbová, E., Kumar, J., Holm, P., Ljungqvist, F. C., Matthews, J. A., Dahl, J. R. & Seim, A. *Nature Food* (2026, in press). https://doi.org/10.1038/s43016-026-01460-8

Data and code to reproduce all figures, tables and numbers in the main text, Extended Data and Supplementary Information.

## Quick start

```bash
conda env create -f environment.yml && conda activate famine-enso   # Python 3.10
Rscript install_R_packages.R                                          # R 4.5.1 packages

python "preprocessing codes/run_all.py"   # optional: rebuilds processed data/ (~1 min)
Rscript analysis/run_all.R                # models, tables, ED and SI figures (~5 min)
python analysis/run_all.py                # machine learning, maps, Figs. 1–5 (~4 min)
```

Outputs are written to `analysis/output/` (`figures/main`, `figures/extended data`, `figures/appendix`, `tables`, `data`). Run time is about 10 minutes on a laptop. All random procedures are seeded.

## Data

All data come from published sources and are included in `data/` (raw) and `processed data/` (analysis-ready), except three large climate fields. Download these into `data/` before running preprocessing or `analysis/run_all.py` (the R stage does not need them):

- `owda.nc`: Old World Drought Atlas (Cook et al. 2015), https://www1.ncdc.noaa.gov/pub/data/paleo/treering/reconstructions/europe/owda.nc
- `ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc` and `ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc`: ModE-RA (Valler et al. 2024), https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/ModE/ModE-RA/ensstat

Sources of the included data:

| Data | Source |
|---|---|
| Famine chronology | Ljungqvist, Seim & Collet (2024), *WIREs Clim. Change* 15, e859 |
| ENSO (NINO3.4) reconstruction | Cook & Cane (2024), *Geophys. Res. Lett.* 51, e2024GL109759 |
| NAO; jet stream latitude | Ortega et al. (2015), *Nature* 523, 71–74; Xu et al. (2024), *Nature* 634, 600–608 |
| Armed conflicts | Brecke (1999), Conflict Catalog |
| Volcanic forcing | Sigl et al. (2015), *Nature* 523, 543–549 |
| Grain prices | Ljungqvist & Seim (2024), *Eur. Econ. Rev.* 170, 104852; Federico, Schulze & Volckart (2021), *J. Econ. Hist.* 81, 276–308 |
| Grain harvests | Ljungqvist et al. (2023), *Clim. Past* 19, 2463–2491, and the sources in Supplementary Table 2 |
| Fish prices | Matthews (2021), PhD thesis, Trinity College Dublin; Allaire & Holm (2022), *Food Hist.* 20, 41–66; Allen & Unger (2019) |
| Fish catches | Holm et al. (2021), *Fish Fish.*; Norfish data, https://doi.org/10.6084/m9.figshare.c.5514351.v1 |
| Alternative wheat prices and yield ratios | Waldinger (2022), *J. Polit. Econ.* 130, 2275–2314 (processed Allen 2001 and Slicher van Bath 1963 data) |

Rights to the data remain with the original authors; please cite them when using the data.

## Code

| Script | Produces |
|---|---|
| `preprocessing codes/process_*.py` | Analysis datasets in `processed data/` |
| `analysis/fig2_onset.R` | Fig. 2a,b; Extended Data Table 1; SI permutation and leave-one-out tests (onset) |
| `analysis/fig3_harvest.R` | Fig. 3c,d; Extended Data Fig. 2; SI harvest robustness figures |
| `analysis/fig4_survival.R` | Fig. 4a,b; SI survival tables and figures |
| `analysis/fig5_prices.R` | Fig. 5; SI price and fish-catch figures |
| `analysis/ed_nao_jsl_enso.R`, `ed_volcanic.R` | Extended Data Figs. 1 and 4 |
| `analysis/si_nonlinearity.R`, `tables.R` | SI non-linearity figure; SI summary and fixed-effects tables |
| `analysis/ml_onset_survival.py` | Fig. 2c–e, Fig. 4c,d and SI machine-learning figures |
| `analysis/assemble_figures.py` | Figs. 1–5 (including the maps in Fig. 3a,b) |
| `analysis/mechanisms.py`, `spatial_maps.py`, `maps_data_coverage.py` | Extended Data Fig. 3; SI teleconnection and data-coverage maps and table |

## Licence and citation

Code: MIT licence (`LICENSE`). Please cite the article:

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

Contact: Emile Esmaili (ee4561@princeton.edu)
