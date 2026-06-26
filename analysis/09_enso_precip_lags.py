"""
09_enso_precip_lags.py  -  Extended-data figure: ENSO x (AMJJ precip / JJA scPDSI)
correlation maps at concurrent (T) and one-year lag (T+1).

Year T = ENSO year (NINO3.4 reconstruction, Cook 2024).
Rows:
  Top    - AMJJ precipitation (ModE-RA)
  Bottom - JJA scPDSI (OWDA)

Columns:
  (T)   concurrent
  (T+1) field one year after ENSO

Style matches the main-figure teleconnection maps in assemble_figures.py
(BrBG colormap, BoundaryNorm, EuroPP projection, stippling for p<=0.10).

Saved to analysis/output/figures/extended data/figED_enso_precip_lags.pdf
"""

from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import cftime
from scipy.stats import pearsonr

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# -- paths --------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).parent
ROOT       = SCRIPT_DIR.parent
DATA_RAW   = ROOT / "data"
DATA_PROC  = ROOT / "processed data"
OUT_DIR    = SCRIPT_DIR / "output" / "figures" / "extended data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_RANGE = slice(1500, 1800)
ALPHA      = 0.10
EXTENT     = [-10, 35, 36, 70]
PROJ       = ccrs.EuroPP()

# -- main-figure rcParams -----------------------------------------------------
FS  = 22
LAB = 24
PAN = 42

mpl.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":         FS,
    "axes.titlesize":    FS,
    "axes.labelsize":    LAB,
    "xtick.labelsize":   FS,
    "ytick.labelsize":   FS,
    "legend.fontsize":   FS,
    "axes.linewidth":    1.4,
    "xtick.major.width": 1.4,
    "ytick.major.width": 1.4,
    "xtick.major.size":  5,
    "ytick.major.size":  5,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
})

_PAN_FP = fm.FontProperties(size=PAN, weight="bold")


def _label(ax, letter, x=-0.10, y=1.05):
    ax.text(x, y, letter, transform=ax.transAxes,
            font_properties=_PAN_FP, va="bottom", ha="left")


# -- ENSO index ---------------------------------------------------------------
def _build_nino34() -> pd.DataFrame:
    enso3   = pd.read_csv(DATA_RAW / "cook2024-R15-ENSO-Rec-1500-2000.txt",
                          delimiter="\t", comment="#", na_values="NA")
    lat_lon = pd.read_csv(DATA_RAW / "cook2024-ENSO-latlon.txt",
                          delimiter="\t", comment="#", na_values="NA")
    em = enso3.melt(id_vars=["Year"], var_name="gridpoint", value_name="enso")
    em["gridpoint"] = em["gridpoint"].astype(int)
    em = em.merge(lat_lon, on="gridpoint")

    enso_xr = em.set_index(["Year", "lat", "lon"])["enso"].to_xarray()
    enso_xr = enso_xr - enso_xr.sel(Year=slice(1801, 1900)).mean(dim="Year")
    nino34  = enso_xr.sel(lat=slice(-5, 5),
                          lon=slice(-170, -120)).mean(dim=["lat", "lon"])
    df = nino34.to_dataframe().reset_index().rename(columns={"enso": "nino34"})
    return df[df["Year"] <= 1800]


# -- AMJJ precip --------------------------------------------------------------
def _load_amjj_precip() -> xr.DataArray:
    """Yearly AMJJ-mean ModE-RA precipitation over Europe (Year, lat, lon)."""
    xds = xr.open_dataset(
        DATA_RAW / "ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc",
        use_cftime=True,
    )
    xds = xds.sel(time=slice("1499-01-01", "1801-12-31"))
    xds["precip"] = (np.sqrt(np.cos(np.deg2rad(xds["latitude"])) + 1e-6)
                     * xds["totprec"] * 86400 * 30)
    eu = xds.sel(latitude=slice(70, 35), longitude=slice(-15, 35))
    eu_amjj = eu.sel(time=eu.time.dt.month.isin([4, 5, 6, 7]))
    annual = eu_amjj.groupby(eu_amjj.time.dt.year).mean(dim="time")
    year_dim = [d for d in annual.dims if d not in ("latitude", "longitude")][0]
    annual = annual.rename({year_dim: "Year"})
    return annual["precip"].transpose("Year", "latitude", "longitude")


# -- AMJJ temp ---------------------------------------------------------------
def _load_amjj_temp() -> xr.DataArray:
    """Yearly AMJJ-mean ModE-RA 2-m temperature anomaly over Europe."""
    xds = xr.open_dataset(
        DATA_RAW / "ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc",
        use_cftime=True,
    )
    xds = xds.sel(time=slice("1499-01-01", "1801-12-31"))
    xds["temp"] = (np.sqrt(np.cos(np.deg2rad(xds["latitude"])) + 1e-6)
                   * xds["temp2"])
    eu = xds.sel(latitude=slice(70, 35), longitude=slice(-15, 35))
    eu_amjj = eu.sel(time=eu.time.dt.month.isin([4, 5, 6, 7]))
    annual = eu_amjj.groupby(eu_amjj.time.dt.year).mean(dim="time")
    year_dim = [d for d in annual.dims if d not in ("latitude", "longitude")][0]
    annual = annual.rename({year_dim: "Year"})
    return annual["temp"].transpose("Year", "latitude", "longitude")


# -- OWDA scPDSI --------------------------------------------------------------
def _load_owda_pdsi() -> xr.DataArray:
    """Annual JJA scPDSI (OWDA) over Europe (Year, lat, lon)."""
    owda = xr.open_dataset(DATA_RAW / "owda.nc")
    owda["PDSI_w"] = (np.sqrt(np.cos(np.deg2rad(owda["lat"])) + 1e-6)
                      * owda["pdsi"])
    owda = (owda.assign_coords(Year=owda.time)
            .swap_dims({"time": "Year"})
            .sel(Year=slice(1499, 1801))
            .sel(lat=slice(35, 70), lon=slice(-15, 35))
            .rename({"lat": "latitude", "lon": "longitude"}))
    return owda["PDSI_w"].transpose("Year", "latitude", "longitude")


# -- correlation grid ---------------------------------------------------------
def _corr_grid(field_da: xr.DataArray, nino_series: xr.DataArray):
    lat = field_da["latitude"].values
    lon = field_da["longitude"].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    nv    = nino_series.values
    corrs = np.full((len(lat), len(lon)), np.nan)
    pvals = np.full_like(corrs, np.nan)
    for i in range(len(lat)):
        for j in range(len(lon)):
            fts  = field_da[:, i, j].values
            mask = np.isfinite(fts) & np.isfinite(nv)
            if mask.sum() > 10:
                r, p = pearsonr(fts[mask], nv[mask])
                corrs[i, j] = r
                pvals[i, j] = p
    return corrs, pvals, lat2d, lon2d


def _compute_lag(field_da: xr.DataArray, nino_df: pd.DataFrame, lag: int):
    """
    Correlate field at year (T + lag) with ENSO at year T.

    lag = -1: field leads ENSO by one year (field at T-1)
    lag =  0: concurrent
    lag = +1: field lags ENSO by one year (field at T+1)
    """
    nino = nino_df.copy()
    nino["MatchYear"] = nino["Year"] + lag   # year of field to fetch

    field_years  = field_da["Year"].values.astype(int)
    common_match = np.intersect1d(field_years, nino["MatchYear"].values)
    nino = nino[nino["MatchYear"].isin(common_match)]

    field = field_da.sel(Year=nino["MatchYear"].values)
    field = field.assign_coords(Year=("Year", nino["Year"].values))
    series_xr = xr.DataArray(nino["nino34"].values,
                             coords={"Year": nino["Year"].values},
                             dims=["Year"])
    series_xr = series_xr.sel(Year=slice(1500, 1800))
    field     = field.sel(Year=slice(1500, 1800))
    common = np.intersect1d(field["Year"].values, series_xr["Year"].values)
    return _corr_grid(field.sel(Year=common), series_xr.sel(Year=common))


# -- drawing ------------------------------------------------------------------
def _draw_panel(ax, corr_vals, p_vals, lat2d, lon2d, title, cmap_name="BrBG"):
    levels = np.linspace(-0.2, 0.2, 17)
    cmap_  = plt.get_cmap(cmap_name, len(levels) + 1)
    norm   = BoundaryNorm(levels, ncolors=cmap_.N, extend="both")

    im = ax.pcolormesh(lon2d, lat2d, corr_vals,
                       cmap=cmap_, norm=norm, shading="auto",
                       transform=ccrs.PlateCarree())
    sig = p_vals <= ALPHA
    ax.scatter(lon2d[sig], lat2d[sig], s=25, color="black",
               transform=ccrs.PlateCarree(), zorder=5, alpha=0.7)

    ax.coastlines(linewidth=0.9)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.5)
    ax.add_feature(cfeature.LAND,  facecolor="whitesmoke", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="aliceblue",  zorder=0)
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.6, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-10, 0, 10, 20, 30])
    gl.ylocator = mticker.FixedLocator([40, 50, 60, 70])
    gl.xlabel_style = {"size": FS - 5}
    gl.ylabel_style = {"size": FS - 5}

    cb = plt.colorbar(im, ax=ax, orientation="vertical",
                      fraction=0.046, pad=0.06, extend="both")
    cb.ax.tick_params(labelsize=FS - 4)
    cb.set_label("Pearson r", fontsize=FS - 2)
    ax.set_title(title, fontsize=FS, fontweight="bold", loc="left")


# -- main ---------------------------------------------------------------------
def main():
    print("Building NINO3.4 index ...")
    nino_df = _build_nino34()

    print("Loading ModE-RA AMJJ precipitation ...")
    precip_da = _load_amjj_precip()
    print("Loading OWDA JJA scPDSI ...")
    pdsi_da = _load_owda_pdsi()

    lags = [0, 1]
    rows = [
        ("AMJJ precip", precip_da,
         [r"$\mathrm{corr}\,$(ENSO$_T$, AMJJ precip$_T$)",
          r"$\mathrm{corr}\,$(ENSO$_T$, AMJJ precip$_{T+1}$)"],
         ["a", "b"]),
        ("JJA scPDSI", pdsi_da,
         [r"$\mathrm{corr}\,$(ENSO$_T$, JJA scPDSI$_T$)",
          r"$\mathrm{corr}\,$(ENSO$_T$, JJA scPDSI$_{T+1}$)"],
         ["c", "d"]),
    ]

    fig = plt.figure(figsize=(20, 18))
    gs  = gridspec.GridSpec(2, 2, figure=fig, wspace=0.20, hspace=0.20)

    for row_i, (row_label, field, titles, letters) in enumerate(rows):
        for k, (lag, title, letter) in enumerate(zip(lags, titles, letters)):
            print(f"  {row_label} lag {lag:+d} ...")
            corr, pval, lat2d, lon2d = _compute_lag(field, nino_df, lag)
            ax = fig.add_subplot(gs[row_i, k], projection=PROJ)
            _draw_panel(ax, corr, pval, lat2d, lon2d, title)
            _label(ax, letter)

    out_path = OUT_DIR / "figED_enso_precip_lags.pdf"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
