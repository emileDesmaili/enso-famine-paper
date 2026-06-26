"""
08_spatial_maps.py  –  Supplementary teleconnection maps for appendix.

Produces 2×2 panel maps (nino3, nino34, nino4, nino12) for each variable:
  • PDSI          (OWDA scPDSI, JJA season)
  • Summer Temp   (ModE-RA, AMJJ)
  • Winter Temp   (ModE-RA, NDJF)
  • Summer Precip (ModE-RA, AMJJ)
  • Winter Precip (ModE-RA, NDJF)

Outputs saved to analysis/output/figures/appendix/:
  PDSI_telecomap.pdf
  SummerTemp_telecomap.pdf
  WinterTemp_telecomap.pdf
  SummerPrecip_telecomap.pdf
  WinterPrecip_telecomap.pdf
"""

from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import cftime
from scipy.stats import pearsonr

# ── paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
ROOT       = SCRIPT_DIR.parent
DATA_RAW   = ROOT / "data"
OUT_DIR    = SCRIPT_DIR / "output" / "figures" / "appendix"
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_RANGE = slice(1500, 1800)
ALPHA      = 0.10
EXTENT     = [-10, 35, 36, 70]   # lon_min, lon_max, lat_min, lat_max
PROJ       = ccrs.EuroPP()

# ── rcParams ───────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       11,
    "axes.titlesize":  11,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})


# ══════════════════════════════════════════════════════════════════════════════
# ENSO index construction
# ══════════════════════════════════════════════════════════════════════════════
def _build_enso_df() -> pd.DataFrame:
    enso3   = pd.read_csv(DATA_RAW / "cook2024-R15-ENSO-Rec-1500-2000.txt",
                          delimiter="\t", comment="#", na_values="NA")
    lat_lon = pd.read_csv(DATA_RAW / "cook2024-ENSO-latlon.txt",
                          delimiter="\t", comment="#", na_values="NA")
    em      = enso3.melt(id_vars=["Year"], var_name="gridpoint", value_name="enso")
    em["gridpoint"] = em["gridpoint"].astype(int)
    em      = em.merge(lat_lon, on="gridpoint")

    enso_xr = em.set_index(["Year", "lat", "lon"])["enso"].to_xarray()
    enso_xr = enso_xr - enso_xr.sel(Year=slice(1801, 1900)).mean(dim="Year")

    nino3  = enso_xr.sel(lat=slice(-5, 5),  lon=slice(-150, -90)).mean(dim=["lat", "lon"])
    nino34 = enso_xr.sel(lat=slice(-5, 5),  lon=slice(-170, -120)).mean(dim=["lat", "lon"])
    nino4  = enso_xr.sel(lat=slice(-5, 5),  lon=slice(-200, -150)).mean(dim=["lat", "lon"])
    nino12 = enso_xr.sel(lat=slice(-10, 0), lon=slice(-90, -80)).mean(dim=["lat", "lon"])

    enso = nino3.to_dataframe().reset_index().rename(columns={"enso": "nino3"})
    for key, da in [("nino34", nino34), ("nino4", nino4), ("nino12", nino12)]:
        enso = enso.merge(da.to_dataframe().reset_index().rename(columns={"enso": key}),
                          on="Year")
    return enso[enso["Year"] <= 1800]


# ══════════════════════════════════════════════════════════════════════════════
# Correlation helpers
# ══════════════════════════════════════════════════════════════════════════════
def _corr_grid_xr(field_da: xr.DataArray, nino_da: xr.DataArray):
    """Compute Pearson r and p-value at each grid point."""
    lat  = field_da["latitude"].values
    lon  = field_da["longitude"].values
    lon2d, lat2d = np.meshgrid(lon, lat)

    nv    = nino_da.values
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


# ══════════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════════
def _draw_panel(ax, corr_vals, p_vals, lat2d, lon2d,
                nino_name: str, cmap: str, levels: np.ndarray):
    cmap_obj = plt.get_cmap(cmap, len(levels) - 1)
    norm     = BoundaryNorm(levels, ncolors=cmap_obj.N)

    im = ax.pcolormesh(lon2d, lat2d, corr_vals,
                       cmap=cmap_obj, norm=norm, shading="auto",
                       transform=ccrs.PlateCarree())
    if p_vals is not None:
        sig = p_vals <= ALPHA
        ax.scatter(lon2d[sig], lat2d[sig], s=4, color="black",
                   transform=ccrs.PlateCarree(), zorder=5)

    ax.coastlines(linewidth=0.7)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.5)
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.set_title(f"Correlation with {nino_name.upper()}", fontsize=10)

    cbar = plt.colorbar(im, ax=ax, orientation="vertical",
                        fraction=0.024, pad=0.04)
    cbar.set_label("Pearson r", fontsize=9)
    cbar.set_ticks(levels[::2])


def _make_4panel(results: dict, title: str, cmap: str,
                 levels: np.ndarray, out_name: str):
    """results: {nino_name: (corr_vals, p_vals, lat2d, lon2d)}"""
    nino_keys = ["nino3", "nino34", "nino4", "nino12"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9),
                             subplot_kw={"projection": PROJ})
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)

    for ax, key in zip(axes.flat, nino_keys):
        corr, pval, lat2d, lon2d = results[key]
        _draw_panel(ax, corr, pval, lat2d, lon2d, key, cmap, levels)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = OUT_DIR / out_name
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path.name}")


# ══════════════════════════════════════════════════════════════════════════════
# Dataset loaders
# ══════════════════════════════════════════════════════════════════════════════
def _load_owda(nino_xr: xr.Dataset) -> xr.Dataset:
    owda = xr.open_dataset(DATA_RAW / "owda.nc")
    owda["PDSI_w"] = (np.sqrt(np.cos(np.deg2rad(owda["lat"])) + 1e-6)
                      * owda["pdsi"])
    owda = (owda.assign_coords(Year=owda.time)
            .swap_dims({"time": "Year"})
            .sel(Year=YEAR_RANGE)
            .sel(lat=slice(35, 70), lon=slice(-15, 35))
            .rename({"lat": "latitude", "lon": "longitude"}))
    return xr.merge([owda, nino_xr])


def _load_temp(nino_xr: xr.Dataset):
    """Returns (winter_merged, summer_merged) for temperature."""
    xds = xr.open_dataset(
        DATA_RAW / "ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc",
        use_cftime=True)
    xds = xds.sel(time=slice("1500-01-01", "1800-12-31"))
    xds["temp"] = (np.sqrt(np.cos(np.deg2rad(xds["latitude"])) + 1e-6)
                   * xds["temp2"])
    eu = xds.sel(latitude=slice(70, 35), longitude=slice(-15, 35))

    def shift_dec(t_arr):
        return [cftime.DatetimeGregorian(t.year + 1, t.month, t.day)
                if t.month in [11, 12] else t for t in t_arr]

    eu_sh = eu.assign_coords(time=("time", shift_dec(eu.time.values)))
    eu_w_sel = eu_sh.sel(time=eu_sh.time.dt.month.isin([11, 12, 1, 2]))
    eu_w = eu_w_sel.groupby(eu_w_sel.time.dt.year).mean(dim="time")
    eu_w = eu_w.rename({[d for d in eu_w.dims if d not in ("latitude", "longitude")][0]: "Year"})
    eu_s_sel = eu_sh.sel(time=eu_sh.time.dt.month.isin([5, 6, 7, 8]))
    eu_s = eu_s_sel.groupby(eu_s_sel.time.dt.year).mean(dim="time")
    eu_s = eu_s.rename({[d for d in eu_s.dims if d not in ("latitude", "longitude")][0]: "Year"})
    return xr.merge([eu_w, nino_xr]), xr.merge([eu_s, nino_xr])


def _load_precip(nino_xr: xr.Dataset):
    """Returns (winter_merged, summer_merged) for precipitation."""
    xds = xr.open_dataset(
        DATA_RAW / "ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc",
        use_cftime=True)
    xds = xds.sel(time=slice("1500-01-01", "1800-12-31"))
    xds["precip"] = (np.sqrt(np.cos(np.deg2rad(xds["latitude"])) + 1e-6)
                     * xds["totprec"] * 86400 * 30)
    eu = xds.sel(latitude=slice(70, 35), longitude=slice(-15, 35))

    def shift_dec(t_arr):
        return [cftime.DatetimeGregorian(t.year + 1, t.month, t.day)
                if t.month in [11, 12] else t for t in t_arr]

    eu_sh = eu.assign_coords(time=("time", shift_dec(eu.time.values)))
    eu_w_sel = eu_sh.sel(time=eu_sh.time.dt.month.isin([11, 12, 1, 2]))
    eu_w = eu_w_sel.groupby(eu_w_sel.time.dt.year).mean(dim="time")
    eu_w = eu_w.rename({[d for d in eu_w.dims if d not in ("latitude", "longitude")][0]: "Year"})
    eu_s_sel = eu_sh.sel(time=eu_sh.time.dt.month.isin([4, 5, 6, 7]))
    eu_s = eu_s_sel.groupby(eu_s_sel.time.dt.year).mean(dim="time")
    eu_s = eu_s.rename({[d for d in eu_s.dims if d not in ("latitude", "longitude")][0]: "Year"})
    return xr.merge([eu_w, nino_xr]), xr.merge([eu_s, nino_xr])


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
def _compute_results(merged: xr.Dataset, var_name: str,
                     nino_keys: list[str]) -> dict:
    da = merged[var_name].transpose("Year", "latitude", "longitude")
    results = {}
    for key in nino_keys:
        nino_s = merged[key]
        common = np.intersect1d(da["Year"].values, nino_s["Year"].values)
        c, p, la, lo = _corr_grid_xr(da.sel(Year=common), nino_s.sel(Year=common))
        results[key] = (c, p, la, lo)
    return results


def main():
    nino_keys = ["nino3", "nino34", "nino4", "nino12"]
    levels_bwg = np.linspace(-0.2, 0.2, 9)
    levels_rbu = np.linspace(-0.2, 0.2, 9)

    print("Building ENSO indices …")
    enso_df = _build_enso_df()
    nino_xr = xr.Dataset.from_dataframe(enso_df.set_index("Year"))

    # ── PDSI ──────────────────────────────────────────────────────────────────
    owda_path = DATA_RAW / "owda.nc"
    if owda_path.exists():
        print("PDSI teleconnections …")
        merged = _load_owda(nino_xr)
        results = _compute_results(merged, "PDSI_w", nino_keys)
        _make_4panel(results,
                     title="NINO × JJA scPDSI (1500–1800)",
                     cmap="BrBG", levels=levels_bwg,
                     out_name="PDSI_telecomap.pdf")
    else:
        print("owda.nc not found – skipping PDSI")

    # ── Temperature ───────────────────────────────────────────────────────────
    temp_path = DATA_RAW / "ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc"
    if temp_path.exists():
        print("Temperature teleconnections …")
        merged_tw, merged_ts = _load_temp(nino_xr)
        for season_label, merged, out_name, title in [
            ("winter", merged_tw, "WinterTemp_telecomap.pdf",
             "NINO × NDJF Temperature (1500–1800)"),
            ("summer", merged_ts, "SummerTemp_telecomap.pdf",
             "NINO × AMJJ Temperature (1500–1800)"),
        ]:
            print(f"  {season_label} …")
            results = _compute_results(merged, "temp", nino_keys)
            _make_4panel(results, title=title, cmap="RdBu_r",
                         levels=levels_rbu, out_name=out_name)
    else:
        print("ModE-RA temp file not found – skipping temperature")

    # ── Precipitation ─────────────────────────────────────────────────────────
    prec_path = DATA_RAW / "ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc"
    if prec_path.exists():
        print("Precipitation teleconnections …")
        merged_pw, merged_ps = _load_precip(nino_xr)
        for season_label, merged, out_name, title in [
            ("winter", merged_pw, "WinterPrecip_telecomap.pdf",
             "NINO × NDJF Precipitation (1500–1800)"),
            ("summer", merged_ps, "SummerPrecip_telecomap.pdf",
             "NINO × AMJJ Precipitation (1500–1800)"),
        ]:
            print(f"  {season_label} …")
            results = _compute_results(merged, "precip", nino_keys)
            _make_4panel(results, title=title, cmap="BrBG",
                         levels=levels_bwg, out_name=out_name)
    else:
        print("ModE-RA precip file not found – skipping precipitation")

    print(f"\nAll supplementary teleconnection maps saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
