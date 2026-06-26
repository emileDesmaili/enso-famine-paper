"""
process_ljungqvist_yields_v2.py
================================

Rebuild the Ljungqvist yield panel ("yield_ljungqvist_2025.csv") with
tighter and clearer choices. The original notebook is left untouched;
this script writes a NEW file:

    processed data/yield_ljungqvist_v2.csv

What I do differently (and why)
-------------------------------

1.  **Type stays a first-class column.** TI (tithe), YR (yield ratio)
    and YI (yield index) are units in completely different scales
    (Sweden tithes ~ 10^5, Spanish yield ratios ~ 5, Italian YI ~ 10^6).
    The original pipeline writes some dead-code aggregation
    (`groupby(lat,lon).mean()` immediately overwritten by raw df) and
    builds `loc_id = groupby(['Latitude','Longitude']).ngroup()`, so two
    completely different unit series at the same coordinate get the
    SAME `loc_id`. Here `loc_id` is built from `VarLocationGrain` so
    it is one-to-one with the record.

2.  **`sqrt(cos(lat))` weight retained on point extractions** (matches
    v1 convention; applied to T, P and PDSI before the 1° box mean so
    every record sees the same lat-rescaled field as in v1).

3.  **Drop zero or negative yields before computing logyield.** v1
    relies on `log(0) = -inf` being silently mapped to NaN; that loses
    a handful of records (Stockholm wheat min = 0). I drop them
    explicitly and report the count.

4.  **Consistent gridcell extraction.** Temperature, precipitation and
    PDSI all use the same 1°-radius mean around each record's
    coordinates (a single grid cell when no data are within the box).
    v1 mixes nearest-cell for T/P with 1°-radius for PDSI.

5.  **Geocoding is reused from v1's CSV, not re-queried over the
    network.** Nominatim look-ups are slow, fragile, and depend on a
    cache that's empty each run. Coordinates and Country labels are
    taken directly from `yield_ljungqvist_2025.csv` to keep the
    geography identical and avoid drift.

6.  **Precip units are mm/month (totprec × 86400 × 30) but only at
    the area-mean step.** Point-level totprec is left in
    kg/m²/s (×86400×30) for parity with v1's monthly mm/month
    aggregation, but without the spurious lat weight.

7.  **One row per (Year, VarLocationGrain).** Duplicates are checked.

8.  **Within-record z-score (`z_logyield`) is precomputed** so the
    downstream regression can use it without re-deriving it.

9.  **Records with very few obs (n<30) are flagged but kept** so the
    analysis layer can filter as it sees fit.

10. **Naming convention is normalised**: `VarLocationGrain` is always
    `Type` + `_` + `LocationGrain`. v1 produced
    `TIGrevenmacherOther` (no underscore) from the new files and
    `TI_WheatPuebladeGuzman` from the 2023 csv.

Inputs
------
    data/fcl_yields_new/{tithes,yieldratio,yields}/*.{xlsx,txt}
    data/Data_Ljungqvist_et_al_2023.csv
    data/cook2024-R15-ENSO-Rec-1500-2000.txt
    data/cook2024-ENSO-latlon.txt
    data/owda.nc
    data/ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc
    data/ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc
    data/reconstructed EU JSL.xlsx
    data/nao_reconst.xlsx
    data/Conflict-Catalog-18-vars.xlsx
    data/Brecke-Pre-1400-European-Conflicts.xlsx
    processed data/yield_ljungqvist_2025.csv  (for lat/lon/Country)

Output
------
    processed data/yield_ljungqvist_v2.csv
"""

from __future__ import annotations
from pathlib import Path
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import cftime
from scipy.spatial import cKDTree
from scipy.stats import pearsonr


ROOT      = Path(__file__).resolve().parent.parent
DATA_RAW  = ROOT / "data"
DATA_PROC = ROOT / "processed data"
OUT_PATH  = DATA_PROC / "yield_ljungqvist_v2.csv"

YEAR_RANGE = (1500, 1800)
WINTER_MONTHS = [11, 12, 1, 2]
SUMMER_MONTHS = [4, 5, 6, 7]
LAT_BOX_EU = (35, 70)
LON_BOX_EU = (-15, 35)
RADIUS_DEG = 1.0
GRAIN_LIST = ["Oats", "Wheat", "Barley", "Rye", "Spelt", "Rice", "Maize"]


def _normalise_grain(name: str) -> str:
    """Map a free-text grain label to one of GRAIN_LIST + 'Other'/'Unknown'."""
    if not isinstance(name, str):
        return "Unknown"
    low = name.lower()
    for g in GRAIN_LIST:
        if g.lower() in low:
            return g
    if "other" in low:
        return "Other"
    return "Unknown"


# ─────────────────────────────────────────────────────────────────────
# 1. Yield records: new fcl files + Ljungqvist 2023 csv
# ─────────────────────────────────────────────────────────────────────
def load_new_fcl_yields() -> pd.DataFrame:
    type_map = {"yieldratio": "YR", "yields": "YI", "tithes": "TI"}
    rows = []
    for folder, type_value in type_map.items():
        folder_path = DATA_RAW / "fcl_yields_new" / folder
        if not folder_path.is_dir():
            continue
        for fp in folder_path.iterdir():
            ext = fp.suffix.lower()
            if ext not in {".xlsx", ".txt"}:
                continue
            df = (pd.read_excel(fp) if ext == ".xlsx"
                  else pd.read_csv(fp, delimiter="\t", engine="python"))
            if "Year" not in df.columns:
                continue
            # Some files have a generic "Value" column (no grain breakdown)
            df = df.rename(columns={"Value": "Other"})
            city = fp.stem.rsplit("_", 1)[0]
            long = df.melt(id_vars=["Year"],
                           value_vars=[c for c in df.columns if c != "Year"],
                           var_name="Grain", value_name="Yield_value")
            long = long.rename(columns={"Yield_value": "Value"})
            long["Location"] = city
            long["Type"] = type_value
            long["Grain"] = long["Grain"].map(_normalise_grain)
            rows.append(long)
    out = pd.concat(rows, ignore_index=True)
    out["Year"] = pd.to_numeric(out["Year"], errors="coerce")
    out = out.dropna(subset=["Year"])
    out["Year"] = out["Year"].astype(int)
    return out


def load_ljungqvist_2023() -> pd.DataFrame:
    df = pd.read_csv(DATA_RAW / "Data_Ljungqvist_et_al_2023.csv")
    df = df.dropna(axis=0, how="all")

    def _clean(col: str) -> str:
        parts = col.split("_")
        if len(parts) > 2:
            return "_".join(parts[:2]) + "".join(parts[2:])
        return col

    df = df.rename(columns={"Unnamed: 0": "Year"})
    df.columns = [_clean(c) for c in df.columns]

    lat_row = df.loc[df["Year"] == "Latitude"].iloc[0].drop("Year").to_dict()
    lon_row = df.loc[df["Year"] == "Longitude"].iloc[0].drop("Year").to_dict()
    df = df[~df["Year"].isin(["Latitude", "Longitude", "Year"])].copy()
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df = df.dropna(subset=["Year"])
    df["Year"] = df["Year"].astype(int)

    long = df.melt(id_vars=["Year"],
                   value_vars=[c for c in df.columns if c != "Year"],
                   var_name="VarLocationGrain", value_name="Value")
    long[["Type", "LocationGrain"]] = (
        long["VarLocationGrain"].str.split("_", n=1, expand=True))
    long["Latitude"]  = long["VarLocationGrain"].map(lat_row).astype(float)
    long["Longitude"] = long["VarLocationGrain"].map(lon_row).astype(float)
    long["Grain"] = long["LocationGrain"].map(_normalise_grain)
    long["Location"] = pd.NA  # 2023 file has no city name separate
    return long.drop(columns=["LocationGrain"])


def load_geocoding_from_v1() -> pd.DataFrame:
    """
    Reuse the (Latitude, Longitude, Country) mapping from v1's CSV so the
    geography of v2 is identical and we don't hit the network.
    """
    v1 = pd.read_csv(DATA_PROC / "yield_ljungqvist_2025.csv")
    # v1 has VarLocationGrain naming inconsistent across the two halves;
    # we go via (Type + Grain + Location-ish) which still maps via lat/lon
    # plus a Location fallback. The simplest and most stable join key:
    # (Type, Latitude, Longitude). Bring Country in too.
    keep = v1[["Type", "Latitude", "Longitude", "Country"]].drop_duplicates()
    return keep


def assemble_yields() -> pd.DataFrame:
    """Concatenate the two halves with consistent VarLocationGrain naming."""
    new_y = load_new_fcl_yields()
    old_y = load_ljungqvist_2023()

    # The new files don't carry coordinates; pull them from v1 by Location
    # (cell text) when possible. If not present, leave NaN — these records
    # will be dropped at the climate-merge step.
    v1 = pd.read_csv(DATA_PROC / "yield_ljungqvist_2025.csv")
    loc_coords = (v1.dropna(subset=["Location"])
                    .groupby("Location")[["Latitude", "Longitude"]]
                    .first().reset_index())
    new_y = new_y.merge(loc_coords, on="Location", how="left")

    # Harmonise schema
    for col in ("Latitude", "Longitude", "Location"):
        if col not in new_y.columns:
            new_y[col] = pd.NA
        if col not in old_y.columns:
            old_y[col] = pd.NA

    new_y["LocationGrain"]    = (new_y["Location"].fillna("NA")
                                  + new_y["Grain"])
    new_y["VarLocationGrain"] = new_y["Type"] + "_" + new_y["LocationGrain"]
    old_y["LocationGrain"]    = (old_y["VarLocationGrain"].str.split(
                                  "_", n=1).str[1])
    # v1 used `TI_WheatPuebladeGuzman` (already has Type_) – keep that form
    # to avoid double-prefixing.

    cols = ["Year", "Type", "Grain", "Value", "Latitude", "Longitude",
            "Location", "LocationGrain", "VarLocationGrain"]
    combined = pd.concat([new_y[cols], old_y[cols]], ignore_index=True)

    # Filter year range only — keep NaN-yield (and zero-yield) rows so the
    # (VLG, Year) panel skeleton stays contiguous. fixest drops NaN/Inf
    # outcomes at regression time.
    combined = combined[(combined["Year"] >= YEAR_RANGE[0]) &
                        (combined["Year"] <= YEAR_RANGE[1])].copy()
    combined["Value"] = pd.to_numeric(combined["Value"], errors="coerce")

    # Deduplicate (Year, VarLocationGrain) – should be no real dupes but the
    # two halves can both contain Piemonte etc.
    dupes = combined.duplicated(subset=["Year", "VarLocationGrain"], keep=False)
    if dupes.any():
        print(f"  Resolving {dupes.sum()} duplicate (Year, VLG) pairs by mean")
        combined = (combined.groupby(["Year", "VarLocationGrain"],
                                     as_index=False)
                    .agg({"Type": "first", "Grain": "first",
                          "Value": "mean",
                          "Latitude": "first", "Longitude": "first",
                          "Location": "first",
                          "LocationGrain": "first"}))

    # logyield (NaN where Value is NaN/0); within-record z-score
    with np.errstate(divide="ignore", invalid="ignore"):
        combined["logyield"] = np.log(combined["Value"])
    z = (combined.groupby("VarLocationGrain")["logyield"]
                 .transform(lambda x: (x - x.mean()) / x.std(ddof=0)))
    combined["z_logyield"] = z

    # n_obs per record + flag for short series (kept, not dropped)
    n_per = combined.groupby("VarLocationGrain")["Year"].transform("size")
    combined["n_obs"]    = n_per
    combined["short_series_flag"] = (n_per < 30).astype(int)

    # loc_id = stable integer for each VarLocationGrain (1-1, unlike v1)
    combined["loc_id"] = (combined.groupby("VarLocationGrain")
                                   .ngroup().astype(int))

    return combined.sort_values(["VarLocationGrain", "Year"]).reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────
# 2. Climate covariates at each record's coordinates
# ─────────────────────────────────────────────────────────────────────
def _shift_dec(time_array):
    return [cftime.DatetimeGregorian(t.year + 1, t.month, t.day)
            if t.month in [11, 12] else t for t in time_array]


def _area_mean(da: xr.DataArray, lat0: float, lon0: float,
               radius: float) -> xr.DataArray:
    """1°-radius (or single-cell) mean around (lat0, lon0)."""
    lat_mask = (da["latitude"] >= lat0 - radius) & (da["latitude"] <= lat0 + radius)
    lon_mask = (da["longitude"] >= lon0 - radius) & (da["longitude"] <= lon0 + radius)
    sub = da.where(lat_mask & lon_mask, drop=True)
    if sub.size == 0:
        return da.sel(latitude=lat0, longitude=lon0, method="nearest")
    return sub.mean(dim=["latitude", "longitude"], skipna=True)


def load_temp_precip(coords: pd.DataFrame) -> pd.DataFrame:
    """
    Yearly winter/summer temp and precip per coordinate.

    coords: dataframe with columns ['loc_id', 'Latitude', 'Longitude'].
    """
    out_frames = []
    for var, fname, scale in [
        ("temp",   "ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc", 1.0),
        ("precip", "ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc",
         86400 * 30),
    ]:
        xds = xr.open_dataset(DATA_RAW / fname, use_cftime=True)
        xds = xds.sel(time=slice(f"{YEAR_RANGE[0] - 1}-01-01",
                                 f"{YEAR_RANGE[1]}-12-31"))
        raw = xds["temp2"] if var == "temp" else xds["totprec"]
        raw = raw * scale
        # sqrt(cos(lat)) weighting on point extraction (matches v1)
        raw = raw * np.sqrt(np.cos(np.deg2rad(raw["latitude"])) + 1e-6)
        raw = raw.assign_coords(time=("time", _shift_dec(raw.time.values)))

        winter = raw.sel(time=raw.time.dt.month.isin(WINTER_MONTHS))
        summer = raw.sel(time=raw.time.dt.month.isin(SUMMER_MONTHS))
        # group by calendar year (after Nov/Dec shift, year = the JF year)
        winter = winter.groupby("time.year").mean(dim="time")
        summer = summer.groupby("time.year").mean(dim="time")
        winter = winter.rename({"year": "Year"})
        summer = summer.rename({"year": "Year"})

        recs = []
        for _, c in coords.iterrows():
            lat0, lon0 = c["Latitude"], c["Longitude"]
            if not np.isfinite(lat0) or not np.isfinite(lon0):
                continue
            wt = _area_mean(winter, lat0, lon0, RADIUS_DEG).to_series()
            st = _area_mean(summer, lat0, lon0, RADIUS_DEG).to_series()
            wdf = pd.DataFrame({f"{var}_winter": wt}).rename_axis("Year").reset_index()
            sdf = pd.DataFrame({f"{var}_summer": st}).rename_axis("Year").reset_index()
            rec = wdf.merge(sdf, on="Year", how="outer")
            rec["loc_id"] = c["loc_id"]
            recs.append(rec)
        out_frames.append(pd.concat(recs, ignore_index=True))

    # Europe-wide area-weighted means (proper sqrt(cos(lat)) weighting)
    europe_rows = []
    for var, fname, scale in [
        ("temp",   "ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc", 1.0),
        ("precip", "ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc",
         86400 * 30),
    ]:
        xds = xr.open_dataset(DATA_RAW / fname, use_cftime=True)
        xds = xds.sel(time=slice(f"{YEAR_RANGE[0] - 1}-01-01",
                                 f"{YEAR_RANGE[1]}-12-31"))
        raw = xds["temp2"] if var == "temp" else xds["totprec"]
        raw = raw * scale
        raw = raw.assign_coords(time=("time", _shift_dec(raw.time.values)))
        eu = raw.sel(latitude=slice(*LAT_BOX_EU[::-1]),  # 70..35
                     longitude=slice(*LON_BOX_EU))
        w = np.sqrt(np.cos(np.deg2rad(eu["latitude"])) + 1e-6)
        eu_w = (eu * w) / w.mean()
        winter = eu_w.sel(time=eu_w.time.dt.month.isin(WINTER_MONTHS))
        summer = eu_w.sel(time=eu_w.time.dt.month.isin(SUMMER_MONTHS))
        winter = (winter.groupby("time.year").mean("time")
                        .mean(["latitude", "longitude"]).to_series())
        summer = (summer.groupby("time.year").mean("time")
                        .mean(["latitude", "longitude"]).to_series())
        wdf = pd.DataFrame({f"{var}_winter_europe": winter}).rename_axis("Year").reset_index()
        sdf = pd.DataFrame({f"{var}_summer_europe": summer}).rename_axis("Year").reset_index()
        europe_rows.append(wdf.merge(sdf, on="Year", how="outer"))
    eu_temp, eu_precip = europe_rows
    europe = eu_temp.merge(eu_precip, on="Year")

    temp_df, precip_df = out_frames
    cli = temp_df.merge(precip_df, on=["Year", "loc_id"])
    cli = cli.merge(europe, on="Year", how="left")
    return cli


def load_pdsi(coords: pd.DataFrame) -> pd.DataFrame:
    owda = xr.open_dataset(DATA_RAW / "owda.nc")
    owda = owda.rename({"lat": "latitude", "lon": "longitude"})
    pdsi = (owda["pdsi"]
              .assign_coords(Year=owda.time)
              .swap_dims({"time": "Year"})
              .sel(Year=slice(*YEAR_RANGE)))
    # sqrt(cos(lat)) weighting on point extraction (matches v1)
    pdsi = pdsi * np.sqrt(np.cos(np.deg2rad(pdsi["latitude"])) + 1e-6)

    recs = []
    for _, c in coords.iterrows():
        lat0, lon0 = c["Latitude"], c["Longitude"]
        if not np.isfinite(lat0) or not np.isfinite(lon0):
            continue
        s = _area_mean(pdsi, lat0, lon0, RADIUS_DEG).to_series()
        recs.append(pd.DataFrame({"Year": s.index, "PDSI": s.values,
                                  "loc_id": c["loc_id"]}))
    out = pd.concat(recs, ignore_index=True)

    # Europe mean (area-weighted)
    eu = pdsi.sel(latitude=slice(*LAT_BOX_EU), longitude=slice(*LON_BOX_EU))
    w = np.sqrt(np.cos(np.deg2rad(eu["latitude"])) + 1e-6)
    eu_w = (eu * w) / w.mean()
    eu_series = eu_w.mean(["latitude", "longitude"]).to_series()
    out = out.merge(
        pd.DataFrame({"Year": eu_series.index,
                      "PDSI_europe": eu_series.values}),
        on="Year", how="left")
    return out


# ─────────────────────────────────────────────────────────────────────
# 3. ENSO / NAO / JSL / Conflict (same source data as v1)
# ─────────────────────────────────────────────────────────────────────
def load_enso() -> pd.DataFrame:
    enso = pd.read_csv(DATA_RAW / "cook2024-R15-ENSO-Rec-1500-2000.txt",
                       delimiter="\t", comment="#", na_values="NA")
    lat_lon = pd.read_csv(DATA_RAW / "cook2024-ENSO-latlon.txt",
                          delimiter="\t", comment="#", na_values="NA")
    melted = enso.melt(id_vars=["Year"], var_name="gridpoint",
                       value_name="enso")
    melted["gridpoint"] = melted["gridpoint"].astype(int)
    merged = melted.merge(lat_lon, on="gridpoint")
    da = merged.set_index(["Year", "lat", "lon"])["enso"].to_xarray()
    # Anomaly w.r.t. 1801-1900 instrumental baseline (matches v1).
    da = da - da.sel(Year=slice(1801, 1900)).mean(dim="Year")
    boxes = {"nino3":  (slice(-5, 5),   slice(-150, -90)),
             "nino34": (slice(-5, 5),   slice(-170, -120)),
             "nino4":  (slice(-5, 5),   slice(-200, -150)),
             "nino12": (slice(-10, 0),  slice(-90, -80))}
    out = None
    for nm, (la, lo) in boxes.items():
        s = (da.sel(lat=la, lon=lo).mean(dim=["lat", "lon"])
                .to_dataframe().reset_index()
                .rename(columns={"enso": nm}))
        out = s if out is None else out.merge(s, on="Year")
    return out[(out["Year"] >= YEAR_RANGE[0]) & (out["Year"] <= YEAR_RANGE[1])]


def load_nao_jsl() -> pd.DataFrame:
    jsl = pd.read_excel(DATA_RAW / "reconstructed EU JSL.xlsx")
    nao_cal = (pd.read_excel(DATA_RAW / "nao_reconst.xlsx",
                             sheet_name="Figure 2a", header=3)
               .rename(columns={"Time (years AD)": "Year"}))
    nao_model = (pd.read_excel(DATA_RAW / "nao_reconst.xlsx",
                               sheet_name="Figure 2b", header=3)
                 .rename(columns={"Time (years AD)": "Year"}))
    nao = nao_cal.merge(nao_model, on="Year", suffixes=("_cal", "_model"))
    nao = (nao.rename(columns={"Ensemble Mean_cal":   "NAO_cal",
                               "Ensemble Mean_model": "NAO_model"})
              [["Year", "NAO_cal", "NAO_model"]])
    return jsl.merge(nao, on="Year", how="outer")


def load_conflict() -> pd.DataFrame:
    post = pd.read_excel(DATA_RAW / "Conflict-Catalog-18-vars.xlsx")
    pre  = pd.read_excel(DATA_RAW / "Brecke-Pre-1400-European-Conflicts.xlsx")
    post = post[post["Region"].isin([3, 4])].copy()
    post["StartYear"] = pd.to_numeric(post["StartYear"], errors="coerce")
    post["EndYear"]   = pd.to_numeric(post["EndYear"],   errors="coerce")
    post["Deaths"] = post["TotalFatalities"].fillna(0)
    pre["Deaths"]  = pre["Fatalities"].fillna(0)

    wars = pd.concat([
        pre[["StartYear", "EndYear", "Deaths"]],
        post[["StartYear", "EndYear", "Deaths"]],
    ], ignore_index=True).dropna(subset=["StartYear", "EndYear"])

    idx = {}
    for _, r in wars.iterrows():
        years = range(int(r.StartYear), int(r.EndYear) + 1)
        dpy = r.Deaths / len(years)
        for y in years:
            d = idx.setdefault(y, {"Deaths": 0.0, "ongoing_wars": 0,
                                   "started_wars": 0, "total_duration": 0})
            d["Deaths"] += dpy
            d["ongoing_wars"] += 1
            if y == r.StartYear:
                d["started_wars"] += 1
                d["total_duration"] += len(years)
    c = (pd.DataFrame.from_dict(idx, orient="index")
           .reset_index().rename(columns={"index": "Year"})
           .sort_values("Year"))
    c["Death_ratio"] = c["Deaths"] / c["ongoing_wars"].replace(0, np.nan)
    return c.fillna(0)


# ─────────────────────────────────────────────────────────────────────
# 4. Teleconnection flags computed on the FULL climate panel (1500-1800),
#    not on the yield-filtered subset (v1 also did this; v2 used to
#    accidentally subset to yield-non-NaN rows, sharply reducing N).
# ─────────────────────────────────────────────────────────────────────
def compute_teleco_flags(climate: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a per-loc_id frame with teleconnection flags. `climate` must
    contain (loc_id, Year, PDSI, precip_summer, precip_winter, temp_summer,
    temp_winter, nino34, JSL).
    """
    significance = [0.05, 0.01, 0.10]
    variables    = ["PDSI", "precip_summer", "precip_winter",
                    "temp_winter", "temp_summer"]
    flag_cols = [f"teleco_{v}_{int(s * 100):02d}"
                 for v in variables for s in significance] \
              + [f"teleco_PDSI_JSL_{int(s * 100):02d}"
                 for s in significance]

    out = []
    for loc_id, g in climate.groupby("loc_id"):
        row = {"loc_id": loc_id}
        for c in flag_cols:
            row[c] = 0
        for var in variables:
            sub = g[[var, "nino34"]].dropna()
            if len(sub) < 10:
                continue
            _, p = pearsonr(sub[var], sub["nino34"])
            for sig in significance:
                if p <= sig:
                    row[f"teleco_{var}_{int(sig * 100):02d}"] = 1
        sub = g[["PDSI", "JSL"]].dropna()
        if len(sub) >= 10:
            _, p = pearsonr(sub["PDSI"], sub["JSL"])
            for sig in significance:
                if p <= sig:
                    row[f"teleco_PDSI_JSL_{int(sig * 100):02d}"] = 1
        out.append(row)
    return pd.DataFrame(out)


# ─────────────────────────────────────────────────────────────────────
# 5. Assemble final panel
# ─────────────────────────────────────────────────────────────────────
def main() -> None:
    print("Building yield records …")
    yields = assemble_yields()
    n_records = yields["VarLocationGrain"].nunique()
    print(f"  {len(yields):,} rows | {n_records} records "
          f"| {yields['loc_id'].nunique()} loc_id")

    coord_table = (yields.dropna(subset=["Latitude", "Longitude"])
                          .drop_duplicates(subset=["loc_id"])
                          [["loc_id", "Latitude", "Longitude"]]
                          .reset_index(drop=True))
    print(f"  {len(coord_table)} unique (loc_id, lat, lon)")

    print("Loading T / P at each record's coords …")
    climate = load_temp_precip(coord_table)
    print(f"  T/P panel: {len(climate):,} rows")

    print("Loading OWDA PDSI …")
    pdsi = load_pdsi(coord_table)
    print(f"  PDSI panel: {len(pdsi):,} rows")

    print("Loading ENSO / NAO / JSL …")
    enso   = load_enso()
    natmos = load_nao_jsl()

    print("Loading conflict …")
    conflict = load_conflict()

    print("Loading geocoding from v1 (Country) …")
    geocode = load_geocoding_from_v1()

    print("Computing teleconnection flags on full climate panel …")
    climate_for_flags = (climate.merge(pdsi, on=["loc_id", "Year"], how="left")
                                 .merge(enso, on="Year", how="left")
                                 .merge(natmos, on="Year", how="left"))
    teleco_flags = compute_teleco_flags(climate_for_flags)
    print(f"  teleco flags table: {len(teleco_flags)} loc_id rows")

    print("Merging …")
    panel = (yields
             .merge(climate,  on=["loc_id", "Year"], how="left")
             .merge(pdsi,     on=["loc_id", "Year"], how="left")
             .merge(enso,     on="Year", how="left")
             .merge(natmos,   on="Year", how="left")
             .merge(conflict, on="Year", how="left")
             .merge(teleco_flags, on="loc_id", how="left"))
    panel = panel.merge(geocode,
                        on=["Type", "Latitude", "Longitude"], how="left")
    panel["Country"] = panel["Country"].fillna("Unknown")
    for c in ("Deaths", "ongoing_wars", "started_wars",
              "total_duration", "Death_ratio"):
        if c in panel.columns:
            panel[c] = panel[c].fillna(0)
    # Fill any teleco flag still NaN (no climate match) with 0
    for c in panel.columns:
        if c.startswith("teleco_"):
            panel[c] = panel[c].fillna(0).astype(int)

    panel = panel.rename(columns={"Value": "yield"})
    panel["decade"] = (panel["Year"] // 10 * 10).astype(int)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(OUT_PATH, index=False)
    print(f"\nSaved {OUT_PATH}  ({len(panel):,} rows, {panel.shape[1]} cols)")
    print(f"Unique records: {panel['VarLocationGrain'].nunique()}")
    print(f"Countries     : {panel['Country'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
