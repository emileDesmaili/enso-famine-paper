"""Replace cells 5 and 6 in maps_data_coverage.ipynb with country-level table."""
import json, uuid, os

NB_PATH = os.path.join(os.path.dirname(__file__), "maps_data_coverage.ipynb")

with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)

def new_cell(source_str):
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": uuid.uuid4().hex[:8],
        "metadata": {},
        "outputs": [],
        "source": source_str,
    }

# ── Cell 5: build country-level rows ─────────────────────────────────────────
SRC_BUILD = """\
# ── Data coverage table – country level within each famine region ────────────
import os

famine = pd.read_csv(f"{ROOT}/processed data/famine_region_data.csv")

# Map price location -> (region, country)
PRICE_LOC = {
    # Great Britain
    "Cambridge":          ("Great Britain",   "England"),
    "Exeter":             ("Great Britain",   "England"),
    "London":             ("Great Britain",   "England"),
    "Oxford":             ("Great Britain",   "England"),
    # France
    "Angers":             ("France",          "France"),
    "Avignon":            ("France",          "France"),
    "Grenoble":           ("France",          "France"),
    "Limoges":            ("France",          "France"),
    "Paris":              ("France",          "France"),
    "Toulouse":           ("France",          "France"),
    "Tours":              ("France",          "France"),
    "Douai":              ("France",          "France"),
    "Strasbourg":         ("France",          "France"),
    # Low Countries
    "Amsterdam":          ("Low Countries",   "Netherlands"),
    "Arnhem":             ("Low Countries",   "Netherlands"),
    "Leiden":             ("Low Countries",   "Netherlands"),
    "Utrecht":            ("Low Countries",   "Netherlands"),
    "Brussels":           ("Low Countries",   "Belgium"),
    "Ghent":              ("Low Countries",   "Belgium"),
    "Xanten":             ("Low Countries",   "Germany"),
    # Central Europe
    "Aachen":             ("Central Europe",  "Germany"),
    "Augsburg":           ("Central Europe",  "Germany"),
    "Cologne":            ("Central Europe",  "Germany"),
    "D\\u00fcren":        ("Central Europe",  "Germany"),
    "Frankfurt":          ("Central Europe",  "Germany"),
    "Leipzig":            ("Central Europe",  "Germany"),
    "Munich":             ("Central Europe",  "Germany"),
    "Speyer":             ("Central Europe",  "Germany"),
    "W\\u00fcrzburg":     ("Central Europe",  "Germany"),
    "Basle":              ("Central Europe",  "Switzerland"),
    "Z\\u00fcrich":       ("Central Europe",  "Switzerland"),
    "Vienna":             ("Central Europe",  "Austria"),
    "Wels":               ("Central Europe",  "Austria"),
    # Nordic Countries
    "Gdansk":             ("Nordic Countries","Poland"),
    # Italy
    "Bassano del Grappa": ("Italy",           "Italy"),
    "Naples":             ("Italy",           "Italy"),
    "Pisa":               ("Italy",           "Italy"),
    "Siena":              ("Italy",           "Italy"),
    # Spain
    "Barcelona":          ("Spain",           "Spain"),
    "Madrid":             ("Spain",           "Spain"),
    "New Castile":        ("Spain",           "Spain"),
    "Valencia":           ("Spain",           "Spain"),
    "Coimbra ":           ("Spain",           "Portugal"),
    "Lisbon":             ("Spain",           "Portugal"),
    "\\u00c9vora":        ("Spain",           "Portugal"),
}

# Map yield country -> famine region
YIELD_REGION_MAP = {
    "Switzerland": "Central Europe",
    "Hungary":     "Central Europe",
    "France":      "France",
    "Spain":       "Spain",
    "Italy":       "Italy",
    "Sweden":      "Nordic Countries",
    "Denmark":     "Nordic Countries",
    # Germany sites are in Luxembourg -> Low Countries
    "Germany":     "Low Countries",
}

famine_years = (
    famine[famine["Famine"] == 1]
    .groupby("Region")["Year"].nunique()
    .rename("n_famine")
)

# Price: one row per location
price_locs = (
    grain.rename(columns={"Latitude": "lat", "Longitude": "lon"})
    .groupby("Location")
    .agg(year_min=("Year", "min"), year_max=("Year", "max"),
         teleco=("teleco_PDSI_10", "max"))
    .reset_index()
)
price_locs["Region"]  = price_locs["Location"].map(lambda x: PRICE_LOC.get(x, (None,None))[0])
price_locs["Country"] = price_locs["Location"].map(lambda x: PRICE_LOC.get(x, (None,None))[1])

# Yield: one row per lat/lon
yield_locs = (
    harvest
    .groupby(["Country", "Latitude", "Longitude"])
    .agg(year_min=("Year", "min"), year_max=("Year", "max"),
         teleco=("teleco_PDSI_10", "max"))
    .reset_index()
    .rename(columns={"Country": "Country"})
)
yield_locs["Region"] = yield_locs["Country"].map(YIELD_REGION_MAP)

# Aggregate price by region+country
price_by_country = (
    price_locs[price_locs["Region"].notna()]
    .groupby(["Region", "Country"])
    .agg(p_locs=("Location", "count"),
         p_ymin=("year_min", "min"),
         p_ymax=("year_max", "max"),
         p_teleco=("teleco", "sum"))
    .reset_index()
)

# Aggregate yield by region+country
yield_by_country = (
    yield_locs[yield_locs["Region"].notna()]
    .groupby(["Region", "Country"])
    .agg(y_locs=("Latitude", "count"),
         y_ymin=("year_min", "min"),
         y_ymax=("year_max", "max"),
         y_teleco=("teleco", "sum"))
    .reset_index()
)

# Merge on region+country (outer so we keep entries with only price or only yield)
merged = price_by_country.merge(yield_by_country, on=["Region","Country"], how="outer")

REGIONS = [
    "Central Europe", "France", "Great Britain", "Ireland",
    "Italy", "Low Countries", "Nordic Countries", "Russia/Ukraine", "Spain",
]
merged["_region_order"] = merged["Region"].map({r: i for i, r in enumerate(REGIONS)})
merged = merged.sort_values(["_region_order", "Country"]).reset_index(drop=True)
merged
"""

# ── Cell 6: write LaTeX ───────────────────────────────────────────────────────
SRC_LATEX = r"""
# ── Write granular LaTeX table ────────────────────────────────────────────────
TAB_DIR = f"{ROOT}/analysis/output/tables"
os.makedirs(TAB_DIR, exist_ok=True)

def fmt(val, fallback="---"):
    if val is None or (hasattr(val, '__float__') and __import__('math').isnan(float(val))):
        return fallback
    return str(int(val))

def yrange(ymin, ymax):
    if val_ok(ymin) and val_ok(ymax):
        return f"{int(ymin)}--{int(ymax)}"
    return "---"

def val_ok(v):
    import math
    try: return not math.isnan(float(v))
    except: return False

lines = []
lines.append(r"\begin{table}[ht]")
lines.append(r"\centering")
lines.append(r"\scriptsize")
lines.append(
    r"\caption{Data coverage by famine region and country (1500--1800). "
    r"Locs = number of distinct geographic locations; "
    r"TC = locations teleconnected to ENSO via PDSI ($p<0.10$); "
    r"Years = earliest--latest observation. "
    r"Famine years refer to the region-level panel. "
    r"Germany price sites in the Low Countries row are Xanten (Germany). "
    r"Germany yield sites (Remich, Grevenmacher) are located in Luxembourg.}"
)
lines.append(r"\label{tab:data_coverage_granular}")
lines.append(r"\begin{tabular}{ll r rrl rrl}")
lines.append(r"\toprule")
lines.append(
    r"Region & Country & Famine & \multicolumn{3}{c}{Grain price} & \multicolumn{3}{c}{Yield} \\"
)
lines.append(r" & & yrs & Locs & TC & Years & Locs & TC & Years \\")
lines.append(r"\midrule")

REGIONS = [
    "Central Europe", "France", "Great Britain", "Ireland",
    "Italy", "Low Countries", "Nordic Countries", "Russia/Ukraine", "Spain",
]

famine_printed = set()

for region in REGIONS:
    sub = merged[merged["Region"] == region]
    n_fam = int(famine_years.get(region, 0))

    if len(sub) == 0:
        # region with no data
        fam_str = str(n_fam) if region not in famine_printed else ""
        famine_printed.add(region)
        lines.append(f"{region} & --- & {fam_str} & --- & --- & --- & --- & --- & --- \\\\")
        continue

    for i, (_, row) in enumerate(sub.iterrows()):
        reg_str = region if i == 0 else ""
        fam_str = str(n_fam) if i == 0 else ""
        country = str(row["Country"]) if val_ok(row.get("Country")) else "---"

        pl  = fmt(row.get("p_locs"))
        pt  = fmt(row.get("p_teleco"))
        pyr = yrange(row.get("p_ymin"), row.get("p_ymax"))

        yl  = fmt(row.get("y_locs"))
        yt  = fmt(row.get("y_teleco"))
        yyr = yrange(row.get("y_ymin"), row.get("y_ymax"))

        lines.append(f"{reg_str} & {country} & {fam_str} & {pl} & {pt} & {pyr} & {yl} & {yt} & {yyr} \\\\")

    lines.append(r"\midrule")

# Remove last \midrule and add \bottomrule
lines[-1] = r"\bottomrule"
lines.append(r"\end{tabular}")
lines.append(r"\end{table}")

out_tex = f"{TAB_DIR}/data_coverage_table.tex"
with open(out_tex, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print(f"Saved -> {out_tex}")
"""

# Replace cells 5 and 6
nb["cells"][5] = new_cell(SRC_BUILD)
nb["cells"][6] = new_cell(SRC_LATEX)

with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
print(f"Done. Notebook has {len(nb['cells'])} cells.")
