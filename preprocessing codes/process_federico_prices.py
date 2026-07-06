#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
import cftime
import cartopy.crs as ccrs
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm


# # Prices

# In[2]:


fsv1 = pd.read_excel('../data/FSV prices.xlsx', sheet_name='1508-1785')
#add _fsv1 to column names
fsv1.rename(columns={"Market" : "Year"}, inplace=True)

#fsv1 = fsv1.drop(columns=['Unnamed: 16', 'Number','CV','Average', 'SD'])
fsv1 = fsv1.drop(columns=['Unnamed: 32', 'Number'])
fsv1 = fsv1.loc[(fsv1['Year'] >= 1500) & (fsv1['Year'] <= 1800)]

fsv1 = fsv1.set_index('Year')

#fsv1 = (fsv1 - fsv1.mean()) / fsv1.std()
print(fsv1.columns)


# In[3]:


prices = fsv1.copy()

start = prices.index.min()
end = prices.index.max()
prices.plot(figsize=(12, 8), title=f'Allen-FSV {start}-{end}',subplots=True, layout=(8, 5), sharex=True)


# In[4]:


from geopy.geocoders import Nominatim
import time

cities = prices.columns.tolist()

geolocator = Nominatim(user_agent="city_geocoder")

city_to_coords = {}

for city in cities:
    try:
        location = geolocator.geocode(city)
        if location:
            city_to_coords[city] = (location.latitude, location.longitude)
            print(f"{city}: ({location.latitude}, {location.longitude})")
        else:
            print(f"{city}: Not found")
    except Exception as e:
        print(f"Error with {city}: {e}")
    time.sleep(1)  # To avoid getting blocked by the server



# # PDSI

# In[5]:


# Load OWDA and compute weighted PDSI
owda = xr.open_dataset("../data/owda.nc")
owda["PDSI_weighted"] = np.sqrt(np.cos(np.deg2rad(owda["lat"])) + 1e-6) * owda["pdsi"]

# Flatten the 3D (time, lat, lon) data into a long table of gridpoints
pdsi_grid = owda["PDSI_weighted"].stack(points=("lat", "lon")).reset_coords()

# Drop points where all values are NaN
valid_points = pdsi_grid.dropna(dim="points", how="all")

# Extract coordinates for KDTree
grid_coords = np.column_stack([valid_points["lat"].values, valid_points["lon"].values])
tree = cKDTree(grid_coords)

# Prepare city coordinates
city_coords = np.array(list(city_to_coords.values()))  # shape (n_cities, 2)

# Find nearest grid point for each city
distances, indices = tree.query(city_coords, k=1, workers = -1)

# For each city, extract the time series from the nearest grid point
selected_data = []
years = owda["time"].values

for i, (city, coord) in enumerate(city_to_coords.items()):
    grid_idx = indices[i]
    lat = valid_points["lat"].values[grid_idx]
    lon = valid_points["lon"].values[grid_idx]

    # Get the PDSI series from that grid point
    city_ts = owda["PDSI_weighted"].sel(lat=lat, lon=lon, method="nearest")
    city_ts = city_ts.assign_coords(Location=city)
    selected_data.append(city_ts)

# Combine the city time series (DataArrays) into a Dataset
pdsi_selected = xr.Dataset({
    'PDSI': xr.concat(selected_data, dim="loc_id")
})

# Rename time to Year
pdsi_selected = pdsi_selected.rename({"time": "Year"})

# Add domain averages as DataArrays (with the same time/Year dimension)
pdsi_selected["PDSI_europe"] = owda["PDSI_weighted"].sel(lat=slice(35, 70), lon=slice(-10, 40)).mean(dim=["lat", "lon"]).rename({"time": "Year"})
pdsi_selected["PDSI_southern"] = owda["PDSI_weighted"].sel(lat=slice(20, 40), lon=slice(-10, 20)).mean(dim=["lat", "lon"]).rename({"time": "Year"})

pdsi_df = pdsi_selected.to_dataframe().reset_index()
pdsi_df.rename(columns={"PDSI_weighted": "PDSI"}, inplace=True)


# # Temp

# In[6]:


# Load dataset
xds = xr.open_dataset('../data/ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc')

# Extract year and filter by given range
xds = xds.assign(Year=xds.time.dt.year)
xds = xds.sel(time=slice(str(1500), str(1800)))

# Apply latitude weighting
xds['temp'] = np.sqrt(np.cos(np.deg2rad(xds['latitude'])) + 1e-6) * (xds['temp2'])

winter_months = [11, 12, 1, 2]
summer_months = [4, 5, 6, 7]

def shift_decembers(time_array):
    return [cftime.DatetimeGregorian(t.year + 1, t.month, t.day) if t.month in [11, 12] else t for t in time_array]

xds = xds.assign_coords(time=("time", shift_decembers(xds.time.values)))



# --- **(1) Compute City-Level Temperatures** ---
def extract_location_data(xds_season):
    return [
        xds_season.sel(latitude=lat, longitude=lon, method='nearest')
        .to_dataframe().reset_index().assign(Location=city)
        for city, (lat, lon) in city_to_coords.items()
    ]

# Compute city-level seasonal means
xds_winter_cities = xds.sel(time=xds.time.dt.month.isin(winter_months)).groupby('Year').mean(dim='time')
xds_summer_cities = xds.sel(time=xds.time.dt.month.isin(summer_months)).groupby('Year').mean(dim='time')

temp_winter = pd.concat(extract_location_data(xds_winter_cities), ignore_index=True)
temp_summer = pd.concat(extract_location_data(xds_summer_cities), ignore_index=True)

# --- **(2) Compute European Domain-Average Temperatures** ---
# Define European bounding box (adjust as needed)
lat_bounds = (70, 35)   # Latitude range for Europe
lon_bounds = (-15, 35)  # Longitude range for Europe

# Subset dataset to Europe only
xds_europe = xds.sel(latitude=slice(*lat_bounds), longitude=slice(*lon_bounds))

# Function to compute spatially averaged seasonal mean
def compute_seasonal_mean(xds, months):
    return (xds.sel(time=xds.time.dt.month.isin(months))
            .groupby('Year').mean(dim='time')
            .mean(dim=['latitude', 'longitude']))

# Compute European domain-average temperatures
xds_winter_europe = compute_seasonal_mean(xds_europe, winter_months)
xds_summer_europe = compute_seasonal_mean(xds_europe, summer_months)

# Convert to DataFrame
europe_temp_df = pd.DataFrame({
    'Year': xds_winter_europe.Year.values,
    'temp_winter_europe': xds_winter_europe['temp'].values,
    'temp_summer_europe': xds_summer_europe['temp'].values
})

# --- **(3) Merge Everything Together** ---
temp_df = pd.merge(temp_winter, temp_summer, on=['Year', 'Location'], suffixes=('_winter', '_summer'))

# Drop unnecessary columns
temp_df = temp_df.drop(columns=['latitude_winter', 'longitude_winter', 'latitude_summer', 'longitude_summer', 'temp2_winter', 'temp2_summer'])

# Merge with European averages
temp_df = pd.merge(temp_df, europe_temp_df, on='Year', how='left')




# # PRecip

# In[7]:


# Load dataset
xds_p = xr.open_dataset('../data/ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc')

# Extract year and filter by given range
xds_p = xds_p.assign(Year=xds_p.time.dt.year)
xds_p = xds_p.sel(time=slice(str(1500), str(1800)))

# Apply latitude weighting
xds_p['precip'] = np.sqrt(np.cos(np.deg2rad(xds_p['latitude'])) + 1e-6) * xds_p['totprec']*86400*30

# Define seasons
winter_months = [11, 12, 1, 2]
summer_months = [4, 5, 6, 7]

def shift_decembers(time_array):
    return [cftime.DatetimeGregorian(t.year + 1, t.month, t.day) if t.month in [11, 12] else t for t in time_array]

xds_p = xds_p.assign_coords(time=("time", shift_decembers(xds_p.time.values)))



# --- **(1) Compute City-Level preciperatures** ---
def extract_location_data(xds_p_season):
    return [
        xds_p_season.sel(latitude=lat, longitude=lon, method='nearest')
        .to_dataframe().reset_index().assign(Location=city)
        for city, (lat, lon) in city_to_coords.items()
    ]

# Compute city-level seasonal means
xds_p_winter_cities = xds_p.sel(time=xds_p.time.dt.month.isin(winter_months)).groupby('Year').mean(dim='time')
xds_p_summer_cities = xds_p.sel(time=xds_p.time.dt.month.isin(summer_months)).groupby('Year').mean(dim='time')

precip_winter = pd.concat(extract_location_data(xds_p_winter_cities), ignore_index=True)
precip_summer = pd.concat(extract_location_data(xds_p_summer_cities), ignore_index=True)

# --- **(2) Compute European Domain-Average preciperatures** ---
# Define European bounding box (adjust as needed)
lat_bounds = (70, 35)   # Latitude range for Europe
lon_bounds = (-15, 35)  # Longitude range for Europe

# Subset dataset to Europe only
xds_p_europe = xds_p.sel(latitude=slice(*lat_bounds), longitude=slice(*lon_bounds))

# Function to compute spatially averaged seasonal mean
def compute_seasonal_mean(xds_p, months):
    return (xds_p.sel(time=xds_p.time.dt.month.isin(months))
            .groupby('Year').mean(dim='time')
            .mean(dim=['latitude', 'longitude']))

# Compute European domain-average preciperatures
xds_p_winter_europe = compute_seasonal_mean(xds_p_europe, winter_months)
xds_p_summer_europe = compute_seasonal_mean(xds_p_europe, summer_months)

# Convert to DataFrame
europe_precip_df = pd.DataFrame({
    'Year': xds_p_winter_europe.Year.values,
    'precip_winter_europe': xds_p_winter_europe['precip'].values,
    'precip_summer_europe': xds_p_summer_europe['precip'].values
})

# --- **(3) Merge Everything Together** ---
precip_df = pd.merge(precip_winter, precip_summer, on=['Year', 'Location'], suffixes=('_winter', '_summer'))

# Drop unnecessary columns
precip_df = precip_df.drop(columns=['latitude_winter', 'longitude_winter', 'latitude_summer', 'longitude_summer', 'totprec_winter', 'totprec_summer'])

# Merge with European averages
precip_df = pd.merge(precip_df, europe_precip_df, on='Year', how='left')


# In[8]:


precip_df


# # ENSO

# In[9]:


enso3 = pd.read_csv("../data/cook2024-R15-ENSO-Rec-1500-2000.txt", delimiter="\t", comment="#", na_values="NA")
lat_lon = pd.read_csv("../data/cook2024-ENSO-latlon.txt", delimiter="\t", comment="#", na_values="NA")


# In[10]:


# Melt the ENSO DataFrame to long format
enso_melted = enso3.melt(id_vars=["Year"], var_name="gridpoint", value_name="enso")
enso_melted["gridpoint"] = enso_melted["gridpoint"].astype(int)

# Merge with the lat/lon DataFrame
enso_merged = enso_melted.merge(lat_lon, on="gridpoint")



# Convert to xarray
enso_xr = enso_merged.set_index(["Year", "lat", "lon"])["enso"].to_xarray()

enso_xr = enso_xr  - enso_xr.sel(Year=slice(1801, 1900)).mean(dim="Year")

nino3 = enso_xr.sel(lat=slice(-5, 5), lon=slice(-150, -90)).mean(dim=["lat", "lon"])
nino34 = enso_xr.sel(lat=slice(-5, 5), lon=slice(-170, -120)).mean(dim=["lat", "lon"])
nino4 = enso_xr.sel(lat=slice(-5, 5), lon = slice(-200, -150)).mean(dim=["lat", "lon"])
nino12 = enso_xr.sel(lat=slice(-10, 0), lon=slice(-90, -80)).mean(dim=["lat", "lon"])


# In[11]:


enso = pd.merge(nino3.to_dataframe().reset_index().rename(columns={"enso": "nino3"}), nino34.to_dataframe().reset_index().rename(columns={"enso": "nino34"}), on="Year")
enso = pd.merge(enso, nino4.to_dataframe().reset_index().rename(columns={"enso": "nino4"}), on="Year")
enso = pd.merge(enso, nino12.to_dataframe().reset_index().rename(columns={"enso": "nino12"}), on="Year")


# In[12]:


enso.plot(x="Year", y=["nino3", "nino34"], figsize=(12, 6))


# In[13]:


#compute correlation of enso
enso.corr()


# # JSL and NAO

# In[14]:


jsl = pd.read_excel("../data/reconstructed EU JSL.xlsx")

nao_cal = pd.read_excel("../data/nao_reconst.xlsx", sheet_name="Figure 2a", header=3)
nao_cal = nao_cal.rename(columns={"Time (years AD)" : "Year"})
nao_model = pd.read_excel("../data/nao_reconst.xlsx", sheet_name="Figure 2b", header=3)
nao_model = nao_model.rename(columns={"Time (years AD)" : "Year"})

nao = pd.merge(nao_cal, nao_model, on="Year", suffixes=('_cal', '_model'))

nao = nao.rename(columns={"Ensemble Mean_cal" : "NAO_cal", "Ensemble Mean_model" : "NAO_model"})[["Year", "NAO_cal", "NAO_model"]]


# # COnflict

# In[15]:


# Load the data
post1400 = pd.read_excel('../data/Conflict-Catalog-18-vars.xlsx')
pre1400 = pd.read_excel('../data/Brecke-Pre-1400-European-Conflicts.xlsx')


# Convert StartYear and EndYear to integer if not NaN
post1400['StartYear'] = post1400['StartYear'].apply(lambda x: int(x) if not pd.isna(x) else np.nan)
post1400['EndYear'] = post1400['EndYear'].apply(lambda x: int(x) if not pd.isna(x) else np.nan)

# Filter Region to be 3 and 4
post1400 = post1400[post1400['Region'].isin([3, 4])]
# Select relevant columns
post1400 = post1400[['TotalFatalities', 'MilFatalities', 'StartYear', 'EndYear']].copy()
pre1400 = pre1400[['Fatalities', 'StartYear', 'EndYear']].copy()

# Fill NaN values in TotalFatalities and MilFatalities with 0
post1400['TotalFatalities'] = post1400['TotalFatalities'].fillna(0)
post1400['MilFatalities'] = post1400['MilFatalities'].fillna(0)
post1400['Deaths'] = post1400['TotalFatalities']

pre1400['Deaths'] = pre1400['Fatalities'].fillna(0)

# Drop TotalFatalities and MilFatalities
post1400 = post1400.drop(columns=['TotalFatalities', 'MilFatalities'])

pre1400 = pre1400.drop(columns=['Fatalities'])

wars = pd.concat([pre1400,post1400])

# Initialize an empty dictionary for tracking deaths per year
death_index = {}

# Populate the index with death counts, ongoing wars, started wars, and total duration
for _, row in wars.iterrows():
    if pd.isna(row['StartYear']) or pd.isna(row['EndYear']):
        continue
    
    year_range = range(int(row['StartYear']), int(row['EndYear']) + 1)
    deaths_per_year = row['Deaths'] / len(year_range)
    duration = len(year_range)
    
    for year in year_range:
        if year in death_index:
            death_index[year]['Deaths'] += deaths_per_year
            death_index[year]['ongoing_wars'] += 1
            if year == row['StartYear']:
                death_index[year]['started_wars'] += 1
                death_index[year]['total_duration'] += duration
        else:
            death_index[year] = {
                'Deaths': deaths_per_year,
                'ongoing_wars': 1,
                'started_wars': 1 if year == row['StartYear'] else 0,
                'total_duration': duration if year == row['StartYear'] else 0
            }

# Convert dictionary to DataFrame
conflict = pd.DataFrame.from_dict(death_index, orient='index').reset_index()

conflict['Death_ratio'] = conflict['Deaths'] / conflict['ongoing_wars']
conflict.rename(columns={'index': 'Year'}, inplace=True)
conflict = conflict.sort_values(by='Year').reset_index(drop=True).fillna(0)

conflict.set_index('Year')['Deaths'].plot()
plt.xlabel('Year')
plt.ylabel('Deaths')
plt.title('Yearly War Deaths')
plt.show()


# # Merge everyone and save

# In[16]:


df_total = pd.merge(temp_df, pdsi_df, on=['Year', 'Location']) 
df_total = pd.merge(df_total, precip_df, on=['Year', 'Location'])
df_total = pd.merge(df_total, enso, on=['Year'])
df_total = pd.merge(df_total, jsl, on=['Year'])
df_total = pd.merge(df_total, nao, on=['Year'])
prices['Year'] = prices.index
prices_long = prices.melt(id_vars=['Year'], var_name='Location', value_name='price')
df_total = pd.merge(df_total, prices_long, on=['Year', 'Location'])
df_total = df_total.sort_values(['Location', 'Year'])

df_total = pd.merge(df_total, conflict, on=['Year'], how='left')
df_total[['Deaths', 'ongoing_wars', 'started_wars', 'total_duration', 'Death_ratio']] = df_total[['Deaths', 'ongoing_wars', 'started_wars', 'total_duration', 'Death_ratio']].fillna(0)

# add coords to df_total
for city, coords in city_to_coords.items():
    df_total.loc[df_total['Location'] == city, 'lat'] = coords[0]
    df_total.loc[df_total['Location'] == city, 'lon'] = coords[1]

df_total['logprice'] = np.log(df_total['price'])
df_total['logprice'] = df_total['logprice'].replace([np.inf, -np.inf], np.nan)


# In[17]:


df_total


# In[18]:


df_total.describe()


# # Teleco

# In[19]:


import pandas as pd
from scipy.stats import pearsonr, kendalltau, spearmanr

# User-specified variables and significance levels
variables = ['PDSI', 'precip_summer','precip_winter','temp_winter','temp_summer']  # <- replace with your desired variables
significance_levels = [0.05, 0.01, 0.10]  # 5%, 1%, 10%

df_climate = pd.merge(temp_df, pdsi_df, on=['Year', 'Location']) 
df_climate = pd.merge(df_climate, precip_df, on=['Year', 'Location'])
df_climate = pd.merge(df_climate, enso, on=['Year'])

# Initialize result columns
for var in variables:
    for sig in significance_levels:
        colname = f"teleco_{var}_{int(sig*100):02d}"
        df_total[colname] = 0  # default value

# Loop over locations and variables
for loc_id, group in df_climate.groupby('Location'):
    for var in variables:
        valid_data = group[[var, 'nino34']].dropna()
        if len(valid_data) >= 2:
            corr, pval = pearsonr(valid_data[var], valid_data['nino34'])
            for sig in significance_levels:
                if pval < sig:
                    colname = f"teleco_{var}_{int(sig*100):02d}"
                    df_total.loc[df_total['Location'] == loc_id, colname] = 1


# In[20]:


df_total.to_csv('../processed data/federico_prices_enso.csv', index=False)


# In[ ]:




