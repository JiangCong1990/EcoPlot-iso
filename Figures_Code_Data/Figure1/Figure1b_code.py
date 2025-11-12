import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sklearn.linear_model import LinearRegression
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# === Load NAO index data ===
nao_path = r'/data/scratch/jiangcong/P/North_Atlantic_Oscillation.csv'
nao_data = pd.read_csv(nao_path, delimiter=",", parse_dates=[["date", "MONTH"]])
nao_data.rename(columns={'date_MONTH': 'date', 'INDEX': 'NAO'}, inplace=True)
nao_data['year'] = nao_data['date'].dt.year
nao_data['month'] = nao_data['date'].dt.month

# === Calculate DJFM average NAO per winter ===
djfm_years = []
djfm_nao_values = []
for year in range(1978, 2020):
    djfm = nao_data[((nao_data['year'] == year) & (nao_data['month'] == 12)) |
                    ((nao_data['year'] == year + 1) & (nao_data['month'].isin([1, 2])))]
    if len(djfm) == 3:
        djfm_years.append(year + 1)
        djfm_nao_values.append(djfm['NAO'].mean())
nao_series = pd.Series(djfm_nao_values, index=djfm_years, name='NAO_DJF')

# === Load precipitation data ===
ds = xr.open_dataset("./Monthly/europe_output.nc")
precip = ds['precipitation']
precip = precip.sel(time=slice('1979-06', '2019-08'))

# === Assign summer year ===
def assign_summer_year(t):
    t = pd.to_datetime(t)
    return t.year if t.month in [6, 7, 8] else None

summer_years = [assign_summer_year(t) for t in precip.time.values]
precip = precip.assign_coords(summer_year=('time', summer_years))
precip_summer = precip.sel(time=precip['time.month'].isin([6, 7, 8]))
precip_summer_mean = precip_summer.groupby('summer_year').mean('time')

# === Align with NAO series (DJFM previous winter) ===
precip_summer_mean = precip_summer_mean.sel(summer_year=slice(1979, 2019))
nao_series = nao_series.loc[1979:2019]

# === Detrend helper ===
def detrend_series(y):
    """Remove linear trend from a 1D array-like, return residuals."""
    x = np.arange(len(y)).reshape(-1, 1)
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return y
    model = LinearRegression().fit(x[mask], y[mask])
    trend = model.predict(x)
    return y - trend

# === Detrend NAO ===
nao_values = detrend_series(nao_series.values)

# === Compute correlation ===
r_map = np.full((precip_summer_mean.sizes['lat'], precip_summer_mean.sizes['lon']), np.nan)
p_map = np.full_like(r_map, np.nan)

for i in range(precip_summer_mean.sizes['lat']):
    for j in range(precip_summer_mean.sizes['lon']):
        y = precip_summer_mean[:, i, j].values
        if not np.any(np.isnan(y)):
            y_detr = detrend_series(y)   # detrend precipitation
            if np.all(np.isnan(y_detr)):
                continue
            r, p = pearsonr(y_detr, nao_values)
            r_map[i, j] = r
            p_map[i, j] = p

# === Build categorical map: 4 classes based on sign and significance ===
# Codes (for cmap index):
#   0 = non-significant negative (lightblue)
#   1 = significant negative (blue)
#   2 = significant positive (red)
#   3 = non-significant positive (lightcoral)

class_map = np.full_like(r_map, np.nan)

for i in range(r_map.shape[0]):
    for j in range(r_map.shape[1]):
        r = r_map[i, j]
        p = p_map[i, j]
        if np.isnan(r) or np.isnan(p):
            continue
        if r > 0:
            if p < 0.05:
                class_map[i, j] = 2
            else:
                class_map[i, j] = 3
        elif r < 0:
            if p < 0.05:
                class_map[i, j] = 1
            else:
                class_map[i, j] = 0

# === Define discrete colormap ===
cmap = ListedColormap(["lightgreen", "darkgreen", "red", "yellow"])

# === Plot ===
plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

im = ax.pcolormesh(
    precip_summer_mean.lon, precip_summer_mean.lat, class_map,
    cmap=cmap, transform=ccrs.PlateCarree()
)

# === Add features ===
ocean_mask = cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='face', facecolor='none')
ax.add_feature(ocean_mask, facecolor='white', zorder=2)
ax.coastlines(zorder=3)
ax.add_feature(cfeature.BORDERS, zorder=3)

# DMC location
dmc_lon, dmc_lat = 14.25, 52.3833
ax.scatter(dmc_lon, dmc_lat, marker='o', facecolors='none', edgecolors='black',
           s=20, transform=ccrs.PlateCarree(), zorder=4)
ax.text(dmc_lon + 1, dmc_lat, 'DMC', transform=ccrs.PlateCarree(), fontsize=10, color='black')

ax.text(0.05, 0.95, "(b)", transform=ax.transAxes, fontsize=10, fontweight='bold')

# Europe extent + ticks
import matplotlib.ticker as mticker
ax.set_extent([-12, 45, 35, 72], crs=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False; gl.right_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(-10, 50, 5))
gl.ylocator = mticker.FixedLocator(np.arange(35, 75, 5))
gl.xlabel_style = {'size': 8}; gl.ylabel_style = {'size': 8}

# === Custom legend ===
legend_elements = [
    Patch(facecolor=cmap.colors[2], edgecolor='k', label="Significant positive (p<0.05)"),
    Patch(facecolor=cmap.colors[3], edgecolor='k', label="Non-significant positive (p≥0.05)"),
    Patch(facecolor=cmap.colors[1], edgecolor='k', label="Significant negative (p<0.05)"),
    Patch(facecolor=cmap.colors[0], edgecolor='k', label="Non-significant negative (p≥0.05)")
]
#ax.legend(handles=legend_elements, ncol=2, loc='lower center', fontsize=8)

ax.legend(
    handles=legend_elements,
    ncol=2,
    loc='upper center',
    bbox_to_anchor=(0.5, -0.05),  # move legend below the figure
    fontsize=10,
    frameon=False
)

plt.tight_layout()
plt.savefig("winterNAO_vs_summerPrecip_discrete_colors.png", dpi=300)
plt.show()




# === Save plotting data (r_map, p_map, class_map) ===
import xarray as xr
import os

# Create output dataset
ds_out = xr.Dataset(
    {
        "correlation": (["lat", "lon"], r_map),
        "p_value": (["lat", "lon"], p_map),
        "class_map": (["lat", "lon"], class_map)
    },
    coords={
        "lat": precip_summer_mean.lat,
        "lon": precip_summer_mean.lon
    },
    attrs={
        "description": "Correlation between winter NAO (DJFM) and summer precipitation (JJA) for 1979–2019",
        "note": "class_map: 0=non-sig negative, 1=sig negative, 2=sig positive, 3=non-sig positive"
    }
)

# Define save directory
output_dir = "./plotting_data"
os.makedirs(output_dir, exist_ok=True)

# Save to NetCDF
nc_path = os.path.join(output_dir, "NAO_summer_precip_plotdata.nc")
ds_out.to_netcdf(nc_path)
print(f"✅ Plotting data saved to: {nc_path}")

# (Optional) Save small CSV for quick viewing of correlation values
df_corr = pd.DataFrame(
    {
        "lat": np.repeat(precip_summer_mean.lat.values, len(precip_summer_mean.lon)),
        "lon": np.tile(precip_summer_mean.lon.values, len(precip_summer_mean.lat)),
        "r": r_map.flatten(),
        "p": p_map.flatten(),
        "class": class_map.flatten()
    }
)
df_corr.to_csv(os.path.join(output_dir, "NAO_summer_precip_plotdata.csv"), index=False)
print("✅ CSV summary also saved.")