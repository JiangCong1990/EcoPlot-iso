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
import matplotlib.ticker as mticker

# === Helper function: detrend ===
def detrend_series(y):
    """Remove linear trend from a 1D array-like, return residuals."""
    x = np.arange(len(y)).reshape(-1, 1)
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return y
    model = LinearRegression().fit(x[mask], y[mask])
    trend = model.predict(x)
    return y - trend

# === Load NAO index data ===
nao_path = r'/data/scratch/jiangcong/P/North_Atlantic_Oscillation.csv'
nao_data = pd.read_csv(nao_path, delimiter=",", parse_dates=[["date", "MONTH"]])
nao_data.rename(columns={'date_MONTH': 'date', 'INDEX': 'NAO'}, inplace=True)
nao_data['year'] = nao_data['date'].dt.year
nao_data['month'] = nao_data['date'].dt.month

# === Compute DJF (winter-year) NAO means ===
djf_years = []
djf_nao_values = []
for year in range(1978, 2021):
    djf = nao_data[((nao_data['year'] == year) & (nao_data['month'] == 12)) |
                   ((nao_data['year'] == year + 1) & (nao_data['month'].isin([1, 2])))]
    if len(djf) == 3:
        djf_years.append(year + 1)
        djf_nao_values.append(djf['NAO'].mean())
nao_series = pd.Series(djf_nao_values, index=djf_years, name='NAO_DJF')

# === Load precipitation data ===
ds = xr.open_dataset("./Monthly/europe_output.nc")
precip = ds['precipitation']

# === Assign summer year (JJA) ===
def assign_summer_year(t):
    t = pd.to_datetime(t)
    return t.year if t.month in [6, 7, 8] else None

summer_years = [assign_summer_year(t) for t in precip.time.values]
precip = precip.assign_coords(summer_year=('time', summer_years))
precip_summer = precip.sel(time=precip['time.month'].isin([6, 7, 8]))
precip_summer_mean = precip_summer.groupby('summer_year').mean('time')

# 🔧 Fix: ensure summer_year is integer for slicing
precip_summer_mean = precip_summer_mean.assign_coords(
    summer_year=precip_summer_mean['summer_year'].astype(int)
)

# === Define periods ===
periods = {
    "1979–1999": (1979, 1999),
    "2000–2020": (2000, 2020)
}

# === Define discrete colormap ===
cmap = ListedColormap(["lightgreen", "darkgreen", "red", "yellow"])

# === Store maps for side-by-side plotting ===
maps = {}

for label, (start, end) in periods.items():
    # Align years
    common_years = np.intersect1d(
        precip_summer_mean['summer_year'].sel(summer_year=slice(start, end)).values,
        nao_series.loc[start:end].index.values
    )
    nao_sub = nao_series.loc[common_years]
    precip_sub = precip_summer_mean.sel(summer_year=common_years)

    # Detrend NAO
    nao_detr = detrend_series(nao_sub.values)

    # Correlation
    r_map = np.full((precip_sub.sizes['lat'], precip_sub.sizes['lon']), np.nan)
    p_map = np.full_like(r_map, np.nan)

    for i in range(precip_sub.sizes['lat']):
        for j in range(precip_sub.sizes['lon']):
            y = precip_sub[:, i, j].values
            if not np.any(np.isnan(y)):
                y_detr = detrend_series(y)
                if np.all(np.isnan(y_detr)) or np.std(y_detr) == 0:
                    continue
                r, p = pearsonr(y_detr, nao_detr)
                r_map[i, j] = r
                p_map[i, j] = p

    # === Build categorical map ===
    class_map = np.full_like(r_map, np.nan)
    for i in range(r_map.shape[0]):
        for j in range(r_map.shape[1]):
            r = r_map[i, j]
            p = p_map[i, j]
            if np.isnan(r) or np.isnan(p):
                continue
            if r > 0:
                class_map[i, j] = 2 if p < 0.05 else 3
            elif r < 0:
                class_map[i, j] = 1 if p < 0.05 else 0

    maps[label] = (class_map, precip_sub.lon, precip_sub.lat)

# === Plot two panels side by side ===
fig, axes = plt.subplots(1, 2, figsize=(12, 6),
                         subplot_kw={'projection': ccrs.PlateCarree()})

titles = list(periods.keys())
# === Panel titles and labels ===
panel_labels = ['(a)', '(b)', '(c)', '(d)'] 

for idx, ax in enumerate(axes):
    label = titles[idx]
    class_map, lon, lat = maps[label]

    im = ax.pcolormesh(
        lon, lat, class_map, cmap=cmap, transform=ccrs.PlateCarree()
    )

    # Decorations
    ocean_mask = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                              edgecolor='face', facecolor='none')
    ax.add_feature(ocean_mask, facecolor='white', zorder=2)
    ax.coastlines(zorder=3)
    ax.add_feature(cfeature.BORDERS, zorder=3)

    # DMC location
    dmc_lon, dmc_lat = 14.25, 52.3833
    ax.scatter(dmc_lon, dmc_lat, marker='o', facecolors='none', edgecolors='black',
               s=25, transform=ccrs.PlateCarree(), zorder=4)
    ax.text(dmc_lon + 1, dmc_lat, 'DMC', transform=ccrs.PlateCarree(),
            fontsize=10, color='black')

    ax.set_extent([-12, 45, 35, 72], crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5,
                      linestyle='--')
    gl.top_labels = False; gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-10, 50, 5))
    gl.ylocator = mticker.FixedLocator(np.arange(35, 75, 5))
    gl.xlabel_style = {'size': 8}; gl.ylabel_style = {'size': 8}
    # Add subplot panel label (upper-left corner)
    ax.text(-11, 71, panel_labels[idx], transform=ccrs.PlateCarree(),
            fontsize=12, fontweight='bold', va='top', ha='left')

    #ax.set_title(f"({chr(98+idx)}) {label}", fontsize=11, loc='left', fontweight='bold')
    ax.set_title(f"Winter NAO vs Summer Precipitation ({label})", 
             fontsize=12, loc='left', fontweight='bold')

# === Common legend below ===
legend_elements = [
    Patch(facecolor=cmap.colors[2], edgecolor='k', label="Significant positive (p<0.05)"),
    Patch(facecolor=cmap.colors[3], edgecolor='k', label="Non-significant positive (p≥0.05)"),
    Patch(facecolor=cmap.colors[1], edgecolor='k', label="Significant negative (p<0.05)"),
    Patch(facecolor=cmap.colors[0], edgecolor='k', label="Non-significant negative (p≥0.05)")
]
#fig.legend(handles=legend_elements, ncol=2, loc='upper center',
#           bbox_to_anchor=(0.5, -0.03), fontsize=9, frameon=False)

plt.tight_layout()
plt.subplots_adjust(bottom=0.12)
plt.savefig("winterNAO_vs_summerPrecip_discrete_1979-1999_2000-2020_detrended.png",
            dpi=300, bbox_inches="tight")
plt.show()

