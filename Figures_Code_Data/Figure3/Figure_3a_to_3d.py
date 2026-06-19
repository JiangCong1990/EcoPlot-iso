import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
import os

# =========================
# Helper function: detrend
# =========================
def detrend_series(y):
    """Remove linear trend from a 1D array-like, return residuals."""
    y = np.asarray(y, dtype=float)
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return y
    x = np.arange(len(y)).reshape(-1, 1)
    model = LinearRegression().fit(x[mask], y[mask])
    trend = model.predict(x)
    return y - trend


# =========================
# Load NAO index data
# =========================
nao_path = r'/data/scratch/jiangcong/P/North_Atlantic_Oscillation.csv'
nao_data = pd.read_csv(nao_path, delimiter=",")

# keep code B style because it is more explicit and robust
nao_data["date"] = pd.to_datetime(
    nao_data["date"].astype(str) + "-" + nao_data["MONTH"].astype(str) + "-15",
    errors="coerce", format="%Y-%m-%d"
)
nao_data.rename(columns={'INDEX': 'NAO'}, inplace=True)
nao_data['year'] = nao_data['date'].dt.year
nao_data['month'] = nao_data['date'].dt.month

# =========================
# Compute DJF (winter-year) NAO means
# =========================
djf_years = []
djf_nao_values = []
for year in range(1978, 2026):
    djf = nao_data[
        ((nao_data['year'] == year) & (nao_data['month'] == 12)) |
        ((nao_data['year'] == year + 1) & (nao_data['month'].isin([1, 2])))
    ]
    if len(djf) == 3:
        djf_years.append(year + 1)
        djf_nao_values.append(djf['NAO'].mean())

nao_series = pd.Series(
    djf_nao_values,
    index=pd.Index(djf_years, name="year"),
    name='NAO_DJF'
)

print(f"[INFO] NAO DJF series prepared: {nao_series.index.min()} → {nao_series.index.max()}")


# =========================
# Common settings
# =========================
cmap = ListedColormap(["lightgreen", "darkgreen", "red", "yellow"])
dmc_lon, dmc_lat = 14.25, 52.3833

save_dir = "plotting_data"
os.makedirs(save_dir, exist_ok=True)

legend_elements = [
    Patch(facecolor=cmap.colors[2], edgecolor='k', label="Significant positive (p<0.05)"),
    Patch(facecolor=cmap.colors[3], edgecolor='k', label="Non-significant positive (p≥0.05)"),
    Patch(facecolor=cmap.colors[1], edgecolor='k', label="Significant negative (p<0.05)"),
    Patch(facecolor=cmap.colors[0], edgecolor='k', label="Non-significant negative (p≥0.05)")
]

# =========================
# PART A: Winter NAO vs Summer Precipitation
# =========================
print("[INFO] Processing precipitation correlation maps ...")

# === Load precipitation data ===
ds_precip = xr.open_dataset("./Monthly/new/Monthly/europe_output_0p25deg.nc")
precip = ds_precip['precipitation']

# === Assign summer year (JJA) ===
def assign_summer_year(t):
    t = pd.to_datetime(t)
    return t.year if t.month in [6, 7, 8] else None

summer_years = [assign_summer_year(t) for t in precip.time.values]
precip = precip.assign_coords(summer_year=('time', summer_years))
precip_summer = precip.sel(time=precip['time.month'].isin([6, 7, 8]))
precip_summer_mean = precip_summer.groupby('summer_year').mean('time')

# ensure summer_year is integer
precip_summer_mean = precip_summer_mean.assign_coords(
    summer_year=precip_summer_mean['summer_year'].astype(int)
)

# === Define precipitation analysis periods ===
periods_precip = {
    "1979–1999": (1979, 1999),
    "2000–2024": (2000, 2024)
}

# === Store correlation maps ===
r_all_precip, p_all_precip, class_all_precip, period_labels_precip = [], [], [], []

for label, (start, end) in periods_precip.items():
    print(f"[INFO] Processing precipitation for {label} ...")

    common_years = np.intersect1d(
        precip_summer_mean['summer_year'].sel(summer_year=slice(start, end)).values,
        nao_series.loc[start:end].index.values
    )

    nao_sub = nao_series.loc[common_years]
    precip_sub = precip_summer_mean.sel(summer_year=common_years).transpose("summer_year", "lat", "lon")

    # Detrend NAO
    nao_detr = detrend_series(nao_sub.values)

    # Prepare arrays
    r_map = np.full((precip_sub.sizes['lat'], precip_sub.sizes['lon']), np.nan)
    p_map = np.full_like(r_map, np.nan)

    # Pixel-wise correlation
    for i in range(precip_sub.sizes['lat']):
        for j in range(precip_sub.sizes['lon']):
            y = precip_sub.isel(lat=i, lon=j).values
            if np.any(np.isnan(y)):
                continue
            y_detr = detrend_series(y)
            if np.all(np.isnan(y_detr)) or np.std(y_detr) == 0:
                continue
            if len(y_detr) != len(nao_detr):
                continue
            r, p = pearsonr(y_detr, nao_detr)
            r_map[i, j] = r
            p_map[i, j] = p

    # Build categorical map
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

    r_all_precip.append(r_map)
    p_all_precip.append(p_map)
    class_all_precip.append(class_map)
    period_labels_precip.append(label)

# === Save precipitation correlations to NetCDF ===
nc_out_path_precip = os.path.join(save_dir, "winterNAO_summerPrecip_correlations.nc")

r_da_precip = xr.DataArray(
    data=np.stack(r_all_precip),
    dims=("period", "lat", "lon"),
    coords={
        "period": period_labels_precip,
        "lat": precip_summer_mean.lat,
        "lon": precip_summer_mean.lon,
    },
    name="correlation_r",
    attrs={"description": "Pearson correlation coefficient between detrended DJF NAO and JJA precipitation"}
)

p_da_precip = xr.DataArray(
    data=np.stack(p_all_precip),
    dims=("period", "lat", "lon"),
    coords={
        "period": period_labels_precip,
        "lat": precip_summer_mean.lat,
        "lon": precip_summer_mean.lon
    },
    name="correlation_p",
    attrs={"description": "Two-tailed p-values for Pearson correlation"}
)

class_da_precip = xr.DataArray(
    data=np.stack(class_all_precip),
    dims=("period", "lat", "lon"),
    coords={
        "period": period_labels_precip,
        "lat": precip_summer_mean.lat,
        "lon": precip_summer_mean.lon
    },
    name="correlation_class",
    attrs={
        "description": "Categorical significance map: 0=ns neg, 1=sig neg, 2=sig pos, 3=ns pos",
        "note": "Significance threshold p<0.05"
    }
)

ds_out_precip = xr.Dataset({"r": r_da_precip, "p": p_da_precip, "class": class_da_precip})
ds_out_precip.attrs["title"] = "Winter NAO vs Summer Precipitation Correlation Maps"
ds_out_precip.attrs["note"] = "Detrended correlation analysis for 1979–1999 and 2000–2020"
ds_out_precip.to_netcdf(nc_out_path_precip)
print(f"✅ Precipitation correlation maps saved to: {nc_out_path_precip}")


# =========================
# PART B: Winter NAO vs SPEI-3
# =========================
print("[INFO] Processing SPEI-3 correlation maps ...")

# === Load SPEI-3 dataset ===
ds_spei = xr.open_dataset("/data/scratch/jiangcong/ERA5/monthly/era5_spi_spei_europe/era5_europe_spei_gamma_03.nc")
spei = ds_spei["spei_gamma_03"].sel(time=slice("1979-01", None))
print(f"[INFO] SPEI dataset loaded, time range: {str(spei.time.values[0])} → {str(spei.time.values[-1])}")

# === Define seasons ===
seasons = {"DJF": [12, 1, 2], "MAM": [3, 4, 5], "JJA": [6, 7, 8], "SON": [9, 10, 11]}
season_names = {"DJF": "Winter", "MAM": "Spring", "JJA": "Summer", "SON": "Autumn"}

# === Representative monthly SPEI-3 for each season ===
# SPEI-3 is already a 3-month accumulated index, so use the ending month
seasonal_means = {
    "DJF": spei.sel(time=spei["time.month"] == 2).groupby("time.year").first("time"),   # Feb -> Dec-Jan-Feb
    "MAM": spei.sel(time=spei["time.month"] == 5).groupby("time.year").first("time"),   # May -> Mar-Apr-May
    "JJA": spei.sel(time=spei["time.month"] == 8).groupby("time.year").first("time"),   # Aug -> Jun-Jul-Aug
    "SON": spei.sel(time=spei["time.month"] == 11).groupby("time.year").first("time")   # Nov -> Sep-Oct-Nov
}

# === Analysis periods ===
periods_spei = {
    "1979–1999": (1979, 1999),
    "2000–2025": (2000, 2025)
}

season_title_map = {
    "DJF": "February SPEI-3",
    "JJA": "August SPEI-3"
}

# === Prepare to save plotting data ===
nc_path_spei = os.path.join(save_dir, "winterNAO_vs_SPEI3_correlations.nc")

r_all_spei, p_all_spei, class_all_spei = [], [], []
season_list_spei, period_list_spei = [], []

# === Compute and collect correlation maps ===
for season in ["JJA"]:
#for season in ["DJF", "JJA"]:
    if season not in seasonal_means:
        continue

    spei_season = seasonal_means[season]

    for label, (start, end) in periods_spei.items():
        print(f"[INFO] Processing {season} for {label} ...")

        spei_years = spei_season["year"].values
        common_years = np.intersect1d(nao_series.index.values, spei_years)
        years_period = [y for y in common_years if (y >= start) and (y <= end)]

        if len(years_period) < 5:
            print(f"[WARN] Not enough overlapping years for {season} {label}")
            continue

        nao_common = nao_series.loc[years_period]
        nao_detr = detrend_series(nao_common.values)
        spei_period = spei_season.sel(year=years_period).transpose("year", "lat", "lon")

        r_map = np.full((spei_period.sizes["lat"], spei_period.sizes["lon"]), np.nan)
        p_map = np.full_like(r_map, np.nan)

        for i in range(spei_period.sizes["lat"]):
            for j in range(spei_period.sizes["lon"]):
                y = spei_period.isel(lat=i, lon=j).values
                y_detr = detrend_series(y)

                if len(y_detr) != len(nao_detr):
                    continue

                valid = (~np.isnan(y_detr)) & (~np.isnan(nao_detr))
                if valid.sum() < 3:
                    continue

                y_valid = y_detr[valid]
                nao_valid = nao_detr[valid]

                if np.std(y_valid) == 0 or np.std(nao_valid) == 0:
                    continue

                r, p = pearsonr(y_valid, nao_valid)
                r_map[i, j], p_map[i, j] = r, p

        # classify by significance
        class_map = np.full_like(r_map, np.nan)
        for i in range(r_map.shape[0]):
            for j in range(r_map.shape[1]):
                r, p = r_map[i, j], p_map[i, j]
                if np.isnan(r) or np.isnan(p):
                    continue
                if r > 0:
                    class_map[i, j] = 2 if p < 0.05 else 3
                elif r < 0:
                    class_map[i, j] = 1 if p < 0.05 else 0

        r_all_spei.append(r_map)
        p_all_spei.append(p_map)
        class_all_spei.append(class_map)
        season_list_spei.append(season)
        period_list_spei.append(label)

# === Save SPEI correlations to NetCDF ===
r_da_spei = xr.DataArray(
    np.stack(r_all_spei),
    dims=("index", "lat", "lon"),
    coords={
        "index": np.arange(len(r_all_spei)),
        "lat": spei.lat,
        "lon": spei.lon,
        "season": ("index", season_list_spei),
        "period": ("index", period_list_spei)
    },
    name="correlation_r",
    attrs={"description": "Pearson correlation between detrended winter NAO (DJF mean) and seasonal SPEI-3"}
)

p_da_spei = xr.DataArray(
    np.stack(p_all_spei),
    dims=("index", "lat", "lon"),
    coords={
        "index": np.arange(len(p_all_spei)),
        "lat": spei.lat,
        "lon": spei.lon,
        "season": ("index", season_list_spei),
        "period": ("index", period_list_spei)
    },
    name="correlation_p",
    attrs={"description": "Two-tailed p-values for correlation"}
)

class_da_spei = xr.DataArray(
    np.stack(class_all_spei),
    dims=("index", "lat", "lon"),
    coords={
        "index": np.arange(len(class_all_spei)),
        "lat": spei.lat,
        "lon": spei.lon,
        "season": ("index", season_list_spei),
        "period": ("index", period_list_spei)
    },
    name="correlation_class",
    attrs={
        "description": "Categorical significance: 0=ns neg, 1=sig neg, 2=sig pos, 3=ns pos",
        "note": "Significance threshold p<0.05"
    }
)

ds_out_spei = xr.Dataset({"r": r_da_spei, "p": p_da_spei, "class": class_da_spei})
ds_out_spei.attrs["title"] = "Winter NAO vs Seasonal SPEI-3 Correlations (1979–2024)"
ds_out_spei.attrs["note"] = "Each index represents a season-period combination."
ds_out_spei.to_netcdf(nc_path_spei)
print(f"✅ SPEI-3 correlation maps saved to: {nc_path_spei}")


# =========================
# Combined plotting
# =========================
print("[INFO] Plotting combined multi-panel figure ...")

fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten()

# ---- first row: precipitation panels (a, b)
titles_precip = list(periods_precip.keys())
panel_labels_precip = ['(a)', '(b)']

for idx in range(len(class_all_precip)):
    ax = axes[idx]
    class_map = class_all_precip[idx]
    label = titles_precip[idx]

    im = ax.pcolormesh(
        precip_summer_mean.lon, precip_summer_mean.lat, class_map,
        cmap=cmap, transform=ccrs.PlateCarree()
    )

    ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder=2)
    ax.add_feature(cfeature.COASTLINE, zorder=3)
    ax.add_feature(cfeature.BORDERS, zorder=3)
    ax.scatter(dmc_lon, dmc_lat, marker='o', facecolors='none', edgecolors='black',
               s=25, transform=ccrs.PlateCarree(), zorder=4)
    ax.text(dmc_lon + 1, dmc_lat, 'DMC', transform=ccrs.PlateCarree(),
            fontsize=9, color='black')

    ax.set_extent([-12, 45, 35, 72], crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-10, 50, 5))
    gl.ylocator = mticker.FixedLocator(np.arange(35, 75, 5))
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

    ax.set_title(panel_labels_precip[idx]+f" Winter NAO vs Summer Precipitation ({label})",
                 fontsize=10, loc='left', fontweight='bold')


# ---- remaining four panels: SPEI panels (c-f)
panel_labels_spei = ['(c)', '(d)', '(e)', '(f)']
start_ax = 2

for idx in range(len(class_all_spei)):
    ax = axes[start_ax + idx]
    class_map = class_all_spei[idx]
    season = season_list_spei[idx]
    period = period_list_spei[idx]

    im = ax.pcolormesh(
        spei.lon, spei.lat, class_map,
        cmap=cmap, transform=ccrs.PlateCarree()
    )

    ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder=2)
    ax.add_feature(cfeature.COASTLINE, zorder=3)
    ax.add_feature(cfeature.BORDERS, zorder=3)
    ax.scatter(dmc_lon, dmc_lat, marker='o', facecolors='none', edgecolors='black',
               s=25, transform=ccrs.PlateCarree(), zorder=4)
    ax.text(dmc_lon + 1, dmc_lat, 'DMC', transform=ccrs.PlateCarree(),
            fontsize=9, color='black')


    ax.set_extent([-12, 45, 35, 72], crs=ccrs.PlateCarree())
    #ax.set_title(panel_labels_spei[idx]+ f" Winter NAO vs {season_names[season]} SPEI-3 ({period})",
    #             fontsize=10, ha='left', fontweight='bold')
    ax.set_title(panel_labels_spei[idx] + f" Winter NAO vs {season_title_map[season]} ({period})",
             fontsize=10, x=0.01, ha='left', fontweight='bold')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-10, 50, 5))
    gl.ylocator = mticker.FixedLocator(np.arange(35, 75, 5))
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

# common legend
fig.legend(handles=legend_elements, ncol=2, loc='lower center',
           bbox_to_anchor=(0.5, 0.02), fontsize=10, frameon=False)

plt.tight_layout()
plt.subplots_adjust(bottom=0.1)
plt.savefig("winterNAO_vs_precip_and_SPEI3_discrete_combined_a_to_d.png",
            dpi=300, bbox_inches="tight")
plt.show()

print("[INFO] Combined multi-panel figure saved successfully.")
