"""
combined_winterNAO_vs_summer_hydroclimate_maps.py

Compute and plot pixel-wise Pearson correlation between detrended winter NAO (DJF mean)
and four summer hydroclimate variables over Europe:

(a) Summer precipitation (JJA mean)
(b) August SPEI-3
(c) Summer satellite soil moisture (JJA mean from monthly files)
(d) Summer GRACE terrestrial water storage / LWE thickness (JJA mean)

Output:
  - One combined 2x2 figure
  - One NetCDF file with r, p, and class maps for all four variables
"""

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

# ==========================================================
# SETTINGS
# ==========================================================
NAO_FILE = "/data/scratch/jiangcong/P/North_Atlantic_Oscillation.csv"

PRECIP_FILE = "./Monthly/new/Monthly/europe_output.nc"
PRECIP_VAR = "precipitation"

SPEI_FILE = "/data/scratch/jiangcong/ERA5/monthly/era5_spi_spei_europe/era5_europe_spei_gamma_03.nc"
SPEI_VAR = "spei_gamma_03"

SM_DIR = "/data/scratch/jiangcong/SM_data/data/SSMV-COMBINED_EU"

GRACE_FILE = "/data/scratch/jiangcong/grace_data/GRCTellus.JPL.200204_202512.GLO.RL06.3M.MSCNv04CRI_eu.nc"
GRACE_VAR = "lwe_thickness"

EU_EXTENT = [-12, 45, 35, 72]
SUMMER_MONTHS = [6, 7, 8]

MONTH_ABBR = {
    1: "JAN", 2: "FEB", 3: "MAR", 4: "APR", 5: "MAY", 6: "JUN",
    7: "JUL", 8: "AUG", 9: "SEP", 10: "OCT", 11: "NOV", 12: "DEC"
}

# 0=ns neg, 1=sig neg, 2=sig pos, 3=ns pos
cmap = ListedColormap(["lightgreen", "darkgreen", "red", "yellow"])

os.makedirs("plotting_data", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ==========================================================
# HELPER: detrend
# ==========================================================
def detrend_series(y):
    """Remove linear trend from a 1D array-like, return residuals."""
    x = np.arange(len(y)).reshape(-1, 1)
    y = np.asarray(y, dtype=float)
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return y
    model = LinearRegression().fit(x[mask], y[mask])
    trend = model.predict(x)
    return y - trend

# ==========================================================
# HELPER: remove low-frequency variability
# ==========================================================
def remove_low_frequency_component(y, window=11, center=True, min_periods=6):
    """
    Remove low-frequency variability from a 1D array-like series
    by subtracting an 11-year running mean.
    """
    y = np.asarray(y, dtype=float)
    y_series = pd.Series(y)
    low_freq = y_series.rolling(window=window, center=center, min_periods=min_periods).mean()
    return (y_series - low_freq).values

# ==========================================================
# HELPER: winter DJF NAO
# ==========================================================
def build_djf_nao_series(nao_file):
    nao_data = pd.read_csv(nao_file, delimiter=",")

    # keep this logic explicit to avoid parse_dates warning
    nao_data["date"] = pd.to_datetime(
        nao_data["date"].astype(str) + "-" + nao_data["MONTH"].astype(str) + "-15",
        format="%Y-%m-%d",
        errors="coerce"
    )

    nao_data.rename(columns={"INDEX": "NAO"}, inplace=True)
    nao_data["year"] = nao_data["date"].dt.year.astype(int)
    nao_data["month"] = nao_data["date"].dt.month.astype(int)

    djf_years = []
    djf_nao_values = []

    min_year = int(nao_data["year"].min())
    max_year = int(nao_data["year"].max())

    for year in range(min_year, max_year):
        djf = nao_data[
            ((nao_data["year"] == year) & (nao_data["month"] == 12)) |
            ((nao_data["year"] == year + 1) & (nao_data["month"].isin([1, 2])))
        ]
        if len(djf) == 3:
            djf_years.append(year + 1)   # winter-year
            djf_nao_values.append(djf["NAO"].mean())

    return pd.Series(djf_nao_values, index=pd.Index(djf_years, name="year"), name="NAO_DJF")

# ==========================================================
# HELPER: classify correlation map
# ==========================================================
def build_class_map(r_map, p_map):
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

    return class_map

# ==========================================================
# HELPER: pixel-wise correlation
# target_da must have dims (year, lat, lon)
# ==========================================================
def compute_corr_maps(nao_series, target_da):
    common_years = np.intersect1d(
        nao_series.index.values.astype(int),
        target_da["year"].values.astype(int)
    )

    if common_years.size < 5:
        raise ValueError(f"Not enough overlapping years for correlation: n={common_years.size}")

    nao_sub = nao_series.loc[common_years]
    target_sub = target_da.sel(year=common_years).transpose("year", "lat", "lon")

    nao_detr = detrend_series(nao_sub.values)

    r_map = np.full((target_sub.sizes["lat"], target_sub.sizes["lon"]), np.nan)
    p_map = np.full_like(r_map, np.nan)

    for i in range(target_sub.sizes["lat"]):
        for j in range(target_sub.sizes["lon"]):
            y = target_sub[:, i, j].values
            if np.all(np.isnan(y)):
                continue

            y_detr = detrend_series(y)
            valid = ~np.isnan(y_detr) & ~np.isnan(nao_detr)

            if valid.sum() < 5:
                continue
            if np.nanstd(y_detr[valid]) == 0 or np.nanstd(nao_detr[valid]) == 0:
                continue

            r, p = pearsonr(y_detr[valid], nao_detr[valid])
            r_map[i, j] = r
            p_map[i, j] = p

    class_map = build_class_map(r_map, p_map)
    return common_years, r_map, p_map, class_map

# ==========================================================
# HELPER: plot one panel
# ==========================================================
def plot_class_panel(ax, lon, lat, class_map, title, panel_label,
                     show_left_labels=True, show_bottom_labels=True):
    ax.pcolormesh(lon, lat, class_map, cmap=cmap, transform=ccrs.PlateCarree())

    ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder=2)
    ax.add_feature(cfeature.COASTLINE, zorder=3, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, zorder=3, linewidth=0.4)

    dmc_lon, dmc_lat = 14.25, 52.3833
    ax.scatter(
        dmc_lon, dmc_lat,
        marker='o', facecolors='none', edgecolors='black',
        s=20, transform=ccrs.PlateCarree(), zorder=4
    )
    ax.text(dmc_lon + 1, dmc_lat, "DMC", transform=ccrs.PlateCarree(),
            fontsize=9, color="black")

    ax.set_extent(EU_EXTENT, crs=ccrs.PlateCarree())
    ax.set_title(f"{panel_label} {title}", fontsize=11, loc="left", fontweight="bold")

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = show_left_labels
    gl.bottom_labels = show_bottom_labels
    gl.xlocator = mticker.FixedLocator(np.arange(-10, 50, 5))
    gl.ylocator = mticker.FixedLocator(np.arange(35, 75, 5))
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

# ==========================================================
# LOAD NAO
# ==========================================================
nao_series = build_djf_nao_series(NAO_FILE)

# ==========================================================
# A. PRECIPITATION: JJA mean
# ==========================================================
ds_precip = xr.open_dataset(PRECIP_FILE)
precip = ds_precip[PRECIP_VAR]

if not np.issubdtype(precip["time"].dtype, np.datetime64):
    precip["time"] = xr.decode_cf(ds_precip).time

precip = precip.sel(time=precip["time.month"].isin(SUMMER_MONTHS))
precip_summer_mean = precip.groupby("time.year").mean("time", skipna=True)
precip_summer_mean = precip_summer_mean.assign_coords(
    year=precip_summer_mean["year"].astype(int)
)

years_precip, precip_r, precip_p, precip_class = compute_corr_maps(
    nao_series, precip_summer_mean
)

# ==========================================================
# B. SPEI-3: AUGUST only
# ==========================================================
ds_spei = xr.open_dataset(SPEI_FILE)
spei = ds_spei[SPEI_VAR]

if not np.issubdtype(spei["time"].dtype, np.datetime64):
    spei["time"] = xr.decode_cf(ds_spei).time

spei_aug = spei.sel(time=spei["time.month"] == 8)
spei_aug = spei_aug.groupby("time.year").mean("time", skipna=True)
spei_aug = spei_aug.assign_coords(
    year=spei_aug["year"].astype(int)
)

years_spei, spei_r, spei_p, spei_class = compute_corr_maps(
    nao_series, spei_aug
)

# ==========================================================
# C. SATELLITE SOIL MOISTURE: JJA mean from JUN/JUL/AUG files
# ==========================================================
sm_monthly_list = []

for m in SUMMER_MONTHS:
    mon = MONTH_ABBR[m]
    sm_file = f"{SM_DIR}/SSMV_EU_{mon}_1978_2024.nc"
    if not os.path.exists(sm_file):
        raise FileNotFoundError(f"Soil moisture file not found: {sm_file}")

    ds_sm = xr.open_dataset(sm_file)
    sm = ds_sm["sm"]
    sm = sm.where(sm != -9999)

    if not np.issubdtype(sm["time"].dtype, np.datetime64):
        sm["time"] = xr.decode_cf(ds_sm).time

    sm_mean = sm.groupby("time.year").mean("time", skipna=True)
    sm_mean = sm_mean.assign_coords(year=sm_mean["year"].astype(int))
    sm_monthly_list.append(sm_mean)

sm_summer_mean = xr.concat(sm_monthly_list, dim="summer_month").mean("summer_month", skipna=True)
sm_summer_mean = sm_summer_mean.transpose("year", "lat", "lon")

years_sm, sm_r, sm_p, sm_class = compute_corr_maps(
    nao_series, sm_summer_mean
)

# ==========================================================
# D. GRACE: JJA mean
# ==========================================================
ds_grace = xr.open_dataset(GRACE_FILE)
grace = ds_grace[GRACE_VAR].where(ds_grace[GRACE_VAR] != -99999)

if not np.issubdtype(grace["time"].dtype, np.datetime64):
    grace["time"] = xr.decode_cf(ds_grace).time

if grace.lon.max() > 180:
    grace = grace.assign_coords(lon=((grace.lon + 180) % 360) - 180).sortby("lon")
    ds_grace = ds_grace.assign_coords(lon=((ds_grace.lon + 180) % 360) - 180).sortby("lon")

if "land_mask" in ds_grace.variables:
    land = ds_grace["land_mask"]
    if land.lon.max() > 180:
        land = land.assign_coords(lon=((land.lon + 180) % 360) - 180).sortby("lon")
    grace = grace.where(land > 0)

grace = grace.sel(time=grace["time.month"].isin(SUMMER_MONTHS))
grace_summer_mean = grace.groupby("time.year").mean("time", skipna=True)
grace_summer_mean = grace_summer_mean.assign_coords(
    year=grace_summer_mean["year"].astype(int)
)

years_grace, grace_r, grace_p, grace_class = compute_corr_maps(
    nao_series, grace_summer_mean
)

# ==========================================================
# SAVE ALL TO ONE NETCDF
# ==========================================================
ds_out = xr.Dataset()

ds_out["precip_r"] = xr.DataArray(
    precip_r,
    dims=("precip_lat", "precip_lon"),
    coords={
        "precip_lat": precip_summer_mean.lat.values,
        "precip_lon": precip_summer_mean.lon.values
    }
)
ds_out["precip_p"] = xr.DataArray(
    precip_p,
    dims=("precip_lat", "precip_lon"),
    coords={
        "precip_lat": precip_summer_mean.lat.values,
        "precip_lon": precip_summer_mean.lon.values
    }
)
ds_out["precip_class"] = xr.DataArray(
    precip_class,
    dims=("precip_lat", "precip_lon"),
    coords={
        "precip_lat": precip_summer_mean.lat.values,
        "precip_lon": precip_summer_mean.lon.values
    }
)

ds_out["spei_r"] = xr.DataArray(
    spei_r,
    dims=("spei_lat", "spei_lon"),
    coords={
        "spei_lat": spei_aug.lat.values,
        "spei_lon": spei_aug.lon.values
    }
)
ds_out["spei_p"] = xr.DataArray(
    spei_p,
    dims=("spei_lat", "spei_lon"),
    coords={
        "spei_lat": spei_aug.lat.values,
        "spei_lon": spei_aug.lon.values
    }
)
ds_out["spei_class"] = xr.DataArray(
    spei_class,
    dims=("spei_lat", "spei_lon"),
    coords={
        "spei_lat": spei_aug.lat.values,
        "spei_lon": spei_aug.lon.values
    }
)

ds_out["sm_r"] = xr.DataArray(
    sm_r,
    dims=("sm_lat", "sm_lon"),
    coords={
        "sm_lat": sm_summer_mean.lat.values,
        "sm_lon": sm_summer_mean.lon.values
    }
)
ds_out["sm_p"] = xr.DataArray(
    sm_p,
    dims=("sm_lat", "sm_lon"),
    coords={
        "sm_lat": sm_summer_mean.lat.values,
        "sm_lon": sm_summer_mean.lon.values
    }
)
ds_out["sm_class"] = xr.DataArray(
    sm_class,
    dims=("sm_lat", "sm_lon"),
    coords={
        "sm_lat": sm_summer_mean.lat.values,
        "sm_lon": sm_summer_mean.lon.values
    }
)

ds_out["grace_r"] = xr.DataArray(
    grace_r,
    dims=("grace_lat", "grace_lon"),
    coords={
        "grace_lat": grace_summer_mean.lat.values,
        "grace_lon": grace_summer_mean.lon.values
    }
)
ds_out["grace_p"] = xr.DataArray(
    grace_p,
    dims=("grace_lat", "grace_lon"),
    coords={
        "grace_lat": grace_summer_mean.lat.values,
        "grace_lon": grace_summer_mean.lon.values
    }
)
ds_out["grace_class"] = xr.DataArray(
    grace_class,
    dims=("grace_lat", "grace_lon"),
    coords={
        "grace_lat": grace_summer_mean.lat.values,
        "grace_lon": grace_summer_mean.lon.values
    }
)

ds_out.attrs["title"] = "Correlation between detrended winter NAO and summer hydroclimate variables"
ds_out.attrs["note"] = (
    "Precipitation = JJA mean; SPEI-3 = August only; "
    "Satellite soil moisture = JJA mean; GRACE LWE = JJA mean. "
    "class: 0=ns neg, 1=sig neg, 2=sig pos, 3=ns pos"
)

nc_path = "plotting_data/winterNAO_vs_summerHydroclimate_combined_correlations.nc"
ds_out.to_netcdf(nc_path)
print(f"✅ Combined plotting data saved to: {nc_path}")

# ==========================================================
# PLOT 2x2 FIGURE
# ==========================================================
fig, axes = plt.subplots(
    2, 2, figsize=(12, 9),
    subplot_kw={"projection": ccrs.PlateCarree()}
)
axes = axes.flatten()

plot_class_panel(
    axes[0],
    precip_summer_mean.lon, precip_summer_mean.lat, precip_class,
    f"Winter NAO vs summer precipitation (JJA, {years_precip.min()}–{years_precip.max()})",
    "(a)",
    show_left_labels=True, show_bottom_labels=False
)

plot_class_panel(
    axes[1],
    spei_aug.lon, spei_aug.lat, spei_class,
    f"Winter NAO vs August SPEI-3 ({years_spei.min()}–{years_spei.max()})",
    "(b)",
    show_left_labels=False, show_bottom_labels=False
)

plot_class_panel(
    axes[2],
    sm_summer_mean.lon, sm_summer_mean.lat, sm_class,
    f"Winter NAO vs summer surface soil moisture (JJA, {years_sm.min()}–{years_sm.max()})",
    "(c)",
    show_left_labels=True, show_bottom_labels=True
)

plot_class_panel(
    axes[3],
    grace_summer_mean.lon, grace_summer_mean.lat, grace_class,
    f"Winter NAO vs summer GRACE water storage (JJA, {years_grace.min()}–{years_grace.max()})",
    "(d)",
    show_left_labels=False, show_bottom_labels=True
)

legend_elements = [
    Patch(facecolor=cmap.colors[2], edgecolor='k', label="Significant positive (p<0.05)"),
    Patch(facecolor=cmap.colors[3], edgecolor='k', label="Non-significant positive (p≥0.05)"),
    Patch(facecolor=cmap.colors[1], edgecolor='k', label="Significant negative (p<0.05)"),
    Patch(facecolor=cmap.colors[0], edgecolor='k', label="Non-significant negative (p≥0.05)")
]

fig.legend(
    handles=legend_elements,
    ncol=2,
    loc='lower center',
    bbox_to_anchor=(0.5, 0.02),
    fontsize=10,
    frameon=False
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.12, wspace=0.08, hspace=0.10)

out_png = "figures/winterNAO_vs_summerHydroclimate_discrete_2x2_fulltime.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()

print(f"✅ Combined figure saved to: {out_png}")
