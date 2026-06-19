"""
winterNAO_vs_ERA5_VIWV_ANOMALIES_DECtoNOV_POSNEG.py

Purpose:
1) Plot Arctic monthly climatology of ERA5 vertically integrated water vapour flux (VIWV) from DEC..NOV
   - shading: VIWV magnitude climatology
   - arrows : viwve / viwvn vectors
2) Plot POS and NEG winter NAO composite anomalies relative to climatology
   - shading: VIWV magnitude anomaly
   - arrows : viwve / viwvn anomalies
3) Use North Polar Stereo projection for better Arctic visualization
"""

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import matplotlib.ticker as mticker

# ==========================================================
# === SETTINGS ===
# ==========================================================
VIWV_FILE = "/data/scratch/jiangcong/ERA5/monthly/download/viwv_1979-2025_renamed_squeezed.nc"

VAR_UNITS = {
    "viwv": r"kg m$^{-1}$ s$^{-1}$"
}

MONTHS = [12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
MONTH_ABBR = {
    1: "JAN", 2: "FEB", 3: "MAR", 4: "APR", 5: "MAY", 6: "JUN",
    7: "JUL", 8: "AUG", 9: "SEP", 10: "OCT", 11: "NOV", 12: "DEC"
}

# DJF-years (Dec belongs to next year)
POS_DJF_YEARS = np.array([1989, 1995, 2000, 2012, 2015, 2016, 2018, 2020], dtype=int)
NEG_DJF_YEARS = np.array([1979, 1985, 1987, 1996, 2010, 2011, 2021], dtype=int)

N_COLS = 3
ARCTIC_LAT_MIN = 40.0

# VIWV quiver settings (from code B)
QUIVER_STEP = 18
QUIVER_WIDTH = 0.003

VIWV_QUIVER_SCALE_CLIM = 3000
VIWV_QUIVER_SCALE_ANOM = 600

VIWV_QUIVER_KEY_VALUE_CLIM = 200
VIWV_QUIVER_KEY_LABEL_CLIM = r"200 kg m$^{-1}$ s$^{-1}$"

VIWV_QUIVER_KEY_VALUE_ANOM = 50
VIWV_QUIVER_KEY_LABEL_ANOM = r"50 kg m$^{-1}$ s$^{-1}$"

os.makedirs("plotting_data", exist_ok=True)
os.makedirs("figures", exist_ok=True)

dmc_lon, dmc_lat = 14.25, 52.3833

# ==========================================================
# === Helper: standardize coordinate names/dims
# ==========================================================
def standardize_coords(ds):
    ds = ds.copy()

    rename_map = {}
    if "valid_time" in ds.coords or "valid_time" in ds.dims:
        rename_map["valid_time"] = "time"
    if "latitude" in ds.coords or "latitude" in ds.dims:
        rename_map["latitude"] = "lat"
    if "longitude" in ds.coords or "longitude" in ds.dims:
        rename_map["longitude"] = "lon"
    if rename_map:
        ds = ds.rename(rename_map)

    for dim in ["number"]:
        if dim in ds.dims and ds.sizes[dim] == 1:
            ds = ds.isel({dim: 0}, drop=True)

    if "lon" in ds.coords:
        lon = ds["lon"]
        if float(lon.max()) > 180:
            ds = ds.assign_coords(lon=(((lon + 180) % 360) - 180)).sortby("lon")

    return ds

# ==========================================================
# === Helper: robust symmetric vmax for anomaly plots
# ==========================================================
def get_robust_symmetric_vmax(data_dict, percentile=98):
    vals = []
    for m in data_dict:
        arr = data_dict[m].values
        arr = arr[np.isfinite(arr)]
        if arr.size > 0:
            vals.append(np.abs(arr).ravel())
    if len(vals) == 0:
        return 1.0
    vals = np.concatenate(vals)
    vmax = np.nanpercentile(vals, percentile)
    if (not np.isfinite(vmax)) or (vmax == 0):
        vmax = 1.0
    return float(vmax)

# ==========================================================
# === Helper: make circular boundary for polar plots
# ==========================================================
def set_circular_boundary(ax):
    theta = np.linspace(0, 2*np.pi, 200)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = plt.matplotlib.path.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

def add_barents_kara_highlight(ax, lon_center=60, lat_center=75,
                               lon_radius=35, lat_radius=8,
                               color='red', lw=1.8):
    t = np.linspace(0, 2 * np.pi, 300)
    lons = lon_center + lon_radius * np.cos(t)
    lats = lat_center + lat_radius * np.sin(t)

    ax.plot(
        lons, lats,
        color=color,
        linewidth=lw,
        transform=ccrs.PlateCarree(),
        zorder=10
    )

# ==========================================================
# === Plot helper: composite anomaly figure (POS or NEG)
# ==========================================================
def plot_composite(tag, anoms_by_month, uv_anoms_by_month, vmax, out_png):
    n = len(MONTHS)
    ncols = int(N_COLS)
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(3.8 * ncols, 3.8 * nrows))

    im = None

    for idx, m in enumerate(MONTHS):
        mon = MONTH_ABBR[m]
        ax = plt.subplot(nrows, ncols, idx + 1, projection=ccrs.NorthPolarStereo())

        da_plot = anoms_by_month[m]
        lon, lat = da_plot.lon, da_plot.lat

        im = ax.pcolormesh(
            lon, lat, da_plot,
            cmap="RdBu_r", vmin=-vmax, vmax=vmax,
            transform=ccrs.PlateCarree()
        )

        # VIWV anomaly vectors
        u_anom, v_anom = uv_anoms_by_month[m]
        u0 = u_anom[::QUIVER_STEP, ::QUIVER_STEP]
        v0 = v_anom[::QUIVER_STEP, ::QUIVER_STEP]

        qu = ax.quiver(
            u0.lon, u0.lat,
            u0.values, v0.values,
            transform=ccrs.PlateCarree(),
            scale=VIWV_QUIVER_SCALE_ANOM,
            width=QUIVER_WIDTH
        )
        ax.quiverkey(
            qu, 0.85, -0.15,
            VIWV_QUIVER_KEY_VALUE_ANOM, VIWV_QUIVER_KEY_LABEL_ANOM,
            transform=ax.transAxes
        )

        ax.scatter(
            dmc_lon, dmc_lat,
            marker='o', facecolors='none', edgecolors='red',
            s=20, transform=ccrs.PlateCarree(), zorder=4
        )
        ax.text(dmc_lon + 1, dmc_lat, 'DMC', transform=ccrs.PlateCarree(),
            fontsize=9, color='black', zorder=4)
        
       # ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='none', zorder=2)
        ax.coastlines(linewidth=0.6, zorder=3)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, zorder=2)
        ax.set_extent([-180, 180, ARCTIC_LAT_MIN, 90], crs=ccrs.PlateCarree())
        add_barents_kara_highlight(ax)

        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=0.5,
            color='gray',
            zorder=5,
            alpha=0.8,
            linestyle='--'
        )
        gl.xlabel_style = {'size': 8, 'rotation': 0}
        gl.ylabel_style = {'size': 8, 'rotation': 0}
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
        gl.ylocator = mticker.FixedLocator([50, 60, 70, 80])

        set_circular_boundary(ax)

        label = f"({chr(97 + idx)}) {mon}"
        ax.text(0.1, 0.98, label, transform=ax.transAxes,
                ha='center', va='top', fontsize=11, fontweight='bold')

    plt.tight_layout()
    fig.subplots_adjust(wspace=0.0, hspace=0.1, bottom=0.10)

    cax = fig.add_axes([0.20, 0.05, 0.60, 0.02])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")

    if tag == "POS":
        nao_label = "NAO+"
    elif tag == "NEG":
        nao_label = "NAO−"
    else:
        nao_label = "NAO"

    unit = VAR_UNITS.get("viwv", "")
    cbar.set_label(f"VIWV anomaly (winter {nao_label} composite, {unit})", fontsize=15)
    cbar.ax.tick_params(labelsize=12)

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show(block=False)

# ==========================================================
# === Plot helper: monthly climatology figure
# ==========================================================
def plot_climatology(clims_by_month, uv_clim_by_month, out_png):
    n = len(MONTHS)
    ncols = int(N_COLS)
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(4.2 * ncols, 3.8 * nrows))

    im = None

    vals = []
    for m in MONTHS:
        vals.append(np.nanmin(clims_by_month[m].values))
        vals.append(np.nanmax(clims_by_month[m].values))
    vmin = float(np.nanmin(vals))
    vmax = float(np.nanmax(vals))

    if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or (vmin == vmax):
        vmin, vmax = 0.0, 1.0

    for idx, m in enumerate(MONTHS):
        mon = MONTH_ABBR[m]
        ax = plt.subplot(nrows, ncols, idx + 1, projection=ccrs.NorthPolarStereo())

        da_plot = clims_by_month[m]
        lon, lat = da_plot.lon, da_plot.lat

        im = ax.pcolormesh(
            lon, lat, da_plot,
            vmin=vmin, vmax=vmax,
            transform=ccrs.PlateCarree()
        )

        # VIWV climatology vectors
        u_clim, v_clim = uv_clim_by_month[m]
        u0 = u_clim[::QUIVER_STEP, ::QUIVER_STEP]
        v0 = v_clim[::QUIVER_STEP, ::QUIVER_STEP]

        qu = ax.quiver(
            u0.lon, u0.lat,
            u0.values, v0.values,
            transform=ccrs.PlateCarree(),
            scale=VIWV_QUIVER_SCALE_CLIM,
            width=QUIVER_WIDTH
        )
        ax.quiverkey(
            qu, 0.85, -0.15,
            VIWV_QUIVER_KEY_VALUE_CLIM, VIWV_QUIVER_KEY_LABEL_CLIM,
            transform=ax.transAxes
        )

        ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='none', zorder=2)
        ax.coastlines(linewidth=0.6, zorder=3)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, zorder=3)

        ax.set_extent([-180, 180, ARCTIC_LAT_MIN, 90], crs=ccrs.PlateCarree())

        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=0.5,
            color='gray',
            alpha=0.4,
            zorder=5,
            linestyle='--'
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))

        lat_labels = np.arange(int(np.ceil(ARCTIC_LAT_MIN / 10) * 10), 90, 10)
        gl.ylocator = mticker.FixedLocator(lat_labels)

        set_circular_boundary(ax)
        ax.set_title(f"CLIM {mon}", fontsize=11, loc='left', fontweight='bold')

    plt.tight_layout()
    fig.subplots_adjust(bottom=0.09)

    cax = fig.add_axes([0.20, 0.06, 0.60, 0.02])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    unit = VAR_UNITS.get("viwv", "")
    cbar.set_label(f"VIWV climatology monthly mean ({unit})")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show(block=False)

# ==========================================================
# === MAIN
# ==========================================================
if not os.path.exists(VIWV_FILE):
    raise FileNotFoundError(f"VIWV file not found: {VIWV_FILE}")

print(f"\n=== Processing VIWV: {VIWV_FILE} ===")

ds = xr.open_dataset(VIWV_FILE)
ds = standardize_coords(ds)

if ("viwve" not in ds.variables) or ("viwvn" not in ds.variables):
    raise KeyError("Variables 'viwve' and 'viwvn' not found in viwv file")

da_u = ds["viwve"]
da_v = ds["viwvn"]
da = np.sqrt(da_u**2 + da_v**2)
da.name = "viwv"

# Ensure time decoded
if not np.issubdtype(da["time"].dtype, np.datetime64):
    da["time"] = xr.decode_cf(ds).time
    da_u["time"] = da["time"]
    da_v["time"] = da["time"]

# ensure latitude is ascending before subsetting
if da["lat"][0] > da["lat"][-1]:
    da = da.sortby("lat")
    da_u = da_u.sortby("lat")
    da_v = da_v.sortby("lat")

# subset Arctic
da = da.sel(lat=slice(ARCTIC_LAT_MIN, 90))
da_u = da_u.sel(lat=slice(ARCTIC_LAT_MIN, 90))
da_v = da_v.sel(lat=slice(ARCTIC_LAT_MIN, 90))

print("Arctic subset shape:", da.shape)
print("Latitude range:", float(da.lat.min()), "to", float(da.lat.max()))

# Create year/month coords
t = pd.to_datetime(da["time"].values)
year = t.year
month = t.month
djf_year = year + (month == 12).astype(int)

da = da.assign_coords(
    year=xr.DataArray(year, dims="time", coords={"time": da["time"]}),
    month=xr.DataArray(month, dims="time", coords={"time": da["time"]}),
    djf_year=xr.DataArray(djf_year, dims="time", coords={"time": da["time"]}),
)
da_u = da_u.assign_coords(year=da["year"], month=da["month"], djf_year=da["djf_year"])
da_v = da_v.assign_coords(year=da["year"], month=da["month"], djf_year=da["djf_year"])

# restrict to analysis period
da = da.sel(time=slice("1979-01-01", "2025-12-31"))
da_u = da_u.sel(time=slice("1979-01-01", "2025-12-31"))
da_v = da_v.sel(time=slice("1979-01-01", "2025-12-31"))

# ----------------------------------------------------------
# Compute climatology and POS/NEG anomalies
# ----------------------------------------------------------
pos_anoms, neg_anoms = {}, {}
clims = {}

pos_uv_anoms, neg_uv_anoms = {}, {}
clim_uv = {}

for m in MONTHS:
    mon = MONTH_ABBR[m]

    da_m = da.where(da["month"] == m, drop=True)
    if da_m.sizes.get("time", 0) < 3:
        raise ValueError(f"[viwv {mon}] Not enough timesteps for this month.")

    da_m_mean = da_m.groupby("djf_year").mean("time", skipna=True)
    y_all = da_m_mean["djf_year"].values.astype(int)

    clim = da_m_mean.mean("djf_year", skipna=True)
    clims[m] = clim

    pos_years = np.intersect1d(y_all, POS_DJF_YEARS)
    neg_years = np.intersect1d(y_all, NEG_DJF_YEARS)

    if pos_years.size < 1:
        raise ValueError(f"[viwv {mon}] No overlapping POS DJF-years found in data.")
    if neg_years.size < 1:
        raise ValueError(f"[viwv {mon}] No overlapping NEG DJF-years found in data.")

    pos_comp = da_m_mean.sel(djf_year=pos_years).mean("djf_year", skipna=True)
    neg_comp = da_m_mean.sel(djf_year=neg_years).mean("djf_year", skipna=True)

    pos_anom = pos_comp - clim
    neg_anom = neg_comp - clim

    pos_anoms[m] = pos_anom
    neg_anoms[m] = neg_anom

    # u/v anomalies for VIWV
    da_u_m = da_u.where(da_u["month"] == m, drop=True)
    da_v_m = da_v.where(da_v["month"] == m, drop=True)

    u_m_mean = da_u_m.groupby("djf_year").mean("time", skipna=True)
    v_m_mean = da_v_m.groupby("djf_year").mean("time", skipna=True)

    u_clim = u_m_mean.mean("djf_year", skipna=True)
    v_clim = v_m_mean.mean("djf_year", skipna=True)
    clim_uv[m] = (u_clim, v_clim)

    u_pos_comp = u_m_mean.sel(djf_year=pos_years).mean("djf_year", skipna=True)
    v_pos_comp = v_m_mean.sel(djf_year=pos_years).mean("djf_year", skipna=True)

    u_neg_comp = u_m_mean.sel(djf_year=neg_years).mean("djf_year", skipna=True)
    v_neg_comp = v_m_mean.sel(djf_year=neg_years).mean("djf_year", skipna=True)

    pos_uv_anoms[m] = (u_pos_comp - u_clim, v_pos_comp - v_clim)
    neg_uv_anoms[m] = (u_neg_comp - u_clim, v_neg_comp - v_clim)

    # save monthly data
    nc_out_path = f"plotting_data/winterNAO_{mon.lower()}_viwv_anomaly_POSNEG.nc"
    ds_out = xr.Dataset({
        "pos_anomaly": pos_anom,
        "neg_anomaly": neg_anom,
        "climatology": clim,
        "pos_composite": pos_comp,
        "neg_composite": neg_comp,
        "pos_u_anom": pos_uv_anoms[m][0],
        "pos_v_anom": pos_uv_anoms[m][1],
        "neg_u_anom": neg_uv_anoms[m][0],
        "neg_v_anom": neg_uv_anoms[m][1],
        "u_climatology": clim_uv[m][0],
        "v_climatology": clim_uv[m][1],
    })
    ds_out.attrs["title"] = f"VIWV composite monthly anomalies (POS/NEG DJF-year NAO) - {mon}"
    ds_out.attrs["viwv_file"] = VIWV_FILE
    ds_out.attrs["pos_djf_years"] = ",".join(map(str, POS_DJF_YEARS.tolist()))
    ds_out.attrs["neg_djf_years"] = ",".join(map(str, NEG_DJF_YEARS.tolist()))
    ds_out.to_netcdf(nc_out_path)

# ----------------------------------------------------------
# Shared robust anomaly range for POS and NEG
# ----------------------------------------------------------
all_anoms = {}
all_anoms.update(pos_anoms)
all_anoms.update(neg_anoms)
anom_vmax = get_robust_symmetric_vmax(all_anoms, percentile=98)

if not np.isfinite(anom_vmax) or anom_vmax == 0:
    anom_vmax = 1.0

# POS figure
out_png_pos = "figures/winterNAO_composite_anomaly_viwv_POS_DECtoNOV_Arctic.png"
plot_composite("POS", pos_anoms, pos_uv_anoms, anom_vmax, out_png_pos)
print(f"✅ POS figure saved: {out_png_pos}")

# NEG figure
out_png_neg = "figures/winterNAO_composite_anomaly_viwv_NEG_DECtoNOV_Arctic.png"
plot_composite("NEG", neg_anoms, neg_uv_anoms, anom_vmax, out_png_neg)
print(f"✅ NEG figure saved: {out_png_neg}")

# CLIM figure
out_png_clim = "figures/monthly_climatology_viwv_DECtoNOV_Arctic.png"
plot_climatology(clims, clim_uv, out_png_clim)
print(f"✅ CLIM figure saved: {out_png_clim}")
