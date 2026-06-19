"""
winterNAO_vs_HadISST_ArcticSeaIce_ANOMALIES_DECtoNOV_POSNEG.py

Purpose:
1) Plot Arctic sea-ice concentration (sic) monthly climatology (DEC..NOV)
2) Plot POS and NEG winter NAO composite anomalies relative to climatology
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
SEAICE_FILE = "/data/scratch/jiangcong/ERA5/monthly/sea_ice_cover/HadISST_ice.nc"   # <- change if needed

MONTHS = [12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
MONTH_ABBR = {
    1: "JAN", 2: "FEB", 3: "MAR", 4: "APR", 5: "MAY", 6: "JUN",
    7: "JUL", 8: "AUG", 9: "SEP", 10: "OCT", 11: "NOV", 12: "DEC"
}

# DJF-years (Dec belongs to next year)
POS_DJF_YEARS = np.array([1989, 1995, 2000, 2012, 2015, 2016, 2018, 2020], dtype=int)
NEG_DJF_YEARS = np.array([1979, 1985, 1987, 1996, 2010, 2011, 2021], dtype=int)

N_COLS = 3
ARCTIC_LAT_MIN = 49.0   # plot 50N to pole

os.makedirs("plotting_data", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ==========================================================
# === Helper: standardize coordinate names/dims
# ==========================================================
def standardize_coords(ds):
    ds = ds.copy()

    rename_map = {}
    if "latitude" in ds.coords or "latitude" in ds.dims:
        rename_map["latitude"] = "lat"
    if "longitude" in ds.coords or "longitude" in ds.dims:
        rename_map["longitude"] = "lon"
    if rename_map:
        ds = ds.rename(rename_map)

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
        return 0.1
    vals = np.concatenate(vals)
    vmax = np.nanpercentile(vals, percentile)
    if (not np.isfinite(vmax)) or (vmax == 0):
        vmax = 0.1
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
def plot_composite(tag, anoms_by_month, vmax, out_png):
    n = len(MONTHS)
    ncols = int(N_COLS)
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(3.8* ncols, 3.8 * nrows))

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

        # optional sea-ice edge contour from anomaly = 0 is not very useful, skip
        #ax.coastlines(linewidth=0.6)
        #ax.add_feature(cfeature.BORDERS, linewidth=0.3)

        ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='none', zorder=2)
        ax.coastlines(linewidth=0.6, zorder=3)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, zorder=2)
        ax.set_extent([-180, 180, ARCTIC_LAT_MIN, 90], crs=ccrs.PlateCarree())
        #ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.4, linestyle='--')
        add_barents_kara_highlight(ax)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', zorder=5, alpha=0.8, linestyle='--')

        gl.xlabel_style = {'size': 8, 'rotation': 0}
        gl.ylabel_style = {'size': 8, 'rotation': 0}

        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
        gl.ylocator = mticker.FixedLocator([75, 85])
        
        set_circular_boundary(ax)
        #ax.set_title(f"({chr(97 + idx)}) {mon}", fontsize=11, loc='left', fontweight='bold')

        label = f"({chr(98 + idx)}) {mon}"

        ax.text(0.5, 0.98, label, transform=ax.transAxes, ha='center', va='top', fontsize=11, fontweight='bold')

        

    plt.tight_layout()
    # manual spacing control: better than mixing with tight_layout()
    fig.subplots_adjust(wspace=0.0, hspace=0.1, bottom=0.09)

    cax = fig.add_axes([0.20, 0.05, 0.60, 0.02])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")


    if tag == "POS":
        nao_label = "NAO+"
    elif tag == "NEG":
        nao_label = "NAO−"
    else:
        nao_label = "NAO"

    cbar.set_label(f"Sea-ice concentration anomaly (winter {nao_label} composite)", fontsize=15)

    cbar.ax.tick_params(labelsize=12)

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show(block=False)

# ==========================================================
# === Plot helper: monthly climatology figure
# ==========================================================
def plot_climatology(clims_by_month, out_png):
    n = len(MONTHS)
    ncols = int(N_COLS)
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(4.2 * ncols, 3.8 * nrows))

    im = None
    vmin, vmax = 0.0, 1.0   # sea-ice concentration fraction

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

        # add 15% sea-ice contour as common sea-ice edge
        cs = ax.contour(
            lon, lat, da_plot,
            levels=[0.15],
            colors='k',
            linewidths=0.7,
            transform=ccrs.PlateCarree()
        )

        #ax.coastlines(linewidth=0.6)
        #ax.add_feature(cfeature.BORDERS, linewidth=0.3)

        ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='none', zorder=2)
        ax.coastlines(linewidth=0.6, zorder=3)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, zorder=3)

        ax.set_extent([-180, 180, ARCTIC_LAT_MIN, 90], crs=ccrs.PlateCarree())
        #ax.gridlines(draw_labels=False, linewidth=0.3, color='gray', alpha=0.4, linestyle='--')

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.4, zorder=5, linestyle='--')
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
    cbar.set_label("Sea-ice concentration climatology")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show(block=False)

# ==========================================================
# === MAIN
# ==========================================================
if not os.path.exists(SEAICE_FILE):
    raise FileNotFoundError(f"Sea-ice file not found: {SEAICE_FILE}")

print(f"\n=== Processing sea ice: {SEAICE_FILE} ===")

ds = xr.open_dataset(SEAICE_FILE)
ds = standardize_coords(ds)

if "sic" not in ds.variables:
    raise KeyError("Variable 'sic' not found in HadISST_ice.nc")

da = ds["sic"]
da.name = "sic"

# decode time if needed
if not np.issubdtype(da["time"].dtype, np.datetime64):
    ds = xr.decode_cf(ds)
    da = ds["sic"]

# replace fill values with NaN
#da = da.where(da > -1e20)

# subset Arctic
#da = da.sel(lat=slice(ARCTIC_LAT_MIN, 90))



# replace fill values with NaN
da = da.where(da > -1e20)

# ensure latitude is ascending before subsetting
if da["lat"][0] > da["lat"][-1]:
    da = da.sortby("lat")

# subset Arctic
da = da.sel(lat=slice(ARCTIC_LAT_MIN, 90))

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


# restrict to analysis period
da = da.sel(time=slice("1979-01-01", "2025-12-31"))

# ----------------------------------------------------------
# Compute climatology and POS/NEG anomalies
# ----------------------------------------------------------
pos_anoms, neg_anoms = {}, {}
clims = {}

for m in MONTHS:
    mon = MONTH_ABBR[m]

    da_m = da.where(da["month"] == m, drop=True)
    if da_m.sizes.get("time", 0) < 3:
        raise ValueError(f"[sic {mon}] Not enough timesteps for this month.")

    da_m_mean = da_m.groupby("djf_year").mean("time", skipna=True)
    y_all = da_m_mean["djf_year"].values.astype(int)

    clim = da_m_mean.mean("djf_year", skipna=True)
    clims[m] = clim

    pos_years = np.intersect1d(y_all, POS_DJF_YEARS)
    neg_years = np.intersect1d(y_all, NEG_DJF_YEARS)

    if pos_years.size < 1:
        raise ValueError(f"[sic {mon}] No overlapping POS DJF-years found in data.")
    if neg_years.size < 1:
        raise ValueError(f"[sic {mon}] No overlapping NEG DJF-years found in data.")

    pos_comp = da_m_mean.sel(djf_year=pos_years).mean("djf_year", skipna=True)
    neg_comp = da_m_mean.sel(djf_year=neg_years).mean("djf_year", skipna=True)

    pos_anom = pos_comp - clim
    neg_anom = neg_comp - clim

    pos_anoms[m] = pos_anom
    neg_anoms[m] = neg_anom

    # save monthly data
    nc_out_path = f"plotting_data/winterNAO_{mon.lower()}_seaice_anomaly_POSNEG.nc"
    ds_out = xr.Dataset({
        "pos_anomaly": pos_anom,
        "neg_anomaly": neg_anom,
        "climatology": clim,
        "pos_composite": pos_comp,
        "neg_composite": neg_comp,
    })
    ds_out.attrs["title"] = f"Sea-ice concentration composite monthly anomalies (POS/NEG DJF-year NAO) - {mon}"
    ds_out.attrs["seaice_file"] = SEAICE_FILE
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
    anom_vmax = 0.1

# POS figure
out_png_pos = "figures/winterNAO_composite_anomaly_seaice_POS_DECtoNOV_Arctic.png"
plot_composite("POS", pos_anoms, anom_vmax, out_png_pos)
print(f"✅ POS figure saved: {out_png_pos}")

# NEG figure
out_png_neg = "figures/winterNAO_composite_anomaly_seaice_NEG_DECtoNOV_Arctic.png"
plot_composite("NEG", neg_anoms, anom_vmax, out_png_neg)
print(f"✅ NEG figure saved: {out_png_neg}")

# CLIM figure
out_png_clim = "figures/monthly_climatology_seaice_DECtoNOV_Arctic.png"
plot_climatology(clims, out_png_clim)
print(f"✅ CLIM figure saved: {out_png_clim}")
