import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, linregress
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sklearn.linear_model import LinearRegression
import matplotlib.ticker as mticker
import os

# ==========================================================
# === SETTINGS
# ==========================================================
WINDOW_LEN = 21
SIG_ALPHA = 0.05

DMC_LON, DMC_LAT = 14.25, 52.3833
EU_EXTENT = [-12, 45, 35, 72]

PLOT_LAND_ONLY = True

# Option: skip recomputation and only plot from saved NetCDF outputs
PLOT_ONLY_FROM_NC = True   # True = read processed outputs and only plot figure

os.makedirs("plotting_data", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ==========================================================
# === Helper function: detrend
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
# === Vectorized correlation map
# ==========================================================
def corr_map_vectorized_valid(Y, x):
    """
    Y: np.ndarray (T, nlat, nlon)
    x: np.ndarray (T,)
    Returns r_map (nlat, nlon)
    Uses valid-mask correlation instead of requiring all years finite.
    """
    Y = np.asarray(Y, dtype=float)
    x = np.asarray(x, dtype=float).ravel()

    T, nlat, nlon = Y.shape
    r_map = np.full((nlat, nlon), np.nan, dtype=float)

    for i in range(nlat):
        for j in range(nlon):
            y = Y[:, i, j]
            valid = np.isfinite(y) & np.isfinite(x)
            if valid.sum() < 3:
                continue

            yv = y[valid]
            xv = x[valid]

            if np.std(yv) == 0 or np.std(xv) == 0:
                continue

            r_map[i, j] = np.corrcoef(yv, xv)[0, 1]

    return r_map

# ==========================================================
# === Plot helper functions
# ==========================================================
def add_common_map_features(ax, extent, dmc_lon, dmc_lat, plot_land_only=True):
    if plot_land_only:
        ocean = cfeature.NaturalEarthFeature(
            "physical", "ocean", "50m", edgecolor="face", facecolor="white"
        )
        ax.add_feature(ocean, zorder=2)

    ax.add_feature(cfeature.COASTLINE, zorder=3, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, zorder=3, linewidth=0.4)

    ax.scatter(dmc_lon, dmc_lat, marker='o', facecolors='none', edgecolors='black',
               s=25, transform=ccrs.PlateCarree(), zorder=4)
    ax.text(dmc_lon + 1, dmc_lat, 'DMC', transform=ccrs.PlateCarree(),
            fontsize=9, color='black', zorder=4)

    ax.set_extent(extent, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-10, 50, 5))
    gl.ylocator = mticker.FixedLocator(np.arange(35, 75, 5))
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    return gl

def plot_map_panel(ax, lon, lat, data2d, title, *,
                   cmap="RdBu_r", vmin=None, vmax=None,
                   sig_mask=None, cbar_label=None):
    im = ax.pcolormesh(lon, lat, data2d, cmap=cmap, vmin=vmin, vmax=vmax,
                       transform=ccrs.PlateCarree(), zorder=1)

    if sig_mask is not None:
        Lon2, Lat2 = np.meshgrid(lon.values, lat.values)
        step = 3
        ax.scatter(Lon2[::step, ::step][sig_mask[::step, ::step]],
                   Lat2[::step, ::step][sig_mask[::step, ::step]],
                   s=4, c="darkgreen", marker="o", linewidths=0,
                   transform=ccrs.PlateCarree(), zorder=2)

    add_common_map_features(ax, EU_EXTENT, DMC_LON, DMC_LAT, plot_land_only=PLOT_LAND_ONLY)

    ax.set_title(title, fontsize=10, loc='left', fontweight='bold')

    return im

def plot_three_panel_figure(trend_precip, sig_precip,
                            trend_feb, sig_feb,
                            trend_aug, sig_aug,
                            lon_precip, lat_precip,
                            lon_spei, lat_spei,
                            out_png):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8),
                             subplot_kw={'projection': ccrs.PlateCarree()})
    axes = axes.flatten()

    vmax1 = np.nanpercentile(np.abs(trend_precip), 95)
    if not np.isfinite(vmax1) or vmax1 == 0:
        vmax1 = 0.01

    vmax3 = np.nanpercentile(np.abs(trend_aug), 95)
    if not np.isfinite(vmax3) or vmax3 == 0:
        vmax3 = 0.01

    # Use one shared color scale for both panels
    vmax_common = max(vmax1, vmax3)

    im1 = plot_map_panel(
        axes[0], lon_precip, lat_precip, trend_precip,
        "(e) Correlation trend: Winter NAO vs summer precipitation",
        cmap="RdBu_r", vmin=-vmax_common, vmax=vmax_common,
        sig_mask=sig_precip
    )

    im2 = plot_map_panel(
        axes[1], lon_spei, lat_spei, trend_aug,
        "(f) Correlation trend: Winter NAO vs August SPEI-3",
        cmap="RdBu_r", vmin=-vmax_common, vmax=vmax_common,
        sig_mask=sig_aug
    )

    # One shared colorbar
    cbar = fig.colorbar(
        im2, ax=axes, orientation="horizontal", shrink=0.5)
    cbar.set_label("Trend in Fisher-z correlation (per year)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.28)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"✅ Figure saved: {out_png}")

# ==========================================================
# === Running-window trend function
# ==========================================================
def compute_running_corr_trend(data_by_year, year_coord_name, nao_series, window_len):
    all_common_years = np.intersect1d(
        data_by_year[year_coord_name].values.astype(int),
        nao_series.index.values.astype(int)
    )
    all_common_years = np.array(sorted(all_common_years), dtype=int)

    if all_common_years.size < window_len:
        raise ValueError(f"Not enough years ({all_common_years.size}) for window length {window_len}")

    data_all = data_by_year.sel({year_coord_name: all_common_years}).transpose(year_coord_name, "lat", "lon").values
    nao_all = nao_series.loc[all_common_years].values

    window_start_indices = np.arange(0, all_common_years.size - window_len + 1)
    n_windows = len(window_start_indices)

    nlat = data_by_year.sizes['lat']
    nlon = data_by_year.sizes['lon']
    r_run = np.full((n_windows, nlat, nlon), np.nan)
    window_center_years = []

    for w, i0 in enumerate(window_start_indices):
        yrs = all_common_years[i0:i0 + window_len]
        window_center_years.append(int(np.round(np.mean(yrs))))

        Yw = data_all[i0:i0 + window_len, :, :]
        xw = nao_all[i0:i0 + window_len]

        xw_detr = detrend_series(xw)

        Tt, nlat, nlon = Yw.shape
        Yd = np.full_like(Yw, np.nan, dtype=float)

        for i in range(nlat):
            for j in range(nlon):
                y = Yw[:, i, j]
                y_detr = detrend_series(y)
                Yd[:, i, j] = y_detr

        r_map = corr_map_vectorized_valid(Yd, xw_detr)
        r_run[w, :, :] = r_map

    window_center_years = np.array(window_center_years, dtype=int)

    trend_map = np.full((nlat, nlon), np.nan)
    trend_p_map = np.full((nlat, nlon), np.nan)

    for i in range(nlat):
        for j in range(nlon):
            r_ts = r_run[:, i, j]
            valid = np.isfinite(r_ts)
            if valid.sum() < 5:
                continue

            r_ts_clip = np.clip(r_ts[valid], -1 + 1e-6, 1 - 1e-6)
            z_ts = np.arctanh(r_ts_clip)
            t_ts = window_center_years[valid]

            lr = linregress(t_ts, z_ts)
            trend_map[i, j] = lr.slope
            trend_p_map[i, j] = lr.pvalue

    trend_sig_raw = np.isfinite(trend_p_map) & (trend_p_map < SIG_ALPHA)

    return trend_map, trend_p_map, trend_sig_raw, r_run, window_center_years

# ==========================================================
# === Output paths
# ==========================================================
out_nc = f"plotting_data/winterNAO_runningCorrTrend_precip_FebSPEI3_AugSPEI3_{WINDOW_LEN}yr.nc"
out_png = f"figures/winterNAO_runningCorrTrend_precip_FebSPEI3_AugSPEI3_{WINDOW_LEN}yr_e_to_f.png"

# ==========================================================
# === Branch 1: plot only from saved NetCDF
# ==========================================================
if PLOT_ONLY_FROM_NC:
    print(f"[INFO] PLOT_ONLY_FROM_NC = True")
    print(f"[INFO] Reading processed file: {out_nc}")

    if not os.path.exists(out_nc):
        raise FileNotFoundError(
            f"Saved NetCDF not found:\n{out_nc}\n"
            f"Please set PLOT_ONLY_FROM_NC = False first to compute and save it."
        )

    ds_plot = xr.open_dataset(out_nc)

    trend_precip = ds_plot["trend_precip"].values
    trend_sig_precip = ds_plot["trend_sig_precip"].values.astype(bool)

    trend_feb = ds_plot["trend_feb"].values
    trend_sig_feb = ds_plot["trend_sig_feb"].values.astype(bool)

    trend_aug = ds_plot["trend_aug"].values
    trend_sig_aug = ds_plot["trend_sig_aug"].values.astype(bool)

    lon_precip = ds_plot["lon_precip"]
    lat_precip = ds_plot["lat_precip"]

    lon_spei = ds_plot["lon_spei"]
    lat_spei = ds_plot["lat_spei"]

    plot_three_panel_figure(
        trend_precip, trend_sig_precip,
        trend_feb, trend_sig_feb,
        trend_aug, trend_sig_aug,
        lon_precip, lat_precip,
        lon_spei, lat_spei,
        out_png
    )

# ==========================================================
# === Branch 2: compute + save + plot
# ==========================================================
else:
    # ==========================================================
    # === Load NAO index data
    # ==========================================================
    nao_path = r'/data/scratch/jiangcong/P/North_Atlantic_Oscillation.csv'
    nao_data = pd.read_csv(nao_path, delimiter=",")

    nao_data["date"] = pd.to_datetime(
        nao_data["date"].astype(str) + "-" + nao_data["MONTH"].astype(str) + "-15",
        errors="coerce", format="%Y-%m-%d"
    )
    nao_data.rename(columns={'INDEX': 'NAO'}, inplace=True)
    nao_data['year'] = nao_data['date'].dt.year
    nao_data['month'] = nao_data['date'].dt.month

    # === Compute DJF (winter-year) NAO means ===
    djf_years = []
    djf_nao_values = []
    for year in range(1978, 2025):
        djf = nao_data[
            ((nao_data['year'] == year) & (nao_data['month'] == 12)) |
            ((nao_data['year'] == year + 1) & (nao_data['month'].isin([1, 2])))
        ]
        if len(djf) == 3:
            djf_years.append(year + 1)
            djf_nao_values.append(djf['NAO'].mean())

    nao_series = pd.Series(djf_nao_values, index=pd.Index(djf_years, name="year"), name='NAO_DJF')

    # ==========================================================
    # === Load precipitation data
    # ==========================================================
    ds_precip = xr.open_dataset("./Monthly/new/Monthly/europe_output_0p25deg.nc")
    precip = ds_precip['precipitation']

    def assign_summer_year(t):
        t = pd.to_datetime(t)
        return t.year if t.month in [6, 7, 8] else None

    summer_years = [assign_summer_year(t) for t in precip.time.values]
    precip = precip.assign_coords(summer_year=('time', summer_years))
    precip_summer = precip.sel(time=precip['time.month'].isin([6, 7, 8]))
    precip_summer_mean = precip_summer.groupby('summer_year').mean('time')
    precip_summer_mean = precip_summer_mean.assign_coords(
        summer_year=precip_summer_mean['summer_year'].astype(int)
    )

    # ==========================================================
    # === Load SPEI-3 data
    # ==========================================================
    ds_spei = xr.open_dataset("/data/scratch/jiangcong/ERA5/monthly/era5_spi_spei_europe/era5_europe_spei_gamma_03.nc")
    spei = ds_spei["spei_gamma_03"].sel(time=slice("1979-01", None))

    # February SPEI-3 = Dec-Jan-Feb
    spei_feb = spei.sel(time=spei["time.month"] == 2).groupby("time.year").first("time")
    spei_feb = spei_feb.transpose("year", "lat", "lon")

    # August SPEI-3 = Jun-Jul-Aug
    spei_aug = spei.sel(time=spei["time.month"] == 8).groupby("time.year").first("time")
    spei_aug = spei_aug.transpose("year", "lat", "lon")

    # ==========================================================
    # === Compute running-correlation trends
    # ==========================================================
    print("[INFO] Computing running-correlation trend: winter NAO vs summer precipitation ...")
    trend_precip, trend_p_precip, trend_sig_precip, rrun_precip, years_precip = compute_running_corr_trend(
        precip_summer_mean, "summer_year", nao_series, WINDOW_LEN
    )

    print("[INFO] Computing running-correlation trend: winter NAO vs February SPEI-3 ...")
    trend_feb, trend_p_feb, trend_sig_feb, rrun_feb, years_feb = compute_running_corr_trend(
        spei_feb, "year", nao_series, WINDOW_LEN
    )

    print("[INFO] Computing running-correlation trend: winter NAO vs August SPEI-3 ...")
    trend_aug, trend_p_aug, trend_sig_aug, rrun_aug, years_aug = compute_running_corr_trend(
        spei_aug, "year", nao_series, WINDOW_LEN
    )

    # ==========================================================
    # === Save output
    # ==========================================================
    ds_out = xr.Dataset({
        "trend_precip": xr.DataArray(
            trend_precip,
            dims=("lat_precip", "lon_precip"),
            coords={
                "lat_precip": ("lat_precip", precip_summer_mean.lat.values),
                "lon_precip": ("lon_precip", precip_summer_mean.lon.values)
            }
        ),
        "trend_p_precip": xr.DataArray(
            trend_p_precip,
            dims=("lat_precip", "lon_precip"),
            coords={
                "lat_precip": ("lat_precip", precip_summer_mean.lat.values),
                "lon_precip": ("lon_precip", precip_summer_mean.lon.values)
            }
        ),
        "trend_sig_precip": xr.DataArray(
            trend_sig_precip.astype(np.int8),
            dims=("lat_precip", "lon_precip"),
            coords={
                "lat_precip": ("lat_precip", precip_summer_mean.lat.values),
                "lon_precip": ("lon_precip", precip_summer_mean.lon.values)
            }
        ),

        "trend_feb": xr.DataArray(
            trend_feb,
            dims=("lat_spei", "lon_spei"),
            coords={
                "lat_spei": ("lat_spei", spei.lat.values),
                "lon_spei": ("lon_spei", spei.lon.values)
            }
        ),
        "trend_p_feb": xr.DataArray(
            trend_p_feb,
            dims=("lat_spei", "lon_spei"),
            coords={
                "lat_spei": ("lat_spei", spei.lat.values),
                "lon_spei": ("lon_spei", spei.lon.values)
            }
        ),
        "trend_sig_feb": xr.DataArray(
            trend_sig_feb.astype(np.int8),
            dims=("lat_spei", "lon_spei"),
            coords={
                "lat_spei": ("lat_spei", spei.lat.values),
                "lon_spei": ("lon_spei", spei.lon.values)
            }
        ),

        "trend_aug": xr.DataArray(
            trend_aug,
            dims=("lat_spei", "lon_spei"),
            coords={
                "lat_spei": ("lat_spei", spei.lat.values),
                "lon_spei": ("lon_spei", spei.lon.values)
            }
        ),
        "trend_p_aug": xr.DataArray(
            trend_p_aug,
            dims=("lat_spei", "lon_spei"),
            coords={
                "lat_spei": ("lat_spei", spei.lat.values),
                "lon_spei": ("lon_spei", spei.lon.values)
            }
        ),
        "trend_sig_aug": xr.DataArray(
            trend_sig_aug.astype(np.int8),
            dims=("lat_spei", "lon_spei"),
            coords={
                "lat_spei": ("lat_spei", spei.lat.values),
                "lon_spei": ("lon_spei", spei.lon.values)
            }
        ),
    })

    ds_out.attrs["title"] = f"Trend in {WINDOW_LEN}-year running correlation with winter NAO"
    ds_out.attrs["note"] = "Panels include summer precipitation, February SPEI-3, and August SPEI-3; trend tested on Fisher-z transformed running correlations"
    ds_out.to_netcdf(out_nc)
    print(f"✅ NetCDF saved: {out_nc}")

    # ==========================================================
    # === Plot final figure
    # ==========================================================
    plot_three_panel_figure(
        trend_precip, trend_sig_precip,
        trend_feb, trend_sig_feb,
        trend_aug, trend_sig_aug,
        precip_summer_mean.lon, precip_summer_mean.lat,
        spei.lon, spei.lat,
        out_png
    )
