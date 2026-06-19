import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

# ==========================================================
# SETTINGS
# ==========================================================
SIC_FILE = "/data/scratch/jiangcong/ERA5/monthly/sea_ice_cover/HadISST_ice.nc"
NAO_FILE = "/data/scratch/jiangcong/P/North_Atlantic_Oscillation.csv"

OUTDIR = "./figures_nao_sic"
os.makedirs(OUTDIR, exist_ok=True)

MONTHS = [12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
MONTH_ABBR = {
    1: "JAN", 2: "FEB", 3: "MAR", 4: "APR", 5: "MAY", 6: "JUN",
    7: "JUL", 8: "AUG", 9: "SEP", 10: "OCT", 11: "NOV", 12: "DEC"
}

ARCTIC_LAT_MIN = 70.0


# ==========================================================
# HELPERS
# ==========================================================
def standardize_time_lat_lon(ds):
    ds = ds.copy()
    rename_map = {}

    if "latitude" in ds.coords or "latitude" in ds.dims:
        rename_map["latitude"] = "lat"
    if "longitude" in ds.coords or "longitude" in ds.dims:
        rename_map["longitude"] = "lon"
    if "valid_time" in ds.coords or "valid_time" in ds.dims:
        rename_map["valid_time"] = "time"

    if rename_map:
        ds = ds.rename(rename_map)

    if "lon" in ds.coords:
        lon = ds["lon"]
        if float(lon.max()) > 180:
            ds = ds.assign_coords(lon=((lon + 180) % 360) - 180).sortby("lon")

    return ds


def ensure_datetime_time(ds):
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        ds = xr.decode_cf(ds)
    return ds


def detrend_1d(y):
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)
    mask = np.isfinite(y)

    if mask.sum() < 3:
        return out

    x = np.arange(len(y), dtype=float)
    slope, intercept, _, _, _ = stats.linregress(x[mask], y[mask])
    out[mask] = y[mask] - (slope * x[mask] + intercept)
    return out


def pearsonr_valid(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan, np.nan, mask.sum()

    xx = x[mask]
    yy = y[mask]

    if np.allclose(np.nanstd(xx), 0.0) or np.allclose(np.nanstd(yy), 0.0):
        return np.nan, np.nan, mask.sum()

    r, p = stats.pearsonr(xx, yy)
    return r, p, mask.sum()


def assign_year_month_djfyear_to_da(da):
    t = pd.to_datetime(da["time"].values)
    year = t.year
    month = t.month
    djf_year = year + (month == 12).astype(int)

    da = da.assign_coords(
        year=("time", year),
        month=("time", month),
        djf_year=("time", djf_year),
    )
    return da


def monthly_series_by_djfyear(index_da, month):
    da_m = index_da.where(index_da["month"] == month, drop=True)
    da_y = da_m.groupby("djf_year").mean("time", skipna=True)
    return da_y


def seasonal_mean_by_djfyear(index_da, months):
    da_s = index_da.where(index_da["month"].isin(months), drop=True)
    da_y = da_s.groupby("djf_year").mean("time", skipna=True)
    return da_y


def correlate_with_nao(nao_djf, target_dict):
    rows = []
    for m in MONTHS:
        y = target_dict[m]

        common_years = np.intersect1d(
            nao_djf.index.values.astype(int),
            y["djf_year"].values.astype(int)
        )

        if common_years.size < 10:
            rows.append({
                "month": m,
                "month_name": MONTH_ABBR[m],
                "n": common_years.size,
                "r": np.nan,
                "p": np.nan
            })
            continue

        x = nao_djf.loc[common_years].values.astype(float)
        yy = y.sel(djf_year=common_years).values.astype(float)

        x_dt = detrend_1d(x)
        y_dt = detrend_1d(yy)

        r, p, n = pearsonr_valid(x_dt, y_dt)

        rows.append({
            "month": m,
            "month_name": MONTH_ABBR[m],
            "n": n,
            "r": r,
            "p": p
        })

    return pd.DataFrame(rows)


def corr_season(nao_djf, target_da, season_name):
    common_years = np.intersect1d(
        nao_djf.index.values.astype(int),
        target_da["djf_year"].values.astype(int)
    )

    x = nao_djf.loc[common_years].values.astype(float)
    y = target_da.sel(djf_year=common_years).values.astype(float)

    x_dt = detrend_1d(x)
    y_dt = detrend_1d(y)

    r, p, n = pearsonr_valid(x_dt, y_dt)

    return {"season": season_name, "n": n, "r": r, "p": p}


def plot_monthly_corr(df, title, out_png):
    x = np.arange(len(df))
    r = df["r"].values
    p = df["p"].values
    labels = df["month_name"].values

    plt.figure(figsize=(8.5*1.2, 2.5*1.2))
    plt.axhline(0, color="k", lw=0.8)
    plt.plot(x, r, marker="o", color="steelblue", label="Correlation (r)")

    sig = np.isfinite(p) & (p < 0.05)
    plt.scatter(x[sig], r[sig], s=70, marker="*", color="red",label="p < 0.05", zorder=3)

    plt.xticks(x, labels)
    plt.ylabel("Pearson r")
    plt.ylim(-0.6, 0)   # as you requested
    plt.title(title)
    plt.legend(loc="lower right")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()


# ==========================================================
# 1. LOAD NAO INDEX DATA
# ==========================================================
nao_data = pd.read_csv(
    NAO_FILE,
    delimiter=",",
    parse_dates=[["date", "MONTH"]]
)

nao_data.rename(columns={"date_MONTH": "date", "INDEX": "NAO"}, inplace=True)
nao_data["year"] = nao_data["date"].dt.year
nao_data["month"] = nao_data["date"].dt.month

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
        djf_years.append(year + 1)
        djf_nao_values.append(djf["NAO"].mean())

nao_series = pd.Series(djf_nao_values, index=djf_years, name="NAO_DJF")
print("\n=== Winter DJF NAO series ===")
print(nao_series.head())
print(nao_series.tail())


# ==========================================================
# 2. READ HADISST SIC
# ==========================================================
ds_sic = xr.open_dataset(SIC_FILE)
ds_sic = standardize_time_lat_lon(ds_sic)
ds_sic = ensure_datetime_time(ds_sic)

sic = ds_sic["sic"].where(ds_sic["sic"] > -1e20)
sic = assign_year_month_djfyear_to_da(sic)

# Robust Arctic selection for ascending/descending latitude
lat_vals = sic["lat"].values
if lat_vals[0] < lat_vals[-1]:
    sic_arctic = sic.sel(lat=slice(ARCTIC_LAT_MIN, 90))
else:
    sic_arctic = sic.sel(lat=slice(90, ARCTIC_LAT_MIN))

print("\n=== SIC diagnostics ===")
print("Original lat first/last:", float(sic["lat"].values[0]), float(sic["lat"].values[-1]))
print("Arctic sic shape:", sic_arctic.shape)
print("Arctic lat min/max:", float(sic_arctic["lat"].min()), float(sic_arctic["lat"].max()))
print("Finite Arctic SIC count:", int(np.isfinite(sic_arctic.values).sum()))
print("Arctic SIC min/max:", float(np.nanmin(sic_arctic.values)), float(np.nanmax(sic_arctic.values)))

# latitude weights
weights_lat = np.cos(np.deg2rad(sic_arctic["lat"]))
weights_lat.name = "weights_lat"

# ==========================================================
# 3. BUILD ARCTIC MEAN SIC INDEX
# ==========================================================
mean_sic_index = sic_arctic.weighted(weights_lat).mean(dim=("lat", "lon"))
mean_sic_index.name = "Arctic_mean_SIC_70_90N"
mean_sic_index = assign_year_month_djfyear_to_da(mean_sic_index)

print("\n=== Index diagnostics ===")
print("mean_sic_index first 5:")
print(mean_sic_index.to_series().head())
print("mean_sic_index std:", float(mean_sic_index.to_series().std()))


# ==========================================================
# 4. MONTHLY SERIES DEC..NOV BY DJF-YEAR
# ==========================================================
monthly_mean_sic = {m: monthly_series_by_djfyear(mean_sic_index, m) for m in MONTHS}

mean_sic_jja = seasonal_mean_by_djfyear(mean_sic_index, [6, 7, 8])
mean_sic_jas = seasonal_mean_by_djfyear(mean_sic_index, [7, 8, 9])
mean_sic_son = seasonal_mean_by_djfyear(mean_sic_index, [9, 10, 11])

# ==========================================================
# 5. RESTRICT TO COMMON PERIOD
# ==========================================================
all_target_years = np.unique(np.concatenate([
    monthly_mean_sic[m]["djf_year"].values.astype(int) for m in MONTHS
]))
common_all = np.intersect1d(nao_series.index.values.astype(int), all_target_years)

nao_series = nao_series.loc[common_all]

for m in MONTHS:
    monthly_mean_sic[m] = monthly_mean_sic[m].sel(
        djf_year=np.intersect1d(monthly_mean_sic[m]["djf_year"], common_all)
    )

mean_sic_jja = mean_sic_jja.sel(djf_year=np.intersect1d(mean_sic_jja["djf_year"], common_all))
mean_sic_jas = mean_sic_jas.sel(djf_year=np.intersect1d(mean_sic_jas["djf_year"], common_all))
mean_sic_son = mean_sic_son.sel(djf_year=np.intersect1d(mean_sic_son["djf_year"], common_all))

print(f"\nCommon DJF years used: {common_all.min()} - {common_all.max()} (n={len(common_all)})")


# ==========================================================
# 6. CORRELATIONS
# ==========================================================
df_monthly_mean_sic = correlate_with_nao(nao_series, monthly_mean_sic)

season_results_mean_sic = pd.DataFrame([
    corr_season(nao_series, mean_sic_jja, "JJA"),
    corr_season(nao_series, mean_sic_jas, "JAS"),
    corr_season(nao_series, mean_sic_son, "SON"),
])

print("\n=== Monthly correlation: detrended winter NAO vs detrended Arctic mean SIC (70-90N) ===")
print(df_monthly_mean_sic.to_string(index=False))

print("\n=== Seasonal correlation: detrended winter NAO vs detrended Arctic mean SIC ===")
print(season_results_mean_sic.to_string(index=False))


# ==========================================================
# 7. PLOTS
# ==========================================================
plot_monthly_corr(
    df_monthly_mean_sic,
    "(a) Monthly correlation: winter NAO vs Arctic mean SIC (70–90°N)",
    os.path.join(OUTDIR, "winterNAO_vs_Arctic_meanSIC_monthly_detrended.png")
)
