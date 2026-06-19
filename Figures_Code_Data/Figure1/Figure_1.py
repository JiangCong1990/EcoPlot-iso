##################################
### Figure 1 
##################################
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from netCDF4 import Dataset
from scipy.stats import pearsonr
from scipy import stats
from sklearn.linear_model import LinearRegression
from matplotlib.gridspec import GridSpec

# =========================
# CONFIG / PATHS
# =========================
base_dir = r"C:\Users\cjiang\work\Ecoplot_four_sites_root\new_root_equation\new"
base_dir = r"C:\Users\cjiang\OneDrive\paper\Second_Paper_IGB\plot_code\data\EcoPlot-iso"
base_path = base_dir

# ---- Heatmap-code paths ----
climate_all_path = os.path.join(base_dir, r"Ecolpt_root_forest_longterm\Demnitz_inpmod_Fo.csv")
streamflow_path = r"C:\Users\cjiang\work\data\longterm_Q_DM26.csv"
nao_path = r"C:\Users\cjiang\work\data\North_Atlantic_Oscillation(NAO)\North_Atlantic_Oscillation(NAO).csv"
sim_data_path = os.path.join(base_dir, r"Ecolpt_root_forest_longterm\All_Sim_Obs_ForestA_Time_Series.csv")
groundwater_path = r"C:\Users\cjiang\work\data\Groundwater\groundwater_DMC.csv"
lai_path = os.path.join(base_dir, r"Ecolpt_root_forest_longterm\Demnitz_inpmod_Fo_LAI.csv")

# ---- Scatter-code paths ----
land_use_paths = {
    "Broadleaf": fr"{base_path}\Ecolpt_root_forest_longterm\All_Sim_Obs_ForestA_Time_Series.csv",
    "Crop":      fr"{base_path}\Ecolpt_root_crop_longterm\All_Sim_Obs_ForestA_Time_Series.csv",
    "Grassland": fr"{base_path}\Ecolpt_root_grassland_longterm\All_Sim_Obs_ForestA_Time_Series.csv",
}

colors = {
    "Agroforest": "#1f77b4",
    "Broadleaf":  "#ff7f0e",
    "Conifer":    "#2ca02c",
    "Crop":       "#d62728",
    "Grassland":  "#9467bd",
}

results_mean_file = fr"{base_path}\Ecolpt_root_forest_longterm\Modelsim_results_mean.nc"
inp_csv = fr"{base_path}\Ecolpt_root_forest_longterm\Demnitz_inpmod_Fo.csv"

corr_method = "pearson"
ALPHA = 0.05
DETREND = True

# =========================
# LOAD DATA
# =========================
# ---- Common / heatmap data ----
climate_data = pd.read_csv(climate_all_path, skiprows=range(2, 367))
streamflow_data = pd.read_csv(streamflow_path, sep=",", parse_dates=["Date"])
nao_data = pd.read_csv(nao_path, delimiter=",")
sim_data = pd.read_csv(sim_data_path, parse_dates=["date"])
groundwater_data = pd.read_csv(groundwater_path, sep=r"\s+", parse_dates=["Time"])
lai_data = pd.read_csv(lai_path, parse_dates=["date"])

# ---- Scatter data ----
with Dataset(results_mean_file, "r") as nc:
    results_mean = nc.variables["mean_results"][:]      # (ages, scalings, forest_types, components, time)
    forest_ages = nc.variables["Forest_Ages"][:]
    scaling_factors = nc.variables["Scaling_Factors"][:]
    time_steps = nc.variables["Time"][:]

scaling_idx = int(np.where(np.isclose(scaling_factors, 1.0))[0][0])
age_idx = int(np.where(np.isclose(forest_ages, 0))[0][0])

inp = pd.read_csv(inp_csv, sep=",")
inp["Date"] = pd.to_datetime(inp["date"])
inp = inp.iloc[366:]  # spin-up skip
inp["Year"] = inp["Date"].dt.year
inp["Month"] = inp["Date"].dt.month

# ---- Scatter dataframes for CSV land uses ----
df_dict = {}
for land_use, file_path in land_use_paths.items():
    df = pd.read_csv(file_path, parse_dates=["date"]).set_index("date")
    df["Year"] = df.index.year
    df["Month"] = df.index.month
    df_dict[land_use] = df

# =========================
# PREPROCESS
# =========================
climate_data = climate_data.dropna(how="all")
climate_data["date"] = pd.to_datetime(climate_data["date"])

streamflow_data = streamflow_data.rename(columns={"Date": "date", "Discharge": "Q"})
streamflow_data["date"] = pd.to_datetime(streamflow_data["date"], dayfirst=True)
streamflow_data = streamflow_data[
    (streamflow_data["date"] >= "2000-01-01") &
    (streamflow_data["date"] <= "2024-12-31")
]

# --- NAO for heatmap panel (a) ---
nao_data.rename(columns={"INDEX": "NAO"}, inplace=True)
nao_data["date"] = pd.to_datetime(
    nao_data["date"].astype(str) + "-" + nao_data["MONTH"].astype(str),
    errors="coerce"
)
nao_data = nao_data[
    (nao_data["date"] >= "1999-12-01") &
    (nao_data["date"] <= "2024-12-31")
]

sim_data = sim_data[
    (sim_data["date"] >= "2000-01-01") &
    (sim_data["date"] <= "2024-12-31")
]
sim_data.set_index("date", inplace=True)

lai_data.set_index("date", inplace=True)

# ---- groundwater ----
groundwater_data.set_index("Time", inplace=True)
elevations = {"GW3": 55.22, "GW4": 57.03, "GW5": 55.46, "GW7": 54.85, "GW8": 57.68}
for col in groundwater_data.columns:
    if col in elevations:
        groundwater_data[col] = elevations[col] - groundwater_data[col]

full_date_range = climate_data["date"]
groundwater_data = groundwater_data.reindex(full_date_range)

# =========================
# HELPERS
# =========================
def seasonal_mean_pandas(series, months, adjust_december=False):
    df = series.copy().to_frame(name="value")
    if adjust_december:
        df["year"] = df.index.year + (df.index.month == 12).astype(int)
    else:
        df["year"] = df.index.year
    df = df[df.index.month.isin(months)]
    return df.groupby("year")["value"].mean()

def monthly_mean(series: pd.Series, month: int, adjust_dec: bool = False) -> pd.Series:
    df = series.copy().to_frame(name="value")
    if adjust_dec:
        df["year"] = df.index.year + (df.index.month == 12).astype(int)
    else:
        df["year"] = df.index.year
    df = df[df.index.month == month]
    return df.groupby("year")["value"].mean()

def compute_anomalies(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.mean()

def detrend_series(y):
    """Remove linear trend, return residuals."""
    if isinstance(y, pd.Series):
        y_arr = y.values.astype(float)
        idx = y.index
    else:
        y_arr = np.asarray(y, dtype=float)
        idx = None

    x = np.arange(len(y_arr)).reshape(-1, 1)
    mask = ~np.isnan(y_arr)

    if mask.sum() < 2:
        return y if isinstance(y, pd.Series) else y_arr

    model = LinearRegression().fit(x[mask], y_arr[mask])
    trend = model.predict(x)
    resid = y_arr - trend

    if isinstance(y, pd.Series):
        return pd.Series(resid, index=idx)
    return resid

def seasonal_mask(months):
    return inp["Month"].isin(months)

def compute_csv_anomalies(df, var_list, season, full_years):
    months = [6, 7, 8] if season == "summer" else [3, 4, 5]
    seasonal_df = df[df["Month"].isin(months)]
    series = seasonal_df[var_list].mean(axis=1)
    seasonal_df = pd.DataFrame({"Year": df.loc[seasonal_df.index, "Year"], "value": series})
    annual_means = seasonal_df.groupby("Year")["value"].mean()
    clim = annual_means.mean()
    anomalies = (annual_means - clim) / (clim if clim != 0 else 1.0)
    return anomalies.reindex(full_years)

def compute_nc_anomalies(results, component, forest_idx, season):
    months = [6, 7, 8] if season == "summer" else [3, 4, 5]
    sm = results[age_idx, scaling_idx, forest_idx, component, :]
    sm_season = sm[seasonal_mask(months).values]
    df = pd.DataFrame({"Date": inp.loc[seasonal_mask(months), "Date"], "val": sm_season})
    df["Year"] = df["Date"].dt.year
    annual = df.groupby("Year")["val"].mean()
    clim = annual.mean()
    anomalies = (annual - clim) / (clim if clim != 0 else 1.0)
    return anomalies.reindex(full_years)

def do_corr(x, y, method="pearson"):
    if method == "spearman":
        r, p = stats.spearmanr(x, y)
    else:
        r, p = stats.pearsonr(x, y)
    return r, p

# =========================
# PANEL (a): MONTHLY HEATMAP FROM CODE B
# =========================
# Winter NAO (DJF mean, Dec assigned to following year)
winter_nao_heat = (
    pd.concat([
        monthly_mean(nao_data.set_index("date")["NAO"], 12, adjust_dec=True),
        monthly_mean(nao_data.set_index("date")["NAO"], 1,  adjust_dec=True),
        monthly_mean(nao_data.set_index("date")["NAO"], 2,  adjust_dec=True),
    ], axis=1)
    .mean(axis=1)
)

def monthly_vars(month: int, adjust_dec: bool = False):
    mname = pd.to_datetime(str(month), format="%m").strftime("%b")
    return {
        f"{mname}_Precip": monthly_mean(climate_data.set_index("date")["P_mm"], month, adjust_dec),
        f"{mname}_Temp": monthly_mean(climate_data.set_index("date")["Air_Temp_oC"], month, adjust_dec),
        f"{mname}_PET": monthly_mean(climate_data.set_index("date")["PET_mm"], month, adjust_dec),
        f"{mname}_SM_upper": monthly_mean(sim_data[[f"STO_{i}" for i in range(1, 101)]].mean(axis=1), month, adjust_dec),
        f"{mname}_SM_lower": monthly_mean(sim_data[[f"GW_{i}" for i in range(1, 101)]].mean(axis=1), month, adjust_dec),
        f"{mname}_SM_deeper": monthly_mean(sim_data[[f"Sdeep_{i}" for i in range(1, 101)]].mean(axis=1), month, adjust_dec),
        f"{mname}_Recharge": monthly_mean(sim_data[[f"Recharge_{i}" for i in range(1, 101)]].mean(axis=1), month, adjust_dec),
        f"{mname}_Streamflow": monthly_mean(streamflow_data.set_index("date")["Q"], month, adjust_dec),
        f"{mname}_GW_Depth": monthly_mean(groundwater_data.mean(axis=1), month, adjust_dec),
        f"{mname}_Broadleaf LAI": monthly_mean(lai_data["LAI_2"], month, adjust_dec),
        f"{mname}_Coniferous LAI": monthly_mean(lai_data["LAI_3"], month, adjust_dec),
    }

desired_order = [
    "Precip",
    "Temp",
    "PET",
    "SM_upper",
    "SM_lower",
    "SM_deeper",
    "Recharge",
    "Streamflow",
    "GW_Depth",
    "Broadleaf LAI",
    "Coniferous LAI"
]

all_month_vars = {}
for m in range(1, 13):
    all_month_vars.update(monthly_vars(m, adjust_dec=(m == 12)))

response_df_a = pd.DataFrame(all_month_vars)

# use Jul_SM_upper as anchor, same logic as code B
anchor = response_df_a.filter(like="Jul_SM_upper").dropna()
valid_years_a = sorted(set(winter_nao_heat.index) & set(anchor.index))

winter_nao_heat = winter_nao_heat.loc[valid_years_a]
driver_df_a = pd.DataFrame({"Winter_NAO": winter_nao_heat})
response_df_a = response_df_a.loc[valid_years_a]

correlation_data, annotation_data = [], []
for col in response_df_a.columns:
    month = col.split("_")[0]
    variable = "_".join(col.split("_")[1:])

    x = driver_df_a["Winter_NAO"]
    y = response_df_a[col]
    valid = x.notna() & y.notna()

    if valid.sum() >= 2:
        x_use = x[valid]
        y_use = y[valid]

        if DETREND:
            x_use = detrend_series(x_use)
            y_use = detrend_series(y_use)

        r, p = pearsonr(x_use, y_use)
        correlation_data.append((month, variable, r))
        annotation_data.append((month, variable, f"{r:.2f}{'*' if p < 0.05 else ''}"))

corr_df_a = pd.DataFrame(correlation_data, columns=["Month", "Variable", "Correlation"])
annot_df_a = pd.DataFrame(annotation_data, columns=["Month", "Variable", "Annot"])

month_order = ["Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"]
heatmap_corr_a = corr_df_a.pivot(index="Month", columns="Variable", values="Correlation").loc[month_order]
heatmap_annot_a = annot_df_a.pivot(index="Month", columns="Variable", values="Annot").loc[month_order]

heatmap_corr_a = heatmap_corr_a[desired_order]
heatmap_annot_a = heatmap_annot_a[desired_order]

# =========================
# SCATTER DATA (same as code A)
# =========================
nao = pd.read_csv(nao_path, delimiter=",", parse_dates=[["date", "MONTH"]])
nao.rename(columns={"date_MONTH": "date", "INDEX": "NAO"}, inplace=True)
nao.set_index("date", inplace=True)
nao["month"] = nao.index.month
nao["year"] = nao.index.year
nao["winter_year"] = np.where(nao["month"] == 12, nao["year"] + 1, nao["year"])

djf = nao[nao["month"].isin([12, 1, 2])]
grp = djf.groupby("winter_year")
complete_djf = grp["NAO"].count() == 3
winter_nao_all = grp["NAO"].mean()[complete_djf]

soil_year_min = int(inp["Year"].min())
soil_year_max = int(inp["Year"].max())
full_years = list(range(soil_year_min, soil_year_max + 1))
winter_nao = winter_nao_all[winter_nao_all.index.isin(full_years)]

print(f"Analysis years: {full_years[0]}–{full_years[-1]}")
print(f"First/last winter NAO years kept: {winter_nao.index.min()}–{winter_nao.index.max()}")

# =========================
# COMBINED FIGURE
# top row: centered, narrower (a)
# bottom row: (b) (c) (d)
# =========================
fig = plt.figure(figsize=(12, 11))
gs = GridSpec(
    2, 3,
    figure=fig,
    height_ratios=[1.3, 1.0],
    width_ratios=[1, 1, 1],
    hspace=0.2,
    wspace=0.25,
)

# narrower, centered panel (a)
gs_top = gs[0, :].subgridspec(1, 12)
ax_heat = fig.add_subplot(gs_top[0, 1:11])

# bottom row unchanged
ax_sc1 = fig.add_subplot(gs[1, 0])
ax_sc2 = fig.add_subplot(gs[1, 1])
ax_sc3 = fig.add_subplot(gs[1, 2])

# keep b c d similar
ax_sc1.set_box_aspect(0.85)
ax_sc2.set_box_aspect(0.85)
ax_sc3.set_box_aspect(0.85)

# ---- (a) monthly heatmap ----
hm = sns.heatmap(
    heatmap_corr_a.astype(float),
    annot=heatmap_annot_a,
    fmt="",
    cmap="coolwarm",
    center=0,
    cbar=True,
    cbar_kws={
        "label": "Pearson correlation",
        "shrink": 0.88,
        "pad": 0.02
    },
    ax=ax_heat
)

ax_heat.text(-0.08, 1.03, "(a)", transform=ax_heat.transAxes,
             fontsize=12, fontweight="bold")
ax_heat.set_xlabel("Variables", fontsize=10)
ax_heat.set_ylabel("Month", fontsize=10)
plt.setp(ax_heat.get_xticklabels(), rotation=25, ha="center")

# make colorbar label/ticks a bit cleaner
cbar = hm.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
cbar.set_label("Pearson correlation", fontsize=10)

# ---- (b)(c)(d) scatter panels ----
scatter_axes = [ax_sc1, ax_sc2, ax_sc3]
soil_layers = ["STO", "GW", "Sdeep"]
titles = [
    "Surface layer (0–10 cm)",
    "Lower layer (10–30 cm)",
    "Deep layer (30–100 cm)",
]
panel_labels = ["(b)", "(c)", "(d)"]
land_uses = ["Agroforest", "Broadleaf", "Conifer", "Crop", "Grassland"]
comp_map = {"STO": 16, "GW": 17, "Sdeep": 18}

nao_series = winter_nao.copy()
if DETREND:
    nao_series = pd.Series(detrend_series(nao_series.values), index=nao_series.index)

for i, (ax, layer, title, plab) in enumerate(zip(scatter_axes, soil_layers, titles, panel_labels)):
    for lu in land_uses:
        if lu in ["Agroforest", "Broadleaf", "Conifer"]:
            forest_idx = ["Agroforest", "Broadleaf", "Conifer"].index(lu)
            anomaly = compute_nc_anomalies(results_mean, comp_map[layer], forest_idx, "summer")
        else:
            var_list = [f"{layer}_{k}" for k in range(1, 101)]
            anomaly = compute_csv_anomalies(df_dict[lu], var_list, "summer", full_years)

        aligned = anomaly[anomaly.index.isin(nao_series.index)]
        nao_aligned = nao_series[nao_series.index.isin(aligned.index)]

        y = aligned.values
        if DETREND:
            y = detrend_series(y)
        x = nao_aligned.values

        mask = (~np.isnan(x)) & (~np.isnan(y))
        n = int(mask.sum())

        if n >= 3:
            r, p = do_corr(x[mask], y[mask], corr_method)
            model = LinearRegression().fit(x[mask].reshape(-1, 1), y[mask])
            slope = float(model.coef_[0])

            sig = (p < ALPHA)
            star = "*" if sig else ""
            label = f"{lu} (k={slope:.2f}, r={r:.2f}{star})"

            point_alpha = 0.95 if sig else 0.6
            line_style = "-" if sig else "--"
            line_alpha = 0.9 if sig else 0.7

            ax.scatter(
                x, y, s=22, alpha=point_alpha,
                color=colors[lu], edgecolor="none", label=label
            )

            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = model.predict(x_line.reshape(-1, 1))
            ax.plot(
                x_line, y_line, linestyle=line_style,
                color=colors[lu], alpha=line_alpha
            )

    ax.axhline(0, color="k", linestyle=":", lw=1)
    if i == 0:
        ax.set_ylabel("Soil moisture relative anomaly (-)", fontsize=10)
    else:
        ax.set_ylabel("")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Winter NAO (DJF mean)", fontsize=10)
    ax.text(-0.16, 1.03, plab, transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=1,
        frameon=True,
        fontsize=9
    )
    ax.tick_params(labelleft=True)

fig.text(0.9, -0.01, "* denotes p<0.05", ha="right", va="bottom", fontsize=10)

plt.savefig("combined_monthlyheatmap_scatter_2row3col_detrended.png", dpi=500, bbox_inches="tight")
plt.show()

# =========================
# SAVE DATA
# =========================
output_dir = os.path.join(base_dir, "plotting_data")
os.makedirs(output_dir, exist_ok=True)

heatmap_corr_a.to_csv(os.path.join(output_dir, "panel_a_monthly_corr_matrix.csv"), float_format="%.3f")
heatmap_annot_a.to_csv(os.path.join(output_dir, "panel_a_monthly_annot_matrix.csv"), index=True)
driver_df_a.to_csv(os.path.join(output_dir, "panel_a_driver_data.csv"), float_format="%.3f")
response_df_a.to_csv(os.path.join(output_dir, "panel_a_response_data.csv"), float_format="%.3f")

print("✅ Combined figure saved as: combined_monthlyheatmap_scatter_2row3col_detrended.png")
print("✅ Plotting data saved to:", output_dir)