import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from sklearn.linear_model import LinearRegression
from scipy import stats

# =========================
# CONFIG / PATHS
# =========================
base_path = r"C:\Users\cjiang\work\Ecoplot_four_sites_root\new_root_equation\new"

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
nao_path = r"C:\Users\cjiang\work\data\North_Atlantic_Oscillation(NAO)\North_Atlantic_Oscillation(NAO).csv"

# 'pearson' or 'spearman'
corr_method = "pearson"
# corr_method = "spearman"

# Significance threshold
ALPHA = 0.05

# Toggle detrending
DETREND = True

# =========================
# LOAD DATA
# =========================
with Dataset(results_mean_file, 'r') as nc:
    results_mean = nc.variables['mean_results'][:]      # (ages, scalings, forest_types, components, time)
    forest_ages = nc.variables['Forest_Ages'][:]
    scaling_factors = nc.variables['Scaling_Factors'][:]
    time_steps = nc.variables['Time'][:]

scaling_idx = int(np.where(np.isclose(scaling_factors, 1.0))[0][0])
age_idx = int(np.where(np.isclose(forest_ages, 0))[0][0])

inp = pd.read_csv(inp_csv, sep=",")
inp["Date"] = pd.to_datetime(inp["date"])
inp = inp.iloc[366:]  # spin-up skip
inp["Year"] = inp["Date"].dt.year
inp["Month"] = inp["Date"].dt.month

# ---- NAO: compute DJF means by winter year ----
nao = pd.read_csv(nao_path, delimiter=",", parse_dates=[["date", "MONTH"]])
nao.rename(columns={'date_MONTH': 'date', 'INDEX': 'NAO'}, inplace=True)
nao.set_index('date', inplace=True)
nao['month'] = nao.index.month
nao['year']  = nao.index.year
nao['winter_year'] = np.where(nao['month'] == 12, nao['year'] + 1, nao['year'])

djf = nao[nao['month'].isin([12, 1, 2])]
grp = djf.groupby('winter_year')
complete_djf = grp['NAO'].count() == 3
winter_nao_all = grp['NAO'].mean()[complete_djf]

# ---- Dataframes for CSV-based land uses ----
df_dict = {}
for land_use, file_path in land_use_paths.items():
    df = pd.read_csv(file_path, parse_dates=['date']).set_index('date')
    df['Year'] = df.index.year
    df['Month'] = df.index.month
    df_dict[land_use] = df

# ---- Define analysis window ----
soil_year_min = int(inp["Year"].min())
soil_year_max = int(inp["Year"].max())
full_years = list(range(soil_year_min, soil_year_max + 1))

winter_nao = winter_nao_all[winter_nao_all.index.isin(full_years)]

print(f"Analysis years: {full_years[0]}–{full_years[-1]}")
print(f"First/last winter NAO years kept: {winter_nao.index.min()}–{winter_nao.index.max()}")

# =========================
# HELPERS
# =========================
def seasonal_mask(months):
    return inp["Month"].isin(months)

def compute_csv_anomalies(df, var_list, season, full_years):
    months = [6, 7, 8] if season == "summer" else [3, 4, 5]
    seasonal_df = df[df["Month"].isin(months)]
    series = seasonal_df[var_list].mean(axis=1)
    seasonal_df = pd.DataFrame({'Year': df.loc[seasonal_df.index, 'Year'], 'value': series})
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

def fmt_p(p, alpha=ALPHA):
    return "p<0.05" if p < alpha else f"p={p:.3f}"

def detrend_series(y):
    """Remove linear trend, return residuals."""
    x = np.arange(len(y)).reshape(-1, 1)
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return y
    model = LinearRegression().fit(x[mask], y[mask])
    trend = model.predict(x)
    return y - trend

# =========================
# PLOTTING — 1×3 LAYOUT
# =========================
def plot_scatter_landuse_vs_nao(season: str):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)
    soil_layers = ["STO", "GW", "Sdeep"]
    titles = [
        f"Surface (0–10 cm) — {season.capitalize()}",
        f"Lower (10–30 cm) — {season.capitalize()}",
        f"Deep (30–100 cm) — {season.capitalize()}",
    ]
    panel_labels = ["(c)", "(d)", "(e)"]
    land_uses = ["Agroforest", "Broadleaf", "Conifer", "Crop", "Grassland"]
    comp_map = {"STO": 16, "GW": 17, "Sdeep": 18}

    nao_series = winter_nao.copy()
    if DETREND:
        nao_series = pd.Series(detrend_series(nao_series.values), index=nao_series.index)

    for i, (layer, title) in enumerate(zip(soil_layers, titles)):
        ax = axes[i]

        for lu in land_uses:
            if lu in ["Agroforest", "Broadleaf", "Conifer"]:
                forest_idx = ["Agroforest", "Broadleaf", "Conifer"].index(lu)
                anomaly = compute_nc_anomalies(results_mean, comp_map[layer], forest_idx, season)
            else:
                var_list = [f"{layer}_{k}" for k in range(1, 101)]
                anomaly = compute_csv_anomalies(df_dict[lu], var_list, season, full_years)

            aligned = anomaly[anomaly.index.isin(nao_series.index)]
            nao_aligned = nao_series[nao_series.index.isin(aligned.index)]
            print(nao_aligned)
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
                line_style  = "-" if sig else "--"
                line_alpha  = 0.9 if sig else 0.7

                ax.scatter(x, y, s=22, alpha=point_alpha,
                           color=colors[lu], edgecolor='none', label=label)

                x_line = np.linspace(x.min(), x.max(), 100)
                y_line = model.predict(x_line.reshape(-1, 1))
                ax.plot(x_line, y_line, linestyle=line_style,
                        color=colors[lu], alpha=line_alpha)

        ax.axhline(0, color="k", linestyle=":", lw=1)
        if i == 0:
            ax.set_ylabel("Relative anomaly (unitless)", fontsize=10)
        ax.set_xlabel("Last-winter NAO (DJF mean)", fontsize=10)
        ax.text(-0.1, 0.98, panel_labels[i], transform=ax.transAxes,
                fontsize=12, fontweight='bold')
        ax.legend(loc='upper center', bbox_to_anchor=(0.2, -0.18),
                  ncol=1, frameon=False, fontsize=10)
        ax.tick_params(labelleft=True)

    fig = axes[0].figure
    fig.text(1.02, 0.09, "* denotes p<0.05, n=25", ha="right", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(f"soil_moisture_scatter_vs_NAO_{season}_{corr_method}_1x3_detrended.png",
                dpi=300, bbox_inches="tight")
    plt.show()

# =========================
# RUN
# =========================
plot_scatter_landuse_vs_nao("summer")