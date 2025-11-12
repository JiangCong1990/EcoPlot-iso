import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.stats import gamma, norm, fisk
from gma.climet.Index import SPI, SPEI
import spei

# ================================
# 1. PATHS
# ================================
base_dir = r"C:\Users\cjiang\work\Ecoplot_four_sites_root\new_root_equation\new"
climate_csv = rf"{base_dir}\Ecolpt_root_forest_longterm\Demnitz_inpmod_Fo.csv"
streamflow_csv = r"C:\Users\cjiang\work\data\longterm_Q_DM26.csv"
gw_csv = r"C:\Users\cjiang\work\data\Groundwater\groundwater_DMC.csv"
sim_forest_csv = rf"{base_dir}\Ecolpt_root_forest_longterm\All_Sim_Obs_ForestA_Time_Series.csv"
sim_crop_csv = rf"{base_dir}\Ecolpt_root_crop_longterm\All_Sim_Obs_ForestA_Time_Series.csv"

# ================================
# 2. HELPER FUNCTIONS
# ================================
def z_index_scale(series, scale):
    rolled = series.rolling(scale, min_periods=scale).mean()
    return (rolled - rolled.mean()) / rolled.std()

def gamma_index(series, scale):
    roll = series.rolling(scale, min_periods=scale).sum()
    valid = roll.dropna()
    shift = abs(valid.min()) + 1e-6
    shape, loc, scale_fit = gamma.fit(valid + shift, floc=0)
    cdf = gamma.cdf(roll + shift, shape, loc=0, scale=scale_fit)
    return pd.Series(norm.ppf(cdf.clip(1e-6, 1-1e-6)), index=series.index)

def loglogistic_index(series, scale):
    roll = series.rolling(scale, min_periods=scale).sum()
    valid = roll.dropna()
    shift = abs(valid.min()) + 1e-6 if valid.min() <= 0 else 0
    shape, loc, scale_fit = fisk.fit(valid + shift, floc=0)
    cdf = fisk.cdf(roll + shift, shape, loc=0, scale=scale_fit)
    return pd.Series(norm.ppf(cdf.clip(1e-6, 1-1e-6)), index=series.index)

def clip_period(series, start='2000-01-01', end='2024-12-31'):
    return series.loc[start:end]

def load_sim_z(path, prefix, scale, anomaly=False):
    df = pd.read_csv(path, parse_dates=['date']).set_index('date')
    cols = [c for c in df.columns if c.startswith(prefix) and c[len(prefix):].isdigit()]
    monthly = df[cols].mean(axis=1).resample('M').mean()
    if prefix.startswith('STO_'):
        monthly = monthly / 10
    elif prefix.startswith('GW_'):
        monthly = monthly / 20
    elif prefix.startswith('Sdeep_'):
        monthly = monthly / 70
    if anomaly:
        monthly = monthly.groupby(monthly.index.month).transform(lambda x: x - x.mean())
    return pd.Series(SPI(clip_period(monthly), Scale=scale, Periodicity=12, Distribution="Gamma"),
                     index=monthly.index)

# ================================
# 3. LOAD & COMPUTE INDICES
# ================================
cl = pd.read_csv(climate_csv, skiprows=range(2, 367))
cl['date'] = pd.to_datetime(cl['date'])
monthly = cl.set_index('date').resample('M').agg({'P_mm':'sum','PET_mm':'sum'})
wb = monthly['P_mm'] - monthly['PET_mm']

# --- Streamflow SSI
q = pd.read_csv(streamflow_csv, parse_dates=['Date'], dayfirst=True)
q = q.rename(columns={'Date': 'date', 'Discharge': 'Q'}).set_index('date').resample('M').mean()
ssi_1 = spei.ssfi(clip_period(q['Q']), timescale=1, fit_freq="ME")

# --- Groundwater SGI
gw = pd.read_csv(gw_csv, sep=r'\s+', parse_dates=['Time']).set_index('Time')
elev = {'GW3':55.22,'GW4':57.03,'GW5':55.46,'GW7':54.85,'GW8':57.68}
for k, v in elev.items():
    gw[k] = v - gw[k]
gw_mean = gw.resample('M').mean().mean(axis=1)
sgi_1 = -1 * spei.sgi(clip_period(gw_mean), timescale=1, fit_freq="ME")

# --- Climate indices
spei1  = pd.Series(SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=1),  index=monthly.index)
spei3  = pd.Series(SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=3),  index=monthly.index)
spei12 = pd.Series(SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=12), index=monthly.index)

# --- Simulated storage indices
ssmi_up_f_1 = load_sim_z(sim_forest_csv, 'STO_', 1)
ssmi_dp_f_1 = load_sim_z(sim_forest_csv, 'Sdeep_', 1)
ssmi_up_c_1 = load_sim_z(sim_crop_csv, 'STO_', 1)
ssmi_dp_c_1 = load_sim_z(sim_crop_csv, 'Sdeep_', 1)

# ================================
# 4. PLOT
# ================================
drought_years = [2018, 2022]
months = np.arange(1, 13)

# Define mapping from short to full names for clarity in labels/legends
label_map = {
    "SPEI-3": "SPEI-3  (Standardized Precipitation–Evapotranspiration Index, 3-month)",
    "SPEI-12": "SPEI-12 (Standardized Precipitation–Evapotranspiration Index, 12-month)",
    "SSMI Upper (Forest)": "SSMI (Standardized Soil Moisture Index, 0–10 cm, Forest)",
    "SSMI Deeper (Forest)": "SSMI (Standardized Soil Moisture Index, 30–100 cm, Forest)",
    "SSMI Upper (Crop)": "SSMI (Standardized Soil Moisture Index, 0–10 cm, Crop)",
    "SSMI Deeper (Crop)": "SSMI (Standardized Soil Moisture Index, 30–100 cm, Crop)",
    "SSI": "SSI (Standardized Streamflow Index)",
    "SGI": "SGI (Standardized Groundwater Index)"
}

series_list = [
  #  (spei1,  "SPEI-1"),
    (spei3,  "SPEI-3"),
    (spei12, "SPEI-12"),
    (ssmi_up_f_1, "SSMI Upper (Forest)"),
    (ssmi_dp_f_1, "SSMI Deeper (Forest)"),
    (ssmi_up_c_1, "SSMI Upper (Crop)"),
    (ssmi_dp_c_1, "SSMI Deeper (Crop)"),
    (ssi_1,  "SSI"),
    (sgi_1,  "SGI"),
]

fig, ax = plt.subplots(
    len(series_list), 2,
    figsize=(12, 1.8 * len(series_list)),
    gridspec_kw={'width_ratios': [4, 1]},
    sharex='col'
)

letters = list("abcdefghijklmnopqrstuvwxyz")  # will continue i, j, k, ...

for idx, (series, label) in enumerate(series_list):
    ax_ts, ax_clim = ax[idx, 0], ax[idx, 1]
    full_label = label_map[label]  # expanded version


    # ---- Time series panel ----
    years = np.arange(series.index.year.min(), series.index.year.max() + 1)
    for y in years:
        if y % 2 == 0:
            ax_ts.axvspan(pd.Timestamp(f"{y}-01-01"), pd.Timestamp(f"{y}-12-31"),
                          color='lightgrey', alpha=0.25, zorder=0)
    ax_ts.plot(series.index, series, lw=1.2, color='black', label=full_label)
    ax_ts.axhline(0, color='gray', linestyle='-', lw=0.8)
    ax_ts.axhline(-1, color='black', linestyle=':', lw=1.0)
    ax_ts.axhline(1, color='black', linestyle=':', lw=1.0)
    ax_ts.axhline(-2, color='black', linestyle=':', lw=1.0)
    ax_ts.axhline(2, color='black', linestyle=':', lw=1.0)
    ax_ts.fill_between(series.index, series, -1, where=series < -1, color='red', alpha=0.4)
    ax_ts.fill_between(series.index, series, 1, where=series > 1, color='blue', alpha=0.3)
    ax_ts.set_xlim(pd.Timestamp("2000-01-01"), pd.Timestamp("2024-12-31"))
    if label == "SSI":   # custom ylim for SSI
        ax_ts.set_ylim(-3.5, 6)
        ax_ts.set_yticks(np.arange(-3, 7, 1))
    else:                # default for others
        ax_ts.set_ylim(-4.5, 3)
        ax_ts.set_yticks(np.arange(-4, 4, 1))
    ax_ts.set_ylabel(label, fontsize=11)
    ax_ts.text(0.01, 0.90, f"({letters[idx]})", transform=ax_ts.transAxes, fontweight='bold')
    ax_ts.legend(loc='lower left')

    # ---- Climatology panel ----
    monthly_grp = series.groupby(series.index.month)
    for yr in series.index.year.unique():
        yearly = series[series.index.year == yr]
        mon_vals = yearly.groupby(yearly.index.month).mean()
        ax_clim.plot(months, mon_vals.reindex(months), color='gray', lw=1.0, alpha=1.0)
    colors = ['red', 'blue', 'orange', 'purple']
    for i, dy in enumerate(drought_years):
        sel = series[series.index.year == dy]
        mon_vals = sel.groupby(sel.index.month).mean()
        ax_clim.plot(months, mon_vals.reindex(months), lw=2.5, color=colors[i % len(colors)], label=f"{dy}")
    ax_clim.set_xlim(1, 12)
    ax_clim.set_xticks(months)
    if label == "SSI":   # custom ylim for SSI
        ax_clim.set_ylim(-3.5, 6)
        ax_clim.set_yticks(np.arange(-3, 7, 1))
    else:                # default for others
        ax_clim.set_ylim(-4.5, 3)
        ax_clim.set_yticks(np.arange(-4, 4, 1))
    #ax_clim.set_ylabel(label, fontsize=11)
    ax_clim.legend(loc='lower left', ncol=2, fontsize=8)
    ax_clim.axhline(-1, color='black', linestyle=':', lw=1.0)
    ax_clim.axhline(1, color='black', linestyle=':', lw=1.0)
    ax_clim.axhline(-2, color='black', linestyle=':', lw=1.0)
    ax_clim.axhline(2, color='black', linestyle=':', lw=1.0)

# X-axis formatting for time series
for ax_ts in ax[:, 0]:
    ax_ts.xaxis.set_major_locator(mdates.YearLocator(1))
    ax_ts.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax_ts.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig("merged_drought_indices.png", dpi=300, bbox_inches='tight')
plt.show()
print(spei1.nsmallest(10))