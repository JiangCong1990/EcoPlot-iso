##################################
### Figure 6 
##################################
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
base_dir = r"C:\Users\cjiang\OneDrive\paper\Second_Paper_IGB\plot_code\data\EcoPlot-iso"

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

    return pd.Series(
        SPI(clip_period(monthly), Scale=scale, Periodicity=12, Distribution="Gamma"),
        index=monthly.index
    )

# ================================
# 3. LOAD & COMPUTE INDICES
# ================================
cl = pd.read_csv(climate_csv, skiprows=range(2, 367))
cl['date'] = pd.to_datetime(cl['date'])
monthly = cl.set_index('date').resample('M').agg({'P_mm': 'sum', 'PET_mm': 'sum'})
wb = monthly['P_mm'] - monthly['PET_mm']

# --- Streamflow SSI
q = pd.read_csv(streamflow_csv, parse_dates=['Date'], dayfirst=True)
q = q.rename(columns={'Date': 'date', 'Discharge': 'Q'}).set_index('date').resample('M').mean()
ssi_1 = spei.ssfi(clip_period(q['Q']), timescale=1, fit_freq="ME")

# --- Groundwater SGI
gw = pd.read_csv(gw_csv, sep=r'\s+', parse_dates=['Time']).set_index('Time')
elev = {'GW3': 55.22, 'GW4': 57.03, 'GW5': 55.46, 'GW7': 54.85, 'GW8': 57.68}
for k, v in elev.items():
    gw[k] = v - gw[k]
gw_mean = gw.resample('M').mean().mean(axis=1)
sgi_1 = -1 * spei.sgi(clip_period(gw_mean), timescale=1, fit_freq="ME")

# --- Climate index
spei6 = pd.Series(
    SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=6),
    index=monthly.index
)
spei6 = spei6.interpolate(limit=1)

spei1 = pd.Series(
    SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=1),
    index=monthly.index
)

spei3 = pd.Series(
    SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=3),
    index=monthly.index
)


spei12 = pd.Series(
    SPEI(clip_period(monthly['P_mm']), monthly['PET_mm'], Axis=0, Scale=12),
    index=monthly.index
)

# --- Simulated storage indices (Forest only)
ssmi_up_f_1 = load_sim_z(sim_forest_csv, 'STO_', 1)
ssmi_dp_f_1 = load_sim_z(sim_forest_csv, 'Sdeep_', 1)
ssmi_up_c_1 = load_sim_z(sim_crop_csv, 'STO_', 1)
ssmi_dp_c_1 = load_sim_z(sim_crop_csv, 'Sdeep_', 1)

# ================================
# 4. PLOT (SIMPLIFIED FIGURE 2)
# ================================
label_map = {
   # "SPEI-1": "SPEI-1  (Standardized Precipitation–Evapotranspiration Index, 1-month)",
   "SPEI-3": "SPEI-3  (Standardized Precipitation–Evapotranspiration Index, 3-month)",
   "SPEI-6": "SPEI-6  (Standardized Precipitation–Evapotranspiration Index, 6-month)",

    "SPEI-12": "SPEI-12 (Standardized Precipitation–Evapotranspiration Index, 12-month)",
    "SSMI Upper (Forest)": "SSMI (Standardized Soil Moisture Index, Forest, 0–10 cm)",
    "SSMI Deeper (Forest)": "SSMI (Standardized Soil Moisture Index, Forest, 30–100 cm)",
  #  "SSMI Upper (Crop)": "SSMI (Standardized Soil Moisture Index, 0–10 cm, Crop)",
 #   "SSMI Deeper (Crop)": "SSMI (Standardized Soil Moisture Index, 30–100 cm, Crop)",
    "SSI": "SSI ((Standardized Streamflow Index)",
    "SGI": "SGI (Standardized Groundwater Index)"
}

series_list = [
 #    (spei1, "SPEI-1"),
    (spei3, "SPEI-3"),
 #    (spei6, "SPEI-6"),
    (spei12, "SPEI-12"),
    (ssmi_up_f_1, "SSMI Upper (Forest)"),
    (ssmi_dp_f_1, "SSMI Deeper (Forest)"),
 #   (ssmi_up_c_1, "SSMI Upper (Crop)"),
  #  (ssmi_dp_c_1, "SSMI Deeper (Crop)"),
    (ssi_1, "SSI"),
    (sgi_1, "SGI"),
]

# Drought event windows to highlight
event_windows = [
    ("2018-03-01", "2019-08-31", "2018–2019 drought"),
    ("2022-03-01", "2023-02-28", "2022–2023 drought"),
]

fig, ax = plt.subplots(
    len(series_list), 1,
    figsize=(12, 2 * len(series_list)),
    sharex=True
)

letters = list("abcdefghijklmnopqrstuvwxyz")

for idx, (series, label) in enumerate(series_list):
    ax_ts = ax[idx]
    full_label = label_map[label]

    # alternating annual background shading
    years = np.arange(series.index.year.min(), series.index.year.max() + 1)
    for y in years:
        if y % 2 == 0:
            ax_ts.axvspan(pd.Timestamp(f"{y}-01-01"), pd.Timestamp(f"{y}-12-31"),
                          color='lightgrey', alpha=0.2, zorder=0)

    # highlight major drought events
    for j, (start, end, txt) in enumerate(event_windows):
        ax_ts.axvspan(pd.Timestamp(start), pd.Timestamp(end),
                      color='darkred', alpha=0.18, zorder=1)
        if idx == 0:
            ax_ts.text(pd.Timestamp(start)-pd.Timedelta(days=300), 1.05, txt, transform=ax_ts.get_xaxis_transform(),fontsize=10, fontweight='bold')

    # time series
    ax_ts.plot(series.index, series, lw=1.5, color='black', label=full_label, zorder=3)

    # reference lines
    ax_ts.axhline(0, color='gray', linestyle='-', lw=0.8)
    for val in [-2, -1, 1, 2]:
        ax_ts.axhline(val, color='black',  alpha=0.5, linestyle='--', lw=0.5)

    # dry / wet shading
    ax_ts.fill_between(series.index, series, -1, where=series < -1,
                       color='red', alpha=0.3, zorder=2)
    ax_ts.fill_between(series.index, series, 1, where=series > 1,
                       color='blue', alpha=0.3, zorder=2)

    ax_ts.set_xlim(pd.Timestamp("2000-01-01"), pd.Timestamp("2024-12-31"))

    # axis range
    if label == "SSI":
        ax_ts.set_ylim(-3.5, 6)
        ax_ts.set_yticks(np.arange(-3, 7, 1))
    else:
        ax_ts.set_ylim(-4.0, 3.0)
        ax_ts.set_yticks(np.arange(-4, 4, 1))

    ax_ts.set_ylabel(label, fontsize=11)
    ax_ts.text(0.01, 0.88, f"({letters[idx]})", transform=ax_ts.transAxes,
               fontweight='bold', fontsize=11)
    ax_ts.legend(loc='lower left', frameon=False, fontsize=10)

ax[-1].set_xlabel("Time", fontsize=11)

# X-axis formatting
ax[-1].xaxis.set_major_locator(mdates.YearLocator(1))
ax[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
for axi in ax:
    axi.tick_params(axis='x', rotation=45)
    
plt.tight_layout()
plt.subplots_adjust(hspace=0.15)
plt.savefig("Figure2_revised_simplified.png", dpi=300, bbox_inches='tight')
plt.show()

# --- Add middle layer (Forest)
ssmi_mid_f_1 = load_sim_z(sim_forest_csv, 'GW_', 1)

# ================================
# 4. PRINT LOWEST 10 VALUES
# ================================
def print_lowest(series, name, n=10):
    s = series.dropna().sort_values().head(n)
    print(f"\n{name} - Lowest {n} values:")
    print(s)

# SSMI upper (forest)
print_lowest(ssmi_up_f_1, "SSMI Upper (Forest)")

print_lowest(ssmi_mid_f_1, "SSMI Middle (Forest)")  

# SSMI deeper (forest)
print_lowest(ssmi_dp_f_1, "SSMI Deeper (Forest)")

# Streamflow index (SSI)
print_lowest(ssi_1, "SSI (Streamflow Index)")

# Groundwater index (SGI)
print_lowest(sgi_1, "SGI (Groundwater Index)")