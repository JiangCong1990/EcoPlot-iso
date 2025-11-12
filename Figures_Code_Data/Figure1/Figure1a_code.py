import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import os

# === Load paths ===
base_dir = r"C:\Users\cjiang\work\Ecoplot_four_sites_root\new_root_equation\new"
climate_all_path = os.path.join(base_dir, r"Ecolpt_root_forest_longterm\Demnitz_inpmod_Fo.csv")
streamflow_path = r"C:\Users\cjiang\work\data\longterm_Q_DM26.csv"
nao_path = r"C:\Users\cjiang\work\data\North_Atlantic_Oscillation(NAO)\North_Atlantic_Oscillation(NAO).csv"
sim_data_path = os.path.join(base_dir, r"Ecolpt_root_forest_longterm\All_Sim_Obs_ForestA_Time_Series.csv")
groundwater_path = r'C:\Users\cjiang\work\data\Groundwater\groundwater_DMC.csv'
lai_path = os.path.join(base_dir, r"Ecolpt_root_forest_longterm\\Demnitz_inpmod_Fo_LAI.csv")

# === Load datasets ===
climate_data = pd.read_csv(climate_all_path, skiprows=range(2, 367))
streamflow_data = pd.read_csv(streamflow_path, sep=',', parse_dates=['Date'])
nao_data = pd.read_csv(nao_path, delimiter=",")
sim_data = pd.read_csv(sim_data_path, parse_dates=['date'])
groundwater_data = pd.read_csv(groundwater_path, sep=r'\s+', parse_dates=['Time'])
lai_data = pd.read_csv(lai_path, parse_dates=['date'])

# === Preprocess ===
climate_data = climate_data.dropna(how='all')
climate_data['date'] = pd.to_datetime(climate_data['date'])
streamflow_data = streamflow_data.rename(columns={'Date': 'date', 'Discharge': 'Q'})
streamflow_data['date'] = pd.to_datetime(streamflow_data['date'], dayfirst=True)
streamflow_data = streamflow_data[(streamflow_data['date'] >= '2000-01-01') & (streamflow_data['date'] <= '2024-12-31')]
nao_data.rename(columns={'INDEX': 'NAO'}, inplace=True)
nao_data['date'] = pd.to_datetime(nao_data['date'].astype(str) + '-' + nao_data['MONTH'].astype(str), errors='coerce')
nao_data = nao_data[(nao_data['date'] >= '1999-01-01') & (nao_data['date'] <= '2024-12-31')]
sim_data.set_index('date', inplace=True)
lai_data.set_index('date', inplace=True)

# === Preprocess groundwater ===
groundwater_data.set_index('Time', inplace=True)
elevations = {'GW3': 55.22, 'GW4': 57.03, 'GW5': 55.46, 'GW7': 54.85, 'GW8': 57.68}
for col in groundwater_data.columns:
    if col in elevations:
        groundwater_data[col] = elevations[col] - groundwater_data[col]
full_date_range = climate_data['date']
groundwater_data = groundwater_data.reindex(full_date_range)

# === Seasonal mean helper ===
def seasonal_mean_pandas(series, months, adjust_december=False):
    df = series.copy().to_frame(name='value')
    if adjust_december:
        df['year'] = df.index.year + (df.index.month == 12).astype(int)
    else:
        df['year'] = df.index.year
    df = df[df.index.month.isin(months)]
    return df.groupby('year')['value'].mean()

def compute_anomalies(series: pd.Series) -> pd.Series:
    """Return anomalies standardized to mean=0, std=1 (z-score)."""
    return (series - series.mean()) / series.mean()
    

# === Climate drivers ===
winter_nao = seasonal_mean_pandas(nao_data.set_index('date')['NAO'], [12, 1, 2], adjust_december=True)
spring_precip = seasonal_mean_pandas(climate_data.set_index('date')['P_mm'], [3, 4, 5])
summer_temp = seasonal_mean_pandas(climate_data.set_index('date')['Air_Temp_oC'], [6, 7, 8])
summer_pet = seasonal_mean_pandas(climate_data.set_index('date')['PET_mm'], [6, 7, 8])
summer_precip = seasonal_mean_pandas(climate_data.set_index('date')['P_mm'], [6, 7, 8])
summer_p_pet = seasonal_mean_pandas(climate_data.set_index('date')['P_mm'] - climate_data.set_index('date')['PET_mm'], [6, 7, 8])

drivers = {
    'Winter NAO': winter_nao,
    'Spring Precip': spring_precip,
    'Summer Temp': summer_temp,
    'Summer PET': summer_pet,
    'Summer Precip': summer_precip,
    'Summer P minus PET': summer_p_pet,
}

# === Responses ===
summer_months = [6, 7, 8]
responses_raw = {
    'SM upper': seasonal_mean_pandas(sim_data[[f'STO_{i}' for i in range(1, 101)]].mean(axis=1), summer_months),
    'SM lower': seasonal_mean_pandas(sim_data[[f'GW_{i}' for i in range(1, 101)]].mean(axis=1), summer_months),
    'SM deeper': seasonal_mean_pandas(sim_data[[f'Sdeep_{i}' for i in range(1, 101)]].mean(axis=1), summer_months),
    'Transpiration': seasonal_mean_pandas(sim_data[[f'Tr_{i}' for i in range(1, 101)]].mean(axis=1), summer_months),
    'GW Recharge': seasonal_mean_pandas(sim_data[[f'Recharge_{i}' for i in range(1, 101)]].mean(axis=1), summer_months),
    'Streamflow': seasonal_mean_pandas(streamflow_data.set_index('date')['Q'], summer_months),
    'GW Depth': seasonal_mean_pandas(groundwater_data.mean(axis=1), summer_months),

}

# Convert only soil moisture to anomalies (z-scores)
responses = responses_raw.copy()
for sm_key in ['SM upper', 'SM lower', 'SM deeper']:
    responses[sm_key] = compute_anomalies(responses_raw[sm_key])

# === Align years ===
years = sorted(set.intersection(*(set(df.index) for df in drivers.values())) & set(responses['SM upper'].index))
driver_df = pd.DataFrame({k: v.loc[years] for k, v in drivers.items()})
response_df = pd.DataFrame({k: v.loc[years] for k, v in responses.items()})

# === Detrend helper ===
def detrend_series(y: pd.Series) -> pd.Series:
    x = np.arange(len(y)).reshape(-1, 1)
    mask = ~y.isna()
    if mask.sum() < 2:
        return y
    model = LinearRegression().fit(x[mask], y[mask])
    trend = model.predict(x)
    return pd.Series(y.values - trend, index=y.index)

# === Correlations (detrended) ===
corr_matrix = pd.DataFrame(index=drivers.keys(), columns=responses.keys())
annot_matrix = pd.DataFrame(index=drivers.keys(), columns=responses.keys())

for d in drivers:
    for r in responses:
        x = driver_df[d]
        y = response_df[r]

        common = x.index.intersection(y.index)
        x_vals = x.loc[common]
        y_vals = y.loc[common]
        valid = x_vals.notna() & y_vals.notna()
        x_vals = x_vals[valid]
        y_vals = y_vals[valid]

        if len(x_vals) >= 2:
            # detrend both before correlation
            x_detr = detrend_series(x_vals)
            y_detr = detrend_series(y_vals)
            r_val, p_val = pearsonr(x_detr, y_detr)
            corr_matrix.loc[d, r] = r_val
            annot_matrix.loc[d, r] = f"{r_val:.2f}{'*' if p_val < 0.05 else ''}"
        else:
            corr_matrix.loc[d, r] = np.nan
            annot_matrix.loc[d, r] = ""

# === Heatmap ===
fig, ax = plt.subplots(figsize=(7, 4))
sns.heatmap(corr_matrix.astype(float), annot=annot_matrix, fmt="", cmap='coolwarm', center=0,
            cbar_kws={'label': 'Pearson Correlation (detrended)'}, ax=ax)

ax.text(-0.3, 0.95, "(a)", transform=ax.transAxes, fontsize=10, fontweight='bold')
ax.set_xlabel("Summer Response Variables")
ax.set_ylabel("Climate Drivers and Variables")
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

fig.tight_layout()
fig.savefig("summer_drought_drivers_correlation_detrended.png", dpi=300)
plt.show()