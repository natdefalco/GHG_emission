# -----------------------------------------------------------
# BLOCK A – Setup and Constants
# ------------------------------------------------------------
# This block defines required libraries, constants for chamber volume and area,
# threshold parameters for flux analysis, and file paths for input and outputs.
# These paths and constants should be edited as needed for each run.

import json
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os
import seaborn as sns
from sklearn.linear_model import LinearRegression

# --- Paths (Update as needed) ---
json_path = r"C:\Users\your_path\.json"
base_name = os.path.splitext(os.path.basename(json_path))[0]
folder = r"C:\Users\your_path\"
os.makedirs(folder, exist_ok=True)

excel_output_path = os.path.join(folder, f"{base_name}_ALLDATA.xlsx")
csv_output_path = os.path.join(folder, f"{base_name}_data_smart_chamber.csv")
save_data_path = os.path.join(folder, "flux_results.csv")
regression_data_path = os.path.join(folder, "regression_parameters.csv")
raw_plot_path = os.path.join(folder, "raw_plot_CH4.png")
flux_plot_path = os.path.join(folder, "flux_plot_CH4.png")
raw_plot_path_co2 = os.path.join(folder, "raw_plot_CO2.png")
flux_plot_path_co2 = os.path.join(folder, "flux_plot_CO2.png")

# --- Constants ---

V = 3672.16 / 1e6  # m³
A = 318 / 1e4   # m²
R_gas = 8.314   # J/mol·K
min_points = 20
min_r2 = 0.80
skip_points = 10

# -----------------------------------------------------------
# BLOCK B – Load and Parse JSON
# ------------------------------------------------------------
# This block loads the JSON file and extracts time series and summary data
# from each repetition. It flattens nested fields and organizes them into rows.

with open(json_path, 'r') as f:
    json_data = json.load(f)

timeseries_rows = []
summary_rows = []
event_counter = 1
original_name = json_data['name']

# --- Parse JSON ---
for dataset in json_data['datasets']:
    for dataset_key, dataset_value in dataset.items():
        reps = dataset_value['reps']
        rep_nums = sorted(reps.keys(), key=lambda x: int(x.split('_')[-1]))
        current_event_name = f"{original_name}_{event_counter}"

        for rep_key in rep_nums:
            rep = reps[rep_key]
            header = rep.get('header', {})
            rep_num = header.get('RepNum')
            footer = rep.get('footer', {})
            fluxes = footer.get('fluxes', [])

            # When RepNum restarts at 1, increment the event counter
            if rep_num == 1 and len(summary_rows) > 0:
                event_counter += 1
                current_event_name = f"{original_name}_{event_counter}"

            # Flatten fluxes
            flux_dict = {}
            for flux in fluxes:
                gas = flux.get('name', '')
                for key, val in flux.items():
                    if key != 'name':
                        flux_dict[f'Flux - {gas} - {key}'] = val

            data = rep['data']
            timestamps = data['timestamp']
            for i in range(len(timestamps)):
                row = {
                    'Dataset Name': current_event_name,
                    'RepNum': rep_num,
                    'Timestamp': timestamps[i],
                    'CH4': data.get('ch4', [None])[i],
                    'CO2': data.get('co2', [None])[i],
                    'Chamber Temperature': data.get('chamber_t', [None])[i],
                    'Chamber Pressure': data.get('chamber_p', [None])[i],
                    'H2O': data.get('h2o', [None])[i],
                    'Soil Temp': data.get('soil_t', [None])[i],
                    'Soil P (c)': data.get('soilp_c', [None])[i],
                    'Soil P (m)': data.get('soilp_m', [None])[i],
                    'Soil P (t)': data.get('soilp_t', [None])[i],
                    'Error': data.get('err', [None])[i],
                }
                timeseries_rows.append(row)

            summary_row = {
                'Dataset Name': current_event_name,
                'RepNum': rep_num,
                'Date': header.get('Date'),
                'Chamber': header.get('Chamber'),
                'Area': header.get('Area'),
                'ChamVolume': header.get('ChamVolume'),
                'IrgaVolume': header.get('IrgaVolume'),
                'TotalVolume': header.get('TotalVolume'),
                'Latitude': header.get('latitude'),
                'Longitude': header.get('longitude'),
                'Altitude': header.get('altitude'),
                'Footer - P_o': footer.get('P_o'),
                'Footer - T_o': footer.get('T_o'),
                'Footer - W_o': footer.get('W_o'),
            }
            summary_row.update(flux_dict)
            summary_rows.append(summary_row)

# -----------------------------------------------------------
# BLOCK C – Create and Save DataFrames
# ------------------------------------------------------------
# Converts parsed lists into DataFrames and saves them to Excel and CSV files.

df_timeseries = pd.DataFrame(timeseries_rows)
df_summary = pd.DataFrame(summary_rows)
df_timeseries.to_csv(csv_output_path, index=False)

with pd.ExcelWriter(excel_output_path, engine='xlsxwriter') as writer:
    df_timeseries.to_excel(writer, sheet_name='timeseries', index=False)
    df_summary.to_excel(writer, sheet_name='summary', index=False)

# -----------------------------------------------------------
# BLOCK D – Organize and Preprocess for Analysis
# ------------------------------------------------------------
# Prepares data by sorting, cleaning column names, and organizing dataset order.

df = df_timeseries.copy()
df.columns = df.columns.str.strip()
df = df.sort_values(by=['Dataset Name', 'RepNum', 'Timestamp']).reset_index(drop=True)

import re

def extract_event_number(name):
    match = re.search(r'_(\d+)$', name)
    return int(match.group(1)) if match else float('inf')

# Group and sort numerically by event number
first_times = df.groupby('Dataset Name')['Timestamp'].min()
dataset_names = sorted(first_times.index.tolist(), key=extract_event_number)

event_lookup = {name: i + 1 for i, name in enumerate(dataset_names)}

flux_data = []
regression_params = []
colors = plt.get_cmap('tab20')

# -----------------------------------------------------------
# BLOCK E – Raw Plots
# ------------------------------------------------------------

raw_plot_paths = {'CH4': raw_plot_path, 'CO2': raw_plot_path_co2}

for gas in ['CH4', 'CO2']:
    value_col = gas
    raw_path = raw_plot_paths[gas]

    plt.figure(figsize=(14, 6))
    time_offset = 0
    for i, name in enumerate(dataset_names):
        event_df = df[df['Dataset Name'] == name]
        reps = sorted(event_df['RepNum'].dropna().unique())
        for rep in reps:
            rep_df = event_df[event_df['RepNum'] == rep].copy()
            rep_df['local_time'] = rep_df['Timestamp'] - rep_df['Timestamp'].min()
            rep_df['shifted_time'] = rep_df['local_time'] + time_offset
            label = name if rep == 1 else None
            plt.plot(rep_df['shifted_time'], rep_df[value_col], color=colors(i), label=label)
            time_offset += rep_df['local_time'].max() + 5
    plt.xlabel('Time (s)')
    plt.ylabel('CH₄ Concentration (ppb)' if gas == 'CH4' else 'CO₂ Concentration (ppm)')
    plt.title(f'{gas} Raw Measurements (Grouped by Events)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(raw_path, dpi=300)
    plt.close()
# -----------------------------------------------------------
# BLOCK F – Flux and Regression Analysis
# ------------------------------------------------------------


for name in dataset_names:
    event_df = df[df['Dataset Name'] == name]
    reps = sorted(event_df['RepNum'].dropna().unique())
    for rep in reps:
        rep_df = event_df[event_df['RepNum'] == rep].copy()
        rep_df['Time'] = rep_df['Timestamp'] - rep_df['Timestamp'].min()
        for gas in ['CH4', 'CO2']:
            value_col = gas
            slope = intercept = r2 = p_val = std_err = flux = None
            T_C = rep_df['Chamber Temperature'].mean()
            T_K = T_C + 273.15
            P_kPa = rep_df['Chamber Pressure'].mean()
            P_Pa = P_kPa*1000  # Convert kPa to Pa

            if len(rep_df) >= skip_points + min_points:
                subset = rep_df.iloc[skip_points:skip_points + min_points]
                if not subset[value_col].isnull().any():
                    slope, intercept, r_val, p_val, std_err = linregress(subset['Time'], subset[value_col])
                    r2 = r_val ** 2
                    mol_per_m3 = P_Pa / (R_gas * T_K)  # Ideal gas law
                    conversion_factor = 3600 * 1e6  # mol → µmol, s → h
                    dc_dt = slope / (1e9 if gas == 'CH4' else 1e6)  # ppb or ppm → mol/mol

                    flux = dc_dt * P_Pa * V / (R_gas * T_K * A) * conversion_factor

            flux_data.append({
                'Dataset Name': name,
                'Event Number': event_lookup[name],
                'RepNum': rep,
                'Gas': gas,
                'R²': r2,
                'Passed R²': r2 is not None and r2 >= min_r2,
                'Chamber Temp (°C)': T_C,
                'Flux (µmol/m²/h)': flux
            })

            if slope is not None and intercept is not None:
                regression_params.append({
                    'Dataset Name': name,
                    'Event Number': event_lookup[name],
                    'RepNum': rep,
                    'Gas': gas,
                    'Slope': slope,
                    'Intercept': intercept,
                    'R²': r2,
                    'p-value': p_val,
                    'Standard Error': std_err
                })

# -----------------------------------------------------------
# BLOCK G – Save Flux and Regression Results
# ------------------------------------------------------------


flux_df = pd.DataFrame(flux_data)
regression_df = pd.DataFrame(regression_params)
flux_df.to_csv(save_data_path, index=False)
regression_df.to_csv(regression_data_path, index=False)

# --- Regression Fits Plots ---
for gas in ['CH4', 'CO2']:
    value_col = gas
    plt.figure(figsize=(14, 6))
    for name in dataset_names:
        event_df = df[df['Dataset Name'] == name]
        reps = sorted(event_df['RepNum'].dropna().unique())
        for rep in reps:
            rep_df = event_df[event_df['RepNum'] == rep].copy()
            rep_df['Time'] = rep_df['Timestamp'] - rep_df['Timestamp'].min()
            # Find regression row for this event, rep, and gas with adequate R²
            row = next((f for f in regression_params if f['Dataset Name'] == name and f['RepNum'] == rep
                        and f['Gas'] == gas and f['R²'] >= min_r2), None)
            if not row:
                continue
            subset = rep_df.iloc[skip_points:skip_points + min_points]
            fit = row['Slope'] * subset['Time'] + row['Intercept']
            plt.scatter(subset['Time'], subset[value_col], s=10, label=f"{name} - Rep {rep}")
            plt.plot(subset['Time'], fit, color='black', linestyle='--', alpha=0.6)
    plt.title(f'{gas} Regression Fits (R² ≥ {min_r2})')
    plt.xlabel("Time (s)")
    plt.ylabel(f"{gas} Concentration (ppm)")
    plt.legend(fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(folder, f"regression_fit_{gas}.png"), dpi=300)
    plt.close()


# -----------------------------------------------------------
# BLOCK E –  Flux Boxplots
# ------------------------------------------------------------

for gas, plot_path in [('CH4', flux_plot_path), ('CO2', flux_plot_path_co2)]:
    plt.figure(figsize=(10, 6))
    sub_df = flux_df[(flux_df['Gas'] == gas) & (flux_df['Passed R²'])]
    event_summary = sub_df.groupby('Dataset Name').agg(
        mean_flux=('Flux (µmol/m²/h)', 'mean'),
        std_flux=('Flux (µmol/m²/h)', 'std'),
        n=('Flux (µmol/m²/h)', 'count')
    ).reset_index()
    x_labels = [f"{name}\n(n={n})" for name, n in zip(event_summary['Dataset Name'], event_summary['n'])]
    plt.errorbar(x_labels, event_summary['mean_flux'], yerr=event_summary['std_flux'],
                 fmt='o', capsize=5, linestyle='None', color='teal')
    plt.title(f'{gas} Flux per Event (± std)')
    plt.xlabel('Event')
    plt.ylabel(f'{gas} Flux (µmol/m²/h)')
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(rotation=0)
    plt.savefig(plot_path, dpi=300)
    plt.close()
summary_df = df_summary.copy()
flux_df_filtered = flux_df.copy()

# -----------------------------------------------------------
# BLOCK F – CH4 Validation
# ------------------------------------------------------------

ch4 = flux_df_filtered[flux_df_filtered["Gas"] == "CH4"].merge(
    summary_df[["Dataset Name", "RepNum", "Flux - ch4 - F_o"]],
    on=["Dataset Name", "RepNum"]
).dropna()

ch4["Your Flux (nmol/m²/s)"] = ch4["Flux (µmol/m²/h)"] * 1e3 / 3600
ch4 = ch4[(ch4["Your Flux (nmol/m²/s)"] < 10) & (ch4["Flux - ch4 - F_o"] < 10)]
X = ch4["Flux - ch4 - F_o"].values.reshape(-1, 1)
y = ch4["Your Flux (nmol/m²/s)"].values
model = LinearRegression().fit(X, y)
r2 = model.score(X, y)

plt.figure(figsize=(6, 6))
sns.regplot(x="Flux - ch4 - F_o", y="Your Flux (nmol/m²/s)", data=ch4, ci=None)
plt.title(f"CH₄ Flux Validation\nR² = {r2:.3f}")
plt.xlabel("Device CH₄ Flux (nmol/m²/s)")
plt.ylabel("Calculated CH₄ Flux (nmol/m²/s)")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(folder, "ch4_flux_validation.png"), dpi=300)
plt.close()

# -----------------------------------------------------------
# BLOCK D – Raw Plots
# ------------------------------------------------------------
# CO₂ VALIDATION
from scipy.stats import zscore
co2 = flux_df_filtered[flux_df_filtered["Gas"] == "CO2"].merge(
    summary_df[["Dataset Name", "RepNum", "Flux - co2 - F_o"]],
    on=["Dataset Name", "RepNum"]
).dropna()

co2["Your Flux (umol/m²/s)"] = co2["Flux (µmol/m²/h)"]/3600

# Optionally filter by R²
co2 = co2[co2["R²"] >= 0.50]

# Z-score filter to remove outliers, only if enough data
if len(co2) >= 5:
    co2_z = co2[["Your Flux (umol/m²/s)", "Flux - co2 - F_o"]].apply(zscore)
    co2 = co2[(co2_z.abs() < 3).all(axis=1)]

# Only run validation if data still remains
if len(co2) >= 2:
    X_co2 = co2["Flux - co2 - F_o"].values.reshape(-1, 1)
    y_co2 = co2["Your Flux (umol/m²/s)"].values
    model_co2 = LinearRegression().fit(X_co2, y_co2)
    r2_co2 = model_co2.score(X_co2, y_co2)

    if r2_co2 < 0.2:
        print(f"⚠️ CO₂ validation weak (R² = {r2_co2:.2f})")

    plt.figure(figsize=(6, 6))
    sns.regplot(x="Flux - co2 - F_o", y="Your Flux (umol/m²/s)", data=co2, ci=None)
    plt.title(f"CO₂ Flux Validation\nR² = {r2_co2:.3f}")
    plt.xlabel("Device CO₂ Flux (µmol/m²/s)")
    plt.ylabel("Calculated CO₂ Flux (µmol/m²/s)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(folder, "co2_flux_validation.png"), dpi=300)
    plt.close()
else:
    print("Not enough CO₂ data points for validation after filtering.")




print(" All steps completed: .json → Excel + Analysis + +Validation + Plots")
