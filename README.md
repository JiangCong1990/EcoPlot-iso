1. Overview

EcoPlot-iso is a tracer-aided ecohydrological modeling framework designed to simulate key ecohydrological and isotopic transformations that characterize water partitioning at the plot scale.
It integrates hydrological, isotopic, and vegetation processes within the Soil–Plant–Atmosphere Continuum (SPAC) to quantify water fluxes, storage changes, and isotopic dynamics under different climatic and land-use conditions.

The model has been applied in diverse climatic and hydrological settings, including:

a one-year simulation in Scotland (Stevenson et al., 2023),

a one-year simulation at the Demnitzer Mill Creek (DMC) catchment in Germany (Landgraf et al., 2023),

a four-year tropical application in Costa Rica (Birkel et al., 2024), and

a long-term (2000–2024) tracer-aided simulation for drought resilience assessment in DMC (Jiang et al., 2024).

Building on these applications, the current version extends EcoPlot-iso for long-term ecohydrological simulations and management scenario analysis, integrating a new depth-dependent root-uptake module.

2. Model Framework and Structure

EcoPlot-iso is a process-based conceptual model that simulates:

Hydrological fluxes: interception, throughfall, infiltration, preferential flow, surface runoff, percolation, groundwater recharge.

Evapotranspiration components: canopy evaporation, soil evaporation, and transpiration.

Isotopic processes: fractionation and mixing in canopy and soil layers.

The vertical structure consists of:

One canopy layer, and

Three soil layers:

0–10 cm (shallow)

10–30 cm (middle)

30–100 cm (deep)

Each layer tracks both water and stable isotope (δ²H or δ¹⁸O) balances.
The isotope module allows explicit separation of evaporation and transpiration losses, improving estimates of water partitioning and residence times.

 <img width="940" height="456" alt="image" src="https://github.com/user-attachments/assets/4e7b24a9-ea3c-473f-9dd0-66702caca42b" />

Figure 2. (a) Schematic representation of the ecohydrological fluxes and water partitioning in the EcoPlot-iso model illustrating major water fluxes and storage components; (b) Conceptual framework and key parameters of the EcoPlot-iso model(Landgraf et al., 2023; Stevenson et al., 2023), highlighting the key ecohydrological processes simulated in this study.


4. Required Input Data
Variable	Unit	Description
Precipitation (P)	mm d⁻¹	Amount of rainfall or snowfall
δ²H / δ¹⁸O of precipitation	‰ VSMOW	Isotopic composition of input water
Air temperature (Tair)	°C	Mean air temperature
Potential evapotranspiration (PET)	mm d⁻¹	Computed e.g. by FAO-56 Penman–Monteith
Relative humidity (RH)	%	Daily mean humidity
Leaf area index (LAI)	m² m⁻²	Vegetation canopy variable
Soil hydraulic parameters	–	Porosity, field capacity, wilting point, Ks, etc.
Vegetation parameters	–	Interception capacity, rooting depth parameter β, rL₁–rL₃, g₁–g₃ (stomatal conductance)


5. Model Outputs
Output	Unit	Description
Ei, Es, Tr	mm d⁻¹	Canopy, soil, and transpiration fluxes
Qs, Recharge	mm d⁻¹	Surface runoff and groundwater recharge
dI, dSTO, dGW, dSdeep	mm	Storage changes
δ²H / δ¹⁸O	‰	Isotopic composition of fluxes and storages
ET/P, Tr/ET	–	Aggregated water-partitioning indices

Outputs are stored in CSV or NetCDF format for analysis and visualization.


7. References

Birkel, C., Jiang, C. et al. (2024). Tracer-aided ecohydrological modeling in humid tropical catchments.

Landgraf, J., Jiang, C. et al. (2023). Long-term ecohydrological simulation in the Demnitz Mill Creek catchment.

Stevenson, E., Jiang, C. et al. (2023). Process-based tracer-aided model for plot-scale hydrological partitioning.

Jiang, C., Tetzlaff, D. et al. (2024). Assessing the drought resilience of different land-management scenarios using a tracer-aided ecohydrological model with variable root-uptake distributions.

Jiang, C., Soulsby, C. et al. (2025). Predicting summer droughts in Central Europe from winter NAO. Nature Water (submitted).
