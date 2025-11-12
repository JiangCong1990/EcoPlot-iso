1. Overview

EcoPlot-iso is a tracer-aided ecohydrological modeling framework designed to simulate key ecohydrological and isotopic transformations that characterize water partitioning at the plot scale.
It integrates hydrological, isotopic, and vegetation processes within the Soil–Plant–Atmosphere Continuum (SPAC) to quantify water fluxes, storage changes, and isotopic dynamics under different climatic and land-use conditions.

The model has been applied in diverse climatic and hydrological settings, including:

** a one-year simulation in Scotland (Stevenson et al., 2023),

** a one-year simulation at the Demnitzer Mill Creek (DMC) catchment in Germany (Landgraf et al., 2023),

** a four-year tropical application in Costa Rica (Birkel et al., 2024), and

** a long-term (2000–2024) tracer-aided simulation for drought resilience assessment in DMC (Jiang et al., 2025).

Building on these applications, the current version extends EcoPlot-iso for long-term ecohydrological simulations and management scenario analysis, integrating a new depth-dependent root-uptake module (Jiang et al., 2025).

2. Model Framework and Structure

EcoPlot-iso is a process-based conceptual model that simulates:

Hydrological fluxes: interception, throughfall, infiltration, preferential flow, surface runoff, percolation, groundwater recharge.

Evapotranspiration components: canopy evaporation, soil evaporation, and transpiration.

Isotopic processes: fractionation and mixing in canopy and soil layers.

The vertical structure consists of:

** One canopy layer, and

Three soil layers:

** 0–10 cm (shallow)

** 10–30 cm (middle)

** 30–100 cm (deep)

Each layer tracks both water and stable isotope (δ²H or δ¹⁸O) balances.
The isotope module allows explicit separation of evaporation and transpiration losses, improving estimates of water partitioning and residence times.

 <img width="940" height="456" alt="image" src="https://github.com/user-attachments/assets/4e7b24a9-ea3c-473f-9dd0-66702caca42b" />

Figure 2. (a) Schematic representation of the ecohydrological fluxes and water partitioning in the EcoPlot-iso model illustrating major water fluxes and storage components; (b) Conceptual framework and key parameters of the EcoPlot-iso model(Landgraf et al., 2023; Stevenson et al., 2023), highlighting the key ecohydrological processes simulated in this study.


3. Required Input Data
Variable	Unit	Description
Precipitation (P)	mm d⁻¹	Amount of rainfall or snowfall
δ²H / δ¹⁸O of precipitation	‰ VSMOW	Isotopic composition of input water
Air temperature (Tair)	°C	Mean air temperature
Potential evapotranspiration (PET)	mm d⁻¹	Computed e.g. by FAO-56 Penman–Monteith
Relative humidity (RH)	%	Daily mean humidity
Leaf area index (LAI)	m² m⁻²	Vegetation canopy variable
Soil hydraulic parameters	–	Porosity, field capacity, wilting point, Ks, etc.
Vegetation parameters	–	Interception capacity, rooting depth parameter β, rL₁–rL₃, g₁–g₃ (stomatal conductance)


4. Model Outputs
Output	Unit	Description
Ei, Es, Tr	mm d⁻¹	Canopy, soil, and transpiration fluxes
Qs, Recharge	mm d⁻¹	Surface runoff and groundwater recharge
dI, dSTO, dGW, dSdeep	mm	Storage changes
δ²H / δ¹⁸O	‰	Isotopic composition of fluxes and storages
ET/P, Tr/ET	–	Aggregated water-partitioning indices

Outputs are stored in CSV or NetCDF format for analysis and visualization.

5. References

** Stevenson, J. L., Birkel, C., Comte, J. C., Tetzlaff, D., Marx, C., Neill, A., Maneta, M., Boll, J., & Soulsby, C. (2023). Quantifying heterogeneity in ecohydrological partitioning in urban green spaces through the integration of empirical and modelling approaches. Environmental Monitoring and Assessment, 195(4). https://doi.org/10.1007/s10661-023-11055-6

** Birkel, C., Tetzlaff, D., Ring, A. M., & Soulsby, C. (2025). Does high-resolution in-situ xylem and atmospheric vapor isotope data help improve modeled estimates of ecohydrological partitioning? Agricultural and Forest Meteorology, 365, 110467. https://doi.org/10.1016/j.agrformet.2025.110467

** Landgraf, J., Tetzlaff, D., Birkel, C., Stevenson, J. L., & Soulsby, C. (2023). Assessing land-use effects on ecohydrological partitioning in the critical zone through isotope-aided modelling. Earth Surface Processes and Landforms, 48(15), 3199–3219. https://doi.org/10.1002/esp.5691

** Jiang, C., Tetzlaff, D., Wu, S., Birkel, C., Laudon, H., & Soulsby, C. (2025). Assessing the drought resilience of different land-management scenarios using a tracer-aided ecohydrological model with variable root-uptake distributions. EGUsphere, 2025, 1–34. https://doi.org/10.5194/egusphere-2025-2533 (under discussion).

** Jiang, C., Soulsby, C., Laudon, H., Wu, S., & Tetzlaff, D. (2025). Predicting summer droughts in Central Europe from winter NAO. Nature Water (submitted).
