**1. Overview**

EcoPlot-iso is a tracer-aided ecohydrological modeling framework designed to simulate key ecohydrological and isotopic transformations that characterize water partitioning at the plot scale.
It integrates hydrological, isotopic, and vegetation processes within the Soil–Plant–Atmosphere Continuum (SPAC) to quantify water fluxes, storage changes, and isotopic dynamics under different climatic and land-use conditions.

The model has been applied in diverse climatic and hydrological settings, including:

** a one-year simulation in Scotland (Stevenson et al., 2023),

** a one-year simulation at the Demnitzer Mill Creek (DMC) catchment in Germany (Landgraf et al., 2023),

** a four-year tropical application in Costa Rica (Birkel et al., 2024), and

** a long-term (2000–2024) tracer-aided simulation for drought resilience assessment in DMC (Jiang et al., 2025).

Building on these applications, the current version extends EcoPlot-iso for long-term ecohydrological simulations and management scenario analysis, integrating a new depth-dependent root-uptake module (Jiang et al., 2025).

**2. Model Framework and Structure**

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

Figure 2. (a) Schematic representation of the ecohydrological fluxes and water partitioning in the EcoPlot-iso model illustrating major water fluxes and storage components; (b) Conceptual framework and key parameters of the EcoPlot-iso model(Landgraf et al., 2023; Stevenson et al., 2023), highlighting the key ecohydrological processes simulated in this study (Jiang et al., 2025).


**3. Required Input Data**
   
| Variable                          | Unit       | Description                                                           |
|-----------------------------------|------------|-----------------------------------------------------------------------|
| Precipitation (P)                 | mm d⁻¹     | Amount of rainfall or snowfall                                        |
| δ²H of precipitation              | ‰ VSMOW    | Isotopic composition of input water                                   |
| Air temperature (Tair)            | °C         | Mean air temperature                                                  |
| Potential evapotranspiration (PET)| mm d⁻¹     | Computed by FAO-56 Penman–Monteith                                    |
| Relative humidity (RH)            | %          | Daily humidity                                                   |
| Leaf area index (LAI)             | m² m⁻²     | Vegetation canopy variable                                            |


**4. Model Outputs**

| Output Variable              | Unit      | Description                                                                          |
|------------------------------|-----------|--------------------------------------------------------------------------------------|
| Canopy Evaporation (Ei)      | mm/d      | Water evaporated from canopy interception                                            |
| Soil Evaporation (Es)        | mm/d      | Direct evaporation from soil surface                                                 |
| Transpiration (Tr)           | mm/d      | Water lost via plant transpiration                                                   |
| Interception Storage (I)     | mm        | Water held on canopy before evaporation or throughfall                               |
| Throughfall (Tf)             | mm/d      | Rainfall passing through canopy to the soil                                          |
| Infiltration (Inf)           | mm/d      | Water entering top soil layer after rainfall                                         |
| Preferential Flow (PF)       | mm/d      | Fast-flow water bypassing matrix via macropores or cracks                            |
| Surface Runoff (Qs)          | mm/d      | Overland flow generated at soil surface                                              |
| Percolation (Perc)           | mm/d      | Downward movement of water between soil layers                                       |
| Groundwater Recharge         | mm/d      | Water flux from deepest soil layer into groundwater                                  |
| Soil Moisture (0–10 cm)      | mm        | Water stored in shallow soil layer                                                   |
| Soil Moisture (10–30 cm)     | mm        | Water stored in middle soil layer                                                    |
| Soil Moisture (30–100 cm)    | mm        | Water stored in deep soil layer                                                      |
| Isotopic Composition (δ²H)   | ‰         | Stable isotope ratios of water in fluxes and storages (for each process/layer)       |

**5. Run the model**
   
   (a) Install required R packages from Conda-forge
   
   conda install -c conda-forge r-data.table r-fme r-tidyverse r-gridextra r-corrplot r-gplots r-rcolorbrewer r-factoextra r-ggplot2 r-ggpubr r-ggsci r-scales r-lubridate r-cowplot r-hydrogof
   
   (b) Execute the EcoPlot-iso R script for different sites
   
    Rscript Script_SWBiso_Forest.R

   (c) Execute the R script for the forest management scenarios

   Rscript Script_SWBiso_Scenario.R

   
**6. References**

** Stevenson, J. L., Birkel, C., Comte, J. C., Tetzlaff, D., Marx, C., Neill, A., Maneta, M., Boll, J., & Soulsby, C. (2023). Quantifying heterogeneity in ecohydrological partitioning in urban green spaces through the integration of empirical and modelling approaches. Environmental Monitoring and Assessment, 195(4). https://doi.org/10.1007/s10661-023-11055-6

** Birkel, C., Arciniega-Esparza, S., Maneta, M. P., Boll, J., Stevenson, J. L., Benegas-Negri, L., Tetzlaff, D., & Soulsby, C. (2024). Importance of measured transpiration fluxes for modelled ecohydrological partitioning in a tropical agroforestry system.
Agricultural and Forest Meteorology, 346, 109870. https://doi.org/10.1016/j.agrformet.2023.109870

** Landgraf, J., Tetzlaff, D., Birkel, C., Stevenson, J. L., & Soulsby, C. (2023). Assessing land-use effects on ecohydrological partitioning in the critical zone through isotope-aided modelling. Earth Surface Processes and Landforms, 48(15), 3199–3219. https://doi.org/10.1002/esp.5691

** Jiang, C., Tetzlaff, D., Wu, S., Birkel, C., Laudon, H., & Soulsby, C. (2025). Assessing the drought resilience of different land-management scenarios using a tracer-aided ecohydrological model with variable root-uptake distributions. EGUsphere, 2025, 1–34. https://doi.org/10.5194/egusphere-2025-2533 (under discussion).

** Jiang, C., Soulsby, C., Laudon, H., Wu, S., & Tetzlaff, D. (2025). Predicting summer droughts in Central Europe from winter NAO (submitted).
