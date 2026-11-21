# Peatland and lagged drought effects on fire activity in Kalimantan

R code for the paper *"Peatland amplifies and prolongs lagged drought effects on tropical fire activity: a Bayesian distributed lag non-linear model"* (F. Ardiansyah, submitted for publication).

## Data (all publicly available, downloaded by the scripts or manually)

| Data | Source |
|---|---|
| VIIRS S-NPP 375 m active fires, 2012–2025 | NASA FIRMS, https://firms.modaps.eosdis.nasa.gov |
| CHIRPS v3.0 monthly rainfall | Climate Hazards Center; Kalimantan subset from the IRI Data Library (`data.nc`, see script 01) |
| District boundaries | GADM 4.1, https://gadm.org |
| Peatland | PEATMAP (Xu et al. 2018), https://archive.researchdata.leeds.ac.uk/251/ (`SEA_Peatland.shp`) |
| Niño 3.4 and DMI indices | NOAA PSL |

## Software

R 4.5, R-INLA 25.10.19, and the CRAN packages `dlnm`, `sf`, `terra`, `geodata`, `ncdf4`, `spdep`, `dplyr`, `tidyr`, `readr`, `lubridate`, `ggplot2`, `MASS`, `sn`.

## Scripts (run in this order)

| Script | Content |
|---|---|
| `01_siapkan_data.R` | Builds the district × month panel: fire counts, rainfall, climate indices, adjacency graph |
| `02_eksplorasi_dlnm.R` | Exploratory analysis, cross-basis, negative binomial DLNM without random effects (M0, M1) |
| `03a_inla_spasial.R` | M2: DLNM + BYM2 (INLA) |
| `03b_inla_musiman.R` | M2, M3a, M3b: adding seasonal and yearly effects; WAIC comparison |
| `04_gambut.R` | Peat proportion per district from PEATMAP |
| `05_inla_gambut.R` | M4: cross-basis × peat interaction |
| `06_diagnostik_validasi.R` | Residual diagnostics and hindcasts of 2015, 2019 and 2023 |
| `07_uji_ketahanan.R` | M5a and M5b: regional and district-specific AR1 terms |
| `08_sensitivitas.R` | Sensitivity of M5b to lag length and spline specification |
| `09_simulasi.R` | Simulation study (30 replicates per scenario; resumable) |
| `10b_gambar_inggris_cepat.R` | Figures for the paper (from saved results, no model refitting) |

Each script writes its results to a folder `hasil_XX/`. Comments in the scripts are in Indonesian.

## Model naming

M2 = DLNM + BYM2; M3a = M2 + season; M3b = M3a + year; M4 = M3b + peat interaction; M5a = M4 + regional AR1; M5b = M4 + district-specific AR1 (final model).
