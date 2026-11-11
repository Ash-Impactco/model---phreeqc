# Ilha Carbon ERW data pack

One file per dataset. Every row has a **status**:
- **measured**: a real lab or supplier value
- **stand-in**: measured, but at a different site or on a different rock (Marques Britas, Sanguinho)
- **assumed**: a modelling choice or literature value
- **TO MEASURE**: an empty row waiting for TECHNOVIA or Lomba da Maia results

The same tables are in `Ilha_Carbon_ERW_Data.xlsx`, one sheet each.

## Rock
| File | Contents |
|---|---|
| 01_rock_major_oxides.csv | Bulk chemistry: Marques Britas Brita 2 (model feedstock), olivine block, polished basalt; empty TECHNOVIA rows |
| 02_rock_physical_properties.csv | Density, porosity, absorption, strength, grading, chloride, ASR; moisture and bulk density to measure |
| 03_rock_grain_size.csv | Measured point (11% < 63 µm), assumed split used in the model, row for laser PSD |
| 04_rock_sem_eds_spots.csv | 3 SEM-EDS spots (semi-quantitative; spot 2 hit a P- and Ti-rich grain; Cr detected) |
| 05_rock_mineralogy.csv | Normative mineral estimate used by the model; glass unknown; row for XRD-Rietveld |
| 06_rock_trace_metals.csv | Isometric-required trace metals, radioactivity, leaching and deportment, all to measure |
| 07_rock_co2_potential.csv | tCO2 per tonne and per hectare from oxides (Renforth method) |

## Soil
| File | Contents |
|---|---|
| 08_soil_sanguinho_existing.csv | All 9 Sanguinho lab results (Jan/Oct 2024); Moreira and Terra da Riba are model soils A and B |
| 09_soil_lomba_da_maia_screening.csv | Template for the 10 composites + 1 duplicate (+ optional forest) on Luis's land, with GPS columns |
| 10_soil_methods_and_requests.csv | Which soil tests to request, why, method and priority |

## Site, climate and model
| File | Contents |
|---|---|
| 11_site_and_logistics.csv | Area, control, rate, tonnage, quarry and field coordinates, trucks, sampling plan |
| 12_climate_and_hydrology.csv | Temperature, rainfall, ET, drainage, soil CO2, water content (several still to source) |
| 13_model_rate_parameters.csv | Kinetic constants per mineral with sources |
| 14_model_settings.csv | Column, roughness, field factor and other model choices |
| 15_model_results_10yr.csv | Latest scenario results (10 years, stand-in data) |
| 16_references.csv | Sources and links |
| 17_isometric_accounting_terms.csv | Each Isometric EWA v1.3 term (Eq. 1, 2, 4), how the model estimates it, and what must be measured for credits |
| 18_emission_factors.csv | Placeholder emission factors for haulage, spreading, milling (replace with real data) |

## Points to check
- **Moreira and Terra da Riba potassium.** Two Notion tables disagree (K2O 1019 vs 286, and 806 vs 101 mg/kg). Check the original lab reports.
- **Mn oxide** is written as Mn2O4 in one Marques Britas source. It is probably Mn3O4.
- **The Sanguinho reports measured exchange acidity at pH 7** (for example 8.23 and 15.68 cmol/kg). That supports the model's "CEC at pH 7" acidity pool. **Ask for KCl exchangeable acidity at Lomba da Maia** as well, so both versions can be tested.
- **Nothing here is from TECHNOVIA or Lomba da Maia yet.** Fill files 01–06 and 09 as results arrive, then copy the values into the model's `data/inputs.json`, or let the GitHub prompt build the loader.
