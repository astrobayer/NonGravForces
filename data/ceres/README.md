# data/ceres — CERES EBAF TOA monthly fluxes (excerpt, test oracle)

Two months (December 2019 and January 2020, the two nodes around the 2020-01-01
GRACE-FO day) of the CERES Energy Balanced and Filled TOA product, Edition 4.2.1,
as published for obs4MIPs:

| File | Variable | Source file (ESGF, dataset version v20260220) | SHA-256 of the source file |
|---|---|---|---|
| `rsut_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc` | TOA outgoing shortwave, W m-2 | `rsut_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc` | `6869b2207d86ec6581a6518f63b65f41fc209dda4ea5e64afbf6f3b6d47ea592` |
| `rlut_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc` | TOA outgoing longwave, W m-2 | `rlut_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc` | `1e611d51416d512a945ee4a61983f5c88a991f08a6aee6d6d7662ea76bd9ab5a` |
| `rsdt_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc` | TOA incident shortwave, W m-2 | `rsdt_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc` | `3e4cd4ad8f1c283369cd2dd385ac3f655c4e10c9ed2a6b954b829097a61125e8` |

Source: `https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/obs4MIPs/NASA-LaRC/CERES-EBAF-4-2-1/mon/<variable>/gn/v20260220/`,
downloaded anonymously on 2026-09-15 and checked against the SHA-256 ESGF publishes.

**How the excerpt was made.** Time indices 238 and 239 of each source file (mid-month
nodes 2019-12-15 06:00 and 2020-01-15 00:00 UTC) were read with MATLAB's `ncread` and
written unchanged (single precision) to netCDF-3 classic files, together with the
axes, their bounds and the attributes `io.importCERES` reads (`variable_id`,
`source_id`, `source_version_number`, `tracking_id`, `doi`, `license`, `units`, ...).
The global attributes `excerpt_of` and `excerpt_of_sha256` name the source. Classic
format was chosen so that `tools/pyxcheck/ceres_check.py` can read the same bytes with
a standard-library parser. Reading back gave the source values bitwise.

The full files are not committed. `setupNonGravForces` fetches them into the data pool,
`<DATAPOOL>/CERES/EBAF/`, when `'FetchData'` allows it (into this folder until 2026-11-08;
downloaded data no longer lands in the toolbox).

**Licence and acknowledgement.** Data in these files were produced by NASA Langley
Research Center and are licensed under the Creative Commons Attribution 4.0
International License (https://creativecommons.org/licenses/by/4.0/); the licence text
is the files' own `license` attribute. Acknowledge their use following
https://ceres.larc.nasa.gov/ — CERES EBAF-TOA Edition 4.2.1,
doi:10.5067/TERRA-AQUA-NOAA20/CERES/EBAF_L3B004.2.1; Loeb et al. (2018), J. Climate 31, 895–918.
The excerpt changes nothing but the time extent and the file format.

**What the grid means** (EBAF Ed4.2 Data Quality Summary): fluxes at a 20 km TOA
reference level; geodetic latitudes (the product's global means use geodetic zone
weights); monthly means with no diurnal cycle.
