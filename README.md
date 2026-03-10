# InSAR-based Landslide Detection — Intag Zone, Cotacachi, Ecuador (December 2023)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![Sentinel-1](https://img.shields.io/badge/Data-Sentinel--1-green.svg)](https://sentinel.esa.int/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18930871.svg)](https://doi.org/10.5281/zenodo.18930871)

## Overview

This repository contains the complete Python and Google Earth Engine processing pipeline used in:

> Pinto Páez, V. (under review). *InSAR-based detection of landslide-induced surface displacement in tropical cloud forest: the December 2023 Intag event, Cotacachi, Ecuador*. Earth Sciences Research Journal (ESRJ), Universidad Nacional de Colombia.

The scripts reproduce all analyses presented in the article, from Sentinel-1 data search to the final publication figures. The pipeline detects and characterizes landslide-induced surface displacement using Differential InSAR (DInSAR) and SBAS time series techniques in a tropical cloud forest environment (84.6% forest cover, Imbabura province, Ecuador).

---

## Study Area

- **Location:** Intag zone, Cotacachi canton, Imbabura province, Ecuador
- **AOI mask:** INTAJ.gpkg — 1,528 km², 6 parishes
- **Event date:** December 19, 2023
- **Affected population:** 5,029 people across 6 parishes (SGR, 2023)
- **Infrastructure damage:** 40 road segments, 1,325 m of road

---

## Key Results

| Metric | Value |
|--------|-------|
| InSAR deformation clusters detected | 1,481 |
| High-confidence clusters (4/4 tests) | 151 |
| Maximum LOS displacement (co-event) | −15.18 cm (P5) |
| Mean LOS displacement (co-event) | −7.13 cm |
| Landslide susceptibility model AUC | 0.670 |
| C-band coherence ≥ 0.3 (forest) | 6.5% of pixels |
| Clusters in High + Very High susceptibility zones | 64.9% |
| Spatial concentration ratio | 17.5× |

---

## Pipeline Structure

| Script | Description |
|--------|-------------|
| `00 CONFIGURACIÓN...` | Central project configuration |
| `P00` | Environment verification |
| `P01` | Environment check |
| `P02` | Sentinel-1 SLC search (from GPKG bounding box) |
| `P03` | Strategic DInSAR pair selection |
| `P04` | Submission of 8 pairs to HyP3 (INSAR_GAMMA) |
| `P04b_v2` | Path 18 selection, submission and download |
| `P05 v2.0` | InSAR product extraction, clipping and analysis |
| `P05b v2.0` | Ascending path fusion (Path 120 + Path 18) |
| `P05c` | File renaming for P06 compatibility |
| `P06 v2.0` | Deformation zone extraction (1,481 clusters) |
| `P06b` | Scientific cluster verification (4-test protocol) |
| `P07` | Optical validation with Sentinel-2 (GEE — JavaScript) |
| `P07b` | Cross-validation InSAR vs. optical (Sentinel-2) |
| `P07c` | Independent optical landslide inventory (Phase 1) |
| `P07c_auxiliar / auxiliar2` | GEE auxiliary layer export (JavaScript) |
| `P07d` | Spatial cross-reference InSAR and optical inventories |
| `P08` | Sentinel-1 coverage verification for SBAS Phase 3 |
| `P09` | SBAS stack design and submission to HyP3 |
| `P09b` | Diagnostic: submitted vs. pending SBAS pairs |
| `P09c` | Submission of 372 pending SBAS pairs (account 2) |
| `P09d` | SBAS processing monitoring (both accounts) |
| `P09e` | Download completed SBAS products (account 1) |
| `P09f` | Download missing SBAS products (account 2) |
| `P10` | SBAS product extraction and preparation for MintPy |
| `P11` | NISAR data availability verification for Intag |
| `P12` | NISAR L-band vs. Sentinel-1 C-band coherence analysis |
| `P12b` | C-band vs. L-band coherence comparative figures |
| `P13` | DEM morphometric variables (Phase 5, Step 5.1) |
| `P14` | Precipitation as landslide trigger (Phase 5, Step 5.2) |
| `P15` | Land cover analysis (Phase 5, Step 5.3) |
| `P15b` | Soil moisture analysis (Phase 5, Step 5.3) |
| `P16` | Land cover and use (Phase 5, Step 5.4) |
| `P17` | Landslide susceptibility model — Random Forest (AUC=0.670) |
| `P18` | Clipping to real AOI and publication figures |
| `P19` | Regeneration of all figures with corrected AOI |
| `P20c` | Susceptibility map over OpenStreetMap |
| `P21` | Final publication figures for ESRJ (7 figures, English) |

---

## Data Sources

| Dataset | Source | Resolution |
|---------|--------|------------|
| Sentinel-1 SLC | NASA ASF / ESA (Path 40 Descending) | 5×20 m |
| Sentinel-2 optical | ESA / Google Earth Engine | 10 m |
| Copernicus DEM | ESA | 30 m |
| CHIRPS precipitation | UCSB Climate Hazards Group | 5 km |
| WorldPop population | WorldPop (Tatem, 2017) | 100 m |
| Land cover | ESA WorldCover 2021 | 10 m |
| AOI mask | INTAJ.gpkg (this study) | 1:50,000 |

---

## Requirements

```
Python >= 3.8
hyp3_sdk
asf_search
geopandas
rasterio
numpy
pandas
scikit-learn
matplotlib
shapely
pyproj
mintpy (for SBAS time series)
```

Install dependencies:
```bash
pip install hyp3_sdk asf_search geopandas rasterio numpy pandas scikit-learn matplotlib shapely pyproj
```

---

## Credentials Setup

Scripts that interact with NASA Earthdata (HyP3, ASF) require credentials. 
**Never hardcode credentials in scripts.** Configure them via environment variables:

```bash
# Windows CMD
set EARTHDATA_USER=your_username
set EARTHDATA_PASS=your_password
```

Or configure a `_netrc` file:
```
machine urs.earthdata.nasa.gov
    login your_username
    password your_password
```

---

## Reproducibility

All scripts print their name, version, execution date, and MD5 hash of the primary input at runtime. The complete processing chain is documented in the Supplementary Methods of the article.

**AOI mask:** Always use `INTAJ.gpkg` (1,528 km², 6 parishes) — NOT the bounding box (4,263 km²).

---

## Author

**Víctor Pinto Páez**  
Independent researcher  
ORCID: [0009-0001-5573-8294](https://orcid.org/0009-0001-5573-8294)  
Email: vpintopaez@hotmail.com  
Location: Ibarra, Ecuador

---

## Citation

If you use these scripts, please cite:

> Pinto Páez, V. (under review). InSAR-based detection of landslide-induced surface displacement in tropical cloud forest: the December 2023 Intag event, Cotacachi, Ecuador. *Earth Sciences Research Journal*. *(DOI to be assigned upon acceptance)*

**Scripts repository:**
> Pinto Páez, V. (2026). *InSAR analysis scripts for Intag landslide detection, Ecuador (v1.0.0)*. Zenodo. https://doi.org/10.5281/zenodo.18930871

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

- NASA ASF HyP3 for InSAR processing
- ESA Copernicus programme for Sentinel-1 and Sentinel-2 data
- Secretaría de Gestión de Riesgos del Ecuador (SGR) for the damage inventory (Report No. 0729)
- INEC Ecuador for census data
- OpenStreetMap contributors (ODbL license)
