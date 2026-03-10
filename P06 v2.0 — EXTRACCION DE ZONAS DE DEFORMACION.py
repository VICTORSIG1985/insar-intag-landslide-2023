# =============================================================================
# P06 v2.0 — EXTRACCION DE ZONAS DE DEFORMACION
# =============================================================================
# Proyecto: InSAR Deslizamientos Zona de Intag, Cotacachi, Ecuador
# Autor:    Víctor Pinto
# Fecha:    2026-02-27
# Version:  2.0
#
# PROPOSITO:
#   Identificar y extraer zonas de deformacion significativa a partir de
#   los productos InSAR procesados por P05 v2 y P05b v2. Aplica:
#     1. Enmascaramiento por coherencia (gamma > umbral)
#     2. Filtrado por desplazamiento (|d| > umbral cm)
#     3. Clustering espacial de pixeles deformados (connected components)
#     4. Extraccion de centroides y poligonos de cada cluster
#     5. Estadisticas por cluster (area, desplazamiento medio/max, coherencia)
#     6. Exportacion a GeoPackage + CSV + visualizaciones
#
# CAMBIOS v1.x -> v2.0:
#   - Clasificacion temporal corregida (P05_FIX): 10 CO, 1 PRE, 1 POST
#   - Mosaicos corregidos: solo M_ASC_CO_12d y M_ASC_CO_24d
#   - Conversion explicita m->cm con verificacion (FACTOR = 100)
#   - Datasets reorganizados: DESC(4) + ASC_MOSAIC(2) + ASC_INDIV(4)
#
# ENTRADA:
#   hyp3_clipped/   -> rasters individuales recortados al AOI
#   hyp3_mosaic/    -> mosaicos CRB ascendentes
#
# SALIDA:
#   hyp3_deformation/
#     clusters_{id}.gpkg, centroids_{id}.gpkg, summary_{id}.csv,
#     map_{id}.png, P06_v2_resumen_completo.json, P06_v2_catalogo_clusters.csv
#
# EJECUTAR EN: Spyder (Run file F5)
# Ref: Protocolo Maestro InSAR Intag v2.0, Fase 2, Paso 2.7
#      Dini et al. (2026), coherence-based landslide detection
# =============================================================================

import json
import os
import sys
import csv
import warnings
from datetime import datetime, date
from collections import Counter

import numpy as np

# --- Verificacion de librerias -----------------------------------------------
print("=" * 80)
print("P06 v2.0 — EXTRACCION DE ZONAS DE DEFORMACION")
print("=" * 80)
print(f"  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("\n--- Verificando librerias ---\n")

_requeridas = {
    "rasterio":    "rasterio",
    "geopandas":   "geopandas",
    "numpy":       "numpy",
    "matplotlib":  "matplotlib",
    "scipy":       "scipy",
    "shapely":     "shapely",
}
_faltantes = []
for nombre, modulo in _requeridas.items():
    try:
        mod = __import__(modulo)
        ver = getattr(mod, "__version__", "?")
        print(f"  ok {nombre:12s} v{ver}")
    except ImportError:
        _faltantes.append(nombre)
        print(f"  XX {nombre:12s} NO instalado")

if _faltantes:
    print(f"\n  [ERROR] Faltan: {', '.join(_faltantes)}")
    print(f"  Instalar: conda install -c conda-forge {' '.join(_faltantes)}")
    sys.exit(1)

import geopandas as gpd
import rasterio
from rasterio.features import shapes as rasterio_shapes
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.ndimage import label as scipy_label
from shapely.geometry import shape, Point

warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)
print("  ok Todas las librerias disponibles")

# =============================================================================
# 1. CONFIGURACION
# =============================================================================
print("\n--- CONFIGURACION ---")

# -- Rutas --------------------------------------------------------------------
SAVE_DIR      = r"D:\POSGRADOS\INTAG\data\sentinel1"
CLIPPED_DIR   = os.path.join(SAVE_DIR, "hyp3_clipped")
MOSAIC_DIR    = os.path.join(SAVE_DIR, "hyp3_mosaic")
DEFORM_DIR    = os.path.join(SAVE_DIR, "hyp3_deformation")
AOI_PATH      = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"

os.makedirs(DEFORM_DIR, exist_ok=True)

print(f"  Clipped:      {CLIPPED_DIR}")
print(f"  Mosaicos:     {MOSAIC_DIR}")
print(f"  Deformacion:  {DEFORM_DIR}")
print(f"  AOI:          {AOI_PATH}")

# -- Parametros de analisis ---------------------------------------------------
#
# COHERENCIA_MIN: 0.3
#   Umbral estandar en InSAR para filtrar pixeles con fase no confiable.
#   En bosque nublado tropical, la coherencia tipica es baja (0.15-0.35)
#   por la densa vegetacion. 0.3 es el minimo para estimaciones de
#   desplazamiento interpretables (Hanssen, 2001).
COHERENCIA_MIN = 0.3

# DESPLAZAMIENTO_MIN_CM: 3.0
#   Umbral minimo de desplazamiento LOS para considerar deformacion
#   significativa. Basado en la precision tipica de Sentinel-1 C-band
#   para interferogramas de 12-48 dias (~1-2 cm de ruido atmosferico).
#   3 cm excede el ruido con margen de seguridad (SNR > 1.5).
DESPLAZAMIENTO_MIN_CM = 3.0

# DESPLAZAMIENTO_FUERTE_CM: 5.0
#   Umbral para deformacion "fuerte" que indica movimiento claramente
#   asociado a deslizamiento activo. Usado para clasificar clusteres.
DESPLAZAMIENTO_FUERTE_CM = 5.0

# MIN_PIXELS_CLUSTER: 4
#   Minimo de pixeles conectados para formar un cluster valido.
#   Con resolucion HyP3 de ~80m (20x4 looks), 4 pixeles = ~25,600 m2
#   (2.56 ha), consistente con el tamano minimo de deslizamientos
#   detectables por Sentinel-1 (Dini et al., 2026).
MIN_PIXELS_CLUSTER = 4

# FACTOR_M_A_CM: 100
#   HyP3 entrega desplazamientos en METROS.
#   Los umbrales estan en CENTIMETROS.
#   CRITICO: Sin esta conversion, los umbrales de 3-5 cm nunca se alcanzan
#   porque los valores raw estan en el rango [-0.22, 0.05] m.
FACTOR_M_A_CM = 100

print(f"\n  Parametros:")
print(f"    Coherencia minima:       {COHERENCIA_MIN}")
print(f"    Desplazamiento minimo:   {DESPLAZAMIENTO_MIN_CM} cm")
print(f"    Desplazamiento fuerte:   {DESPLAZAMIENTO_FUERTE_CM} cm")
print(f"    Pixeles min. por cluster: {MIN_PIXELS_CLUSTER}")
print(f"    Factor m->cm:            x{FACTOR_M_A_CM}")

# -- Fecha del evento ---------------------------------------------------------
EVENTO = date(2023, 12, 19)

# =============================================================================
# 2. INVENTARIO DE DATASETS (clasificacion temporal corregida v2.0)
# =============================================================================
print("\n--- INVENTARIO DE DATASETS ---")


def clasificar_temporal(pre_str, post_str, evento=EVENTO):
    # PRE-EVENTO:  ambas fechas < evento
    # CO-EVENTO:   master < evento, secondary >= evento
    # POST-EVENTO: ambas fechas > evento
    pre  = date(int(pre_str[:4]), int(pre_str[4:6]), int(pre_str[6:8]))
    post = date(int(post_str[:4]), int(post_str[4:6]), int(post_str[6:8]))
    if pre < evento and post < evento:
        return "PRE-EVENTO"
    elif pre >= evento and post > evento:
        return "POST-EVENTO"
    else:
        return "CO-EVENTO"


# -- PARES MAESTRO CORREGIDO (fuente: P05_FIX / pares_maestro_corregido_v2.json)
PARES = {
    # Descendente Path 40 — 100% AOI
    "D1": {"track": "P40",  "dir": "DESC", "pre": "20231208", "post": "20231220", "bl": 12},
    "D2": {"track": "P40",  "dir": "DESC", "pre": "20231126", "post": "20231220", "bl": 24},
    "D3": {"track": "P40",  "dir": "DESC", "pre": "20231208", "post": "20240113", "bl": 36},
    "D4": {"track": "P40",  "dir": "DESC", "pre": "20231208", "post": "20240125", "bl": 48},
    # Ascendente Path 120 — 69.2% AOI
    "A1": {"track": "P120", "dir": "ASC",  "pre": "20231213", "post": "20231225", "bl": 12},
    "A2": {"track": "P120", "dir": "ASC",  "pre": "20231201", "post": "20231225", "bl": 24},
    "A3": {"track": "P120", "dir": "ASC",  "pre": "20231213", "post": "20240106", "bl": 24},
    "A4": {"track": "P120", "dir": "ASC",  "pre": "20231213", "post": "20240130", "bl": 48},
    # Ascendente Path 18 — 69.8% AOI
    "P18_A3": {"track": "P18", "dir": "ASC", "pre": "20231218", "post": "20240111", "bl": 24},
    "P18_B1": {"track": "P18", "dir": "ASC", "pre": "20231124", "post": "20231218", "bl": 24},
    "P18_B2": {"track": "P18", "dir": "ASC", "pre": "20231218", "post": "20231230", "bl": 12},
    "P18_B3": {"track": "P18", "dir": "ASC", "pre": "20231230", "post": "20240123", "bl": 24},
}

# Agregar clasificacion temporal
for pid, info in PARES.items():
    info["cat"] = clasificar_temporal(info["pre"], info["post"])

# -- Definicion de DATASETS para analisis -------------------------------------
DATASETS = []

# 2a. Descendentes individuales (100% AOI)
for pid in ["D1", "D2", "D3", "D4"]:
    p = PARES[pid]
    DATASETS.append({
        "id":       pid,
        "tipo":     "individual",
        "track":    p["track"],
        "dir":      p["dir"],
        "cat":      p["cat"],
        "baseline": p["bl"],
        "pre":      p["pre"],
        "post":     p["post"],
        "cobertura": "100%",
        "corr_path":      os.path.join(CLIPPED_DIR, pid, f"{pid}_corr.tif"),
        "los_disp_path":  os.path.join(CLIPPED_DIR, pid, f"{pid}_los_disp.tif"),
        "vert_disp_path": os.path.join(CLIPPED_DIR, pid, f"{pid}_vert_disp.tif"),
        "dem_path":       os.path.join(CLIPPED_DIR, pid, f"{pid}_dem.tif"),
    })

# 2b. Mosaicos ascendentes CRB (100% AOI combinado)
MOSAICOS = {
    "M_ASC_CO_12d": {
        "p120": "A1", "p18": "P18_B2",
        "cat": "CO-EVENTO", "baseline": 12,
        "label": "Co-evento corto (~12d)",
    },
    "M_ASC_CO_24d": {
        "p120": "A3", "p18": "P18_A3",
        "cat": "CO-EVENTO", "baseline": 24,
        "label": "Co-evento amplio (~24d)",
    },
}

for mid, minfo in MOSAICOS.items():
    mos_dir = os.path.join(MOSAIC_DIR, mid)
    DATASETS.append({
        "id":       mid,
        "tipo":     "mosaic",
        "track":    "P120+P18",
        "dir":      "ASC",
        "cat":      minfo["cat"],
        "baseline": minfo["baseline"],
        "pre":      "multi",
        "post":     "multi",
        "cobertura": "~100%",
        "label":    minfo["label"],
        "corr_path":      os.path.join(mos_dir, f"{mid}_corr.tif"),
        "los_disp_path":  os.path.join(mos_dir, f"{mid}_los_disp.tif"),
        "vert_disp_path": os.path.join(mos_dir, f"{mid}_vert_disp.tif"),
        "dem_path":       os.path.join(mos_dir, f"{mid}_dem.tif"),
    })

# 2c. Ascendentes individuales sin mosaico (69.2% o 69.8% AOI)
for pid in ["A2", "A4", "P18_B1", "P18_B3"]:
    p = PARES[pid]
    DATASETS.append({
        "id":       pid,
        "tipo":     "individual",
        "track":    p["track"],
        "dir":      p["dir"],
        "cat":      p["cat"],
        "baseline": p["bl"],
        "pre":      p["pre"],
        "post":     p["post"],
        "cobertura": "69.2%" if p["track"] == "P120" else "69.8%",
        "corr_path":      os.path.join(CLIPPED_DIR, pid, f"{pid}_corr.tif"),
        "los_disp_path":  os.path.join(CLIPPED_DIR, pid, f"{pid}_los_disp.tif"),
        "vert_disp_path": os.path.join(CLIPPED_DIR, pid, f"{pid}_vert_disp.tif"),
        "dem_path":       os.path.join(CLIPPED_DIR, pid, f"{pid}_dem.tif"),
    })

# Mostrar inventario
print(f"\n  Total datasets: {len(DATASETS)}")
print(f"  {'─'*78}")
print(f"  {'ID':<16} {'Tipo':<10} {'Track':<8} {'Cat':<14} {'BL':>3} {'Cob':>6}")
print(f"  {'─'*78}")
for ds in DATASETS:
    print(f"  {ds['id']:<16} {ds['tipo']:<10} {ds['track']:<8} "
          f"{ds['cat']:<14} {ds['baseline']:>3}d {ds['cobertura']:>6}")
print(f"  {'─'*78}")

# =============================================================================
# 3. FUNCIONES DE ANALISIS
# =============================================================================


def cargar_raster(path, banda=1):
    # Carga un raster GeoTIFF y retorna (data, transform, crs, nodata).
    # Convierte nodata a NaN para operaciones numericas.
    if not os.path.exists(path):
        return None, None, None, None
    with rasterio.open(path) as src:
        data = src.read(banda).astype(np.float64)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs
    if nodata is not None:
        data[data == nodata] = np.nan
    return data, transform, crs, nodata


def extraer_clusters(corr, los_cm, vert_cm, transform, crs,
                     coh_min, disp_min_cm, min_pix):
    # Pipeline de extraccion de clusteres de deformacion:
    #   1. Mascara de coherencia: gamma > coh_min
    #   2. Mascara de desplazamiento: |LOS| > disp_min_cm
    #   3. Mascara combinada: coherencia AND desplazamiento
    #   4. Connected component labeling (scipy.ndimage.label, 8-conectividad)
    #   5. Filtrado por tamano minimo (min_pix pixeles)
    #   6. Estadisticas por cluster
    #   7. Vectorizacion a poligonos y centroides
    #
    # Retorna: clusters_gdf, centroids_gdf, stats_list, mask_combined, labeled

    # Paso 1: Mascara de coherencia
    mask_coh = (corr > coh_min) & np.isfinite(corr)

    # Paso 2: Mascara de desplazamiento (valor absoluto: subsidencia + uplift)
    mask_disp = (np.abs(los_cm) > disp_min_cm) & np.isfinite(los_cm)

    # Paso 3: Combinacion
    mask_combined = mask_coh & mask_disp
    n_combined = int(np.sum(mask_combined))

    if n_combined == 0:
        return None, None, [], mask_combined, None

    # Paso 4: Connected components (8-conectividad)
    structure = np.ones((3, 3), dtype=int)
    labeled, n_features = scipy_label(mask_combined.astype(int), structure=structure)

    # Paso 5+6: Filtrar por tamano minimo + estadisticas
    cluster_stats = []

    for label_id in range(1, n_features + 1):
        mask_label = labeled == label_id
        n_pix = int(np.sum(mask_label))

        if n_pix < min_pix:
            labeled[mask_label] = 0
            continue

        los_vals = los_cm[mask_label]
        vert_vals = vert_cm[mask_label] if vert_cm is not None else None
        coh_vals = corr[mask_label]

        # Centroide (fila, columna promedio -> coordenadas)
        rows, cols = np.where(mask_label)
        mean_row = float(np.mean(rows))
        mean_col = float(np.mean(cols))
        cx, cy = rasterio.transform.xy(transform, mean_row, mean_col)

        # Area
        pixel_area_m2 = abs(transform.a * transform.e)
        area_m2 = n_pix * pixel_area_m2
        area_ha = area_m2 / 10000.0

        stats = {
            "cluster_id":    label_id,
            "n_pixels":      n_pix,
            "area_m2":       round(area_m2, 1),
            "area_ha":       round(area_ha, 3),
            "centroid_x":    round(cx, 6),
            "centroid_y":    round(cy, 6),
            "los_mean_cm":   round(float(np.nanmean(los_vals)), 3),
            "los_min_cm":    round(float(np.nanmin(los_vals)), 3),
            "los_max_cm":    round(float(np.nanmax(los_vals)), 3),
            "los_std_cm":    round(float(np.nanstd(los_vals)), 3),
            "coh_mean":      round(float(np.nanmean(coh_vals)), 3),
            "coh_min":       round(float(np.nanmin(coh_vals)), 3),
        }

        if vert_vals is not None and np.any(np.isfinite(vert_vals)):
            stats["vert_mean_cm"] = round(float(np.nanmean(vert_vals)), 3)
            stats["vert_min_cm"]  = round(float(np.nanmin(vert_vals)), 3)
            stats["vert_max_cm"]  = round(float(np.nanmax(vert_vals)), 3)

        # Clasificacion de intensidad
        max_abs_los = max(abs(stats["los_min_cm"]), abs(stats["los_max_cm"]))
        if max_abs_los >= DESPLAZAMIENTO_FUERTE_CM:
            stats["intensidad"] = "FUERTE"
        else:
            stats["intensidad"] = "MODERADO"

        cluster_stats.append(stats)

    if not cluster_stats:
        return None, None, [], mask_combined, labeled

    # Paso 7a: Vectorizacion — poligonos
    poly_geoms = []
    poly_labels = []

    for stats in cluster_stats:
        lid = stats["cluster_id"]
        mask_bin = (labeled == lid).astype(np.uint8)

        for geom_dict, value in rasterio_shapes(
            mask_bin, mask=(mask_bin == 1), transform=transform
        ):
            poly_geoms.append(shape(geom_dict))
            poly_labels.append(lid)

    if poly_geoms:
        clusters_gdf = gpd.GeoDataFrame(
            {"cluster_id": poly_labels},
            geometry=poly_geoms,
            crs=crs,
        )
        # Disolver por cluster_id (unir fragmentos del mismo cluster)
        clusters_gdf = clusters_gdf.dissolve(by="cluster_id").reset_index()

        # Agregar estadisticas
        stats_df = {s["cluster_id"]: s for s in cluster_stats}
        for col in ["n_pixels", "area_ha", "los_mean_cm", "los_min_cm",
                     "los_max_cm", "coh_mean", "intensidad"]:
            clusters_gdf[col] = clusters_gdf["cluster_id"].map(
                lambda lid, c=col: stats_df.get(lid, {}).get(c, None)
            )
    else:
        clusters_gdf = None

    # Paso 7b: Centroides
    centroid_geoms = [Point(s["centroid_x"], s["centroid_y"]) for s in cluster_stats]
    centroids_gdf = gpd.GeoDataFrame(
        cluster_stats,
        geometry=centroid_geoms,
        crs=crs,
    )

    return clusters_gdf, centroids_gdf, cluster_stats, mask_combined, labeled


def generar_mapa(los_cm, corr, mask, labeled, transform, crs,
                 dataset_id, cat, baseline, output_path):
    # Genera mapa de deformacion con dos paneles:
    #   Izquierdo: desplazamiento LOS enmascarado por coherencia
    #   Derecho: clusteres identificados sobre coherencia

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))

    # -- Panel 1: Desplazamiento LOS ------------------------------------------
    los_masked = np.where(corr > COHERENCIA_MIN, los_cm, np.nan)

    valid_vals = los_masked[np.isfinite(los_masked)]
    if len(valid_vals) > 0:
        vmax = max(abs(np.nanpercentile(valid_vals, 2)),
                   abs(np.nanpercentile(valid_vals, 98)))
    else:
        vmax = 10

    norm1 = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    extent = rasterio.transform.array_bounds(
        los_cm.shape[0], los_cm.shape[1], transform)
    extent_plot = [extent[0], extent[2], extent[1], extent[3]]

    im1 = ax1.imshow(los_masked, cmap="RdBu", norm=norm1,
                     extent=extent_plot, origin="upper")
    plt.colorbar(im1, ax=ax1, label="Desplazamiento LOS (cm)", shrink=0.7)
    ax1.set_title(f"{dataset_id} — Desplazamiento LOS\n"
                  f"[{cat}, {baseline}d, gamma>{COHERENCIA_MIN}]",
                  fontsize=11, fontweight="bold")
    ax1.set_xlabel("Longitud")
    ax1.set_ylabel("Latitud")

    # -- Panel 2: Clusteres ----------------------------------------------------
    ax2.imshow(corr, cmap="gray", vmin=0, vmax=1,
               extent=extent_plot, origin="upper", alpha=0.5)

    if labeled is not None and np.any(labeled > 0):
        cluster_display = np.where(labeled > 0, labeled, np.nan)
        ax2.imshow(cluster_display, cmap="tab20", interpolation="nearest",
                   extent=extent_plot, origin="upper", alpha=0.8)

    n_clusters = len(set(labeled[labeled > 0])) if labeled is not None else 0
    ax2.set_title(f"{dataset_id} — {n_clusters} clusteres detectados\n"
                  f"[|LOS|>{DESPLAZAMIENTO_MIN_CM}cm, >={MIN_PIXELS_CLUSTER}px]",
                  fontsize=11, fontweight="bold")
    ax2.set_xlabel("Longitud")
    ax2.set_ylabel("Latitud")

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


# =============================================================================
# 4. PROCESAMIENTO DE TODOS LOS DATASETS
# =============================================================================
print(f"\n{'='*80}")
print("PROCESAMIENTO DE DATASETS")
print(f"{'='*80}")

resultados_globales = []
catalogo_clusters = []
figuras_generadas = []

for i, ds in enumerate(DATASETS, 1):
    ds_id = ds["id"]
    print(f"\n{'─'*80}")
    print(f"  [{i}/{len(DATASETS)}] {ds_id}  ({ds['cat']}, {ds['baseline']}d, "
          f"{ds['track']} {ds['dir']})")
    print(f"{'─'*80}")

    # -- 4a. Verificar existencia de archivos ---------------------------------
    archivos_ok = True
    for key in ["corr_path", "los_disp_path"]:
        if not os.path.exists(ds[key]):
            print(f"  XX No encontrado: {ds[key]}")
            archivos_ok = False

    if not archivos_ok:
        print(f"  !! Saltando {ds_id} — archivos faltantes")
        resultados_globales.append({
            "dataset_id": ds_id, "status": "ARCHIVOS_FALTANTES",
            "cat": ds["cat"], "baseline": ds["baseline"],
        })
        continue

    # -- 4b. Cargar rasters ---------------------------------------------------
    print(f"  Cargando rasters...")
    corr, transform, crs, _ = cargar_raster(ds["corr_path"])
    los_raw, _, _, _ = cargar_raster(ds["los_disp_path"])

    vert_raw = None
    if os.path.exists(ds.get("vert_disp_path", "")):
        vert_raw, _, _, _ = cargar_raster(ds["vert_disp_path"])

    if corr is None or los_raw is None:
        print(f"  XX Error al cargar rasters")
        resultados_globales.append({
            "dataset_id": ds_id, "status": "ERROR_LECTURA",
        })
        continue

    # -- 4c. CONVERSION METROS -> CENTIMETROS ---------------------------------
    # CRITICO: HyP3 entrega desplazamientos en METROS
    # Los umbrales de analisis estan en CENTIMETROS
    los_cm = los_raw * FACTOR_M_A_CM
    vert_cm = vert_raw * FACTOR_M_A_CM if vert_raw is not None else None

    # Verificacion de conversion
    los_finite = los_raw[np.isfinite(los_raw)]
    if len(los_finite) > 0:
        print(f"  -- Verificacion de unidades --")
        print(f"     Raw (metros):   min={np.nanmin(los_finite):.4f}  "
              f"max={np.nanmax(los_finite):.4f}  mean={np.nanmean(los_finite):.4f}")
        print(f"     x{FACTOR_M_A_CM} (cm):      "
              f"min={np.nanmin(los_finite)*FACTOR_M_A_CM:.2f}  "
              f"max={np.nanmax(los_finite)*FACTOR_M_A_CM:.2f}  "
              f"mean={np.nanmean(los_finite)*FACTOR_M_A_CM:.2f}")

    # -- 4d. Estadisticas generales -------------------------------------------
    corr_finite = corr[np.isfinite(corr) & (corr > 0)]
    los_cm_finite = los_cm[np.isfinite(los_cm)]

    coh_mean = float(np.nanmean(corr_finite)) if len(corr_finite) > 0 else 0
    coh_gt03 = float(np.sum(corr_finite > 0.3) / len(corr_finite) * 100) \
               if len(corr_finite) > 0 else 0
    los_mean = float(np.nanmean(los_cm_finite)) if len(los_cm_finite) > 0 else 0
    los_p5   = float(np.nanpercentile(los_cm_finite, 5)) if len(los_cm_finite) > 0 else 0
    los_p95  = float(np.nanpercentile(los_cm_finite, 95)) if len(los_cm_finite) > 0 else 0

    print(f"\n  Estadisticas generales:")
    print(f"    Coherencia: mu={coh_mean:.3f}, >{COHERENCIA_MIN}={coh_gt03:.1f}%")
    print(f"    LOS (cm):   mu={los_mean:.2f}, P5={los_p5:.2f}, P95={los_p95:.2f}")

    # -- 4e. Extraccion de clusteres ------------------------------------------
    print(f"\n  Extrayendo clusteres (gamma>{COHERENCIA_MIN}, "
          f"|LOS|>{DESPLAZAMIENTO_MIN_CM}cm)...")

    clusters_gdf, centroids_gdf, cluster_stats, mask, labeled = extraer_clusters(
        corr, los_cm, vert_cm, transform, crs,
        COHERENCIA_MIN, DESPLAZAMIENTO_MIN_CM, MIN_PIXELS_CLUSTER
    )

    n_clusters = len(cluster_stats)
    n_fuerte = sum(1 for s in cluster_stats if s["intensidad"] == "FUERTE")
    n_moderado = n_clusters - n_fuerte

    print(f"    Clusteres encontrados: {n_clusters}")
    print(f"      FUERTE (|LOS|>{DESPLAZAMIENTO_FUERTE_CM}cm): {n_fuerte}")
    print(f"      MODERADO ({DESPLAZAMIENTO_MIN_CM}-"
          f"{DESPLAZAMIENTO_FUERTE_CM}cm): {n_moderado}")

    # -- 4f. Mostrar detalle de clusteres -------------------------------------
    if n_clusters > 0:
        print(f"\n    {'CL':>4} {'PIXELS':>6} {'AREA_ha':>8} "
              f"{'LOS_mu':>8} {'LOS_min':>8} {'LOS_max':>8} "
              f"{'COH_mu':>6} {'INT':<9} {'LON':>10} {'LAT':>9}")
        print(f"    {'─'*90}")
        for s in sorted(cluster_stats, key=lambda x: x["los_min_cm"]):
            print(f"    {s['cluster_id']:>4} {s['n_pixels']:>6} "
                  f"{s['area_ha']:>8.3f} "
                  f"{s['los_mean_cm']:>8.2f} {s['los_min_cm']:>8.2f} "
                  f"{s['los_max_cm']:>8.2f} "
                  f"{s['coh_mean']:>6.3f} {s['intensidad']:<9} "
                  f"{s['centroid_x']:>10.5f} {s['centroid_y']:>9.5f}")
        print(f"    {'─'*90}")

    # -- 4g. Exportar resultados ----------------------------------------------
    print(f"\n  Exportando...")

    # GeoPackage: poligonos
    gpkg_poly = None
    if clusters_gdf is not None and len(clusters_gdf) > 0:
        gpkg_poly = os.path.join(DEFORM_DIR, f"clusters_{ds_id}.gpkg")
        clusters_gdf.to_file(gpkg_poly, driver="GPKG")
        print(f"    ok {os.path.basename(gpkg_poly)}")

    # GeoPackage: centroides
    gpkg_cent = None
    if centroids_gdf is not None and len(centroids_gdf) > 0:
        gpkg_cent = os.path.join(DEFORM_DIR, f"centroids_{ds_id}.gpkg")
        centroids_gdf.to_file(gpkg_cent, driver="GPKG")
        print(f"    ok {os.path.basename(gpkg_cent)}")

    # CSV: estadisticas
    if cluster_stats:
        csv_path = os.path.join(DEFORM_DIR, f"summary_{ds_id}.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=cluster_stats[0].keys())
            writer.writeheader()
            writer.writerows(cluster_stats)
        print(f"    ok {os.path.basename(csv_path)}")

        # Agregar al catalogo global
        for s in cluster_stats:
            s_copy = dict(s)
            s_copy["dataset_id"] = ds_id
            s_copy["cat_temporal"] = ds["cat"]
            s_copy["baseline_d"] = ds["baseline"]
            s_copy["track"] = ds["track"]
            catalogo_clusters.append(s_copy)

    # Mapa
    map_path = os.path.join(DEFORM_DIR, f"map_{ds_id}.png")
    try:
        generar_mapa(los_cm, corr, mask, labeled, transform, crs,
                     ds_id, ds["cat"], ds["baseline"], map_path)
        figuras_generadas.append(map_path)
        print(f"    ok {os.path.basename(map_path)}")
    except Exception as e:
        print(f"    XX Error en mapa: {e}")

    # -- 4h. Registrar resultado global ---------------------------------------
    resultados_globales.append({
        "dataset_id":     ds_id,
        "status":         "OK",
        "tipo":           ds["tipo"],
        "track":          ds["track"],
        "dir":            ds["dir"],
        "cat":            ds["cat"],
        "baseline":       ds["baseline"],
        "cobertura":      ds["cobertura"],
        "coh_mean":       round(coh_mean, 3),
        "coh_pct_gt03":   round(coh_gt03, 1),
        "los_mean_cm":    round(los_mean, 2),
        "los_p5_cm":      round(los_p5, 2),
        "los_p95_cm":     round(los_p95, 2),
        "n_clusters":     n_clusters,
        "n_fuerte":       n_fuerte,
        "n_moderado":     n_moderado,
        "gpkg_poly":      gpkg_poly,
        "gpkg_cent":      gpkg_cent,
    })

# =============================================================================
# 5. CATALOGO UNIFICADO DE CLUSTERES
# =============================================================================
print(f"\n{'='*80}")
print("CATALOGO UNIFICADO")
print(f"{'='*80}")

if catalogo_clusters:
    cat_path = os.path.join(DEFORM_DIR, "P06_v2_catalogo_clusters.csv")
    with open(cat_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=catalogo_clusters[0].keys())
        writer.writeheader()
        writer.writerows(catalogo_clusters)
    print(f"\n  ok {cat_path}")
    print(f"    Total clusteres en catalogo: {len(catalogo_clusters)}")

    # Resumen por dataset
    ds_counts = Counter(c["dataset_id"] for c in catalogo_clusters)
    for ds_id, count in ds_counts.most_common():
        fuerte = sum(1 for c in catalogo_clusters
                     if c["dataset_id"] == ds_id and c["intensidad"] == "FUERTE")
        print(f"    {ds_id:<16}: {count} clusteres ({fuerte} fuertes)")
else:
    print("\n  !! No se detectaron clusteres en ningun dataset")

# =============================================================================
# 6. RESUMEN JSON
# =============================================================================
print(f"\n{'='*80}")
print("RESUMEN FINAL")
print(f"{'='*80}")

resumen = {
    "generado":      datetime.now().isoformat(),
    "version":       "P06 v2.0",
    "evento":        "2023-12-19",
    "parametros": {
        "coherencia_min":           COHERENCIA_MIN,
        "desplazamiento_min_cm":    DESPLAZAMIENTO_MIN_CM,
        "desplazamiento_fuerte_cm": DESPLAZAMIENTO_FUERTE_CM,
        "min_pixels_cluster":       MIN_PIXELS_CLUSTER,
        "factor_m_a_cm":            FACTOR_M_A_CM,
    },
    "clasificacion_temporal": {
        "PRE-EVENTO":  [p for p, i in PARES.items() if i["cat"] == "PRE-EVENTO"],
        "CO-EVENTO":   [p for p, i in PARES.items() if i["cat"] == "CO-EVENTO"],
        "POST-EVENTO": [p for p, i in PARES.items() if i["cat"] == "POST-EVENTO"],
    },
    "mosaicos_validos": list(MOSAICOS.keys()),
    "datasets_procesados": len([r for r in resultados_globales if r["status"] == "OK"]),
    "total_clusters":      len(catalogo_clusters),
    "resultados":          resultados_globales,
}

json_path = os.path.join(DEFORM_DIR, "P06_v2_resumen_completo.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(resumen, f, indent=2, ensure_ascii=False, default=str)
print(f"\n  ok Resumen: {json_path}")

# -- Tabla resumen ------------------------------------------------------------
print(f"\n  {'─'*85}")
print(f"  {'DATASET':<16} {'CAT':<14} {'BL':>3} {'COH_mu':>6} {'LOS_mu':>8} "
      f"{'LOS_P5':>8} {'LOS_P95':>8} {'CL':>4} {'F':>3}")
print(f"  {'─'*85}")
for r in resultados_globales:
    if r["status"] != "OK":
        print(f"  {r['dataset_id']:<16} {'-- SALTADO --'}")
        continue
    print(f"  {r['dataset_id']:<16} {r['cat']:<14} {r['baseline']:>3}d "
          f"{r['coh_mean']:>6.3f} {r['los_mean_cm']:>8.2f} "
          f"{r.get('los_p5_cm',0):>8.2f} {r.get('los_p95_cm',0):>8.2f} "
          f"{r['n_clusters']:>4} {r['n_fuerte']:>3}")
print(f"  {'─'*85}")
print(f"  CL = clusteres totales | F = clusteres FUERTES "
      f"(|LOS|>{DESPLAZAMIENTO_FUERTE_CM}cm)")

# -- Totales ------------------------------------------------------------------
total_ok = sum(1 for r in resultados_globales if r["status"] == "OK")
total_cl = sum(r.get("n_clusters", 0) for r in resultados_globales)
total_fu = sum(r.get("n_fuerte", 0) for r in resultados_globales)

print(f"\n  RESUMEN:")
print(f"    Datasets procesados:   {total_ok}/{len(DATASETS)}")
print(f"    Total clusteres:       {total_cl}")
print(f"    Clusteres FUERTES:     {total_fu}")
print(f"    Figuras generadas:     {len(figuras_generadas)}")

print(f"\n  ARCHIVOS DE SALIDA:")
print(f"    {DEFORM_DIR}/")
print(f"      clusters_*.gpkg     -> poligonos para QGIS")
print(f"      centroids_*.gpkg    -> centroides para QGIS")
print(f"      summary_*.csv       -> estadisticas por cluster")
print(f"      map_*.png           -> mapas de deformacion")
print(f"      P06_v2_catalogo_clusters.csv  -> catalogo unificado")
print(f"      P06_v2_resumen_completo.json  -> metadatos completos")

print(f"\n  SIGUIENTE PASO:")
print(f"    1. Abrir clusters_*.gpkg en QGIS sobre imagen Sentinel-2")
print(f"    2. Comparar centroides con inventario de deslizamientos")
print(f"    3. Validar con fotos aereas / noticias / reportes GAD")
print(f"    4. Integracion multi-fuente en P07")

print(f"\n{'='*80}")