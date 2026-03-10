# =============================================================================
# P06b_VERIFICACION — VERIFICACION CIENTIFICA INTEGRAL DE CLUSTERES InSAR
# =============================================================================
# Proyecto: InSAR Deslizamientos Zona de Intag, Cotacachi, Ecuador
# Autor:    Víctor Pinto
# Fecha:    2026-02-27
#
# PROPOSITO:
#   Aplicar 4 pruebas independientes de verificacion a cada cluster de
#   deformacion detectado por P06 v2.0, asignando una puntuacion de
#   confiabilidad (0-4) basada en multiples lineas de evidencia convergentes.
#
# PRUEBAS:
#   1. Comparacion temporal PRE vs CO vs POST (hipotesis nula)
#   2. Correlacion topografica con pendiente (DEM)
#   3. Analisis de patron espacial (tamano de cluster)
#   4. Sensibilidad al umbral de coherencia
#
# SISTEMA DE PUNTUACION:
#   4/4 = ALTA CONFIANZA    → Inventario final del articulo
#   3/4 = CONFIANZA MODERADA → Incluir con nota
#   2/4 = BAJA CONFIANZA    → Material suplementario
#   0-1 = FALSO POSITIVO    → Descartar
#
# ENTRADA:
#   hyp3_deformation/  → clusters y centroides de P06 v2.0
#   hyp3_clipped/      → rasters originales recortados
#   hyp3_mosaic/       → mosaicos CRB
#
# SALIDA:
#   hyp3_deformation/verificacion/
#     clusters_verificados_M_ASC_CO_12d.gpkg  → con puntuacion
#     reporte_verificacion.json               → estadisticas completas
#     fig_verificacion_*.png                  → figuras diagnosticas
#
# EJECUTAR EN: Spyder (Run file F5)
# PREREQUISITO: P06 v2.0 ejecutado con todos los datasets
# Ref: Protocolo de Verificacion Cientifica (P06_Protocolo_Verificacion.docx)
# =============================================================================

import json
import os
import sys
import warnings
from datetime import datetime

import numpy as np

print("=" * 80)
print("P06_VERIFICACION — VERIFICACION CIENTIFICA INTEGRAL")
print("=" * 80)
print(f"  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("\n--- Verificando librerias ---\n")

_requeridas = {
    "rasterio":   "rasterio",
    "geopandas":  "geopandas",
    "numpy":      "numpy",
    "matplotlib": "matplotlib",
    "scipy":      "scipy",
    "shapely":    "shapely",
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
    sys.exit(1)

import geopandas as gpd
import rasterio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.ndimage import label as scipy_label
from shapely.geometry import Point

warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
print("  ok Todas las librerias disponibles")

# =============================================================================
# 1. CONFIGURACION
# =============================================================================
print(f"\n{'='*80}")
print("1. CONFIGURACION")
print(f"{'='*80}")

SAVE_DIR     = r"D:\POSGRADOS\INTAG\data\sentinel1"
DEFORM_DIR   = os.path.join(SAVE_DIR, "hyp3_deformation")
CLIPPED_DIR  = os.path.join(SAVE_DIR, "hyp3_clipped")
MOSAIC_DIR   = os.path.join(SAVE_DIR, "hyp3_mosaic")
VERIF_DIR    = os.path.join(DEFORM_DIR, "verificacion")

os.makedirs(VERIF_DIR, exist_ok=True)

# Dataset principal a verificar
TARGET = "M_ASC_CO_12d"

# Rutas de rasters del mosaico co-evento 12d
TARGET_CORR = os.path.join(MOSAIC_DIR, TARGET, f"{TARGET}_corr.tif")
TARGET_LOS  = os.path.join(MOSAIC_DIR, TARGET, f"{TARGET}_los_disp.tif")
TARGET_DEM  = os.path.join(MOSAIC_DIR, TARGET, f"{TARGET}_dem.tif")

# Rutas de clusters (GeoPackage de centroides)
CENT_CO   = os.path.join(DEFORM_DIR, f"centroids_{TARGET}.gpkg")
CENT_PRE  = os.path.join(DEFORM_DIR, "centroids_P18_B1.gpkg")
CENT_POST = os.path.join(DEFORM_DIR, "centroids_P18_B3.gpkg")

# Parametros (identicos a P06 v2.0)
COHERENCIA_MIN           = 0.3
DESPLAZAMIENTO_MIN_CM    = 3.0
DESPLAZAMIENTO_FUERTE_CM = 5.0
MIN_PIXELS_CLUSTER       = 4
FACTOR_M_A_CM            = 100

# Parametros de verificacion
BUFFER_COINCIDENCIA_M    = 500     # metros para Test 1
PENDIENTE_MIN_GRADOS     = 5.0     # grados para Test 2
AREA_MAX_HA              = 500.0   # hectareas para Test 3
COHERENCIA_TEST_ALTA     = 0.4     # umbral alto para Test 4
COHERENCIA_TEST_MUY_ALTA = 0.5     # umbral muy alto para Test 4
RETENCION_MIN_PCT        = 50.0    # % pixeles retenidos para Test 4

print(f"  Target:              {TARGET}")
print(f"  Buffer coincidencia: {BUFFER_COINCIDENCIA_M} m")
print(f"  Pendiente minima:    {PENDIENTE_MIN_GRADOS}°")
print(f"  Area maxima:         {AREA_MAX_HA} ha")
print(f"  Coherencia test:     {COHERENCIA_TEST_ALTA} / {COHERENCIA_TEST_MUY_ALTA}")
print(f"  Retencion minima:    {RETENCION_MIN_PCT}%")
print(f"  Salida:              {VERIF_DIR}")

# Verificar archivos
archivos_requeridos = {
    "Raster coherencia": TARGET_CORR,
    "Raster LOS":        TARGET_LOS,
    "Raster DEM":        TARGET_DEM,
    "Centroides CO":     CENT_CO,
    "Centroides PRE":    CENT_PRE,
    "Centroides POST":   CENT_POST,
}

print(f"\n  Verificando archivos:")
todos_ok = True
for nombre, ruta in archivos_requeridos.items():
    existe = os.path.exists(ruta)
    print(f"    {'ok' if existe else 'XX'} {nombre}: {os.path.basename(ruta)}")
    if not existe:
        todos_ok = False

if not todos_ok:
    print(f"\n  [ERROR] Archivos faltantes. Ejecutar P06 v2.0 primero.")
    sys.exit(1)

# =============================================================================
# 2. CARGAR DATOS
# =============================================================================
print(f"\n{'='*80}")
print("2. CARGANDO DATOS")
print(f"{'='*80}")

# 2a. Centroides
print(f"\n  Cargando centroides...")
gdf_co   = gpd.read_file(CENT_CO)
gdf_pre  = gpd.read_file(CENT_PRE)
gdf_post = gpd.read_file(CENT_POST)

print(f"    CO-EVENTO ({TARGET}):  {len(gdf_co)} clusteres")
print(f"    PRE-EVENTO (P18_B1):   {len(gdf_pre)} clusteres")
print(f"    POST-EVENTO (P18_B3):  {len(gdf_post)} clusteres")
print(f"    CRS: {gdf_co.crs}")

# 2b. Rasters
print(f"\n  Cargando rasters...")


def cargar_raster(path, banda=1):
    """Carga un raster GeoTIFF y retorna (data, transform, crs, nodata)."""
    with rasterio.open(path) as src:
        data = src.read(banda).astype(np.float64)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs
    if nodata is not None:
        data[data == nodata] = np.nan
    return data, transform, crs, nodata


corr_data, corr_transform, corr_crs, _ = cargar_raster(TARGET_CORR)
los_raw, _, _, _ = cargar_raster(TARGET_LOS)
dem_data, _, _, _ = cargar_raster(TARGET_DEM)

los_cm = los_raw * FACTOR_M_A_CM

print(f"    Coherencia: {corr_data.shape}")
print(f"    LOS (cm):   min={np.nanmin(los_cm):.2f}, max={np.nanmax(los_cm):.2f}")
print(f"    DEM:         min={np.nanmin(dem_data):.0f}m, max={np.nanmax(dem_data):.0f}m")

# =============================================================================
# 3. CALCULAR PENDIENTE DESDE DEM
# =============================================================================
print(f"\n{'='*80}")
print("3. CALCULANDO PENDIENTE")
print(f"{'='*80}")

# Tamano de pixel en metros
pixel_x = abs(corr_transform.a)
pixel_y = abs(corr_transform.e)

# Gradiente con numpy
dy, dx = np.gradient(dem_data, pixel_y, pixel_x)
slope_rad = np.arctan(np.sqrt(dx**2 + dy**2))
slope_deg = np.degrees(slope_rad)

# Limpiar donde DEM es NaN
slope_deg[~np.isfinite(dem_data)] = np.nan

print(f"  Pixel size: {pixel_x:.1f} x {pixel_y:.1f} m")
print(f"  Pendiente: min={np.nanmin(slope_deg):.1f}°, "
      f"max={np.nanmax(slope_deg):.1f}°, "
      f"mean={np.nanmean(slope_deg):.1f}°")

# =============================================================================
# 4. RE-GENERAR MASCARA Y LABELED PARA TEST 4
# =============================================================================
print(f"\n{'='*80}")
print("4. RE-GENERANDO CLUSTERES PARA SENSIBILIDAD DE COHERENCIA")
print(f"{'='*80}")

# Mascara original (gamma > 0.3 AND |LOS| > 3 cm)
mask_coh_03 = (corr_data > COHERENCIA_MIN) & np.isfinite(corr_data)
mask_disp = (np.abs(los_cm) > DESPLAZAMIENTO_MIN_CM) & np.isfinite(los_cm)
mask_original = mask_coh_03 & mask_disp

structure = np.ones((3, 3), dtype=int)
labeled_03, n_features_03 = scipy_label(mask_original.astype(int), structure=structure)
print(f"  Mascara original (gamma>0.3): {n_features_03} features antes de filtro de tamano")

# Mascara con gamma > 0.4
mask_coh_04 = (corr_data > COHERENCIA_TEST_ALTA) & np.isfinite(corr_data)
mask_04 = mask_coh_04 & mask_disp
labeled_04, n_04 = scipy_label(mask_04.astype(int), structure=structure)
print(f"  Mascara gamma>0.4: {n_04} features")

# Mascara con gamma > 0.5
mask_coh_05 = (corr_data > COHERENCIA_TEST_MUY_ALTA) & np.isfinite(corr_data)
mask_05 = mask_coh_05 & mask_disp
labeled_05, n_05 = scipy_label(mask_05.astype(int), structure=structure)
print(f"  Mascara gamma>0.5: {n_05} features")

# =============================================================================
# 5. APLICAR 4 PRUEBAS A CADA CLUSTER
# =============================================================================
print(f"\n{'='*80}")
print("5. APLICANDO 4 PRUEBAS A CADA CLUSTER")
print(f"{'='*80}")

resultados = []

for idx, row in gdf_co.iterrows():
    cl_id = row["cluster_id"]
    cx = row.geometry.x
    cy = row.geometry.y
    area_ha = row.get("area_ha", 0)
    n_pixels = row.get("n_pixels", 0)
    los_mean = row.get("los_mean_cm", 0)
    los_min = row.get("los_min_cm", 0)
    los_max = row.get("los_max_cm", 0)
    coh_mean = row.get("coh_mean", 0)
    intensidad = row.get("intensidad", "?")

    resultado = {
        "cluster_id": cl_id,
        "n_pixels":   n_pixels,
        "area_ha":    area_ha,
        "los_mean_cm": los_mean,
        "los_min_cm":  los_min,
        "los_max_cm":  los_max,
        "coh_mean":    coh_mean,
        "intensidad":  intensidad,
        "cx": cx, "cy": cy,
    }

    # ── PRUEBA 1: Comparacion temporal ──────────────────────────────────
    # ¿Existe un cluster PRE o POST dentro de BUFFER_COINCIDENCIA_M?
    punto_co = Point(cx, cy)
    buffer_co = punto_co.buffer(BUFFER_COINCIDENCIA_M)

    # Buscar coincidencia con PRE
    coincide_pre = False
    for _, row_pre in gdf_pre.iterrows():
        if buffer_co.contains(row_pre.geometry):
            coincide_pre = True
            break

    # Buscar coincidencia con POST
    coincide_post = False
    for _, row_post in gdf_post.iterrows():
        if buffer_co.contains(row_post.geometry):
            coincide_post = True
            break

    # PASA si NO coincide con PRE (señal unica del co-evento)
    test1_pasa = not coincide_pre
    resultado["T1_coincide_pre"]  = coincide_pre
    resultado["T1_coincide_post"] = coincide_post
    resultado["T1_pasa"] = test1_pasa

    # ── PRUEBA 2: Pendiente ─────────────────────────────────────────────
    # Extraer pendiente en el centroide
    try:
        col_idx = int((cx - corr_transform.c) / corr_transform.a)
        row_idx = int((cy - corr_transform.f) / corr_transform.e)

        # Asegurar que estamos dentro de los limites
        row_idx = max(0, min(row_idx, slope_deg.shape[0] - 1))
        col_idx = max(0, min(col_idx, slope_deg.shape[1] - 1))

        # Pendiente en ventana 3x3 alrededor del centroide (mas robusto)
        r_min = max(0, row_idx - 1)
        r_max = min(slope_deg.shape[0], row_idx + 2)
        c_min = max(0, col_idx - 1)
        c_max = min(slope_deg.shape[1], col_idx + 2)
        ventana = slope_deg[r_min:r_max, c_min:c_max]
        slope_val = float(np.nanmean(ventana)) if np.any(np.isfinite(ventana)) else 0.0

        # Elevacion
        ventana_dem = dem_data[r_min:r_max, c_min:c_max]
        elev_val = float(np.nanmean(ventana_dem)) if np.any(np.isfinite(ventana_dem)) else 0.0
    except Exception:
        slope_val = 0.0
        elev_val = 0.0

    test2_pasa = slope_val >= PENDIENTE_MIN_GRADOS
    resultado["T2_pendiente_deg"] = round(slope_val, 2)
    resultado["T2_elevacion_m"]   = round(elev_val, 1)
    resultado["T2_pasa"] = test2_pasa

    # ── PRUEBA 3: Tamano de cluster ─────────────────────────────────────
    test3_pasa = area_ha <= AREA_MAX_HA
    if area_ha <= 100:
        t3_clase = "COMPATIBLE"
    elif area_ha <= 500:
        t3_clase = "SOSPECHOSO"
    else:
        t3_clase = "ARTEFACTO"

    resultado["T3_area_clase"] = t3_clase
    resultado["T3_pasa"] = test3_pasa

    # ── PRUEBA 4: Sensibilidad de coherencia ────────────────────────────
    # Verificar cuantos pixeles del cluster sobreviven gamma > 0.4
    try:
        # Ubicar pixeles del cluster en labeled_03
        # Usar el cluster_id del labeled original
        # Necesitamos encontrar el label correspondiente en labeled_03
        # Buscar el label en la posicion del centroide
        label_at_centroid = labeled_03[row_idx, col_idx]

        if label_at_centroid > 0:
            # Mascara de este cluster
            mask_this = (labeled_03 == label_at_centroid)
            n_pix_03 = int(np.sum(mask_this))

            # Cuantos de estos pixeles tienen gamma > 0.4
            n_pix_04 = int(np.sum(mask_this & mask_coh_04))
            pct_04 = (n_pix_04 / n_pix_03 * 100) if n_pix_03 > 0 else 0

            # Cuantos tienen gamma > 0.5
            n_pix_05 = int(np.sum(mask_this & mask_coh_05))
            pct_05 = (n_pix_05 / n_pix_03 * 100) if n_pix_03 > 0 else 0
        else:
            n_pix_03 = n_pixels
            pct_04 = 0.0
            pct_05 = 0.0
    except Exception:
        n_pix_03 = n_pixels
        pct_04 = 0.0
        pct_05 = 0.0

    test4_pasa = pct_04 >= RETENCION_MIN_PCT
    resultado["T4_pix_03"]  = n_pix_03
    resultado["T4_pct_04"]  = round(pct_04, 1)
    resultado["T4_pct_05"]  = round(pct_05, 1)
    resultado["T4_pasa"] = test4_pasa

    # ── PUNTUACION FINAL ────────────────────────────────────────────────
    score = sum([test1_pasa, test2_pasa, test3_pasa, test4_pasa])

    if score == 4:
        confianza = "ALTA"
    elif score == 3:
        confianza = "MODERADA"
    elif score == 2:
        confianza = "BAJA"
    else:
        confianza = "FALSO_POSITIVO"

    resultado["score"] = score
    resultado["confianza"] = confianza

    resultados.append(resultado)

print(f"  Clusteres evaluados: {len(resultados)}")

# =============================================================================
# 6. RESUMEN DE RESULTADOS
# =============================================================================
print(f"\n{'='*80}")
print("6. RESUMEN DE VERIFICACION")
print(f"{'='*80}")

# Conteo por nivel de confianza
from collections import Counter
confianza_counts = Counter(r["confianza"] for r in resultados)
score_counts = Counter(r["score"] for r in resultados)

print(f"\n  PUNTUACION DE CONFIABILIDAD:")
print(f"  {'─'*60}")
for nivel, etiqueta in [(4, "ALTA (4/4)"), (3, "MODERADA (3/4)"),
                         (2, "BAJA (2/4)"), (1, "FALSO POSITIVO (1/4)"),
                         (0, "FALSO POSITIVO (0/4)")]:
    n = score_counts.get(nivel, 0)
    pct = n / len(resultados) * 100 if resultados else 0
    barra = "#" * int(pct / 2)
    print(f"    {etiqueta:25s} {n:5d} ({pct:5.1f}%) {barra}")
print(f"  {'─'*60}")
print(f"    TOTAL:                    {len(resultados):5d}")

# Resumen de cada prueba individual
print(f"\n  RESULTADOS POR PRUEBA:")
print(f"  {'─'*60}")
for prueba, nombre in [("T1_pasa", "T1 Temporal (no coincide PRE)"),
                        ("T2_pasa", "T2 Pendiente (>= 5°)"),
                        ("T3_pasa", "T3 Tamano (<= 500 ha)"),
                        ("T4_pasa", "T4 Coherencia (>= 50% en 0.4)")]:
    n_pasa = sum(1 for r in resultados if r[prueba])
    pct = n_pasa / len(resultados) * 100 if resultados else 0
    print(f"    {nombre:40s} PASA: {n_pasa:5d} ({pct:5.1f}%)")
print(f"  {'─'*60}")

# Estadisticas de clusters de ALTA CONFIANZA
alta = [r for r in resultados if r["confianza"] == "ALTA"]
print(f"\n  CLUSTERS DE ALTA CONFIANZA ({len(alta)}):")
if alta:
    los_vals = [r["los_mean_cm"] for r in alta]
    area_vals = [r["area_ha"] for r in alta]
    slope_vals = [r["T2_pendiente_deg"] for r in alta]
    elev_vals = [r["T2_elevacion_m"] for r in alta]
    print(f"    LOS media:     {np.mean(los_vals):.2f} cm "
          f"(min={np.min(los_vals):.2f}, max={np.max(los_vals):.2f})")
    print(f"    Area media:    {np.mean(area_vals):.2f} ha "
          f"(min={np.min(area_vals):.2f}, max={np.max(area_vals):.2f})")
    print(f"    Pendiente:     {np.mean(slope_vals):.1f}° "
          f"(min={np.min(slope_vals):.1f}°, max={np.max(slope_vals):.1f}°)")
    print(f"    Elevacion:     {np.mean(elev_vals):.0f} m "
          f"(min={np.min(elev_vals):.0f} m, max={np.max(elev_vals):.0f} m)")

    # Desglose por intensidad
    fuerte_alta = sum(1 for r in alta if r["intensidad"] == "FUERTE")
    print(f"    FUERTES:       {fuerte_alta}")
    print(f"    MODERADOS:     {len(alta) - fuerte_alta}")

# Estadisticas del test 1 (temporal)
print(f"\n  DETALLE TEST 1 (Comparacion temporal):")
n_coincide_pre = sum(1 for r in resultados if r["T1_coincide_pre"])
n_coincide_post = sum(1 for r in resultados if r["T1_coincide_post"])
n_unico = sum(1 for r in resultados
              if not r["T1_coincide_pre"] and not r["T1_coincide_post"])
print(f"    Coincide con PRE:   {n_coincide_pre:5d} "
      f"({n_coincide_pre/len(resultados)*100:.1f}%)")
print(f"    Coincide con POST:  {n_coincide_post:5d} "
      f"({n_coincide_post/len(resultados)*100:.1f}%)")
print(f"    Unico (no PRE, no POST): {n_unico:5d} "
      f"({n_unico/len(resultados)*100:.1f}%)")

# Estadisticas de pendiente
print(f"\n  DETALLE TEST 2 (Pendiente):")
slopes = [r["T2_pendiente_deg"] for r in resultados]
for rango, lo, hi in [("< 5° (sospechoso)", 0, 5),
                       ("5-15° (suave)", 5, 15),
                       ("15-30° (moderada)", 15, 30),
                       ("> 30° (pronunciada)", 30, 90)]:
    n = sum(1 for s in slopes if lo <= s < hi)
    print(f"    {rango:25s} {n:5d} ({n/len(resultados)*100:.1f}%)")

# Estadisticas de tamano
print(f"\n  DETALLE TEST 3 (Tamano):")
for clase in ["COMPATIBLE", "SOSPECHOSO", "ARTEFACTO"]:
    n = sum(1 for r in resultados if r["T3_area_clase"] == clase)
    print(f"    {clase:20s} {n:5d} ({n/len(resultados)*100:.1f}%)")

# Estadisticas de coherencia
print(f"\n  DETALLE TEST 4 (Sensibilidad coherencia):")
pcts_04 = [r["T4_pct_04"] for r in resultados]
print(f"    % pixeles con gamma>0.4: "
      f"mean={np.mean(pcts_04):.1f}%, median={np.median(pcts_04):.1f}%")
pcts_05 = [r["T4_pct_05"] for r in resultados]
print(f"    % pixeles con gamma>0.5: "
      f"mean={np.mean(pcts_05):.1f}%, median={np.median(pcts_05):.1f}%")

# =============================================================================
# 7. EXPORTAR RESULTADOS
# =============================================================================
print(f"\n{'='*80}")
print("7. EXPORTANDO RESULTADOS")
print(f"{'='*80}")

# 7a. GeoPackage con puntuacion
print(f"\n  Creando GeoPackage verificado...")
geoms = [Point(r["cx"], r["cy"]) for r in resultados]
gdf_verif = gpd.GeoDataFrame(resultados, geometry=geoms, crs=gdf_co.crs)

# Eliminar columnas de geometria redundantes
for col in ["cx", "cy"]:
    if col in gdf_verif.columns:
        gdf_verif = gdf_verif.drop(columns=[col])

gpkg_path = os.path.join(VERIF_DIR, f"clusters_verificados_{TARGET}.gpkg")
gdf_verif.to_file(gpkg_path, driver="GPKG")
print(f"    ok {os.path.basename(gpkg_path)}")

# 7b. CSV
import csv
csv_path = os.path.join(VERIF_DIR, f"verificacion_{TARGET}.csv")
keys = resultados[0].keys() if resultados else []
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=keys)
    writer.writeheader()
    writer.writerows(resultados)
print(f"    ok {os.path.basename(csv_path)}")

# 7c. JSON reporte completo
reporte = {
    "generado": datetime.now().isoformat(),
    "version": "P06_VERIFICACION v1.0",
    "target": TARGET,
    "parametros_p06": {
        "coherencia_min": COHERENCIA_MIN,
        "desplazamiento_min_cm": DESPLAZAMIENTO_MIN_CM,
        "desplazamiento_fuerte_cm": DESPLAZAMIENTO_FUERTE_CM,
        "min_pixels_cluster": MIN_PIXELS_CLUSTER,
    },
    "parametros_verificacion": {
        "buffer_coincidencia_m": BUFFER_COINCIDENCIA_M,
        "pendiente_min_grados": PENDIENTE_MIN_GRADOS,
        "area_max_ha": AREA_MAX_HA,
        "coherencia_test_alta": COHERENCIA_TEST_ALTA,
        "coherencia_test_muy_alta": COHERENCIA_TEST_MUY_ALTA,
        "retencion_min_pct": RETENCION_MIN_PCT,
    },
    "datasets_comparacion": {
        "PRE": {"id": "P18_B1", "n_clusters": len(gdf_pre), "baseline": 24},
        "CO":  {"id": TARGET, "n_clusters": len(gdf_co), "baseline": 12},
        "POST": {"id": "P18_B3", "n_clusters": len(gdf_post), "baseline": 24},
    },
    "resumen_confianza": {
        "ALTA": confianza_counts.get("ALTA", 0),
        "MODERADA": confianza_counts.get("MODERADA", 0),
        "BAJA": confianza_counts.get("BAJA", 0),
        "FALSO_POSITIVO": confianza_counts.get("FALSO_POSITIVO", 0),
    },
    "resumen_pruebas": {
        "T1_pasa": sum(1 for r in resultados if r["T1_pasa"]),
        "T2_pasa": sum(1 for r in resultados if r["T2_pasa"]),
        "T3_pasa": sum(1 for r in resultados if r["T3_pasa"]),
        "T4_pasa": sum(1 for r in resultados if r["T4_pasa"]),
    },
}
json_path = os.path.join(VERIF_DIR, f"reporte_verificacion_{TARGET}.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(reporte, f, indent=2, ensure_ascii=False, default=str)
print(f"    ok {os.path.basename(json_path)}")

# =============================================================================
# 8. FIGURAS DIAGNOSTICAS
# =============================================================================
print(f"\n{'='*80}")
print("8. GENERANDO FIGURAS")
print(f"{'='*80}")

# 8a. Mapa de confianza
print(f"\n  Generando mapa de confianza...")
fig, ax = plt.subplots(1, 1, figsize=(12, 10))

# Fondo: coherencia
extent = rasterio.transform.array_bounds(
    corr_data.shape[0], corr_data.shape[1], corr_transform)
extent_plot = [extent[0], extent[2], extent[1], extent[3]]
ax.imshow(corr_data, cmap="gray", vmin=0, vmax=1,
          extent=extent_plot, origin="upper", alpha=0.4)

# Clusters coloreados por confianza
colores = {"ALTA": "#27AE60", "MODERADA": "#2471A3",
           "BAJA": "#E67E22", "FALSO_POSITIVO": "#C0392B"}
for nivel in ["FALSO_POSITIVO", "BAJA", "MODERADA", "ALTA"]:
    sub = [r for r in resultados if r["confianza"] == nivel]
    if sub:
        xs = [r["geometry"].x if hasattr(r.get("geometry", None), "x")
              else r.get("cx", 0) for r in sub]
        ys = [r["geometry"].y if hasattr(r.get("geometry", None), "y")
              else r.get("cy", 0) for r in sub]
        # Fallback: usar gdf_verif
        sub_gdf = gdf_verif[gdf_verif["confianza"] == nivel]
        ax.scatter(sub_gdf.geometry.x, sub_gdf.geometry.y,
                   c=colores[nivel], s=8, alpha=0.7, label=f"{nivel} ({len(sub_gdf)})",
                   edgecolors="none")

ax.legend(loc="upper right", fontsize=9)
ax.set_title(f"Verificación Científica — {TARGET}\n"
             f"Confiabilidad de clústeres (4 pruebas convergentes)",
             fontsize=12, fontweight="bold")
ax.set_xlabel("Este (m)")
ax.set_ylabel("Norte (m)")

fig.tight_layout()
map_path = os.path.join(VERIF_DIR, f"fig_mapa_confianza_{TARGET}.png")
fig.savefig(map_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"    ok {os.path.basename(map_path)}")

# 8b. Dashboard de 4 pruebas
print(f"  Generando dashboard de pruebas...")
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Panel 1: Comparacion temporal (barras PRE vs CO vs POST)
ax1 = axes[0, 0]
cats = ["PRE\n(P18_B1)", f"CO\n({TARGET})", "POST\n(P18_B3)"]
vals = [len(gdf_pre), len(gdf_co), len(gdf_post)]
bars = ax1.bar(cats, vals, color=["#3498DB", "#E74C3C", "#F39C12"], edgecolor="black")
for bar, val in zip(bars, vals):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
             str(val), ha="center", fontweight="bold")
ax1.set_ylabel("Número de clústeres")
ax1.set_title("T1: Comparación temporal\n(PRE vs CO vs POST)", fontweight="bold")

# Panel 2: Distribucion de pendiente
ax2 = axes[0, 1]
slopes_arr = np.array([r["T2_pendiente_deg"] for r in resultados])
ax2.hist(slopes_arr[slopes_arr < 50], bins=30, color="#2ECC71",
         edgecolor="black", alpha=0.8)
ax2.axvline(PENDIENTE_MIN_GRADOS, color="red", linestyle="--",
            label=f"Umbral = {PENDIENTE_MIN_GRADOS}°")
ax2.set_xlabel("Pendiente (grados)")
ax2.set_ylabel("Número de clústeres")
ax2.set_title("T2: Distribución de pendiente\nen centroides", fontweight="bold")
ax2.legend()

# Panel 3: Distribucion de tamano (log)
ax3 = axes[1, 0]
areas = np.array([r["area_ha"] for r in resultados])
areas_log = areas[areas > 0]
ax3.hist(np.log10(areas_log), bins=30, color="#9B59B6",
         edgecolor="black", alpha=0.8)
ax3.axvline(np.log10(100), color="orange", linestyle="--", label="100 ha")
ax3.axvline(np.log10(500), color="red", linestyle="--", label="500 ha")
ax3.set_xlabel("log₁₀(Área en ha)")
ax3.set_ylabel("Número de clústeres")
ax3.set_title("T3: Distribución de tamaño\n(escala logarítmica)", fontweight="bold")
ax3.legend()

# Panel 4: Sensibilidad de coherencia
ax4 = axes[1, 1]
pcts = np.array([r["T4_pct_04"] for r in resultados])
ax4.hist(pcts, bins=30, color="#E67E22", edgecolor="black", alpha=0.8)
ax4.axvline(RETENCION_MIN_PCT, color="red", linestyle="--",
            label=f"Umbral = {RETENCION_MIN_PCT}%")
ax4.set_xlabel("% píxeles retenidos con γ > 0.4")
ax4.set_ylabel("Número de clústeres")
ax4.set_title("T4: Sensibilidad de coherencia\n(retención con γ > 0.4)", fontweight="bold")
ax4.legend()

fig.suptitle(f"Verificación Científica — {TARGET} — {len(resultados)} clústeres",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
dash_path = os.path.join(VERIF_DIR, f"fig_dashboard_verificacion_{TARGET}.png")
fig.savefig(dash_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"    ok {os.path.basename(dash_path)}")

# 8c. Distribucion de puntuaciones (pie chart)
print(f"  Generando grafico de puntuaciones...")
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
labels = []
sizes = []
colors_pie = []
for nivel, color in [("ALTA", "#27AE60"), ("MODERADA", "#2471A3"),
                      ("BAJA", "#E67E22"), ("FALSO_POSITIVO", "#C0392B")]:
    n = confianza_counts.get(nivel, 0)
    if n > 0:
        pct = n / len(resultados) * 100
        labels.append(f"{nivel}\n({n}, {pct:.1f}%)")
        sizes.append(n)
        colors_pie.append(color)

ax.pie(sizes, labels=labels, colors=colors_pie, autopct="",
       startangle=90, textprops={"fontsize": 11})
ax.set_title(f"Distribución de confiabilidad — {TARGET}\n"
             f"Total: {len(resultados)} clústeres", fontsize=13, fontweight="bold")

fig.tight_layout()
pie_path = os.path.join(VERIF_DIR, f"fig_confianza_pie_{TARGET}.png")
fig.savefig(pie_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"    ok {os.path.basename(pie_path)}")

# =============================================================================
# 9. RESUMEN FINAL
# =============================================================================
print(f"\n{'='*80}")
print("RESUMEN FINAL")
print(f"{'='*80}")

print(f"""
  TARGET: {TARGET} ({len(resultados)} clusteres evaluados)

  CONFIABILIDAD:
    ALTA (4/4):           {confianza_counts.get('ALTA', 0):5d} clusteres
    MODERADA (3/4):       {confianza_counts.get('MODERADA', 0):5d} clusteres
    BAJA (2/4):           {confianza_counts.get('BAJA', 0):5d} clusteres
    FALSO POSITIVO (0-1): {confianza_counts.get('FALSO_POSITIVO', 0):5d} clusteres

  ARCHIVOS DE SALIDA:
    {VERIF_DIR}/
      clusters_verificados_{TARGET}.gpkg   → centroides con puntuacion
      verificacion_{TARGET}.csv            → tabla completa
      reporte_verificacion_{TARGET}.json   → metadatos
      fig_mapa_confianza_{TARGET}.png      → mapa por confianza
      fig_dashboard_verificacion_{TARGET}.png → dashboard 4 pruebas
      fig_confianza_pie_{TARGET}.png       → distribucion de confianza

  SIGUIENTE PASO:
    1. Abrir clusters_verificados_*.gpkg en QGIS
    2. Filtrar por confianza = 'ALTA' para inventario final
    3. Comparar ALTA CONFIANZA con inventario GAD / fotos aereas
    4. Documentar en articulo cientifico
""")

print(f"{'='*80}")
print("  VERIFICACION COMPLETADA")
print(f"{'='*80}")