# ============================================================================
# P07c — INVENTARIO ÓPTICO INDEPENDIENTE DE DESLIZAMIENTOS (Fase 1)
# ============================================================================
# Proyecto: InSAR Deslizamientos — Zona de Intag, Cotacachi, Ecuador
# Protocolo: Maestro InSAR Intag v2.0 — Fase 1, Paso 1.1
# Fecha: 28 de febrero de 2026
# Autor: Víctor Pinto
#
# PROPÓSITO:
#   Generar el inventario óptico independiente de deslizamientos usando
#   exclusivamente datos Sentinel-2 + DEM + cobertura del suelo.
#   Este inventario es un PRODUCTO AUTÓNOMO de la Fase 1, independiente
#   del análisis InSAR. Posteriormente se cruzará con el inventario InSAR
#   en el Paso 1.3 para asignar niveles de confianza.
#
# METODOLOGÍA (Notti et al. 2023, NHESS 23:2625-2648):
#   Candidato óptico = dNDVI < umbral AND pendiente ≥ 15°
#   Filtros: excluir agua (WorldCover 80) y urbano (WorldCover 50)
#   Evidencia adicional: dNBR > umbral, dBSI > umbral
#   Análisis de sensibilidad: umbrales -0.10, -0.15, -0.20
#
# INPUTS:
#   De P07 (GEE):
#     - P07_dNDVI.tif, P07_dNBR.tif, P07_dBSI.tif
#     - P07_NDVI_pre.tif, P07_n_obs.tif
#   De P07c_auxiliar (GEE):
#     - P07c_slope_deg.tif (Copernicus GLO-30)
#     - P07c_worldcover.tif (ESA WorldCover 2021)
#
# OUTPUTS:
#   - P07c_inventario_optico.gpkg (polígonos con atributos)
#   - P07c_inventario_optico.csv (tabla exportable)
#   - P07c_sensibilidad_umbrales.csv (análisis de sensibilidad)
#   - P07c_mapa_inventario.png (mapa del inventario)
#   - P07c_histograma_atributos.png (distribución de atributos)
#   - P07c_log.txt (log de ejecución)
#
# REFERENCIAS:
#   Notti, D. et al. (2023). NHESS 23:2625-2648.
#   Behling, R. et al. (2014). Remote Sensing, 6(9), 8026-8055.
#   ESA WorldCover 2021: doi.org/10.5281/zenodo.7254221
#   Copernicus GLO-30: doi.org/10.5270/ESA-c5d3d65
# ============================================================================

import os
import sys
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import shapes
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
import pandas as pd
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# SECCIÓN 0: CONFIGURACIÓN
# ============================================================================

# --- RUTAS ---
DIR_OPTICO = r"D:\POSGRADOS\INTAG\data\sentinel2\VALIDACION\INTAG_InSAR_P07"
DIR_OUTPUT = os.path.join(DIR_OPTICO, "inventario_optico")
os.makedirs(DIR_OUTPUT, exist_ok=True)

# --- ARCHIVOS DE ENTRADA ---
RASTERS = {
    "dNDVI": os.path.join(DIR_OPTICO, "P07_dNDVI.tif"),
    "dNBR": os.path.join(DIR_OPTICO, "P07_dNBR.tif"),
    "dBSI": os.path.join(DIR_OPTICO, "P07_dBSI.tif"),
    "NDVI_pre": os.path.join(DIR_OPTICO, "P07_NDVI_pre.tif"),
    "n_obs": os.path.join(DIR_OPTICO, "P07_n_obs.tif"),
    "slope": os.path.join(DIR_OPTICO, "P07c_slope_deg.tif"),
    "worldcover": os.path.join(DIR_OPTICO, "P07c_worldcover.tif"),
}

# --- PARÁMETROS CIENTÍFICOS ---
# Umbral principal de dNDVI (Notti et al. 2023, NHESS 23:2625-2648)
DNDVI_THRESHOLD = -0.15

# Pendiente mínima para candidato a deslizamiento (Notti et al. 2023)
SLOPE_MIN_DEG = 15.0

# Umbrales de índices adicionales (mismos que P07b)
DNBR_THRESHOLD = 0.05   # Behling et al. (2014)
DBSI_THRESHOLD = 0.02   # Ariza & Dávila (2021)

# Calidad mínima del composite (observaciones S2)
MIN_OBS = 3

# Área mínima de candidato (píxeles, a 10m = 100 m² por píxel)
MIN_PIXELS = 10  # = 1,000 m² = 0.1 ha

# Umbrales para análisis de sensibilidad
SENSITIVITY_THRESHOLDS = [-0.05, -0.10, -0.15, -0.20, -0.25]

# Clases WorldCover a excluir
EXCLUDE_CLASSES = [50, 80]  # 50=urbano, 80=agua

# --- LOG ---
log_lines = []
def log(msg):
    timestamp = datetime.now().strftime("%H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    log_lines.append(line)

# ============================================================================
# SECCIÓN 1: VERIFICACIÓN DE INPUTS
# ============================================================================

log("=" * 70)
log("P07c — INVENTARIO ÓPTICO INDEPENDIENTE (Fase 1 del Protocolo)")
log(f"Fecha de ejecución: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)
log("")
log("--- SECCIÓN 1: VERIFICACIÓN DE INPUTS ---")

ref_crs = None
ref_shape = None
ref_transform = None
all_ok = True

for name, path in RASTERS.items():
    if os.path.exists(path):
        with rasterio.open(path) as src:
            crs = src.crs
            h, w = src.height, src.width
            res = src.res[0]
            log(f"  OK: {name} — {h}x{w} px, CRS={crs}, res={res}m")
            if ref_crs is None:
                ref_crs = crs
                ref_shape = (h, w)
                ref_transform = src.transform
            # Verificar compatibilidad
            if crs != ref_crs:
                log(f"  ALERTA: {name} tiene CRS diferente ({crs} vs {ref_crs})")
    else:
        log(f"  FALTA: {path}")
        all_ok = False

if not all_ok:
    log("")
    log("ERROR: Faltan archivos. Ejecutar P07c_auxiliar en GEE primero.")
    log("Archivos necesarios de GEE: P07c_slope_deg.tif, P07c_worldcover.tif")
    sys.exit(1)

log("")

# ============================================================================
# SECCIÓN 2: CARGA Y ALINEACIÓN DE RASTERS
# ============================================================================

log("--- SECCIÓN 2: CARGA Y ALINEACIÓN DE RASTERS ---")

# Cargar raster de referencia (dNDVI)
with rasterio.open(RASTERS["dNDVI"]) as src_ref:
    ref_profile = src_ref.profile.copy()
    ref_crs = src_ref.crs
    ref_transform = src_ref.transform
    ref_shape = (src_ref.height, src_ref.width)
    ref_bounds = src_ref.bounds
    dNDVI = src_ref.read(1)

log(f"  Referencia: dNDVI {ref_shape[0]}x{ref_shape[1]} px, CRS={ref_crs}")

# Función para cargar y alinear un raster al de referencia
def load_aligned(path, name, method='bilinear'):
    with rasterio.open(path) as src:
        if src.crs == ref_crs and src.shape == ref_shape and src.transform == ref_transform:
            data = src.read(1)
            log(f"  {name}: cargado directo (misma grilla)")
            return data
        else:
            # Reproyectar/alinear
            resamp = Resampling.bilinear if method == 'bilinear' else Resampling.nearest
            dst = np.empty(ref_shape, dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1),
                destination=dst,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=ref_transform,
                dst_crs=ref_crs,
                resampling=resamp
            )
            log(f"  {name}: reproyectado/alineado ({src.shape} → {ref_shape})")
            return dst

# Cargar todos los rasters alineados
dNBR = load_aligned(RASTERS["dNBR"], "dNBR")
dBSI = load_aligned(RASTERS["dBSI"], "dBSI")
NDVI_pre = load_aligned(RASTERS["NDVI_pre"], "NDVI_pre")
n_obs = load_aligned(RASTERS["n_obs"], "n_obs")
slope = load_aligned(RASTERS["slope"], "slope", method='bilinear')
worldcover = load_aligned(RASTERS["worldcover"], "worldcover", method='nearest')

log("")

# ============================================================================
# SECCIÓN 3: CREACIÓN DE MÁSCARAS
# ============================================================================

log("--- SECCIÓN 3: CREACIÓN DE MÁSCARAS ---")

# 3.1 Máscara de datos válidos (dNDVI no es NaN/NoData)
valid_mask = np.isfinite(dNDVI) & (dNDVI != 0)
n_valid = np.sum(valid_mask)
log(f"  Píxeles válidos (dNDVI): {n_valid:,} de {dNDVI.size:,} ({100*n_valid/dNDVI.size:.1f}%)")

# 3.2 Máscara de calidad del composite (≥ MIN_OBS observaciones)
quality_mask = n_obs >= MIN_OBS
n_quality = np.sum(valid_mask & quality_mask)
log(f"  Píxeles con calidad suficiente (≥{MIN_OBS} obs): {n_quality:,} ({100*n_quality/max(n_valid,1):.1f}% de válidos)")

# 3.3 Máscara de pendiente (≥ 15° según Notti et al. 2023)
slope_mask = slope >= SLOPE_MIN_DEG
n_slope = np.sum(valid_mask & slope_mask)
log(f"  Píxeles con pendiente ≥ {SLOPE_MIN_DEG}°: {n_slope:,} ({100*n_slope/max(n_valid,1):.1f}% de válidos)")

# 3.4 Máscara de exclusión (agua + urbano de WorldCover)
exclusion_mask = np.zeros_like(worldcover, dtype=bool)
for cls in EXCLUDE_CLASSES:
    n_cls = np.sum(worldcover == cls)
    exclusion_mask |= (worldcover == cls)
    cls_name = {50: "urbano", 80: "agua"}.get(cls, str(cls))
    log(f"  Píxeles {cls_name} (WorldCover {cls}): {n_cls:,}")

land_mask = ~exclusion_mask
n_land = np.sum(valid_mask & land_mask)
log(f"  Píxeles de tierra (excluidos agua+urbano): {n_land:,}")

# 3.5 Máscara combinada base
base_mask = valid_mask & quality_mask & slope_mask & land_mask
n_base = np.sum(base_mask)
log(f"  Píxeles que pasan todos los filtros base: {n_base:,}")

log("")

# ============================================================================
# SECCIÓN 4: DETECCIÓN MULTI-CRITERIO
# ============================================================================

log("--- SECCIÓN 4: DETECCIÓN MULTI-CRITERIO ---")

# 4.1 Criterio principal: dNDVI < umbral (protocolo Notti et al. 2023)
dndvi_detect = base_mask & (dNDVI < DNDVI_THRESHOLD)
n_dndvi = np.sum(dndvi_detect)
log(f"  Criterio A (dNDVI < {DNDVI_THRESHOLD}): {n_dndvi:,} píxeles")

# 4.2 Criterio secundario B: dNBR > umbral
dnbr_detect = base_mask & np.isfinite(dNBR) & (dNBR > DNBR_THRESHOLD)
n_dnbr = np.sum(dnbr_detect)
log(f"  Criterio B (dNBR > {DNBR_THRESHOLD}): {n_dnbr:,} píxeles")

# 4.3 Criterio secundario C: dBSI > umbral
dbsi_detect = base_mask & np.isfinite(dBSI) & (dBSI > DBSI_THRESHOLD)
n_dbsi = np.sum(dbsi_detect)
log(f"  Criterio C (dBSI > {DBSI_THRESHOLD}): {n_dbsi:,} píxeles")

# 4.4 Candidatos: OBLIGATORIO cumplir criterio A (dNDVI + pendiente)
# Score de evidencia: cuántos criterios se cumplen (A es obligatorio)
candidates = dndvi_detect.copy()
score_map = np.zeros(ref_shape, dtype=np.int8)
score_map[dndvi_detect] += 1   # A cumplido
score_map[dnbr_detect] += 1    # B cumplido
score_map[dbsi_detect] += 1    # C cumplido

# Solo donde A se cumple (dNDVI < umbral AND pendiente ≥ 15°)
score_map[~dndvi_detect] = 0

n_score1 = np.sum(score_map == 1)  # solo dNDVI
n_score2 = np.sum(score_map == 2)  # dNDVI + 1 secundario
n_score3 = np.sum(score_map == 3)  # los 3 índices

log(f"")
log(f"  Píxeles candidatos (A obligatorio):")
log(f"    Score 1/3 (solo dNDVI+pendiente): {n_score1:,}")
log(f"    Score 2/3 (dNDVI + 1 secundario): {n_score2:,}")
log(f"    Score 3/3 (los 3 índices):        {n_score3:,}")
log(f"    TOTAL CANDIDATOS:                 {n_dndvi:,}")

log("")

# ============================================================================
# SECCIÓN 5: VECTORIZACIÓN Y ATRIBUTOS
# ============================================================================

log("--- SECCIÓN 5: VECTORIZACIÓN Y ATRIBUTOS ---")

# 5.1 Etiquetar componentes conectados
from scipy.ndimage import label as ndlabel

labeled, n_features = ndlabel(candidates.astype(np.int32))
log(f"  Componentes conectados encontrados: {n_features}")

# 5.2 Filtrar por tamaño mínimo
# Contar píxeles por etiqueta
labels_flat = labeled.ravel()
counts = np.bincount(labels_flat)
# Máscara de etiquetas que superan el mínimo (excluir etiqueta 0 = fondo)
valid_labels = np.where(counts >= MIN_PIXELS)[0]
valid_labels = valid_labels[valid_labels > 0]  # excluir fondo
n_valid_clusters = len(valid_labels)
log(f"  Clusters con ≥ {MIN_PIXELS} píxeles (≥ {MIN_PIXELS*100/10000:.2f} ha): {n_valid_clusters}")

# Crear máscara filtrada
filtered = np.isin(labeled, valid_labels)
n_filtered_px = np.sum(filtered)
log(f"  Píxeles en clusters válidos: {n_filtered_px:,}")

# 5.3 Vectorizar
log(f"  Vectorizando...")

polygons = []
values = []
for geom_dict, val in shapes(labeled.astype(np.int32), mask=filtered, transform=ref_transform):
    if val > 0:
        polygons.append(shape(geom_dict))
        values.append(int(val))

log(f"  Polígonos extraídos: {len(polygons)}")

if len(polygons) == 0:
    log("")
    log("ALERTA: No se detectaron candidatos con los umbrales actuales.")
    log("Considerar ajustar DNDVI_THRESHOLD o SLOPE_MIN_DEG.")
    # Guardar log y salir
    log_path = os.path.join(DIR_OUTPUT, "P07c_log.txt")
    with open(log_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(log_lines))
    sys.exit(0)

# 5.4 Disolver polígonos con la misma etiqueta (multi-parte → uni-parte)
log(f"  Disolviendo polígonos por cluster_id...")
cluster_polys = {}
for poly, val in zip(polygons, values):
    if val not in cluster_polys:
        cluster_polys[val] = []
    cluster_polys[val].append(poly)

dissolved = []
cluster_ids = []
for cid, polys in cluster_polys.items():
    merged = unary_union(polys)
    dissolved.append(merged)
    cluster_ids.append(cid)

log(f"  Clusters disueltos: {len(dissolved)}")

# 5.5 Calcular atributos
log(f"  Calculando atributos...")

records = []
for i, (cid, geom) in enumerate(zip(cluster_ids, dissolved)):
    centroid = geom.centroid

    # Máscara del cluster en el raster
    cluster_mask = (labeled == cid)
    n_px = np.sum(cluster_mask)

    # Área en hectáreas
    area_m2 = n_px * 100  # 10m × 10m = 100 m²
    area_ha = area_m2 / 10000

    # Estadísticas espectrales dentro del cluster
    dndvi_vals = dNDVI[cluster_mask]
    dndvi_vals = dndvi_vals[np.isfinite(dndvi_vals)]
    dndvi_mean = np.mean(dndvi_vals) if len(dndvi_vals) > 0 else np.nan
    dndvi_min = np.min(dndvi_vals) if len(dndvi_vals) > 0 else np.nan

    dnbr_vals = dNBR[cluster_mask]
    dnbr_vals = dnbr_vals[np.isfinite(dnbr_vals)]
    dnbr_mean = np.mean(dnbr_vals) if len(dnbr_vals) > 0 else np.nan

    dbsi_vals = dBSI[cluster_mask]
    dbsi_vals = dbsi_vals[np.isfinite(dbsi_vals)]
    dbsi_mean = np.mean(dbsi_vals) if len(dbsi_vals) > 0 else np.nan

    ndvi_pre_vals = NDVI_pre[cluster_mask]
    ndvi_pre_vals = ndvi_pre_vals[np.isfinite(ndvi_pre_vals)]
    ndvi_pre_mean = np.mean(ndvi_pre_vals) if len(ndvi_pre_vals) > 0 else np.nan

    slope_vals = slope[cluster_mask]
    slope_vals = slope_vals[np.isfinite(slope_vals)]
    slope_mean = np.mean(slope_vals) if len(slope_vals) > 0 else np.nan

    # Score de evidencia espectral (media del cluster)
    score = 0
    if dndvi_mean < DNDVI_THRESHOLD:
        score += 1
    if dnbr_mean > DNBR_THRESHOLD:
        score += 1
    if dbsi_mean > DBSI_THRESHOLD:
        score += 1

    # Clasificación
    if score == 3:
        evidencia = "fuerte"
    elif score == 2:
        evidencia = "moderada"
    else:
        evidencia = "debil"

    # Cobertura pre-evento dominante
    wc_vals = worldcover[cluster_mask]
    wc_vals = wc_vals[wc_vals > 0]
    if len(wc_vals) > 0:
        wc_dominant = int(np.bincount(wc_vals.astype(int)).argmax())
    else:
        wc_dominant = 0

    wc_names = {10:"bosque", 20:"matorral", 30:"pastizal", 40:"cultivo",
                60:"suelo_desnudo", 70:"nieve", 90:"humedal", 95:"manglar"}
    wc_name = wc_names.get(wc_dominant, f"clase_{wc_dominant}")

    records.append({
        "opt_id": i + 1,
        "n_pixels": n_px,
        "area_ha": round(area_ha, 3),
        "centroid_x": round(centroid.x, 1),
        "centroid_y": round(centroid.y, 1),
        "dNDVI_mean": round(dndvi_mean, 4),
        "dNDVI_min": round(dndvi_min, 4),
        "dNBR_mean": round(dnbr_mean, 4),
        "dBSI_mean": round(dbsi_mean, 4),
        "NDVI_pre_mean": round(ndvi_pre_mean, 4),
        "slope_mean_deg": round(slope_mean, 1),
        "score_espectral": score,
        "evidencia": evidencia,
        "cobertura_pre": wc_name,
        "worldcover_class": wc_dominant,
    })

gdf = gpd.GeoDataFrame(records, geometry=dissolved, crs=ref_crs)
log(f"  Inventario construido: {len(gdf)} candidatos ópticos")

# Resumen por evidencia
for ev in ["fuerte", "moderada", "debil"]:
    n = len(gdf[gdf["evidencia"] == ev])
    log(f"    Evidencia {ev}: {n}")

log("")

# ============================================================================
# SECCIÓN 6: ANÁLISIS DE SENSIBILIDAD
# ============================================================================

log("--- SECCIÓN 6: ANÁLISIS DE SENSIBILIDAD ---")
log(f"  Umbrales evaluados: {SENSITIVITY_THRESHOLDS}")

sensitivity_records = []
for thresh in SENSITIVITY_THRESHOLDS:
    detect = base_mask & (dNDVI < thresh)
    n_px = np.sum(detect)

    # Contar clusters
    lab_s, n_feat_s = ndlabel(detect.astype(np.int32))
    counts_s = np.bincount(lab_s.ravel())
    valid_s = np.sum(counts_s[1:] >= MIN_PIXELS)

    sensitivity_records.append({
        "dNDVI_threshold": thresh,
        "pixels_detected": int(n_px),
        "area_km2": round(n_px * 100 / 1e6, 2),
        "clusters_total": int(n_feat_s),
        "clusters_valid": int(valid_s),
    })
    log(f"    dNDVI < {thresh}: {n_px:,} px, {valid_s} clusters válidos")

df_sens = pd.DataFrame(sensitivity_records)
log("")

# ============================================================================
# SECCIÓN 7: ESTADÍSTICAS RESUMEN
# ============================================================================

log("--- SECCIÓN 7: ESTADÍSTICAS RESUMEN ---")

log(f"  Total candidatos ópticos: {len(gdf)}")
log(f"  Área total: {gdf['area_ha'].sum():.1f} ha ({gdf['area_ha'].sum()/100:.2f} km²)")
log(f"  Área media: {gdf['area_ha'].mean():.2f} ha")
log(f"  Área mediana: {gdf['area_ha'].median():.2f} ha")
log(f"  Rango de área: [{gdf['area_ha'].min():.3f}, {gdf['area_ha'].max():.2f}] ha")
log(f"  dNDVI medio: {gdf['dNDVI_mean'].mean():.4f}")
log(f"  Pendiente media: {gdf['slope_mean_deg'].mean():.1f}°")
log(f"  NDVI pre-evento medio: {gdf['NDVI_pre_mean'].mean():.4f}")

# Distribución por cobertura pre-evento
log(f"  Distribución por cobertura pre-evento:")
for cob, n in gdf['cobertura_pre'].value_counts().items():
    log(f"    {cob}: {n} ({100*n/len(gdf):.1f}%)")

log("")

# ============================================================================
# SECCIÓN 8: EXPORTACIÓN
# ============================================================================

log("--- SECCIÓN 8: EXPORTACIÓN ---")

# 8.1 GeoPackage
gpkg_path = os.path.join(DIR_OUTPUT, "P07c_inventario_optico.gpkg")
gdf.to_file(gpkg_path, driver="GPKG")
log(f"  Exportado: {gpkg_path}")

# 8.2 CSV
csv_path = os.path.join(DIR_OUTPUT, "P07c_inventario_optico.csv")
gdf.drop(columns='geometry').to_csv(csv_path, index=False)
log(f"  Exportado: {csv_path}")

# 8.3 Sensibilidad
sens_path = os.path.join(DIR_OUTPUT, "P07c_sensibilidad_umbrales.csv")
df_sens.to_csv(sens_path, index=False)
log(f"  Exportado: {sens_path}")

# 8.4 Raster binario de detección (para QGIS)
raster_path = os.path.join(DIR_OUTPUT, "P07c_deteccion_binaria.tif")
out_profile = ref_profile.copy()
out_profile.update(dtype='uint8', count=1, nodata=255, compress='LZW')
out_raster = np.zeros(ref_shape, dtype=np.uint8)
out_raster[filtered] = score_map[filtered]
out_raster[~valid_mask] = 255
with rasterio.open(raster_path, 'w', **out_profile) as dst:
    dst.write(out_raster, 1)
log(f"  Exportado: {raster_path}")

log("")

# ============================================================================
# SECCIÓN 9: VISUALIZACIONES
# ============================================================================

log("--- SECCIÓN 9: VISUALIZACIONES ---")

# 9.1 Mapa del inventario
log(f"  Generando mapa del inventario...")
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
colors = {"fuerte": "#d62728", "moderada": "#ff7f0e", "debil": "#2ca02c"}
for ev in ["debil", "moderada", "fuerte"]:
    subset = gdf[gdf["evidencia"] == ev]
    if len(subset) > 0:
        subset.plot(ax=ax, color=colors[ev], alpha=0.7,
                    label=f"{ev.capitalize()} ({len(subset)})", edgecolor='black', linewidth=0.3)

ax.set_title(f"P07c — Inventario Óptico Independiente\n"
             f"dNDVI < {DNDVI_THRESHOLD}, pendiente ≥ {SLOPE_MIN_DEG}° | "
             f"n = {len(gdf)} candidatos", fontsize=13)
ax.set_xlabel("Este (m) — EPSG:32617")
ax.set_ylabel("Norte (m)")
ax.legend(title="Evidencia espectral", loc="upper right")
ax.ticklabel_format(style='plain')
plt.tight_layout()
map_path = os.path.join(DIR_OUTPUT, "P07c_mapa_inventario.png")
plt.savefig(map_path, dpi=200)
plt.close()
log(f"  Exportado: {map_path}")

# 9.2 Histogramas de atributos
log(f"  Generando histogramas de atributos...")
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

# Área
axes[0,0].hist(gdf['area_ha'], bins=50, color='steelblue', edgecolor='black', linewidth=0.5)
axes[0,0].set_xlabel('Área (ha)')
axes[0,0].set_ylabel('Frecuencia')
axes[0,0].set_title('Distribución de área')
axes[0,0].axvline(gdf['area_ha'].median(), color='red', linestyle='--', label=f'Mediana={gdf["area_ha"].median():.2f} ha')
axes[0,0].legend(fontsize=8)

# dNDVI
axes[0,1].hist(gdf['dNDVI_mean'], bins=50, color='darkgreen', edgecolor='black', linewidth=0.5)
axes[0,1].set_xlabel('dNDVI medio')
axes[0,1].set_title('dNDVI en candidatos')
axes[0,1].axvline(DNDVI_THRESHOLD, color='red', linestyle='--', label=f'Umbral={DNDVI_THRESHOLD}')
axes[0,1].legend(fontsize=8)

# Pendiente
axes[0,2].hist(gdf['slope_mean_deg'], bins=50, color='sienna', edgecolor='black', linewidth=0.5)
axes[0,2].set_xlabel('Pendiente (°)')
axes[0,2].set_title('Pendiente media')
axes[0,2].axvline(SLOPE_MIN_DEG, color='red', linestyle='--', label=f'Mínimo={SLOPE_MIN_DEG}°')
axes[0,2].legend(fontsize=8)

# dNBR
axes[1,0].hist(gdf['dNBR_mean'].dropna(), bins=50, color='darkorange', edgecolor='black', linewidth=0.5)
axes[1,0].set_xlabel('dNBR medio')
axes[1,0].set_title('dNBR en candidatos')
axes[1,0].axvline(DNBR_THRESHOLD, color='red', linestyle='--', label=f'Umbral={DNBR_THRESHOLD}')
axes[1,0].legend(fontsize=8)

# dBSI
axes[1,1].hist(gdf['dBSI_mean'].dropna(), bins=50, color='goldenrod', edgecolor='black', linewidth=0.5)
axes[1,1].set_xlabel('dBSI medio')
axes[1,1].set_title('dBSI en candidatos')
axes[1,1].axvline(DBSI_THRESHOLD, color='red', linestyle='--', label=f'Umbral={DBSI_THRESHOLD}')
axes[1,1].legend(fontsize=8)

# Sensibilidad
axes[1,2].bar([str(t) for t in df_sens['dNDVI_threshold']],
              df_sens['clusters_valid'], color='mediumpurple', edgecolor='black')
axes[1,2].set_xlabel('Umbral dNDVI')
axes[1,2].set_ylabel('Clusters válidos')
axes[1,2].set_title('Sensibilidad al umbral')
axes[1,2].axhline(len(gdf), color='red', linestyle='--', label=f'Seleccionado: {DNDVI_THRESHOLD}')
axes[1,2].legend(fontsize=8)

plt.suptitle(f"P07c — Atributos del inventario óptico (n={len(gdf)})", fontsize=14, y=1.02)
plt.tight_layout()
hist_path = os.path.join(DIR_OUTPUT, "P07c_histograma_atributos.png")
plt.savefig(hist_path, dpi=200, bbox_inches='tight')
plt.close()
log(f"  Exportado: {hist_path}")

# 9.3 NDVI pre-evento vs dNDVI
log(f"  Generando scatter NDVI_pre vs dNDVI...")
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(gdf['NDVI_pre_mean'], gdf['dNDVI_mean'],
                     c=gdf['score_espectral'], cmap='RdYlGn_r',
                     s=gdf['area_ha']*5, alpha=0.6, edgecolors='black', linewidths=0.3,
                     vmin=1, vmax=3)
ax.set_xlabel('NDVI pre-evento')
ax.set_ylabel('dNDVI (post - pre)')
ax.set_title('P07c — NDVI pre-evento vs pérdida de vegetación')
ax.axhline(DNDVI_THRESHOLD, color='red', linestyle='--', alpha=0.5)
cbar = plt.colorbar(scatter, ax=ax, label='Score espectral (1-3)')
cbar.set_ticks([1, 2, 3])
plt.tight_layout()
scatter_path = os.path.join(DIR_OUTPUT, "P07c_scatter_ndvi_pre_vs_dndvi.png")
plt.savefig(scatter_path, dpi=200)
plt.close()
log(f"  Exportado: {scatter_path}")

log("")

# ============================================================================
# SECCIÓN 10: RESUMEN FINAL
# ============================================================================

log("=" * 70)
log("RESUMEN FINAL P07c — INVENTARIO ÓPTICO INDEPENDIENTE")
log("=" * 70)
log(f"  Umbral dNDVI:              {DNDVI_THRESHOLD} (Notti et al. 2023)")
log(f"  Pendiente mínima:          {SLOPE_MIN_DEG}° (protocolo)")
log(f"  Área mínima:               {MIN_PIXELS} px = {MIN_PIXELS*100/10000:.2f} ha")
log(f"  Exclusiones:               agua (WC 80) + urbano (WC 50)")
log(f"")
log(f"  CANDIDATOS ÓPTICOS:        {len(gdf)}")
log(f"    Evidencia fuerte (3/3):  {len(gdf[gdf['evidencia']=='fuerte'])}")
log(f"    Evidencia moderada (2/3):{len(gdf[gdf['evidencia']=='moderada'])}")
log(f"    Evidencia débil (1/3):   {len(gdf[gdf['evidencia']=='debil'])}")
log(f"")
log(f"  Directorio de salida: {DIR_OUTPUT}")
log("=" * 70)

# Guardar log
log_path = os.path.join(DIR_OUTPUT, "P07c_log.txt")
with open(log_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines))
log(f"")
log(f"Log guardado: {log_path}")
log(f"")
log(f"P07c completado exitosamente.")
log(f"")
log(f"SIGUIENTE PASO: Paso 1.2 — Interpretación visual en Google Earth Pro")
log(f"LUEGO: Paso 1.3 — Cruzar inventario óptico (P07c) con inventario InSAR (P06b)")