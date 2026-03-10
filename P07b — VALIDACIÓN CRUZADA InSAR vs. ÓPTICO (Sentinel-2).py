# =============================================================================
# P07b — VALIDACIÓN CRUZADA InSAR vs. ÓPTICO (Sentinel-2)
# Proyecto: Detección de deslizamientos por InSAR — Intag, Cotacachi, Ecuador
# Evento: 19 de diciembre de 2023
# Autor: Víctor Pinto
# Fecha: 28 de febrero de 2026
# Python: 3.x / Spyder IDE
# =============================================================================
#
# OBJETIVO:
# Validar independientemente los clústeres de deformación detectados por
# DInSAR (Sentinel-1) mediante cruce espacial con índices de cambio óptico
# derivados de composites Sentinel-2 pre/post evento (P07 GEE Cloud Score+).
#
# FUNDAMENTO CIENTÍFICO:
# Un deslizamiento causa simultáneamente:
#   - Desplazamiento del terreno → detectable por InSAR (deformación LOS)
#   - Remoción de vegetación → detectable por dNDVI negativo
#   - Exposición de suelo desnudo → detectable por dNBR positivo, dBSI positivo
# La concordancia espacial entre ambas fuentes independientes constituye
# evidencia convergente que fortalece la identificación (Xu et al., 2024;
# Dini et al., 2026).
#
# INPUTS:
#   1. Centroides de clústeres InSAR (GeoPackage, de P06 v2.0)
#   2. Clústeres InSAR como polígonos (GeoPackage, de P06 v2.0)
#   3. Catálogo de verificación con scores de confianza (de P06b, si existe)
#   4. Rasters Sentinel-2 de P07 GEE:
#      - P07_dNDVI.tif  (NDVI_post - NDVI_pre; negativo = pérdida vegetación)
#      - P07_dNBR.tif   (NBR_pre - NBR_post; positivo = degradación)
#      - P07_dBSI.tif   (BSI_post - BSI_pre; positivo = más suelo desnudo)
#      - P07_n_obs.tif  (observaciones válidas por píxel)
#      - P07_NDVI_pre.tif  (NDVI del composite pre-evento)
#      - P07_NDVI_post.tif (NDVI del composite post-evento)
#
# OUTPUTS:
#   - P07b_validacion_cruzada.gpkg  (centroides con valores ópticos)
#   - P07b_validacion_cruzada.csv   (tabla completa)
#   - P07b_resumen_validacion.csv   (estadísticas agregadas)
#   - P07b_mapa_concordancia.png    (mapa de concordancia espacial)
#   - P07b_histogramas.png          (distribución de índices en clústeres)
#   - P07b_scatter_deformacion_vs_dNDVI.png
#   - P07b_log_validacion.txt       (log completo para trazabilidad)
#
# CONVENCIONES DE SIGNO (verificadas empíricamente en el script):
#   dNDVI < 0 → pérdida de vegetación (cicatriz de deslizamiento)
#   dNBR  > 0 → degradación/perturbación del terreno
#   dBSI  > 0 → aumento de suelo desnudo expuesto
#
# CRITERIO DE CONCORDANCIA:
#   Un clúster InSAR se considera "validado ópticamente" si cumple AL MENOS
#   2 de 3 condiciones:
#     (a) dNDVI mediano < -0.05 (pérdida significativa de vegetación)
#     (b) dNBR mediano > 0.05  (degradación significativa)
#     (c) dBSI mediano > 0.02  (aumento de suelo desnudo)
#   Umbrales basados en literatura (Xu et al., 2024; Behling et al., 2014).
# =============================================================================

import os
import sys
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask as rio_mask
from shapely.geometry import mapping
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# SECCIÓN 0: CONFIGURACIÓN DE RUTAS Y PARÁMETROS
# =============================================================================

# --- Rutas principales ---
DIR_INSAR = r"D:\POSGRADOS\INTAG\data\sentinel1\hyp3_deformation"
DIR_OPTICO = r"D:\POSGRADOS\INTAG\data\sentinel2\VALIDACION\INTAG_InSAR_P07"
DIR_VERIF = os.path.join(DIR_INSAR, "verificacion")
DIR_OUTPUT = os.path.join(DIR_OPTICO, "validacion_cruzada")

# Crear directorio de salida
os.makedirs(DIR_OUTPUT, exist_ok=True)

# --- Dataset principal para validación ---
# M_ASC_CO_12d = mosaico ascendente co-evento 12 días (el más relevante)
DATASET_PRINCIPAL = "M_ASC_CO_12d"

# --- Todos los datasets disponibles ---
DATASETS_TODOS = [
    "M_ASC_CO_12d", "M_ASC_CO_24d",
    "D1", "D2", "D3", "D4",
    "A2", "A4",
    "P18_B1", "P18_B3"
]

# --- Archivos ópticos P07 ---
RASTERS_OPTICOS = {
    "dNDVI": os.path.join(DIR_OPTICO, "P07_dNDVI.tif"),
    "dNBR": os.path.join(DIR_OPTICO, "P07_dNBR.tif"),
    "dBSI": os.path.join(DIR_OPTICO, "P07_dBSI.tif"),
    "n_obs": os.path.join(DIR_OPTICO, "P07_n_obs.tif"),
    "NDVI_pre": os.path.join(DIR_OPTICO, "P07_NDVI_pre.tif"),
    "NDVI_post": os.path.join(DIR_OPTICO, "P07_NDVI_post.tif"),
}

# --- Umbrales de concordancia óptica ---
# Basados en literatura científica para detección de cicatrices de deslizamiento
UMBRAL_dNDVI = -0.05   # dNDVI < este valor = pérdida significativa de vegetación
UMBRAL_dNBR = 0.05     # dNBR > este valor = degradación significativa
UMBRAL_dBSI = 0.02     # dBSI > este valor = aumento de suelo desnudo
MIN_CRITERIOS = 2       # Mínimo de criterios ópticos para validar un clúster
MIN_OBS_CALIDAD = 3     # Mínimo de observaciones S2 válidas para considerar fiable

# --- Configuración de visualización ---
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# --- Log ---
LOG_PATH = os.path.join(DIR_OUTPUT, "P07b_log_validacion.txt")
log_lines = []

def log(msg):
    # Registra mensaje en consola y en archivo de log
    timestamp = datetime.now().strftime("%H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    log_lines.append(line)

log("=" * 70)
log("P07b — VALIDACIÓN CRUZADA InSAR vs. ÓPTICO (Sentinel-2)")
log(f"Fecha de ejecución: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

# =============================================================================
# SECCIÓN 1: VERIFICACIÓN DE INPUTS
# =============================================================================

log("\n--- SECCIÓN 1: VERIFICACIÓN DE INPUTS ---")

# 1.1 Verificar rasters ópticos
log("\n[1.1] Verificando rasters ópticos P07...")
rasters_ok = True
for nombre, ruta in RASTERS_OPTICOS.items():
    if os.path.exists(ruta):
        with rasterio.open(ruta) as src:
            log(f"  OK: {nombre} — {src.width}x{src.height} px, "
                f"CRS={src.crs}, res={src.res[0]:.1f}m")
    else:
        log(f"  FALTA: {nombre} — {ruta}")
        rasters_ok = False

if not rasters_ok:
    log("ERROR: Faltan rasters ópticos. Ejecutar P07 en GEE primero.")
    sys.exit(1)

# 1.2 Verificar archivos InSAR
log("\n[1.2] Verificando archivos InSAR P06...")

# Buscar centroides del dataset principal
centroids_path = os.path.join(DIR_INSAR, f"centroids_{DATASET_PRINCIPAL}.gpkg")
clusters_path = os.path.join(DIR_INSAR, f"clusters_{DATASET_PRINCIPAL}.gpkg")

if os.path.exists(centroids_path):
    gdf_centroids = gpd.read_file(centroids_path)
    log(f"  OK: Centroides {DATASET_PRINCIPAL} — {len(gdf_centroids)} clústeres")
else:
    log(f"  FALTA: {centroids_path}")
    sys.exit(1)

if os.path.exists(clusters_path):
    gdf_clusters = gpd.read_file(clusters_path)
    log(f"  OK: Clústeres {DATASET_PRINCIPAL} — {len(gdf_clusters)} polígonos")
else:
    log(f"  FALTA: {clusters_path}")
    sys.exit(1)

# 1.3 Buscar catálogo de verificación P06b (scores de confianza)
log("\n[1.3] Buscando catálogo de verificación P06b...")
catalogo_path = None
posibles_catalogo = [
    os.path.join(DIR_VERIF, "catalogo_verificado.gpkg"),
    os.path.join(DIR_VERIF, "catalogo_unificado.gpkg"),
    os.path.join(DIR_VERIF, "P06b_catalogo_verificado.gpkg"),
    os.path.join(DIR_VERIF, "clusters_verificados.gpkg"),
]

# Buscar también cualquier gpkg en la carpeta verificacion
if os.path.isdir(DIR_VERIF):
    gpkgs_verif = glob.glob(os.path.join(DIR_VERIF, "*.gpkg"))
    posibles_catalogo.extend(gpkgs_verif)

for ruta in posibles_catalogo:
    if os.path.exists(ruta):
        catalogo_path = ruta
        gdf_catalogo = gpd.read_file(catalogo_path)
        log(f"  OK: Catálogo P06b encontrado — {catalogo_path}")
        log(f"      {len(gdf_catalogo)} registros, columnas: {list(gdf_catalogo.columns)}")
        break

if catalogo_path is None:
    log("  INFO: No se encontró catálogo P06b. Se usarán centroides de P06 directamente.")
    log("        (La columna 'confidence_score' no estará disponible)")

# 1.4 Verificar compatibilidad CRS
log("\n[1.4] Verificando compatibilidad de CRS...")
with rasterio.open(RASTERS_OPTICOS["dNDVI"]) as src:
    crs_optico = src.crs

crs_insar = gdf_centroids.crs
log(f"  CRS óptico (S2): {crs_optico}")
log(f"  CRS InSAR (S1):  {crs_insar}")

if crs_optico != crs_insar:
    log("  Reproyectando centroides InSAR al CRS del raster óptico...")
    gdf_centroids = gdf_centroids.to_crs(crs_optico)
    gdf_clusters = gdf_clusters.to_crs(crs_optico)
    if catalogo_path:
        gdf_catalogo = gdf_catalogo.to_crs(crs_optico)
    log(f"  Reproyección completada a {crs_optico}")
else:
    log("  CRS compatibles — sin reproyección necesaria")

# =============================================================================
# SECCIÓN 2: EXTRACCIÓN DE VALORES ÓPTICOS EN CENTROIDES
# =============================================================================

log("\n--- SECCIÓN 2: EXTRACCIÓN DE VALORES ÓPTICOS EN CENTROIDES ---")

# 2.1 Extracción puntual en centroides (valor del píxel donde cae el centroide)
log("\n[2.1] Extrayendo valores ópticos puntuales en centroides...")

# Usar el dataframe base: catálogo P06b si existe, si no centroides P06
if catalogo_path:
    gdf_base = gdf_catalogo.copy()
    log(f"  Base: catálogo P06b ({len(gdf_base)} clústeres con scores de confianza)")
else:
    gdf_base = gdf_centroids.copy()
    log(f"  Base: centroides P06 ({len(gdf_base)} clústeres)")

# Extraer coordenadas de centroides
if gdf_base.geometry.geom_type.iloc[0] == 'Point':
    coords = [(geom.x, geom.y) for geom in gdf_base.geometry]
else:
    # Si son polígonos, calcular centroides
    coords = [(geom.centroid.x, geom.centroid.y) for geom in gdf_base.geometry]

# Extraer valor de cada raster en cada centroide
for nombre_idx, ruta_raster in RASTERS_OPTICOS.items():
    valores = []
    with rasterio.open(ruta_raster) as src:
        for x, y in coords:
            try:
                # Convertir coordenada a fila/columna del raster
                row, col = src.index(x, y)
                # Verificar que está dentro de los límites
                if 0 <= row < src.height and 0 <= col < src.width:
                    val = src.read(1)[row, col]
                    # Verificar NoData
                    if src.nodata is not None and val == src.nodata:
                        valores.append(np.nan)
                    else:
                        valores.append(float(val))
                else:
                    valores.append(np.nan)
            except Exception:
                valores.append(np.nan)

    gdf_base[f"opt_{nombre_idx}"] = valores
    n_valid = np.sum(~np.isnan(valores))
    log(f"  {nombre_idx}: {n_valid}/{len(valores)} centroides con valor válido "
        f"(media={np.nanmean(valores):.4f}, rango=[{np.nanmin(valores):.4f}, "
        f"{np.nanmax(valores):.4f}])")

# =============================================================================
# SECCIÓN 3: ESTADÍSTICAS ZONALES EN POLÍGONOS DE CLÚSTERES
# =============================================================================

log("\n--- SECCIÓN 3: ESTADÍSTICAS ZONALES EN POLÍGONOS ---")
log(f"Procesando {len(gdf_clusters)} polígonos de clústeres...")

# Preparar columnas para estadísticas zonales
indices_zonales = ["dNDVI", "dNBR", "dBSI"]
stats_nombres = ["median", "mean", "std", "min", "max", "n_pixels"]

for idx_name in indices_zonales:
    for stat in stats_nombres:
        gdf_clusters[f"zonal_{idx_name}_{stat}"] = np.nan

# Extraer estadísticas zonales para cada clúster
n_procesados = 0
n_errores = 0

for i, row in gdf_clusters.iterrows():
    geom = row.geometry
    geom_json = [mapping(geom)]

    for idx_name in indices_zonales:
        ruta_raster = RASTERS_OPTICOS[idx_name]
        try:
            with rasterio.open(ruta_raster) as src:
                out_image, out_transform = rio_mask(src, geom_json, crop=True,
                                                     nodata=np.nan, filled=True)
                datos = out_image[0]  # Banda 1
                # Filtrar NoData
                if src.nodata is not None:
                    datos = np.where(datos == src.nodata, np.nan, datos)
                datos_valid = datos[~np.isnan(datos)]

                if len(datos_valid) > 0:
                    gdf_clusters.at[i, f"zonal_{idx_name}_median"] = np.median(datos_valid)
                    gdf_clusters.at[i, f"zonal_{idx_name}_mean"] = np.mean(datos_valid)
                    gdf_clusters.at[i, f"zonal_{idx_name}_std"] = np.std(datos_valid)
                    gdf_clusters.at[i, f"zonal_{idx_name}_min"] = np.min(datos_valid)
                    gdf_clusters.at[i, f"zonal_{idx_name}_max"] = np.max(datos_valid)
                    gdf_clusters.at[i, f"zonal_{idx_name}_n_pixels"] = len(datos_valid)
        except Exception as e:
            n_errores += 1

    n_procesados += 1
    if n_procesados % 200 == 0:
        log(f"  Procesados {n_procesados}/{len(gdf_clusters)} clústeres...")

log(f"  Completado: {n_procesados} clústeres procesados, {n_errores} errores")

# Estadísticas resumen de las zonales
for idx_name in indices_zonales:
    col = f"zonal_{idx_name}_median"
    vals = gdf_clusters[col].dropna()
    log(f"  {idx_name} zonal mediana: n={len(vals)}, "
        f"media={vals.mean():.4f}, rango=[{vals.min():.4f}, {vals.max():.4f}]")

# =============================================================================
# SECCIÓN 4: VERIFICACIÓN DE CONVENCIÓN DE SIGNOS
# =============================================================================

log("\n--- SECCIÓN 4: VERIFICACIÓN EMPÍRICA DE CONVENCIONES DE SIGNO ---")

# Si hay deformación real (deslizamientos), esperamos:
#   dNDVI negativo (pérdida veg), dNBR positivo (degradación), dBSI positivo (suelo)
# Verificamos que la mayoría de clústeres con deformación fuerte cumplan esto

# Usar valores puntuales del centroide
dndvi_vals = gdf_base["opt_dNDVI"].dropna()
dnbr_vals = gdf_base["opt_dNBR"].dropna()
dbsi_vals = gdf_base["opt_dBSI"].dropna()

pct_dndvi_neg = (dndvi_vals < 0).sum() / len(dndvi_vals) * 100 if len(dndvi_vals) > 0 else 0
pct_dnbr_pos = (dnbr_vals > 0).sum() / len(dnbr_vals) * 100 if len(dnbr_vals) > 0 else 0
pct_dbsi_pos = (dbsi_vals > 0).sum() / len(dbsi_vals) * 100 if len(dbsi_vals) > 0 else 0

log(f"  dNDVI < 0 (pérdida veg):    {pct_dndvi_neg:.1f}% de clústeres")
log(f"  dNBR  > 0 (degradación):    {pct_dnbr_pos:.1f}% de clústeres")
log(f"  dBSI  > 0 (suelo desnudo):  {pct_dbsi_pos:.1f}% de clústeres")

# Evaluar si la convención es coherente
convencion_ok = True
if pct_dndvi_neg < 40:
    log("  ALERTA: Mayoría de dNDVI son positivos — posible convención invertida")
    log("          Esperado: dNDVI = NDVI_post - NDVI_pre (negativo = pérdida)")
    convencion_ok = False
if pct_dnbr_pos < 40:
    log("  ALERTA: Mayoría de dNBR son negativos — posible convención invertida")
    log("          Esperado: dNBR = NBR_pre - NBR_post (positivo = degradación)")
    convencion_ok = False

if convencion_ok:
    log("  CONVENCIONES VERIFICADAS: signos coherentes con deslizamientos")
else:
    log("  REVISAR convenciones de signo del script GEE P07")

# =============================================================================
# SECCIÓN 5: CLASIFICACIÓN DE CONCORDANCIA InSAR-ÓPTICO
# =============================================================================

log("\n--- SECCIÓN 5: CLASIFICACIÓN DE CONCORDANCIA ---")

# 5.1 Clasificar cada centroide
log("\n[5.1] Aplicando criterios de concordancia en centroides...")

# Criterio A: dNDVI < umbral (pérdida significativa de vegetación)
gdf_base["crit_dNDVI"] = gdf_base["opt_dNDVI"] < UMBRAL_dNDVI

# Criterio B: dNBR > umbral (degradación significativa)
gdf_base["crit_dNBR"] = gdf_base["opt_dNBR"] > UMBRAL_dNBR

# Criterio C: dBSI > umbral (aumento de suelo desnudo)
gdf_base["crit_dBSI"] = gdf_base["opt_dBSI"] > UMBRAL_dBSI

# Criterio D: calidad del dato óptico (suficientes observaciones S2)
gdf_base["crit_calidad"] = gdf_base["opt_n_obs"] >= MIN_OBS_CALIDAD

# Contar criterios ópticos cumplidos (A, B, C — sin contar calidad)
gdf_base["n_criterios_opticos"] = (
    gdf_base["crit_dNDVI"].astype(int) +
    gdf_base["crit_dNBR"].astype(int) +
    gdf_base["crit_dBSI"].astype(int)
)

# Clasificación de concordancia
def clasificar_concordancia(row):
    # Si no hay dato óptico válido
    if pd.isna(row["opt_dNDVI"]) and pd.isna(row["opt_dNBR"]):
        return "sin_dato_optico"
    # Si la calidad óptica es baja
    if not row["crit_calidad"]:
        return "baja_calidad_optica"
    # Clasificar por número de criterios
    n = row["n_criterios_opticos"]
    if n >= 3:
        return "concordancia_fuerte"     # 3/3 criterios
    elif n >= MIN_CRITERIOS:
        return "concordancia_moderada"   # 2/3 criterios
    elif n == 1:
        return "concordancia_debil"      # 1/3 criterios
    else:
        return "sin_concordancia"        # 0/3 criterios

gdf_base["concordancia"] = gdf_base.apply(clasificar_concordancia, axis=1)

# Estadísticas de concordancia
log("\n[5.2] Resultados de concordancia InSAR-Óptico:")
log(f"  Umbrales aplicados: dNDVI<{UMBRAL_dNDVI}, dNBR>{UMBRAL_dNBR}, dBSI>{UMBRAL_dBSI}")
log(f"  Mínimo criterios para validación: {MIN_CRITERIOS}/3")
log(f"  Mínimo observaciones S2 para calidad: {MIN_OBS_CALIDAD}")

conteo = gdf_base["concordancia"].value_counts()
total = len(gdf_base)
log(f"\n  Total clústeres analizados: {total}")
for cat in ["concordancia_fuerte", "concordancia_moderada", "concordancia_debil",
            "sin_concordancia", "baja_calidad_optica", "sin_dato_optico"]:
    n = conteo.get(cat, 0)
    pct = n / total * 100
    log(f"  {cat:30s}: {n:5d} ({pct:5.1f}%)")

# Tasa de validación (concordancia fuerte + moderada / total con dato válido)
n_con_dato = total - conteo.get("sin_dato_optico", 0) - conteo.get("baja_calidad_optica", 0)
n_validados = conteo.get("concordancia_fuerte", 0) + conteo.get("concordancia_moderada", 0)
if n_con_dato > 0:
    tasa_validacion = n_validados / n_con_dato * 100
    log(f"\n  TASA DE VALIDACIÓN ÓPTICA: {n_validados}/{n_con_dato} = {tasa_validacion:.1f}%")
    log(f"  (Clústeres InSAR confirmados por >=2 índices ópticos independientes)")
else:
    tasa_validacion = 0
    log("\n  ADVERTENCIA: No hay suficientes datos ópticos válidos para calcular tasa")

# 5.3 Si existe score de confianza P06b, cruzar con concordancia óptica
if catalogo_path and "confidence_score" in gdf_base.columns:
    log("\n[5.3] Cruce concordancia óptica vs. score de confianza P06b:")
    for score in sorted(gdf_base["confidence_score"].unique()):
        subset = gdf_base[gdf_base["confidence_score"] == score]
        n_sub = len(subset)
        n_val = subset["concordancia"].isin(
            ["concordancia_fuerte", "concordancia_moderada"]).sum()
        pct = n_val / n_sub * 100 if n_sub > 0 else 0
        log(f"  Score {score}/4: {n_val}/{n_sub} validados ópticamente ({pct:.1f}%)")
elif "score" in gdf_base.columns:
    log("\n[5.3] Cruce concordancia óptica vs. score de confianza P06b:")
    col_score = "score"
    for score in sorted(gdf_base[col_score].unique()):
        subset = gdf_base[gdf_base[col_score] == score]
        n_sub = len(subset)
        n_val = subset["concordancia"].isin(
            ["concordancia_fuerte", "concordancia_moderada"]).sum()
        pct = n_val / n_sub * 100 if n_sub > 0 else 0
        log(f"  Score {score}/4: {n_val}/{n_sub} validados ópticamente ({pct:.1f}%)")

# =============================================================================
# SECCIÓN 6: EXPORTACIÓN DE RESULTADOS
# =============================================================================

log("\n--- SECCIÓN 6: EXPORTACIÓN DE RESULTADOS ---")

# 6.1 GeoPackage con centroides y valores ópticos
gpkg_path = os.path.join(DIR_OUTPUT, "P07b_validacion_cruzada.gpkg")
gdf_base.to_file(gpkg_path, driver="GPKG")
log(f"  Exportado: {gpkg_path}")

# 6.2 CSV completo
csv_path = os.path.join(DIR_OUTPUT, "P07b_validacion_cruzada.csv")
df_csv = gdf_base.drop(columns=["geometry"]).copy()
df_csv.to_csv(csv_path, index=False, float_format="%.6f")
log(f"  Exportado: {csv_path}")

# 6.3 GeoPackage con polígonos y estadísticas zonales
gpkg_zonas = os.path.join(DIR_OUTPUT, "P07b_clusters_estadisticas_zonales.gpkg")
gdf_clusters.to_file(gpkg_zonas, driver="GPKG")
log(f"  Exportado: {gpkg_zonas}")

# 6.4 Resumen estadístico
resumen = {
    "metrica": [],
    "valor": []
}

resumen["metrica"].append("total_clusters_analizados")
resumen["valor"].append(total)
resumen["metrica"].append("clusters_con_dato_optico_valido")
resumen["valor"].append(n_con_dato)
resumen["metrica"].append("clusters_validados_opticamente")
resumen["valor"].append(n_validados)
resumen["metrica"].append("tasa_validacion_pct")
resumen["valor"].append(round(tasa_validacion, 2))

for cat in ["concordancia_fuerte", "concordancia_moderada", "concordancia_debil",
            "sin_concordancia", "baja_calidad_optica", "sin_dato_optico"]:
    resumen["metrica"].append(f"n_{cat}")
    resumen["valor"].append(conteo.get(cat, 0))

resumen["metrica"].append("umbral_dNDVI")
resumen["valor"].append(UMBRAL_dNDVI)
resumen["metrica"].append("umbral_dNBR")
resumen["valor"].append(UMBRAL_dNBR)
resumen["metrica"].append("umbral_dBSI")
resumen["valor"].append(UMBRAL_dBSI)
resumen["metrica"].append("min_criterios")
resumen["valor"].append(MIN_CRITERIOS)
resumen["metrica"].append("min_obs_calidad")
resumen["valor"].append(MIN_OBS_CALIDAD)

df_resumen = pd.DataFrame(resumen)
resumen_path = os.path.join(DIR_OUTPUT, "P07b_resumen_validacion.csv")
df_resumen.to_csv(resumen_path, index=False)
log(f"  Exportado: {resumen_path}")

# =============================================================================
# SECCIÓN 7: VISUALIZACIONES
# =============================================================================

log("\n--- SECCIÓN 7: VISUALIZACIONES ---")

# 7.1 Mapa de concordancia espacial
log("\n[7.1] Generando mapa de concordancia...")

fig, ax = plt.subplots(1, 1, figsize=(12, 10))

# Fondo: NDVI pre-evento
with rasterio.open(RASTERS_OPTICOS["NDVI_pre"]) as src:
    ndvi_pre = src.read(1)
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
    ndvi_pre = np.where(ndvi_pre == src.nodata, np.nan, ndvi_pre) if src.nodata else ndvi_pre

ax.imshow(ndvi_pre, extent=extent, cmap="Greens", alpha=0.5,
          vmin=0, vmax=0.8, origin="upper")

# Plotear centroides coloreados por concordancia
colores_conc = {
    "concordancia_fuerte": "#d62728",       # rojo
    "concordancia_moderada": "#ff7f0e",     # naranja
    "concordancia_debil": "#ffbb78",        # naranja claro
    "sin_concordancia": "#1f77b4",          # azul
    "baja_calidad_optica": "#7f7f7f",       # gris
    "sin_dato_optico": "#d3d3d3",           # gris claro
}

for cat, color in colores_conc.items():
    subset = gdf_base[gdf_base["concordancia"] == cat]
    if len(subset) > 0:
        if subset.geometry.geom_type.iloc[0] == 'Point':
            xs = [g.x for g in subset.geometry]
            ys = [g.y for g in subset.geometry]
        else:
            xs = [g.centroid.x for g in subset.geometry]
            ys = [g.centroid.y for g in subset.geometry]
        n_cat = len(subset)
        label = f"{cat.replace('_', ' ').title()} (n={n_cat})"
        ax.scatter(xs, ys, c=color, s=8, alpha=0.7, label=label, edgecolors='none')

ax.set_title(f"P07b — Concordancia InSAR-Óptico\n"
             f"{DATASET_PRINCIPAL} | Tasa validación: {tasa_validacion:.1f}%",
             fontsize=13, fontweight='bold')
ax.set_xlabel("Longitud (m)")
ax.set_ylabel("Latitud (m)")
ax.legend(loc="lower left", fontsize=8, markerscale=2)

plt.tight_layout()
mapa_path = os.path.join(DIR_OUTPUT, "P07b_mapa_concordancia.png")
plt.savefig(mapa_path, dpi=300, bbox_inches='tight')
plt.close()
log(f"  Exportado: {mapa_path}")

# 7.2 Histogramas de índices ópticos en clústeres InSAR
log("\n[7.2] Generando histogramas de distribución...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

indices_plot = [
    ("opt_dNDVI", "dNDVI", UMBRAL_dNDVI, "Pérdida vegetación", "tab:green"),
    ("opt_dNBR", "dNBR", UMBRAL_dNBR, "Degradación", "tab:red"),
    ("opt_dBSI", "dBSI", UMBRAL_dBSI, "Suelo desnudo", "tab:brown"),
]

for ax, (col, titulo, umbral, desc, color) in zip(axes, indices_plot):
    vals = gdf_base[col].dropna()
    if len(vals) == 0:
        ax.text(0.5, 0.5, "Sin datos", ha='center', va='center', transform=ax.transAxes)
        continue

    ax.hist(vals, bins=50, color=color, alpha=0.7, edgecolor='black', linewidth=0.3)
    ax.axvline(umbral, color='red', linestyle='--', linewidth=2,
               label=f"Umbral: {umbral}")
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)

    # Anotar porcentaje que cumple el criterio
    if "dNDVI" in col:
        pct_cumple = (vals < umbral).sum() / len(vals) * 100
        ax.fill_betweenx([0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 100],
                         vals.min(), umbral, alpha=0.1, color='red')
    else:
        pct_cumple = (vals > umbral).sum() / len(vals) * 100

    ax.set_title(f"{titulo}\n{desc}: {pct_cumple:.1f}% cumple umbral", fontsize=11)
    ax.set_xlabel(titulo)
    ax.set_ylabel("Frecuencia")
    ax.legend(fontsize=8)

plt.suptitle(f"Distribución de índices ópticos Sentinel-2 en clústeres InSAR\n"
             f"({DATASET_PRINCIPAL}, n={len(gdf_base)})", fontsize=13, fontweight='bold')
plt.tight_layout()
hist_path = os.path.join(DIR_OUTPUT, "P07b_histogramas.png")
plt.savefig(hist_path, dpi=300, bbox_inches='tight')
plt.close()
log(f"  Exportado: {hist_path}")

# 7.3 Scatter: deformación InSAR vs dNDVI
log("\n[7.3] Generando scatter deformación vs dNDVI...")

# Buscar columna de desplazamiento en los datos
cols_desp = [c for c in gdf_base.columns if any(
    x in c.lower() for x in ["desp", "disp", "mean_cm", "median_cm", "displacement"])]

if len(cols_desp) > 0:
    col_desp = cols_desp[0]
    log(f"  Columna de desplazamiento encontrada: {col_desp}")

    fig, ax = plt.subplots(figsize=(8, 8))

    # Datos válidos
    mask = ~(gdf_base[col_desp].isna() | gdf_base["opt_dNDVI"].isna())
    x_desp = gdf_base.loc[mask, col_desp]
    y_dndvi = gdf_base.loc[mask, "opt_dNDVI"]

    if len(x_desp) > 0:
        # Colorear por concordancia
        conc = gdf_base.loc[mask, "concordancia"]
        for cat, color in colores_conc.items():
            m = conc == cat
            if m.sum() > 0:
                ax.scatter(x_desp[m], y_dndvi[m], c=color, s=15, alpha=0.6,
                           label=cat.replace('_', ' ').title(), edgecolors='none')

        # Líneas de referencia
        ax.axhline(UMBRAL_dNDVI, color='green', linestyle='--', alpha=0.5,
                    label=f"Umbral dNDVI={UMBRAL_dNDVI}")
        ax.axhline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

        # Correlación
        from scipy import stats as sp_stats
        r, p = sp_stats.pearsonr(x_desp, y_dndvi)
        ax.set_title(f"Deformación InSAR vs. Cambio de Vegetación (dNDVI)\n"
                     f"r = {r:.3f}, p = {p:.2e}, n = {len(x_desp)}", fontsize=12)
        ax.set_xlabel(f"Desplazamiento InSAR ({col_desp}) [cm]")
        ax.set_ylabel("dNDVI (post - pre)")
        ax.legend(fontsize=7, loc="upper left")
        log(f"  Correlación Pearson: r={r:.3f}, p={p:.2e}")

    plt.tight_layout()
    scatter_path = os.path.join(DIR_OUTPUT, "P07b_scatter_deformacion_vs_dNDVI.png")
    plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  Exportado: {scatter_path}")
else:
    log("  INFO: No se encontró columna de desplazamiento para scatter plot")
    log(f"  Columnas disponibles: {list(gdf_base.columns)}")

# 7.4 Diagrama de barras: concordancia por categoría de confianza P06b
if catalogo_path and any(c in gdf_base.columns for c in ["confidence_score", "score"]):
    log("\n[7.4] Generando gráfico concordancia vs confianza P06b...")
    col_score = "confidence_score" if "confidence_score" in gdf_base.columns else "score"

    fig, ax = plt.subplots(figsize=(10, 6))

    scores_unicos = sorted(gdf_base[col_score].unique())
    categorias = ["concordancia_fuerte", "concordancia_moderada",
                  "concordancia_debil", "sin_concordancia"]
    colores_cat = ["#d62728", "#ff7f0e", "#ffbb78", "#1f77b4"]

    x = np.arange(len(scores_unicos))
    width = 0.2

    for j, (cat, color) in enumerate(zip(categorias, colores_cat)):
        pcts = []
        for s in scores_unicos:
            subset = gdf_base[(gdf_base[col_score] == s) & gdf_base["crit_calidad"]]
            if len(subset) > 0:
                pcts.append((subset["concordancia"] == cat).sum() / len(subset) * 100)
            else:
                pcts.append(0)
        ax.bar(x + j * width, pcts, width, label=cat.replace('_', ' ').title(),
               color=color)

    ax.set_xlabel("Score de confianza P06b")
    ax.set_ylabel("Porcentaje de clústeres (%)")
    ax.set_title("Concordancia óptica por nivel de confianza InSAR")
    ax.set_xticks(x + 1.5 * width)
    ax.set_xticklabels([f"Score {s}/4" for s in scores_unicos])
    ax.legend()

    plt.tight_layout()
    bar_path = os.path.join(DIR_OUTPUT, "P07b_concordancia_vs_confianza.png")
    plt.savefig(bar_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  Exportado: {bar_path}")

# =============================================================================
# SECCIÓN 8: RESUMEN FINAL Y LOG
# =============================================================================

log("\n" + "=" * 70)
log("RESUMEN FINAL P07b — VALIDACIÓN CRUZADA InSAR vs. ÓPTICO")
log("=" * 70)
log(f"  Dataset InSAR:              {DATASET_PRINCIPAL}")
log(f"  Total clústeres:            {total}")
log(f"  Con dato óptico válido:     {n_con_dato}")
log(f"  Validados ópticamente:      {n_validados}")
log(f"  TASA DE VALIDACIÓN:         {tasa_validacion:.1f}%")
log(f"")
log(f"  Concordancia fuerte (3/3):  {conteo.get('concordancia_fuerte', 0)}")
log(f"  Concordancia moderada (2/3):{conteo.get('concordancia_moderada', 0)}")
log(f"  Concordancia débil (1/3):   {conteo.get('concordancia_debil', 0)}")
log(f"  Sin concordancia (0/3):     {conteo.get('sin_concordancia', 0)}")
log(f"")
log(f"  Directorio de salida: {DIR_OUTPUT}")
log("=" * 70)

# Guardar log completo
with open(LOG_PATH, 'w', encoding='utf-8') as f:
    f.write("\n".join(log_lines))
log(f"\nLog guardado: {LOG_PATH}")

log("\nP07b completado exitosamente.")