# ============================================================================
# P07d — CRUCE ESPACIAL DE INVENTARIOS InSAR Y ÓPTICO (Paso 1.3)
# ============================================================================
# Proyecto: InSAR Deslizamientos — Zona de Intag, Cotacachi, Ecuador
# Protocolo: Maestro InSAR Intag v2.0 — Fase 1, Paso 1.3
# Fecha: 28 de febrero de 2026
# Autor: Víctor Pinto
#
# PROPÓSITO:
#   Cruzar espacialmente los dos inventarios independientes:
#     - InSAR (P06b): 1,481 clusters de deformación del terreno (centroides)
#     - Óptico (P07c): 4,122 candidatos de pérdida de vegetación (polígonos)
#   para generar un inventario unificado con clasificación de confianza
#   basada en evidencia convergente multi-sensor.
#
# LÓGICA DEL CRUCE:
#   Nivel 1 — INTERSECCIÓN DIRECTA:
#     Centroide InSAR cae dentro de polígono óptico.
#   Nivel 2 — PROXIMIDAD:
#     Centroide InSAR dentro de buffer de X metros del polígono óptico
#     más cercano (captura desplazamiento por diferencia de resolución
#     y desfase temporal entre adquisiciones radar/óptica).
#
# CLASIFICACIÓN DE CONFIANZA:
#   CONCORDANTE_DIRECTO: centroide InSAR dentro de polígono óptico
#   CONCORDANTE_PROXIMO: centroide InSAR dentro de buffer del polígono
#   SOLO_INSAR:  cluster InSAR sin correspondencia óptica
#   SOLO_OPTICO: candidato óptico sin correspondencia InSAR
#
# ENTRADA:
#   P06b: clusters_verificados_M_ASC_CO_12d.gpkg  (puntos, 1481)
#   P07c: P07c_inventario_optico.gpkg              (polígonos, 4122)
#
# SALIDA:
#   inventario_unificado_P07d.gpkg  — inventario cruzado completo
#   P07d_cruce_insar.gpkg           — clusters InSAR con resultado del cruce
#   P07d_cruce_optico.gpkg          — polígonos ópticos con resultado del cruce
#   P07d_resumen_cruce.csv          — tabla resumen
#   P07d_estadisticas_cruce.json    — metadatos y estadísticas
#   P07d_mapa_cruce.png             — mapa de clasificación
#   P07d_dashboard_cruce.png        — dashboard analítico
#   P07d_log.txt                    — log de ejecución
#
# EJECUTAR EN: Spyder (Run file F5)
# PREREQUISITO: P06b y P07c ejecutados
# Ref: Protocolo Maestro InSAR Intag v2.0, Fase 1 Paso 1.3
# ============================================================================

import os
import sys
import json
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# SECCIÓN 0: CONFIGURACIÓN
# ============================================================================

# --- RUTAS ---
DIR_INSAR  = r"D:\POSGRADOS\INTAG\data\sentinel1\hyp3_deformation\verificacion"
DIR_OPTICO = r"D:\POSGRADOS\INTAG\data\sentinel2\VALIDACION\INTAG_InSAR_P07\inventario_optico"
DIR_OUTPUT = r"D:\POSGRADOS\INTAG\data\sentinel2\VALIDACION\INTAG_InSAR_P07\validacion_cruzada"

os.makedirs(DIR_OUTPUT, exist_ok=True)

# --- ARCHIVOS DE ENTRADA ---
INSAR_GPKG  = os.path.join(DIR_INSAR,  "clusters_verificados_M_ASC_CO_12d.gpkg")
OPTICO_GPKG = os.path.join(DIR_OPTICO, "P07c_inventario_optico.gpkg")

# --- PARÁMETROS ---
# Buffer de proximidad en metros para Nivel 2
# Justificación: Sentinel-1 píxel ~14m, Sentinel-2 píxel 10m,
# más incertidumbre de geolocalización y desfase temporal.
# Un buffer conservador de 50m cubre ~3-4 píxeles de tolerancia.
BUFFER_PROXIMIDAD_M = 50

# Buffer extendido para análisis de sensibilidad
BUFFERS_SENSIBILIDAD = [25, 50, 75, 100, 150, 200]

# Dataset InSAR de referencia
TARGET_INSAR = "M_ASC_CO_12d"

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
log("P07d — CRUCE ESPACIAL DE INVENTARIOS InSAR Y ÓPTICO")
log(f"Protocolo Maestro v2.0 — Fase 1, Paso 1.3")
log(f"Fecha de ejecución: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)
log("")
log("--- SECCIÓN 1: VERIFICACIÓN DE INPUTS ---")

archivos = {
    "InSAR (P06b)":  INSAR_GPKG,
    "Óptico (P07c)": OPTICO_GPKG,
}

all_ok = True
for nombre, ruta in archivos.items():
    if os.path.exists(ruta):
        log(f"  OK: {nombre} — {os.path.basename(ruta)}")
    else:
        log(f"  FALTA: {nombre} — {ruta}")
        all_ok = False

if not all_ok:
    log("")
    log("ERROR: Faltan archivos de entrada.")
    log("Ejecutar P06b (InSAR) y P07c (óptico) primero.")
    sys.exit(1)

log("")

# ============================================================================
# SECCIÓN 2: CARGA DE DATOS
# ============================================================================

log("--- SECCIÓN 2: CARGA DE DATOS ---")

# 2.1 Inventario InSAR (centroides = puntos)
gdf_insar = gpd.read_file(INSAR_GPKG)
log(f"  InSAR (P06b): {len(gdf_insar)} clusters")
log(f"    CRS: {gdf_insar.crs}")
log(f"    Tipo geometría: {gdf_insar.geom_type.unique()}")
log(f"    Columnas: {list(gdf_insar.columns)}")

# Verificar que son puntos
if not all(gdf_insar.geom_type == "Point"):
    log("  ALERTA: InSAR contiene geometrías que no son puntos.")
    log("    Extrayendo centroides...")
    gdf_insar["geometry"] = gdf_insar.geometry.centroid

# 2.2 Inventario óptico (polígonos)
gdf_optico = gpd.read_file(OPTICO_GPKG)
log(f"  Óptico (P07c): {len(gdf_optico)} candidatos")
log(f"    CRS: {gdf_optico.crs}")
log(f"    Tipo geometría: {gdf_optico.geom_type.unique()}")
log(f"    Columnas: {list(gdf_optico.columns)}")

# 2.3 Verificar CRS compatible
if gdf_insar.crs != gdf_optico.crs:
    log(f"  ALERTA: CRS diferentes — reproyectando óptico a {gdf_insar.crs}")
    gdf_optico = gdf_optico.to_crs(gdf_insar.crs)
    log(f"  Reproyectado a {gdf_optico.crs}")

crs_trabajo = gdf_insar.crs
log(f"  CRS de trabajo: {crs_trabajo}")

# 2.4 Resumen de atributos clave
log("")
log("  Resumen InSAR:")
if "confianza" in gdf_insar.columns:
    for nivel in ["ALTA", "MODERADA", "BAJA", "FALSO_POSITIVO"]:
        n = len(gdf_insar[gdf_insar["confianza"] == nivel])
        if n > 0:
            log(f"    {nivel}: {n}")

log("  Resumen óptico:")
if "evidencia" in gdf_optico.columns:
    for nivel in ["fuerte", "moderada", "debil"]:
        n = len(gdf_optico[gdf_optico["evidencia"] == nivel])
        if n > 0:
            log(f"    {nivel}: {n}")

log("")

# ============================================================================
# SECCIÓN 3: CRUCE ESPACIAL — NIVEL 1 (INTERSECCIÓN DIRECTA)
# ============================================================================

log("--- SECCIÓN 3: CRUCE ESPACIAL — NIVEL 1 (INTERSECCIÓN DIRECTA) ---")
log(f"  Método: centroide InSAR dentro de polígono óptico")

# Spatial join: puntos InSAR dentro de polígonos ópticos
sjoin_direct = gpd.sjoin(
    gdf_insar,
    gdf_optico,
    how="left",
    predicate="within"
)

# Identificar clusters InSAR con correspondencia directa
# Un cluster puede caer dentro de múltiples polígonos (poco probable pero posible)
# Tomar el polígono con mayor score espectral en caso de duplicados
insar_con_optico_directo = sjoin_direct[sjoin_direct["index_right"].notna()]
insar_sin_optico = sjoin_direct[sjoin_direct["index_right"].isna()]

# Eliminar duplicados — conservar el match con mayor score espectral
if "score_espectral" in insar_con_optico_directo.columns:
    insar_con_optico_directo = (
        insar_con_optico_directo
        .sort_values("score_espectral", ascending=False)
        .drop_duplicates(subset="cluster_id", keep="first")
    )

n_directo = len(insar_con_optico_directo)
n_sin = len(gdf_insar) - n_directo

log(f"  Resultado Nivel 1:")
log(f"    InSAR con correspondencia DIRECTA: {n_directo} ({100*n_directo/len(gdf_insar):.1f}%)")
log(f"    InSAR sin correspondencia directa: {n_sin} ({100*n_sin/len(gdf_insar):.1f}%)")

# Set de IDs con match directo
ids_insar_directo = set(insar_con_optico_directo["cluster_id"].values)
# Set de IDs ópticos con match directo
ids_optico_directo = set(insar_con_optico_directo["index_right"].dropna().astype(int).values)

log("")

# ============================================================================
# SECCIÓN 4: CRUCE ESPACIAL — NIVEL 2 (PROXIMIDAD CON BUFFER)
# ============================================================================

log("--- SECCIÓN 4: CRUCE ESPACIAL — NIVEL 2 (PROXIMIDAD) ---")
log(f"  Buffer de proximidad: {BUFFER_PROXIMIDAD_M} m")
log(f"  Solo para clusters InSAR sin match directo")

# Clusters InSAR que no tuvieron match directo
ids_sin_directo = set(gdf_insar["cluster_id"].values) - ids_insar_directo
gdf_insar_pendiente = gdf_insar[gdf_insar["cluster_id"].isin(ids_sin_directo)].copy()

log(f"  Clusters InSAR pendientes: {len(gdf_insar_pendiente)}")

# Crear buffer alrededor de polígonos ópticos
log(f"  Creando buffer de {BUFFER_PROXIMIDAD_M}m alrededor de polígonos ópticos...")
gdf_optico_buffer = gdf_optico.copy()
gdf_optico_buffer["geometry"] = gdf_optico.geometry.buffer(BUFFER_PROXIMIDAD_M)

# Spatial join: puntos pendientes dentro de buffer de polígonos ópticos
sjoin_buffer = gpd.sjoin(
    gdf_insar_pendiente,
    gdf_optico_buffer,
    how="left",
    predicate="within"
)

# Identificar matches por proximidad
insar_con_optico_proximo = sjoin_buffer[sjoin_buffer["index_right"].notna()]

# Eliminar duplicados
if len(insar_con_optico_proximo) > 0 and "score_espectral" in insar_con_optico_proximo.columns:
    insar_con_optico_proximo = (
        insar_con_optico_proximo
        .sort_values("score_espectral", ascending=False)
        .drop_duplicates(subset="cluster_id", keep="first")
    )

n_proximo = len(insar_con_optico_proximo)

log(f"  Resultado Nivel 2:")
log(f"    InSAR con correspondencia PRÓXIMA ({BUFFER_PROXIMIDAD_M}m): {n_proximo}")

ids_insar_proximo = set(insar_con_optico_proximo["cluster_id"].values)
ids_optico_proximo = set(insar_con_optico_proximo["index_right"].dropna().astype(int).values)

log("")

# ============================================================================
# SECCIÓN 5: CLASIFICACIÓN FINAL
# ============================================================================

log("--- SECCIÓN 5: CLASIFICACIÓN FINAL ---")

# 5.1 Clasificar cada cluster InSAR
insar_resultado = []
for idx, row in gdf_insar.iterrows():
    cid = row["cluster_id"]
    resultado = {
        "cluster_id": cid,
        "geometry": row.geometry,
    }

    # Copiar atributos InSAR relevantes
    for col in ["n_pixels", "area_ha", "los_mean_cm", "los_min_cm",
                "los_max_cm", "coh_mean", "intensidad", "score",
                "confianza", "T2_pendiente_deg", "T2_elevacion_m"]:
        if col in gdf_insar.columns:
            resultado[f"insar_{col}"] = row.get(col, None)

    # Asignar clasificación del cruce
    if cid in ids_insar_directo:
        resultado["cruce_tipo"] = "CONCORDANTE_DIRECTO"
        # Obtener info del polígono óptico match
        match_row = insar_con_optico_directo[
            insar_con_optico_directo["cluster_id"] == cid
        ].iloc[0]
        resultado["opt_id_match"] = match_row.get("opt_id", None)
        resultado["opt_dNDVI"] = match_row.get("dNDVI_mean", None)
        resultado["opt_score"] = match_row.get("score_espectral", None)
        resultado["opt_evidencia"] = match_row.get("evidencia", None)
        resultado["opt_area_ha"] = match_row.get("area_ha_right", match_row.get("area_ha", None))
        resultado["distancia_m"] = 0.0
    elif cid in ids_insar_proximo:
        resultado["cruce_tipo"] = "CONCORDANTE_PROXIMO"
        match_row = insar_con_optico_proximo[
            insar_con_optico_proximo["cluster_id"] == cid
        ].iloc[0]
        resultado["opt_id_match"] = match_row.get("opt_id", None)
        resultado["opt_dNDVI"] = match_row.get("dNDVI_mean", None)
        resultado["opt_score"] = match_row.get("score_espectral", None)
        resultado["opt_evidencia"] = match_row.get("evidencia", None)
        resultado["opt_area_ha"] = match_row.get("area_ha_right", match_row.get("area_ha", None))
        # Calcular distancia al polígono óptico más cercano
        opt_idx = int(match_row["index_right"])
        opt_geom = gdf_optico.loc[opt_idx, "geometry"]
        resultado["distancia_m"] = round(row.geometry.distance(opt_geom), 1)
    else:
        resultado["cruce_tipo"] = "SOLO_INSAR"
        resultado["opt_id_match"] = None
        resultado["opt_dNDVI"] = None
        resultado["opt_score"] = None
        resultado["opt_evidencia"] = None
        resultado["opt_area_ha"] = None
        resultado["distancia_m"] = None

    insar_resultado.append(resultado)

gdf_cruce_insar = gpd.GeoDataFrame(insar_resultado, crs=crs_trabajo)

# 5.2 Clasificar cada candidato óptico
ids_optico_con_match = ids_optico_directo | ids_optico_proximo

optico_resultado = []
for idx, row in gdf_optico.iterrows():
    resultado = {
        "opt_id": row["opt_id"],
        "geometry": row.geometry,
    }

    # Copiar atributos ópticos relevantes
    for col in ["n_pixels", "area_ha", "dNDVI_mean", "dNDVI_min",
                "dNBR_mean", "dBSI_mean", "NDVI_pre_mean",
                "slope_mean_deg", "score_espectral", "evidencia",
                "cobertura_pre", "worldcover_class"]:
        if col in gdf_optico.columns:
            resultado[f"opt_{col}"] = row.get(col, None)

    # Asignar clasificación
    if idx in ids_optico_directo:
        resultado["cruce_tipo"] = "CON_INSAR_DIRECTO"
    elif idx in ids_optico_proximo:
        resultado["cruce_tipo"] = "CON_INSAR_PROXIMO"
    else:
        resultado["cruce_tipo"] = "SOLO_OPTICO"

    optico_resultado.append(resultado)

gdf_cruce_optico = gpd.GeoDataFrame(optico_resultado, crs=crs_trabajo)

# 5.3 Conteos finales
n_concordante_directo = len(gdf_cruce_insar[gdf_cruce_insar["cruce_tipo"] == "CONCORDANTE_DIRECTO"])
n_concordante_proximo = len(gdf_cruce_insar[gdf_cruce_insar["cruce_tipo"] == "CONCORDANTE_PROXIMO"])
n_concordante_total = n_concordante_directo + n_concordante_proximo
n_solo_insar = len(gdf_cruce_insar[gdf_cruce_insar["cruce_tipo"] == "SOLO_INSAR"])

n_opt_con_insar = len(gdf_cruce_optico[gdf_cruce_optico["cruce_tipo"].isin(
    ["CON_INSAR_DIRECTO", "CON_INSAR_PROXIMO"])])
n_solo_optico = len(gdf_cruce_optico[gdf_cruce_optico["cruce_tipo"] == "SOLO_OPTICO"])

log(f"  RESULTADO FINAL DEL CRUCE:")
log(f"  {'─'*60}")
log(f"  PERSPECTIVA InSAR ({len(gdf_insar)} clusters):")
log(f"    CONCORDANTE DIRECTO:  {n_concordante_directo:5d} ({100*n_concordante_directo/len(gdf_insar):.1f}%)")
log(f"    CONCORDANTE PRÓXIMO:  {n_concordante_proximo:5d} ({100*n_concordante_proximo/len(gdf_insar):.1f}%)")
log(f"    TOTAL CONCORDANTE:    {n_concordante_total:5d} ({100*n_concordante_total/len(gdf_insar):.1f}%)")
log(f"    SOLO INSAR:           {n_solo_insar:5d} ({100*n_solo_insar/len(gdf_insar):.1f}%)")
log(f"  {'─'*60}")
log(f"  PERSPECTIVA ÓPTICO ({len(gdf_optico)} candidatos):")
log(f"    CON CORRESPONDENCIA InSAR: {n_opt_con_insar:5d} ({100*n_opt_con_insar/len(gdf_optico):.1f}%)")
log(f"    SOLO ÓPTICO:               {n_solo_optico:5d} ({100*n_solo_optico/len(gdf_optico):.1f}%)")
log(f"  {'─'*60}")

log("")

# ============================================================================
# SECCIÓN 6: ANÁLISIS CRUZADO DE ATRIBUTOS
# ============================================================================

log("--- SECCIÓN 6: ANÁLISIS CRUZADO DE ATRIBUTOS ---")

# 6.1 Comparación de atributos entre clases InSAR
for clase in ["CONCORDANTE_DIRECTO", "CONCORDANTE_PROXIMO", "SOLO_INSAR"]:
    subset = gdf_cruce_insar[gdf_cruce_insar["cruce_tipo"] == clase]
    if len(subset) > 0:
        log(f"  {clase} (n={len(subset)}):")
        for col, label in [("insar_los_mean_cm", "LOS media (cm)"),
                           ("insar_area_ha", "Área InSAR (ha)"),
                           ("insar_coh_mean", "Coherencia media"),
                           ("insar_T2_pendiente_deg", "Pendiente (°)")]:
            if col in subset.columns:
                vals = pd.to_numeric(subset[col], errors="coerce").dropna()
                if len(vals) > 0:
                    log(f"    {label}: media={vals.mean():.2f}, "
                        f"mediana={vals.median():.2f}, "
                        f"[{vals.min():.2f}, {vals.max():.2f}]")

        # Distribución de confianza InSAR dentro de esta clase
        if "insar_confianza" in subset.columns:
            log(f"    Confianza InSAR:")
            for nivel in ["ALTA", "MODERADA", "BAJA", "FALSO_POSITIVO"]:
                n = len(subset[subset["insar_confianza"] == nivel])
                if n > 0:
                    log(f"      {nivel}: {n} ({100*n/len(subset):.1f}%)")
        log("")

# 6.2 Comparación de atributos entre clases ópticas
for clase in ["CON_INSAR_DIRECTO", "CON_INSAR_PROXIMO", "SOLO_OPTICO"]:
    subset = gdf_cruce_optico[gdf_cruce_optico["cruce_tipo"] == clase]
    if len(subset) > 0:
        log(f"  {clase} (n={len(subset)}):")
        for col, label in [("opt_dNDVI_mean", "dNDVI medio"),
                           ("opt_area_ha", "Área óptica (ha)"),
                           ("opt_slope_mean_deg", "Pendiente (°)"),
                           ("opt_NDVI_pre_mean", "NDVI pre-evento")]:
            if col in subset.columns:
                vals = pd.to_numeric(subset[col], errors="coerce").dropna()
                if len(vals) > 0:
                    log(f"    {label}: media={vals.mean():.4f}, "
                        f"mediana={vals.median():.4f}")

        # Distribución de evidencia óptica
        if "opt_evidencia" in subset.columns:
            log(f"    Evidencia espectral:")
            for nivel in ["fuerte", "moderada", "debil"]:
                n = len(subset[subset["opt_evidencia"] == nivel])
                if n > 0:
                    log(f"      {nivel}: {n} ({100*n/len(subset):.1f}%)")
        log("")

# ============================================================================
# SECCIÓN 7: ANÁLISIS DE SENSIBILIDAD AL BUFFER
# ============================================================================

log("--- SECCIÓN 7: ANÁLISIS DE SENSIBILIDAD AL BUFFER ---")
log(f"  Buffers evaluados: {BUFFERS_SENSIBILIDAD} metros")

sensibilidad = []
for buf in BUFFERS_SENSIBILIDAD:
    # Buffer directo: sjoin ya hecho para buf=0 (directo)
    # Para cada buffer, crear polígonos bufferizados y hacer sjoin
    gdf_buf = gdf_optico.copy()
    if buf > 0:
        gdf_buf["geometry"] = gdf_optico.geometry.buffer(buf)

    sj = gpd.sjoin(gdf_insar, gdf_buf, how="left", predicate="within")
    n_match = sj["index_right"].notna().sum()
    # Clusters únicos con match
    n_unique = sj[sj["index_right"].notna()]["cluster_id"].nunique()

    sensibilidad.append({
        "buffer_m": buf,
        "insar_con_match": n_unique,
        "pct_insar": round(100 * n_unique / len(gdf_insar), 1),
    })
    log(f"    Buffer {buf:4d}m: {n_unique:5d} clusters InSAR con match "
        f"({100*n_unique/len(gdf_insar):.1f}%)")

# Agregar punto sin buffer (intersección directa pura)
sensibilidad_df = pd.DataFrame(sensibilidad)

log("")

# ============================================================================
# SECCIÓN 8: EXPORTACIÓN
# ============================================================================

log("--- SECCIÓN 8: EXPORTACIÓN ---")

# 8.1 GeoPackage — Clusters InSAR con resultado del cruce
gpkg_insar = os.path.join(DIR_OUTPUT, "P07d_cruce_insar.gpkg")
gdf_cruce_insar.to_file(gpkg_insar, driver="GPKG")
log(f"  Exportado: {gpkg_insar}")

# 8.2 GeoPackage — Polígonos ópticos con resultado del cruce
gpkg_optico = os.path.join(DIR_OUTPUT, "P07d_cruce_optico.gpkg")
gdf_cruce_optico.to_file(gpkg_optico, driver="GPKG")
log(f"  Exportado: {gpkg_optico}")

# 8.3 CSV — Resumen del cruce InSAR
csv_insar = os.path.join(DIR_OUTPUT, "P07d_cruce_insar.csv")
gdf_cruce_insar.drop(columns="geometry").to_csv(csv_insar, index=False)
log(f"  Exportado: {csv_insar}")

# 8.4 CSV — Resumen del cruce óptico
csv_optico = os.path.join(DIR_OUTPUT, "P07d_cruce_optico.csv")
gdf_cruce_optico.drop(columns="geometry").to_csv(csv_optico, index=False)
log(f"  Exportado: {csv_optico}")

# 8.5 CSV — Sensibilidad al buffer
csv_sens = os.path.join(DIR_OUTPUT, "P07d_sensibilidad_buffer.csv")
sensibilidad_df.to_csv(csv_sens, index=False)
log(f"  Exportado: {csv_sens}")

# 8.6 JSON — Estadísticas completas
stats = {
    "generado": datetime.now().isoformat(),
    "version": "P07d v1.0",
    "paso_protocolo": "Fase 1, Paso 1.3 — Cruce espacial",
    "parametros": {
        "buffer_proximidad_m": BUFFER_PROXIMIDAD_M,
        "target_insar": TARGET_INSAR,
        "crs": str(crs_trabajo),
    },
    "inventarios_entrada": {
        "insar": {
            "archivo": os.path.basename(INSAR_GPKG),
            "n_clusters": len(gdf_insar),
            "geometria": "puntos (centroides)",
        },
        "optico": {
            "archivo": os.path.basename(OPTICO_GPKG),
            "n_candidatos": len(gdf_optico),
            "geometria": "polígonos",
        },
    },
    "resultados_cruce": {
        "perspectiva_insar": {
            "concordante_directo": n_concordante_directo,
            "concordante_proximo": n_concordante_proximo,
            "concordante_total": n_concordante_total,
            "solo_insar": n_solo_insar,
            "pct_concordante": round(100 * n_concordante_total / len(gdf_insar), 1),
        },
        "perspectiva_optico": {
            "con_insar": n_opt_con_insar,
            "solo_optico": n_solo_optico,
            "pct_con_insar": round(100 * n_opt_con_insar / len(gdf_optico), 1),
        },
    },
    "sensibilidad_buffer": sensibilidad,
}

json_path = os.path.join(DIR_OUTPUT, "P07d_estadisticas_cruce.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(stats, f, indent=2, ensure_ascii=False, default=str)
log(f"  Exportado: {json_path}")

log("")

# ============================================================================
# SECCIÓN 9: VISUALIZACIONES
# ============================================================================

log("--- SECCIÓN 9: VISUALIZACIONES ---")

# 9.1 Mapa de clasificación del cruce
log(f"  Generando mapa de clasificación...")

fig, ax = plt.subplots(1, 1, figsize=(14, 11))

# Capa base: polígonos ópticos SOLO_OPTICO (gris claro)
solo_opt = gdf_cruce_optico[gdf_cruce_optico["cruce_tipo"] == "SOLO_OPTICO"]
if len(solo_opt) > 0:
    solo_opt.plot(ax=ax, color="#D5D8DC", alpha=0.4, edgecolor="#ABB2B9",
                  linewidth=0.2, label=f"Solo óptico ({len(solo_opt)})")

# Polígonos ópticos con InSAR (azul claro)
con_insar = gdf_cruce_optico[gdf_cruce_optico["cruce_tipo"].isin(
    ["CON_INSAR_DIRECTO", "CON_INSAR_PROXIMO"])]
if len(con_insar) > 0:
    con_insar.plot(ax=ax, color="#85C1E9", alpha=0.5, edgecolor="#2E86C1",
                   linewidth=0.3, label=f"Óptico + InSAR ({len(con_insar)})")

# Puntos InSAR por clase
colores_insar = {
    "CONCORDANTE_DIRECTO": "#27AE60",
    "CONCORDANTE_PROXIMO": "#F39C12",
    "SOLO_INSAR":          "#E74C3C",
}

for clase, color in colores_insar.items():
    subset = gdf_cruce_insar[gdf_cruce_insar["cruce_tipo"] == clase]
    if len(subset) > 0:
        subset.plot(ax=ax, color=color, markersize=8, alpha=0.8,
                    edgecolors="black", linewidth=0.3)

# Leyenda manual
legend_elements = [
    Patch(facecolor="#D5D8DC", edgecolor="#ABB2B9", alpha=0.4,
          label=f"Solo óptico ({len(solo_opt)})"),
    Patch(facecolor="#85C1E9", edgecolor="#2E86C1", alpha=0.5,
          label=f"Óptico + InSAR ({len(con_insar)})"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#27AE60",
               markersize=8, label=f"Concordante directo ({n_concordante_directo})"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#F39C12",
               markersize=8, label=f"Concordante próximo ({n_concordante_proximo})"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#E74C3C",
               markersize=8, label=f"Solo InSAR ({n_solo_insar})"),
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=9,
          title="Clasificación del cruce", title_fontsize=10)

ax.set_title(f"P07d — Cruce Espacial de Inventarios InSAR ({len(gdf_insar)}) vs Óptico ({len(gdf_optico)})\n"
             f"Buffer de proximidad: {BUFFER_PROXIMIDAD_M} m | "
             f"Concordancia InSAR: {100*n_concordante_total/len(gdf_insar):.1f}%",
             fontsize=12, fontweight="bold")
ax.set_xlabel("Este (m) — EPSG:32617")
ax.set_ylabel("Norte (m)")
ax.ticklabel_format(style="plain")
plt.tight_layout()

map_path = os.path.join(DIR_OUTPUT, "P07d_mapa_cruce.png")
fig.savefig(map_path, dpi=200, bbox_inches="tight")
plt.close(fig)
log(f"  Exportado: {map_path}")

# 9.2 Dashboard analítico (4 paneles)
log(f"  Generando dashboard analítico...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Panel 1: Barras — Resultado del cruce InSAR
ax1 = axes[0, 0]
cats = ["Concordante\nDirecto", "Concordante\nPróximo", "Solo\nInSAR"]
vals = [n_concordante_directo, n_concordante_proximo, n_solo_insar]
colors_bar = ["#27AE60", "#F39C12", "#E74C3C"]
bars = ax1.bar(cats, vals, color=colors_bar, edgecolor="black", linewidth=0.5)
for bar, val in zip(bars, vals):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 15,
             f"{val}\n({100*val/len(gdf_insar):.1f}%)",
             ha="center", fontsize=9, fontweight="bold")
ax1.set_ylabel("Clusters InSAR")
ax1.set_title("Resultado del cruce — Perspectiva InSAR", fontweight="bold")

# Panel 2: Sensibilidad al buffer
ax2 = axes[0, 1]
ax2.plot(sensibilidad_df["buffer_m"], sensibilidad_df["pct_insar"],
         "o-", color="#2E86C1", linewidth=2, markersize=6)
ax2.axhline(100*n_concordante_directo/len(gdf_insar), color="#27AE60",
            linestyle="--", alpha=0.7, label=f"Sin buffer ({n_concordante_directo})")
ax2.axvline(BUFFER_PROXIMIDAD_M, color="#E74C3C", linestyle=":",
            alpha=0.7, label=f"Buffer seleccionado ({BUFFER_PROXIMIDAD_M}m)")
ax2.set_xlabel("Buffer (metros)")
ax2.set_ylabel("% clusters InSAR con match óptico")
ax2.set_title("Sensibilidad al buffer de proximidad", fontweight="bold")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# Panel 3: Distribución LOS por clase
ax3 = axes[1, 0]
for clase, color in [("CONCORDANTE_DIRECTO", "#27AE60"),
                      ("CONCORDANTE_PROXIMO", "#F39C12"),
                      ("SOLO_INSAR", "#E74C3C")]:
    subset = gdf_cruce_insar[gdf_cruce_insar["cruce_tipo"] == clase]
    if len(subset) > 0 and "insar_los_mean_cm" in subset.columns:
        vals = pd.to_numeric(subset["insar_los_mean_cm"], errors="coerce").dropna()
        if len(vals) > 0:
            ax3.hist(vals, bins=30, alpha=0.6, color=color,
                     edgecolor="black", linewidth=0.3,
                     label=f"{clase.split('_')[-1]} (n={len(subset)})")
ax3.set_xlabel("Desplazamiento LOS medio (cm)")
ax3.set_ylabel("Frecuencia")
ax3.set_title("Distribución de desplazamiento por clase", fontweight="bold")
ax3.legend(fontsize=8)

# Panel 4: Comparación cruce vs confianza InSAR
ax4 = axes[1, 1]
if "insar_confianza" in gdf_cruce_insar.columns:
    # Tabla cruzada: confianza vs tipo de cruce
    cross = pd.crosstab(
        gdf_cruce_insar["insar_confianza"],
        gdf_cruce_insar["cruce_tipo"]
    )
    # Reordenar
    orden_conf = ["ALTA", "MODERADA", "BAJA", "FALSO_POSITIVO"]
    orden_cruce = ["CONCORDANTE_DIRECTO", "CONCORDANTE_PROXIMO", "SOLO_INSAR"]
    cross = cross.reindex(index=[c for c in orden_conf if c in cross.index],
                          columns=[c for c in orden_cruce if c in cross.columns],
                          fill_value=0)
    cross.plot(kind="bar", ax=ax4,
               color=[colores_insar.get(c, "gray") for c in cross.columns],
               edgecolor="black", linewidth=0.3)
    ax4.set_xlabel("Confianza InSAR (P06b)")
    ax4.set_ylabel("Clusters")
    ax4.set_title("Confianza InSAR vs resultado del cruce", fontweight="bold")
    ax4.legend(fontsize=7, title="Cruce", title_fontsize=8)
    ax4.tick_params(axis="x", rotation=0)

fig.suptitle(f"P07d — Dashboard del Cruce Espacial\n"
             f"InSAR: {len(gdf_insar)} clusters | Óptico: {len(gdf_optico)} candidatos | "
             f"Buffer: {BUFFER_PROXIMIDAD_M}m",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()

dash_path = os.path.join(DIR_OUTPUT, "P07d_dashboard_cruce.png")
fig.savefig(dash_path, dpi=200, bbox_inches="tight")
plt.close(fig)
log(f"  Exportado: {dash_path}")

log("")

# ============================================================================
# SECCIÓN 10: RESUMEN FINAL
# ============================================================================

log("=" * 70)
log("RESUMEN FINAL P07d — CRUCE ESPACIAL DE INVENTARIOS")
log("=" * 70)
log(f"")
log(f"  ENTRADAS:")
log(f"    InSAR (P06b):  {len(gdf_insar)} clusters de deformación (centroides)")
log(f"    Óptico (P07c): {len(gdf_optico)} candidatos de pérdida de vegetación (polígonos)")
log(f"    Buffer:        {BUFFER_PROXIMIDAD_M} m")
log(f"")
log(f"  RESULTADO DEL CRUCE — PERSPECTIVA InSAR:")
log(f"    CONCORDANTE DIRECTO:  {n_concordante_directo:5d} ({100*n_concordante_directo/len(gdf_insar):.1f}%)")
log(f"    CONCORDANTE PRÓXIMO:  {n_concordante_proximo:5d} ({100*n_concordante_proximo/len(gdf_insar):.1f}%)")
log(f"    ─────────────────────────────")
log(f"    TOTAL CONCORDANTE:    {n_concordante_total:5d} ({100*n_concordante_total/len(gdf_insar):.1f}%)")
log(f"    SOLO INSAR:           {n_solo_insar:5d} ({100*n_solo_insar/len(gdf_insar):.1f}%)")
log(f"")
log(f"  RESULTADO DEL CRUCE — PERSPECTIVA ÓPTICA:")
log(f"    CON CORRESPONDENCIA InSAR: {n_opt_con_insar:5d} ({100*n_opt_con_insar/len(gdf_optico):.1f}%)")
log(f"    SOLO ÓPTICO:               {n_solo_optico:5d} ({100*n_solo_optico/len(gdf_optico):.1f}%)")
log(f"")
log(f"  INTERPRETACIÓN:")
log(f"    Los clusters CONCORDANTES representan deslizamientos confirmados")
log(f"    por ambos sensores (deformación + pérdida de vegetación).")
log(f"    SOLO INSAR = probable movimiento subsuperficial bajo dosel intacto.")
log(f"    SOLO ÓPTICO = pérdida de vegetación sin deformación InSAR detectable")
log(f"    (erosión superficial, cárcavas, deforestación, o movimientos < 3 cm).")
log(f"")
log(f"  ARCHIVOS DE SALIDA:")
log(f"    {DIR_OUTPUT}/")
log(f"      P07d_cruce_insar.gpkg           — clusters InSAR clasificados")
log(f"      P07d_cruce_optico.gpkg          — polígonos ópticos clasificados")
log(f"      P07d_cruce_insar.csv            — tabla InSAR")
log(f"      P07d_cruce_optico.csv           — tabla óptica")
log(f"      P07d_sensibilidad_buffer.csv    — sensibilidad al buffer")
log(f"      P07d_estadisticas_cruce.json    — metadatos completos")
log(f"      P07d_mapa_cruce.png             — mapa de clasificación")
log(f"      P07d_dashboard_cruce.png        — dashboard analítico")
log(f"")
log("=" * 70)

# Guardar log
log_path = os.path.join(DIR_OUTPUT, "P07d_log.txt")
with open(log_path, "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))
log(f"Log guardado: {log_path}")
log(f"")
log(f"P07d completado exitosamente.")
log(f"")
log(f"SIGUIENTE PASO: Fase 2 — DInSAR co-evento (ya completada)")
log(f"  o Fase 3 — SBAS time series para clusters SOLO_INSAR")