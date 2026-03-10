# -*- coding: utf-8 -*-
# =============================================================================
# P05 v2.0 — EXTRACCIÓN, RECORTE Y ANÁLISIS InSAR (CORREGIDO)
# =============================================================================
# Proyecto: InSAR Deslizamientos Zona de Intag, Cotacachi, Ecuador
# Autor:    Víctor Pinto
# Fecha:    2026-02-27
# Versión:  2.0 (corrige errores de v1.0 documentados en P05_FIX)
#
# CORRECCIONES RESPECTO A v1.0:
# 1. PARES_MAESTRO: fechas verificadas contra P03 (fuente de verdad).
# En v1.0 las fechas de D1-D4 y A1-A4 estaban intercambiadas.
# 2. Clasificación temporal automática basada en fecha del evento.
# En v1.0 las etiquetas eran manuales e incorrectas.
# 3. Sin mosaicos: la fusión CRB se hace en P05b_v2 por separado,
# garantizando que solo se fusionen pares de la misma categoría.
#
# FLUJO:
# 1. Verificar librerías
# 2. Extraer 12 .zip de hyp3_products/
# 3. Cargar AOI (INTAJ.gpkg)
# 4. Identificar cada carpeta HyP3 → par_id (con fechas corregidas)
# 5. Recortar todos los rasters al AOI
# 6. Estadísticas de coherencia y desplazamiento por par
# 7. Visualizaciones individuales y comparativas
# 8. Exportar inventario + estadísticas + clasificación temporal
#
# PRODUCTOS HyP3 INSAR_GAMMA por par:
# *_unw_phase.tif     → Fase unwrapped (rad)
# *_corr.tif          → Coherencia (0-1)
# *_los_disp.tif      → Desplazamiento LOS (m)  [HyP3 entrega en metros]
# *_vert_disp.tif     → Desplazamiento vertical (m)
# *_wrapped_phase.tif → Fase envuelta (rad)
# *_dem.tif           → DEM Copernicus GLO-30 (m)
# *_inc_map.tif       → Ángulo de incidencia (°)
# *_lv_theta.tif      → Look vector theta
# *_lv_phi.tif        → Look vector phi
# *_amp.tif           → Amplitud
#
# ENTRADA: hyp3_products/*.zip + INTAJ.gpkg
# SALIDA:  hyp3_extracted/  → GeoTIFF originales
# hyp3_clipped/    → recortados al AOI (por par_id/)
# hyp3_analysis/   → figuras + estadísticas JSON
#
# PREREQUISITOS: P05_LIMPIEZA ejecutado (carpetas vacías listas)
# SIGUIENTE:     P05b_v2 (fusión CRB ascendente)
#
# Ref: Protocolo Maestro InSAR Intag v2.0, Fase 2, Paso 2.6
# https://hyp3-docs.asf.alaska.edu/guides/insar_product_guide/
# =============================================================================

import json, os, sys, glob, zipfile
from datetime import datetime, date

print("=" * 70)
print("P05 v2.0 — EXTRACCIÓN, RECORTE Y ANÁLISIS InSAR (CORREGIDO)")
print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# ============================================================
# 1. VERIFICACIÓN DE LIBRERÍAS
# ============================================================
print("\n--- VERIFICANDO LIBRERÍAS ---\n")

libs = {}
for nombre, modulo in [
    ("rasterio",   "rasterio"),
    ("geopandas",  "geopandas"),
    ("numpy",      "numpy"),
    ("matplotlib", "matplotlib"),
    ("fiona",      "fiona"),
]:
    try:
        mod = __import__(modulo)
        version = getattr(mod, "__version__", "?")
        libs[nombre] = True
        print(f"  ✓ {nombre:12s} v{version}")
    except ImportError:
        libs[nombre] = False
        print(f"  ✗ {nombre:12s} NO instalado")

faltantes = [n for n, ok in libs.items()
             if not ok and n in ("rasterio", "geopandas", "numpy", "matplotlib")]
if faltantes:
    print(f"\n  [ERROR] Librerías faltantes: {', '.join(faltantes)}")
    print(f"  Instalar con: conda install -c conda-forge {' '.join(faltantes)}")
    sys.exit(1)

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.mask import mask as rasterio_mask
import matplotlib
matplotlib.use("Agg")  # Backend sin GUI para evitar bloqueos
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

print("  ✓ Todas las librerías disponibles")

# ============================================================
# 2. CONFIGURACIÓN
# ============================================================
SAVE_DIR       = r"D:\POSGRADOS\INTAG\data\sentinel1"
PRODUCTS_DIR   = os.path.join(SAVE_DIR, "hyp3_products")
EXTRACTED_DIR  = os.path.join(SAVE_DIR, "hyp3_extracted")
CLIPPED_DIR    = os.path.join(SAVE_DIR, "hyp3_clipped")
ANALYSIS_DIR   = os.path.join(SAVE_DIR, "hyp3_analysis")
AOI_PATH       = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"

for d in [EXTRACTED_DIR, CLIPPED_DIR, ANALYSIS_DIR]:
    os.makedirs(d, exist_ok=True)

print(f"\n  Productos:  {PRODUCTS_DIR}")
print(f"  AOI:        {AOI_PATH}")
print(f"  Extraídos:  {EXTRACTED_DIR}")
print(f"  Recortados: {CLIPPED_DIR}")
print(f"  Análisis:   {ANALYSIS_DIR}")

# ============================================================
# 3. FECHA DEL EVENTO Y CLASIFICACIÓN TEMPORAL
# ============================================================
EVENTO = date(2023, 12, 19)
EVENTO_STR = "2023-12-19"


def clasificar_par(pre_date_str, post_date_str, evento=EVENTO):
    # Clasifica un par interferométrico según relación temporal con evento.
    #
    # Criterio estándar DInSAR:
    # PRE-EVENTO:  master < evento AND secondary < evento
    # CO-EVENTO:   master < evento AND secondary >= evento
    # POST-EVENTO: master >= evento AND secondary > evento
    #
    # Parámetros:
    # pre_date_str:  fecha master YYYYMMDD
    # post_date_str: fecha secondary YYYYMMDD
    # evento:        datetime.date del evento
    #
    # Retorna: "PRE-EVENTO", "CO-EVENTO" o "POST-EVENTO"
    pre  = date(int(pre_date_str[:4]), int(pre_date_str[4:6]), int(pre_date_str[6:8]))
    post = date(int(post_date_str[:4]), int(post_date_str[4:6]), int(post_date_str[6:8]))

    if pre < evento and post < evento:
        return "PRE-EVENTO"
    elif pre >= evento and post > evento:
        return "POST-EVENTO"
    else:
        return "CO-EVENTO"


# ============================================================
# 4. INVENTARIO MAESTRO DE LOS 12 PARES (CORREGIDO)
# ============================================================
# FUENTE DE VERDAD: P03 (pares_dinsar_optimos.json) para D1-D4 y A1-A4
#                   P04b_v2 para P18_A3, P18_B1, P18_B2, P18_B3
#
# CORRECCIÓN v2.0: En v1.0 las fechas de D1-D4 contenían las fechas
# de A1-A4 y viceversa. Error de transcripción documentado en P05_FIX.
#
# Las fechas DEBEN coincidir con los timestamps en los nombres de
# carpeta HyP3 para que identificar_par() funcione correctamente.

PARES_MAESTRO = {
    # ── DESCENDENTE Path 40, Frame 590 (100% AOI) ──
    # Fuente: P03, escenas verificadas con ASF DAAC
    "D1": {
        "track": "P40", "direction": "DESCENDING",
        "path": 40, "frame": 590,
        "pre_date": "20231208", "post_date": "20231220",
        "baseline": 12,
        "cobertura": "6/6 parroquias (100%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231208T110106_20231208T110136_051562_063984_C1B0",
        "post_scene": "S1A_IW_SLC__1SDV_20231220T110105_20231220T110134_051737_063FA1_0CB8",
    },
    "D2": {
        "track": "P40", "direction": "DESCENDING",
        "path": 40, "frame": 590,
        "pre_date": "20231126", "post_date": "20231220",
        "baseline": 24,
        "cobertura": "6/6 parroquias (100%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231126T110106_20231126T110136_051387_06337A_9395",
        "post_scene": "S1A_IW_SLC__1SDV_20231220T110105_20231220T110134_051737_063FA1_0CB8",
    },
    "D3": {
        "track": "P40", "direction": "DESCENDING",
        "path": 40, "frame": 590,
        "pre_date": "20231208", "post_date": "20240113",
        "baseline": 36,
        "cobertura": "6/6 parroquias (100%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231208T110106_20231208T110136_051562_063984_C1B0",
        "post_scene": "S1A_IW_SLC__1SDV_20240113T110104_20240113T110133_052087_064B9D_5E75",
    },
    "D4": {
        "track": "P40", "direction": "DESCENDING",
        "path": 40, "frame": 590,
        "pre_date": "20231208", "post_date": "20240125",
        "baseline": 48,
        "cobertura": "6/6 parroquias (100%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231208T110106_20231208T110136_051562_063984_C1B0",
        "post_scene": "S1A_IW_SLC__1SDV_20240125T110103_20240125T110133_052262_06517A_A2E4",
    },

    # ── ASCENDENTE Path 120, Frame 1182 (69.2% AOI) ──
    # Fuente: P03
    "A1": {
        "track": "P120", "direction": "ASCENDING",
        "path": 120, "frame": 1182,
        "pre_date": "20231213", "post_date": "20231225",
        "baseline": 12,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231213T232947_20231213T233023_051642_063C64_14B1",
        "post_scene": "S1A_IW_SLC__1SDV_20231225T232946_20231225T233022_051817_06426C_168C",
    },
    "A2": {
        "track": "P120", "direction": "ASCENDING",
        "path": 120, "frame": 1182,
        "pre_date": "20231201", "post_date": "20231225",
        "baseline": 24,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231201T232947_20231201T233023_051467_06363D_8F6A",
        "post_scene": "S1A_IW_SLC__1SDV_20231225T232946_20231225T233022_051817_06426C_168C",
    },
    "A3": {
        "track": "P120", "direction": "ASCENDING",
        "path": 120, "frame": 1182,
        "pre_date": "20231213", "post_date": "20240106",
        "baseline": 24,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231213T232947_20231213T233023_051642_063C64_14B1",
        "post_scene": "S1A_IW_SLC__1SDV_20240106T232945_20240106T233021_051992_06486B_CF9A",
    },
    "A4": {
        "track": "P120", "direction": "ASCENDING",
        "path": 120, "frame": 1182,
        "pre_date": "20231213", "post_date": "20240130",
        "baseline": 48,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_scene": "S1A_IW_SLC__1SDV_20231213T232947_20231213T233023_051642_063C64_14B1",
        "post_scene": "S1A_IW_SLC__1SDV_20240130T232944_20240130T233020_052342_065447_619B",
    },

    # ── ASCENDENTE Path 18, Frame 1180 (69.8% AOI) ──
    # Fuente: P04b_v2 (pares alternativos, evitan granule 2023-12-06)
    "P18_A3": {
        "track": "P18", "direction": "ASCENDING",
        "path": 18, "frame": 1180,
        "pre_date": "20231218", "post_date": "20240111",
        "baseline": 24,
        "cobertura": "~69.8% AOI",
    },
    "P18_B1": {
        "track": "P18", "direction": "ASCENDING",
        "path": 18, "frame": 1180,
        "pre_date": "20231124", "post_date": "20231218",
        "baseline": 24,
        "cobertura": "~69.8% AOI",
    },
    "P18_B2": {
        "track": "P18", "direction": "ASCENDING",
        "path": 18, "frame": 1180,
        "pre_date": "20231218", "post_date": "20231230",
        "baseline": 12,
        "cobertura": "~69.8% AOI",
    },
    "P18_B3": {
        "track": "P18", "direction": "ASCENDING",
        "path": 18, "frame": 1180,
        "pre_date": "20231230", "post_date": "20240123",
        "baseline": 24,
        "cobertura": "~69.8% AOI",
    },
}

# Aplicar clasificación temporal automática a cada par
for par_id, info in PARES_MAESTRO.items():
    info["categoria"] = clasificar_par(info["pre_date"], info["post_date"])

# Mostrar clasificación
print(f"\n--- CLASIFICACIÓN TEMPORAL (evento: {EVENTO_STR}) ---\n")
print(f"  {'Par':<8} {'Path':>5} {'Dir':>6}  {'Master':>12}  {'Secondary':>12}  "
      f"{'Días':>4}  {'Categoría':<14}")
print(f"  {'─' * 80}")
for par_id, info in PARES_MAESTRO.items():
    pre_fmt = f"{info['pre_date'][:4]}-{info['pre_date'][4:6]}-{info['pre_date'][6:]}"
    post_fmt = f"{info['post_date'][:4]}-{info['post_date'][4:6]}-{info['post_date'][6:]}"
    dir_s = "DESC" if info["direction"] == "DESCENDING" else "ASC"
    print(f"  {par_id:<8} {info['path']:>5} {dir_s:>6}  {pre_fmt:>12}  "
          f"{post_fmt:>12}  {info['baseline']:>4}  {info['categoria']:<14}")
print(f"  {'─' * 80}")

# Resumen por categoría
cats = {}
for par_id, info in PARES_MAESTRO.items():
    cats.setdefault(info["categoria"], []).append(par_id)
for cat in ["PRE-EVENTO", "CO-EVENTO", "POST-EVENTO"]:
    pares = cats.get(cat, [])
    print(f"  {cat}: {len(pares)} pares → {', '.join(pares)}")

# Sufijos de productos HyP3
SUFIJOS = {
    "unw_phase":     "_unw_phase.tif",
    "corr":          "_corr.tif",
    "los_disp":      "_los_disp.tif",
    "vert_disp":     "_vert_disp.tif",
    "wrapped_phase": "_wrapped_phase.tif",
    "dem":           "_dem.tif",
    "inc_map":       "_inc_map.tif",
    "amp":           "_amp.tif",
    "lv_theta":      "_lv_theta.tif",
    "lv_phi":        "_lv_phi.tif",
}

# ============================================================
# 5. CARGAR AOI
# ============================================================
print("\n--- CARGANDO AOI ---")

if not os.path.exists(AOI_PATH):
    print(f"  [ERROR] No se encontró: {AOI_PATH}")
    sys.exit(1)

aoi = gpd.read_file(AOI_PATH)
print(f"  ✓ AOI cargado: {AOI_PATH}")
print(f"    CRS: {aoi.crs}")
print(f"    Bounds: {aoi.total_bounds}")
print(f"    Parroquias: {len(aoi)}")

if aoi.crs and aoi.crs.to_epsg() != 4326:
    aoi_4326 = aoi.to_crs(epsg=4326)
    print(f"    Reproyectado a EPSG:4326")
else:
    aoi_4326 = aoi

aoi_geom = aoi_4326.geometry.values

# ============================================================
# 6. EXTRAER ARCHIVOS ZIP
# ============================================================
print(f"\n--- EXTRAYENDO PRODUCTOS ---")

zips = sorted(glob.glob(os.path.join(PRODUCTS_DIR, "*.zip")))
print(f"  Archivos .zip encontrados: {len(zips)}")

if len(zips) == 0:
    print(f"  [ERROR] No hay archivos .zip en {PRODUCTS_DIR}")
    sys.exit(1)

carpetas_extraidas = []

for i, zpath in enumerate(zips, 1):
    zname = os.path.basename(zpath)
    carpeta = os.path.join(EXTRACTED_DIR, os.path.splitext(zname)[0])

    if os.path.exists(carpeta) and len(os.listdir(carpeta)) > 0:
        print(f"  [{i}/{len(zips)}] {zname[:60]}... → ya extraído")
        carpetas_extraidas.append(carpeta)
        continue

    print(f"  [{i}/{len(zips)}] Extrayendo {zname[:60]}...")
    try:
        with zipfile.ZipFile(zpath, "r") as zf:
            zf.extractall(EXTRACTED_DIR)
        carpetas_extraidas.append(carpeta)
        n_archivos = len(os.listdir(carpeta)) if os.path.exists(carpeta) else 0
        print(f"    ✓ {n_archivos} archivos")
    except Exception as e:
        print(f"    ✗ Error: {e}")

print(f"\n  Carpetas extraídas: {len(carpetas_extraidas)}")

# ============================================================
# 7. IDENTIFICAR PARES Y SUS PRODUCTOS
# ============================================================
print(f"\n--- IDENTIFICANDO PRODUCTOS ---")


def identificar_par(nombre_carpeta):
    # Identifica el par_id a partir del nombre de carpeta HyP3.
    #
    # Los nombres HyP3 contienen timestamps YYYYMMDDTHHMMSS de ambas
    # escenas. Buscamos pre_date y post_date del PARES_MAESTRO corregido.
    #
    # Retorna par_id o None si no se identifica.
    for par_id, info in PARES_MAESTRO.items():
        pre = info["pre_date"]
        post = info["post_date"]
        if pre in nombre_carpeta and post in nombre_carpeta:
            return par_id
    return None


inventario = {}

for carpeta in sorted(carpetas_extraidas):
    nombre = os.path.basename(carpeta)
    if not os.path.isdir(carpeta):
        continue

    par_id = identificar_par(nombre)

    # Buscar productos por sufijo
    productos = {}
    for clave, sufijo in SUFIJOS.items():
        encontrados = glob.glob(os.path.join(carpeta, f"*{sufijo}"))
        if encontrados:
            productos[clave] = encontrados[0]

    if par_id:
        inventario[par_id] = {
            "carpeta":   carpeta,
            "nombre":    nombre,
            "productos": productos,
            **{k: v for k, v in PARES_MAESTRO[par_id].items()},
        }
        info = PARES_MAESTRO[par_id]
        print(f"  ✓ {par_id:8s} ({info['track']}) "
              f"{info['pre_date']}→{info['post_date']} "
              f"[{info['categoria']}] → {len(productos)} productos")
    else:
        print(f"  ? {nombre[:60]}... → no identificado, {len(productos)} productos")

print(f"\n  Pares identificados: {len(inventario)}/{len(PARES_MAESTRO)}")

# Verificar 12 de 12
faltantes = [p for p in PARES_MAESTRO if p not in inventario]
if faltantes:
    print(f"  ⚠ Pares faltantes: {faltantes}")
    print(f"    Verificar nombres de carpeta vs fechas en PARES_MAESTRO")
    print(f"    Carpetas disponibles:")
    for c in sorted(carpetas_extraidas):
        print(f"      {os.path.basename(c)}")

# ============================================================
# 8. RECORTAR RASTERS AL AOI
# ============================================================
print(f"\n--- RECORTANDO AL AOI ---\n")


def recortar_raster(src_path, dst_path, aoi_gdf, aoi_geom_4326):
    # Recorta un raster GeoTIFF a las geometrías del AOI.
    #
    # Maneja reproyección automática si el raster no está en EPSG:4326
    # (frecuente: HyP3 entrega en UTM 17N o 18N según la zona).
    #
    # Parámetros:
    # src_path:       ruta del GeoTIFF fuente
    # dst_path:       ruta del GeoTIFF recortado
    # aoi_gdf:        GeoDataFrame del AOI (para reproyectar)
    # aoi_geom_4326:  geometrías en WGS84
    #
    # Retorna: True si exitoso, False si error
    try:
        with rasterio.open(src_path) as src:
            src_crs = src.crs
            if src_crs and str(src_crs) != "EPSG:4326":
                aoi_reproj = aoi_gdf.to_crs(src_crs)
                geoms = aoi_reproj.geometry.values
            else:
                geoms = aoi_geom_4326

            out_image, out_transform = rasterio_mask(
                src, geoms, crop=True, nodata=src.nodata
            )

            out_meta = src.meta.copy()
            out_meta.update({
                "driver":    "GTiff",
                "height":    out_image.shape[1],
                "width":     out_image.shape[2],
                "transform": out_transform,
                "compress":  "lzw",
            })

            os.makedirs(os.path.dirname(dst_path), exist_ok=True)
            with rasterio.open(dst_path, "w", **out_meta) as dst:
                dst.write(out_image)
        return True
    except Exception as e:
        print(f"      ✗ Error recortando: {e}")
        return False


total_recortes = 0
total_errores = 0

for par_id, item in sorted(inventario.items()):
    dir_salida = os.path.join(CLIPPED_DIR, par_id)
    n_prod = len(item["productos"])

    # Verificar si ya están recortados
    existentes = 0
    for clave, src_path in item["productos"].items():
        dst_path = os.path.join(dir_salida, os.path.basename(src_path))
        if os.path.exists(dst_path):
            existentes += 1

    if existentes == n_prod:
        print(f"  {par_id}: {n_prod} productos ya recortados, saltando")
        total_recortes += n_prod
        continue

    print(f"  {par_id}: recortando {n_prod} productos...")

    for clave, src_path in item["productos"].items():
        nombre_src = os.path.basename(src_path)
        dst_path = os.path.join(dir_salida, nombre_src)

        if os.path.exists(dst_path):
            total_recortes += 1
            continue

        ok = recortar_raster(src_path, dst_path, aoi, aoi_geom)
        if ok:
            total_recortes += 1
        else:
            total_errores += 1

print(f"\n  ✓ Recortes completados: {total_recortes}")
if total_errores > 0:
    print(f"  ⚠ Errores: {total_errores}")

# ============================================================
# 9. ESTADÍSTICAS DE COHERENCIA Y DESPLAZAMIENTO
# ============================================================
print(f"\n{'=' * 70}")
print("ANÁLISIS ESTADÍSTICO")
print(f"{'=' * 70}\n")


def calcular_stats(raster_path, nombre, factor=1.0, valid_range=None):
    # Calcula estadísticas descriptivas de un raster.
    #
    # Parámetros:
    # raster_path: ruta del GeoTIFF
    # nombre:      prefijo para las claves del dict
    # factor:      multiplicador (ej. 100 para m→cm)
    # valid_range:  (min, max) para filtrar valores fuera de rango
    #
    # Retorna: dict con mean, median, std, min, max, p5, p95, n_valid
    try:
        with rasterio.open(raster_path) as src:
            data = src.read(1).astype(float)
            nodata = src.nodata
            if nodata is not None:
                data[data == nodata] = np.nan
            if valid_range:
                data[(data < valid_range[0]) | (data > valid_range[1])] = np.nan
            data = data * factor
            valid = data[~np.isnan(data)]
            if len(valid) == 0:
                return None
            return {
                f"{nombre}_mean":    float(np.mean(valid)),
                f"{nombre}_median":  float(np.median(valid)),
                f"{nombre}_std":     float(np.std(valid)),
                f"{nombre}_min":     float(np.min(valid)),
                f"{nombre}_max":     float(np.max(valid)),
                f"{nombre}_p5":      float(np.percentile(valid, 5)),
                f"{nombre}_p95":     float(np.percentile(valid, 95)),
                f"{nombre}_n_valid": int(len(valid)),
            }
    except Exception:
        return None


print("  Estadísticas por par (12 pares)\n")

estadisticas_pares = []

for par_id, item in sorted(inventario.items()):
    dir_clip = os.path.join(CLIPPED_DIR, par_id)

    stats = {
        "par_id":    par_id,
        "track":     item["track"],
        "direction": item["direction"],
        "path":      item["path"],
        "frame":     item["frame"],
        "pre_date":  item["pre_date"],
        "post_date": item["post_date"],
        "baseline":  item["baseline"],
        "categoria": item["categoria"],
        "cobertura": item["cobertura"],
    }

    # Coherencia (0-1, sin factor)
    corr_files = glob.glob(os.path.join(dir_clip, "*_corr.tif"))
    if corr_files:
        s = calcular_stats(corr_files[0], "corr", valid_range=(0, 1))
        if s:
            stats.update(s)
            # Porcentaje de píxeles con coherencia > 0.3
            with rasterio.open(corr_files[0]) as src:
                cdata = src.read(1).astype(float)
                nodata = src.nodata
                if nodata is not None:
                    cdata[cdata == nodata] = np.nan
                valid_c = cdata[~np.isnan(cdata)]
                if len(valid_c) > 0:
                    stats["corr_pct_above_03"] = float(
                        np.sum(valid_c > 0.3) / len(valid_c) * 100
                    )

    # Desplazamiento LOS (HyP3 entrega en metros → multiplicar x100 para cm)
    los_files = glob.glob(os.path.join(dir_clip, "*_los_disp.tif"))
    if los_files:
        s = calcular_stats(los_files[0], "los_cm", factor=100)
        if s:
            stats.update(s)

    # Desplazamiento vertical (m → cm)
    vert_files = glob.glob(os.path.join(dir_clip, "*_vert_disp.tif"))
    if vert_files:
        s = calcular_stats(vert_files[0], "vert_cm", factor=100)
        if s:
            stats.update(s)

    estadisticas_pares.append(stats)

    # Imprimir resumen
    pre_fmt = f"{item['pre_date'][4:6]}/{item['pre_date'][6:]}"
    post_fmt = f"{item['post_date'][4:6]}/{item['post_date'][6:]}"
    print(f"  {par_id:8s} {item['track']:5s} {pre_fmt}→{post_fmt} "
          f"({item['baseline']}d) [{item['categoria']}]")
    if "corr_mean" in stats:
        print(f"           Coh: μ={stats['corr_mean']:.3f} "
              f"med={stats['corr_median']:.3f} "
              f">0.3={stats.get('corr_pct_above_03', 0):.0f}%")
    if "los_cm_mean" in stats:
        print(f"           LOS: μ={stats['los_cm_mean']:.2f}cm "
              f"[{stats['los_cm_p5']:.2f}, {stats['los_cm_p95']:.2f}]cm")
    if "vert_cm_mean" in stats:
        print(f"           Vert: μ={stats['vert_cm_mean']:.2f}cm "
              f"[{stats['vert_cm_p5']:.2f}, {stats['vert_cm_p95']:.2f}]cm")

# ============================================================
# 10. VISUALIZACIONES
# ============================================================
print(f"\n{'=' * 70}")
print("VISUALIZACIONES")
print(f"{'=' * 70}\n")


def crear_figura_par(par_id, dir_clip, dir_salida, title_extra=""):
    # Genera figura de 4 paneles: coherencia, fase envuelta, LOS, vertical.
    #
    # Parámetros:
    # par_id:      identificador del par
    # dir_clip:    carpeta con rasters recortados
    # dir_salida:  carpeta donde guardar la figura
    # title_extra: texto adicional para el título
    #
    # Retorna: ruta de la figura guardada
    corr_f = glob.glob(os.path.join(dir_clip, "*_corr.tif"))
    wrap_f = glob.glob(os.path.join(dir_clip, "*_wrapped_phase.tif"))
    los_f  = glob.glob(os.path.join(dir_clip, "*_los_disp.tif"))
    vert_f = glob.glob(os.path.join(dir_clip, "*_vert_disp.tif"))

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.suptitle(f"InSAR {par_id} {title_extra}",
                 fontsize=14, fontweight="bold")

    configs = [
        (axes[0, 0], corr_f,  "Coherencia",               "RdYlGn",  0, 1, 1),
        (axes[0, 1], wrap_f,  "Fase envuelta (rad)",       "hsv",     -np.pi, np.pi, 1),
        (axes[1, 0], los_f,   "Desplazamiento LOS (cm)",   "RdBu_r",  None, None, 100),
        (axes[1, 1], vert_f,  "Desplazamiento Vert (cm)",  "RdBu_r",  None, None, 100),
    ]

    for ax, archivos, titulo, cmap, vmin, vmax, factor in configs:
        if archivos:
            with rasterio.open(archivos[0]) as src:
                data = src.read(1).astype(float)
                nodata = src.nodata
                if nodata is not None:
                    data[data == nodata] = np.nan
                data = data * factor

                if vmin is None:
                    valid = data[~np.isnan(data)]
                    if len(valid) > 0:
                        lim = np.percentile(np.abs(valid), 98)
                        vmin, vmax = -lim, lim
                    else:
                        vmin, vmax = -1, 1

                extent = [src.bounds.left, src.bounds.right,
                          src.bounds.bottom, src.bounds.top]

                im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax,
                               extent=extent, interpolation="nearest")
                plt.colorbar(im, ax=ax, shrink=0.8)
        else:
            ax.text(0.5, 0.5, "No disponible", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")

        ax.set_title(titulo, fontsize=11)
        ax.set_xlabel("Longitud")
        ax.set_ylabel("Latitud")

    plt.tight_layout()
    fig_path = os.path.join(dir_salida, f"insar_{par_id}.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig_path


# 10a. Figuras de pares individuales
print("  A) Pares individuales...")
figuras = []
for par_id, item in sorted(inventario.items()):
    dir_clip = os.path.join(CLIPPED_DIR, par_id)
    if not os.path.isdir(dir_clip):
        continue
    try:
        title_extra = (f"— {item['direction']} {item['track']} | "
                       f"{item['pre_date']}→{item['post_date']} "
                       f"({item['baseline']}d) | {item['categoria']}")
        fig_path = crear_figura_par(par_id, dir_clip, ANALYSIS_DIR, title_extra)
        figuras.append(fig_path)
        print(f"    ✓ {par_id} [{item['categoria']}]")
    except Exception as e:
        print(f"    ✗ {par_id}: {e}")

# 10b. Comparativa de coherencia
print("\n  B) Comparativa de coherencia...")
try:
    fig, ax = plt.subplots(figsize=(18, 7))

    labels = []
    means = []
    colors = []
    cat_colors = {
        "PRE-EVENTO":  "#FFC107",  # amarillo
        "CO-EVENTO":   "#F44336",  # rojo
        "POST-EVENTO": "#4CAF50",  # verde
    }

    for s in estadisticas_pares:
        if "corr_mean" not in s:
            continue
        pre_fmt = f"{s['pre_date'][4:6]}/{s['pre_date'][6:]}"
        post_fmt = f"{s['post_date'][4:6]}/{s['post_date'][6:]}"
        labels.append(f"{s['par_id']}\n{s['track']}\n{pre_fmt}→{post_fmt}\n"
                      f"({s['baseline']}d)\n{s['categoria'][:3]}")
        means.append(s["corr_mean"])

        # Color por track, borde por categoría
        if s["direction"] == "DESCENDING":
            colors.append("#FF5722")
        elif s["track"] == "P120":
            colors.append("#2196F3")
        else:
            colors.append("#4CAF50")

    x = np.arange(len(labels))
    bars = ax.bar(x, means, color=colors, alpha=0.85,
                  edgecolor="black", linewidth=0.5)

    ax.axhline(y=0.3, color="red", linestyle="--", linewidth=1,
               label="Umbral mínimo (0.3)")
    ax.axhline(y=0.6, color="green", linestyle="--", linewidth=0.8,
               alpha=0.5, label="Alta confianza (0.6)")

    for xi, bar in zip(x, bars):
        ax.text(xi, bar.get_height() + 0.008, f"{bar.get_height():.3f}",
                ha="center", va="bottom", fontsize=7, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Coherencia media")
    ax.set_ylim(0, max(means) * 1.3 if means else 1)
    ax.set_title("Coherencia media por par interferométrico\n"
                 "12 pares individuales — Clasificación temporal corregida v2.0")

    legend_elements = [
        Patch(facecolor="#FF5722", label="Descending P40 (100% AOI)"),
        Patch(facecolor="#2196F3", label="Ascending P120 (69.2% AOI)"),
        Patch(facecolor="#4CAF50", label="Ascending P18 (69.8% AOI)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    comp_path = os.path.join(ANALYSIS_DIR, "comparativa_coherencia_v2.png")
    fig.savefig(comp_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figuras.append(comp_path)
    print(f"    ✓ {os.path.basename(comp_path)}")
except Exception as e:
    print(f"    ✗ Error: {e}")

# 10c. Comparativa de desplazamiento LOS
print("\n  C) Comparativa de desplazamiento LOS...")
try:
    fig, ax = plt.subplots(figsize=(18, 7))

    labels_d = []
    los_means = []
    los_p5 = []
    los_p95 = []
    colors_d = []

    for s in estadisticas_pares:
        if "los_cm_mean" not in s:
            continue
        pre_fmt = f"{s['pre_date'][4:6]}/{s['pre_date'][6:]}"
        post_fmt = f"{s['post_date'][4:6]}/{s['post_date'][6:]}"
        labels_d.append(f"{s['par_id']}\n{s['track']}\n({s['baseline']}d)\n"
                        f"{s['categoria'][:3]}")
        los_means.append(s["los_cm_mean"])
        los_p5.append(s["los_cm_p5"])
        los_p95.append(s["los_cm_p95"])

        if s["direction"] == "DESCENDING":
            colors_d.append("#FF5722")
        elif s["track"] == "P120":
            colors_d.append("#2196F3")
        else:
            colors_d.append("#4CAF50")

    x = np.arange(len(labels_d))
    ax.bar(x, los_means, color=colors_d, alpha=0.85,
           edgecolor="black", linewidth=0.5)

    # Barras de error (P5 → P95)
    for xi, mean, p5, p95 in zip(x, los_means, los_p5, los_p95):
        ax.plot([xi, xi], [p5, p95], color="black", linewidth=1.5)
        ax.plot([xi - 0.1, xi + 0.1], [p5, p5], color="black", linewidth=1)
        ax.plot([xi - 0.1, xi + 0.1], [p95, p95], color="black", linewidth=1)

    ax.axhline(y=0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_d, fontsize=8)
    ax.set_ylabel("Desplazamiento LOS (cm)")
    ax.set_title("Desplazamiento LOS medio por par interferométrico\n"
                 "Barras: percentiles 5-95 | Positivo = acercamiento al satélite | "
                 "HyP3 → m × 100 = cm")

    ax.legend(handles=[
        Patch(facecolor="#FF5722", label="DESC P40 (100% AOI)"),
        Patch(facecolor="#2196F3", label="ASC P120 (69.2% AOI)"),
        Patch(facecolor="#4CAF50", label="ASC P18 (69.8% AOI)"),
    ], loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    disp_path = os.path.join(ANALYSIS_DIR, "comparativa_desplazamiento_LOS_v2.png")
    fig.savefig(disp_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figuras.append(disp_path)
    print(f"    ✓ {os.path.basename(disp_path)}")
except Exception as e:
    print(f"    ✗ Error: {e}")

# 10d. Comparativa de desplazamiento vertical
print("\n  D) Comparativa de desplazamiento vertical...")
try:
    fig, ax = plt.subplots(figsize=(18, 7))

    labels_v = []
    vert_means = []
    vert_p5 = []
    vert_p95 = []
    colors_v = []

    for s in estadisticas_pares:
        if "vert_cm_mean" not in s:
            continue
        labels_v.append(f"{s['par_id']}\n{s['track']}\n({s['baseline']}d)\n"
                        f"{s['categoria'][:3]}")
        vert_means.append(s["vert_cm_mean"])
        vert_p5.append(s["vert_cm_p5"])
        vert_p95.append(s["vert_cm_p95"])

        if s["direction"] == "DESCENDING":
            colors_v.append("#FF5722")
        elif s["track"] == "P120":
            colors_v.append("#2196F3")
        else:
            colors_v.append("#4CAF50")

    x = np.arange(len(labels_v))
    ax.bar(x, vert_means, color=colors_v, alpha=0.85,
           edgecolor="black", linewidth=0.5)

    for xi, mean, p5, p95 in zip(x, vert_means, vert_p5, vert_p95):
        ax.plot([xi, xi], [p5, p95], color="black", linewidth=1.5)
        ax.plot([xi - 0.1, xi + 0.1], [p5, p5], color="black", linewidth=1)
        ax.plot([xi - 0.1, xi + 0.1], [p95, p95], color="black", linewidth=1)

    ax.axhline(y=0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_v, fontsize=8)
    ax.set_ylabel("Desplazamiento Vertical (cm)")
    ax.set_title("Desplazamiento vertical medio por par interferométrico\n"
                 "Positivo = uplift | Negativo = subsidencia | "
                 "HyP3 → m × 100 = cm")

    ax.legend(handles=[
        Patch(facecolor="#FF5722", label="DESC P40"),
        Patch(facecolor="#2196F3", label="ASC P120"),
        Patch(facecolor="#4CAF50", label="ASC P18"),
    ], loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    vert_path = os.path.join(ANALYSIS_DIR, "comparativa_desplazamiento_VERT_v2.png")
    fig.savefig(vert_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figuras.append(vert_path)
    print(f"    ✓ {os.path.basename(vert_path)}")
except Exception as e:
    print(f"    ✗ Error: {e}")

# ============================================================
# 11. GUARDAR ESTADÍSTICAS Y METADATOS
# ============================================================
print(f"\n--- GUARDANDO RESULTADOS ---")

output = {
    "generado":       datetime.now().isoformat(),
    "script":         "P05_v2.0",
    "version":        "2.0 — Clasificación temporal corregida",
    "evento":         EVENTO_STR,
    "aoi":            AOI_PATH,
    "n_pares":        len(inventario),
    "clasificacion":  {cat: pares for cat, pares in cats.items()},
    "correcciones_v2": [
        "Fechas D1-D4 y A1-A4 corregidas (estaban intercambiadas en v1.0)",
        "Clasificación temporal automática (no manual)",
        "Mosaicos removidos de P05 (se hacen en P05b_v2 por separado)",
    ],
    "pares_maestro":  {k: {kk: vv for kk, vv in v.items()
                           if kk not in ("productos", "carpeta", "nombre")}
                       for k, v in inventario.items()},
    "estadisticas_pares": estadisticas_pares,
}

stats_path = os.path.join(ANALYSIS_DIR, "estadisticas_insar_v2.json")
with open(stats_path, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"  ✓ Estadísticas: {stats_path}")

# También guardar el PARES_MAESTRO como JSON independiente para P05b y P06
maestro_path = os.path.join(SAVE_DIR, "pares_maestro_v2.json")
with open(maestro_path, "w", encoding="utf-8") as f:
    json.dump({
        "generado": datetime.now().isoformat(),
        "evento": EVENTO_STR,
        "version": "2.0",
        "pares": {k: {kk: vv for kk, vv in v.items()}
                  for k, v in PARES_MAESTRO.items()},
    }, f, indent=2, ensure_ascii=False)
print(f"  ✓ Pares maestro: {maestro_path}")

# ============================================================
# 12. RESUMEN FINAL
# ============================================================
print(f"\n{'=' * 70}")
print("RESUMEN P05 v2.0")
print(f"{'=' * 70}")

print(f"\n  INVENTARIO:")
print(f"    Productos .zip:        {len(zips)}")
print(f"    Pares identificados:   {len(inventario)}/12")
print(f"    Rasters recortados:    {total_recortes}")
print(f"    Figuras generadas:     {len(figuras)}")

print(f"\n  CLASIFICACIÓN TEMPORAL (evento: {EVENTO_STR}):")
for cat in ["PRE-EVENTO", "CO-EVENTO", "POST-EVENTO"]:
    pares = cats.get(cat, [])
    print(f"    {cat:14s}: {len(pares):2d} pares → {', '.join(pares)}")

print(f"\n  COBERTURA POR TRACK:")
print(f"    Descending P40:   4 pares (D1-D4), 100% AOI")
print(f"    Ascending  P120:  4 pares (A1-A4), 69.2% AOI")
print(f"    Ascending  P18:   4 pares (P18_*), 69.8% AOI")

print(f"\n  ESTADÍSTICAS (pares individuales):")
print(f"  {'─' * 80}")
print(f"  {'PAR':<8} {'TRACK':5} {'BASE':>5}  {'CATEG':14s}  {'COH μ':>6}  "
      f"{'COH>0.3':>7}  {'LOS cm':>8}  {'VERT cm':>8}")
print(f"  {'─' * 80}")
for s in estadisticas_pares:
    coh  = f"{s.get('corr_mean', 0):.3f}" if "corr_mean" in s else "N/A"
    pct  = f"{s.get('corr_pct_above_03', 0):.0f}%" if "corr_pct_above_03" in s else "N/A"
    los  = f"{s.get('los_cm_mean', 0):.2f}" if "los_cm_mean" in s else "N/A"
    vrt  = f"{s.get('vert_cm_mean', 0):.2f}" if "vert_cm_mean" in s else "N/A"
    print(f"  {s['par_id']:<8} {s['track']:5} {s['baseline']:>4}d  "
          f"{s['categoria']:14s}  {coh:>6}  {pct:>7}  {los:>8}  {vrt:>8}")
print(f"  {'─' * 80}")

print(f"\n  ARCHIVOS DE SALIDA:")
print(f"    Recortados:     {CLIPPED_DIR}")
print(f"    Análisis:       {ANALYSIS_DIR}")
print(f"    Pares maestro:  {maestro_path}")
print(f"    Estadísticas:   {stats_path}")

print(f"\n  INTERPRETACIÓN:")
print(f"    Coherencia > 0.3: zona útil para análisis InSAR")
print(f"    Coherencia > 0.6: alta confianza en desplazamiento")
print(f"    Desp. LOS: positivo = acercamiento, negativo = alejamiento")
print(f"    Desp. Vert: positivo = uplift, negativo = subsidencia")
print(f"    IMPORTANTE: HyP3 entrega desplazamiento en metros.")
print(f"    Factor ×100 aplicado en estadísticas para reportar en cm.")

print(f"\n  SIGUIENTE PASO:")
print(f"    Ejecutar P05b_v2 → Fusión CRB ascendente (solo 2 mosaicos válidos)")
print(f"    Mosaicos: M_ASC_CO_12d (A1+P18_B2) y M_ASC_CO_24d (A3+P18_A3)")
print(f"{'=' * 70}")