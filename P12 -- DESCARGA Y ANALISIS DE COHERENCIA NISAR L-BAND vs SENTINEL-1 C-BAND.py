# -*- coding: utf-8 -*-
# =============================================================================
# P12 -- DESCARGA Y ANALISIS DE COHERENCIA NISAR L-BAND vs SENTINEL-1 C-BAND
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 4 NISAR
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.0
#
# PROPOSITO:
#   1. Descargar los productos NISAR GUNW (interferogramas geocodificados
#      desenrollados) disponibles sobre la zona de Intag
#   2. Extraer mapas de coherencia L-band
#   3. Comparar cuantitativamente con coherencia C-band de Sentinel-1
#   4. Generar estadisticas y figuras para publicacion
#
# CONTEXTO:
#   La Fase 3 (SBAS) demostro que la coherencia temporal con banda C
#   es ~0 en bosque tropical humedo de Intag (coherencia espacial media
#   0.18). Esta fase evalua si la banda L (NISAR, lambda=24 cm) mejora
#   la coherencia, como predice la teoria (Zebker & Villasenor, 1992;
#   Wei & Sandwell, 2010).
#
# DATOS NISAR DISPONIBLES (verificado 2026-03-08):
#   - 58 productos sobre Intag (liberacion 27/02/2026)
#   - 10 GUNW (interferogramas con coherencia)
#   - 12 RSLC, 12 GSLC, 12 GCOV, 12 SME2
#   - Orbita descendente, 6 fechas unicas (2025-10-29 a 2026-01-09)
#   - NOTA: Datos pre-calibracion. La coherencia NO depende de
#     calibracion radiometrica absoluta (Touzi et al., 1999).
#
# ESTRUCTURA NISAR GUNW (HDF5):
#   /science/LSAR/GUNW/grids/frequencyA/interferogram/
#     - coherenceMagnitude   (lo que necesitamos)
#     - unwrappedPhase
#     - connectedComponents
#   /science/LSAR/GUNW/grids/frequencyA/pixelOffsets/
#   /science/LSAR/GUNW/metadata/
#
# REFERENCIA:
#   - NISAR Data User Guide: https://nisar-docs.asf.alaska.edu/
#   - Wei & Sandwell (2010): Decorrelation L vs C, IEEE TGRS
#   - Zebker & Villasenor (1992): Decorrelacion, IEEE TGRS 30(5)
#   - Brancato et al. (2025): NISAR GUNW Product Specification
# =============================================================================

import os
import sys
import glob
from datetime import datetime

try:
    import asf_search as asf
    import h5py
    import numpy as np
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    print("  Ejecutar: pip install asf_search h5py numpy")
    sys.exit(1)

# ============================================================
# CONFIGURACION
# ============================================================
# Directorio de trabajo Fase 4
FASE4_DIR = r"D:\POSGRADOS\INTAG\data\nisar"
GUNW_DIR = os.path.join(FASE4_DIR, "gunw")
RESULTS_DIR = os.path.join(FASE4_DIR, "resultados")
os.makedirs(GUNW_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# AOI de Intag
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718
AOI_WKT = (f"POLYGON(({AOI_WEST} {AOI_SOUTH}, {AOI_EAST} {AOI_SOUTH}, "
           f"{AOI_EAST} {AOI_NORTH}, {AOI_WEST} {AOI_NORTH}, "
           f"{AOI_WEST} {AOI_SOUTH}))")

# Coherencia C-band de referencia (Sentinel-1, Fase 3)
SENTINEL1_COH_FILE = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3\mintpy\avgSpatialCoh.h5"

# Credenciales NASA Earthdata (las mismas de HyP3)
# asf_search usa las credenciales de .netrc o las pide interactivamente

print("=" * 70)
print("P12 -- DESCARGA Y ANALISIS NISAR L-BAND vs SENTINEL-1 C-BAND")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"Directorio: {FASE4_DIR}")
print(f"AOI: W={AOI_WEST}, E={AOI_EAST}, S={AOI_SOUTH}, N={AOI_NORTH}")

# ============================================================
# PASO 1: BUSCAR PRODUCTOS GUNW SOBRE INTAG
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Busqueda de productos NISAR GUNW")
print(f"{'-' * 70}")

results = asf.search(
    platform=asf.PLATFORM.NISAR,
    intersectsWith=AOI_WKT,
    processingLevel="GUNW",
)
print(f"  Productos GUNW encontrados: {len(results)}")

if len(results) == 0:
    print("  ERROR: No se encontraron productos GUNW.")
    print("  Verificar en https://search.asf.alaska.edu con dataset NISAR")
    sys.exit(1)

# Listar productos
print(f"\n  {'Nombre':<65} {'Fecha':>12}")
print(f"  {'-'*65} {'-'*12}")
for r in results:
    nombre = r.properties.get("sceneName", "N/A")[:63]
    fecha = r.properties.get("startTime", "N/A")[:10]
    print(f"  {nombre:<65} {fecha:>12}")

# ============================================================
# PASO 2: DESCARGAR PRODUCTOS GUNW
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Descarga de productos GUNW")
print(f"{'-' * 70}")

# Verificar cuales ya estan descargados
existentes = set(os.listdir(GUNW_DIR))
por_descargar = []
for r in results:
    nombre = r.properties.get("fileName", "")
    if nombre and nombre not in existentes:
        por_descargar.append(r)

print(f"  Total GUNW: {len(results)}")
print(f"  Ya descargados: {len(results) - len(por_descargar)}")
print(f"  Por descargar: {len(por_descargar)}")

if len(por_descargar) > 0:
    print(f"\n  Iniciando descarga de {len(por_descargar)} productos...")
    print(f"  Destino: {GUNW_DIR}")
    print(f"  NOTA: Requiere autenticacion NASA Earthdata.")
    print(f"  Si no tienes .netrc configurado, se pediran credenciales.")

    try:
        # Crear sesion con autenticacion
        session = asf.ASFSession()
        # Intentar autenticacion automatica (usa .netrc si existe)
        # Si no, se usa autenticacion interactiva
        try:
            session.auth_with_creds("vpintopaez@hotmail.com", "")
        except Exception:
            pass

        for i, r in enumerate(por_descargar, 1):
            nombre = r.properties.get("fileName", "producto")
            print(f"  [{i}/{len(por_descargar)}] Descargando {nombre[:60]}...")
            try:
                r.download(path=GUNW_DIR, session=session)
                print(f"    OK")
            except Exception as e:
                print(f"    ERROR: {str(e)[:80]}")
                print(f"    Intentar descarga manual desde ASF Vertex")

    except Exception as e:
        print(f"\n  Error de autenticacion: {e}")
        print(f"  SOLUCION: Descargar manualmente desde ASF Vertex:")
        print(f"    1. Ir a https://search.asf.alaska.edu")
        print(f"    2. Dataset: NISAR")
        print(f"    3. Dibujar AOI sobre Intag")
        print(f"    4. Filtrar por 'GUNW'")
        print(f"    5. Descargar los {len(results)} productos a: {GUNW_DIR}")
        print(f"\n  Cuando los archivos esten descargados, ejecutar este")
        print(f"  script nuevamente para continuar con el analisis.")
        sys.exit(0)

# ============================================================
# PASO 3: EXPLORAR ESTRUCTURA DE LOS GUNW
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Exploracion de estructura GUNW")
print(f"{'-' * 70}")

gunw_files = sorted(glob.glob(os.path.join(GUNW_DIR, "*.h5")))
if len(gunw_files) == 0:
    gunw_files = sorted(glob.glob(os.path.join(GUNW_DIR, "NISAR*GUNW*")))

print(f"  Archivos GUNW en disco: {len(gunw_files)}")

if len(gunw_files) == 0:
    print(f"  No hay archivos GUNW descargados aun.")
    print(f"  Descargar manualmente a: {GUNW_DIR}")
    print(f"  Luego reejecutar este script.")
    sys.exit(0)

# Explorar el primer archivo para entender la estructura
print(f"\n  Explorando estructura del primer GUNW:")
f = h5py.File(gunw_files[0], "r")

def print_hdf5_tree(g, prefix=""):
    for key in g.keys():
        item = g[key]
        if isinstance(item, h5py.Group):
            print(f"  {prefix}{key}/")
            if prefix.count("/") < 4:
                print_hdf5_tree(item, prefix + "  ")
        elif isinstance(item, h5py.Dataset):
            print(f"  {prefix}{key}: {item.shape} {item.dtype}")

print_hdf5_tree(f)
f.close()

# ============================================================
# PASO 4: EXTRAER COHERENCIA L-BAND DE CADA GUNW
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Extraccion de coherencia L-band")
print(f"{'-' * 70}")

# Posibles rutas de coherencia en NISAR GUNW
COH_PATHS = [
    "/science/LSAR/GUNW/grids/frequencyA/interferogram/coherenceMagnitude",
    "/science/LSAR/GUNW/grids/frequencyA/interferogram/coherence",
    "//science/LSAR/GUNW/data/coherenceMagnitude",
]

lband_stats = []
for i, gf in enumerate(gunw_files):
    nombre = os.path.basename(gf)
    f = h5py.File(gf, "r")

    # Buscar la coherencia
    coh_data = None
    coh_path_found = None
    for cp in COH_PATHS:
        try:
            coh_data = f[cp][:]
            coh_path_found = cp
            break
        except (KeyError, Exception):
            continue

    # Si no encontramos en rutas predefinidas, buscar recursivamente
    if coh_data is None:
        def find_coherence(g, path=""):
            for key in g.keys():
                full = f"{path}/{key}"
                item = g[key]
                if isinstance(item, h5py.Dataset):
                    if "coherence" in key.lower():
                        return full
                elif isinstance(item, h5py.Group):
                    result = find_coherence(item, full)
                    if result:
                        return result
            return None

        coh_path_found = find_coherence(f)
        if coh_path_found:
            coh_data = f[coh_path_found][:]

    if coh_data is not None:
        # Filtrar valores validos (0-1)
        valid = coh_data[(coh_data >= 0) & (coh_data <= 1) & (~np.isnan(coh_data))]
        stats = {
            "archivo": nombre[:60],
            "path_coherencia": coh_path_found,
            "shape": coh_data.shape,
            "min": np.min(valid) if len(valid) > 0 else np.nan,
            "max": np.max(valid) if len(valid) > 0 else np.nan,
            "media": np.mean(valid) if len(valid) > 0 else np.nan,
            "mediana": np.median(valid) if len(valid) > 0 else np.nan,
            "pix_validos": len(valid),
            "pix_total": coh_data.size,
        }
        lband_stats.append(stats)
        print(f"  [{i+1}/{len(gunw_files)}] {nombre[:50]}")
        print(f"    Coherencia L-band: media={stats['media']:.4f}, "
              f"mediana={stats['mediana']:.4f}, "
              f"pixeles={stats['pix_validos']}")
    else:
        print(f"  [{i+1}/{len(gunw_files)}] {nombre[:50]}")
        print(f"    ADVERTENCIA: No se encontro dataset de coherencia")
        # Listar datasets disponibles para diagnostico
        def list_datasets(g, path=""):
            datasets = []
            for key in g.keys():
                full = f"{path}/{key}"
                item = g[key]
                if isinstance(item, h5py.Dataset):
                    datasets.append(full)
                elif isinstance(item, h5py.Group):
                    datasets.extend(list_datasets(item, full))
            return datasets
        ds = list_datasets(f)
        print(f"    Datasets disponibles ({len(ds)}):")
        for d in ds[:15]:
            print(f"      {d}")
        if len(ds) > 15:
            print(f"      ... y {len(ds)-15} mas")

    f.close()

# ============================================================
# PASO 5: CARGAR COHERENCIA C-BAND (SENTINEL-1)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Carga de coherencia C-band (Sentinel-1)")
print(f"{'-' * 70}")

if os.path.exists(SENTINEL1_COH_FILE):
    f = h5py.File(SENTINEL1_COH_FILE, "r")
    cband_coh = f["coherence"][:]
    f.close()
    cband_valid = cband_coh[(cband_coh >= 0) & (cband_coh <= 1) & (~np.isnan(cband_coh))]
    print(f"  Archivo: {SENTINEL1_COH_FILE}")
    print(f"  Shape: {cband_coh.shape}")
    print(f"  Coherencia C-band: media={np.mean(cband_valid):.4f}, "
          f"mediana={np.median(cband_valid):.4f}")
else:
    print(f"  ADVERTENCIA: No se encontro {SENTINEL1_COH_FILE}")
    print(f"  La comparacion C vs L se hara solo con estadisticas.")
    cband_coh = None
    cband_valid = None

# ============================================================
# PASO 6: COMPARACION CUANTITATIVA C-BAND vs L-BAND
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Comparacion cuantitativa C-band vs L-band")
print(f"{'-' * 70}")

if len(lband_stats) > 0:
    # Estadisticas agregadas L-band
    medias_l = [s["media"] for s in lband_stats if not np.isnan(s["media"])]
    medianas_l = [s["mediana"] for s in lband_stats if not np.isnan(s["mediana"])]

    lband_media_global = np.mean(medias_l) if medias_l else np.nan
    lband_mediana_global = np.mean(medianas_l) if medianas_l else np.nan

    # C-band
    cband_media = np.mean(cband_valid) if cband_valid is not None else 0.1891
    cband_mediana = np.median(cband_valid) if cband_valid is not None else 0.1601

    print(f"\n  {'METRICA':<30} {'C-BAND (S1)':>15} {'L-BAND (NISAR)':>15} {'MEJORA':>10}")
    print(f"  {'-'*30} {'-'*15} {'-'*15} {'-'*10}")
    print(f"  {'Coherencia media':<30} {cband_media:>15.4f} {lband_media_global:>15.4f} "
          f"{(lband_media_global/cband_media - 1)*100:>9.1f}%")
    print(f"  {'Coherencia mediana':<30} {cband_mediana:>15.4f} {lband_mediana_global:>15.4f} "
          f"{(lband_mediana_global/cband_mediana - 1)*100:>9.1f}%")
    print(f"  {'Numero de interferogramas':<30} {'1,116':>15} {len(lband_stats):>15}")
    print(f"  {'Longitud de onda':<30} {'5.6 cm':>15} {'24 cm':>15}")
    print(f"  {'Sensor':<30} {'Sentinel-1':>15} {'NISAR':>15}")

    # Pixeles con coherencia util
    if cband_valid is not None:
        for umbral in [0.3, 0.5, 0.7]:
            pct_c = np.sum(cband_valid >= umbral) / len(cband_valid) * 100
            # Para L-band, usar estadisticas de los GUNW
            print(f"  {'Pixeles coh >= ' + str(umbral):<30} {pct_c:>14.1f}%")

else:
    print("  No hay estadisticas L-band disponibles.")
    print("  Descargar los GUNW primero.")

# ============================================================
# PASO 7: REPORTE
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Reporte")
print(f"{'-' * 70}")

reporte_path = os.path.join(RESULTS_DIR, "P12_reporte_coherencia_CvsL.txt")
with open(reporte_path, "w", encoding="utf-8") as rpt:
    rpt.write("=" * 70 + "\n")
    rpt.write("P12 -- COMPARACION COHERENCIA C-BAND vs L-BAND\n")
    rpt.write("=" * 70 + "\n")
    rpt.write(f"Fecha: {datetime.now()}\n\n")
    rpt.write("DATOS:\n")
    rpt.write(f"  C-band: Sentinel-1 Path 40 DESC, 1116 pares, 2017-2024\n")
    rpt.write(f"  L-band: NISAR, {len(gunw_files)} GUNW, 2025-2026\n")
    rpt.write(f"  AOI: Intag, Cotacachi, Imbabura, Ecuador\n\n")

    if len(lband_stats) > 0:
        rpt.write("RESULTADOS POR INTERFEROGRAMA NISAR:\n")
        for s in lband_stats:
            rpt.write(f"  {s['archivo']}\n")
            rpt.write(f"    Coherencia: media={s['media']:.4f}, "
                      f"mediana={s['mediana']:.4f}, "
                      f"pixeles={s['pix_validos']}\n")

        rpt.write(f"\nCOMPARACION GLOBAL:\n")
        rpt.write(f"  C-band media: {cband_media:.4f}\n")
        rpt.write(f"  L-band media: {lband_media_global:.4f}\n")
        if lband_media_global > 0 and cband_media > 0:
            rpt.write(f"  Mejora: {(lband_media_global/cband_media - 1)*100:.1f}%\n")

print(f"  Reporte: {reporte_path}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P12 COMPLETADO")
print(f"{'=' * 70}")
print(f"  GUNW descargados: {len(gunw_files)}")
print(f"  Con coherencia extraida: {len(lband_stats)}")
if len(lband_stats) > 0:
    print(f"  Coherencia L-band media: {lband_media_global:.4f}")
    print(f"  Coherencia C-band media: {cband_media:.4f}")
    if lband_media_global > cband_media:
        print(f"  RESULTADO: L-band supera C-band en "
              f"{(lband_media_global/cband_media - 1)*100:.1f}%")
    print(f"\n  SIGUIENTE PASO: Generar figuras comparativas y")
    print(f"  documentar resultados de Fase 4.")
else:
    print(f"\n  Los archivos GUNW necesitan ser descargados.")
    print(f"  Descargar desde ASF Vertex a: {GUNW_DIR}")
    print(f"  Luego reejecutar este script.")
print(f"{'=' * 70}")