# -*- coding: utf-8 -*-
# =============================================================================
# P10 — EXTRACCIÓN Y PREPARACIÓN DE PRODUCTOS SBAS PARA MintPy
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag — Fase 3 SBAS
# Autor    : Víctor Pinto Páez
# Fecha    : 2026-03-03
# Versión  : 1.0
#
# PROPÓSITO:
#   1. Extraer los 1,140 archivos ZIP de HyP3 en subcarpetas individuales
#   2. Recortar (clip) todos los GeoTIFF al AOI de Intag
#   3. Generar el archivo template de MintPy (IntagSenDT40.txt)
#   4. Verificar completitud de la estructura
#
# ENTRADA:
#   - 1,140 archivos ZIP en hyp3_products/
#   - AOI: INTAJ.gpkg (D:\POSGRADOS\INTAG\INTAJ.gpkg) + buffer 0.01°
#
# SALIDA:
#   - Estructura de carpetas MintPy-compatible en mintpy_prep/
#   - IntagSenDT40.txt (template MintPy)
#   - P10_reporte_preparacion.txt
#
# REFERENCIA:
#   MintPy HyP3 directory structure:
#   https://mintpy.readthedocs.io/en/latest/dir_structure/
# =============================================================================

import os
import sys
import glob
import zipfile
import time
import csv
from datetime import datetime

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    print("[ERROR] GDAL no instalado. Ejecutar: pip install GDAL")
    print("  O instalar via conda: conda install -c conda-forge gdal")
    sys.exit(1)

# ============================================================
# CONFIGURACIÓN
# ============================================================
SBAS_DIR     = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"
PRODUCTS_DIR = os.path.join(SBAS_DIR, "hyp3_products")
MINTPY_DIR   = os.path.join(SBAS_DIR, "mintpy_prep")
TEMPLATE_DIR = os.path.join(SBAS_DIR, "mintpy")
os.makedirs(MINTPY_DIR, exist_ok=True)
os.makedirs(TEMPLATE_DIR, exist_ok=True)

# AOI de Intag (EPSG:4326 — lon/lat)
# Extraído de INTAJ.gpkg (6 parroquias) + buffer 0.01° (~1 km)
# Extraídos de INTAJ.gpkg + buffer 0.01° (~1 km)
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718

# Archivos de salida
OUT_TEMPLATE = os.path.join(TEMPLATE_DIR, "IntagSenDT40.txt")
OUT_REPORTE  = os.path.join(SBAS_DIR, "P10_reporte_preparacion.txt")
OUT_CSV      = os.path.join(SBAS_DIR, "P10_inventario_mintpy.csv")

# Sufijos de productos HyP3 que MintPy necesita
SUFFIXES_NEEDED = [
    "unw_phase",
    "corr",
    "dem",
    "lv_theta",
    "lv_phi",
    "water_mask",
]

print("=" * 70)
print("P10 — EXTRACCIÓN Y PREPARACIÓN PARA MintPy")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"ZIP origen: {PRODUCTS_DIR}")
print(f"Destino MintPy: {MINTPY_DIR}")
print(f"AOI: W={AOI_WEST}, E={AOI_EAST}, S={AOI_SOUTH}, N={AOI_NORTH}")

# ============================================================
# PASO 1: Verificar ZIPs disponibles
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Verificación de ZIPs")
print(f"{'-' * 70}")

zips = sorted(glob.glob(os.path.join(PRODUCTS_DIR, "*.zip")))
print(f"  ZIPs encontrados: {len(zips)}")

if len(zips) == 0:
    print(f"  ERROR: No se encontraron ZIPs en {PRODUCTS_DIR}")
    sys.exit(1)

# ============================================================
# PASO 2: Extraer ZIPs
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Extracción de ZIPs")
print(f"{'-' * 70}")

# Verificar cuántos ya están extraídos
carpetas_existentes = set()
for d in os.listdir(MINTPY_DIR):
    if os.path.isdir(os.path.join(MINTPY_DIR, d)):
        carpetas_existentes.add(d)

print(f"  Carpetas ya extraídas: {len(carpetas_existentes)}")

extraidos = 0
errores_ext = 0
errores_ext_detalle = []
t_inicio = time.time()

for i, zpath in enumerate(zips, 1):
    zname = os.path.splitext(os.path.basename(zpath))[0]
    
    # Saltar si ya extraído
    if zname in carpetas_existentes:
        extraidos += 1
        if i % 200 == 0:
            print(f"  [{i}/{len(zips)}] {extraidos} extraídos (saltando existentes)")
        continue
    
    try:
        with zipfile.ZipFile(zpath, 'r') as zf:
            zf.extractall(MINTPY_DIR)
        extraidos += 1
    except Exception as e:
        errores_ext += 1
        errores_ext_detalle.append({"zip": zname, "error": str(e)})
        if errores_ext <= 10:
            print(f"  ERROR extrayendo {zname}: {str(e)[:60]}")
    
    if i % 100 == 0 or i == len(zips):
        t_el = (time.time() - t_inicio) / 60
        print(f"  [{i}/{len(zips)}] {extraidos} OK, {errores_ext} err | {t_el:.1f} min")

print(f"\n  Extracción completada: {extraidos} OK, {errores_ext} errores")

# ============================================================
# PASO 3: Identificar estructura extraída
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Análisis de estructura extraída")
print(f"{'-' * 70}")

# Listar todas las subcarpetas (cada una = 1 interferograma)
subcarpetas = sorted([
    d for d in os.listdir(MINTPY_DIR)
    if os.path.isdir(os.path.join(MINTPY_DIR, d))
])
print(f"  Subcarpetas (interferogramas): {len(subcarpetas)}")

# Analizar contenido de la primera carpeta para entender la estructura
if len(subcarpetas) > 0:
    muestra = subcarpetas[0]
    muestra_path = os.path.join(MINTPY_DIR, muestra)
    archivos_muestra = os.listdir(muestra_path)
    print(f"\n  Muestra de estructura ({muestra}):")
    for f in sorted(archivos_muestra):
        size = os.path.getsize(os.path.join(muestra_path, f)) / (1024**2)
        print(f"    {f} ({size:.1f} MB)")

# ============================================================
# PASO 4: Recortar (clip) GeoTIFFs al AOI
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Recorte de GeoTIFFs al AOI")
print(f"{'-' * 70}")

# Determinar si los productos ya tienen sufijo _clip
# (HyP3 NO recorta — debemos hacerlo nosotros)
# Verificar si ya están recortados (si _clip.tif existe)
muestra_clips = glob.glob(os.path.join(MINTPY_DIR, subcarpetas[0], "*_clip.tif")) if subcarpetas else []
ya_recortados = len(muestra_clips) > 0

if ya_recortados:
    print(f"  Los productos ya parecen estar recortados (_clip.tif encontrados)")
    print(f"  Saltando recorte.")
else:
    print(f"  Recortando {len(subcarpetas)} interferogramas al AOI...")
    print(f"  Esto puede tomar varias horas para 1,140 productos.\n")
    
    # Determinar si los datos están en UTM o WGS84
    # Los productos HyP3 están en UTM
    # Necesitamos usar las coordenadas del AOI en el SRC de los datos
    muestra_tifs = glob.glob(os.path.join(MINTPY_DIR, subcarpetas[0], "*.tif"))
    src_epsg = None
    if muestra_tifs:
        ds = gdal.Open(muestra_tifs[0])
        if ds:
            from osgeo import osr
            srs = osr.SpatialReference(wkt=ds.GetProjection())
            src_epsg = srs.GetAuthorityCode(None)
            print(f"  SRC detectado: EPSG:{src_epsg}")
            ds = None
    
    recortados = 0
    errores_clip = 0
    t_clip_inicio = time.time()
    
    for i, carpeta in enumerate(subcarpetas, 1):
        carpeta_path = os.path.join(MINTPY_DIR, carpeta)
        tifs = glob.glob(os.path.join(carpeta_path, "*.tif"))
        
        for tif in tifs:
            basename = os.path.basename(tif)
            # No recortar si ya es _clip
            if "_clip.tif" in basename:
                continue
            
            # Solo recortar los productos que MintPy necesita
            es_necesario = False
            for suffix in SUFFIXES_NEEDED:
                if suffix in basename:
                    es_necesario = True
                    break
            
            if not es_necesario:
                continue
            
            # Generar nombre de salida con _clip
            out_name = basename.replace(".tif", "_clip.tif")
            out_path = os.path.join(carpeta_path, out_name)
            
            if os.path.exists(out_path):
                continue
            
            try:
                # gdalwarp con reproyección a WGS84 y recorte al AOI
                # -t_srs EPSG:4326 asegura coordenadas geográficas
                # -te define el bounding box de recorte
                gdal.Warp(
                    out_path, tif,
                    dstSRS="EPSG:4326",
                    outputBounds=[AOI_WEST, AOI_SOUTH, AOI_EAST, AOI_NORTH],
                    resampleAlg="nearest",
                    format="GTiff",
                    creationOptions=["COMPRESS=LZW", "TILED=YES"]
                )
                recortados += 1
            except Exception as e:
                errores_clip += 1
                if errores_clip <= 10:
                    print(f"  ERROR clip {basename}: {str(e)[:60]}")
        
        if i % 50 == 0 or i == len(subcarpetas):
            t_el = (time.time() - t_clip_inicio) / 60
            t_rest = (t_el / i) * (len(subcarpetas) - i) if i > 0 else 0
            print(f"  [{i}/{len(subcarpetas)}] "
                  f"{recortados} tifs recortados | "
                  f"{t_el:.1f} min | ~{t_rest:.0f} min rest.")
    
    print(f"\n  Recorte completado: {recortados} archivos, {errores_clip} errores")

# ============================================================
# PASO 5: Verificar completitud para MintPy
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Verificación de completitud MintPy")
print(f"{'-' * 70}")

completos = 0
incompletos = 0
inventario = []

for carpeta in subcarpetas:
    carpeta_path = os.path.join(MINTPY_DIR, carpeta)
    clips = glob.glob(os.path.join(carpeta_path, "*_clip.tif"))
    clip_names = [os.path.basename(c) for c in clips]
    
    # Verificar que cada sufijo necesario tiene su _clip
    sufijos_ok = []
    for suffix in SUFFIXES_NEEDED:
        encontrado = any(suffix in c for c in clip_names)
        sufijos_ok.append(encontrado)
    
    if all(sufijos_ok):
        completos += 1
        estado = "COMPLETO"
    else:
        incompletos += 1
        faltantes = [s for s, ok in zip(SUFFIXES_NEEDED, sufijos_ok) if not ok]
        estado = f"INCOMPLETO: faltan {','.join(faltantes)}"
    
    inventario.append({
        "carpeta": carpeta,
        "n_clips": len(clips),
        "estado": estado,
    })

print(f"  Interferogramas completos: {completos}/{len(subcarpetas)}")
print(f"  Interferogramas incompletos: {incompletos}")

if completos >= 1100:
    print(f"  ✓ Estructura MintPy viable")
else:
    print(f"  ⚠ Revisar productos incompletos")

# ============================================================
# PASO 6: Generar template MintPy
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Generación de template MintPy")
print(f"{'-' * 70}")

# Rutas relativas para MintPy (usar formato Unix con *)
mintpy_data_path = MINTPY_DIR.replace("\\", "/")

template_content = f"""# vim: set filetype=cfg:
# =============================================================================
# MintPy Template — InSAR Deslizamientos Intag
# =============================================================================
# Path 40 Descendente, Frame 590
# Periodo: 2017-01-19 a 2024-06-29
# Escenas: 242 | Pares: 1,140
# Evento: Deslizamiento 2023-12-19
# =============================================================================

########## 1. load_data
mintpy.load.processor      = hyp3

##---------interferogram datasets:
mintpy.load.unwFile        = {mintpy_data_path}/*/*unw_phase_clip.tif
mintpy.load.corFile         = {mintpy_data_path}/*/*corr_clip.tif

##---------geometry datasets:
mintpy.load.demFile         = {mintpy_data_path}/*/*dem_clip.tif
mintpy.load.incAngleFile    = {mintpy_data_path}/*/*lv_theta_clip.tif
mintpy.load.azAngleFile     = {mintpy_data_path}/*/*lv_phi_clip.tif
mintpy.load.waterMaskFile   = {mintpy_data_path}/*/*water_mask_clip.tif

########## 2. modify_network
# Umbral de coherencia temporal para excluir interferogramas ruidosos
mintpy.network.coherenceBased  = yes
mintpy.network.minCoherence    = 0.3

########## 3. reference_point
# Punto de referencia estable (ajustar tras inspección visual)
# Coordenadas aproximadas de zona estable al este del AOI
mintpy.reference.lalo      = auto

########## 4. correct_unwrap_error
mintpy.unwrapError.method  = phase_closure

########## 5. correct_troposphere
# Corrección troposférica con ERA5
# NOTA: Requiere cuenta CDS (Climate Data Store) y .cdsapirc configurado
# Si no está disponible, desactivar con: mintpy.troposphericDelay.method = no
mintpy.troposphericDelay.method = ERA5

########## 6. deramp
# Remover rampa residual (gradiente orbital)
mintpy.deramp               = linear

########## 7. network_inversion
# Coherencia temporal mínima para la inversión
mintpy.networkInversion.minTempCoh = 0.4

########## 8. velocity
# Estimación de velocidad media
mintpy.velocity.startDate   = auto
mintpy.velocity.endDate     = auto
"""

with open(OUT_TEMPLATE, "w", encoding="utf-8") as f:
    f.write(template_content)

print(f"  ✓ Template guardado: {OUT_TEMPLATE}")

# ============================================================
# PASO 7: Inventario y reporte
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Inventario y reporte")
print(f"{'-' * 70}")

with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["carpeta", "n_clips", "estado"])
    for item in inventario:
        writer.writerow([item["carpeta"], item["n_clips"], item["estado"]])
print(f"  ✓ {OUT_CSV}")

# Calcular tamaño total de clips
size_clips = 0
n_clips_total = 0
for carpeta in subcarpetas:
    clips = glob.glob(os.path.join(MINTPY_DIR, carpeta, "*_clip.tif"))
    for c in clips:
        size_clips += os.path.getsize(c)
        n_clips_total += 1

with open(OUT_REPORTE, "w", encoding="utf-8") as f:
    f.write("=" * 70 + "\n")
    f.write("P10 — REPORTE DE PREPARACIÓN PARA MintPy\n")
    f.write("=" * 70 + "\n")
    f.write(f"Fecha: {datetime.now()}\n\n")
    f.write(f"ZIPs procesados: {len(zips)}\n")
    f.write(f"Interferogramas extraídos: {len(subcarpetas)}\n")
    f.write(f"Interferogramas completos (6 _clip.tif): {completos}\n")
    f.write(f"Interferogramas incompletos: {incompletos}\n")
    f.write(f"Total archivos _clip.tif: {n_clips_total}\n")
    f.write(f"Tamaño _clip.tif: {size_clips / (1024**3):.1f} GB\n")
    f.write(f"Template MintPy: {OUT_TEMPLATE}\n")
    f.write(f"\nAOI: W={AOI_WEST}, E={AOI_EAST}, S={AOI_SOUTH}, N={AOI_NORTH}\n")
    f.write(f"SRC salida: EPSG:4326 (WGS84)\n")
print(f"  ✓ {OUT_REPORTE}")

# ============================================================
# RESUMEN FINAL
# ============================================================
print(f"\n{'=' * 70}")
print("P10 COMPLETADO")
print(f"{'=' * 70}")
print(f"  Interferogramas: {len(subcarpetas)}")
print(f"  Completos para MintPy: {completos}")
print(f"  Archivos _clip.tif: {n_clips_total}")
print(f"  Tamaño clips: {size_clips / (1024**3):.1f} GB")
print(f"  Template: {OUT_TEMPLATE}")
print(f"\n  SIGUIENTE PASO:")
print(f"    Instalar MintPy: conda install -c conda-forge mintpy")
print(f"    Ejecutar: smallbaselineApp.py {OUT_TEMPLATE}")
print(f"{'=' * 70}")