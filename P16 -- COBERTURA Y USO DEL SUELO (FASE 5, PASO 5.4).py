# -*- coding: utf-8 -*-
# =============================================================================
# P16 -- COBERTURA Y USO DEL SUELO (FASE 5, PASO 5.4)
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 5 Integracion
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.0
#
# PROPOSITO:
#   Caracterizar la cobertura del suelo en el AOI de Intag y analizar
#   la distribucion de deslizamientos (inventarios InSAR y optico)
#   por tipo de cobertura. Esto establece la relacion entre cobertura
#   vegetal y susceptibilidad a deslizamientos.
#
# FUENTES DE DATOS:
#   1. ESA WorldCover 2021 (10 m): ESA/WorldCover/v200
#      - 11 clases de cobertura global
#      - Referencia: Zanaga et al. (2022), ESA
#
#   2. Sentinel-2 NDVI pre-evento (de P07)
#
# REFERENCIA:
#   - Protocolo Maestro InSAR Intag v2.0, Paso 5.4
#   - Guzzetti et al. (2012): Landslide inventory maps, Earth-Sci Rev
# =============================================================================

import os
import sys
import numpy as np
from datetime import datetime

try:
    import ee
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from osgeo import gdal, osr
    gdal.UseExceptions()
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    sys.exit(1)

ee.Initialize()

# ============================================================
# CONFIGURACION
# ============================================================
FASE5_DIR = r"D:\POSGRADOS\INTAG\data\fase5_integracion"
LC_DIR = os.path.join(FASE5_DIR, "cobertura")
os.makedirs(LC_DIR, exist_ok=True)

AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718
aoi = ee.Geometry.Rectangle([AOI_WEST, AOI_SOUTH, AOI_EAST, AOI_NORTH])

# Clases ESA WorldCover v200
WC_CLASSES = {
    10: {"name": "Bosque cerrado", "color": "#006400"},
    20: {"name": "Arbustos", "color": "#FFBB22"},
    30: {"name": "Pastizal", "color": "#FFFF4C"},
    40: {"name": "Cultivo", "color": "#F096FF"},
    50: {"name": "Urbanizado", "color": "#FA0000"},
    60: {"name": "Suelo desnudo", "color": "#B4B4B4"},
    70: {"name": "Nieve/Hielo", "color": "#F0F0F0"},
    80: {"name": "Agua", "color": "#0064C8"},
    90: {"name": "Humedal herbaceo", "color": "#0096A0"},
    95: {"name": "Manglar", "color": "#00CF75"},
    100: {"name": "Musgo/Liquen", "color": "#FAE6A0"},
}

print("=" * 70)
print("P16 -- COBERTURA Y USO DEL SUELO")
print("=" * 70)
print(f"Fecha: {datetime.now()}")

# ============================================================
# PASO 1: EXTRAER ESA WORLDCOVER 2021
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Extraccion de ESA WorldCover 2021 (10 m)")
print(f"{'-' * 70}")

wc = ee.ImageCollection("ESA/WorldCover/v200").first().clip(aoi)

# Estadisticas por clase
print("  Calculando area por clase...")
class_areas = {}
total_pixels = 0

for code, info in WC_CLASSES.items():
    mask = wc.eq(code)
    area = mask.multiply(ee.Image.pixelArea()).reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=aoi,
        scale=10,
        maxPixels=1e10
    ).get("Map").getInfo()

    if area is not None and area > 0:
        area_km2 = area / 1e6
        class_areas[code] = area_km2
        print(f"    {info['name']:<25} {area_km2:>8.2f} km2")

total_area = sum(class_areas.values())
print(f"\n  Area total AOI: {total_area:.2f} km2")

print(f"\n  Distribucion porcentual:")
print(f"  {'Clase':<25} {'Area (km2)':>10} {'Porcentaje':>10}")
print(f"  {'-'*25} {'-'*10} {'-'*10}")
for code in sorted(class_areas.keys()):
    info = WC_CLASSES[code]
    area = class_areas[code]
    pct = area / total_area * 100
    print(f"  {info['name']:<25} {area:>10.2f} {pct:>9.1f}%")

# ============================================================
# PASO 2: DESCARGAR RASTER DE COBERTURA
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Descarga de raster WorldCover")
print(f"{'-' * 70}")

# Extraer como array via sampleRectangle
# Usar escala de 80m para compatibilidad con DEM/InSAR
wc_resampled = wc.reproject(crs="EPSG:4326", scale=80)

print("  Extrayendo raster (escala 80m para compatibilidad con DEM)...")
try:
    wc_arr_raw = wc_resampled.sampleRectangle(
        region=aoi,
        defaultValue=0
    ).get("Map").getInfo()
    wc_arr = np.array(wc_arr_raw, dtype=np.uint8)
    print(f"  Shape: {wc_arr.shape}")
    print(f"  Clases presentes: {np.unique(wc_arr)}")
except Exception as e:
    print(f"  Error extrayendo raster: {e}")
    print(f"  Usando escala mayor...")
    wc_resampled = wc.reproject(crs="EPSG:4326", scale=200)
    wc_arr_raw = wc_resampled.sampleRectangle(
        region=aoi,
        defaultValue=0
    ).get("Map").getInfo()
    wc_arr = np.array(wc_arr_raw, dtype=np.uint8)
    print(f"  Shape: {wc_arr.shape}")

# Guardar como GeoTIFF
ny, nx = wc_arr.shape
lon_step_wc = (AOI_EAST - AOI_WEST) / nx
lat_step_wc = -(AOI_NORTH - AOI_SOUTH) / ny

tif_path = os.path.join(LC_DIR, "P16_worldcover_2021.tif")
driver = gdal.GetDriverByName("GTiff")
ds = driver.Create(tif_path, nx, ny, 1, gdal.GDT_Byte,
                   options=["COMPRESS=LZW"])
ds.SetGeoTransform([AOI_WEST, lon_step_wc, 0, AOI_NORTH, 0, lat_step_wc])
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)
ds.SetProjection(srs.ExportToWkt())
ds.GetRasterBand(1).WriteArray(wc_arr)
ds.FlushCache()
ds = None
print(f"  GeoTIFF: {tif_path}")

# ============================================================
# PASO 3: NDVI PRE-EVENTO (Sentinel-2)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: NDVI pre-evento (Sentinel-2)")
print(f"{'-' * 70}")

s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED") \
    .filterDate("2023-10-01", "2023-12-18") \
    .filterBounds(aoi) \
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30))

print(f"  Imagenes Sentinel-2 disponibles: {s2.size().getInfo()}")

if s2.size().getInfo() > 0:
    ndvi_pre = s2.median().normalizedDifference(["B8", "B4"]).rename("NDVI").clip(aoi)

    # Estadisticas NDVI por clase de cobertura
    print(f"\n  NDVI medio pre-evento por cobertura:")
    print(f"  {'Clase':<25} {'NDVI medio':>10}")
    print(f"  {'-'*25} {'-'*10}")

    for code in sorted(class_areas.keys()):
        info = WC_CLASSES[code]
        mask = wc.eq(code)
        ndvi_masked = ndvi_pre.updateMask(mask)
        ndvi_mean = ndvi_masked.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=aoi,
            scale=80,
            maxPixels=1e9
        ).get("NDVI").getInfo()

        if ndvi_mean is not None:
            print(f"  {info['name']:<25} {ndvi_mean:>10.4f}")

    # Extraer raster NDVI
    ndvi_resampled = ndvi_pre.reproject(crs="EPSG:4326", scale=80)
    try:
        ndvi_arr_raw = ndvi_resampled.sampleRectangle(
            region=aoi,
            defaultValue=-1
        ).get("NDVI").getInfo()
        ndvi_arr = np.array(ndvi_arr_raw, dtype=np.float32)
        print(f"\n  NDVI raster shape: {ndvi_arr.shape}")
        ndvi_valid = ndvi_arr[ndvi_arr > -0.5]
        print(f"  NDVI: min={np.min(ndvi_valid):.4f}, max={np.max(ndvi_valid):.4f}, "
              f"media={np.mean(ndvi_valid):.4f}")
    except Exception as e:
        print(f"  Error extrayendo NDVI: {e}")
        ndvi_arr = None
else:
    print(f"  Sin imagenes Sentinel-2 disponibles.")
    ndvi_arr = None

# ============================================================
# PASO 4: FIGURAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Generacion de figuras")
print(f"{'-' * 70}")

extent = [AOI_WEST, AOI_EAST, AOI_SOUTH, AOI_NORTH]

# FIGURA 1: Mapa de cobertura + grafico de barras
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle("Cobertura del Suelo - Zona de Intag, Cotacachi\n"
             "ESA WorldCover 2021 (10 m)",
             fontsize=14, fontweight="bold")

# Panel 1: Mapa
# Crear colormap personalizado
from matplotlib.colors import ListedColormap, BoundaryNorm
codes_present = sorted([c for c in class_areas.keys()])
colors_map = [WC_CLASSES[c]["color"] for c in codes_present]
cmap_wc = ListedColormap(colors_map)
bounds = codes_present + [codes_present[-1] + 10]
norm_wc = BoundaryNorm(bounds, cmap_wc.N)

# Mapear valores a indices para colormap
wc_display = np.zeros_like(wc_arr, dtype=np.float32)
for i, code in enumerate(codes_present):
    wc_display[wc_arr == code] = i

cmap_idx = ListedColormap(colors_map)
im = axes[0].imshow(wc_display, cmap=cmap_idx, extent=extent, aspect="auto",
                     vmin=-0.5, vmax=len(codes_present)-0.5,
                     interpolation="nearest")
axes[0].set_xlabel("Longitud")
axes[0].set_ylabel("Latitud")
axes[0].set_title("Mapa de Cobertura")

# Leyenda
patches = [mpatches.Patch(color=WC_CLASSES[c]["color"],
           label=f"{WC_CLASSES[c]['name']} ({class_areas[c]:.0f} km2)")
           for c in codes_present if class_areas[c] > 0.5]
axes[0].legend(handles=patches, loc="lower left", fontsize=7,
               framealpha=0.9)

# Panel 2: Barras
nombres = [WC_CLASSES[c]["name"] for c in codes_present]
areas = [class_areas[c] for c in codes_present]
colores = [WC_CLASSES[c]["color"] for c in codes_present]
pcts = [a / total_area * 100 for a in areas]

y_pos = np.arange(len(nombres))
bars = axes[1].barh(y_pos, pcts, color=colores, edgecolor="gray", linewidth=0.5)
axes[1].set_yticks(y_pos)
axes[1].set_yticklabels(nombres, fontsize=9)
axes[1].set_xlabel("Porcentaje del AOI (%)")
axes[1].set_title("Distribucion de Cobertura")
axes[1].grid(True, alpha=0.3, axis="x")

# Etiquetas en barras
for bar, pct, area in zip(bars, pcts, areas):
    if pct > 2:
        axes[1].text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                     f"{pct:.1f}% ({area:.0f} km2)",
                     va="center", fontsize=8)

plt.tight_layout()
fig1_path = os.path.join(LC_DIR, "P16_fig1_cobertura.png")
plt.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 1: {fig1_path}")

# FIGURA 2: NDVI + Cobertura
if ndvi_arr is not None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle("NDVI Pre-Evento y Cobertura del Suelo - Intag",
                 fontsize=13, fontweight="bold")

    im1 = axes[0].imshow(ndvi_arr, cmap="RdYlGn", extent=extent, aspect="auto",
                          vmin=-0.2, vmax=0.9)
    axes[0].set_xlabel("Longitud")
    axes[0].set_ylabel("Latitud")
    axes[0].set_title("NDVI Pre-Evento (Oct-Dic 2023)")
    plt.colorbar(im1, ax=axes[0], label="NDVI", shrink=0.8)

    im2 = axes[1].imshow(wc_display, cmap=cmap_idx, extent=extent, aspect="auto",
                          vmin=-0.5, vmax=len(codes_present)-0.5,
                          interpolation="nearest")
    axes[1].set_xlabel("Longitud")
    axes[1].set_ylabel("Latitud")
    axes[1].set_title("Cobertura WorldCover 2021")
    axes[1].legend(handles=patches, loc="lower left", fontsize=7,
                   framealpha=0.9)

    plt.tight_layout()
    fig2_path = os.path.join(LC_DIR, "P16_fig2_ndvi_cobertura.png")
    plt.savefig(fig2_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Figura 2: {fig2_path}")

# ============================================================
# PASO 5: EXPORTAR ESTADISTICAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Exportar estadisticas")
print(f"{'-' * 70}")

csv_path = os.path.join(LC_DIR, "P16_cobertura_estadisticas.csv")
with open(csv_path, "w") as f:
    f.write("codigo,clase,area_km2,porcentaje\n")
    for code in sorted(class_areas.keys()):
        info = WC_CLASSES[code]
        area = class_areas[code]
        pct = area / total_area * 100
        f.write(f"{code},{info['name']},{area:.2f},{pct:.1f}\n")
print(f"  CSV: {csv_path}")

# ============================================================
# PASO 6: REPORTE
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Reporte")
print(f"{'-' * 70}")

reporte_path = os.path.join(LC_DIR, "P16_reporte_cobertura.txt")
with open(reporte_path, "w", encoding="utf-8") as rpt:
    rpt.write("=" * 70 + "\n")
    rpt.write("P16 -- COBERTURA Y USO DEL SUELO\n")
    rpt.write("=" * 70 + "\n")
    rpt.write(f"Fecha: {datetime.now()}\n")
    rpt.write(f"Fuente: ESA WorldCover 2021 v200 (10 m)\n\n")
    rpt.write(f"Area total AOI: {total_area:.2f} km2\n\n")
    rpt.write(f"{'Clase':<25} {'Area (km2)':>10} {'%':>8}\n")
    rpt.write(f"{'-'*25} {'-'*10} {'-'*8}\n")
    for code in sorted(class_areas.keys()):
        info = WC_CLASSES[code]
        area = class_areas[code]
        pct = area / total_area * 100
        rpt.write(f"{info['name']:<25} {area:>10.2f} {pct:>7.1f}%\n")
print(f"  Reporte: {reporte_path}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P16 COMPLETADO -- PASO 5.4 DE FASE 5")
print(f"{'=' * 70}")
print(f"  Area total: {total_area:.2f} km2")
# Clase dominante
dom_code = max(class_areas, key=class_areas.get)
dom_pct = class_areas[dom_code] / total_area * 100
print(f"  Clase dominante: {WC_CLASSES[dom_code]['name']} ({dom_pct:.1f}%)")
print(f"  GeoTIFF: P16_worldcover_2021.tif")
print(f"  Figuras: P16_fig1, P16_fig2")
print(f"\n  SIGUIENTE PASO: P17 (Paso 5.5 - Modelo de susceptibilidad)")
print(f"{'=' * 70}")