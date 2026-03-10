# -*- coding: utf-8 -*-
# =============================================================================
# P20c -- MAPA DE SUSCEPTIBILIDAD SOBRE OpenStreetMap
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-09
# Version  : 1.0
#
# PROPOSITO:
#   Superponer el mapa de susceptibilidad sobre OpenStreetMap para que
#   se vean carreteras, pueblos, rios y topografia con las zonas de
#   riesgo claramente identificables.
# =============================================================================

import os
import sys
import numpy as np
from datetime import datetime

try:
    import geopandas as gpd
    import contextily as cx
    from osgeo import gdal, osr, ogr
    gdal.UseExceptions()
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import ListedColormap
    from pyproj import Transformer
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    sys.exit(1)

# ============================================================
# CONFIGURACION
# ============================================================
AOI_GPKG = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"
SUSC_PATH = r"D:\POSGRADOS\INTAG\data\fase5_integracion\susceptibilidad\P17_susceptibilidad_probabilidad.tif"
CLASS_PATH = r"D:\POSGRADOS\INTAG\data\fase5_integracion\susceptibilidad\P17_susceptibilidad_clases.tif"
INV_PATH = r"D:\POSGRADOS\INTAG\data\sentinel2\VALIDACION\INTAG_InSAR_P07\validacion_cruzada\P07d_cruce_insar.gpkg"
PUB_DIR = r"D:\POSGRADOS\INTAG\data\fase5_integracion\publicacion"
os.makedirs(PUB_DIR, exist_ok=True)

print("=" * 70)
print("P20c -- MAPA SUSCEPTIBILIDAD SOBRE OpenStreetMap")
print("=" * 70)
print(f"Fecha: {datetime.now()}")

# ============================================================
# PASO 1: CARGAR DATOS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Carga de datos")
print(f"{'-' * 70}")

gdf = gpd.read_file(AOI_GPKG)
if gdf.crs.to_epsg() != 4326:
    gdf = gdf.to_crs(4326)
gdf_3857 = gdf.to_crs(3857)

# Centroides
parroquias = {}
for _, row in gdf_3857.iterrows():
    parroquias[row["DPA_DESPAR"]] = (row.geometry.centroid.x, row.geometry.centroid.y)

# Susceptibilidad
ds = gdal.Open(CLASS_PATH)
classes = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
gt = ds.GetGeoTransform()
ny, nx = classes.shape
ds = None

ds = gdal.Open(SUSC_PATH)
prob = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
ds = None

# Mascara poligono
driver_mem = gdal.GetDriverByName("MEM")
ds_mask = driver_mem.Create("", nx, ny, 1, gdal.GDT_Byte)
ds_mask.SetGeoTransform(gt)
srs_mask = osr.SpatialReference()
srs_mask.ImportFromEPSG(4326)
ds_mask.SetProjection(srs_mask.ExportToWkt())
ds_mask.GetRasterBand(1).Fill(0)
ogr_drv = ogr.GetDriverByName("Memory")
ogr_ds = ogr_drv.CreateDataSource("")
ogr_lyr = ogr_ds.CreateLayer("mask", srs_mask, ogr.wkbPolygon)
for _, row in gdf.iterrows():
    feat_ogr = ogr.Feature(ogr_lyr.GetLayerDefn())
    geom_ogr = ogr.CreateGeometryFromWkt(row.geometry.wkt)
    feat_ogr.SetGeometry(geom_ogr)
    ogr_lyr.CreateFeature(feat_ogr)
gdal.RasterizeLayer(ds_mask, [1], ogr_lyr, burn_values=[1])
mask = ds_mask.GetRasterBand(1).ReadAsArray().astype(bool)
ds_mask = None
ogr_ds = None

classes_masked = np.where(mask, classes, np.nan)

# Inventario de deslizamientos
inv = gpd.read_file(INV_PATH)
inv = inv.to_crs(3857)
print(f"  Deslizamientos: {len(inv)}")

# Extent en Web Mercator
transformer = Transformer.from_crs(4326, 3857, always_xy=True)
x_min_3857, y_min_3857 = transformer.transform(gt[0], gt[3] + ny * gt[5])
x_max_3857, y_max_3857 = transformer.transform(gt[0] + nx * gt[1], gt[3])
raster_extent_3857 = [x_min_3857, x_max_3857, y_min_3857, y_max_3857]

total_bounds = gdf_3857.total_bounds
margin = 2000

print(f"  Datos cargados OK")

# ============================================================
# PASO 2: MAPA 1 - SUSCEPTIBILIDAD SOBRE OSM
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Mapa susceptibilidad sobre OSM")
print(f"{'-' * 70}")

fig, ax = plt.subplots(figsize=(16, 14))

# Susceptibilidad con transparencia
susc_display = np.full_like(classes_masked, np.nan)
susc_display[classes_masked == 2] = 2  # Media
susc_display[classes_masked == 3] = 3  # Alta
susc_display[classes_masked == 4] = 4  # Muy alta

cmap_3 = ListedColormap(["#FFD700", "#FF6600", "#CC0000"])
ax.imshow(susc_display, cmap=cmap_3, extent=raster_extent_3857,
          aspect="auto", alpha=0.5, vmin=1.5, vmax=4.5, zorder=2,
          interpolation="nearest")

# Contorno parroquias
gdf_3857.boundary.plot(ax=ax, color="black", linewidth=2, zorder=3)

# Nombres
for nombre, (x, y) in parroquias.items():
    nombre_display = nombre.replace("SEIS DE JULIO DE ", "")
    nombre_display = nombre_display.replace("GARCIA MORENO", "GARCÍA MORENO")
    nombre_display = nombre_display.replace("PENAHERRERA", "PEÑAHERRERA")
    nombre_display = nombre_display.replace("GUTIERREZ", "GUTIÉRREZ")
    ax.text(x, y, nombre_display, fontsize=11, fontweight="bold",
            color="black", ha="center", va="center", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      alpha=0.85, edgecolor="black", linewidth=1))

ax.set_xlim(total_bounds[0] - margin, total_bounds[2] + margin)
ax.set_ylim(total_bounds[1] - margin, total_bounds[3] + margin)

print("  Descargando basemap OSM...")
cx.add_basemap(ax, source=cx.providers.OpenStreetMap.Mapnik, zoom=12,
               attribution="© OpenStreetMap contributors")

ax.set_title("MAPA DE SUSCEPTIBILIDAD A DESLIZAMIENTOS\n"
             "Zona de Intag, Cantón Cotacachi, Imbabura, Ecuador\n"
             "Basado en evento del 19 de diciembre de 2023",
             fontsize=14, fontweight="bold")

patches = [
    mpatches.Patch(color="#FFD700", alpha=0.6, label="Susceptibilidad MEDIA (25.2%)"),
    mpatches.Patch(color="#FF6600", alpha=0.6, label="Susceptibilidad ALTA (3.7%)"),
    mpatches.Patch(color="#CC0000", alpha=0.6, label="Susceptibilidad MUY ALTA (0.1%)"),
    mpatches.Patch(facecolor="none", edgecolor="black", linewidth=2,
                   label="Límites parroquiales"),
]
ax.legend(handles=patches, loc="lower left", fontsize=10, framealpha=0.95)

info = (f"Área: 1,528 km2 | 6 parroquias\n"
        f"Alta+Muy Alta: 3.7% (~57 km2)\n"
        f"Modelo: Random Forest AUC=0.67\n"
        f"Fuente: Sentinel-1 InSAR + DEM\n"
        f"Marzo 2026")
ax.text(0.98, 0.02, info, transform=ax.transAxes, fontsize=8,
        va="bottom", ha="right", family="monospace",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9,
                  edgecolor="black"), zorder=5)

ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
path1 = os.path.join(PUB_DIR, "MAPA_OSM_SUSCEPTIBILIDAD.png")
plt.savefig(path1, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Mapa 1: {path1}")

# ============================================================
# PASO 3: MAPA 2 - VALIDACION (deslizamientos + susceptibilidad)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Mapa validacion deslizamientos vs susceptibilidad")
print(f"{'-' * 70}")

fig, ax = plt.subplots(figsize=(16, 14))

# Susceptibilidad 5 clases con transparencia
cmap_5 = ListedColormap(["#2166AC", "#67A9CF", "#FDDBC7", "#EF8A62", "#B2182B"])
ax.imshow(classes_masked, cmap=cmap_5, extent=raster_extent_3857,
          aspect="auto", alpha=0.45, vmin=-0.5, vmax=4.5, zorder=2,
          interpolation="nearest")

# Deslizamientos como puntos
inv.plot(ax=ax, color="black", markersize=3, alpha=0.6, zorder=4,
         label=f"Deslizamientos (n={len(inv)})")

# Contorno
gdf_3857.boundary.plot(ax=ax, color="black", linewidth=1.5, zorder=3)

# Nombres
for nombre, (x, y) in parroquias.items():
    nombre_display = nombre.replace("SEIS DE JULIO DE ", "")
    nombre_display = nombre_display.replace("GARCIA MORENO", "GARCÍA MORENO")
    nombre_display = nombre_display.replace("PENAHERRERA", "PEÑAHERRERA")
    nombre_display = nombre_display.replace("GUTIERREZ", "GUTIÉRREZ")
    ax.text(x, y, nombre_display, fontsize=9, fontweight="bold",
            color="black", ha="center", va="center", zorder=5,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      alpha=0.8, edgecolor="gray", linewidth=0.5))

ax.set_xlim(total_bounds[0] - margin, total_bounds[2] + margin)
ax.set_ylim(total_bounds[1] - margin, total_bounds[3] + margin)

print("  Descargando basemap OSM...")
cx.add_basemap(ax, source=cx.providers.OpenStreetMap.Mapnik, zoom=12,
               attribution="© OpenStreetMap contributors")

ax.set_title("VALIDACIÓN: DEFORMACIÓN DETECTADA vs SUSCEPTIBILIDAD\n"
             "1,481 zonas con movimiento del terreno detectado - 19/12/2023 (puntos negros)\n"
             "82.2% cayeron en zonas de susceptibilidad Media, Alta o Muy Alta",
             fontsize=13, fontweight="bold")

patches_v = [
    mpatches.Patch(color="#2166AC", alpha=0.5, label="Muy baja (4.9% de eventos)"),
    mpatches.Patch(color="#67A9CF", alpha=0.5, label="Baja (12.8%)"),
    mpatches.Patch(color="#FDDBC7", alpha=0.5, label="Media (17.4%)"),
    mpatches.Patch(color="#EF8A62", alpha=0.5, label="Alta (57.0%)"),
    mpatches.Patch(color="#B2182B", alpha=0.5, label="Muy alta (7.9%)"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
               markersize=6, label=f"Zonas de deformación (n=1,481)"),
]
ax.legend(handles=patches_v, loc="lower left", fontsize=9, framealpha=0.95,
          title="Clase de susceptibilidad (% de eventos reales)",
          title_fontsize=10)

# Tabla resumen
tabla_text = (
    "RESUMEN DE VALIDACIÓN\n"
    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    "Clase         Eventos    %\n"
    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    "Muy baja         73    4.9%\n"
    "Baja            190   12.8%\n"
    "Media           257   17.4%\n"
    "Alta            844   57.0%\n"
    "Muy alta        117    7.9%\n"
    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    "Media+Alta+MAlta: 82.2%\n"
    "Alta+Muy alta: 64.9%\n"
    "(en solo 3.7% del área)\n"
    "Marzo 2026"
)
ax.text(0.98, 0.02, tabla_text, transform=ax.transAxes, fontsize=8,
        va="bottom", ha="right", family="monospace",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.95,
                  edgecolor="black"), zorder=5)

ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
path2 = os.path.join(PUB_DIR, "MAPA_OSM_VALIDACION_DESLIZAMIENTOS.png")
plt.savefig(path2, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Mapa 2: {path2}")

# ============================================================
# PASO 4: MAPA 3 - ZONAS CRITICAS SOBRE OSM
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Mapa zonas criticas sobre OSM")
print(f"{'-' * 70}")

fig, ax = plt.subplots(figsize=(16, 14))

# Solo alta + muy alta
critico = np.full_like(classes_masked, np.nan)
critico[classes_masked >= 3] = classes_masked[classes_masked >= 3]

cmap_crit = ListedColormap(["#FF6600", "#CC0000"])
ax.imshow(critico, cmap=cmap_crit, extent=raster_extent_3857,
          aspect="auto", alpha=0.6, vmin=2.5, vmax=4.5, zorder=2,
          interpolation="nearest")

gdf_3857.boundary.plot(ax=ax, color="black", linewidth=2, zorder=3)

for nombre, (x, y) in parroquias.items():
    nombre_display = nombre.replace("SEIS DE JULIO DE ", "")
    nombre_display = nombre_display.replace("GARCIA MORENO", "GARCÍA MORENO")
    nombre_display = nombre_display.replace("PENAHERRERA", "PEÑAHERRERA")
    nombre_display = nombre_display.replace("GUTIERREZ", "GUTIÉRREZ")
    ax.text(x, y, nombre_display, fontsize=11, fontweight="bold",
            color="black", ha="center", va="center", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      alpha=0.85, edgecolor="red", linewidth=1.5))

ax.set_xlim(total_bounds[0] - margin, total_bounds[2] + margin)
ax.set_ylim(total_bounds[1] - margin, total_bounds[3] + margin)

print("  Descargando basemap OSM...")
cx.add_basemap(ax, source=cx.providers.OpenStreetMap.Mapnik, zoom=12,
               attribution="© OpenStreetMap contributors")

ax.set_title("ZONAS DE RIESGO PRIORITARIO\n"
             "Susceptibilidad ALTA y MUY ALTA a movimientos del terreno\n"
             "3.7% del área (~57 km2) concentra el 64.9% de las deformaciones detectadas",
             fontsize=14, fontweight="bold")

patches_c = [
    mpatches.Patch(color="#FF6600", alpha=0.7, label="ALTA - Monitoreo recomendado"),
    mpatches.Patch(color="#CC0000", alpha=0.7, label="MUY ALTA - Acción prioritaria"),
]
ax.legend(handles=patches_c, loc="lower left", fontsize=11, framealpha=0.95)

info_c = "Marzo 2026"
ax.text(0.98, 0.02, info_c, transform=ax.transAxes, fontsize=9,
        va="bottom", ha="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9,
                  edgecolor="black"), zorder=5)

ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
path3 = os.path.join(PUB_DIR, "MAPA_OSM_ZONAS_CRITICAS.png")
plt.savefig(path3, dpi=300, bbox_inches="tight")
plt.close()
print(f"  Mapa 3: {path3}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P20c COMPLETADO")
print(f"{'=' * 70}")
print(f"  1. MAPA_OSM_SUSCEPTIBILIDAD.png (general sobre OSM)")
print(f"  2. MAPA_OSM_VALIDACION_DESLIZAMIENTOS.png (eventos reales vs modelo)")
print(f"  3. MAPA_OSM_ZONAS_CRITICAS.png (zonas prioritarias para autoridades)")
print(f"{'=' * 70}")
