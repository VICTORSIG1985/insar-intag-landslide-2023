# -*- coding: utf-8 -*-
# =============================================================================
# P13 -- VARIABLES MORFOMETRICAS DEL DEM (FASE 5, PASO 5.1)
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 5 Integracion
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.0
#
# PROPOSITO:
#   Derivar variables morfometricas del DEM Copernicus GLO-30 para la
#   zona de Intag. Estas variables son predictores fundamentales en el
#   modelo de susceptibilidad a deslizamientos (Paso 5.5).
#
# VARIABLES A GENERAR:
#   1. Pendiente (grados)
#   2. Aspecto/Orientacion (grados)
#   3. Curvatura plana (plan curvature)
#   4. Curvatura de perfil (profile curvature)
#   5. TWI (Topographic Wetness Index)
#   6. TPI (Topographic Position Index)
#   7. Rugosidad del terreno (Terrain Roughness Index)
#
# DATOS DE ENTRADA:
#   - DEM de geometryGeo.h5 (Fase 3, generado por MintPy desde GLO-30)
#   - Alternativa: descarga directa de Copernicus GLO-30
#
# REFERENCIA:
#   - Guzzetti et al. (2012): Landslide inventory maps, Earth-Sci Rev
#   - Beven & Kirkby (1979): TWI, Hydrological Sciences Journal
#   - Wilson & Gallant (2000): Terrain Analysis, Wiley
# =============================================================================

import os
import sys
import numpy as np
from datetime import datetime

try:
    import h5py
    from scipy import ndimage
    import matplotlib.pyplot as plt
    from osgeo import gdal, osr
    gdal.UseExceptions()
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    sys.exit(1)

# ============================================================
# CONFIGURACION
# ============================================================
FASE5_DIR = r"D:\POSGRADOS\INTAG\data\fase5_integracion"
DEM_DIR = os.path.join(FASE5_DIR, "morfometria")
os.makedirs(DEM_DIR, exist_ok=True)

# DEM de MintPy (Fase 3)
GEOM_FILE = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3\mintpy\inputs\geometryGeo.h5"

# Inventarios para superposicion
INSAR_INVENTORY = r"D:\POSGRADOS\INTAG\data\sentinel1\dinsar_fase2\resultados\clusters_verificados_M_ASC_CO_12d.gpkg"
OPTICAL_INVENTORY = r"D:\POSGRADOS\INTAG\data\sentinel1\dinsar_fase2\resultados\P07c_inventario_optico.gpkg"

# AOI
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718

print("=" * 70)
print("P13 -- VARIABLES MORFOMETRICAS DEL DEM")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"Directorio: {DEM_DIR}")

# ============================================================
# PASO 1: CARGAR DEM
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Carga del DEM")
print(f"{'-' * 70}")

f = h5py.File(GEOM_FILE, "r")
dem = f["height"][:]
# Coordenadas
lat_first = 0.555718
lat_step = -0.0007210577689243029
lon_first = -79.285032
lon_step = 0.0007205173192771056
f.close()

ny, nx = dem.shape
lats = lat_first + np.arange(ny) * lat_step
lons = lon_first + np.arange(nx) * lon_step

# Tamano de pixel en metros (aproximado en el ecuador)
dx = lon_step * 111320 * np.cos(np.radians(np.mean([AOI_SOUTH, AOI_NORTH])))
dy = abs(lat_step) * 111320
print(f"  DEM shape: {ny} x {nx}")
print(f"  Elevacion: min={np.nanmin(dem):.0f} m, max={np.nanmax(dem):.0f} m")
print(f"  Tamano pixel: {dx:.1f} x {dy:.1f} m")

# ============================================================
# PASO 2: PENDIENTE (grados)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Pendiente")
print(f"{'-' * 70}")

# Gradiente en Y y X (en metros)
gy, gx = np.gradient(dem, dy, dx)
slope_rad = np.arctan(np.sqrt(gx**2 + gy**2))
slope_deg = np.degrees(slope_rad)

print(f"  Min: {np.nanmin(slope_deg):.2f} grados")
print(f"  Max: {np.nanmax(slope_deg):.2f} grados")
print(f"  Media: {np.nanmean(slope_deg):.2f} grados")
print(f"  Mediana: {np.nanmedian(slope_deg):.2f} grados")

# ============================================================
# PASO 3: ASPECTO (grados, 0=N, 90=E, 180=S, 270=W)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Aspecto (orientacion)")
print(f"{'-' * 70}")

aspect = np.degrees(np.arctan2(-gx, gy))
aspect = np.where(aspect < 0, aspect + 360, aspect)
print(f"  Min: {np.nanmin(aspect):.2f}")
print(f"  Max: {np.nanmax(aspect):.2f}")
print(f"  Media: {np.nanmean(aspect):.2f}")

# ============================================================
# PASO 4: CURVATURA PLANA Y DE PERFIL
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Curvatura plana y de perfil")
print(f"{'-' * 70}")

# Segundas derivadas
gyy, gyx = np.gradient(gy, dy, dx)
gxy, gxx = np.gradient(gx, dy, dx)

# Curvatura plana (plan curvature) - convergencia/divergencia del flujo
p2 = gx**2
q2 = gy**2
pq = p2 + q2
pq[pq == 0] = 1e-10

plan_curv = -(gxx * q2 - 2 * gx * gy * gxy + gyy * p2) / (pq**1.5)
# Limitar valores extremos
plan_curv = np.clip(plan_curv, -0.1, 0.1)

# Curvatura de perfil (profile curvature) - aceleracion del flujo
prof_curv = -(gxx * p2 + 2 * gx * gy * gxy + gyy * q2) / (pq**1.5)
prof_curv = np.clip(prof_curv, -0.1, 0.1)

print(f"  Curvatura plana: min={np.nanmin(plan_curv):.6f}, max={np.nanmax(plan_curv):.6f}")
print(f"  Curvatura perfil: min={np.nanmin(prof_curv):.6f}, max={np.nanmax(prof_curv):.6f}")

# ============================================================
# PASO 5: TWI (Topographic Wetness Index)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: TWI (Topographic Wetness Index)")
print(f"{'-' * 70}")

# TWI = ln(a / tan(b))
# a = area de contribucion (aproximada con filtro acumulativo)
# b = pendiente

# Aproximacion del area de contribucion usando flujo D8
# Simplificacion: usar filtro uniforme como proxy
# (Para TWI preciso se necesaria SAGA GIS o similar)
# Usamos una aproximacion basada en el tamano de celda
cell_area = dx * dy
slope_tan = np.tan(slope_rad)
slope_tan[slope_tan < 0.001] = 0.001  # Evitar log(inf)

# Area de contribucion aproximada con filtro acumulativo simple
# Usamos un kernel grande para simular la acumulacion
flow_acc = ndimage.uniform_filter(np.ones_like(dem), size=15) * (15**2) * cell_area

twi = np.log(flow_acc / slope_tan)
twi = np.clip(twi, 0, 25)

print(f"  Min: {np.nanmin(twi):.2f}")
print(f"  Max: {np.nanmax(twi):.2f}")
print(f"  Media: {np.nanmean(twi):.2f}")
print(f"  NOTA: TWI aproximado. Para precision, usar SAGA GIS con D8/D-inf.")

# ============================================================
# PASO 6: TPI (Topographic Position Index)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: TPI (Topographic Position Index)")
print(f"{'-' * 70}")

# TPI = elevacion - elevacion media en ventana circular
# Ventana de 500m (~6 pixeles)
radius_px = max(1, int(500 / dx))
tpi = dem - ndimage.uniform_filter(dem, size=2*radius_px+1)

print(f"  Radio: {radius_px} pixeles (~{radius_px*dx:.0f} m)")
print(f"  Min: {np.nanmin(tpi):.2f} m")
print(f"  Max: {np.nanmax(tpi):.2f} m")
print(f"  Media: {np.nanmean(tpi):.2f} m")

# ============================================================
# PASO 7: RUGOSIDAD (Terrain Roughness Index)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Rugosidad del terreno (TRI)")
print(f"{'-' * 70}")

# TRI = desviacion estandar de la elevacion en ventana 3x3
# Implementacion: sqrt(mean((dem - dem_mean)^2)) en ventana 3x3
dem_mean = ndimage.uniform_filter(dem, size=3)
dem_sq_mean = ndimage.uniform_filter(dem**2, size=3)
tri = np.sqrt(np.maximum(dem_sq_mean - dem_mean**2, 0))

print(f"  Min: {np.nanmin(tri):.2f} m")
print(f"  Max: {np.nanmax(tri):.2f} m")
print(f"  Media: {np.nanmean(tri):.2f} m")

# ============================================================
# PASO 8: GUARDAR COMO GeoTIFF
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 8: Exportar variables como GeoTIFF")
print(f"{'-' * 70}")

def save_geotiff(data, filename, description):
    filepath = os.path.join(DEM_DIR, filename)
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(filepath, nx, ny, 1, gdal.GDT_Float32,
                       options=["COMPRESS=LZW", "TILED=YES"])
    ds.SetGeoTransform([lon_first, lon_step, 0, lat_first, 0, lat_step])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(data.astype(np.float32))
    band.SetNoDataValue(-9999)
    band.SetDescription(description)
    ds.FlushCache()
    ds = None
    size_mb = os.path.getsize(filepath) / (1024**2)
    print(f"  {filename}: {size_mb:.1f} MB")
    return filepath

variables = [
    (dem, "P13_dem.tif", "Elevacion (m)"),
    (slope_deg, "P13_pendiente.tif", "Pendiente (grados)"),
    (aspect, "P13_aspecto.tif", "Aspecto (grados, 0=N)"),
    (plan_curv, "P13_curvatura_plana.tif", "Curvatura plana (1/m)"),
    (prof_curv, "P13_curvatura_perfil.tif", "Curvatura de perfil (1/m)"),
    (twi, "P13_twi.tif", "TWI - Topographic Wetness Index"),
    (tpi, "P13_tpi.tif", "TPI - Topographic Position Index (m)"),
    (tri, "P13_tri.tif", "TRI - Terrain Roughness Index (m)"),
]

for data, fname, desc in variables:
    save_geotiff(data, fname, desc)

# ============================================================
# PASO 9: FIGURAS DIAGNOSTICAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 9: Figuras diagnosticas")
print(f"{'-' * 70}")

extent = [lons[0], lons[-1], lats[-1], lats[0]]

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle("Variables Morfometricas del DEM - Zona de Intag, Cotacachi",
             fontsize=14, fontweight="bold")

configs = [
    (dem, "Elevacion (m)", "terrain", None, None),
    (slope_deg, "Pendiente (grados)", "YlOrRd", 0, 50),
    (aspect, "Aspecto (grados)", "hsv", 0, 360),
    (plan_curv, "Curvatura plana", "RdBu_r", -0.01, 0.01),
    (prof_curv, "Curvatura perfil", "RdBu_r", -0.01, 0.01),
    (twi, "TWI", "Blues", 5, 20),
    (tpi, "TPI (m)", "RdBu_r", -200, 200),
    (tri, "Rugosidad TRI (m)", "YlOrRd", 0, 50),
]

for ax, (data, title, cmap, vmin, vmax) in zip(axes.flat, configs):
    im = ax.imshow(data, cmap=cmap, extent=extent, aspect="auto",
                   vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("Longitud", fontsize=8)
    ax.set_ylabel("Latitud", fontsize=8)
    ax.tick_params(labelsize=7)
    plt.colorbar(im, ax=ax, shrink=0.8)

plt.tight_layout()
fig_path = os.path.join(DEM_DIR, "P13_variables_morfometricas.png")
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Figura: {fig_path}")

# ============================================================
# PASO 10: ESTADISTICAS RESUMEN
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 10: Estadisticas resumen")
print(f"{'-' * 70}")

print(f"\n  {'Variable':<25} {'Min':>10} {'Max':>10} {'Media':>10} {'Mediana':>10} {'Std':>10}")
print(f"  {'-'*25} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
stats_data = [
    ("Elevacion (m)", dem),
    ("Pendiente (grados)", slope_deg),
    ("Aspecto (grados)", aspect),
    ("Curvatura plana", plan_curv),
    ("Curvatura perfil", prof_curv),
    ("TWI", twi),
    ("TPI (m)", tpi),
    ("Rugosidad TRI (m)", tri),
]
for name, data in stats_data:
    v = data[~np.isnan(data)]
    print(f"  {name:<25} {np.min(v):>10.2f} {np.max(v):>10.2f} "
          f"{np.mean(v):>10.2f} {np.median(v):>10.2f} {np.std(v):>10.2f}")

# Distribucion de pendiente por clases
print(f"\n  Distribucion de pendiente:")
for low, high, label in [(0,5,"Plano"), (5,15,"Suave"), (15,25,"Moderado"),
                          (25,35,"Empinado"), (35,45,"Muy empinado"), (45,90,"Escarpado")]:
    n = np.sum((slope_deg >= low) & (slope_deg < high))
    pct = n / slope_deg.size * 100
    print(f"    {label} ({low}-{high} grados): {n:,} pixeles ({pct:.1f}%)")

# ============================================================
# PASO 11: REPORTE
# ============================================================
reporte_path = os.path.join(DEM_DIR, "P13_reporte_morfometria.txt")
with open(reporte_path, "w", encoding="utf-8") as rpt:
    rpt.write("=" * 70 + "\n")
    rpt.write("P13 -- VARIABLES MORFOMETRICAS DEL DEM\n")
    rpt.write("=" * 70 + "\n")
    rpt.write(f"Fecha: {datetime.now()}\n\n")
    rpt.write(f"DEM: Copernicus GLO-30 (via geometryGeo.h5)\n")
    rpt.write(f"Dimensiones: {ny} x {nx} pixeles\n")
    rpt.write(f"Tamano pixel: {dx:.1f} x {dy:.1f} m\n")
    rpt.write(f"Elevacion: {np.nanmin(dem):.0f} - {np.nanmax(dem):.0f} m\n\n")
    rpt.write("Variables generadas:\n")
    for _, fname, desc in variables:
        rpt.write(f"  {fname}: {desc}\n")
    rpt.write(f"\nFigura: P13_variables_morfometricas.png\n")
print(f"\n  Reporte: {reporte_path}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P13 COMPLETADO -- PASO 5.1 DE FASE 5")
print(f"{'=' * 70}")
print(f"  Variables generadas: {len(variables)}")
print(f"  GeoTIFFs en: {DEM_DIR}")
print(f"  Figura: P13_variables_morfometricas.png")
print(f"\n  SIGUIENTE PASO: P14 (Paso 5.2 - Precipitacion como disparador)")
print(f"  Requiere: Google Earth Engine (CHIRPS / ERA5-Land)")
print(f"{'=' * 70}")