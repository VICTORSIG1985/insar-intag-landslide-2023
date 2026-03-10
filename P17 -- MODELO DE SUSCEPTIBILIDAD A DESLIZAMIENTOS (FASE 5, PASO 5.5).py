# -*- coding: utf-8 -*-
# =============================================================================
# P17 -- MODELO DE SUSCEPTIBILIDAD A DESLIZAMIENTOS (FASE 5, PASO 5.5)
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 5 Integracion
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.2 (filtrado al poligono real INTAJ.gpkg)
#
# PROPOSITO:
#   Generar un mapa de susceptibilidad a deslizamientos integrando
#   variables morfometricas, cobertura del suelo, y los inventarios
#   InSAR y optico mediante un modelo Random Forest.
#
# REFERENCIA:
#   - Protocolo Maestro InSAR Intag v2.0, Paso 5.5
#   - Breiman (2001): Random Forests, Machine Learning 45:5-32
#   - Reichenbach et al. (2018): Statistical landslide susceptibility
# =============================================================================

import os
import sys
import numpy as np
from datetime import datetime

try:
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import (roc_auc_score, roc_curve,
                                  classification_report, confusion_matrix)
    import matplotlib.pyplot as plt
    import geopandas as gpd
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    print("  Instalar: pip install scikit-learn")
    sys.exit(1)

# ============================================================
# CONFIGURACION
# ============================================================
FASE5_DIR = r"D:\POSGRADOS\INTAG\data\fase5_integracion"
MORF_DIR = os.path.join(FASE5_DIR, "morfometria")
LC_DIR = os.path.join(FASE5_DIR, "cobertura")
MODEL_DIR = os.path.join(FASE5_DIR, "susceptibilidad")
os.makedirs(MODEL_DIR, exist_ok=True)

# Inventarios de deslizamientos (rutas corregidas)
INSAR_CRUCE = r"D:\POSGRADOS\INTAG\data\sentinel2\VALIDACION\INTAG_InSAR_P07\validacion_cruzada\P07d_cruce_insar.gpkg"
INSAR_VERIF = r"D:\POSGRADOS\INTAG\data\sentinel1\hyp3_deformation\verificacion\clusters_verificados_M_ASC_CO_12d.gpkg"

# AOI
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718

# Rasters predictores
RASTERS = {
    "elevacion": os.path.join(MORF_DIR, "P13_dem.tif"),
    "pendiente": os.path.join(MORF_DIR, "P13_pendiente.tif"),
    "aspecto": os.path.join(MORF_DIR, "P13_aspecto.tif"),
    "curv_plana": os.path.join(MORF_DIR, "P13_curvatura_plana.tif"),
    "curv_perfil": os.path.join(MORF_DIR, "P13_curvatura_perfil.tif"),
    "twi": os.path.join(MORF_DIR, "P13_twi.tif"),
    "tpi": os.path.join(MORF_DIR, "P13_tpi.tif"),
    "rugosidad": os.path.join(MORF_DIR, "P13_tri.tif"),
    "cobertura": os.path.join(LC_DIR, "P16_worldcover_2021.tif"),
}

# Poligono real del AOI (6 parroquias de Intag)
AOI_GPKG = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"

# Random Forest
N_TREES = 500
TEST_SIZE = 0.3
RANDOM_STATE = 42
N_NEGATIVE_RATIO = 3

print("=" * 70)
print("P17 -- MODELO DE SUSCEPTIBILIDAD A DESLIZAMIENTOS")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"Metodo: Random Forest ({N_TREES} arboles)")
print(f"Validacion: {int((1-TEST_SIZE)*100)}/{int(TEST_SIZE*100)} train/test")

# ============================================================
# PASO 1: CARGAR RASTERS PREDICTORES
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Carga de rasters predictores")
print(f"{'-' * 70}")

raster_data = {}
ref_shape = None
ref_gt = None

for name, path in RASTERS.items():
    if not os.path.exists(path):
        print(f"  ADVERTENCIA: {name} no encontrado en {path}")
        continue
    ds = gdal.Open(path)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray().astype(np.float32)
    gt = ds.GetGeoTransform()
    ds = None

    if ref_shape is None:
        ref_shape = data.shape
        ref_gt = gt
        print(f"  Referencia: {ref_shape}, GT={gt[:2]}...{gt[4:6]}")

    if data.shape != ref_shape:
        from scipy.ndimage import zoom
        zoom_y = ref_shape[0] / data.shape[0]
        zoom_x = ref_shape[1] / data.shape[1]
        if name == "cobertura":
            data = zoom(data, (zoom_y, zoom_x), order=0)
        else:
            data = zoom(data, (zoom_y, zoom_x), order=1)
        print(f"  {name}: redimensionado a {ref_shape}")

    raster_data[name] = data
    print(f"  {name}: {data.shape}, min={np.nanmin(data):.2f}, max={np.nanmax(data):.2f}")

n_predictores = len(raster_data)
print(f"\n  Total predictores: {n_predictores}")

# Crear mascara del poligono real
print(f"\n  Creando mascara del poligono INTAJ.gpkg...")
gdf_aoi = gpd.read_file(AOI_GPKG)
if gdf_aoi.crs.to_epsg() != 4326:
    gdf_aoi = gdf_aoi.to_crs(4326)

# Rasterizar poligono
driver_mem = gdal.GetDriverByName("MEM")
ds_mask = driver_mem.Create("", ref_shape[1], ref_shape[0], 1, gdal.GDT_Byte)
ds_mask.SetGeoTransform(ref_gt)
srs_mask = osr.SpatialReference()
srs_mask.ImportFromEPSG(4326)
ds_mask.SetProjection(srs_mask.ExportToWkt())
ds_mask.GetRasterBand(1).Fill(0)

ogr_drv = ogr.GetDriverByName("Memory")
ogr_ds = ogr_drv.CreateDataSource("")
ogr_lyr = ogr_ds.CreateLayer("mask", srs_mask, ogr.wkbPolygon)
for _, row in gdf_aoi.iterrows():
    feat_ogr = ogr.Feature(ogr_lyr.GetLayerDefn())
    geom_ogr = ogr.CreateGeometryFromWkt(row.geometry.wkt)
    feat_ogr.SetGeometry(geom_ogr)
    ogr_lyr.CreateFeature(feat_ogr)
gdal.RasterizeLayer(ds_mask, [1], ogr_lyr, burn_values=[1])
polygon_mask = ds_mask.GetRasterBand(1).ReadAsArray().astype(bool)
ds_mask = None
ogr_ds = None

print(f"  Pixeles dentro del poligono: {np.sum(polygon_mask):,} ({np.sum(polygon_mask)/polygon_mask.size:.1%})")
area_km2 = gdf_aoi.to_crs(32618).area.sum() / 1e6
print(f"  Area real: {area_km2:.1f} km2")

# ============================================================
# PASO 2: CARGAR INVENTARIO DE DESLIZAMIENTOS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Carga del inventario de deslizamientos")
print(f"{'-' * 70}")

inv_path = None
for candidate in [INSAR_CRUCE, INSAR_VERIF]:
    if os.path.exists(candidate):
        inv_path = candidate
        break

if inv_path is None:
    print(f"  ERROR: No se encontro inventario de deslizamientos.")
    print(f"  Buscados: {INSAR_CRUCE}")
    print(f"            {INSAR_VERIF}")
    sys.exit(1)

print(f"  Inventario: {os.path.basename(inv_path)}")

ds_inv = ogr.Open(inv_path)
layer = ds_inv.GetLayer()
n_features = layer.GetFeatureCount()
srs_layer = layer.GetSpatialRef()
print(f"  Features: {n_features}")
if srs_layer:
    print(f"  CRS: {srs_layer.GetAuthorityCode(None)}")

slide_coords = []
for feat in layer:
    geom = feat.GetGeometryRef()
    if geom is not None:
        if geom.GetGeometryType() in [ogr.wkbPoint, ogr.wkbPoint25D]:
            x, y = geom.GetX(), geom.GetY()
        else:
            centroid = geom.Centroid()
            x, y = centroid.GetX(), centroid.GetY()

        if abs(x) > 360:
            src_srs = srs_layer if srs_layer else osr.SpatialReference()
            if not srs_layer:
                src_srs.ImportFromEPSG(32618)
            dst_srs = osr.SpatialReference()
            dst_srs.ImportFromEPSG(4326)
            dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            transform = osr.CoordinateTransformation(src_srs, dst_srs)
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(x, y)
            point.Transform(transform)
            x, y = point.GetX(), point.GetY()

        if AOI_WEST <= x <= AOI_EAST and AOI_SOUTH <= y <= AOI_NORTH:
            # Verificar que esta dentro del poligono real
            row_px, col_px = int((y - ref_gt[3]) / ref_gt[5]), int((x - ref_gt[0]) / ref_gt[1])
            if 0 <= row_px < ref_shape[0] and 0 <= col_px < ref_shape[1]:
                if polygon_mask[row_px, col_px]:
                    slide_coords.append((x, y))

ds_inv = None
print(f"  Deslizamientos en AOI: {len(slide_coords)}")

# ============================================================
# PASO 3: CONVERTIR COORDENADAS A PIXELES
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Extraccion de valores en puntos de deslizamiento")
print(f"{'-' * 70}")

ny, nx = ref_shape

def coord_to_pixel(x, y, gt):
    col = int((x - gt[0]) / gt[1])
    row = int((y - gt[3]) / gt[5])
    return row, col

slide_values = []
for x, y in slide_coords:
    row, col = coord_to_pixel(x, y, ref_gt)
    if 0 <= row < ny and 0 <= col < nx:
        vals = {}
        valid = True
        for name, data in raster_data.items():
            val = data[row, col]
            if np.isnan(val) or val == -9999:
                valid = False
                break
            vals[name] = val
        if valid:
            slide_values.append(vals)

print(f"  Puntos de deslizamiento con valores validos: {len(slide_values)}")

n_negative = len(slide_values) * N_NEGATIVE_RATIO
print(f"  Generando {n_negative} puntos no-deslizamiento (ratio {N_NEGATIVE_RATIO}:1)...")

np.random.seed(RANDOM_STATE)
neg_values = []
max_attempts = n_negative * 20
attempts = 0

slide_mask = np.zeros((ny, nx), dtype=bool)
for x, y in slide_coords:
    row, col = coord_to_pixel(x, y, ref_gt)
    if 0 <= row < ny and 0 <= col < nx:
        r_min = max(0, row - 3)
        r_max = min(ny, row + 4)
        c_min = max(0, col - 3)
        c_max = min(nx, col + 4)
        slide_mask[r_min:r_max, c_min:c_max] = True

while len(neg_values) < n_negative and attempts < max_attempts:
    row = np.random.randint(0, ny)
    col = np.random.randint(0, nx)
    attempts += 1

    # Solo dentro del poligono real
    if not polygon_mask[row, col]:
        continue

    if slide_mask[row, col]:
        continue

    vals = {}
    valid = True
    for name, data in raster_data.items():
        val = data[row, col]
        if np.isnan(val) or val == -9999 or val == 0:
            valid = False
            break
        vals[name] = val
    if valid:
        neg_values.append(vals)

print(f"  Puntos no-deslizamiento generados: {len(neg_values)}")

# ============================================================
# PASO 4: PREPARAR DATOS PARA RANDOM FOREST
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Preparacion de datos para Random Forest")
print(f"{'-' * 70}")

var_names = sorted(raster_data.keys())
print(f"  Variables: {var_names}")

X_pos = np.array([[v[name] for name in var_names] for v in slide_values])
X_neg = np.array([[v[name] for name in var_names] for v in neg_values])

X = np.vstack([X_pos, X_neg])
y = np.hstack([np.ones(len(X_pos)), np.zeros(len(X_neg))])

print(f"  Positivos (deslizamiento): {len(X_pos)}")
print(f"  Negativos (no-deslizamiento): {len(X_neg)}")
print(f"  Total muestras: {len(X)}")

asp_idx = var_names.index("aspecto")
asp_rad = np.radians(X[:, asp_idx])
X_final = np.column_stack([
    X[:, :asp_idx],
    np.sin(asp_rad),
    np.cos(asp_rad),
    X[:, asp_idx+1:]
])
var_names_final = (var_names[:asp_idx] +
                   ["aspecto_sin", "aspecto_cos"] +
                   var_names[asp_idx+1:])
print(f"  Variables finales ({len(var_names_final)}): {var_names_final}")

X_train, X_test, y_train, y_test = train_test_split(
    X_final, y, test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y
)
print(f"  Train: {len(X_train)} | Test: {len(X_test)}")

# ============================================================
# PASO 5: ENTRENAR RANDOM FOREST
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Entrenamiento Random Forest")
print(f"{'-' * 70}")

rf = RandomForestClassifier(
    n_estimators=N_TREES,
    max_depth=None,
    min_samples_split=5,
    min_samples_leaf=2,
    class_weight="balanced",
    random_state=RANDOM_STATE,
    n_jobs=-1,
    oob_score=True,
)

print(f"  Entrenando con {N_TREES} arboles...")
rf.fit(X_train, y_train)
print(f"  OOB Score: {rf.oob_score_:.4f}")

importances = rf.feature_importances_
sorted_idx = np.argsort(importances)[::-1]
print(f"\n  Importancia de variables:")
for i in sorted_idx:
    print(f"    {var_names_final[i]:<20} {importances[i]:.4f}")

# ============================================================
# PASO 6: VALIDACION
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Validacion")
print(f"{'-' * 70}")

y_pred = rf.predict(X_test)
y_prob = rf.predict_proba(X_test)[:, 1]

auc = roc_auc_score(y_test, y_prob)
print(f"  AUC-ROC: {auc:.4f}")

print(f"\n  Reporte de clasificacion:")
report = classification_report(y_test, y_pred,
                                target_names=["No desliz.", "Desliz."])
print(report)

cm = confusion_matrix(y_test, y_pred)
print(f"  Matriz de confusion:")
print(f"    {'':>15} {'Pred No':>10} {'Pred Si':>10}")
print(f"    {'Real No':<15} {cm[0,0]:>10} {cm[0,1]:>10}")
print(f"    {'Real Si':<15} {cm[1,0]:>10} {cm[1,1]:>10}")

# ============================================================
# PASO 7: GENERAR MAPA DE SUSCEPTIBILIDAD
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Generacion del mapa de susceptibilidad")
print(f"{'-' * 70}")

print("  Prediciendo probabilidad pixel por pixel...")

asp_data = raster_data["aspecto"]
asp_sin = np.sin(np.radians(asp_data))
asp_cos = np.cos(np.radians(asp_data))

stack_layers = []
for name in var_names_final:
    if name == "aspecto_sin":
        stack_layers.append(asp_sin)
    elif name == "aspecto_cos":
        stack_layers.append(asp_cos)
    else:
        stack_layers.append(raster_data[name])

stack = np.stack(stack_layers, axis=-1)
stack_flat = stack.reshape(-1, len(var_names_final))

valid_mask = np.all(~np.isnan(stack_flat), axis=1) & np.all(stack_flat != -9999, axis=1)

susceptibility_flat = np.full(stack_flat.shape[0], np.nan)
if np.sum(valid_mask) > 0:
    print(f"  Pixeles a predecir: {np.sum(valid_mask):,}")
    block_size = 100000
    valid_indices = np.where(valid_mask)[0]
    for i in range(0, len(valid_indices), block_size):
        block_idx = valid_indices[i:i+block_size]
        block_data = stack_flat[block_idx]
        block_data = np.nan_to_num(block_data, nan=0, posinf=0, neginf=0)
        probs = rf.predict_proba(block_data)[:, 1]
        susceptibility_flat[block_idx] = probs
        if (i // block_size) % 5 == 0:
            print(f"    Bloque {i//block_size + 1}...")

susceptibility_map = susceptibility_flat.reshape(ny, nx)
# Aplicar mascara del poligono real
susceptibility_map = np.where(polygon_mask, susceptibility_map, np.nan)
print(f"  Mapa generado: {ny} x {nx} (enmascarado al poligono real)")
print(f"  Susceptibilidad: min={np.nanmin(susceptibility_map):.4f}, "
      f"max={np.nanmax(susceptibility_map):.4f}, "
      f"media={np.nanmean(susceptibility_map):.4f}")

print(f"\n  Clasificacion en 5 niveles:")
classes = np.full_like(susceptibility_map, np.nan)
thresholds = [0.2, 0.4, 0.6, 0.8]
labels = ["Muy baja", "Baja", "Media", "Alta", "Muy alta"]

classes[susceptibility_map < thresholds[0]] = 0
classes[(susceptibility_map >= thresholds[0]) & (susceptibility_map < thresholds[1])] = 1
classes[(susceptibility_map >= thresholds[1]) & (susceptibility_map < thresholds[2])] = 2
classes[(susceptibility_map >= thresholds[2]) & (susceptibility_map < thresholds[3])] = 3
classes[susceptibility_map >= thresholds[3]] = 4

for i, label in enumerate(labels):
    n = np.sum(classes == i)
    pct = n / np.sum(~np.isnan(classes)) * 100 if np.sum(~np.isnan(classes)) > 0 else 0
    print(f"    {label:<12} {n:>8,} pixeles ({pct:.1f}%)")

# ============================================================
# PASO 8: GUARDAR RASTERS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 8: Exportar rasters")
print(f"{'-' * 70}")

def save_raster(data, filename, dtype=gdal.GDT_Float32):
    filepath = os.path.join(MODEL_DIR, filename)
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(filepath, nx, ny, 1, dtype,
                       options=["COMPRESS=LZW", "TILED=YES"])
    ds.SetGeoTransform(ref_gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(data.astype(np.float32))
    band.SetNoDataValue(-9999)
    ds.FlushCache()
    ds = None
    size = os.path.getsize(filepath) / (1024**2)
    print(f"  {filename}: {size:.1f} MB")

save_raster(susceptibility_map, "P17_susceptibilidad_probabilidad.tif")
save_raster(classes, "P17_susceptibilidad_clases.tif", gdal.GDT_Byte)

# ============================================================
# PASO 9: FIGURAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 9: Generacion de figuras")
print(f"{'-' * 70}")

extent = [AOI_WEST, AOI_EAST, AOI_SOUTH, AOI_NORTH]

# Funcion para dibujar contorno del poligono
def plot_polygon_border(ax, gdf_local, color="black", lw=1):
    for _, row in gdf_local.iterrows():
        if row.geometry.geom_type == "MultiPolygon":
            for poly in row.geometry.geoms:
                x, y = poly.exterior.xy
                ax.plot(x, y, color=color, linewidth=lw)
        elif row.geometry.geom_type == "Polygon":
            x, y = row.geometry.exterior.xy
            ax.plot(x, y, color=color, linewidth=lw)

# FIGURA 1: Mapa + ROC
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle("Modelo de Susceptibilidad a Deslizamientos - Intag, Cotacachi\n"
             f"Random Forest ({N_TREES} arboles) | AUC = {auc:.4f} | AOI: {area_km2:.0f} km2",
             fontsize=13, fontweight="bold")

from matplotlib.colors import ListedColormap
cmap_susc = ListedColormap(["#2166AC", "#67A9CF", "#FDDBC7", "#EF8A62", "#B2182B"])
im = axes[0].imshow(classes, cmap=cmap_susc, extent=extent, aspect="auto",
                     vmin=-0.5, vmax=4.5, interpolation="nearest")
plot_polygon_border(axes[0], gdf_aoi, color="black", lw=1.5)
axes[0].set_xlabel("Longitud")
axes[0].set_ylabel("Latitud")
axes[0].set_title("Mapa de Susceptibilidad")
cbar = plt.colorbar(im, ax=axes[0], ticks=[0, 1, 2, 3, 4], shrink=0.8)
cbar.set_ticklabels(labels)

fpr, tpr, _ = roc_curve(y_test, y_prob)
axes[1].plot(fpr, tpr, color="#B2182B", linewidth=2,
             label=f"Random Forest (AUC={auc:.4f})")
axes[1].plot([0, 1], [0, 1], "k--", linewidth=1, label="Aleatorio (AUC=0.5)")
axes[1].fill_between(fpr, tpr, alpha=0.1, color="#B2182B")
axes[1].set_xlabel("Tasa de Falsos Positivos", fontsize=11)
axes[1].set_ylabel("Tasa de Verdaderos Positivos", fontsize=11)
axes[1].set_title("Curva ROC")
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.3)
axes[1].set_xlim(0, 1)
axes[1].set_ylim(0, 1)

plt.tight_layout()
fig1_path = os.path.join(MODEL_DIR, "P17_fig1_susceptibilidad_ROC.png")
plt.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 1: {fig1_path}")

# FIGURA 2: Importancia de variables
fig, ax = plt.subplots(figsize=(10, 6))
y_pos = np.arange(len(var_names_final))
sorted_imp = importances[sorted_idx]
sorted_names = [var_names_final[i] for i in sorted_idx]

ax.barh(y_pos, sorted_imp[::-1], color="#2166AC", edgecolor="gray", linewidth=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(sorted_names[::-1], fontsize=10)
ax.set_xlabel("Importancia (Gini)", fontsize=11)
ax.set_title("Importancia de Variables Predictoras\n"
             f"Random Forest ({N_TREES} arboles)",
             fontsize=12, fontweight="bold")
ax.grid(True, alpha=0.3, axis="x")

plt.tight_layout()
fig2_path = os.path.join(MODEL_DIR, "P17_fig2_importancia_variables.png")
plt.savefig(fig2_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 2: {fig2_path}")

# FIGURA 3: Dashboard
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

ax1 = fig.add_subplot(gs[0, 0:2])
im1 = ax1.imshow(susceptibility_map, cmap="RdYlGn_r", extent=extent,
                  aspect="auto", vmin=0, vmax=1)
plot_polygon_border(ax1, gdf_aoi, color="black", lw=1.5)
ax1.set_title("Probabilidad de Deslizamiento", fontweight="bold")
ax1.set_xlabel("Longitud")
ax1.set_ylabel("Latitud")
plt.colorbar(im1, ax=ax1, label="Probabilidad", shrink=0.8)

ax2 = fig.add_subplot(gs[0, 2])
ax2.plot(fpr, tpr, color="#B2182B", linewidth=2, label=f"AUC={auc:.4f}")
ax2.plot([0, 1], [0, 1], "k--", linewidth=1)
ax2.fill_between(fpr, tpr, alpha=0.1, color="#B2182B")
ax2.set_xlabel("FPR")
ax2.set_ylabel("TPR")
ax2.set_title("Curva ROC", fontweight="bold")
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

ax3 = fig.add_subplot(gs[1, 0])
ax3.barh(y_pos, sorted_imp[::-1], color="#2166AC")
ax3.set_yticks(y_pos)
ax3.set_yticklabels(sorted_names[::-1], fontsize=8)
ax3.set_xlabel("Importancia")
ax3.set_title("Variables Predictoras", fontweight="bold")
ax3.grid(True, alpha=0.3, axis="x")

ax4 = fig.add_subplot(gs[1, 1])
im4 = ax4.imshow(classes, cmap=cmap_susc, extent=extent, aspect="auto",
                  vmin=-0.5, vmax=4.5, interpolation="nearest")
plot_polygon_border(ax4, gdf_aoi, color="black", lw=1.5)
ax4.set_title("Clases de Susceptibilidad", fontweight="bold")
ax4.set_xlabel("Longitud")
ax4.set_ylabel("Latitud")
cbar4 = plt.colorbar(im4, ax=ax4, ticks=[0, 1, 2, 3, 4], shrink=0.8)
cbar4.set_ticklabels(labels)

ax5 = fig.add_subplot(gs[1, 2])
ax5.axis("off")
tabla = [
    ["Metrica", "Valor"],
    ["AUC-ROC", f"{auc:.4f}"],
    ["OOB Score", f"{rf.oob_score_:.4f}"],
    ["Arboles", f"{N_TREES}"],
    ["Muestras train", f"{len(X_train)}"],
    ["Muestras test", f"{len(X_test)}"],
    ["Predictores", f"{len(var_names_final)}"],
    ["Deslizamientos", f"{len(slide_values)}"],
    ["Variable top", f"{sorted_names[0]}"],
]
table = ax5.table(cellText=tabla, loc="center", cellLoc="center")
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.5)
table[0, 0].set_facecolor("#D5E8F0")
table[0, 1].set_facecolor("#D5E8F0")
ax5.set_title("Resumen del Modelo", fontweight="bold", pad=20)

fig.suptitle("Modelo de Susceptibilidad a Deslizamientos\n"
             f"Zona de Intag (6 parroquias, {area_km2:.0f} km2), Cotacachi, Imbabura, Ecuador",
             fontsize=15, fontweight="bold", y=0.98)

fig3_path = os.path.join(MODEL_DIR, "P17_fig3_dashboard_susceptibilidad.png")
plt.savefig(fig3_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 3: {fig3_path}")

# ============================================================
# PASO 10: REPORTE
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 10: Reporte")
print(f"{'-' * 70}")

reporte_path = os.path.join(MODEL_DIR, "P17_reporte_susceptibilidad.txt")
with open(reporte_path, "w", encoding="utf-8") as rpt:
    rpt.write("=" * 70 + "\n")
    rpt.write("P17 -- MODELO DE SUSCEPTIBILIDAD A DESLIZAMIENTOS\n")
    rpt.write("=" * 70 + "\n")
    rpt.write(f"Fecha: {datetime.now()}\n\n")
    rpt.write(f"METODO: Random Forest ({N_TREES} arboles)\n")
    rpt.write(f"AOI: Poligono real INTAJ.gpkg ({area_km2:.0f} km2, 6 parroquias)\n")
    rpt.write(f"VALIDACION: {int((1-TEST_SIZE)*100)}/{int(TEST_SIZE*100)} train/test\n\n")
    rpt.write(f"RESULTADOS:\n")
    rpt.write(f"  AUC-ROC: {auc:.4f}\n")
    rpt.write(f"  OOB Score: {rf.oob_score_:.4f}\n\n")
    rpt.write(f"VARIABLES PREDICTORAS (por importancia):\n")
    for i in sorted_idx:
        rpt.write(f"  {var_names_final[i]:<20} {importances[i]:.4f}\n")
    rpt.write(f"\nDATOS:\n")
    rpt.write(f"  Inventario: {os.path.basename(inv_path)}\n")
    rpt.write(f"  Deslizamientos: {len(slide_values)}\n")
    rpt.write(f"  No-deslizamiento: {len(neg_values)}\n")
    rpt.write(f"  Train: {len(X_train)} | Test: {len(X_test)}\n\n")
    rpt.write(f"CLASES DE SUSCEPTIBILIDAD:\n")
    for i, label in enumerate(labels):
        n = np.sum(classes == i)
        pct = n / np.sum(~np.isnan(classes)) * 100 if np.sum(~np.isnan(classes)) > 0 else 0
        rpt.write(f"  {label:<12} {n:>8,} pixeles ({pct:.1f}%)\n")
    rpt.write(f"\nCLASIFICACION:\n{report}\n")
    rpt.write(f"\nMATRIZ DE CONFUSION:\n")
    rpt.write(f"  Pred No: {cm[0,0]:>6} {cm[0,1]:>6}\n")
    rpt.write(f"  Pred Si: {cm[1,0]:>6} {cm[1,1]:>6}\n")
print(f"  Reporte: {reporte_path}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P17 COMPLETADO -- PASO 5.5 DE FASE 5")
print(f"{'=' * 70}")
print(f"  AUC-ROC: {auc:.4f}")
print(f"  OOB Score: {rf.oob_score_:.4f}")
print(f"  Variable mas importante: {sorted_names[0]} ({sorted_imp[0]:.4f})")
print(f"  Deslizamientos usados: {len(slide_values)}")
print(f"\n  Rasters: P17_susceptibilidad_probabilidad.tif")
print(f"           P17_susceptibilidad_clases.tif")
print(f"  Figuras: P17_fig1, P17_fig2, P17_fig3")
print(f"\n  FASE 5 COMPLETADA.")
print(f"  SIGUIENTE: Documentacion de Fase 5 y Fase 6 (Publicacion)")
print(f"{'=' * 70}")