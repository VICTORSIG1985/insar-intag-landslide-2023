# -*- coding: utf-8 -*-
# =============================================================================
# P15 -- ANALISIS DE HUMEDAD DEL SUELO (FASE 5, PASO 5.3)
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 5 Integracion
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.1 (corregido: SMAP v008 + manejo dias sin datos)
#
# PROPOSITO:
#   Analizar la humedad del suelo como predictor de los deslizamientos
#   del 19 de diciembre de 2023. La hipotesis es que la saturacion del
#   suelo redujo la resistencia al corte, provocando el colapso masivo
#   al combinarse con la precipitacion extrema (+2.3 sigma, P14).
#
# FUENTES DE DATOS:
#   1. NASA SMAP L4 v008 (9 km, 3 horas): NASA/SMAP/SPL4SMGP/008
#      - Humedad del suelo superficial y de zona raiz
#      - Referencia: Reichle et al. (2017), J. Geophys. Res.
#
#   2. ERA5-Land Soil Moisture (9 km, diario): ECMWF/ERA5_LAND/DAILY_AGGR
#      - Humedad volumetrica en 4 capas (0-7, 7-28, 28-100, 100-289 cm)
#      - Referencia: Munoz-Sabater et al. (2021), ESSD
#
# REFERENCIA:
#   - Protocolo Maestro InSAR Intag v2.0, Paso 5.3
#   - Iverson (2000): Landslide triggering by rain infiltration, WRR
#   - Segoni et al. (2018): Soil moisture and landslides, Landslides
# =============================================================================

import os
import sys
import numpy as np
from datetime import datetime, timedelta

try:
    import ee
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    sys.exit(1)

ee.Initialize()

# ============================================================
# CONFIGURACION
# ============================================================
FASE5_DIR = r"D:\POSGRADOS\INTAG\data\fase5_integracion"
SM_DIR = os.path.join(FASE5_DIR, "humedad_suelo")
os.makedirs(SM_DIR, exist_ok=True)

AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718
aoi = ee.Geometry.Rectangle([AOI_WEST, AOI_SOUTH, AOI_EAST, AOI_NORTH])

EVENTO = datetime(2023, 12, 19)
EVENTO_STR = "2023-12-19"

INICIO = "2023-06-01"
FIN = "2024-03-31"

print("=" * 70)
print("P15 -- ANALISIS DE HUMEDAD DEL SUELO")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"Evento: {EVENTO_STR}")
print(f"Periodo: {INICIO} a {FIN}")

# ============================================================
# PASO 1: SMAP L4 v008 - HUMEDAD SUPERFICIAL Y ZONA RAIZ
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: SMAP L4 v008 (9 km, 3 horas)")
print(f"{'-' * 70}")

smap = ee.ImageCollection("NASA/SMAP/SPL4SMGP/008") \
    .filterDate(INICIO, FIN) \
    .filterBounds(aoi)

# Agrupar por dia y promediar (SMAP tiene multiples obs/dia)
days = ee.List.sequence(
    ee.Date(INICIO).millis(),
    ee.Date(FIN).millis(),
    86400000
)

def smap_daily_mean(day_millis):
    day = ee.Date(day_millis)
    day_end = day.advance(1, "day")
    daily = smap.filterDate(day, day_end)
    count = daily.size()
    mean_img = ee.Algorithms.If(
        count.gt(0),
        daily.mean(),
        ee.Image.constant([0, 0]).rename(["sm_surface", "sm_rootzone"])
    )
    mean_img = ee.Image(mean_img)
    date_str = day.format("YYYY-MM-dd")
    surface = ee.Algorithms.If(
        count.gt(0),
        mean_img.select("sm_surface").reduceRegion(
            reducer=ee.Reducer.mean(), geometry=aoi, scale=9000
        ).get("sm_surface"),
        None
    )
    rootzone = ee.Algorithms.If(
        count.gt(0),
        mean_img.select("sm_rootzone").reduceRegion(
            reducer=ee.Reducer.mean(), geometry=aoi, scale=9000
        ).get("sm_rootzone"),
        None
    )
    return ee.Feature(None, {
        "date": date_str,
        "sm_surface": surface,
        "sm_rootzone": rootzone
    })

print("  Extrayendo datos SMAP diarios (puede tomar 2-3 min)...")
smap_fc = ee.FeatureCollection(days.map(smap_daily_mean))
smap_list = smap_fc.getInfo()["features"]

smap_dates = []
smap_surface = []
smap_rootzone = []
for feat in smap_list:
    props = feat["properties"]
    if props["sm_surface"] is not None and props["sm_rootzone"] is not None:
        smap_dates.append(datetime.strptime(props["date"], "%Y-%m-%d"))
        smap_surface.append(props["sm_surface"])
        smap_rootzone.append(props["sm_rootzone"])

smap_dates = np.array(smap_dates)
smap_surface = np.array(smap_surface)
smap_rootzone = np.array(smap_rootzone)

print(f"  Dias con datos: {len(smap_dates)}")
if len(smap_dates) > 0:
    print(f"  Periodo: {smap_dates[0].strftime('%Y-%m-%d')} a {smap_dates[-1].strftime('%Y-%m-%d')}")
    print(f"  Humedad superficial (0-5 cm):")
    print(f"    Media: {np.mean(smap_surface):.4f} m3/m3")
    print(f"    Min: {np.min(smap_surface):.4f}, Max: {np.max(smap_surface):.4f}")

    idx_ev = np.argmin(np.abs(smap_dates - EVENTO))
    print(f"    Dia del evento: {smap_surface[idx_ev]:.4f} m3/m3")

    print(f"  Humedad zona raiz (0-100 cm):")
    print(f"    Media: {np.mean(smap_rootzone):.4f} m3/m3")
    print(f"    Min: {np.min(smap_rootzone):.4f}, Max: {np.max(smap_rootzone):.4f}")
    print(f"    Dia del evento: {smap_rootzone[idx_ev]:.4f} m3/m3")

    pct_surf = np.sum(smap_surface < smap_surface[idx_ev]) / len(smap_surface) * 100
    pct_root = np.sum(smap_rootzone < smap_rootzone[idx_ev]) / len(smap_rootzone) * 100
    print(f"    Percentil superficie dia evento: {pct_surf:.0f}%")
    print(f"    Percentil zona raiz dia evento: {pct_root:.0f}%")
else:
    print("  ADVERTENCIA: Sin datos SMAP para el periodo.")
    idx_ev = 0
    pct_surf = 0
    pct_root = 0

# ============================================================
# PASO 2: ERA5-Land - HUMEDAD EN 4 CAPAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: ERA5-Land humedad del suelo (4 capas)")
print(f"{'-' * 70}")

era5_sm = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR") \
    .filterDate(INICIO, FIN) \
    .filterBounds(aoi)

capas_era5 = [
    ("volumetric_soil_water_layer_1", "0-7 cm"),
    ("volumetric_soil_water_layer_2", "7-28 cm"),
    ("volumetric_soil_water_layer_3", "28-100 cm"),
    ("volumetric_soil_water_layer_4", "100-289 cm"),
]

def get_era5_sm(image):
    date = image.date().format("YYYY-MM-dd")
    props = {"date": date}
    for band, _ in capas_era5:
        val = image.select(band).reduceRegion(
            reducer=ee.Reducer.mean(), geometry=aoi, scale=9000
        ).get(band)
        props[band] = val
    return ee.Feature(None, props)

print("  Extrayendo datos ERA5-Land (puede tomar 1-2 min)...")
era5_sm_fc = era5_sm.map(get_era5_sm)
era5_sm_list = era5_sm_fc.getInfo()["features"]

era5_sm_dates = []
era5_layers = {band: [] for band, _ in capas_era5}

for feat in era5_sm_list:
    props = feat["properties"]
    if props[capas_era5[0][0]] is not None:
        era5_sm_dates.append(datetime.strptime(props["date"], "%Y-%m-%d"))
        for band, _ in capas_era5:
            era5_layers[band].append(props[band])

era5_sm_dates = np.array(era5_sm_dates)
for band in era5_layers:
    era5_layers[band] = np.array(era5_layers[band])

print(f"  Dias con datos: {len(era5_sm_dates)}")
idx_ev_era5 = 0
if len(era5_sm_dates) > 0:
    idx_ev_era5 = np.argmin(np.abs(era5_sm_dates - EVENTO))
    print(f"\n  Humedad del suelo dia del evento ({EVENTO_STR}):")
    for band, label in capas_era5:
        val = era5_layers[band][idx_ev_era5]
        media = np.mean(era5_layers[band])
        pct = np.sum(era5_layers[band] < val) / len(era5_layers[band]) * 100
        print(f"    {label}: {val:.4f} m3/m3 (media={media:.4f}, percentil={pct:.0f}%)")

# ============================================================
# PASO 3: ANALISIS TEMPORAL - EVOLUCION PRE-EVENTO
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Evolucion pre-evento (30 dias)")
print(f"{'-' * 70}")

coef_surf = [0, 0]
coef_root = [0, 0]
if len(smap_dates) > 0:
    inicio_30d = EVENTO - timedelta(days=30)
    mask_pre = (smap_dates >= inicio_30d) & (smap_dates <= EVENTO)

    if np.sum(mask_pre) > 2:
        print(f"  SMAP 30 dias pre-evento:")
        print(f"    Superficie: de {np.min(smap_surface[mask_pre]):.4f} a {np.max(smap_surface[mask_pre]):.4f} m3/m3")
        print(f"    Zona raiz: de {np.min(smap_rootzone[mask_pre]):.4f} a {np.max(smap_rootzone[mask_pre]):.4f} m3/m3")

        dias_pre = np.array([(d - inicio_30d).days for d in smap_dates[mask_pre]])
        coef_surf = np.polyfit(dias_pre, smap_surface[mask_pre], 1)
        coef_root = np.polyfit(dias_pre, smap_rootzone[mask_pre], 1)
        print(f"    Tendencia superficie: {coef_surf[0]*30:+.4f} m3/m3 en 30 dias")
        print(f"    Tendencia zona raiz: {coef_root[0]*30:+.4f} m3/m3 en 30 dias")
    else:
        print(f"  Insuficientes datos SMAP en ventana de 30 dias.")

# ============================================================
# PASO 4: CARGAR PRECIPITACION CHIRPS (P14)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Carga de precipitacion CHIRPS (de P14)")
print(f"{'-' * 70}")

chirps_csv = os.path.join(FASE5_DIR, "precipitacion", "P14_chirps_diario.csv")
chirps_dates_p = []
chirps_precip = []
if os.path.exists(chirps_csv):
    with open(chirps_csv, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split(",")
            chirps_dates_p.append(datetime.strptime(parts[0], "%Y-%m-%d"))
            chirps_precip.append(float(parts[1]))
    chirps_dates_p = np.array(chirps_dates_p)
    chirps_precip = np.array(chirps_precip)
    print(f"  CHIRPS cargado: {len(chirps_dates_p)} dias")
else:
    print(f"  ADVERTENCIA: No se encontro {chirps_csv}")
    print(f"  Ejecutar P14 primero.")

# ============================================================
# PASO 5: FIGURAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Generacion de figuras")
print(f"{'-' * 70}")

# FIGURA 1: Serie temporal completa
fig, axes = plt.subplots(3, 1, figsize=(14, 14), sharex=True)
fig.suptitle("Humedad del Suelo y Precipitacion - Deslizamientos de Intag\n"
             "Evento: 19 de diciembre de 2023",
             fontsize=14, fontweight="bold")

# Panel 1: Precipitacion
if len(chirps_dates_p) > 0:
    axes[0].bar(chirps_dates_p, chirps_precip, width=1, color="#2166AC",
                alpha=0.7, label="CHIRPS diario")
    axes[0].axvline(EVENTO, color="red", linestyle="--", linewidth=2)
    axes[0].set_ylabel("Precipitacion (mm/dia)")
    axes[0].set_title("Precipitacion Diaria (CHIRPS)", fontweight="bold")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)
    axes[0].invert_yaxis()

# Panel 2: SMAP
if len(smap_dates) > 0:
    axes[1].plot(smap_dates, smap_surface, "-", color="#D95F02", linewidth=1.5,
                 label="Superficial (0-5 cm)")
    axes[1].plot(smap_dates, smap_rootzone, "-", color="#1B9E77", linewidth=1.5,
                 label="Zona raiz (0-100 cm)")
    axes[1].axvline(EVENTO, color="red", linestyle="--", linewidth=2,
                    label=f"Evento ({EVENTO_STR})")
    axes[1].set_ylabel("Humedad (m$^3$/m$^3$)")
    axes[1].set_title("Humedad del Suelo SMAP L4 v008", fontweight="bold")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)
else:
    axes[1].text(0.5, 0.5, "Sin datos SMAP disponibles",
                 ha="center", va="center", transform=axes[1].transAxes, fontsize=14)
    axes[1].set_title("SMAP L4 - Sin datos", fontweight="bold")

# Panel 3: ERA5-Land 4 capas
if len(era5_sm_dates) > 0:
    colores = ["#E66101", "#FDB863", "#B2ABD2", "#5E3C99"]
    for (band, label), color in zip(capas_era5, colores):
        axes[2].plot(era5_sm_dates, era5_layers[band], "-", color=color,
                     linewidth=1.2, label=label)
    axes[2].axvline(EVENTO, color="red", linestyle="--", linewidth=2,
                    label="Evento")
    axes[2].set_ylabel("Humedad (m$^3$/m$^3$)")
    axes[2].set_xlabel("Fecha")
    axes[2].set_title("Humedad del Suelo ERA5-Land (4 capas)", fontweight="bold")
    axes[2].legend(fontsize=8, ncol=2)
    axes[2].grid(True, alpha=0.3)
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axes[2].xaxis.set_major_locator(mdates.MonthLocator())

plt.tight_layout()
fig1_path = os.path.join(SM_DIR, "P15_fig1_humedad_temporal.png")
plt.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 1: {fig1_path}")

# FIGURA 2: Zoom 60 dias
fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
fig.suptitle("Humedad del Suelo - Zoom 60 Dias Alrededor del Evento",
             fontsize=13, fontweight="bold")

zoom_ini = EVENTO - timedelta(days=30)
zoom_fin = EVENTO + timedelta(days=30)

if len(chirps_dates_p) > 0:
    mask_z = (chirps_dates_p >= zoom_ini) & (chirps_dates_p <= zoom_fin)
    axes[0].bar(chirps_dates_p[mask_z], chirps_precip[mask_z], width=1,
                color="#2166AC", alpha=0.7)
    axes[0].axvline(EVENTO, color="red", linestyle="--", linewidth=2)
    axes[0].set_ylabel("Precipitacion (mm/dia)")
    axes[0].set_title("Precipitacion Diaria")
    axes[0].grid(True, alpha=0.3)
    axes[0].invert_yaxis()

if len(smap_dates) > 0:
    mask_z = (smap_dates >= zoom_ini) & (smap_dates <= zoom_fin)
    if np.sum(mask_z) > 0:
        axes[1].plot(smap_dates[mask_z], smap_surface[mask_z], "o-", color="#D95F02",
                     linewidth=1.5, markersize=4, label="Superficial (0-5 cm)")
        axes[1].plot(smap_dates[mask_z], smap_rootzone[mask_z], "s-", color="#1B9E77",
                     linewidth=1.5, markersize=4, label="Zona raiz (0-100 cm)")
        axes[1].axvline(EVENTO, color="red", linestyle="--", linewidth=2,
                        label="Evento 19-dic-2023")
        axes[1].annotate(f"{smap_surface[idx_ev]:.3f}",
                         xy=(EVENTO, smap_surface[idx_ev]),
                         xytext=(EVENTO + timedelta(days=3), smap_surface[idx_ev] + 0.01),
                         arrowprops=dict(arrowstyle="->", color="#D95F02"),
                         fontsize=10, fontweight="bold", color="#D95F02")
        axes[1].annotate(f"{smap_rootzone[idx_ev]:.3f}",
                         xy=(EVENTO, smap_rootzone[idx_ev]),
                         xytext=(EVENTO + timedelta(days=3), smap_rootzone[idx_ev] - 0.01),
                         arrowprops=dict(arrowstyle="->", color="#1B9E77"),
                         fontsize=10, fontweight="bold", color="#1B9E77")
        axes[1].set_ylabel("Humedad (m$^3$/m$^3$)")
        axes[1].set_title("Humedad del Suelo SMAP L4")
        axes[1].legend(fontsize=9)
        axes[1].grid(True, alpha=0.3)
elif len(era5_sm_dates) > 0:
    # Si no hay SMAP, usar ERA5 en el zoom
    mask_z = (era5_sm_dates >= zoom_ini) & (era5_sm_dates <= zoom_fin)
    colores = ["#E66101", "#FDB863", "#B2ABD2", "#5E3C99"]
    for (band, label), color in zip(capas_era5, colores):
        axes[1].plot(era5_sm_dates[mask_z], era5_layers[band][mask_z], "o-",
                     color=color, linewidth=1.2, markersize=3, label=label)
    axes[1].axvline(EVENTO, color="red", linestyle="--", linewidth=2)
    axes[1].set_ylabel("Humedad (m$^3$/m$^3$)")
    axes[1].set_title("Humedad ERA5-Land (zoom)")
    axes[1].legend(fontsize=8, ncol=2)
    axes[1].grid(True, alpha=0.3)

axes[1].set_xlabel("Fecha")
axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
axes[1].xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))

plt.tight_layout()
fig2_path = os.path.join(SM_DIR, "P15_fig2_humedad_zoom.png")
plt.savefig(fig2_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 2: {fig2_path}")

# FIGURA 3: Correlacion precipitacion-humedad
r_surf = 0
r_root = 0
# Usar SMAP si disponible, si no ERA5 capa 1
use_smap = len(smap_dates) > 5
use_era5 = len(era5_sm_dates) > 5

if len(chirps_dates_p) > 0 and (use_smap or use_era5):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Correlacion Precipitacion - Humedad del Suelo",
                 fontsize=13, fontweight="bold")

    if use_smap:
        sm_dates_corr = smap_dates
        sm_surf_corr = smap_surface
        sm_root_corr = smap_rootzone
        src_label = "SMAP"
    else:
        sm_dates_corr = era5_sm_dates
        sm_surf_corr = era5_layers[capas_era5[0][0]]
        sm_root_corr = era5_layers[capas_era5[2][0]]
        src_label = "ERA5-Land"

    precip_interp = np.interp(
        [d.toordinal() for d in sm_dates_corr],
        [d.toordinal() for d in chirps_dates_p],
        chirps_precip
    )
    precip_acum7 = np.convolve(precip_interp, np.ones(7), mode="same")

    axes[0].scatter(precip_acum7, sm_surf_corr, c="#D95F02", alpha=0.4, s=15)
    axes[0].set_xlabel("Precipitacion acumulada 7 dias (mm)")
    axes[0].set_ylabel("Humedad superficial (m$^3$/m$^3$)")
    axes[0].set_title(f"Superficie ({src_label})")
    r_surf = np.corrcoef(precip_acum7, sm_surf_corr)[0, 1]
    axes[0].text(0.05, 0.95, f"r = {r_surf:.3f}",
                 transform=axes[0].transAxes, fontsize=12, fontweight="bold",
                 verticalalignment="top")
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(precip_acum7, sm_root_corr, c="#1B9E77", alpha=0.4, s=15)
    axes[1].set_xlabel("Precipitacion acumulada 7 dias (mm)")
    axes[1].set_ylabel("Humedad profunda (m$^3$/m$^3$)")
    axes[1].set_title(f"Profundidad ({src_label})")
    r_root = np.corrcoef(precip_acum7, sm_root_corr)[0, 1]
    axes[1].text(0.05, 0.95, f"r = {r_root:.3f}",
                 transform=axes[1].transAxes, fontsize=12, fontweight="bold",
                 verticalalignment="top")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    fig3_path = os.path.join(SM_DIR, "P15_fig3_correlacion.png")
    plt.savefig(fig3_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Figura 3: {fig3_path}")

# ============================================================
# PASO 6: EXPORTAR DATOS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Exportar datos")
print(f"{'-' * 70}")

if len(smap_dates) > 0:
    csv_smap = os.path.join(SM_DIR, "P15_smap_diario.csv")
    with open(csv_smap, "w") as f:
        f.write("fecha,sm_surface_m3m3,sm_rootzone_m3m3\n")
        for d, s, rz in zip(smap_dates, smap_surface, smap_rootzone):
            f.write(f"{d.strftime('%Y-%m-%d')},{s:.6f},{rz:.6f}\n")
    print(f"  CSV SMAP: {csv_smap}")

if len(era5_sm_dates) > 0:
    csv_era5 = os.path.join(SM_DIR, "P15_era5land_humedad.csv")
    with open(csv_era5, "w") as f:
        header = "fecha," + ",".join([f"swvl_{label.replace(' ','_').replace('-','_')}" for _, label in capas_era5])
        f.write(header + "\n")
        for i, d in enumerate(era5_sm_dates):
            vals = ",".join([f"{era5_layers[band][i]:.6f}" for band, _ in capas_era5])
            f.write(f"{d.strftime('%Y-%m-%d')},{vals}\n")
    print(f"  CSV ERA5: {csv_era5}")

# ============================================================
# PASO 7: REPORTE
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Reporte")
print(f"{'-' * 70}")

reporte_path = os.path.join(SM_DIR, "P15_reporte_humedad.txt")
with open(reporte_path, "w", encoding="utf-8") as rpt:
    rpt.write("=" * 70 + "\n")
    rpt.write("P15 -- ANALISIS DE HUMEDAD DEL SUELO\n")
    rpt.write("=" * 70 + "\n")
    rpt.write(f"Fecha: {datetime.now()}\n")
    rpt.write(f"Evento: {EVENTO_STR}\n\n")

    if len(smap_dates) > 0:
        rpt.write("SMAP L4 v008 - DIA DEL EVENTO:\n")
        rpt.write(f"  Superficie (0-5 cm): {smap_surface[idx_ev]:.4f} m3/m3 "
                  f"(percentil {pct_surf:.0f}%)\n")
        rpt.write(f"  Zona raiz (0-100 cm): {smap_rootzone[idx_ev]:.4f} m3/m3 "
                  f"(percentil {pct_root:.0f}%)\n\n")

    if len(era5_sm_dates) > 0:
        rpt.write("ERA5-Land - DIA DEL EVENTO:\n")
        for band, label in capas_era5:
            val = era5_layers[band][idx_ev_era5]
            pct = np.sum(era5_layers[band] < val) / len(era5_layers[band]) * 100
            rpt.write(f"  {label}: {val:.4f} m3/m3 (percentil {pct:.0f}%)\n")

    rpt.write(f"\nCORRELACION PRECIPITACION-HUMEDAD:\n")
    rpt.write(f"  r(precip_7d, superficie): {r_surf:.3f}\n")
    rpt.write(f"  r(precip_7d, profundidad): {r_root:.3f}\n")

print(f"  Reporte: {reporte_path}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P15 COMPLETADO -- PASO 5.3 DE FASE 5")
print(f"{'=' * 70}")
if len(smap_dates) > 0:
    print(f"  SMAP superficie dia evento: {smap_surface[idx_ev]:.4f} m3/m3 (P{pct_surf:.0f})")
    print(f"  SMAP zona raiz dia evento: {smap_rootzone[idx_ev]:.4f} m3/m3 (P{pct_root:.0f})")
if len(era5_sm_dates) > 0:
    print(f"  ERA5 capa 0-7cm dia evento: {era5_layers[capas_era5[0][0]][idx_ev_era5]:.4f} m3/m3")
    print(f"  ERA5 capa 28-100cm dia evento: {era5_layers[capas_era5[2][0]][idx_ev_era5]:.4f} m3/m3")
print(f"  Correlacion precip-superficie: r={r_surf:.3f}")
print(f"  Correlacion precip-profundidad: r={r_root:.3f}")
print(f"\n  Figuras: P15_fig1, P15_fig2, P15_fig3")
print(f"\n  SIGUIENTE PASO: P16 (Paso 5.4 - Cobertura y uso del suelo)")
print(f"{'=' * 70}")