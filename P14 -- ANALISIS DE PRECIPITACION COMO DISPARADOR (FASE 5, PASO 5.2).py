# -*- coding: utf-8 -*-
# =============================================================================
# P14 -- ANALISIS DE PRECIPITACION COMO DISPARADOR (FASE 5, PASO 5.2)
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 5 Integracion
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.0
#
# PROPOSITO:
#   Analizar la precipitacion como disparador de los deslizamientos del
#   19 de diciembre de 2023 en Intag. Se extraen datos de CHIRPS y
#   ERA5-Land para calcular precipitacion acumulada en ventanas de
#   1, 3, 7, 14, 30 dias previos al evento.
#
# FUENTES DE DATOS:
#   1. CHIRPS Daily (5 km): UCSB-CHG/CHIRPS/DAILY
#      - Precipitacion diaria, 1981-presente
#      - Referencia: Funk et al. (2015), Scientific Data 2:150066
#
#   2. ERA5-Land Daily (9 km): ECMWF/ERA5_LAND/DAILY_AGGR
#      - Precipitacion diaria agregada, 1950-presente
#      - Referencia: Munoz-Sabater et al. (2021), ESSD 13:4349-4383
#
# REFERENCIA:
#   - Notti et al. (2023): Semi-automatic landslide mapping, NHESS
#   - Guzzetti et al. (2008): Rainfall thresholds, Landslides 5:3-17
#   - Protocolo Maestro InSAR Intag v2.0, Paso 5.2
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

# Inicializar GEE
ee.Initialize()

# ============================================================
# CONFIGURACION
# ============================================================
FASE5_DIR = r"D:\POSGRADOS\INTAG\data\fase5_integracion"
PRECIP_DIR = os.path.join(FASE5_DIR, "precipitacion")
os.makedirs(PRECIP_DIR, exist_ok=True)

# AOI de Intag
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718

# Fecha del evento
EVENTO = datetime(2023, 12, 19)
EVENTO_STR = "2023-12-19"

# AOI como geometria GEE
aoi = ee.Geometry.Rectangle([AOI_WEST, AOI_SOUTH, AOI_EAST, AOI_NORTH])

# Ventanas de acumulacion (dias antes del evento)
VENTANAS = [1, 3, 7, 14, 30]

# Periodo de contexto (3 meses antes y 1 mes despues)
INICIO_CONTEXTO = "2023-09-01"
FIN_CONTEXTO = "2024-01-31"

print("=" * 70)
print("P14 -- ANALISIS DE PRECIPITACION COMO DISPARADOR")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"Evento: {EVENTO_STR}")
print(f"AOI: Intag, Cotacachi")

# ============================================================
# PASO 1: CHIRPS - SERIE TEMPORAL DIARIA
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Extraccion de CHIRPS (5 km, diario)")
print(f"{'-' * 70}")

chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY") \
    .filterDate(INICIO_CONTEXTO, FIN_CONTEXTO) \
    .filterBounds(aoi) \
    .select("precipitation")

# Calcular precipitacion media diaria sobre el AOI
def get_chirps_daily(image):
    date = image.date().format("YYYY-MM-dd")
    mean_val = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi,
        scale=5000
    ).get("precipitation")
    return ee.Feature(None, {"date": date, "precip_mm": mean_val})

chirps_fc = chirps.map(get_chirps_daily)
chirps_list = chirps_fc.getInfo()["features"]

chirps_dates = []
chirps_precip = []
for feat in chirps_list:
    props = feat["properties"]
    if props["precip_mm"] is not None:
        chirps_dates.append(datetime.strptime(props["date"], "%Y-%m-%d"))
        chirps_precip.append(props["precip_mm"])

chirps_dates = np.array(chirps_dates)
chirps_precip = np.array(chirps_precip)

print(f"  Dias con datos: {len(chirps_dates)}")
print(f"  Periodo: {chirps_dates[0].strftime('%Y-%m-%d')} a {chirps_dates[-1].strftime('%Y-%m-%d')}")
print(f"  Precipitacion total periodo: {np.sum(chirps_precip):.1f} mm")
print(f"  Precipitacion media diaria: {np.mean(chirps_precip):.1f} mm/dia")

# Precipitacion del dia del evento
idx_evento = np.argmin(np.abs(chirps_dates - EVENTO))
print(f"  Precipitacion dia del evento ({EVENTO_STR}): {chirps_precip[idx_evento]:.1f} mm")

# ============================================================
# PASO 2: ERA5-Land - SERIE TEMPORAL DIARIA
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Extraccion de ERA5-Land (9 km, diario)")
print(f"{'-' * 70}")

era5 = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR") \
    .filterDate(INICIO_CONTEXTO, FIN_CONTEXTO) \
    .filterBounds(aoi) \
    .select("total_precipitation_sum")

def get_era5_daily(image):
    date = image.date().format("YYYY-MM-dd")
    mean_val = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi,
        scale=9000
    ).get("total_precipitation_sum")
    return ee.Feature(None, {"date": date, "precip_m": mean_val})

era5_fc = era5.map(get_era5_daily)
era5_list = era5_fc.getInfo()["features"]

era5_dates = []
era5_precip = []
for feat in era5_list:
    props = feat["properties"]
    if props["precip_m"] is not None:
        era5_dates.append(datetime.strptime(props["date"], "%Y-%m-%d"))
        era5_precip.append(props["precip_m"] * 1000)  # metros a mm

era5_dates = np.array(era5_dates)
era5_precip = np.array(era5_precip)

print(f"  Dias con datos: {len(era5_dates)}")
if len(era5_dates) > 0:
    print(f"  Periodo: {era5_dates[0].strftime('%Y-%m-%d')} a {era5_dates[-1].strftime('%Y-%m-%d')}")
    print(f"  Precipitacion total periodo: {np.sum(era5_precip):.1f} mm")
    print(f"  Precipitacion media diaria: {np.mean(era5_precip):.1f} mm/dia")
    idx_evento_era5 = np.argmin(np.abs(era5_dates - EVENTO))
    print(f"  Precipitacion dia del evento: {era5_precip[idx_evento_era5]:.1f} mm")
else:
    print(f"  ADVERTENCIA: Sin datos ERA5-Land para el periodo.")

# ============================================================
# PASO 3: PRECIPITACION ACUMULADA PRE-EVENTO
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Precipitacion acumulada pre-evento")
print(f"{'-' * 70}")

print(f"\n  CHIRPS - Precipitacion acumulada antes del {EVENTO_STR}:")
print(f"  {'Ventana':>10} {'Acumulada (mm)':>15} {'Media diaria':>15}")
print(f"  {'-'*10} {'-'*15} {'-'*15}")

chirps_acum = {}
for v in VENTANAS:
    inicio = EVENTO - timedelta(days=v)
    mask = (chirps_dates >= inicio) & (chirps_dates < EVENTO)
    acum = np.sum(chirps_precip[mask])
    media = acum / v if v > 0 else 0
    chirps_acum[v] = acum
    print(f"  {v:>7} dias {acum:>15.1f} {media:>15.1f}")

if len(era5_dates) > 0:
    print(f"\n  ERA5-Land - Precipitacion acumulada antes del {EVENTO_STR}:")
    print(f"  {'Ventana':>10} {'Acumulada (mm)':>15} {'Media diaria':>15}")
    print(f"  {'-'*10} {'-'*15} {'-'*15}")

    era5_acum = {}
    for v in VENTANAS:
        inicio = EVENTO - timedelta(days=v)
        mask = (era5_dates >= inicio) & (era5_dates < EVENTO)
        acum = np.sum(era5_precip[mask])
        media = acum / v if v > 0 else 0
        era5_acum[v] = acum
        print(f"  {v:>7} dias {acum:>15.1f} {media:>15.1f}")

# ============================================================
# PASO 4: CONTEXTO HISTORICO (percentiles)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Contexto historico - CHIRPS 2017-2023")
print(f"{'-' * 70}")

# Descargar datos historicos diciembre de cada anio
print("  Descargando datos historicos de diciembre (2017-2022)...")
hist_dic_precip = []
for year in range(2017, 2023):
    chirps_hist = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY") \
        .filterDate(f"{year}-12-01", f"{year}-12-31") \
        .filterBounds(aoi) \
        .select("precipitation")

    total = chirps_hist.sum().reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi,
        scale=5000
    ).get("precipitation").getInfo()

    if total is not None:
        hist_dic_precip.append({"year": year, "total_mm": total})
        print(f"    Diciembre {year}: {total:.1f} mm")

# Diciembre 2023
chirps_dic2023 = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY") \
    .filterDate("2023-12-01", "2023-12-31") \
    .filterBounds(aoi) \
    .select("precipitation")

total_dic2023 = chirps_dic2023.sum().reduceRegion(
    reducer=ee.Reducer.mean(),
    geometry=aoi,
    scale=5000
).get("precipitation").getInfo()

print(f"    Diciembre 2023 (evento): {total_dic2023:.1f} mm")

hist_vals = [h["total_mm"] for h in hist_dic_precip]
if len(hist_vals) > 0:
    media_hist = np.mean(hist_vals)
    std_hist = np.std(hist_vals)
    percentil = sum(1 for v in hist_vals if v < total_dic2023) / len(hist_vals) * 100
    anomalia = (total_dic2023 - media_hist) / std_hist if std_hist > 0 else 0

    print(f"\n  Media historica diciembre: {media_hist:.1f} mm")
    print(f"  Desviacion estandar: {std_hist:.1f} mm")
    print(f"  Diciembre 2023: {total_dic2023:.1f} mm")
    print(f"  Anomalia: {anomalia:+.1f} desviaciones estandar")
    print(f"  Percentil: {percentil:.0f}%")

# ============================================================
# PASO 5: MAPA DE PRECIPITACION ACUMULADA 7 DIAS PRE-EVENTO
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Mapa espacial de precipitacion acumulada (7 dias)")
print(f"{'-' * 70}")

inicio_7d = (EVENTO - timedelta(days=7)).strftime("%Y-%m-%d")
chirps_7d = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY") \
    .filterDate(inicio_7d, EVENTO_STR) \
    .filterBounds(aoi) \
    .select("precipitation") \
    .sum()

# Extraer como array
precip_7d_arr = chirps_7d.sampleRectangle(
    region=aoi,
    defaultValue=0
).get("precipitation").getInfo()

precip_7d = np.array(precip_7d_arr)
print(f"  Shape mapa: {precip_7d.shape}")
print(f"  Precipitacion 7d: min={np.min(precip_7d):.1f}, max={np.max(precip_7d):.1f}, "
      f"media={np.mean(precip_7d):.1f} mm")

# ============================================================
# PASO 6: FIGURAS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Generacion de figuras")
print(f"{'-' * 70}")

# FIGURA 1: Serie temporal diaria CHIRPS
fig, axes = plt.subplots(3, 1, figsize=(14, 14))
fig.suptitle("Analisis de Precipitacion - Deslizamientos de Intag\n"
             "Evento: 19 de diciembre de 2023",
             fontsize=14, fontweight="bold")

# Panel 1: Serie temporal diaria
ax1 = axes[0]
ax1.bar(chirps_dates, chirps_precip, width=1, color="#2166AC", alpha=0.7,
        label="CHIRPS diario")
if len(era5_dates) > 0:
    ax1.plot(era5_dates, era5_precip, "r-", linewidth=0.8, alpha=0.6,
             label="ERA5-Land diario")
ax1.axvline(EVENTO, color="red", linestyle="--", linewidth=2,
            label=f"Evento ({EVENTO_STR})")
ax1.set_ylabel("Precipitacion (mm/dia)", fontsize=11)
ax1.set_title("Precipitacion Diaria", fontweight="bold")
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax1.xaxis.set_major_locator(mdates.MonthLocator())

# Panel 2: Precipitacion acumulada movil (7 dias)
ax2 = axes[1]
# Calcular acumulada movil 7 dias
acum_7d = np.convolve(chirps_precip, np.ones(7), mode="same")
ax2.fill_between(chirps_dates, 0, acum_7d, color="#2166AC", alpha=0.3)
ax2.plot(chirps_dates, acum_7d, color="#2166AC", linewidth=1.5,
         label="Acumulada 7 dias (CHIRPS)")
ax2.axvline(EVENTO, color="red", linestyle="--", linewidth=2,
            label=f"Evento ({EVENTO_STR})")

# Marcar el valor del evento
idx_ev = np.argmin(np.abs(chirps_dates - EVENTO))
ax2.annotate(f"{acum_7d[idx_ev]:.0f} mm",
             xy=(EVENTO, acum_7d[idx_ev]),
             xytext=(EVENTO + timedelta(days=5), acum_7d[idx_ev] + 20),
             arrowprops=dict(arrowstyle="->", color="red"),
             fontsize=11, fontweight="bold", color="red")

ax2.set_ylabel("Precipitacion acumulada 7d (mm)", fontsize=11)
ax2.set_title("Precipitacion Acumulada Movil (7 dias)", fontweight="bold")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax2.xaxis.set_major_locator(mdates.MonthLocator())

# Panel 3: Barras de acumulacion pre-evento
ax3 = axes[2]
x_pos = np.arange(len(VENTANAS))
width = 0.35
bars_c = ax3.bar(x_pos - width/2, [chirps_acum[v] for v in VENTANAS],
                 width, color="#2166AC", alpha=0.8, label="CHIRPS")
if len(era5_dates) > 0:
    bars_e = ax3.bar(x_pos + width/2, [era5_acum[v] for v in VENTANAS],
                     width, color="#B2182B", alpha=0.8, label="ERA5-Land")

# Etiquetas en las barras
for bar in bars_c:
    h = bar.get_height()
    ax3.text(bar.get_x() + bar.get_width()/2., h + 1,
             f"{h:.0f}", ha="center", va="bottom", fontsize=9)

ax3.set_xlabel("Ventana temporal pre-evento", fontsize=11)
ax3.set_ylabel("Precipitacion acumulada (mm)", fontsize=11)
ax3.set_title("Precipitacion Acumulada Pre-Evento por Ventana Temporal",
              fontweight="bold")
ax3.set_xticks(x_pos)
ax3.set_xticklabels([f"{v} dias" for v in VENTANAS])
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3, axis="y")

plt.tight_layout()
fig1_path = os.path.join(PRECIP_DIR, "P14_fig1_precipitacion_temporal.png")
plt.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 1: {fig1_path}")

# FIGURA 2: Mapa espacial + contexto historico
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Precipitacion Espacial y Contexto Historico - Intag",
             fontsize=13, fontweight="bold")

# Panel 1: Mapa 7 dias
extent = [AOI_WEST, AOI_EAST, AOI_SOUTH, AOI_NORTH]
im = axes[0].imshow(precip_7d, cmap="Blues", extent=extent, aspect="auto",
                     vmin=0)
axes[0].set_xlabel("Longitud")
axes[0].set_ylabel("Latitud")
axes[0].set_title(f"Precipitacion acumulada\n{inicio_7d} a {EVENTO_STR} (7 dias)")
plt.colorbar(im, ax=axes[0], label="mm", shrink=0.8)

# Panel 2: Contexto historico diciembre
years = [h["year"] for h in hist_dic_precip] + [2023]
totals = [h["total_mm"] for h in hist_dic_precip] + [total_dic2023]
colors = ["#2166AC"] * len(hist_dic_precip) + ["#B2182B"]
axes[1].bar(years, totals, color=colors)
axes[1].axhline(media_hist, color="gray", linestyle="--",
                label=f"Media historica ({media_hist:.0f} mm)")
axes[1].set_xlabel("Anio")
axes[1].set_ylabel("Precipitacion diciembre (mm)")
axes[1].set_title("Contexto Historico - Diciembre")
axes[1].legend()
axes[1].grid(True, alpha=0.3, axis="y")
# Etiqueta 2023
axes[1].annotate(f"{total_dic2023:.0f} mm", xy=(2023, total_dic2023),
                 xytext=(2023, total_dic2023 + 20),
                 ha="center", fontsize=10, fontweight="bold", color="red")

plt.tight_layout()
fig2_path = os.path.join(PRECIP_DIR, "P14_fig2_mapa_contexto.png")
plt.savefig(fig2_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Figura 2: {fig2_path}")

# ============================================================
# PASO 7: EXPORTAR DATOS COMO CSV
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Exportar datos")
print(f"{'-' * 70}")

# CSV de serie temporal CHIRPS
csv_path = os.path.join(PRECIP_DIR, "P14_chirps_diario.csv")
with open(csv_path, "w") as f:
    f.write("fecha,precipitacion_mm\n")
    for d, p in zip(chirps_dates, chirps_precip):
        f.write(f"{d.strftime('%Y-%m-%d')},{p:.2f}\n")
print(f"  CSV CHIRPS: {csv_path}")

# CSV de ERA5
if len(era5_dates) > 0:
    csv_era5 = os.path.join(PRECIP_DIR, "P14_era5land_diario.csv")
    with open(csv_era5, "w") as f:
        f.write("fecha,precipitacion_mm\n")
        for d, p in zip(era5_dates, era5_precip):
            f.write(f"{d.strftime('%Y-%m-%d')},{p:.2f}\n")
    print(f"  CSV ERA5: {csv_era5}")

# ============================================================
# PASO 8: REPORTE
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 8: Reporte")
print(f"{'-' * 70}")

reporte_path = os.path.join(PRECIP_DIR, "P14_reporte_precipitacion.txt")
with open(reporte_path, "w", encoding="utf-8") as rpt:
    rpt.write("=" * 70 + "\n")
    rpt.write("P14 -- ANALISIS DE PRECIPITACION COMO DISPARADOR\n")
    rpt.write("=" * 70 + "\n")
    rpt.write(f"Fecha: {datetime.now()}\n")
    rpt.write(f"Evento: {EVENTO_STR}\n\n")
    rpt.write("PRECIPITACION ACUMULADA PRE-EVENTO (CHIRPS):\n")
    for v in VENTANAS:
        rpt.write(f"  {v:2d} dias: {chirps_acum[v]:.1f} mm\n")
    rpt.write(f"\nDia del evento: {chirps_precip[idx_evento]:.1f} mm\n")
    rpt.write(f"\nCONTEXTO HISTORICO DICIEMBRE:\n")
    rpt.write(f"  Media 2017-2022: {media_hist:.1f} mm\n")
    rpt.write(f"  Diciembre 2023: {total_dic2023:.1f} mm\n")
    rpt.write(f"  Anomalia: {anomalia:+.1f} sigma\n")
    rpt.write(f"  Percentil: {percentil:.0f}%\n")
print(f"  Reporte: {reporte_path}")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P14 COMPLETADO -- PASO 5.2 DE FASE 5")
print(f"{'=' * 70}")
print(f"  Evento: {EVENTO_STR}")
print(f"  Precipitacion dia del evento (CHIRPS): {chirps_precip[idx_evento]:.1f} mm")
print(f"  Acumulada 7 dias: {chirps_acum[7]:.1f} mm")
print(f"  Acumulada 30 dias: {chirps_acum[30]:.1f} mm")
print(f"  Diciembre 2023: {total_dic2023:.1f} mm")
if len(hist_vals) > 0:
    print(f"  Anomalia: {anomalia:+.1f} sigma vs media historica")
print(f"\n  Figuras: P14_fig1, P14_fig2")
print(f"  CSVs: P14_chirps_diario.csv, P14_era5land_diario.csv")
print(f"\n  SIGUIENTE PASO: P15 (Paso 5.3 - Humedad del suelo)")
print(f"{'=' * 70}")