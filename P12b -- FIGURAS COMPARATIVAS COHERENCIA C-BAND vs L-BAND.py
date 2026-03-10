# -*- coding: utf-8 -*-
# =============================================================================
# P12b -- FIGURAS COMPARATIVAS COHERENCIA C-BAND vs L-BAND
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 4 NISAR
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.0
#
# PROPOSITO:
#   Generar figuras de calidad publicable comparando la coherencia
#   interferometrica de banda C (Sentinel-1) vs banda L (NISAR)
#   sobre la zona de Intag, Cotacachi, Ecuador.
#
# RESULTADOS CLAVE (P12):
#   C-band media: 0.1812 | L-band media: 0.3738 | Mejora: +106%
#
# REFERENCIA:
#   - Wei & Sandwell (2010): Decorrelation L vs C, IEEE TGRS
#   - Zebker & Villasenor (1992): Decorrelacion, IEEE TGRS 30(5)
# =============================================================================

import os
import sys
import glob
from datetime import datetime

try:
    import h5py
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LinearSegmentedColormap
except ImportError as e:
    print(f"[ERROR] Dependencia faltante: {e}")
    sys.exit(1)

# ============================================================
# CONFIGURACION
# ============================================================
GUNW_DIR = r"D:\POSGRADOS\INTAG\data\nisar\gunw"
RESULTS_DIR = r"D:\POSGRADOS\INTAG\data\nisar\resultados"
SENTINEL1_COH_FILE = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3\mintpy\avgSpatialCoh.h5"
os.makedirs(RESULTS_DIR, exist_ok=True)

# AOI de Intag para recorte
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718

print("=" * 70)
print("P12b -- FIGURAS COMPARATIVAS C-BAND vs L-BAND")
print("=" * 70)
print(f"Fecha: {datetime.now()}")

# ============================================================
# PASO 1: CARGAR COHERENCIA C-BAND (SENTINEL-1)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Carga de coherencia C-band")
print(f"{'-' * 70}")

f = h5py.File(SENTINEL1_COH_FILE, "r")
cband_coh = f["coherence"][:]
# Coordenadas del raster C-band
cb_lat_first = float(f.attrs.get("Y_FIRST", 0.555718))
cb_lat_step = float(f.attrs.get("Y_STEP", -0.000721))
cb_lon_first = float(f.attrs.get("X_FIRST", -79.285032))
cb_lon_step = float(f.attrs.get("X_STEP", 0.000721))
f.close()

cb_ny, cb_nx = cband_coh.shape
cb_lats = cb_lat_first + np.arange(cb_ny) * cb_lat_step
cb_lons = cb_lon_first + np.arange(cb_nx) * cb_lon_step
cb_extent = [cb_lons[0], cb_lons[-1], cb_lats[-1], cb_lats[0]]

cband_valid = cband_coh[(cband_coh >= 0) & (cband_coh <= 1) & (~np.isnan(cband_coh))]
print(f"  Shape: {cband_coh.shape}")
print(f"  Media: {np.mean(cband_valid):.4f}, Mediana: {np.median(cband_valid):.4f}")

# ============================================================
# PASO 2: CARGAR COHERENCIA L-BAND (NISAR GUNW)
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Carga de coherencia L-band")
print(f"{'-' * 70}")

gunw_files = sorted(glob.glob(os.path.join(GUNW_DIR, "*.h5")))
print(f"  Archivos GUNW: {len(gunw_files)}")

# Ruta de coherencia en GUNW
COH_PATH = "/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/coherenceMagnitude"

# Cargar todas las coherencias L-band y calcular promedio
lband_arrays = []
lband_extents = []
lband_info = []

for i, gf in enumerate(gunw_files):
    f = h5py.File(gf, "r")
    try:
        coh = f[COH_PATH][:]
        # Coordenadas
        grp = f["/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH"]
        x_coords = grp["xCoordinates"][:]
        y_coords = grp["yCoordinates"][:]

        # Verificar si son UTM o lat/lon
        proj = int(grp["projection"][()])

        lband_arrays.append(coh)
        lband_extents.append({
            "x_min": x_coords[0], "x_max": x_coords[-1],
            "y_min": y_coords[-1], "y_max": y_coords[0],
            "proj": proj, "x": x_coords, "y": y_coords
        })

        valid = coh[(coh >= 0) & (coh <= 1) & (~np.isnan(coh))]
        info = {
            "archivo": os.path.basename(gf)[:55],
            "media": np.mean(valid),
            "mediana": np.median(valid),
            "shape": coh.shape,
            "proj": proj,
        }
        lband_info.append(info)
        print(f"  [{i+1}] {info['archivo'][:45]} coh={info['media']:.4f} EPSG:{proj}")
    except Exception as e:
        print(f"  [{i+1}] Error: {e}")
    f.close()

# Seleccionar el GUNW con mayor coherencia para la figura de mapa
best_idx = np.argmax([info["media"] for info in lband_info])
print(f"\n  Mejor GUNW para figuras: [{best_idx+1}] coh={lband_info[best_idx]['media']:.4f}")

# ============================================================
# PASO 3: FIGURA 1 -- HISTOGRAMAS COMPARATIVOS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Figura 1 -- Histogramas comparativos")
print(f"{'-' * 70}")

# Recopilar todos los valores L-band validos
all_lband_valid = []
for arr in lband_arrays:
    v = arr[(arr >= 0) & (arr <= 1) & (~np.isnan(arr))]
    all_lband_valid.append(v)
all_lband = np.concatenate(all_lband_valid)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle("Comparacion de Coherencia Interferometrica: Banda C vs Banda L\n"
             "Zona de Intag, Cotacachi, Ecuador", fontsize=13, fontweight="bold")

# Histograma C-band
bins = np.arange(0, 1.02, 0.02)
axes[0].hist(cband_valid, bins=bins, color="#2166AC", alpha=0.8,
             edgecolor="white", linewidth=0.3, density=True)
axes[0].axvline(np.mean(cband_valid), color="red", linestyle="--", linewidth=1.5,
                label=f"Media = {np.mean(cband_valid):.3f}")
axes[0].axvline(np.median(cband_valid), color="orange", linestyle="-.", linewidth=1.5,
                label=f"Mediana = {np.median(cband_valid):.3f}")
axes[0].set_xlabel("Coherencia", fontsize=11)
axes[0].set_ylabel("Densidad", fontsize=11)
axes[0].set_title(f"Sentinel-1 Banda C ($\\lambda$=5.6 cm)\n1,116 interferogramas, 2017-2024",
                   fontsize=10)
axes[0].set_xlim(0, 1)
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# Histograma L-band
axes[1].hist(all_lband, bins=bins, color="#B2182B", alpha=0.8,
             edgecolor="white", linewidth=0.3, density=True)
axes[1].axvline(np.mean(all_lband), color="red", linestyle="--", linewidth=1.5,
                label=f"Media = {np.mean(all_lband):.3f}")
axes[1].axvline(np.median(all_lband), color="orange", linestyle="-.", linewidth=1.5,
                label=f"Mediana = {np.median(all_lband):.3f}")
axes[1].set_xlabel("Coherencia", fontsize=11)
axes[1].set_ylabel("Densidad", fontsize=11)
axes[1].set_title(f"NISAR Banda L ($\\lambda$=24 cm)\n{len(gunw_files)} interferogramas, 2025-2026",
                   fontsize=10)
axes[1].set_xlim(0, 1)
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
fig1_path = os.path.join(RESULTS_DIR, "P12b_fig1_histogramas_CvsL.png")
plt.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Guardado: {fig1_path}")

# ============================================================
# PASO 4: FIGURA 2 -- HISTOGRAMAS SUPERPUESTOS
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Figura 2 -- Histogramas superpuestos")
print(f"{'-' * 70}")

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(cband_valid, bins=bins, color="#2166AC", alpha=0.5,
        edgecolor="none", density=True, label=f"Banda C (media={np.mean(cband_valid):.3f})")
ax.hist(all_lband, bins=bins, color="#B2182B", alpha=0.5,
        edgecolor="none", density=True, label=f"Banda L (media={np.mean(all_lband):.3f})")

ax.axvline(np.mean(cband_valid), color="#2166AC", linestyle="--", linewidth=2)
ax.axvline(np.mean(all_lband), color="#B2182B", linestyle="--", linewidth=2)

# Zona de coherencia util para SBAS (>0.3)
ax.axvspan(0.3, 1.0, alpha=0.08, color="green")
ax.text(0.65, ax.get_ylim()[1]*0.9, "Zona viable\npara SBAS",
        ha="center", fontsize=9, color="green", fontstyle="italic")

ax.set_xlabel("Coherencia Interferometrica", fontsize=12)
ax.set_ylabel("Densidad de Probabilidad", fontsize=12)
ax.set_title("Coherencia C-band vs L-band en Bosque Tropical Humedo\n"
             "Zona de Intag, Cotacachi, Imbabura, Ecuador",
             fontsize=13, fontweight="bold")
ax.set_xlim(0, 1)
ax.legend(fontsize=11, loc="upper right")
ax.grid(True, alpha=0.3)

# Texto de mejora
mejora = (np.mean(all_lband) / np.mean(cband_valid) - 1) * 100
ax.text(0.02, 0.95, f"Mejora L/C: +{mejora:.0f}%",
        transform=ax.transAxes, fontsize=12, fontweight="bold",
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

plt.tight_layout()
fig2_path = os.path.join(RESULTS_DIR, "P12b_fig2_histogramas_superpuestos.png")
plt.savefig(fig2_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Guardado: {fig2_path}")

# ============================================================
# PASO 5: FIGURA 3 -- COHERENCIA POR INTERFEROGRAMA
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Figura 3 -- Coherencia por interferograma")
print(f"{'-' * 70}")

fig, ax = plt.subplots(figsize=(12, 6))

# L-band por interferograma
nombres_l = [info["archivo"][-30:-15] for info in lband_info]
medias_l = [info["media"] for info in lband_info]
medianas_l = [info["mediana"] for info in lband_info]

x_pos = np.arange(len(lband_info))
width = 0.35
bars1 = ax.bar(x_pos - width/2, medias_l, width, color="#B2182B",
               alpha=0.8, label="Media L-band")
bars2 = ax.bar(x_pos + width/2, medianas_l, width, color="#EF8A62",
               alpha=0.8, label="Mediana L-band")

# Linea de referencia C-band
ax.axhline(np.mean(cband_valid), color="#2166AC", linestyle="--",
           linewidth=2, label=f"Media C-band ({np.mean(cband_valid):.3f})")
ax.axhline(np.median(cband_valid), color="#67A9CF", linestyle="-.",
           linewidth=2, label=f"Mediana C-band ({np.median(cband_valid):.3f})")

ax.set_xlabel("Interferograma NISAR", fontsize=11)
ax.set_ylabel("Coherencia", fontsize=11)
ax.set_title("Coherencia por Interferograma NISAR (L-band)\nvs Promedio Sentinel-1 (C-band)",
             fontsize=12, fontweight="bold")
ax.set_xticks(x_pos)
ax.set_xticklabels([f"GUNW\n{i+1}" for i in range(len(lband_info))],
                    fontsize=8)
ax.set_ylim(0, 0.6)
ax.legend(fontsize=9, loc="upper right")
ax.grid(True, alpha=0.3, axis="y")

plt.tight_layout()
fig3_path = os.path.join(RESULTS_DIR, "P12b_fig3_coherencia_por_ifg.png")
plt.savefig(fig3_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Guardado: {fig3_path}")

# ============================================================
# PASO 6: FIGURA 4 -- MAPA DE COHERENCIA C vs L
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Figura 4 -- Mapas de coherencia")
print(f"{'-' * 70}")

# Usar el mejor GUNW
best_coh = lband_arrays[best_idx]
best_ext = lband_extents[best_idx]

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle("Mapas de Coherencia Interferometrica: Banda C vs Banda L\n"
             "Zona de Intag, Cotacachi, Ecuador",
             fontsize=13, fontweight="bold")

# Colormap de coherencia
coh_cmap = "inferno"

# C-band
im1 = axes[0].imshow(cband_coh, cmap=coh_cmap, vmin=0, vmax=0.8,
                      extent=cb_extent, aspect="auto")
axes[0].set_xlabel("Longitud", fontsize=10)
axes[0].set_ylabel("Latitud", fontsize=10)
axes[0].set_title(f"Sentinel-1 Banda C ($\\lambda$=5.6 cm)\n"
                   f"Media={np.mean(cband_valid):.3f}",
                   fontsize=11)
plt.colorbar(im1, ax=axes[0], label="Coherencia", shrink=0.8)

# L-band
# Mascara valores fuera de rango
best_coh_plot = np.where((best_coh >= 0) & (best_coh <= 1), best_coh, np.nan)
lb_extent = [best_ext["x_min"], best_ext["x_max"],
             best_ext["y_min"], best_ext["y_max"]]

im2 = axes[1].imshow(best_coh_plot, cmap=coh_cmap, vmin=0, vmax=0.8,
                      extent=lb_extent, aspect="auto")
axes[1].set_xlabel("Este (m)" if best_ext["proj"] != 4326 else "Longitud", fontsize=10)
axes[1].set_ylabel("Norte (m)" if best_ext["proj"] != 4326 else "Latitud", fontsize=10)
lb_valid = best_coh[(best_coh >= 0) & (best_coh <= 1) & (~np.isnan(best_coh))]
axes[1].set_title(f"NISAR Banda L ($\\lambda$=24 cm)\n"
                   f"Media={np.mean(lb_valid):.3f}",
                   fontsize=11)
plt.colorbar(im2, ax=axes[1], label="Coherencia", shrink=0.8)

plt.tight_layout()
fig4_path = os.path.join(RESULTS_DIR, "P12b_fig4_mapas_coherencia.png")
plt.savefig(fig4_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Guardado: {fig4_path}")

# ============================================================
# PASO 7: FIGURA 5 -- DASHBOARD RESUMEN
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 7: Figura 5 -- Dashboard resumen")
print(f"{'-' * 70}")

fig = plt.figure(figsize=(16, 12))
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3)

# Panel 1: Histogramas superpuestos
ax1 = fig.add_subplot(gs[0, 0:2])
ax1.hist(cband_valid, bins=bins, color="#2166AC", alpha=0.5,
         density=True, label=f"C-band (media={np.mean(cband_valid):.3f})")
ax1.hist(all_lband, bins=bins, color="#B2182B", alpha=0.5,
         density=True, label=f"L-band (media={np.mean(all_lband):.3f})")
ax1.axvline(np.mean(cband_valid), color="#2166AC", linestyle="--", linewidth=1.5)
ax1.axvline(np.mean(all_lband), color="#B2182B", linestyle="--", linewidth=1.5)
ax1.set_xlabel("Coherencia")
ax1.set_ylabel("Densidad")
ax1.set_title("Distribucion de Coherencia C vs L", fontweight="bold")
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# Panel 2: Tabla resumen
ax2 = fig.add_subplot(gs[0, 2])
ax2.axis("off")
tabla_data = [
    ["Metrica", "C-band", "L-band"],
    ["Sensor", "Sentinel-1", "NISAR"],
    ["Longitud onda", "5.6 cm", "24 cm"],
    ["Coh. media", f"{np.mean(cband_valid):.4f}", f"{np.mean(all_lband):.4f}"],
    ["Coh. mediana", f"{np.median(cband_valid):.4f}", f"{np.median(all_lband):.4f}"],
    ["N interferogramas", "1,116", f"{len(gunw_files)}"],
    ["Periodo", "2017-2024", "2025-2026"],
    ["Mejora media", "", f"+{mejora:.0f}%"],
    ["SBAS viable?", "NO", "POSIBLE"],
]
table = ax2.table(cellText=tabla_data, loc="center", cellLoc="center")
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.4)
# Formato de cabecera
for j in range(3):
    table[0, j].set_facecolor("#D5E8F0")
    table[0, j].set_text_props(fontweight="bold")
# Fila de mejora
table[7, 2].set_text_props(fontweight="bold", color="darkgreen")
table[8, 1].set_text_props(color="red", fontweight="bold")
table[8, 2].set_text_props(color="darkgreen", fontweight="bold")
ax2.set_title("Resumen Comparativo", fontweight="bold", pad=20)

# Panel 3: Mapa C-band
ax3 = fig.add_subplot(gs[1, 0])
im3 = ax3.imshow(cband_coh, cmap="inferno", vmin=0, vmax=0.7,
                  extent=cb_extent, aspect="auto")
ax3.set_title(f"C-band (media={np.mean(cband_valid):.3f})", fontweight="bold")
ax3.set_xlabel("Longitud")
ax3.set_ylabel("Latitud")
plt.colorbar(im3, ax=ax3, shrink=0.8)

# Panel 4: Mapa L-band
ax4 = fig.add_subplot(gs[1, 1])
im4 = ax4.imshow(best_coh_plot, cmap="inferno", vmin=0, vmax=0.7,
                  extent=lb_extent, aspect="auto")
ax4.set_title(f"L-band (media={np.mean(lb_valid):.3f})", fontweight="bold")
ax4.set_xlabel("Este (m)" if best_ext["proj"] != 4326 else "Longitud")
ax4.set_ylabel("Norte (m)" if best_ext["proj"] != 4326 else "Latitud")
plt.colorbar(im4, ax=ax4, shrink=0.8)

# Panel 5: Percentiles comparativos
ax5 = fig.add_subplot(gs[1, 2])
percentiles = [10, 25, 50, 75, 90, 95, 99]
pct_c = [np.percentile(cband_valid, p) for p in percentiles]
pct_l = [np.percentile(all_lband, p) for p in percentiles]
ax5.plot(percentiles, pct_c, "o-", color="#2166AC", linewidth=2,
         markersize=6, label="C-band")
ax5.plot(percentiles, pct_l, "s-", color="#B2182B", linewidth=2,
         markersize=6, label="L-band")
ax5.axhline(0.3, color="green", linestyle=":", alpha=0.5, label="Umbral SBAS (0.3)")
ax5.set_xlabel("Percentil")
ax5.set_ylabel("Coherencia")
ax5.set_title("Percentiles de Coherencia", fontweight="bold")
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.3)
ax5.set_ylim(0, 1)

fig.suptitle("Analisis Comparativo de Coherencia: Sentinel-1 (C) vs NISAR (L)\n"
             "Bosque Tropical Humedo, Intag, Cotacachi, Ecuador",
             fontsize=14, fontweight="bold", y=0.98)

fig5_path = os.path.join(RESULTS_DIR, "P12b_fig5_dashboard_CvsL.png")
plt.savefig(fig5_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"  Guardado: {fig5_path}")

# ============================================================
# PASO 8: ESTADISTICAS POR PERCENTIL
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 8: Tabla de percentiles")
print(f"{'-' * 70}")

print(f"\n  {'Percentil':>10} {'C-band':>10} {'L-band':>10} {'Mejora':>10}")
print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
for p in percentiles:
    vc = np.percentile(cband_valid, p)
    vl = np.percentile(all_lband, p)
    mej = (vl/vc - 1)*100 if vc > 0 else 0
    print(f"  {p:>10} {vc:>10.4f} {vl:>10.4f} {mej:>9.1f}%")

# Pixeles utiles para SBAS
print(f"\n  Pixeles con coherencia >= umbral:")
print(f"  {'Umbral':>10} {'C-band':>15} {'L-band':>15}")
print(f"  {'-'*10} {'-'*15} {'-'*15}")
for u in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
    pc = np.sum(cband_valid >= u) / len(cband_valid) * 100
    pl = np.sum(all_lband >= u) / len(all_lband) * 100
    print(f"  {u:>10.1f} {pc:>14.1f}% {pl:>14.1f}%")

# ============================================================
# RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("P12b COMPLETADO")
print(f"{'=' * 70}")
print(f"  Figuras generadas:")
print(f"    1. {fig1_path}")
print(f"    2. {fig2_path}")
print(f"    3. {fig3_path}")
print(f"    4. {fig4_path}")
print(f"    5. {fig5_path}")
print(f"\n  RESULTADO PRINCIPAL:")
print(f"    Banda L (NISAR) mejora coherencia en +{mejora:.0f}%")
print(f"    C-band: {np.mean(cband_valid):.4f} -> L-band: {np.mean(all_lband):.4f}")
print(f"    Esto confirma la viabilidad de monitoreo InSAR en")
print(f"    bosque tropical humedo con banda L.")
print(f"\n  SIGUIENTE PASO: Documentar Fase 4 completa")
print(f"{'=' * 70}")