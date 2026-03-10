# -*- coding: utf-8 -*-
# ============================================================
# P05b v2.0 - FUSIÓN ASCENDENTE Path 120 + Path 18
#              PONDERACIÓN POR VARIANZA DE CRAMÉR-RAO
# ============================================================
#
# CORRECCIONES RESPECTO A v1.0:
#   1. Solo 2 mosaicos válidos (misma categoría temporal):
#      - M_ASC_CO_12d: A1 (co-evento 12d) + P18_B2 (co-evento 12d)
#      - M_ASC_CO_24d: A3 (co-evento 24d) + P18_A3 (co-evento 24d)
#   2. Eliminados M_ASC_02 y M_ASC_04 que mezclaban categorías:
#      - M_ASC_02: A2 (co-evento) + P18_B1 (pre-evento) ✗
#      - M_ASC_04: A4 (co-evento) + P18_B3 (post-evento) ✗
#   3. Nomenclatura descriptiva: incluye categoría y baseline
#
# FUNDAMENTO CIENTÍFICO:
#
# La varianza de la fase interferométrica sigue la cota de
# Cramér-Rao (Rodriguez & Martin, 1992; Hanssen, 2001):
#
#   σ²_φ = (1 - γ²) / (2·N·γ²)
#
# donde γ = coherencia [0,1] y N = número equivalente de looks.
#
# HyP3 procesó ambos tracks con los mismos parámetros:
# 20×4 looks → N = 80. Como N es idéntico, se cancela al
# normalizar los pesos. El peso óptimo pixel a pixel es:
#
#   w(i,j) = γ²(i,j) / (1 - γ²(i,j))   [Ec. 1]
#
# que es proporcional a la inversa de la varianza de fase.
# Este peso es el estimador de varianza mínima (MVU) para
# promediar observaciones con ruido heteroscedástico.
#
# ESTRATEGIA POR ZONA:
#   - Zona exclusiva P120: dato P120 directamente
#   - Zona exclusiva P18:  dato P18 directamente
#   - Zona de traslape (~39% del AOI):
#       d_fused = (w_120·d_120 + w_18·d_18) / (w_120 + w_18)
#
# PRODUCTO DE COHERENCIA FUSIONADO:
#   - max(γ_120, γ_18) en traslape como cota conservadora
#     (la combinación ponderada tiene varianza menor que
#      cualquier componente individual)
#
# NOTA SOBRE LOS DISPLACEMENT:
#   Path 120 (Frame 1182) y Path 18 (Frame 1180) son ambos
#   ascendentes con ángulos de incidencia similares (~39°).
#   La fusión directa de vert_disp es exacta (ya descompuesto).
#   La fusión de los_disp es aproximada pero válida dado que
#   Δθ_inc entre ambos tracks es < 2° (Pepe & Calò, 2017).
#
# REFERENCIAS:
#   Rodriguez, E. & Martin, J.M. (1992). Theory and Design
#     of InSAR. IEE Proc. F, 139(2), 147-159.
#   Hanssen, R.F. (2001). Radar Interferometry: Data
#     Interpretation and Error Analysis. Kluwer Academic.
#   Pepe, A. & Calò, F. (2017). A Review of InSAR Multi-Track
#     Approaches. Appl. Sci. 7, 1264.
#   Agram, P.S. & Simons, M. (2015). A noise model for InSAR
#     time series. J. Geophys. Res. 120, 3752-3766.
#
# Entrada: hyp3_clipped/ (recortados al AOI por P05 v2)
# Salida:  hyp3_mosaic/  (mosaicos ASC Path120+Path18)
#          hyp3_analysis/ (diagnósticos CRB)
#
# PREREQUISITOS: P05 v2.0 ejecutado exitosamente (12/12 pares)
# SIGUIENTE:     P06 v2 (análisis de deformación)
# ============================================================
import os
import sys
import glob
import json
import numpy as np
from datetime import datetime

print("=" * 70)
print("P05b v2.0 - FUSIÓN ASCENDENTE CRB (Cramér-Rao Bound)")
print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# ============================================================
# 1. VERIFICACIÓN DE LIBRERÍAS
# ============================================================
print("\n--- VERIFICANDO LIBRERÍAS ---\n")

try:
    import rasterio
    from rasterio.merge import merge as rasterio_merge
    from rasterio.warp import reproject, Resampling, calculate_default_transform
    from rasterio.transform import array_bounds
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    print("  ✓ Todas las librerías disponibles")
except ImportError as e:
    print(f"  [ERROR] Librería faltante: {e}")
    sys.exit(1)

# ============================================================
# 2. CONFIGURACIÓN
# ============================================================
SAVE_DIR     = r"D:\POSGRADOS\INTAG\data\sentinel1"
CLIPPED_DIR  = os.path.join(SAVE_DIR, "hyp3_clipped")
MOSAIC_DIR   = os.path.join(SAVE_DIR, "hyp3_mosaic")
ANALYSIS_DIR = os.path.join(SAVE_DIR, "hyp3_analysis")

for d in [MOSAIC_DIR, ANALYSIS_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================
# 3. EMPAREJAMIENTO TEMPORAL CORREGIDO v2.0
# ============================================================
# REGLA: Solo fusionar pares de la MISMA categoría temporal.
# Evento: 19-dic-2023
#
# Mosaicos válidos (ambos CO-EVENTO):
#   M_ASC_CO_12d: A1 (13dic-25dic, 12d) + P18_B2 (18dic-30dic, 12d)
#   M_ASC_CO_24d: A3 (13dic-06ene, 24d) + P18_A3 (18dic-11ene, 24d)
#
# Pares SIN contraparte compatible (se analizan individualmente):
#   A2 (co-evento 24d) → no tiene par P18 co-evento 24d compatible
#   A4 (co-evento 48d) → no tiene par P18 co-evento 48d
#   P18_B1 (PRE-EVENTO) → único par pre-evento, sin par P120
#   P18_B3 (POST-EVENTO) → único par post-evento, sin par P120
#
# Mosaicos ELIMINADOS (mezclaban categorías):
#   M_ASC_02: A2 (co) + P18_B1 (pre) ✗ señal contaminada
#   M_ASC_04: A4 (co) + P18_B3 (post) ✗ señal contaminada

MOSAIC_PAIRS = {
    "M_ASC_CO_12d": {
        "p120": "A1", "p18": "P18_B2",
        "categoria": "CO-EVENTO",
        "baseline": "12d",
        "label": "Co-evento ~12d (ventana más estrecha al evento)",
        "ventana_120": "13dic→25dic (12d)",
        "ventana_18":  "18dic→30dic (12d)",
        "justificacion": "Ambos pares co-evento con baseline 12d. "
                         "Capturan el desplazamiento más cercano al evento.",
    },
    "M_ASC_CO_24d": {
        "p120": "A3", "p18": "P18_A3",
        "categoria": "CO-EVENTO",
        "baseline": "24d",
        "label": "Co-evento ~24d (ventana extendida)",
        "ventana_120": "13dic→06ene (24d)",
        "ventana_18":  "18dic→11ene (24d)",
        "justificacion": "Ambos pares co-evento con baseline 24d. "
                         "Mayor acumulación de desplazamiento post-evento.",
    },
}

# Productos a fusionar (todos excepto coherencia que tiene merge propio)
PRODUCTOS_FUSION = {
    "los_disp":      "_los_disp.tif",
    "vert_disp":     "_vert_disp.tif",
    "unw_phase":     "_unw_phase.tif",
    "wrapped_phase": "_wrapped_phase.tif",
    "amp":           "_amp.tif",
    "dem":           "_dem.tif",
    "inc_map":       "_inc_map.tif",
    "lv_theta":      "_lv_theta.tif",
    "lv_phi":        "_lv_phi.tif",
}

SUFIJO_CORR = "_corr.tif"

# Parámetros CRB
GAMMA_MIN = 0.05   # Coherencia mínima para incluir en fusión
GAMMA_CLIP = 0.999  # Clamp superior para evitar w → ∞

print(f"\n  Directorio recortados: {CLIPPED_DIR}")
print(f"  Directorio mosaicos:  {MOSAIC_DIR}")
print(f"  Mosaicos a generar:   {len(MOSAIC_PAIRS)}")
print(f"  γ mínima:             {GAMMA_MIN}")
print(f"  γ clamp:              {GAMMA_CLIP}")

for mos_id, info in MOSAIC_PAIRS.items():
    print(f"\n  {mos_id}: {info['p120']} + {info['p18']}")
    print(f"    Categoría: {info['categoria']} | Baseline: {info['baseline']}")
    print(f"    P120: {info['ventana_120']}")
    print(f"    P18:  {info['ventana_18']}")

# ============================================================
# 4. FUNCIONES DE FUSIÓN CIENTÍFICA
# ============================================================

def peso_cramer_rao(gamma, gamma_min=GAMMA_MIN, gamma_clip=GAMMA_CLIP):
    # Calcula peso pixel a pixel desde coherencia usando la
    # cota de Cramér-Rao (Rodriguez & Martin, 1992; Hanssen, 2001).
    #
    # w = γ² / (1 - γ²)   proporcional a 1/σ²_φ
    #
    # Se aplica clamp para estabilidad numérica:
    #   - γ < gamma_min → w = 0 (datos no confiables)
    #   - γ > gamma_clip → γ = gamma_clip (evitar w → ∞)
    w = np.zeros_like(gamma, dtype=np.float64)
    valid = (gamma >= gamma_min) & np.isfinite(gamma)
    g = np.clip(gamma[valid], gamma_min, gamma_clip)
    w[valid] = (g ** 2) / (1.0 - g ** 2)
    return w


def alinear_a_grilla_comun(src_path, dst_shape, dst_transform, dst_crs):
    # Reproyecta un raster a una grilla común definida por
    # shape, transform y CRS destino.
    # Necesario porque Path 120 y Path 18 pueden estar en
    # distintas zonas UTM (17N vs 18N) con grillas no coincidentes.
    with rasterio.open(src_path) as src:
        data = np.full((1, dst_shape[0], dst_shape[1]),
                       np.nan, dtype=np.float64)
        reproject(
            source=rasterio.band(src, 1),
            destination=data[0],
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear,
            src_nodata=src.nodata,
            dst_nodata=np.nan,
        )
        if src.nodata is not None:
            data[data == src.nodata] = np.nan
    return data[0]


def calcular_grilla_comun(path_a, path_b):
    # Calcula la grilla unión (extensión combinada) de dos rasters.
    # Retorna shape, transform, crs para la grilla destino.
    #
    # Si los CRS son diferentes (ej. UTM 17N vs 18N, frecuente
    # en zonas limítrofes como Intag), se reproyecta todo al
    # CRS del primer raster (path_a = P120).
    # Justificación: reproyectar a uno de los dos UTM es más
    # preciso que usar WGS84 geográfico, y las distorsiones
    # entre zonas UTM adyacentes son < 0.04% a esta latitud.
    with rasterio.open(path_a) as a, rasterio.open(path_b) as b:
        dst_crs = a.crs

        if a.crs != b.crs:
            print(f"      [INFO] CRS diferentes: {a.crs} vs {b.crs}")
            print(f"      [INFO] Reproyectando a {dst_crs}")
            from rasterio.warp import transform_bounds
            bounds_b_reproj = transform_bounds(
                b.crs, dst_crs,
                b.bounds.left, b.bounds.bottom,
                b.bounds.right, b.bounds.top
            )
            bounds_b_left, bounds_b_bottom, bounds_b_right, bounds_b_top = bounds_b_reproj
        else:
            bounds_b_left   = b.bounds.left
            bounds_b_bottom = b.bounds.bottom
            bounds_b_right  = b.bounds.right
            bounds_b_top    = b.bounds.top

        # Resolución: usar la del primer raster (ambos HyP3
        # con 20x4 looks → resolución casi idéntica ~80m)
        res_x = abs(a.res[0])
        res_y = abs(a.res[1])

        # Extensión unión
        left   = min(a.bounds.left,   bounds_b_left)
        bottom = min(a.bounds.bottom, bounds_b_bottom)
        right  = max(a.bounds.right,  bounds_b_right)
        top    = max(a.bounds.top,    bounds_b_top)

        # Calcular transform y shape
        width  = int(np.ceil((right - left) / res_x))
        height = int(np.ceil((top - bottom) / res_y))
        transform = rasterio.transform.from_bounds(
            left, bottom, right, top, width, height
        )
        return (height, width), transform, dst_crs


def fusionar_crb(raster_120, raster_18, corr_120, corr_18,
                  dst_path):
    # Fusión pixel a pixel con pesos Cramér-Rao.
    #
    # d_fused = (w_120·d_120 + w_18·d_18) / (w_120 + w_18)
    #
    # Zonas exclusivas: dato directo sin interpolación.
    # Retorna: (True, métricas) o (False, None)

    try:
        # 1. Calcular grilla común
        shape, transform, crs = calcular_grilla_comun(raster_120, raster_18)

        # 2. Alinear todos los rasters a la grilla común
        d120 = alinear_a_grilla_comun(raster_120, shape, transform, crs)
        d18  = alinear_a_grilla_comun(raster_18,  shape, transform, crs)
        g120 = alinear_a_grilla_comun(corr_120,   shape, transform, crs)
        g18  = alinear_a_grilla_comun(corr_18,    shape, transform, crs)

        # 3. Calcular pesos CRB  [Ec. 1]
        w120 = peso_cramer_rao(g120)
        w18  = peso_cramer_rao(g18)

        # 4. Máscaras de validez
        valid_120 = np.isfinite(d120) & (w120 > 0)
        valid_18  = np.isfinite(d18)  & (w18 > 0)
        overlap   = valid_120 & valid_18
        solo_120  = valid_120 & ~valid_18
        solo_18   = ~valid_120 & valid_18

        # 5. Fusión
        result = np.full(shape, np.nan, dtype=np.float64)

        # Zona exclusiva P120: dato directo
        result[solo_120] = d120[solo_120]

        # Zona exclusiva P18: dato directo
        result[solo_18] = d18[solo_18]

        # Zona de traslape: promedio ponderado CRB
        w_sum = w120[overlap] + w18[overlap]
        result[overlap] = (
            (w120[overlap] * d120[overlap] + w18[overlap] * d18[overlap])
            / w_sum
        )

        # 6. Coherencia fusionada (max en traslape)
        corr_fused = np.full(shape, np.nan, dtype=np.float64)
        corr_fused[solo_120] = g120[solo_120]
        corr_fused[solo_18]  = g18[solo_18]
        corr_fused[overlap]  = np.maximum(g120[overlap], g18[overlap])

        # 7. Mapa de contribución relativa de P120
        contrib_120 = np.full(shape, np.nan, dtype=np.float64)
        contrib_120[solo_120] = 1.0
        contrib_120[solo_18]  = 0.0
        contrib_120[overlap]  = w120[overlap] / w_sum

        # 8. Guardar resultado
        nodata_out = 0.0
        result_out = np.where(np.isfinite(result), result, nodata_out)

        out_meta = {
            "driver": "GTiff",
            "dtype": "float32",
            "height": shape[0],
            "width": shape[1],
            "count": 1,
            "crs": crs,
            "transform": transform,
            "compress": "lzw",
            "nodata": nodata_out,
        }

        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        with rasterio.open(dst_path, "w", **out_meta) as dst:
            dst.write(result_out.astype(np.float32), 1)

        # 9. Métricas
        n_total  = np.sum(np.isfinite(result))
        n_120    = np.sum(solo_120)
        n_18     = np.sum(solo_18)
        n_ov     = np.sum(overlap)
        pct_120  = 100.0 * n_120 / max(n_total, 1)
        pct_18   = 100.0 * n_18  / max(n_total, 1)
        pct_ov   = 100.0 * n_ov  / max(n_total, 1)
        w120_mean = np.nanmean(contrib_120[overlap]) if n_ov > 0 else 0

        metricas = {
            "n_total": int(n_total),
            "n_solo_120": int(n_120), "pct_solo_120": pct_120,
            "n_solo_18": int(n_18),   "pct_solo_18": pct_18,
            "n_overlap": int(n_ov),   "pct_overlap": pct_ov,
            "w120_mean_overlap": w120_mean,
            "w18_mean_overlap": 1.0 - w120_mean,
            "corr_fused": corr_fused,
            "contrib_120": contrib_120,
            "result": result,
        }

        return True, metricas

    except Exception as e:
        print(f"      ✗ Error en fusión CRB: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def fusionar_coherencia(corr_120, corr_18, dst_path):
    # Para el producto de coherencia, la fusión usa max() en
    # traslape. Esto es conservador: la coherencia efectiva de
    # la combinación ponderada es >= max(componentes).
    try:
        shape, transform, crs = calcular_grilla_comun(corr_120, corr_18)
        g120 = alinear_a_grilla_comun(corr_120, shape, transform, crs)
        g18  = alinear_a_grilla_comun(corr_18,  shape, transform, crs)

        valid_120 = np.isfinite(g120) & (g120 > 0)
        valid_18  = np.isfinite(g18)  & (g18 > 0)
        overlap   = valid_120 & valid_18
        solo_120  = valid_120 & ~valid_18
        solo_18   = ~valid_120 & valid_18

        result = np.full(shape, np.nan, dtype=np.float64)
        result[solo_120] = g120[solo_120]
        result[solo_18]  = g18[solo_18]
        result[overlap]  = np.maximum(g120[overlap], g18[overlap])

        nodata_out = 0.0
        result_out = np.where(np.isfinite(result), result, nodata_out)

        out_meta = {
            "driver": "GTiff", "dtype": "float32",
            "height": shape[0], "width": shape[1], "count": 1,
            "crs": crs, "transform": transform,
            "compress": "lzw", "nodata": nodata_out,
        }
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        with rasterio.open(dst_path, "w", **out_meta) as dst:
            dst.write(result_out.astype(np.float32), 1)
        return True
    except Exception as e:
        print(f"      ✗ Error en fusión coherencia: {e}")
        return False


def generar_mapa_diagnostico(mos_id, metricas, dir_out, label=""):
    # Genera figura diagnóstica con 3 paneles:
    # (a) Mapa de contribución relativa P120 vs P18
    # (b) Coherencia fusionada
    # (c) Resultado fusionado (desplazamiento)
    try:
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        # (a) Contribución P120
        contrib = metricas["contrib_120"]
        im0 = axes[0].imshow(contrib, cmap="RdYlBu", vmin=0, vmax=1)
        axes[0].set_title(f"{mos_id}: Contribución P120\n"
                          f"(1.0=solo P120, 0.0=solo P18)")
        plt.colorbar(im0, ax=axes[0], shrink=0.8, label="Peso P120")

        # (b) Coherencia fusionada
        corr = metricas["corr_fused"]
        im1 = axes[1].imshow(corr, cmap="magma", vmin=0, vmax=1)
        axes[1].set_title("Coherencia fusionada")
        plt.colorbar(im1, ax=axes[1], shrink=0.8, label="Coherencia")

        # (c) Resultado (displacement)
        data = metricas["result"]
        valid = data[np.isfinite(data)]
        if len(valid) > 0:
            vmax = np.percentile(np.abs(valid), 98)
            norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            im2 = axes[2].imshow(data, cmap="RdBu_r", norm=norm)
        else:
            im2 = axes[2].imshow(data, cmap="RdBu_r")
        axes[2].set_title("Desplazamiento fusionado (m)")
        plt.colorbar(im2, ax=axes[2], shrink=0.8, label="m")

        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])

        plt.suptitle(
            f"Fusión CRB: {mos_id} — {label}\n"
            f"Solo P120: {metricas['pct_solo_120']:.1f}% | "
            f"Traslape: {metricas['pct_overlap']:.1f}% "
            f"(w120={metricas['w120_mean_overlap']:.2f}) | "
            f"Solo P18: {metricas['pct_solo_18']:.1f}%",
            fontsize=11
        )
        plt.tight_layout()

        fig_path = os.path.join(dir_out, f"{mos_id}_diagnostico_CRB.png")
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        return fig_path
    except Exception as e:
        print(f"      [AVISO] No se pudo generar figura: {e}")
        return None

# ============================================================
# 5. VERIFICAR DATOS DE ENTRADA
# ============================================================
print(f"\n--- VERIFICANDO DATOS RECORTADOS ---\n")

if not os.path.isdir(CLIPPED_DIR):
    print(f"  [ERROR] No existe: {CLIPPED_DIR}")
    print(f"  Ejecute P05 v2 primero.")
    sys.exit(1)

# Verificar que existen los 4 pares ascendentes necesarios
pares_requeridos = set()
for mos_info in MOSAIC_PAIRS.values():
    pares_requeridos.add(mos_info["p120"])
    pares_requeridos.add(mos_info["p18"])

disponibles = set()
for d in os.listdir(CLIPPED_DIR):
    if os.path.isdir(os.path.join(CLIPPED_DIR, d)):
        disponibles.add(d)

faltantes = pares_requeridos - disponibles
if faltantes:
    print(f"  [ERROR] Faltan pares recortados: {faltantes}")
    print(f"  Disponibles: {sorted(disponibles)}")
    sys.exit(1)

print(f"  ✓ {len(pares_requeridos)} pares requeridos presentes")
for p in sorted(pares_requeridos):
    n_files = len(glob.glob(os.path.join(CLIPPED_DIR, p, "*.tif")))
    print(f"    {p}: {n_files} archivos .tif")

# ============================================================
# 6. EJECUTAR FUSIÓN POR PAR TEMPORAL
# ============================================================
print(f"\n{'=' * 70}")
print("FUSIÓN ASCENDENTE - Path 120 + Path 18")
print(f"{'=' * 70}")
print(f"\n  Método: Ponderación por varianza de Cramér-Rao (Ec. 1)")
print(f"  Peso pixel: w(i,j) = γ² / (1 - γ²)")
print(f"  Ref: Rodriguez & Martin (1992), Hanssen (2001)")
print(f"  γ mínima: {GAMMA_MIN} | γ clamp: {GAMMA_CLIP}")
print(f"\n  {len(MOSAIC_PAIRS)} mosaicos a generar\n")

resumen_fusion = {}

for mos_id, mos_info in MOSAIC_PAIRS.items():
    p120_id = mos_info["p120"]
    p18_id  = mos_info["p18"]

    print(f"\n  {'─' * 60}")
    print(f"  {mos_id}: {p120_id} (P120) + {p18_id} (P18)")
    print(f"  Categoría: {mos_info['categoria']} | Baseline: {mos_info['baseline']}")
    print(f"  P120: {mos_info['ventana_120']}")
    print(f"  P18:  {mos_info['ventana_18']}")
    print(f"  Justificación: {mos_info['justificacion']}")
    print(f"  {'─' * 60}")

    dir_p120 = os.path.join(CLIPPED_DIR, p120_id)
    dir_p18  = os.path.join(CLIPPED_DIR, p18_id)
    dir_mos  = os.path.join(MOSAIC_DIR, mos_id)

    # Localizar coherencia de cada track
    corr_120 = glob.glob(os.path.join(dir_p120, f"*{SUFIJO_CORR}"))
    corr_18  = glob.glob(os.path.join(dir_p18, f"*{SUFIJO_CORR}"))

    if not corr_120 or not corr_18:
        print(f"    [ERROR] No se encontró coherencia para este par")
        continue

    corr_120_path = corr_120[0]
    corr_18_path  = corr_18[0]

    productos_mos = {}
    metricas_principal = None

    # ── Fusionar coherencia (max en traslape) ──
    corr_dst = os.path.join(dir_mos, f"{mos_id}_corr.tif")
    if not os.path.exists(corr_dst):
        ok = fusionar_coherencia(corr_120_path, corr_18_path, corr_dst)
        if ok:
            productos_mos["corr"] = corr_dst
            print(f"    ✓ corr: max() en traslape")
    else:
        productos_mos["corr"] = corr_dst
        print(f"    ✓ corr: ya existe")

    # ── Fusionar cada producto con pesos CRB ──
    for clave, sufijo in PRODUCTOS_FUSION.items():
        f120 = glob.glob(os.path.join(dir_p120, f"*{sufijo}"))
        f18  = glob.glob(os.path.join(dir_p18, f"*{sufijo}"))

        if not f120 or not f18:
            continue

        dst_path = os.path.join(dir_mos, f"{mos_id}{sufijo}")

        if os.path.exists(dst_path):
            productos_mos[clave] = dst_path
            print(f"    ✓ {clave}: ya existe")
            continue

        ok, metricas = fusionar_crb(
            f120[0], f18[0], corr_120_path, corr_18_path, dst_path
        )

        if ok:
            productos_mos[clave] = dst_path

            # Guardar métricas del vert_disp para diagnóstico
            if clave == "vert_disp":
                metricas_principal = metricas

            print(f"    ✓ {clave}: CRB fusión OK"
                  f" | traslape={metricas['pct_overlap']:.1f}%"
                  f" (w120={metricas['w120_mean_overlap']:.2f},"
                  f" w18={metricas['w18_mean_overlap']:.2f})")
        else:
            print(f"    ✗ {clave}: FALLO")

    # ── Diagnóstico visual ──
    if metricas_principal is not None:
        fig_path = generar_mapa_diagnostico(
            mos_id, metricas_principal, ANALYSIS_DIR,
            label=mos_info["label"]
        )
        if fig_path:
            print(f"    ✓ Diagnóstico: {os.path.basename(fig_path)}")

    resumen_fusion[mos_id] = {
        "p120": p120_id, "p18": p18_id,
        "categoria": mos_info["categoria"],
        "baseline": mos_info["baseline"],
        "label": mos_info["label"],
        "productos": {k: v for k, v in productos_mos.items()},
        "n_productos": len(productos_mos),
    }

# ============================================================
# 7. ESTADÍSTICAS DE MOSAICOS
# ============================================================
print(f"\n{'=' * 70}")
print("ESTADÍSTICAS DE MOSAICOS")
print(f"{'=' * 70}\n")

estadisticas_mosaicos = []

for mos_id, info in resumen_fusion.items():
    stats = {"mosaic_id": mos_id, "label": info["label"]}

    # Coherencia
    corr_mos = info["productos"].get("corr")
    if corr_mos and os.path.exists(corr_mos):
        with rasterio.open(corr_mos) as src:
            data = src.read(1).astype(float)
            nd = src.nodata
            if nd is not None:
                data[data == nd] = np.nan
            data[(data < 0) | (data > 1)] = np.nan
            valid = data[~np.isnan(data)]
            if len(valid) > 0:
                stats["corr_mean"] = float(np.mean(valid))
                stats["corr_median"] = float(np.median(valid))
                stats["corr_pct_above_03"] = float(np.sum(valid > 0.3) / len(valid) * 100)

    # LOS (m → cm)
    los_mos = info["productos"].get("los_disp")
    if los_mos and os.path.exists(los_mos):
        with rasterio.open(los_mos) as src:
            data = src.read(1).astype(float) * 100  # m → cm
            nd = src.nodata
            if nd is not None:
                data[data == nd * 100] = np.nan
            valid = data[~np.isnan(data)]
            if len(valid) > 0:
                stats["los_cm_mean"] = float(np.mean(valid))
                stats["los_cm_p5"] = float(np.percentile(valid, 5))
                stats["los_cm_p95"] = float(np.percentile(valid, 95))

    # Vertical (m → cm)
    vert_mos = info["productos"].get("vert_disp")
    if vert_mos and os.path.exists(vert_mos):
        with rasterio.open(vert_mos) as src:
            data = src.read(1).astype(float) * 100
            nd = src.nodata
            if nd is not None:
                data[data == nd * 100] = np.nan
            valid = data[~np.isnan(data)]
            if len(valid) > 0:
                stats["vert_cm_mean"] = float(np.mean(valid))
                stats["vert_cm_p5"] = float(np.percentile(valid, 5))
                stats["vert_cm_p95"] = float(np.percentile(valid, 95))

    estadisticas_mosaicos.append(stats)

    print(f"  {mos_id}: {info['label']}")
    if "corr_mean" in stats:
        print(f"    Coh: μ={stats['corr_mean']:.3f} "
              f"med={stats['corr_median']:.3f} "
              f">0.3={stats.get('corr_pct_above_03', 0):.0f}%")
    if "los_cm_mean" in stats:
        print(f"    LOS: μ={stats['los_cm_mean']:.2f}cm "
              f"[{stats['los_cm_p5']:.2f}, {stats['los_cm_p95']:.2f}]cm")
    if "vert_cm_mean" in stats:
        print(f"    Vert: μ={stats['vert_cm_mean']:.2f}cm "
              f"[{stats['vert_cm_p5']:.2f}, {stats['vert_cm_p95']:.2f}]cm")

# ============================================================
# 8. GUARDAR METADATOS
# ============================================================
print(f"\n--- GUARDANDO METADATOS ---")

output = {
    "generado": datetime.now().isoformat(),
    "script": "P05b_v2.0",
    "version": "2.0 — Solo mosaicos de misma categoría temporal",
    "metodo": "Cramér-Rao Bound (Rodriguez & Martin, 1992)",
    "ecuacion": "w(i,j) = gamma^2 / (1 - gamma^2)",
    "gamma_min": GAMMA_MIN,
    "gamma_clip": GAMMA_CLIP,
    "mosaicos": {k: {kk: vv for kk, vv in v.items() if kk != "productos"}
                 for k, v in resumen_fusion.items()},
    "estadisticas": estadisticas_mosaicos,
    "pares_sin_fusion": {
        "A2": "CO-EVENTO 24d, sin contraparte P18 compatible",
        "A4": "CO-EVENTO 48d, sin contraparte P18 compatible",
        "P18_B1": "PRE-EVENTO, único par pre, sin contraparte P120",
        "P18_B3": "POST-EVENTO, único par post, sin contraparte P120",
    },
    "mosaicos_eliminados": {
        "M_ASC_02": "A2(co) + P18_B1(pre) — categorías incompatibles",
        "M_ASC_04": "A4(co) + P18_B3(post) — categorías incompatibles",
    },
}

meta_path = os.path.join(ANALYSIS_DIR, "fusion_crb_v2.json")
with open(meta_path, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"  ✓ Metadatos: {meta_path}")

# ============================================================
# 9. RESUMEN FINAL
# ============================================================
print(f"\n{'=' * 70}")
print("RESUMEN P05b v2.0")
print(f"{'=' * 70}")

for mos_id, info in resumen_fusion.items():
    print(f"\n  {mos_id}: {info['p120']}+{info['p18']} | "
          f"{info['categoria']} {info['baseline']}")
    print(f"    {info['label']}")
    print(f"    Productos fusionados: {info['n_productos']}")
    for clave in sorted(info["productos"].keys()):
        print(f"      ✓ {clave}")

print(f"\n  PARES ANALIZADOS INDIVIDUALMENTE (sin fusión):")
print(f"    A2     → CO-EVENTO 24d (solo P120, 69.2% AOI)")
print(f"    A4     → CO-EVENTO 48d (solo P120, 69.2% AOI)")
print(f"    P18_B1 → PRE-EVENTO 24d (solo P18, 69.8% AOI)")
print(f"    P18_B3 → POST-EVENTO 24d (solo P18, 69.8% AOI)")

print(f"\n  VERIFICACIÓN CIENTÍFICA:")
print(f"    Método: Cramér-Rao Bound (MVU estimator)")
print(f"    Ecuación: w = γ² / (1 - γ²)")
print(f"    Derivación: σ²_φ = (1-γ²)/(2Nγ²), Rodriguez & Martin 1992")
print(f"    N = 80 (20×4 looks HyP3), idéntico → se cancela")
print(f"    Coherencia fusionada: max(γ_120, γ_18) en traslape")
print(f"    Limitación: Δθ_inc < 2° entre tracks → fusión LOS aproximada")

print(f"\n  ARCHIVOS DE SALIDA:")
print(f"    Mosaicos:    {MOSAIC_DIR}")
print(f"    Diagnósticos: {ANALYSIS_DIR}")
print(f"    Metadatos:    {meta_path}")

print(f"\n  SIGUIENTE PASO:")
print(f"    Ejecutar P06 v2 → Análisis de deformación")
print(f"    Productos disponibles para P06:")
print(f"      Descending: D1, D2, D3, D4 (100% AOI)")
print(f"      Mosaicos ASC: M_ASC_CO_12d, M_ASC_CO_24d (~100% AOI)")
print(f"      Individuales: A2, A4 (69.2%), P18_B1, P18_B3 (69.8%)")
print(f"{'=' * 70}")