# =============================================================================
# P09 — DISEÑO Y ENVÍO DEL STACK SBAS A HyP3
# Versión: 1.0
# Fecha: 2026-02-28
# Autor: Víctor Pinto
# Proyecto: Detección de deslizamientos mediante InSAR — Zona de Intag, Cotacachi
# Fase 3: Series Temporales SBAS — Protocolo Maestro v1.0
#
# OBJETIVO:
# Diseñar la red de pares interferométricos SBAS a partir de las 242 escenas
# SLC del Path 40 Descendente verificadas por P08, y enviar los jobs a HyP3
# para procesamiento.
#
# JUSTIFICACIÓN DE PARÁMETROS:
# - Umbral temporal: 48 días (tau ~4 días en bosque tropical C-band;
#   pares >48d tendrán coherencia ~0 sobre vegetación)
# - Umbral perpendicular: 200 m (valor por defecto HyP3, minimiza
#   decorrelación geométrica)
# - Conexiones: mínimo 2 vecinos por imagen para garantizar red conectada
#   (Berardino et al. 2002)
#
# INPUTS:
# - D:\POSGRADOS\INTAG\data\sentinel1\diagnostico_cobertura\fase3\
#   P08_escenas_Path_40_fase3.csv (catálogo de 242 escenas)
# - Credenciales NASA Earthdata (para hyp3_sdk)
#
# OUTPUTS (en D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3\):
# - P09_red_sbas_pares.csv (lista de pares con baselines)
# - P09_red_sbas_network.png (gráfico temporal-perpendicular)
# - P09_red_sbas_conectividad.png (verificación de conectividad)
# - P09_red_sbas_histograma.png (distribución de baselines)
# - P09_reporte_red_sbas.txt (reporte completo)
# - P09_jobs_enviados.json (registro de jobs HyP3, si se envían)
#
# DEPENDENCIAS:
# pip install asf_search hyp3_sdk pandas numpy matplotlib networkx
#
# NOTA SOBRE ENVÍO A HyP3:
# El script tiene dos modos:
# - MODO DISEÑO (por defecto): Solo diseña la red y genera diagnósticos.
#   No requiere credenciales. Usar para verificar antes de enviar.
# - MODO ENVÍO: Envía los pares a HyP3. Requiere credenciales NASA
#   Earthdata. Activar con ENVIAR_A_HYP3 = True.
# =============================================================================

import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from itertools import combinations

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

# Ruta de entrada (salida de P08)
ESCENAS_CSV = r"D:\POSGRADOS\INTAG\data\sentinel1\diagnostico_cobertura\fase3\P08_escenas_Path_40_fase3.csv"

# Ruta de salida
OUTPUT_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Parámetros de la red SBAS
TEMPORAL_BASELINE_MAX = 48      # días (justificado por tau ~4d en bosque tropical)
PERPENDICULAR_BASELINE_MAX = 200  # metros (valor por defecto HyP3)
MIN_CONNECTIONS_PER_IMAGE = 2    # mínimo para red conectada

# Modo de operación
ENVIAR_A_HYP3 = True  # MODO ENVÍO ACTIVO
HYP3_PROJECT_NAME = "INTAG_SBAS_Fase3_Path40"

# Parámetros de procesamiento HyP3 (consistentes con Fase 2)
HYP3_PARAMS = {
    "looks": "20x4",
    "include_dem": True,                    # DEM recortado (consistente con Fase 2)
    "include_inc_map": True,                # Mapa ángulo incidencia (consistente con Fase 2)
    "include_look_vectors": True,
    "include_displacement_maps": True,
    "include_wrapped_phase": True,
    "apply_water_mask": False,              # False = consistente con Fase 2 (ríos en AOI)
    "phase_filter_parameter": 0.6,          # 0.6 = consistente con Fase 2 (P04/P04b)
}

# Evento de referencia
FECHA_EVENTO = pd.Timestamp("2023-12-19")

print("=" * 70)
print("P09 — DISEÑO Y ENVÍO DEL STACK SBAS A HyP3")
print("=" * 70)
print(f"Fecha de ejecución: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Modo: {'ENVÍO a HyP3' if ENVIAR_A_HYP3 else 'DISEÑO (sin envío)'}")
print(f"Umbral temporal: {TEMPORAL_BASELINE_MAX} días")
print(f"Umbral perpendicular: {PERPENDICULAR_BASELINE_MAX} m")
print()

# =============================================================================
# PASO 1: CARGAR CATÁLOGO DE ESCENAS
# =============================================================================

print("-" * 70)
print("PASO 1: Carga del catálogo de escenas Path 40")
print("-" * 70)

df_escenas = pd.read_csv(ESCENAS_CSV)
print(f"Escenas cargadas: {len(df_escenas)}")

# Asegurar que la columna de fecha esté en formato correcto
df_escenas['date'] = pd.to_datetime(df_escenas['date'])
df_escenas = df_escenas.sort_values('date').reset_index(drop=True)

print(f"Periodo: {df_escenas['date'].min().strftime('%Y-%m-%d')} a "
      f"{df_escenas['date'].max().strftime('%Y-%m-%d')}")
print(f"Escenas únicas: {df_escenas['date'].nunique()}")

# Verificar duplicados de fecha
duplicados = df_escenas[df_escenas.duplicated(subset='date', keep=False)]
if len(duplicados) > 0:
    print(f"  ADVERTENCIA: {len(duplicados)} escenas con fechas duplicadas.")
    print(f"  Manteniendo primera ocurrencia por fecha.")
    df_escenas = df_escenas.drop_duplicates(subset='date', keep='first').reset_index(drop=True)
    print(f"  Escenas después de deduplicación: {len(df_escenas)}")

print()

# =============================================================================
# PASO 2: OBTENER BASELINES PERPENDICULARES
# =============================================================================

print("-" * 70)
print("PASO 2: Cálculo de baselines perpendiculares")
print("-" * 70)

# Para obtener baselines perpendiculares precisas, necesitamos consultar ASF
# Intentamos usar asf_search; si no está disponible, estimamos a partir
# de los metadatos orbitales

try:
    import asf_search as asf

    print("Obteniendo baselines perpendiculares desde ASF...")

    # Buscar la escena de referencia (centro temporal del stack)
    mid_idx = len(df_escenas) // 2
    ref_scene_id = df_escenas.iloc[mid_idx]['scene_id']
    ref_date = df_escenas.iloc[mid_idx]['date']
    print(f"  Escena de referencia: {ref_scene_id}")
    print(f"  Fecha de referencia: {ref_date.strftime('%Y-%m-%d')}")

    # Buscar baselines usando la escena de referencia
    # asf_search puede calcular baselines relativas a una referencia
    ref_result = asf.search(granule_list=[ref_scene_id])

    if len(ref_result) > 0:
        ref_granule = ref_result[0]
        # Obtener stack de baselines
        stack = ref_granule.stack()

        # Crear diccionario de baselines perpendiculares por fecha
        bperp_dict = {}
        for item in stack:
            props = item.properties
            scene_date = pd.to_datetime(props.get('startTime', '')).normalize()
            bperp = props.get('perpendicularBaseline', None)
            if bperp is not None:
                bperp_dict[scene_date] = float(bperp)

        # Asignar baselines a nuestras escenas
        df_escenas['bperp'] = df_escenas['date'].map(
            lambda d: bperp_dict.get(d.normalize(), np.nan)
        )

        n_con_bperp = df_escenas['bperp'].notna().sum()
        print(f"  Baselines perpendiculares obtenidas: {n_con_bperp}/{len(df_escenas)}")

        if n_con_bperp < len(df_escenas) * 0.5:
            print("  ADVERTENCIA: <50% de escenas con baseline perpendicular.")
            print("  Usando baseline perpendicular = 0 para escenas sin dato.")
            print("  (El filtrado real por Bperp se hará en HyP3/MintPy)")
            df_escenas['bperp'] = df_escenas['bperp'].fillna(0)
    else:
        print("  No se encontró la escena de referencia en ASF.")
        print("  Usando baseline perpendicular = 0 (estimación).")
        df_escenas['bperp'] = 0

except ImportError:
    print("asf_search no disponible. Usando baseline perpendicular = 0.")
    print("(El filtrado real por Bperp se realizará durante el procesamiento)")
    df_escenas['bperp'] = 0

except Exception as e:
    print(f"Error obteniendo baselines: {e}")
    print("Usando baseline perpendicular = 0 (estimación).")
    df_escenas['bperp'] = 0

print()

# =============================================================================
# PASO 3: GENERAR RED DE PARES SBAS
# =============================================================================

print("-" * 70)
print("PASO 3: Generación de la red SBAS")
print("-" * 70)

dates = df_escenas['date'].values
scene_ids = df_escenas['scene_id'].values
bperps = df_escenas['bperp'].values
n_scenes = len(dates)

print(f"Generando pares con:")
print(f"  Temporal baseline máximo: {TEMPORAL_BASELINE_MAX} días")
print(f"  Perpendicular baseline máximo: {PERPENDICULAR_BASELINE_MAX} m")

# Generar todos los pares que cumplen los umbrales
pares = []
for i in range(n_scenes):
    for j in range(i + 1, n_scenes):
        dt = (pd.Timestamp(dates[j]) - pd.Timestamp(dates[i])).days
        dbperp = abs(bperps[j] - bperps[i])

        if dt <= TEMPORAL_BASELINE_MAX and dbperp <= PERPENDICULAR_BASELINE_MAX:
            # Clasificar el par respecto al evento
            d_i = pd.Timestamp(dates[i])
            d_j = pd.Timestamp(dates[j])

            if d_j < FECHA_EVENTO:
                categoria = "PRE-EVENTO"
            elif d_i > FECHA_EVENTO:
                categoria = "POST-EVENTO"
            else:
                categoria = "CO-EVENTO"

            pares.append({
                'idx_ref': i,
                'idx_sec': j,
                'ref_scene': scene_ids[i],
                'sec_scene': scene_ids[j],
                'ref_date': pd.Timestamp(dates[i]).strftime('%Y-%m-%d'),
                'sec_date': pd.Timestamp(dates[j]).strftime('%Y-%m-%d'),
                'temporal_baseline': dt,
                'perp_baseline': round(dbperp, 1),
                'categoria': categoria
            })

df_pares = pd.DataFrame(pares)
print(f"\nPares generados: {len(df_pares)}")
print(f"  Pre-evento: {(df_pares['categoria'] == 'PRE-EVENTO').sum()}")
print(f"  Co-evento: {(df_pares['categoria'] == 'CO-EVENTO').sum()}")
print(f"  Post-evento: {(df_pares['categoria'] == 'POST-EVENTO').sum()}")

# Estadísticas de baselines
print(f"\nEstadísticas de baseline temporal:")
print(f"  Media: {df_pares['temporal_baseline'].mean():.1f} días")
print(f"  Mediana: {df_pares['temporal_baseline'].median():.0f} días")
print(f"  Mín: {df_pares['temporal_baseline'].min()} días")
print(f"  Máx: {df_pares['temporal_baseline'].max()} días")

# Distribución por intervalo temporal
print(f"\nDistribución por intervalo temporal:")
for dt_val in sorted(df_pares['temporal_baseline'].unique()):
    count = (df_pares['temporal_baseline'] == dt_val).sum()
    if count > 0:
        print(f"  {dt_val:3d} días: {count:4d} pares")

print()

# =============================================================================
# PASO 4: VERIFICAR CONECTIVIDAD DE LA RED
# =============================================================================

print("-" * 70)
print("PASO 4: Verificación de conectividad de la red")
print("-" * 70)

try:
    import networkx as nx

    # Construir grafo
    G = nx.Graph()
    for i in range(n_scenes):
        G.add_node(i, date=pd.Timestamp(dates[i]).strftime('%Y-%m-%d'))

    for _, par in df_pares.iterrows():
        G.add_edge(par['idx_ref'], par['idx_sec'],
                    temporal=par['temporal_baseline'])

    # Verificar conectividad
    is_connected = nx.is_connected(G)
    n_components = nx.number_connected_components(G)

    print(f"Red conectada: {'SÍ' if is_connected else 'NO'}")
    print(f"Componentes conexos: {n_components}")

    if not is_connected:
        components = list(nx.connected_components(G))
        print(f"\n  ADVERTENCIA: Red desconectada. Esto impide estimación")
        print(f"  no sesgada de velocidades (ASF Product Guide).")

        for ci, comp in enumerate(components):
            comp_dates = [pd.Timestamp(dates[idx]) for idx in comp]
            print(f"\n  Componente {ci + 1}: {len(comp)} escenas")
            print(f"    Periodo: {min(comp_dates).strftime('%Y-%m-%d')} a "
                  f"{max(comp_dates).strftime('%Y-%m-%d')}")

        # Intentar conectar añadiendo pares con baseline mayor
        print(f"\n  Intentando conectar con pares adicionales...")
        pares_adicionales = []

        for ci in range(len(components) - 1):
            comp_a = list(components[ci])
            comp_b = list(components[ci + 1])

            # Encontrar el par con menor baseline temporal entre componentes
            best_dt = 999
            best_pair = None

            for ia in comp_a:
                for ib in comp_b:
                    dt = abs((pd.Timestamp(dates[ia]) -
                              pd.Timestamp(dates[ib])).days)
                    if dt < best_dt:
                        best_dt = dt
                        best_pair = (ia, ib, dt)

            if best_pair:
                ia, ib, dt = best_pair
                if ia > ib:
                    ia, ib = ib, ia

                d_i = pd.Timestamp(dates[ia])
                d_j = pd.Timestamp(dates[ib])

                if d_j < FECHA_EVENTO:
                    cat = "PRE-EVENTO"
                elif d_i > FECHA_EVENTO:
                    cat = "POST-EVENTO"
                else:
                    cat = "CO-EVENTO"

                dbperp = abs(bperps[ib] - bperps[ia])

                par_adicional = {
                    'idx_ref': ia,
                    'idx_sec': ib,
                    'ref_scene': scene_ids[ia],
                    'sec_scene': scene_ids[ib],
                    'ref_date': d_i.strftime('%Y-%m-%d'),
                    'sec_date': d_j.strftime('%Y-%m-%d'),
                    'temporal_baseline': dt,
                    'perp_baseline': round(dbperp, 1),
                    'categoria': cat
                }
                pares_adicionales.append(par_adicional)
                print(f"    Par puente: {d_i.strftime('%Y-%m-%d')} - "
                      f"{d_j.strftime('%Y-%m-%d')} ({dt} días)")

        if pares_adicionales:
            df_pares_extra = pd.DataFrame(pares_adicionales)
            df_pares = pd.concat([df_pares, df_pares_extra],
                                 ignore_index=True)
            print(f"\n  {len(pares_adicionales)} pares puente añadidos.")
            print(f"  Total pares actualizado: {len(df_pares)}")

            # Verificar nueva conectividad
            G2 = nx.Graph()
            for i in range(n_scenes):
                G2.add_node(i)
            for _, par in df_pares.iterrows():
                G2.add_edge(par['idx_ref'], par['idx_sec'])

            is_connected = nx.is_connected(G2)
            print(f"  Red conectada tras pares puente: "
                  f"{'SÍ' if is_connected else 'NO'}")

            if not is_connected:
                print("  ADVERTENCIA: La red sigue desconectada.")
                print("  Puede ser necesario añadir pares manualmente en Vertex.")
                print("  MintPy puede manejar subconjuntos desconectados pero")
                print("  con menor precisión en la inversión.")

    # Estadísticas de conexiones por nodo
    degrees = dict(G.degree())
    deg_values = list(degrees.values())
    print(f"\nConexiones por escena:")
    print(f"  Media: {np.mean(deg_values):.1f}")
    print(f"  Mín: {min(deg_values)} (escena más aislada)")
    print(f"  Máx: {max(deg_values)}")

    # Identificar escenas con pocas conexiones
    escenas_debiles = [i for i, d in degrees.items() if d < MIN_CONNECTIONS_PER_IMAGE]
    if escenas_debiles:
        print(f"\n  ADVERTENCIA: {len(escenas_debiles)} escenas con "
              f"<{MIN_CONNECTIONS_PER_IMAGE} conexiones:")
        for idx in escenas_debiles[:10]:
            print(f"    {pd.Timestamp(dates[idx]).strftime('%Y-%m-%d')} "
                  f"({degrees[idx]} conexiones)")

except ImportError:
    print("networkx no disponible. Verificación de conectividad omitida.")
    print("Instalar con: pip install networkx")
    print("RECOMENDACIÓN: Verificar conectividad en ASF Vertex antes de enviar.")
    is_connected = None

print()

# =============================================================================
# PASO 5: GENERAR FIGURAS DE DIAGNÓSTICO
# =============================================================================

print("-" * 70)
print("PASO 5: Generación de figuras de diagnóstico")
print("-" * 70)

# --- Figura 1: Red temporal-perpendicular ---
fig, ax = plt.subplots(1, 1, figsize=(16, 8))

# Dibujar conexiones
for _, par in df_pares.iterrows():
    d1 = pd.Timestamp(par['ref_date'])
    d2 = pd.Timestamp(par['sec_date'])
    b1 = bperps[par['idx_ref']]
    b2 = bperps[par['idx_sec']]

    color = '#CCCCCC'
    alpha = 0.3
    lw = 0.5

    if par['temporal_baseline'] <= 12:
        color = '#2196F3'
        alpha = 0.5
        lw = 0.8
    elif par['temporal_baseline'] <= 24:
        color = '#4CAF50'
        alpha = 0.4
        lw = 0.6

    ax.plot([d1, d2], [b1, b2], color=color, alpha=alpha, linewidth=lw)

# Dibujar escenas como puntos
scene_dates_ts = [pd.Timestamp(d) for d in dates]
colors_scatter = ['red' if d >= FECHA_EVENTO else '#1565C0' for d in scene_dates_ts]
ax.scatter(scene_dates_ts, bperps, c=colors_scatter, s=15, zorder=5, edgecolors='white', linewidth=0.5)

# Línea del evento
ax.axvline(x=FECHA_EVENTO, color='red', linestyle='--', linewidth=2,
           label=f'Evento {FECHA_EVENTO.strftime("%Y-%m-%d")}', zorder=4)

ax.set_title(f"Red SBAS — Path 40 Descendente\n"
             f"{len(df_escenas)} escenas, {len(df_pares)} pares "
             f"(Btemp \u2264 {TEMPORAL_BASELINE_MAX}d, Bperp \u2264 {PERPENDICULAR_BASELINE_MAX}m)",
             fontsize=13, fontweight='bold')
ax.set_xlabel("Fecha de adquisición")
ax.set_ylabel("Baseline perpendicular (m)")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.xaxis.set_major_locator(mdates.YearLocator())

fig_path = os.path.join(OUTPUT_DIR, "P09_red_sbas_network.png")
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Guardado: {fig_path}")

# --- Figura 2: Histograma de baselines ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Temporal
ax = axes[0]
ax.hist(df_pares['temporal_baseline'], bins=range(0, TEMPORAL_BASELINE_MAX + 6, 6),
        color='#2196F3', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axvline(x=12, color='green', linestyle='--', linewidth=1.5,
           label='Revisita nominal (12d)')
ax.set_title("Distribución de baseline temporal", fontsize=11, fontweight='bold')
ax.set_xlabel("Baseline temporal (días)")
ax.set_ylabel("Número de pares")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Perpendicular
ax = axes[1]
if df_pares['perp_baseline'].max() > 0:
    ax.hist(df_pares['perp_baseline'], bins=30,
            color='#FF9800', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.set_title("Distribución de baseline perpendicular", fontsize=11, fontweight='bold')
    ax.set_xlabel("Baseline perpendicular (m)")
    ax.set_ylabel("Número de pares")
else:
    ax.text(0.5, 0.5, "Baselines perpendiculares\nno disponibles\n(estimadas como 0)",
            ha='center', va='center', fontsize=12, transform=ax.transAxes)
    ax.set_title("Baseline perpendicular (no disponible)", fontsize=11)
ax.grid(True, alpha=0.3)

plt.suptitle(f"Estadísticas de baselines — {len(df_pares)} pares SBAS",
             fontsize=13, fontweight='bold')
plt.tight_layout()

fig_path = os.path.join(OUTPUT_DIR, "P09_red_sbas_histograma.png")
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Guardado: {fig_path}")

# --- Figura 3: Conexiones por escena ---
fig, ax = plt.subplots(1, 1, figsize=(16, 5))

if 'degrees' in dir():
    deg_by_date = [(pd.Timestamp(dates[i]), degrees.get(i, 0)) for i in range(n_scenes)]
    deg_by_date.sort(key=lambda x: x[0])

    ax.bar([d[0] for d in deg_by_date], [d[1] for d in deg_by_date],
           width=8, color='#2196F3', alpha=0.7)
    ax.axhline(y=MIN_CONNECTIONS_PER_IMAGE, color='red', linestyle='--',
               linewidth=1.5, label=f'Mínimo requerido ({MIN_CONNECTIONS_PER_IMAGE})')
    ax.axvline(x=FECHA_EVENTO, color='red', linestyle='-', linewidth=2,
               label='Evento', alpha=0.7)

    ax.set_title(f"Conexiones por escena en la red SBAS\n"
                 f"Media: {np.mean(deg_values):.1f} | "
                 f"Mín: {min(deg_values)} | Máx: {max(deg_values)}",
                 fontsize=12, fontweight='bold')
    ax.set_xlabel("Fecha")
    ax.set_ylabel("Número de conexiones")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

fig_path = os.path.join(OUTPUT_DIR, "P09_red_sbas_conectividad.png")
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Guardado: {fig_path}")

print()

# =============================================================================
# PASO 6: EXPORTAR RED DE PARES
# =============================================================================

print("-" * 70)
print("PASO 6: Exportación de la red de pares")
print("-" * 70)

# CSV completo de pares
csv_path = os.path.join(OUTPUT_DIR, "P09_red_sbas_pares.csv")
df_pares.to_csv(csv_path, index=False)
print(f"  Guardado: {csv_path} ({len(df_pares)} pares)")

# Reporte textual
reporte_lines = []
reporte_lines.append("=" * 70)
reporte_lines.append("P09 — REPORTE DE DISEÑO DE RED SBAS")
reporte_lines.append("=" * 70)
reporte_lines.append(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
reporte_lines.append(f"Track: Path 40 Descendente (Frame 590)")
reporte_lines.append(f"Escenas: {len(df_escenas)}")
reporte_lines.append(f"Periodo: {df_escenas['date'].min().strftime('%Y-%m-%d')} a "
                     f"{df_escenas['date'].max().strftime('%Y-%m-%d')}")
reporte_lines.append(f"Evento: {FECHA_EVENTO.strftime('%Y-%m-%d')}")
reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("PARÁMETROS DE LA RED")
reporte_lines.append("-" * 70)
reporte_lines.append(f"  Umbral temporal máximo: {TEMPORAL_BASELINE_MAX} días")
reporte_lines.append(f"  Umbral perpendicular máximo: {PERPENDICULAR_BASELINE_MAX} m")
reporte_lines.append(f"  Conexiones mínimas por imagen: {MIN_CONNECTIONS_PER_IMAGE}")
reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("RESULTADO")
reporte_lines.append("-" * 70)
reporte_lines.append(f"  Total de pares: {len(df_pares)}")
reporte_lines.append(f"    Pre-evento: {(df_pares['categoria'] == 'PRE-EVENTO').sum()}")
reporte_lines.append(f"    Co-evento: {(df_pares['categoria'] == 'CO-EVENTO').sum()}")
reporte_lines.append(f"    Post-evento: {(df_pares['categoria'] == 'POST-EVENTO').sum()}")
reporte_lines.append(f"  Baseline temporal media: {df_pares['temporal_baseline'].mean():.1f} días")
reporte_lines.append(f"  Baseline temporal mediana: {df_pares['temporal_baseline'].median():.0f} días")
reporte_lines.append(f"  Red conectada: {'SÍ' if is_connected else 'NO' if is_connected is not None else 'No verificado'}")

if is_connected is not None and not is_connected:
    reporte_lines.append(f"  Componentes conexos: {n_components}")
    reporte_lines.append(f"  Pares puente añadidos: {len(pares_adicionales) if 'pares_adicionales' in dir() else 0}")

reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("DISTRIBUCIÓN POR INTERVALO TEMPORAL")
reporte_lines.append("-" * 70)
for dt_val in sorted(df_pares['temporal_baseline'].unique()):
    count = (df_pares['temporal_baseline'] == dt_val).sum()
    if count > 0:
        bar = "#" * min(count // 5, 40)
        reporte_lines.append(f"  {dt_val:3d} días: {count:4d} pares  {bar}")

reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("PARÁMETROS HyP3 DE PROCESAMIENTO")
reporte_lines.append("-" * 70)
for k, v in HYP3_PARAMS.items():
    reporte_lines.append(f"  {k}: {v}")

reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("ESTIMACIÓN DE RECURSOS")
reporte_lines.append("-" * 70)
reporte_lines.append(f"  Interferogramas a procesar: {len(df_pares)}")
reporte_lines.append(f"  Tiempo estimado HyP3: {len(df_pares) // 80 + 1}-{len(df_pares) // 50 + 1} días")
reporte_lines.append(f"  Espacio en disco estimado: {len(df_pares) * 0.2:.0f}-{len(df_pares) * 0.35:.0f} GB")
reporte_lines.append(f"  Costo: Gratuito (NASA Earthdata / ASF On Demand)")

reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("SIGUIENTE PASO")
reporte_lines.append("-" * 70)
reporte_lines.append("  1. Revisar la red en las figuras de diagnóstico.")
reporte_lines.append("  2. Si la red es satisfactoria, cambiar ENVIAR_A_HYP3 = True")
reporte_lines.append("     y ejecutar de nuevo para enviar los jobs.")
reporte_lines.append("  3. Alternativamente, usar ASF Vertex SBAS Tool con los mismos")
reporte_lines.append("     parámetros para enviar manualmente.")
reporte_lines.append("")
reporte_lines.append("=" * 70)
reporte_lines.append("FIN DEL REPORTE P09")
reporte_lines.append("=" * 70)

reporte_text = "\n".join(reporte_lines)
reporte_path = os.path.join(OUTPUT_DIR, "P09_reporte_red_sbas.txt")
with open(reporte_path, 'w', encoding='utf-8') as f:
    f.write(reporte_text)
print(f"  Guardado: {reporte_path}")

print()
print(reporte_text)

# =============================================================================
# PASO 7: ENVÍO A HyP3 (solo si ENVIAR_A_HYP3 = True)
# =============================================================================

if ENVIAR_A_HYP3:
    print()
    print("-" * 70)
    print("PASO 7: Envío de jobs a HyP3")
    print("-" * 70)

    try:
        import hyp3_sdk

        print("Conectando a HyP3...")
        print("  (Requiere credenciales NASA Earthdata)")
        print("  Si no tienes cuenta: https://urs.earthdata.nasa.gov/users/new")
        print()

        # Conectar a HyP3 (pedirá credenciales interactivamente)
        hyp3 = hyp3_sdk.HyP3()
        print(f"  Conectado exitosamente.")

        # Enviar jobs uno por uno (misma API que P04/P04b)
        # hyp3.submit_insar_job() retorna Batch; job = batch.jobs[0]
        import time as _time
        all_submitted = []
        total_pares = len(df_pares)
        t_start = _time.time()

        print(f"  Enviando {total_pares} jobs a HyP3...")
        print(f"  (Esto tomará ~15-30 minutos para 1126 pares)")
        print()

        for i, (_, par) in enumerate(df_pares.iterrows(), 1):
            try:
                batch = hyp3.submit_insar_job(
                    granule1=par['ref_scene'],
                    granule2=par['sec_scene'],
                    name=HYP3_PROJECT_NAME,
                    **HYP3_PARAMS,
                )
                job = batch.jobs[0]
                all_submitted.append({
                    'job_id': job.job_id,
                    'ref_scene': par['ref_scene'],
                    'sec_scene': par['sec_scene'],
                    'ref_date': str(par['ref_date']),
                    'sec_date': str(par['sec_date']),
                    'temporal_baseline': int(par['temporal_baseline']),
                    'status': job.status_code,
                })

            except Exception as e:
                print(f"  ERROR par {i} ({par['ref_date']}-{par['sec_date']}): {e}")
                all_submitted.append({
                    'ref_scene': par['ref_scene'],
                    'sec_scene': par['sec_scene'],
                    'ref_date': str(par['ref_date']),
                    'sec_date': str(par['sec_date']),
                    'error': str(e),
                })

            # Progreso cada 25 jobs
            if i % 25 == 0 or i == total_pares:
                elapsed = _time.time() - t_start
                rate = i / elapsed if elapsed > 0 else 0
                eta = (total_pares - i) / rate / 60 if rate > 0 else 0
                n_ok = len([j for j in all_submitted if 'job_id' in j])
                n_err = len([j for j in all_submitted if 'error' in j])
                print(f"  [{i}/{total_pares}] {n_ok} OK, {n_err} errores "
                      f"| {elapsed/60:.1f} min transcurridos "
                      f"| ~{eta:.0f} min restantes")

            # Pausa para evitar rate-limiting (0.5s entre llamadas)
            _time.sleep(0.5)

            # Guardado incremental cada 100 jobs (seguridad ante desconexión)
            if i % 100 == 0:
                _backup = {
                    'proyecto': HYP3_PROJECT_NAME,
                    'parcial': True,
                    'enviados_hasta': i,
                    'total_pares': total_pares,
                    'fecha': datetime.now().isoformat(),
                    'jobs': all_submitted
                }
                _bkp_path = os.path.join(OUTPUT_DIR, "P09_jobs_parcial.json")
                with open(_bkp_path, 'w') as _f:
                    json.dump(_backup, _f, indent=2)
                print(f"    (backup guardado: {i} jobs)")

        # Guardar registro de jobs
        jobs_record = {
            'proyecto': HYP3_PROJECT_NAME,
            'fecha_envio': datetime.now().isoformat(),
            'total_jobs': len([j for j in all_submitted if 'job_id' in j]),
            'total_errores': len([j for j in all_submitted if 'error' in j]),
            'parametros': HYP3_PARAMS,
            'jobs': all_submitted
        }

        jobs_path = os.path.join(OUTPUT_DIR, "P09_jobs_enviados.json")
        with open(jobs_path, 'w') as f:
            json.dump(jobs_record, f, indent=2)
        print(f"\n  Registro guardado: {jobs_path}")

        n_ok = len([j for j in all_submitted if 'job_id' in j])
        n_err = len([j for j in all_submitted if 'error' in j])
        print(f"\n  ENVÍO COMPLETADO: {n_ok} jobs enviados, {n_err} errores.")
        print(f"  Monitorear en: https://hyp3-api.asf.alaska.edu/")
        print(f"  Proyecto: {HYP3_PROJECT_NAME}")

    except ImportError:
        print("ERROR: hyp3_sdk no instalada.")
        print("Instalar con: pip install hyp3_sdk")

    except Exception as e:
        print(f"ERROR durante envío: {e}")
        print("Verificar credenciales y conectividad.")

else:
    print("-" * 70)
    print("PASO 7: Envío a HyP3 OMITIDO (modo diseño)")
    print("-" * 70)
    print("  Para enviar, cambiar ENVIAR_A_HYP3 = True y ejecutar de nuevo.")
    print("  O usar ASF Vertex SBAS Tool con parámetros:")
    print(f"    - Temporal baseline: {TEMPORAL_BASELINE_MAX} días")
    print(f"    - Perpendicular baseline: {PERPENDICULAR_BASELINE_MAX} m")
    print(f"    - Track: Path 40 Descendente")
    print(f"    - Periodo: 2017-01-01 a 2024-06-30")

# =============================================================================
# RESUMEN FINAL
# =============================================================================

print()
print("=" * 70)
print("P09 COMPLETADO")
print("=" * 70)
print(f"Archivos generados en: {OUTPUT_DIR}")
print()
for f_name in sorted(os.listdir(OUTPUT_DIR)):
    f_path = os.path.join(OUTPUT_DIR, f_name)
    if f_name.startswith("P09"):
        size = os.path.getsize(f_path)
        print(f"  {f_name} ({size/1024:.1f} KB)")