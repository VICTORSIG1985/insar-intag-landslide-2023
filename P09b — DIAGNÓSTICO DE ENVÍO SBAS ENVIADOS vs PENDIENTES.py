# ============================================================
# P09b — DIAGNÓSTICO DE ENVÍO SBAS: ENVIADOS vs PENDIENTES
# ============================================================
# Analiza P09_jobs_enviados.json para determinar:
#   1. Qué pares se enviaron exitosamente (754)
#   2. Qué pares quedaron pendientes (372)
#   3. Cobertura temporal de cada grupo
#   4. Impacto en conectividad de la red
#   5. Exporta CSV de pares pendientes para reenvío
#
# Entrada: P09_jobs_enviados.json
# Salida:  P09c_diagnostico.txt
#          P09c_pares_pendientes.csv (para P09d)
#          P09c_pares_enviados.csv
#          P09c_diagnostico_temporal.png
# ============================================================
import json, os, sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import numpy as np

print("=" * 70)
print("P09c — DIAGNÓSTICO DE ENVÍO SBAS")
print("=" * 70)

# ============================================================
# 1. CONFIGURACIÓN
# ============================================================
SBAS_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"
JOBS_FILE = os.path.join(SBAS_DIR, "P09_jobs_enviados.json")
PARES_FILE = os.path.join(SBAS_DIR, "P09_red_sbas_pares.csv")
FECHA_EVENTO = pd.Timestamp("2023-12-19")

if not os.path.exists(JOBS_FILE):
    print(f"[ERROR] No se encuentra: {JOBS_FILE}")
    sys.exit(1)

# ============================================================
# 2. CARGAR DATOS
# ============================================================
print("\n--- Cargando datos ---")

with open(JOBS_FILE, 'r') as f:
    data = json.load(f)

jobs = data['jobs']
print(f"  Total de registros: {len(jobs)}")

# Separar enviados vs pendientes
enviados = [j for j in jobs if 'job_id' in j]
pendientes = [j for j in jobs if 'error' in j]

print(f"  Enviados exitosamente: {len(enviados)}")
print(f"  Pendientes (sin créditos): {len(pendientes)}")

# ============================================================
# 3. ANÁLISIS TEMPORAL
# ============================================================
print(f"\n{'─' * 70}")
print("ANÁLISIS TEMPORAL")
print(f"{'─' * 70}")

# Convertir a DataFrames
df_env = pd.DataFrame(enviados)
df_pend = pd.DataFrame(pendientes)

df_env['ref_date'] = pd.to_datetime(df_env['ref_date'])
df_env['sec_date'] = pd.to_datetime(df_env['sec_date'])
df_pend['ref_date'] = pd.to_datetime(df_pend['ref_date'])
df_pend['sec_date'] = pd.to_datetime(df_pend['sec_date'])

# Clasificar por relación con el evento
def clasificar(row):
    if row['ref_date'] < FECHA_EVENTO and row['sec_date'] > FECHA_EVENTO:
        return 'CO-EVENTO'
    elif row['sec_date'] <= FECHA_EVENTO:
        return 'PRE-EVENTO'
    else:
        return 'POST-EVENTO'

df_env['categoria'] = df_env.apply(clasificar, axis=1)
df_pend['categoria'] = df_pend.apply(clasificar, axis=1)

# Estadísticas enviados
print(f"\n  ENVIADOS (754 pares):")
print(f"    Periodo: {df_env['ref_date'].min().strftime('%Y-%m-%d')} a "
      f"{df_env['sec_date'].max().strftime('%Y-%m-%d')}")
env_cat = df_env['categoria'].value_counts()
for cat in ['PRE-EVENTO', 'CO-EVENTO', 'POST-EVENTO']:
    print(f"    {cat}: {env_cat.get(cat, 0)} pares")

# Estadísticas pendientes
print(f"\n  PENDIENTES (372 pares):")
print(f"    Periodo: {df_pend['ref_date'].min().strftime('%Y-%m-%d')} a "
      f"{df_pend['sec_date'].max().strftime('%Y-%m-%d')}")
pend_cat = df_pend['categoria'].value_counts()
for cat in ['PRE-EVENTO', 'CO-EVENTO', 'POST-EVENTO']:
    print(f"    {cat}: {pend_cat.get(cat, 0)} pares")

# ============================================================
# 4. ANÁLISIS DE CRITICIDAD
# ============================================================
print(f"\n{'─' * 70}")
print("ANÁLISIS DE CRITICIDAD DE PARES PENDIENTES")
print(f"{'─' * 70}")

# Pares co-evento pendientes (CRÍTICOS)
co_evento_pend = df_pend[df_pend['categoria'] == 'CO-EVENTO']
print(f"\n  PARES CO-EVENTO PENDIENTES: {len(co_evento_pend)}")
if len(co_evento_pend) > 0:
    print(f"  ⚠ CRÍTICO: Sin estos pares no se captura el evento del 19-dic-2023")
    for _, row in co_evento_pend.iterrows():
        print(f"    {row['ref_date'].strftime('%Y-%m-%d')} → "
              f"{row['sec_date'].strftime('%Y-%m-%d')}")

# Pares post-evento pendientes
post_pend = df_pend[df_pend['categoria'] == 'POST-EVENTO']
print(f"\n  PARES POST-EVENTO PENDIENTES: {len(post_pend)}")
if len(post_pend) > 0:
    print(f"  ⚠ IMPORTANTE: Necesarios para serie temporal post-deslizamiento")

# Pares puente pendientes
puente_pend = df_pend[
    ((df_pend['ref_date'].dt.strftime('%Y-%m-%d') == '2018-05-26') & 
     (df_pend['sec_date'].dt.strftime('%Y-%m-%d') == '2018-07-31')) |
    ((df_pend['ref_date'].dt.strftime('%Y-%m-%d') == '2019-05-03') & 
     (df_pend['sec_date'].dt.strftime('%Y-%m-%d') == '2019-07-02')) |
    ((df_pend['ref_date'].dt.strftime('%Y-%m-%d') == '2020-04-27') & 
     (df_pend['sec_date'].dt.strftime('%Y-%m-%d') == '2020-08-13'))
]
print(f"\n  PARES PUENTE PENDIENTES: {len(puente_pend)}")
if len(puente_pend) > 0:
    print(f"  ⚠ CRÍTICO: Sin pares puente la red queda DESCONECTADA")
    for _, row in puente_pend.iterrows():
        dt = (row['sec_date'] - row['ref_date']).days
        print(f"    {row['ref_date'].strftime('%Y-%m-%d')} → "
              f"{row['sec_date'].strftime('%Y-%m-%d')} ({dt} días)")

# ============================================================
# 5. CONECTIVIDAD DE RED CON SOLO LOS ENVIADOS
# ============================================================
print(f"\n{'─' * 70}")
print("VERIFICACIÓN DE CONECTIVIDAD (solo pares enviados)")
print(f"{'─' * 70}")

try:
    import networkx as nx
    
    G_env = nx.Graph()
    for _, row in df_env.iterrows():
        G_env.add_edge(row['ref_date'].strftime('%Y-%m-%d'),
                       row['sec_date'].strftime('%Y-%m-%d'))
    
    n_comp_env = nx.number_connected_components(G_env)
    conectada_env = n_comp_env == 1
    
    print(f"  Red con 754 pares enviados:")
    print(f"    Componentes conexos: {n_comp_env}")
    print(f"    Red conectada: {'SÍ' if conectada_env else 'NO'}")
    
    if not conectada_env:
        print(f"\n  Detalle de componentes:")
        for i, comp in enumerate(nx.connected_components(G_env), 1):
            fechas = sorted(comp)
            print(f"    Componente {i}: {len(comp)} escenas "
                  f"({fechas[0]} a {fechas[-1]})")
    
    # Red completa (enviados + pendientes)
    G_full = G_env.copy()
    for _, row in df_pend.iterrows():
        G_full.add_edge(row['ref_date'].strftime('%Y-%m-%d'),
                        row['sec_date'].strftime('%Y-%m-%d'))
    
    n_comp_full = nx.number_connected_components(G_full)
    print(f"\n  Red completa (754 + 372 = 1126 pares):")
    print(f"    Componentes conexos: {n_comp_full}")
    print(f"    Red conectada: {'SÍ' if n_comp_full == 1 else 'NO'}")

except ImportError:
    print("  networkx no disponible, saltando análisis de conectividad")

# ============================================================
# 6. RESUMEN CUANTITATIVO
# ============================================================
print(f"\n{'─' * 70}")
print("RESUMEN CUANTITATIVO")
print(f"{'─' * 70}")

print(f"""
  ┌────────────────────────────────────────────────┐
  │  ESTADO ACTUAL DEL ENVÍO SBAS                  │
  ├────────────────────────────────────────────────┤
  │  Total red SBAS:          1,126 pares          │
  │  ✓ Enviados (cuenta 1):     754 pares (67%)    │
  │  ✗ Pendientes:              372 pares (33%)    │
  ├────────────────────────────────────────────────┤
  │  ENVIADOS por categoría:                       │
  │    Pre-evento:   {env_cat.get('PRE-EVENTO', 0):>4} pares                    │
  │    Co-evento:    {env_cat.get('CO-EVENTO', 0):>4} pares                    │
  │    Post-evento:  {env_cat.get('POST-EVENTO', 0):>4} pares                    │
  ├────────────────────────────────────────────────┤
  │  PENDIENTES por categoría:                     │
  │    Pre-evento:   {pend_cat.get('PRE-EVENTO', 0):>4} pares                    │
  │    Co-evento:    {pend_cat.get('CO-EVENTO', 0):>4} pares                    │
  │    Post-evento:  {pend_cat.get('POST-EVENTO', 0):>4} pares                    │
  │    Puente:       {len(puente_pend):>4} pares                    │
  └────────────────────────────────────────────────┘
""")

# ============================================================
# 7. FIGURA DE DIAGNÓSTICO
# ============================================================
print("--- Generando figura de diagnóstico ---")

fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# --- Panel 1: Timeline de pares ---
ax = axes[0]
ax.set_title("Estado del envío SBAS — Pares enviados vs pendientes",
             fontsize=13, fontweight='bold')

# Pares enviados (verde)
for _, row in df_env.iterrows():
    mid = row['ref_date'] + (row['sec_date'] - row['ref_date']) / 2
    dt = (row['sec_date'] - row['ref_date']).days
    ax.plot([row['ref_date'], row['sec_date']], [dt, dt],
            color='#2ecc71', alpha=0.15, linewidth=0.5)

# Pares pendientes (rojo)
for _, row in df_pend.iterrows():
    mid = row['ref_date'] + (row['sec_date'] - row['ref_date']) / 2
    dt = (row['sec_date'] - row['ref_date']).days
    ax.plot([row['ref_date'], row['sec_date']], [dt, dt],
            color='#e74c3c', alpha=0.3, linewidth=0.8)

# Línea del evento
ax.axvline(FECHA_EVENTO, color='red', linewidth=2, linestyle='--', label='Evento 19-dic-2023')
ax.set_xlabel("Fecha")
ax.set_ylabel("Baseline temporal (días)")
ax.legend(loc='upper left')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.xaxis.set_major_locator(mdates.YearLocator())

# Leyenda manual
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='#2ecc71', linewidth=2, label=f'Enviados ({len(enviados)})'),
    Line2D([0], [0], color='#e74c3c', linewidth=2, label=f'Pendientes ({len(pendientes)})'),
    Line2D([0], [0], color='red', linewidth=2, linestyle='--', label='Evento 19-dic-2023'),
]
ax.legend(handles=legend_elements, loc='upper right')

# --- Panel 2: Distribución temporal ---
ax2 = axes[1]
ax2.set_title("Cobertura temporal: escenas cubiertas por cuenta", fontsize=13, fontweight='bold')

# Todas las fechas únicas
all_dates_env = sorted(set(df_env['ref_date'].tolist() + df_env['sec_date'].tolist()))
all_dates_pend = sorted(set(df_pend['ref_date'].tolist() + df_pend['sec_date'].tolist()))
all_dates = sorted(set(all_dates_env + all_dates_pend))

# Contar conexiones por fecha
env_counts = {}
for _, row in df_env.iterrows():
    d1 = row['ref_date'].strftime('%Y-%m-%d')
    d2 = row['sec_date'].strftime('%Y-%m-%d')
    env_counts[d1] = env_counts.get(d1, 0) + 1
    env_counts[d2] = env_counts.get(d2, 0) + 1

pend_counts = {}
for _, row in df_pend.iterrows():
    d1 = row['ref_date'].strftime('%Y-%m-%d')
    d2 = row['sec_date'].strftime('%Y-%m-%d')
    pend_counts[d1] = pend_counts.get(d1, 0) + 1
    pend_counts[d2] = pend_counts.get(d2, 0) + 1

dates_plot = [pd.Timestamp(d) for d in sorted(set(list(env_counts.keys()) + list(pend_counts.keys())))]
env_vals = [env_counts.get(d.strftime('%Y-%m-%d'), 0) for d in dates_plot]
pend_vals = [pend_counts.get(d.strftime('%Y-%m-%d'), 0) for d in dates_plot]

ax2.bar(dates_plot, env_vals, width=8, color='#2ecc71', alpha=0.7, label=f'Cuenta 1 ({len(enviados)} pares)')
ax2.bar(dates_plot, pend_vals, width=8, bottom=env_vals, color='#e74c3c', alpha=0.7, label=f'Cuenta 2 ({len(pendientes)} pares)')
ax2.axvline(FECHA_EVENTO, color='red', linewidth=2, linestyle='--')
ax2.set_xlabel("Fecha")
ax2.set_ylabel("Conexiones por escena")
ax2.legend()
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax2.xaxis.set_major_locator(mdates.YearLocator())

plt.tight_layout()
fig_path = os.path.join(SBAS_DIR, "P09c_diagnostico_temporal.png")
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Guardado: {fig_path}")

# ============================================================
# 8. EXPORTAR CSVs
# ============================================================
print(f"\n--- Exportando CSVs ---")

# CSV de pares pendientes (para P09d)
df_pend_export = df_pend[['ref_scene', 'sec_scene', 'ref_date', 'sec_date', 'categoria']].copy()
df_pend_export['ref_date'] = df_pend_export['ref_date'].dt.strftime('%Y-%m-%d')
df_pend_export['sec_date'] = df_pend_export['sec_date'].dt.strftime('%Y-%m-%d')
pend_csv = os.path.join(SBAS_DIR, "P09c_pares_pendientes.csv")
df_pend_export.to_csv(pend_csv, index=False)
print(f"  ✓ Pares pendientes: {pend_csv} ({len(df_pend_export)} pares)")

# CSV de pares enviados
df_env_export = df_env[['job_id', 'ref_scene', 'sec_scene', 'ref_date', 'sec_date', 
                         'temporal_baseline', 'status', 'categoria']].copy()
df_env_export['ref_date'] = df_env_export['ref_date'].dt.strftime('%Y-%m-%d')
df_env_export['sec_date'] = df_env_export['sec_date'].dt.strftime('%Y-%m-%d')
env_csv = os.path.join(SBAS_DIR, "P09c_pares_enviados.csv")
df_env_export.to_csv(env_csv, index=False)
print(f"  ✓ Pares enviados: {env_csv} ({len(df_env_export)} pares)")

# ============================================================
# 9. REPORTE COMPLETO
# ============================================================
reporte = []
reporte.append("=" * 70)
reporte.append("P09c — DIAGNÓSTICO DE ENVÍO SBAS")
reporte.append("=" * 70)
reporte.append(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
reporte.append(f"")
reporte.append(f"RED SBAS TOTAL: 1,126 pares interferométricos")
reporte.append(f"")
reporte.append(f"CUENTA 1 (enviados): {len(enviados)} pares")
reporte.append(f"  Periodo: {df_env['ref_date'].min().strftime('%Y-%m-%d')} a "
               f"{df_env['sec_date'].max().strftime('%Y-%m-%d')}")
reporte.append(f"  Pre-evento:  {env_cat.get('PRE-EVENTO', 0)}")
reporte.append(f"  Co-evento:   {env_cat.get('CO-EVENTO', 0)}")
reporte.append(f"  Post-evento: {env_cat.get('POST-EVENTO', 0)}")
reporte.append(f"")
reporte.append(f"CUENTA 2 (pendientes): {len(pendientes)} pares")
reporte.append(f"  Periodo: {df_pend['ref_date'].min().strftime('%Y-%m-%d')} a "
               f"{df_pend['sec_date'].max().strftime('%Y-%m-%d')}")
reporte.append(f"  Pre-evento:  {pend_cat.get('PRE-EVENTO', 0)}")
reporte.append(f"  Co-evento:   {pend_cat.get('CO-EVENTO', 0)}")
reporte.append(f"  Post-evento: {pend_cat.get('POST-EVENTO', 0)}")
reporte.append(f"  Puente:      {len(puente_pend)}")
reporte.append(f"")
reporte.append(f"PARES CRÍTICOS PENDIENTES:")
reporte.append(f"  Co-evento:   {len(co_evento_pend)} (capturan el deslizamiento)")
reporte.append(f"  Post-evento: {len(post_pend)} (serie temporal post-evento)")
reporte.append(f"  Puente:      {len(puente_pend)} (conectividad de red)")
reporte.append(f"")
reporte.append(f"ACCIÓN REQUERIDA:")
reporte.append(f"  Crear segunda cuenta NASA Earthdata")
reporte.append(f"  Ejecutar P09d con P09c_pares_pendientes.csv")
reporte.append(f"  Los 372 pares se procesan en paralelo con los 754")

reporte_path = os.path.join(SBAS_DIR, "P09c_diagnostico.txt")
with open(reporte_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(reporte))
print(f"\n  ✓ Reporte: {reporte_path}")

# ============================================================
# 10. RESUMEN FINAL
# ============================================================
print(f"\n{'=' * 70}")
print("P09c COMPLETADO")
print(f"{'=' * 70}")
print(f"\n  Archivos generados:")
for fname in [pend_csv, env_csv, fig_path, reporte_path]:
    sz = os.path.getsize(fname) / 1024
    print(f"    {os.path.basename(fname)} ({sz:.1f} KB)")

print(f"\n  SIGUIENTE PASO:")
print(f"  1. Crear nueva cuenta: https://urs.earthdata.nasa.gov/users/new")
print(f"  2. Configurar credenciales en _netrc")
print(f"  3. Ejecutar P09d_envio_pendientes.py")
print(f"{'=' * 70}")