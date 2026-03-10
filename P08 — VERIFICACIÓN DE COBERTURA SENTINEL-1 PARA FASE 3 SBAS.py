# =============================================================================
# P08 — VERIFICACIÓN DE COBERTURA SENTINEL-1 PARA FASE 3 SBAS
# Versión: 1.0
# Fecha: 2026-02-28
# Autor: Víctor Pinto
# Proyecto: Detección de deslizamientos mediante InSAR — Zona de Intag, Cotacachi
# Fase: Preparación para Fase 3 (SBAS Time Series)
#
# OBJETIVO:
# Verificar geométricamente que los tracks Sentinel-1 usados en Fase 2
# (DInSAR co-evento) cubren el 100% del AOI para la Fase 3 (SBAS).
# Se verifica Path 40 (descendente, cobertura primaria) y Path 120
# (ascendente, cobertura secundaria) usando footprints reales de ASF.
#
# JUSTIFICACIÓN CIENTÍFICA:
# El SBAS debe usar los mismos tracks de Fase 2 para garantizar:
# - Misma geometría de observación (ángulo de incidencia, LOS)
# - Sensibilidad equivalente al desplazamiento
# - Comparabilidad directa entre magnitudes DInSAR y velocidades SBAS
# (Referencia: Berardino et al. 2002; Lanari et al. 2007)
#
# INPUTS:
# - D:\POSGRADOS\INTAG\INTAJ.gpkg (AOI: 6 parroquias de Intag)
# - D:\POSGRADOS\INTAG\data\sentinel1\pares_maestro_v2.json
# - API de ASF (asf_search) para footprints Sentinel-1 SLC
#
# OUTPUTS (en D:\POSGRADOS\INTAG\data\sentinel1\diagnostico_cobertura\fase3\):
# - P08_cobertura_path40_fase3.png (mapa de cobertura)
# - P08_cobertura_path120_fase3.png (mapa de cobertura)
# - P08_disponibilidad_temporal.png (gráfico de escenas por mes)
# - P08_gaps_temporales.png (gráfico de gaps entre adquisiciones)
# - P08_resumen_cobertura_fase3.csv (estadísticas por track)
# - P08_escenas_path40_fase3.csv (catálogo completo de escenas)
# - P08_escenas_path120_fase3.csv (catálogo completo de escenas)
# - P08_reporte_fase3.txt (reporte textual completo)
#
# DEPENDENCIAS:
# pip install asf_search geopandas shapely matplotlib pandas numpy
# =============================================================================

import os
import json
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from shapely.geometry import shape, box, Polygon
from shapely.ops import unary_union
from datetime import datetime, timedelta

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

# Rutas de entrada
AOI_PATH = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"
PARES_PATH = r"D:\POSGRADOS\INTAG\data\sentinel1\pares_maestro_v2.json"
BBOX_PATH = r"D:\POSGRADOS\INTAG\data\sentinel1\bbox_area_estudio.json"

# Ruta de salida
OUTPUT_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\diagnostico_cobertura\fase3"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Periodo de análisis SBAS (justificado por Dini et al. 2026: 5+ años pre-evento)
FECHA_EVENTO = "2023-12-19"
SBAS_START = "2017-01-01"
SBAS_END = "2024-06-30"

# Tracks de Fase 2 (extraídos de pares_maestro_v2.json)
TRACKS = {
    "Path_40": {
        "path": 40,
        "direction": "DESCENDING",
        "frame": 590,
        "pares_fase2": ["D1", "D2", "D3", "D4"],
        "cobertura_fase2": "100% (6/6 parroquias)"
    },
    "Path_120": {
        "path": 120,
        "direction": "ASCENDING",
        "frame": 1182,
        "pares_fase2": ["A1", "A2", "A3", "A4"],
        "cobertura_fase2": "5/6 parroquias (García Moreno parcial 46.8%)"
    },
    "Path_18": {
        "path": 18,
        "direction": "ASCENDING",
        "frame": 1180,
        "pares_fase2": ["P18_A3", "P18_B1", "P18_B2", "P18_B3"],
        "cobertura_fase2": "~69.8% AOI"
    }
}

print("=" * 70)
print("P08 — VERIFICACIÓN DE COBERTURA SENTINEL-1 PARA FASE 3 SBAS")
print("=" * 70)
print(f"Fecha de ejecución: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Periodo SBAS: {SBAS_START} a {SBAS_END}")
print(f"Evento de referencia: {FECHA_EVENTO}")
print()

# =============================================================================
# PASO 1: CARGAR AOI
# =============================================================================

print("-" * 70)
print("PASO 1: Carga del AOI")
print("-" * 70)

# Intentar cargar desde INTAJ.gpkg
try:
    aoi_gdf = gpd.read_file(AOI_PATH)
    print(f"AOI cargado desde: {AOI_PATH}")
    print(f"  CRS: {aoi_gdf.crs}")
    print(f"  Geometrías: {len(aoi_gdf)}")
    
    # Reproyectar a WGS84 si es necesario
    if aoi_gdf.crs and aoi_gdf.crs.to_epsg() != 4326:
        aoi_gdf_wgs84 = aoi_gdf.to_crs(epsg=4326)
        print(f"  Reproyectado a EPSG:4326")
    else:
        aoi_gdf_wgs84 = aoi_gdf.copy()
    
    # Unión de todas las geometrías para obtener AOI completo
    aoi_union = unary_union(aoi_gdf_wgs84.geometry)
    aoi_area_km2 = aoi_gdf.to_crs(epsg=32617).geometry.area.sum() / 1e6
    
    print(f"  Área total AOI: {aoi_area_km2:.2f} km²")
    print(f"  Bounds (WGS84): {aoi_union.bounds}")
    
    # Listar parroquias si hay columna de nombre
    name_cols = [c for c in aoi_gdf.columns if c.lower() in 
                 ['nombre', 'name', 'parroquia', 'dpa_despar', 'dpa_parroq', 'nam']]
    if name_cols:
        print(f"  Parroquias:")
        for _, row in aoi_gdf.iterrows():
            print(f"    - {row[name_cols[0]]}")
    
    aoi_source = "INTAJ.gpkg"

except Exception as e:
    print(f"  ADVERTENCIA: No se pudo cargar INTAJ.gpkg: {e}")
    print(f"  Usando bbox_area_estudio.json como respaldo...")
    
    with open(BBOX_PATH, 'r') as f:
        bbox_data = json.load(f)
    
    ext = bbox_data["extent_original"]
    aoi_union = box(ext["lon_min"], ext["lat_min"], ext["lon_max"], ext["lat_max"])
    aoi_gdf_wgs84 = gpd.GeoDataFrame(geometry=[aoi_union], crs="EPSG:4326")
    aoi_area_km2 = aoi_gdf_wgs84.to_crs(epsg=32617).geometry.area.sum() / 1e6
    
    print(f"  Área total AOI (bbox): {aoi_area_km2:.2f} km²")
    print(f"  Bounds: {aoi_union.bounds}")
    
    aoi_source = "bbox_area_estudio.json"

print()

# =============================================================================
# PASO 2: CONSULTAR ASF PARA CADA TRACK
# =============================================================================

print("-" * 70)
print("PASO 2: Consulta de escenas Sentinel-1 en ASF")
print("-" * 70)

try:
    import asf_search as asf
    print("Librería asf_search disponible.")
except ImportError:
    print("ERROR: asf_search no instalada.")
    print("Instalar con: pip install asf_search")
    print()
    print("ALTERNATIVA: Usando búsqueda por API REST de ASF...")

# Función para buscar escenas en ASF
def buscar_escenas_asf(path_number, direction, start_date, end_date, aoi_wkt):
    """
    Busca escenas Sentinel-1 SLC en ASF para un track específico.
    Retorna DataFrame con metadatos y footprints.
    """
    try:
        import asf_search as asf
        
        # Construir búsqueda
        results = asf.search(
            platform=asf.PLATFORM.SENTINEL1,
            processingLevel=asf.PRODUCT_TYPE.SLC,
            relativeOrbit=path_number,
            flightDirection=direction.upper(),
            intersectsWith=aoi_wkt,
            start=start_date,
            end=end_date,
            beamMode=asf.BEAMMODE.IW
        )
        
        print(f"  Encontradas {len(results)} escenas para Path {path_number} ({direction})")
        
        if len(results) == 0:
            return pd.DataFrame()
        
        # Extraer metadatos
        records = []
        for r in results:
            props = r.properties
            geom = r.geojson()['geometry']
            
            records.append({
                'scene_id': props.get('sceneName', ''),
                'start_time': props.get('startTime', ''),
                'stop_time': props.get('stopTime', ''),
                'path': props.get('pathNumber', path_number),
                'frame': props.get('frameNumber', ''),
                'direction': direction,
                'beam_mode': props.get('beamModeType', ''),
                'polarization': props.get('polarization', ''),
                'geometry_wkt': shape(geom).wkt,
                'footprint': shape(geom)
            })
        
        df = pd.DataFrame(records)
        df['date'] = pd.to_datetime(df['start_time']).dt.date
        df['year'] = pd.to_datetime(df['start_time']).dt.year
        df['month'] = pd.to_datetime(df['start_time']).dt.month
        
        # Eliminar duplicados por fecha (misma órbita)
        df = df.drop_duplicates(subset=['date']).sort_values('date').reset_index(drop=True)
        
        return df
    
    except ImportError:
        # Alternativa usando API REST
        return buscar_escenas_rest(path_number, direction, start_date, end_date, aoi_wkt)

def buscar_escenas_rest(path_number, direction, start_date, end_date, aoi_wkt):
    """
    Alternativa usando la API REST de ASF cuando asf_search no está disponible.
    """
    import urllib.request
    import urllib.parse
    
    base_url = "https://api.daac.asf.alaska.edu/services/search/param"
    
    params = {
        'platform': 'SENTINEL-1',
        'processingLevel': 'SLC',
        'beamMode': 'IW',
        'relativeOrbit': str(path_number),
        'flightDirection': direction.upper(),
        'intersectsWith': aoi_wkt,
        'start': start_date,
        'end': end_date,
        'output': 'geojson',
        'maxResults': 5000
    }
    
    url = base_url + '?' + urllib.parse.urlencode(params)
    
    print(f"  Consultando API REST de ASF para Path {path_number}...")
    
    try:
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'P08_SBAS_Coverage_Check/1.0')
        
        with urllib.request.urlopen(req, timeout=120) as response:
            data = json.loads(response.read().decode())
        
        features = data.get('features', [])
        print(f"  Encontradas {len(features)} escenas para Path {path_number} ({direction})")
        
        if len(features) == 0:
            return pd.DataFrame()
        
        records = []
        for f in features:
            props = f['properties']
            geom = f['geometry']
            
            records.append({
                'scene_id': props.get('sceneName', props.get('fileID', '')),
                'start_time': props.get('startTime', ''),
                'stop_time': props.get('stopTime', ''),
                'path': path_number,
                'frame': props.get('frameNumber', props.get('frame', '')),
                'direction': direction,
                'beam_mode': props.get('beamModeType', 'IW'),
                'polarization': props.get('polarization', ''),
                'geometry_wkt': shape(geom).wkt,
                'footprint': shape(geom)
            })
        
        df = pd.DataFrame(records)
        df['date'] = pd.to_datetime(df['start_time']).dt.date
        df['year'] = pd.to_datetime(df['start_time']).dt.year
        df['month'] = pd.to_datetime(df['start_time']).dt.month
        
        df = df.drop_duplicates(subset=['date']).sort_values('date').reset_index(drop=True)
        
        return df
    
    except Exception as e:
        print(f"  ERROR en consulta REST: {e}")
        return pd.DataFrame()

# Generar WKT del AOI para la búsqueda
aoi_wkt = aoi_union.wkt

# Buscar escenas para cada track
resultados_tracks = {}

for track_name, track_info in TRACKS.items():
    print(f"\n  Buscando {track_name} (Path {track_info['path']}, {track_info['direction']})...")
    
    df = buscar_escenas_asf(
        path_number=track_info['path'],
        direction=track_info['direction'],
        start_date=SBAS_START,
        end_date=SBAS_END,
        aoi_wkt=aoi_wkt
    )
    
    if len(df) > 0:
        resultados_tracks[track_name] = {
            'info': track_info,
            'escenas': df,
            'n_escenas': len(df)
        }
        
        # Resumen por año
        by_year = df.groupby('year').size()
        print(f"  Distribución por año:")
        for year, count in by_year.items():
            print(f"    {year}: {count} escenas")
    else:
        print(f"  SIN DATOS para {track_name}")
        resultados_tracks[track_name] = {
            'info': track_info,
            'escenas': pd.DataFrame(),
            'n_escenas': 0
        }

print()

# =============================================================================
# PASO 3: VERIFICACIÓN GEOMÉTRICA DE COBERTURA
# =============================================================================

print("-" * 70)
print("PASO 3: Verificación geométrica de cobertura sobre AOI")
print("-" * 70)

cobertura_resultados = {}

for track_name, data in resultados_tracks.items():
    df = data['escenas']
    
    if len(df) == 0:
        cobertura_resultados[track_name] = {
            'cobertura_pct': 0,
            'area_cubierta_km2': 0,
            'area_no_cubierta_km2': aoi_area_km2,
            'n_escenas': 0
        }
        continue
    
    print(f"\n  {track_name}:")
    
    # Unión de todos los footprints
    all_footprints = unary_union(df['footprint'].tolist())
    
    # Intersección con AOI
    interseccion = aoi_union.intersection(all_footprints)
    area_interseccion = gpd.GeoDataFrame(
        geometry=[interseccion], crs="EPSG:4326"
    ).to_crs(epsg=32617).geometry.area.values[0] / 1e6
    
    # Zona no cubierta
    zona_no_cubierta = aoi_union.difference(all_footprints)
    area_no_cubierta = gpd.GeoDataFrame(
        geometry=[zona_no_cubierta], crs="EPSG:4326"
    ).to_crs(epsg=32617).geometry.area.values[0] / 1e6
    
    cobertura_pct = (area_interseccion / aoi_area_km2) * 100
    
    cobertura_resultados[track_name] = {
        'cobertura_pct': cobertura_pct,
        'area_cubierta_km2': area_interseccion,
        'area_no_cubierta_km2': area_no_cubierta,
        'n_escenas': len(df),
        'footprint_union': all_footprints,
        'interseccion': interseccion,
        'zona_no_cubierta': zona_no_cubierta
    }
    
    print(f"  Cobertura: {cobertura_pct:.2f}%")
    print(f"  Área cubierta: {area_interseccion:.2f} km²")
    print(f"  Área NO cubierta: {area_no_cubierta:.2f} km²")
    print(f"  Escenas totales: {len(df)}")
    
    # Verificar por parroquia si tenemos geometrías individuales
    if aoi_source == "INTAJ.gpkg" and name_cols:
        print(f"  Cobertura por parroquia:")
        for _, row in aoi_gdf_wgs84.iterrows():
            parr_geom = row.geometry
            parr_inter = parr_geom.intersection(all_footprints)
            parr_area = gpd.GeoDataFrame(
                geometry=[parr_geom], crs="EPSG:4326"
            ).to_crs(epsg=32617).geometry.area.values[0] / 1e6
            parr_cov = gpd.GeoDataFrame(
                geometry=[parr_inter], crs="EPSG:4326"
            ).to_crs(epsg=32617).geometry.area.values[0] / 1e6
            pct = (parr_cov / parr_area) * 100 if parr_area > 0 else 0
            print(f"    {row[name_cols[0]]}: {pct:.1f}% ({parr_cov:.1f}/{parr_area:.1f} km²)")

print()

# =============================================================================
# PASO 4: ANÁLISIS TEMPORAL (GAPS, CONTINUIDAD)
# =============================================================================

print("-" * 70)
print("PASO 4: Análisis temporal de continuidad")
print("-" * 70)

fecha_evento = pd.Timestamp(FECHA_EVENTO)

for track_name, data in resultados_tracks.items():
    df = data['escenas']
    
    if len(df) == 0:
        continue
    
    print(f"\n  {track_name}:")
    
    # Calcular gaps entre adquisiciones consecutivas
    dates = pd.to_datetime(df['date']).sort_values().reset_index(drop=True)
    gaps = dates.diff().dt.days.dropna()
    
    print(f"  Periodo: {dates.min().strftime('%Y-%m-%d')} a {dates.max().strftime('%Y-%m-%d')}")
    print(f"  Revisita mediana: {gaps.median():.0f} días")
    print(f"  Revisita media: {gaps.mean():.1f} días")
    print(f"  Gap máximo: {gaps.max():.0f} días")
    print(f"  Gap mínimo: {gaps.min():.0f} días")
    
    # Identificar gaps > 24 días (más de 2 ciclos de revisita)
    gaps_grandes = gaps[gaps > 24]
    if len(gaps_grandes) > 0:
        print(f"  Gaps > 24 días: {len(gaps_grandes)}")
        # Mostrar los 5 gaps más grandes
        gap_indices = gaps_grandes.nlargest(5).index
        for idx in gap_indices:
            fecha_antes = dates.iloc[idx - 1].strftime('%Y-%m-%d')
            fecha_despues = dates.iloc[idx].strftime('%Y-%m-%d')
            print(f"    {fecha_antes} → {fecha_despues} ({gaps.iloc[idx]:.0f} días)")
    else:
        print(f"  Sin gaps > 24 días")
    
    # Escenas pre-evento y post-evento
    n_pre = (dates < fecha_evento).sum()
    n_post = (dates >= fecha_evento).sum()
    print(f"  Escenas pre-evento (<{FECHA_EVENTO}): {n_pre}")
    print(f"  Escenas post-evento (>={FECHA_EVENTO}): {n_post}")
    
    # Estimación de interferogramas SBAS (umbral temporal 60 días)
    # Cada imagen se conecta con vecinos dentro de 60 días
    n_ifg_estimate = 0
    for i in range(len(dates)):
        for j in range(i + 1, len(dates)):
            delta = (dates.iloc[j] - dates.iloc[i]).days
            if delta <= 60:
                n_ifg_estimate += 1
            elif delta > 60:
                break
    
    print(f"  Interferogramas estimados (umbral 60 días): ~{n_ifg_estimate}")
    
    # Guardar estadísticas
    data['gaps'] = gaps
    data['dates'] = dates
    data['n_pre'] = n_pre
    data['n_post'] = n_post
    data['n_ifg_estimate'] = n_ifg_estimate

print()

# =============================================================================
# PASO 5: GENERAR FIGURAS
# =============================================================================

print("-" * 70)
print("PASO 5: Generación de figuras")
print("-" * 70)

# --- Figura 1: Mapa de cobertura por track ---
for track_name, cob in cobertura_resultados.items():
    if cob['n_escenas'] == 0:
        continue
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Dibujar AOI
    if aoi_source == "INTAJ.gpkg":
        aoi_gdf_wgs84.plot(ax=ax, facecolor='lightblue', edgecolor='black', 
                           linewidth=1.5, alpha=0.5, label='AOI (parroquias)')
    else:
        gpd.GeoDataFrame(geometry=[aoi_union], crs="EPSG:4326").plot(
            ax=ax, facecolor='lightblue', edgecolor='black', 
            linewidth=1.5, alpha=0.5, label='AOI')
    
    # Dibujar intersección
    if not cob['interseccion'].is_empty:
        gpd.GeoDataFrame(geometry=[cob['interseccion']], crs="EPSG:4326").plot(
            ax=ax, facecolor='green', edgecolor='green', 
            linewidth=0.5, alpha=0.3, label=f'Cubierto ({cob["cobertura_pct"]:.1f}%)')
    
    # Dibujar zona no cubierta
    if not cob['zona_no_cubierta'].is_empty and cob['area_no_cubierta_km2'] > 0.01:
        gpd.GeoDataFrame(geometry=[cob['zona_no_cubierta']], crs="EPSG:4326").plot(
            ax=ax, facecolor='red', edgecolor='red', 
            linewidth=0.5, alpha=0.5, label=f'NO cubierto ({cob["area_no_cubierta_km2"]:.1f} km²)')
    
    # Dibujar footprint de una escena representativa
    df_track = resultados_tracks[track_name]['escenas']
    if len(df_track) > 0:
        mid_idx = len(df_track) // 2
        sample_fp = df_track.iloc[mid_idx]['footprint']
        gpd.GeoDataFrame(geometry=[sample_fp], crs="EPSG:4326").plot(
            ax=ax, facecolor='none', edgecolor='blue', 
            linewidth=1, linestyle='--', alpha=0.7, label='Footprint ejemplo')
    
    track_info = TRACKS[track_name]
    ax.set_title(f"Cobertura {track_name} (Path {track_info['path']}, "
                 f"{track_info['direction']})\n"
                 f"Cobertura AOI: {cob['cobertura_pct']:.1f}% | "
                 f"Escenas: {cob['n_escenas']} ({SBAS_START[:4]}–{SBAS_END[:4]})",
                 fontsize=13, fontweight='bold')
    ax.set_xlabel("Longitud")
    ax.set_ylabel("Latitud")
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    fig_path = os.path.join(OUTPUT_DIR, f"P08_cobertura_{track_name}_fase3.png")
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Guardado: {fig_path}")

# --- Figura 2: Disponibilidad temporal (escenas por mes) ---
fig, axes = plt.subplots(len([t for t in resultados_tracks if resultados_tracks[t]['n_escenas'] > 0]), 
                         1, figsize=(14, 4 * len([t for t in resultados_tracks if resultados_tracks[t]['n_escenas'] > 0])),
                         sharex=True)
if not hasattr(axes, '__len__'):
    axes = [axes]

plot_idx = 0
for track_name, data in resultados_tracks.items():
    if data['n_escenas'] == 0:
        continue
    
    ax = axes[plot_idx]
    df = data['escenas']
    
    # Agrupar por mes
    df_temp = df.copy()
    df_temp['year_month'] = pd.to_datetime(df_temp['date']).dt.to_period('M')
    monthly = df_temp.groupby('year_month').size()
    
    # Crear rango completo de meses
    all_months = pd.period_range(start=SBAS_START, end=SBAS_END, freq='M')
    monthly_full = monthly.reindex(all_months, fill_value=0)
    
    colors = ['#2196F3' if m.to_timestamp() < fecha_evento else '#FF5722' 
              for m in monthly_full.index]
    
    ax.bar(range(len(monthly_full)), monthly_full.values, color=colors, alpha=0.7, width=0.8)
    
    # Línea vertical del evento
    evento_idx = None
    for i, m in enumerate(monthly_full.index):
        if m.to_timestamp().year == 2023 and m.to_timestamp().month == 12:
            evento_idx = i
            break
    
    if evento_idx:
        ax.axvline(x=evento_idx, color='red', linestyle='--', linewidth=2, 
                   label=f'Evento {FECHA_EVENTO}')
    
    track_info = TRACKS[track_name]
    ax.set_title(f"{track_name} (Path {track_info['path']}, {track_info['direction']}) — "
                 f"{data['n_escenas']} escenas", fontsize=11, fontweight='bold')
    ax.set_ylabel("Escenas/mes")
    ax.legend(fontsize=9)
    
    # Etiquetas de año
    tick_positions = []
    tick_labels = []
    for i, m in enumerate(monthly_full.index):
        if m.month == 1:
            tick_positions.append(i)
            tick_labels.append(str(m.year))
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.grid(True, axis='y', alpha=0.3)
    
    plot_idx += 1

plt.suptitle("Disponibilidad temporal Sentinel-1 SLC para SBAS\n"
             f"Periodo: {SBAS_START} a {SBAS_END}", fontsize=14, fontweight='bold')
plt.tight_layout()

fig_path = os.path.join(OUTPUT_DIR, "P08_disponibilidad_temporal.png")
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Guardado: {fig_path}")

# --- Figura 3: Distribución de gaps ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for idx, track_name in enumerate(["Path_40", "Path_120"]):
    if track_name not in resultados_tracks or resultados_tracks[track_name]['n_escenas'] == 0:
        continue
    
    data = resultados_tracks[track_name]
    if 'gaps' not in data:
        continue
    
    ax = axes[idx]
    gaps = data['gaps']
    
    ax.hist(gaps, bins=range(0, int(gaps.max()) + 15, 6), 
            color='#2196F3', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(x=12, color='green', linestyle='--', linewidth=1.5, 
               label='Revisita nominal (12d)')
    ax.axvline(x=24, color='orange', linestyle='--', linewidth=1.5, 
               label='2 ciclos (24d)')
    ax.axvline(x=gaps.median(), color='red', linestyle='-', linewidth=2, 
               label=f'Mediana ({gaps.median():.0f}d)')
    
    track_info = TRACKS[track_name]
    ax.set_title(f"{track_name} (Path {track_info['path']})\n"
                 f"Mediana: {gaps.median():.0f}d | Máx: {gaps.max():.0f}d", 
                 fontsize=11, fontweight='bold')
    ax.set_xlabel("Gap entre adquisiciones (días)")
    ax.set_ylabel("Frecuencia")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.suptitle("Distribución de gaps temporales entre adquisiciones", 
             fontsize=13, fontweight='bold')
plt.tight_layout()

fig_path = os.path.join(OUTPUT_DIR, "P08_gaps_temporales.png")
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Guardado: {fig_path}")

print()

# =============================================================================
# PASO 6: EXPORTAR DATOS
# =============================================================================

print("-" * 70)
print("PASO 6: Exportación de datos")
print("-" * 70)

# CSV de resumen de cobertura
resumen_rows = []
for track_name, cob in cobertura_resultados.items():
    track_info = TRACKS[track_name]
    data = resultados_tracks[track_name]
    
    row = {
        'track': track_name,
        'path_number': track_info['path'],
        'direction': track_info['direction'],
        'frame': track_info['frame'],
        'pares_fase2': ', '.join(track_info['pares_fase2']),
        'cobertura_fase2': track_info['cobertura_fase2'],
        'n_escenas_sbas': cob['n_escenas'],
        'cobertura_aoi_pct': round(cob['cobertura_pct'], 2),
        'area_cubierta_km2': round(cob['area_cubierta_km2'], 2),
        'area_no_cubierta_km2': round(cob['area_no_cubierta_km2'], 2),
        'n_pre_evento': data.get('n_pre', 0),
        'n_post_evento': data.get('n_post', 0),
        'n_ifg_estimados': data.get('n_ifg_estimate', 0),
        'apto_sbas': 'SI' if cob['cobertura_pct'] >= 95 and cob['n_escenas'] >= 50 else 'NO'
    }
    resumen_rows.append(row)

resumen_df = pd.DataFrame(resumen_rows)
resumen_path = os.path.join(OUTPUT_DIR, "P08_resumen_cobertura_fase3.csv")
resumen_df.to_csv(resumen_path, index=False)
print(f"  Guardado: {resumen_path}")

# CSV de escenas por track
for track_name, data in resultados_tracks.items():
    df = data['escenas']
    if len(df) == 0:
        continue
    
    # Exportar sin la columna footprint (geometría pesada)
    export_cols = [c for c in df.columns if c != 'footprint']
    csv_path = os.path.join(OUTPUT_DIR, f"P08_escenas_{track_name}_fase3.csv")
    df[export_cols].to_csv(csv_path, index=False)
    print(f"  Guardado: {csv_path}")

print()

# =============================================================================
# PASO 7: REPORTE FINAL
# =============================================================================

print("-" * 70)
print("PASO 7: Generación de reporte")
print("-" * 70)

reporte_lines = []
reporte_lines.append("=" * 70)
reporte_lines.append("P08 — REPORTE DE VERIFICACIÓN DE COBERTURA PARA FASE 3 SBAS")
reporte_lines.append("=" * 70)
reporte_lines.append(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
reporte_lines.append(f"AOI: {AOI_PATH} (fuente: {aoi_source})")
reporte_lines.append(f"Área AOI: {aoi_area_km2:.2f} km²")
reporte_lines.append(f"Periodo SBAS: {SBAS_START} a {SBAS_END}")
reporte_lines.append(f"Evento: {FECHA_EVENTO}")
reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("RESUMEN DE COBERTURA POR TRACK")
reporte_lines.append("-" * 70)

for track_name in ["Path_40", "Path_120", "Path_18"]:
    cob = cobertura_resultados[track_name]
    track_info = TRACKS[track_name]
    data = resultados_tracks[track_name]
    
    reporte_lines.append(f"\n  {track_name} (Path {track_info['path']}, {track_info['direction']}, Frame {track_info['frame']}):")
    reporte_lines.append(f"    Pares Fase 2: {', '.join(track_info['pares_fase2'])}")
    reporte_lines.append(f"    Cobertura Fase 2: {track_info['cobertura_fase2']}")
    reporte_lines.append(f"    Escenas disponibles (SBAS): {cob['n_escenas']}")
    reporte_lines.append(f"    Cobertura AOI: {cob['cobertura_pct']:.2f}%")
    reporte_lines.append(f"    Área cubierta: {cob['area_cubierta_km2']:.2f} km²")
    reporte_lines.append(f"    Área NO cubierta: {cob['area_no_cubierta_km2']:.2f} km²")
    
    if data['n_escenas'] > 0:
        reporte_lines.append(f"    Pre-evento: {data.get('n_pre', 0)} escenas")
        reporte_lines.append(f"    Post-evento: {data.get('n_post', 0)} escenas")
        reporte_lines.append(f"    Interferogramas estimados (60d): ~{data.get('n_ifg_estimate', 0)}")
    
    apto = 'SI' if cob['cobertura_pct'] >= 95 and cob['n_escenas'] >= 50 else 'NO'
    reporte_lines.append(f"    APTO PARA SBAS: {apto}")

reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("DECISIÓN CIENTÍFICA")
reporte_lines.append("-" * 70)

# Determinar track principal
best_track = max(cobertura_resultados.items(), 
                 key=lambda x: (x[1]['cobertura_pct'], x[1]['n_escenas']))

reporte_lines.append(f"\n  Track principal recomendado: {best_track[0]}")
reporte_lines.append(f"  Justificación:")
reporte_lines.append(f"    1. Mayor cobertura del AOI ({best_track[1]['cobertura_pct']:.1f}%)")
reporte_lines.append(f"    2. Mismo track que pares D1-D4 de Fase 2 (consistencia geométrica)")
reporte_lines.append(f"    3. {best_track[1]['n_escenas']} escenas en rango de estudios publicados")
reporte_lines.append(f"    4. Geometría descendente: sensibilidad a desplazamiento E-W")
reporte_lines.append("")
reporte_lines.append("  Track secundario (validación cross-geometry): Path_120")
reporte_lines.append("  Justificación: geometría ascendente complementaria,")
reporte_lines.append("  permite descomposición del vector de desplazamiento.")
reporte_lines.append("")
reporte_lines.append("-" * 70)
reporte_lines.append("PARÁMETROS RECOMENDADOS PARA SBAS (HyP3 + MintPy)")
reporte_lines.append("-" * 70)
reporte_lines.append("  Umbral temporal máximo: 60 días")
reporte_lines.append("  Umbral perpendicular máximo: 200 m (por defecto HyP3)")
reporte_lines.append("  Coherencia mínima: 0.3 (consistente con Fase 2)")
reporte_lines.append("  Resolución: ~80 m (looks: 20x4, consistente con HyP3 INSAR_GAMMA)")
reporte_lines.append("  Referencia espacial: WGS84 / UTM zona 17N (EPSG:32617)")
reporte_lines.append("  Procesador: HyP3 INSAR_GAMMA → MintPy smallbaselineApp.py")
reporte_lines.append("")
reporte_lines.append("=" * 70)
reporte_lines.append("FIN DEL REPORTE P08")
reporte_lines.append("=" * 70)

reporte_text = "\n".join(reporte_lines)

reporte_path = os.path.join(OUTPUT_DIR, "P08_reporte_fase3.txt")
with open(reporte_path, 'w', encoding='utf-8') as f:
    f.write(reporte_text)
print(f"  Guardado: {reporte_path}")

# Imprimir reporte en consola
print()
print(reporte_text)

# =============================================================================
# RESUMEN FINAL
# =============================================================================

print()
print("=" * 70)
print("P08 COMPLETADO")
print("=" * 70)
print(f"Archivos generados en: {OUTPUT_DIR}")
print()
for f_name in sorted(os.listdir(OUTPUT_DIR)):
    f_path = os.path.join(OUTPUT_DIR, f_name)
    if f_name.startswith("P08"):
        size = os.path.getsize(f_path)
        print(f"  {f_name} ({size/1024:.1f} KB)")
print()
print("SIGUIENTE PASO:")
print("  Si Path 40 confirma cobertura >= 95%, proceder a:")
print("  P09 — Diseño del stack SBAS en HyP3 (envío de jobs)")
print("  Protocolo Maestro Fase 3 se establece tras confirmar cobertura.")