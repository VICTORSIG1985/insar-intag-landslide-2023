# -*- coding: utf-8 -*-
"""
=============================================================================
CONFIGURACIÓN CENTRAL DEL PROYECTO InSAR - INTAG
=============================================================================
Archivo: config.py
Proyecto: Análisis InSAR de deslizamientos en Zona de Intag, Cotacachi, Ecuador
Autor: Víctor [Apellido]
Fecha: 2026-02-25
Versión: 1.0

INSTRUCCIONES:
    1. Editar las variables marcadas con ← EDITAR según tu sistema
    2. Ejecutar este archivo una vez para verificar que las rutas existen
    3. Todos los demás scripts importan este archivo
=============================================================================
"""

import os
from datetime import datetime

# =============================================================================
# 1. RUTAS DEL PROYECTO (← EDITAR SEGÚN TU SISTEMA)
# =============================================================================
# Ruta raíz del proyecto
PROJECT_ROOT = r"D:\POSGRADOS\INTAG"  # ← EDITAR si tu ruta es diferente

# Subrutas (se generan automáticamente)
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, "scripts", "python")
GEE_DIR = os.path.join(PROJECT_ROOT, "scripts", "gee")
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
DATA_S1 = os.path.join(DATA_DIR, "sentinel1")
DATA_S2 = os.path.join(DATA_DIR, "sentinel2")
DATA_DEM = os.path.join(DATA_DIR, "dem")
DATA_CHIRPS = os.path.join(DATA_DIR, "chirps")
DATA_LICSAR = os.path.join(DATA_DIR, "licsar")
DATA_NISAR = os.path.join(DATA_DIR, "nisar")
DATA_INV = os.path.join(DATA_DIR, "inventario")
DATA_INV_GE = os.path.join(DATA_INV, "digitalizacion_GE")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "resultados")
RES_FASE1 = os.path.join(RESULTS_DIR, "fase1_inventario")
RES_FASE2 = os.path.join(RESULTS_DIR, "fase2_dinsar")
RES_FASE3 = os.path.join(RESULTS_DIR, "fase3_sbas")
RES_FASE4 = os.path.join(RESULTS_DIR, "fase4_nisar")
RES_FASE5 = os.path.join(RESULTS_DIR, "fase5_integracion")
RES_FASE6 = os.path.join(RESULTS_DIR, "fase6_publicacion")
FIGS_DIR = os.path.join(PROJECT_ROOT, "figuras")
LOGS_DIR = os.path.join(PROJECT_ROOT, "logs")
DOCS_DIR = os.path.join(PROJECT_ROOT, "docs")

# =============================================================================
# 2. CREDENCIALES NASA EARTHDATA (← EDITAR)
# =============================================================================
# Registrarse en: https://urs.earthdata.nasa.gov/users/new
EARTHDATA_USER = ""   # no subir a GitHub
EARTHDATA_PASS = ""   # no subir a GitHub
# =============================================================================
# 3. ÁREA DE ESTUDIO - ZONA DE INTAG
# =============================================================================
# Bounding box del área de estudio (NO modificar salvo ajuste fino)
STUDY_AREA = {
    "name": "Zona de Intag, Cotacachi, Imbabura, Ecuador",
    "lat_min": 0.20,    # Límite sur
    "lat_max": 0.40,    # Límite norte
    "lon_min": -78.65,  # Límite oeste
    "lon_max": -78.45,  # Límite este
    "center_lat": 0.30,
    "center_lon": -78.55,
}

# WKT del bounding box (para ASF search)
STUDY_AREA_WKT = (
    f"POLYGON(("
    f"{STUDY_AREA['lon_min']} {STUDY_AREA['lat_min']}, "
    f"{STUDY_AREA['lon_max']} {STUDY_AREA['lat_min']}, "
    f"{STUDY_AREA['lon_max']} {STUDY_AREA['lat_max']}, "
    f"{STUDY_AREA['lon_min']} {STUDY_AREA['lat_max']}, "
    f"{STUDY_AREA['lon_min']} {STUDY_AREA['lat_min']}))"
)

# Coordenadas aproximadas por parroquia (centroide)
PARROQUIAS = {
    "Apuela":              {"lat": 0.3517, "lon": -78.5117, "alt_range": "1600-1800"},
    "Peñaherrera":         {"lat": 0.3333, "lon": -78.5500, "alt_range": "1800-2200"},
    "Plaza Gutiérrez":     {"lat": 0.3167, "lon": -78.5333, "alt_range": "1400-1800"},
    "6 de Julio Cuellaje": {"lat": 0.3667, "lon": -78.5667, "alt_range": "1800-2400"},
    "Vacas Galindo":       {"lat": 0.3333, "lon": -78.5167, "alt_range": "1200-1600"},
    "García Moreno":       {"lat": 0.2500, "lon": -78.6167, "alt_range": "800-1400"},
}

# =============================================================================
# 4. PARÁMETROS DEL EVENTO
# =============================================================================
EVENT_DATE = datetime(2023, 12, 19)  # Fecha del deslizamiento principal
EVENT_DATE_STR = "2023-12-19"

# Ventanas temporales para búsqueda de imágenes
PRE_EVENT_START = "2023-10-01"   # Inicio búsqueda pre-evento
PRE_EVENT_END = "2023-12-18"     # Fin búsqueda pre-evento (día antes)
POST_EVENT_START = "2023-12-20"  # Inicio búsqueda post-evento (día después)
POST_EVENT_END = "2024-03-31"    # Fin búsqueda post-evento

# Para series temporales SBAS
TS_START = "2017-01-01"  # Inicio serie temporal
TS_END = "2024-12-31"    # Fin serie temporal

# =============================================================================
# 5. PARÁMETROS DE PROCESAMIENTO
# =============================================================================
# Sentinel-1 InSAR
S1_PARAMS = {
    "platform": "Sentinel-1",
    "beam_mode": "IW",           # Interferometric Wide Swath
    "product_type": "SLC",       # Single Look Complex
    "polarization": "VV",        # Vertical-Vertical (preferida para deformación)
    "flight_direction": None,    # None = ambas; "DESCENDING" o "ASCENDING"
    "max_baseline_perp": 200,    # metros (máximo recomendado para InSAR)
}

# ASF HyP3 InSAR On Demand
HYP3_PARAMS = {
    "dem_name": "copernicus",    # Copernicus GLO-30
    "looks": "20x4",             # Looks (azimuth x range) → ~80m resolución
    "include_dem": True,
    "include_inc_map": True,
    "include_look_vectors": True,
    "apply_water_mask": False,   # No enmascarar agua (hay ríos en la zona)
    "phase_filter_parameter": 0.6,  # Goldstein filter strength
}

# NDVI Change Detection (Sentinel-2)
NDVI_PARAMS = {
    "delta_ndvi_threshold": -0.15,  # Umbral ΔNDVI para candidato a deslizamiento
    "slope_threshold_deg": 15,       # Pendiente mínima (grados)
    "min_area_m2": 100,              # Área mínima de cicatriz (m²)
    "cloud_cover_max": 30,           # Máximo % nubes por imagen S-2
}

# =============================================================================
# 6. PARÁMETROS PARA FIGURAS DEL ARTÍCULO
# =============================================================================
FIG_PARAMS = {
    "dpi": 300,
    "format": "png",  # png para draft, tiff para publicación
    "font_family": "Arial",
    "font_size": 10,
    "cmap_deformation": "RdBu_r",   # Colormap para deformación
    "cmap_coherence": "viridis",     # Colormap para coherencia
    "cmap_velocity": "jet",          # Colormap para velocidad
    "cmap_susceptibility": "RdYlGn_r",  # Colormap para susceptibilidad
}

# =============================================================================
# 7. DATOS OFICIALES DEL EVENTO (SGR Reporte 0729)
# =============================================================================
SGR_DATA = {
    "reporte_numero": "0729",
    "fecha_reporte": "2023-12-31",
    "fuente": "Servicio Nacional de Gestión de Riesgos y Emergencias (SGR), Ecuador",
    "afectaciones": {
        "Plaza Gutiérrez":     {"familias": 112, "personas": 328,   "via_m": 120, "tramos": 4},
        "Apuela":              {"familias": 627, "personas": 1939,  "via_m": 205, "tramos": 6},
        "Peñaherrera":         {"familias": 589, "personas": 1870,  "via_m": 315, "tramos": 5},
        "6 de Julio Cuellaje": {"familias": 23,  "personas": 200,   "via_m": 175, "tramos": 8},
        "Vacas Galindo":       {"familias": 194, "personas": 592,   "via_m": 110, "tramos": 4},
        "García Moreno":       {"familias": 25,  "personas": 100,   "via_m": 400, "tramos": 13},
    },
    "totales": {"familias": 1570, "personas": 5029, "via_m": 1325, "tramos": 40},
}

# =============================================================================
# 8. FUNCIÓN DE INICIALIZACIÓN
# =============================================================================
def crear_estructura():
    """Crea todas las carpetas del proyecto si no existen."""
    carpetas = [
        SCRIPTS_DIR, GEE_DIR, DATA_S1, DATA_S2, DATA_DEM, DATA_CHIRPS,
        DATA_LICSAR, DATA_NISAR, DATA_INV, DATA_INV_GE,
        RES_FASE1, RES_FASE2, RES_FASE3, RES_FASE4, RES_FASE5, RES_FASE6,
        FIGS_DIR, LOGS_DIR, DOCS_DIR,
    ]
    creadas = 0
    for carpeta in carpetas:
        if not os.path.exists(carpeta):
            os.makedirs(carpeta, exist_ok=True)
            creadas += 1
    return creadas

def verificar_config():
    """Verifica que la configuración es válida."""
    print("=" * 70)
    print("VERIFICACIÓN DE CONFIGURACIÓN - Proyecto InSAR Intag")
    print("=" * 70)
    
    # Verificar ruta raíz
    if os.path.exists(PROJECT_ROOT):
        print(f"[OK] Ruta del proyecto existe: {PROJECT_ROOT}")
    else:
        print(f"[!!] Ruta del proyecto NO existe: {PROJECT_ROOT}")
        print(f"     Creando estructura...")
    
    # Crear carpetas
    n = crear_estructura()
    print(f"[OK] Estructura de carpetas verificada ({n} nuevas creadas)")
    
    # Verificar credenciales
    if EARTHDATA_USER == "TU_USUARIO_AQUI":
        print("[!!] CREDENCIALES: Debes editar EARTHDATA_USER y EARTHDATA_PASS en config.py")
    else:
        print(f"[OK] Credenciales Earthdata configuradas (usuario: {EARTHDATA_USER})")
    
    # Mostrar área de estudio
    print(f"\n[INFO] Área de estudio: {STUDY_AREA['name']}")
    print(f"       Bounding box: {STUDY_AREA['lat_min']}°N - {STUDY_AREA['lat_max']}°N / "
          f"{STUDY_AREA['lon_min']}°W - {STUDY_AREA['lon_max']}°W")
    print(f"       Fecha del evento: {EVENT_DATE_STR}")
    
    print("\n" + "=" * 70)
    print("Configuración verificada. Proceder con P00_verificar_entorno.py")
    print("=" * 70)


# Ejecutar verificación al importar directamente
if __name__ == "__main__":
    verificar_config()