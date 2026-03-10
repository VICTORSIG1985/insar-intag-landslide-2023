# -*- coding: utf-8 -*-
"""
=============================================================================
P00 - VERIFICACIÓN DEL ENTORNO DE TRABAJO
=============================================================================
Proyecto: InSAR Deslizamientos Zona de Intag
Fase: Preparación
Autor: Víctor [Apellido]
Fecha: 2026-02-25

PROPÓSITO:
    Verificar que todas las librerías necesarias están instaladas
    y que la configuración es correcta antes de empezar el procesamiento.

EJECUTAR EN: Spyder (Run file F5)
=============================================================================
"""

import sys
import importlib
from datetime import datetime

def check_library(name, import_name=None, pip_name=None):
    """Verifica si una librería está instalada y retorna versión."""
    if import_name is None:
        import_name = name
    if pip_name is None:
        pip_name = name
    try:
        mod = importlib.import_module(import_name)
        version = getattr(mod, "__version__", "sin versión")
        return True, version
    except ImportError:
        return False, pip_name

# =============================================================================
# LISTA DE LIBRERÍAS REQUERIDAS
# =============================================================================
print("=" * 70)
print(f"VERIFICACIÓN DEL ENTORNO - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Python: {sys.version}")
print(f"Ejecutando desde: {sys.executable}")
print("=" * 70)

libraries = [
    # (Nombre display, nombre import, nombre pip)
    ("NumPy", "numpy", "numpy"),
    ("Pandas", "pandas", "pandas"),
    ("Matplotlib", "matplotlib", "matplotlib"),
    ("SciPy", "scipy", "scipy"),
    ("Scikit-learn", "sklearn", "scikit-learn"),
    ("Rasterio", "rasterio", "rasterio"),
    ("GeoPandas", "geopandas", "geopandas"),
    ("Shapely", "shapely", "shapely"),
    ("Fiona", "fiona", "fiona"),
    ("ASF Search", "asf_search", "asf_search"),
    ("HyP3 SDK", "hyp3_sdk", "hyp3_sdk"),
    ("xarray", "xarray", "xarray"),
    ("netCDF4", "netCDF4", "netCDF4"),
    ("h5py", "h5py", "h5py"),
    ("RichDEM", "richdem", "richdem"),
]

optional_libraries = [
    ("LiCSBAS", "LiCSBAS", "licsbas"),
    ("MintPy", "mintpy", "mintpy"),
    ("GDAL", "osgeo.gdal", "gdal"),
]

print("\n--- LIBRERÍAS REQUERIDAS ---")
missing_required = []
for display_name, import_name, pip_name in libraries:
    ok, info = check_library(display_name, import_name, pip_name)
    if ok:
        print(f"  [OK]  {display_name:20s} v{info}")
    else:
        print(f"  [!!]  {display_name:20s} NO INSTALADA → pip install {info}")
        missing_required.append(pip_name)

print("\n--- LIBRERÍAS OPCIONALES ---")
for display_name, import_name, pip_name in optional_libraries:
    ok, info = check_library(display_name, import_name, pip_name)
    if ok:
        print(f"  [OK]  {display_name:20s} v{info}")
    else:
        print(f"  [--]  {display_name:20s} No instalada (opcional) → pip install {info}")

# =============================================================================
# VERIFICAR CONFIGURACIÓN
# =============================================================================
print("\n--- CONFIGURACIÓN DEL PROYECTO ---")
try:
    # Ajustar path para importar config
    import os
    config_dir = os.path.dirname(os.path.abspath(__file__))
    if config_dir not in sys.path:
        sys.path.insert(0, config_dir)
    
    import config
    config.verificar_config()
    print("  [OK] config.py importado correctamente")
except Exception as e:
    print(f"  [!!] Error importando config.py: {e}")
    print(f"       Asegúrate de que config.py está en la misma carpeta")

# =============================================================================
# RESUMEN
# =============================================================================
print("\n" + "=" * 70)
if missing_required:
    print("ACCIÓN REQUERIDA: Instalar librerías faltantes ejecutando:")
    print(f"  pip install {' '.join(missing_required)}")
    print("\nEn Anaconda Prompt o terminal de Spyder (Tools > Open Command Prompt)")
else:
    print("TODO CORRECTO. El entorno está listo para trabajar.")
    print("Siguiente paso: Ejecutar P01_buscar_sentinel1.py")
print("=" * 70)
