import sys
sys.path.insert(0, r"D:\POSGRADOS\INTAG\INTAG_PROJECT\scripts\python")
from config import *

print("=" * 70)
print("VERIFICACIÓN DEL ENTORNO")
print("=" * 70)

# Verificar librerías
libs = {
    "numpy": "NumPy", "pandas": "Pandas", "matplotlib": "Matplotlib",
    "scipy": "SciPy", "sklearn": "Scikit-learn", "rasterio": "Rasterio",
    "geopandas": "GeoPandas", "shapely": "Shapely", "fiona": "Fiona",
    "asf_search": "ASF Search", "hyp3_sdk": "HyP3 SDK",
    "xarray": "xarray", "netCDF4": "netCDF4", "h5py": "h5py",
    "osgeo.gdal": "GDAL",
}

for mod, nombre in libs.items():
    try:
        m = __import__(mod)
        v = getattr(m, '__version__', 'OK')
        print(f"  [OK]  {nombre:20s} v{v}")
    except:
        print(f"  [!!]  {nombre:20s} NO INSTALADA")

# Verificar config
print(f"\n[OK] Credenciales: {EARTHDATA_USER}")
print(f"[OK] Área: {STUDY_AREA['name']}")
print(f"[OK] Evento: {EVENT_DATE_STR}")
print(f"[OK] Ruta proyecto: {PROJECT_ROOT}")
print("=" * 70)