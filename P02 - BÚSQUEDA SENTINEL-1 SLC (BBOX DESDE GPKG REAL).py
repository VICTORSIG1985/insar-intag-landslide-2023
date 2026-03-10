# ============================================================
# P02 - BÚSQUEDA SENTINEL-1 SLC (BBOX DESDE GPKG REAL)
# ============================================================
# Lee el extent real de las 6 parroquias afectadas desde
# el GeoPackage exportado de QGIS y busca en ASF DAAC.
# ============================================================
import asf_search as asf
import json, os, csv
import geopandas as gpd
from datetime import datetime
from shapely.ops import unary_union

GPKG_PATH = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"
SAVE_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1"
os.makedirs(SAVE_DIR, exist_ok=True)

print("=" * 70)
print("P02 - BÚSQUEDA SENTINEL-1 SLC (BBOX CORREGIDO)")
print("=" * 70)

# ============================================================
# 1. LEER EXTENT REAL DESDE GPKG
# ============================================================
print("\n--- LEYENDO ÁREA DE ESTUDIO DESDE GPKG ---")
gdf = gpd.read_file(GPKG_PATH)
print(f"  Archivo: {GPKG_PATH}")
print(f"  Parroquias cargadas: {len(gdf)}")
for _, row in gdf.iterrows():
    nombre = row.get("DPA_DESPAR", row.get("dpa_despar", "?"))
    codigo = row.get("DPA_PARROQ", row.get("dpa_parroq", "?"))
    print(f"    - {nombre} (código: {codigo})")

# Obtener bbox real
bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
lon_min, lat_min, lon_max, lat_max = bounds

# Agregar buffer de 0.02° (~2.2 km) para capturar bordes de deslizamientos
BUFFER = 0.02
lon_min -= BUFFER
lat_min -= BUFFER
lon_max += BUFFER
lat_max += BUFFER

print(f"\n  Extent real de parroquias:")
print(f"    Lon: {bounds[0]:.4f} a {bounds[2]:.4f}")
print(f"    Lat: {bounds[1]:.4f} a {bounds[3]:.4f}")
print(f"  Extent con buffer ({BUFFER}°):")
print(f"    Lon: {lon_min:.4f} a {lon_max:.4f}")
print(f"    Lat: {lat_min:.4f} a {lat_max:.4f}")

# Generar WKT
BBOX_WKT = (f"POLYGON(({lon_min} {lat_min}, {lon_max} {lat_min}, "
            f"{lon_max} {lat_max}, {lon_min} {lat_max}, {lon_min} {lat_min}))")
print(f"\n  WKT: {BBOX_WKT}")

# Guardar bbox para referencia
bbox_info = {
    "generado": datetime.now().isoformat(),
    "fuente": GPKG_PATH,
    "parroquias": len(gdf),
    "extent_original": {"lon_min": bounds[0], "lat_min": bounds[1], "lon_max": bounds[2], "lat_max": bounds[3]},
    "buffer_grados": BUFFER,
    "extent_con_buffer": {"lon_min": lon_min, "lat_min": lat_min, "lon_max": lon_max, "lat_max": lat_max},
    "wkt": BBOX_WKT,
}
bbox_path = os.path.join(SAVE_DIR, "bbox_area_estudio.json")
with open(bbox_path, "w", encoding="utf-8") as f:
    json.dump(bbox_info, f, indent=2, ensure_ascii=False)
print(f"  ✓ Guardado: {bbox_path}")

# ============================================================
# 2. BÚSQUEDA EN ASF
# ============================================================
# --- PRE-EVENTO ---
print("\n--- BÚSQUEDA PRE-EVENTO ---")
pre = asf.search(
    platform=asf.PLATFORM.SENTINEL1,
    processingLevel=asf.PRODUCT_TYPE.SLC,
    beamMode=asf.BEAMMODE.IW,
    intersectsWith=BBOX_WKT,
    start="2023-10-01",
    end="2023-12-18",
)
print(f"  PRE-EVENTO: {len(pre)} escenas")

# --- POST-EVENTO ---
print("\n--- BÚSQUEDA POST-EVENTO ---")
post = asf.search(
    platform=asf.PLATFORM.SENTINEL1,
    processingLevel=asf.PRODUCT_TYPE.SLC,
    beamMode=asf.BEAMMODE.IW,
    intersectsWith=BBOX_WKT,
    start="2023-12-20",
    end="2024-03-31",
)
print(f"  POST-EVENTO: {len(post)} escenas")

# --- SERIE TEMPORAL ---
print("\n--- BÚSQUEDA SERIE TEMPORAL ---")
ts = asf.search(
    platform=asf.PLATFORM.SENTINEL1,
    processingLevel=asf.PRODUCT_TYPE.SLC,
    beamMode=asf.BEAMMODE.IW,
    intersectsWith=BBOX_WKT,
    start="2017-01-01",
    end="2024-12-31",
)
print(f"  SERIE TEMPORAL: {len(ts)} escenas")

# ============================================================
# 3. GUARDAR RESULTADOS
# ============================================================
def escenas_a_lista(resultados):
    lista = []
    for r in resultados:
        p = r.properties
        lista.append({
            "sceneName": p.get("sceneName", ""),
            "date": p.get("startTime", "")[:10],
            "flightDirection": p.get("flightDirection", ""),
            "pathNumber": p.get("pathNumber"),
            "frameNumber": p.get("frameNumber"),
            "polarization": p.get("polarization", ""),
            "url": p.get("url", ""),
            "bytes": p.get("bytes", 0),
            "sizeMB": round(p.get("bytes", 0) / 1e6, 1),
        })
    return lista

for nombre, datos in [("pre_evento", pre), ("post_evento", post), ("serie_temporal", ts)]:
    catalogo = {
        "generado": datetime.now().isoformat(),
        "busqueda": nombre,
        "bbox_wkt": BBOX_WKT,
        "bbox": {"lon_min": lon_min, "lat_min": lat_min, "lon_max": lon_max, "lat_max": lat_max},
        "fuente_bbox": "GeoPackage 6 parroquias afectadas + buffer 0.02°",
        "total_escenas": len(datos),
        "escenas": escenas_a_lista(datos),
    }
    ruta = os.path.join(SAVE_DIR, f"catalogo_{nombre}.json")
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(catalogo, f, indent=2, ensure_ascii=False)
    print(f"  ✓ Guardado: {ruta}")

# CSV consolidado
csv_path = os.path.join(SAVE_DIR, "resumen_todas_escenas.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "grupo", "sceneName", "date", "flightDirection",
        "pathNumber", "frameNumber", "polarization", "sizeMB", "url"
    ])
    writer.writeheader()
    for grupo, datos in [("pre_evento", pre), ("post_evento", post), ("serie_temporal", ts)]:
        for r in datos:
            p = r.properties
            writer.writerow({
                "grupo": grupo, "sceneName": p.get("sceneName", ""),
                "date": p.get("startTime", "")[:10],
                "flightDirection": p.get("flightDirection", ""),
                "pathNumber": p.get("pathNumber"), "frameNumber": p.get("frameNumber"),
                "polarization": p.get("polarization", ""),
                "sizeMB": round(p.get("bytes", 0) / 1e6, 1), "url": p.get("url", ""),
            })
print(f"  ✓ Guardado: {csv_path}")

# --- RESUMEN ---
asc = [r for r in ts if r.properties['flightDirection'] == 'ASCENDING']
desc = [r for r in ts if r.properties['flightDirection'] == 'DESCENDING']

print("\n" + "=" * 70)
print("RESUMEN P02")
print("=" * 70)
print(f"  Área: 6 parroquias desde {GPKG_PATH}")
print(f"  BBox: [{lon_min:.4f}, {lat_min:.4f}, {lon_max:.4f}, {lat_max:.4f}]")
print(f"  Pre-evento:  {len(pre)} escenas")
print(f"  Post-evento: {len(post)} escenas")
print(f"  Serie total: {len(ts)} ({len(asc)} ASC + {len(desc)} DESC)")
print(f"\n  Archivos en: {SAVE_DIR}")
print("=" * 70)