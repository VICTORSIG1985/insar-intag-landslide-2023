# ============================================================
# P04b_v2 - PATH 18: SELECCIÓN, ENVÍO, ESPERA Y DESCARGA
# ============================================================
# Script completo que hace TODO de una vez:
#   1. Busca granules Path 18 en ASF
#   2. Selecciona 3 pares alternativos (evita 2023-12-06)
#   3. Envía a HyP3
#   4. Espera hasta que terminen
#   5. Descarga los productos
#
# P18_A3 (2023-12-18 → 2024-01-11) ya fue exitoso y descargado.
# Este script procesa los 3 pares restantes.
#
# Pares alternativos:
#   P18_B1: 2023-11-24 → 2023-12-18 (24d) pre→post evento
#   P18_B2: 2023-12-18 → 2023-12-30 (12d) post evento corto
#   P18_B3: 2023-12-30 → 2024-01-23 (24d) post extendido
#
# hyp3_sdk v7.7.6: submit_insar_job() retorna Batch
# ============================================================
import json, os, sys, time
from datetime import datetime

print("=" * 70)
print("P04b_v2 - PATH 18 COMPLETO (selección → descarga)")
print("=" * 70)

# --- Librerías ---
try:
    import asf_search as asf
    import hyp3_sdk
    import geopandas as gpd
    from shapely.geometry import box
    from shapely.ops import unary_union
    print(f"  asf_search v{asf.__version__}")
    print(f"  hyp3_sdk v{hyp3_sdk.__version__}")
except ImportError as e:
    print(f"[ERROR] Librería faltante: {e}")
    sys.exit(1)

# ============================================================
# 1. CONFIGURACIÓN
# ============================================================
SAVE_DIR     = r"D:\POSGRADOS\INTAG\data\sentinel1"
JOBS_DIR     = os.path.join(SAVE_DIR, "hyp3_jobs")
PRODUCTS_DIR = os.path.join(SAVE_DIR, "hyp3_products")
AOI_PATH     = r"D:\POSGRADOS\INTAG\INTAJ.gpkg"
os.makedirs(JOBS_DIR, exist_ok=True)
os.makedirs(PRODUCTS_DIR, exist_ok=True)

PROJECT_NAME = "Intag_DInSAR_2023_P18v2"

# Parámetros idénticos a P04 (reproducibilidad)
HYP3_PARAMS = {
    "looks":                     "20x4",
    "include_dem":               True,
    "include_inc_map":           True,
    "include_look_vectors":      True,
    "include_displacement_maps": True,
    "include_wrapped_phase":     True,
    "apply_water_mask":          False,
    "phase_filter_parameter":    0.6,
}

# Pares alternativos (evitan granule 2023-12-06 que causó fallos)
PARES = [
    ("P18_B1", "2023-11-24", "2023-12-18", 24, "Pre→Post evento"),
    ("P18_B2", "2023-12-18", "2023-12-30", 12, "Post evento (12d)"),
    ("P18_B3", "2023-12-30", "2024-01-23", 24, "Post extendido"),
]

# Intervalo de consulta (segundos)
POLL_INTERVAL = 120

print(f"\n  Proyecto: {PROJECT_NAME}")
print(f"  Pares a procesar: {len(PARES)}")
print(f"  (P18_A3 ya exitoso, no se reenvía)")

# ============================================================
# 2. BUSCAR GRANULES PATH 18
# ============================================================
print(f"\n{'─' * 70}")
print("PASO 1/5: BUSCANDO GRANULES PATH 18 EN ASF")
print(f"{'─' * 70}")

aoi = gpd.read_file(AOI_PATH)
if aoi.crs and aoi.crs.to_epsg() != 4326:
    aoi = aoi.to_crs(epsg=4326)
bbox_wkt = box(*unary_union(aoi.geometry).bounds).wkt

try:
    results = asf.search(
        intersectsWith=bbox_wkt,
        dataset=asf.DATASET.SENTINEL1,
        processingLevel=[asf.PRODUCT_TYPE.SLC],
        flightDirection=asf.FLIGHT_DIRECTION.ASCENDING,
        start="2023-11-01",
        end="2024-02-28",
    )
except Exception:
    results = asf.geo_search(
        intersectsWith=bbox_wkt,
        platform=[asf.PLATFORM.SENTINEL1],
        processingLevel=[asf.PRODUCT_TYPE.SLC],
        flightDirection=asf.FLIGHT_DIRECTION.ASCENDING,
        start="2023-11-01",
        end="2024-02-28",
    )

# Filtrar Path 18
escenas = {}
for r in results:
    props = r.properties
    if props.get("pathNumber") == 18:
        fecha = str(props.get("startTime", ""))[:10]
        escenas[fecha] = props.get("sceneName", "?")

print(f"  Escenas Path 18: {len(escenas)}")
for f in sorted(escenas.keys()):
    marca = " ← EVITAR" if f == "2023-12-06" else ""
    print(f"    {f}: {escenas[f]}{marca}")

# ============================================================
# 3. CONSTRUIR PARES
# ============================================================
print(f"\n{'─' * 70}")
print("PASO 2/5: CONSTRUYENDO PARES")
print(f"{'─' * 70}")

pares_listos = []

for par_id, pre_date, post_date, baseline, proposito in PARES:
    pre_g = escenas.get(pre_date)
    post_g = escenas.get(post_date)
    
    if pre_g and post_g:
        pares_listos.append({
            "id": par_id,
            "pre_scene": pre_g,
            "post_scene": post_g,
            "pre_date": pre_date,
            "post_date": post_date,
            "baseline_dias": baseline,
            "flightDirection": "ASCENDING",
            "pathNumber": 18,
            "frameNumber": 1180,
            "calidad": "★★★" if baseline <= 24 else "★★",
            "proposito": proposito,
        })
        print(f"  ✓ {par_id}: {pre_date} → {post_date} ({baseline}d)")
    else:
        print(f"  ✗ {par_id}: granule no encontrado (pre={pre_g}, post={post_g})")

if not pares_listos:
    print("\n  [ERROR] No se pudo construir ningún par.")
    sys.exit(1)

print(f"\n  Pares listos: {len(pares_listos)}/{len(PARES)}")

# ============================================================
# 4. CONECTAR Y ENVIAR A HyP3
# ============================================================
print(f"\n{'─' * 70}")
print("PASO 3/5: ENVIANDO A HyP3")
print(f"{'─' * 70}")

try:
    hyp3 = hyp3_sdk.HyP3()
    credits_info = hyp3.check_credits()
    print(f"  ✓ Conectado. Créditos: {credits_info}")
except Exception as e:
    print(f"  ✗ Error conexión: {e}")
    sys.exit(1)

jobs_enviados = []

for i, par in enumerate(pares_listos, 1):
    job_name = f"{PROJECT_NAME}_{par['id']}_{par['pre_date']}_{par['post_date']}"
    print(f"\n  [{i}/{len(pares_listos)}] {par['id']}...")
    
    try:
        batch = hyp3.submit_insar_job(
            granule1=par["pre_scene"],
            granule2=par["post_scene"],
            name=job_name,
            **HYP3_PARAMS,
        )
        job = batch.jobs[0]
        
        par["job_id"] = job.job_id
        par["status"] = job.status_code
        par["job_name"] = job_name
        par["submitted_at"] = datetime.now().isoformat()
        jobs_enviados.append(par)
        
        print(f"    ✓ Job ID: {job.job_id} → {job.status_code}")
    except Exception as e:
        print(f"    ✗ Error: {e}")
        par["job_id"] = None
        par["error"] = str(e)
        par["submitted_at"] = datetime.now().isoformat()
        jobs_enviados.append(par)

exitosos = [j for j in jobs_enviados if j.get("job_id")]
fallidos_envio = [j for j in jobs_enviados if not j.get("job_id")]

print(f"\n  Enviados: {len(exitosos)}/{len(pares_listos)}")
if fallidos_envio:
    print(f"  ⚠ Fallidos en envío: {len(fallidos_envio)}")
    for j in fallidos_envio:
        print(f"    - {j['id']}: {j.get('error','?')}")

if not exitosos:
    print("\n  [ERROR] Ningún job se envió. Abortando.")
    sys.exit(1)

# Guardar metadatos
meta = {
    "generado": datetime.now().isoformat(),
    "proyecto": PROJECT_NAME,
    "procesador": "INSAR_GAMMA",
    "dem": "Copernicus GLO-30 (automático)",
    "track": "Path 18 / Frame 1180 (ASCENDING)",
    "nota": "Pares alternativos. Evita granule 2023-12-06 (causó fallos en v1). "
            "P18_A3 ya exitoso, no incluido aquí.",
    "parametros": HYP3_PARAMS,
    "hyp3_sdk_version": hyp3_sdk.__version__,
    "jobs": jobs_enviados,
}
meta_path = os.path.join(JOBS_DIR, "hyp3_jobs_path18v2_enviados.json")
with open(meta_path, "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2, ensure_ascii=False)
print(f"  ✓ Metadatos: {meta_path}")

# ============================================================
# 5. ESPERAR HASTA QUE TERMINEN
# ============================================================
print(f"\n{'─' * 70}")
print("PASO 4/5: ESPERANDO PROCESAMIENTO")
print(f"{'─' * 70}")

job_ids = {j["id"]: j["job_id"] for j in exitosos}
total = len(job_ids)

while True:
    # Consultar estado
    estados = {}
    for par_id, jid in job_ids.items():
        try:
            job = hyp3.get_job_by_id(jid)
            estados[par_id] = {
                "status": job.status_code,
                "job": job,
            }
        except Exception as e:
            estados[par_id] = {"status": "ERROR", "error": str(e)}
    
    # Mostrar
    ts = datetime.now().strftime("%H:%M:%S")
    simbolos = {"SUCCEEDED": "✓", "FAILED": "✗", "RUNNING": "⟳", "PENDING": "◌"}
    
    n_done = sum(1 for v in estados.values() if v["status"] in ("SUCCEEDED", "FAILED"))
    n_ok = sum(1 for v in estados.values() if v["status"] == "SUCCEEDED")
    n_fail = sum(1 for v in estados.values() if v["status"] == "FAILED")
    
    print(f"\n  [{ts}] ", end="")
    for par_id, info in estados.items():
        s = simbolos.get(info["status"], "?")
        print(f"{s}{par_id} ", end="")
    print(f"  ({n_ok}✓ {n_fail}✗ {total - n_done}⏳)")
    
    # ¿Todos terminaron?
    if n_done >= total:
        break
    
    # Esperar
    print(f"  Próxima consulta en {POLL_INTERVAL}s...", end="", flush=True)
    time.sleep(POLL_INTERVAL)
    print(" consultando...")

# ============================================================
# 6. DESCARGAR PRODUCTOS
# ============================================================
print(f"\n{'─' * 70}")
print("PASO 5/5: DESCARGANDO PRODUCTOS")
print(f"{'─' * 70}")

descargados = []
fallidos_proc = []

for par_id, info in estados.items():
    if info["status"] == "SUCCEEDED":
        print(f"\n  Descargando {par_id}...")
        try:
            job = info["job"]
            archivos = job.download_files(PRODUCTS_DIR)
            print(f"    ✓ {len(archivos)} archivo(s)")
            for a in archivos:
                print(f"      → {os.path.basename(str(a))}")
            descargados.append(par_id)
        except Exception as e:
            print(f"    ✗ Error descarga: {e}")
            fallidos_proc.append(par_id)
    elif info["status"] == "FAILED":
        fallidos_proc.append(par_id)
        print(f"\n  ✗ {par_id}: FAILED en HyP3")

# ============================================================
# 7. GUARDAR REPORTE FINAL
# ============================================================
reporte = {
    "generado": datetime.now().isoformat(),
    "exitosos": descargados,
    "fallidos": fallidos_proc,
    "estados": {k: v["status"] for k, v in estados.items()},
    "nota": "P18_A3 procesado previamente (exitoso). "
            f"Este script procesó {len(job_ids)} pares alternativos.",
}
reporte_path = os.path.join(JOBS_DIR, "hyp3_jobs_path18v2_final.json")
with open(reporte_path, "w", encoding="utf-8") as f:
    json.dump(reporte, f, indent=2, ensure_ascii=False)

# ============================================================
# 8. RESUMEN FINAL
# ============================================================
print(f"\n{'=' * 70}")
print("RESUMEN FINAL P04b_v2")
print(f"{'=' * 70}")
print(f"  Pares enviados:    {len(exitosos)}")
print(f"  Descargados:       {len(descargados)}")
print(f"  Fallidos:          {len(fallidos_proc)}")

if fallidos_proc:
    print(f"\n  ⚠ PARES FALLIDOS:")
    for fp in fallidos_proc:
        print(f"    - {fp}")

print(f"\n  INVENTARIO PATH 18 COMPLETO:")
print(f"    P18_A3 (12-18→01-11, 24d): ✓ ya descargado previamente")
for par_id in job_ids:
    estado = estados[par_id]["status"]
    s = "✓" if estado == "SUCCEEDED" else "✗"
    par = next(p for p in pares_listos if p["id"] == par_id)
    print(f"    {par_id} ({par['pre_date']}→{par['post_date']}, "
          f"{par['baseline_dias']}d): {s} {estado}")

print(f"\n  INVENTARIO COMPLETO DEL PROYECTO:")
print(f"  {'─' * 50}")
print(f"  Descending Path 40  (D1-D4):       ✓ 4 pares")
print(f"  Ascending  Path 120 (A1-A4):       ✓ 4 pares")
print(f"  Ascending  Path 18  (P18_A3+B1-3): ", end="")

total_p18 = 1 + len(descargados)  # P18_A3 + nuevos descargados
if total_p18 >= 4:
    print(f"✓ {total_p18} pares")
    print(f"\n  ¡TODOS LOS DATOS COMPLETOS!")
    print(f"  Siguiente paso: P06 DEFINITIVO")
else:
    print(f"◐ {total_p18} pares ({4 - total_p18} faltantes)")

print(f"\n  Productos en: {PRODUCTS_DIR}")
print(f"  Reportes en:  {JOBS_DIR}")
print(f"{'=' * 70}")