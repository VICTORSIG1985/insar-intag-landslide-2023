# ============================================================
# P04 - ENVÍO DE 8 PARES A HyP3 PARA PROCESAMIENTO InSAR
# ============================================================
# Lee los 8 pares de pares_dinsar_optimos.json y los envía
# al servicio HyP3 (ASF) para procesamiento interferométrico
# usando el procesador INSAR_GAMMA.
#
# DEM: Copernicus GLO-30 (usado automáticamente por HyP3)
# Procesador: GAMMA
#
# Parámetros válidos para submit_insar_job (hyp3_sdk v7.7.6):
#   - looks: "20x4" o "10x2"
#   - include_look_vectors: bool
#   - include_dem: bool
#   - include_inc_map: bool
#   - include_displacement_maps: bool
#   - include_wrapped_phase: bool
#   - apply_water_mask: bool
#   - phase_filter_parameter: float (0.0 - 1.0)
#
# NOTA IMPORTANTE (v7.7.6):
#   submit_insar_job() retorna un objeto Batch, NO un Job.
#   Para obtener el Job individual: batch.jobs[0]
#
# Ref: https://hyp3-docs.asf.alaska.edu/using/sdk_api/
#      https://github.com/ASFHyP3/hyp3-sdk/blob/develop/CHANGELOG.md
#
# Guardado: D:\POSGRADOS\INTAG\data\sentinel1\hyp3_jobs\
# ============================================================
import json, os, sys
from datetime import datetime

try:
    import hyp3_sdk
except ImportError:
    print("[ERROR] hyp3_sdk no instalado. Ejecutar:")
    print("  pip install hyp3_sdk")
    sys.exit(1)

print(f"  hyp3_sdk versión: {hyp3_sdk.__version__}")

# ============================================================
# 1. CONFIGURACIÓN
# ============================================================
SAVE_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1"
JOBS_DIR = os.path.join(SAVE_DIR, "hyp3_jobs")
os.makedirs(JOBS_DIR, exist_ok=True)

# Parámetros VÁLIDOS para INSAR_GAMMA (hyp3_sdk v7.7.6)
# NOTA: dem_name NO existe para InSAR (solo para RTC).
#       HyP3 usa Copernicus GLO-30 automáticamente.
HYP3_PARAMS = {
    "looks":                     "20x4",   # 20 azimuth × 4 range → ~80m
    "include_dem":               True,      # DEM recortado al área
    "include_inc_map":           True,      # Mapa de ángulo de incidencia
    "include_look_vectors":      True,      # Vectores de look (descomp. 2D)
    "include_displacement_maps": True,      # Desplazamiento LOS + vertical
    "include_wrapped_phase":     True,      # Fase envuelta (wrapped)
    "apply_water_mask":          False,     # No enmascarar agua (ríos)
    "phase_filter_parameter":    0.6,       # Goldstein filter strength
}

PROJECT_NAME = "Intag_DInSAR_2023"

print("=" * 70)
print("P04 - ENVÍO DE PARES A HyP3")
print("=" * 70)
print(f"\n  Proyecto HyP3: {PROJECT_NAME}")
print(f"  DEM: Copernicus GLO-30 (automático en INSAR_GAMMA)")
print(f"  Parámetros de procesamiento:")
for k, v in HYP3_PARAMS.items():
    print(f"    {k}: {v}")

# ============================================================
# 2. CARGAR PARES SELECCIONADOS
# ============================================================
pares_path = os.path.join(SAVE_DIR, "pares_dinsar_optimos.json")
with open(pares_path, encoding="utf-8") as f:
    data = json.load(f)

pares = data["pares"]
print(f"\n  Pares a procesar: {len(pares)}")
for p in pares:
    print(f"    {p['id']} {p['calidad']} {p['flightDirection']:10s} "
          f"{p['pre_date']} → {p['post_date']} ({p['baseline_dias']}d)")

# ============================================================
# 3. CONECTAR A HyP3
# ============================================================
print(f"\n--- CONECTANDO A HyP3 ---")
try:
    hyp3 = hyp3_sdk.HyP3()
    print(f"  ✓ Conectado a HyP3 (autenticación via .netrc)")
except Exception as e:
    print(f"  ✗ Error de conexión: {e}")
    print(f"\n  Verificar C:\\Users\\Usuario\\_netrc:")
    print(f"    machine urs.earthdata.nasa.gov")
    print(f"    login TU_USUARIO")
    print(f"    password TU_PASSWORD")
    sys.exit(1)

# Verificar créditos (método actual, no deprecado)
try:
    credits_info = hyp3.check_credits()
    print(f"  Créditos disponibles: {credits_info}")
    if isinstance(credits_info, (int, float)) and credits_info < len(pares):
        print(f"  ⚠ Créditos insuficientes para {len(pares)} pares.")
        sys.exit(1)
except Exception as e:
    print(f"  ⚠ No se pudo verificar créditos: {e}")
    print(f"    Continuando de todos modos...")

# ============================================================
# 4. ENVIAR JOBS
# ============================================================
# CORRECCIÓN v7.7.6: submit_insar_job() retorna un objeto Batch,
# NO un Job individual. Se debe extraer el Job con batch.jobs[0].
# Ref: https://github.com/ASFHyP3/hyp3-sdk (SDK example notebook)
# ============================================================
print(f"\n--- ENVIANDO {len(pares)} PARES A HyP3 ---\n")

jobs_enviados = []

for i, par in enumerate(pares, 1):
    job_name = f"{PROJECT_NAME}_{par['id']}_{par['pre_date']}_{par['post_date']}"
    
    print(f"  [{i}/{len(pares)}] Enviando {par['id']}...")
    print(f"    PRE:  {par['pre_scene']}")
    print(f"    POST: {par['post_scene']}")
    
    try:
        # submit_insar_job retorna Batch (v7.7.6), no Job
        batch = hyp3.submit_insar_job(
            granule1=par["pre_scene"],
            granule2=par["post_scene"],
            name=job_name,
            **HYP3_PARAMS,
        )
        
        # Extraer el Job individual del Batch (1 par = 1 job)
        job = batch.jobs[0]
        
        job_info = {
            "par_id":          par["id"],
            "job_id":          job.job_id,
            "job_name":        job_name,
            "status":          job.status_code,
            "pre_scene":       par["pre_scene"],
            "post_scene":      par["post_scene"],
            "pre_date":        par["pre_date"],
            "post_date":       par["post_date"],
            "baseline_dias":   par["baseline_dias"],
            "flightDirection": par["flightDirection"],
            "calidad":         par["calidad"],
            "cobertura":       par["cobertura"],
            "proposito":       par["proposito"],
            "submitted_at":    datetime.now().isoformat(),
            "hyp3_params":     HYP3_PARAMS,
        }
        jobs_enviados.append(job_info)
        
        print(f"    ✓ Enviado! Job ID: {job.job_id}")
        print(f"    Estado: {job.status_code}\n")
        
    except Exception as e:
        print(f"    ✗ Error: {e}\n")
        jobs_enviados.append({
            "par_id":       par["id"],
            "job_id":       None,
            "error":        str(e),
            "pre_scene":    par["pre_scene"],
            "post_scene":   par["post_scene"],
            "submitted_at": datetime.now().isoformat(),
        })

# ============================================================
# 5. GUARDAR METADATOS DE JOBS
# ============================================================
jobs_exitosos = [j for j in jobs_enviados if j.get("job_id")]
jobs_fallidos = [j for j in jobs_enviados if not j.get("job_id")]

output = {
    "generado":        datetime.now().isoformat(),
    "proyecto_hyp3":   PROJECT_NAME,
    "procesador":      "INSAR_GAMMA",
    "dem":             "Copernicus GLO-30 (automático)",
    "hyp3_sdk_version": hyp3_sdk.__version__,
    "parametros":      HYP3_PARAMS,
    "total_enviados":  len(jobs_exitosos),
    "total_fallidos":  len(jobs_fallidos),
    "jobs":            jobs_enviados,
}

jobs_path = os.path.join(JOBS_DIR, "hyp3_jobs_enviados.json")
with open(jobs_path, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"  ✓ Guardado: {jobs_path}")

ids_path = os.path.join(JOBS_DIR, "job_ids.txt")
with open(ids_path, "w") as f:
    for j in jobs_exitosos:
        f.write(f"{j['par_id']}\t{j['job_id']}\t{j['pre_date']}→{j['post_date']}\n")
print(f"  ✓ Guardado: {ids_path}")

# ============================================================
# 6. RESUMEN
# ============================================================
print(f"\n{'=' * 70}")
print("RESUMEN P04")
print(f"{'=' * 70}")
print(f"  Jobs enviados: {len(jobs_exitosos)}/{len(pares)}")
if jobs_fallidos:
    print(f"  ⚠ Jobs fallidos: {len(jobs_fallidos)}")
    for j in jobs_fallidos:
        print(f"    - {j['par_id']}: {j.get('error', 'desconocido')}")
print(f"\n  Archivos en: {JOBS_DIR}")
print(f"\n  TIEMPO ESTIMADO DE PROCESAMIENTO:")
print(f"  ~30-60 minutos por par × {len(jobs_exitosos)} pares")
print(f"  Total estimado: {len(jobs_exitosos) * 0.5:.0f}–{len(jobs_exitosos) * 1:.0f} horas")
print(f"\n  Para verificar estado, ejecutar P05 (Monitoreo)")
print(f"{'=' * 70}")