# -*- coding: utf-8 -*-
# =============================================================================
# P09f — DESCARGA DE PRODUCTOS SBAS FALTANTES (CUENTA 2)
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag — Fase 3 SBAS
# Autor    : Víctor Pinto Páez
# Fecha    : 2026-03-02
# Versión  : 1.0
#
# ⚠ SEGURIDAD: Borrar las credenciales de este archivo después de ejecutar
# =============================================================================

import os
import sys
import time
import csv
from datetime import datetime

try:
    import hyp3_sdk
except ImportError:
    print("[ERROR] pip install hyp3_sdk")
    sys.exit(1)

# ============================================================
# CONFIGURACIÓN
# ============================================================
SBAS_DIR     = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"
PRODUCTS_DIR = os.path.join(SBAS_DIR, "hyp3_products")
os.makedirs(PRODUCTS_DIR, exist_ok=True)

PROJECT_C2 = "Intag_SBAS_Fase3"

# ⚠ CREDENCIALES CUENTA 2 — BORRAR DESPUÉS DE EJECUTAR ⚠
EARTHDATA_USER = ""   # no subir a GitHub
EARTHDATA_PASS = ""   # no subir a GitHub

OUT_CSV     = os.path.join(SBAS_DIR, "P09f_inventario_descargas.csv")
OUT_REPORTE = os.path.join(SBAS_DIR, "P09f_reporte_descarga.txt")

print("=" * 70)
print("P09f — DESCARGA DE PRODUCTOS SBAS FALTANTES (CUENTA 2)")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"Directorio: {PRODUCTS_DIR}")

# ============================================================
# PASO 1: Inventario existente
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 1: Inventario de productos ya descargados")
print(f"{'-' * 70}")

zips_existentes = set()
for f in os.listdir(PRODUCTS_DIR):
    if f.endswith(".zip"):
        zips_existentes.add(f)

size_existente = sum(
    os.path.getsize(os.path.join(PRODUCTS_DIR, f)) for f in zips_existentes
)
print(f"  Zips existentes: {len(zips_existentes)}")
print(f"  Tamaño: {size_existente / (1024**3):.1f} GB")

# ============================================================
# PASO 2: Conectar a HyP3
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 2: Conexión a HyP3 (cuenta 2)")
print(f"{'-' * 70}")

try:
    hyp3 = hyp3_sdk.HyP3(username=EARTHDATA_USER, password=EARTHDATA_PASS)
    print(f"  ✓ Conectado como: {EARTHDATA_USER}")
except Exception as e:
    print(f"  ERROR: {e}")
    sys.exit(1)

# ============================================================
# PASO 3: Obtener jobs SUCCEEDED
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 3: Jobs completados")
print(f"{'-' * 70}")

jobs_c2 = hyp3.find_jobs(name=PROJECT_C2, status_code="SUCCEEDED")
print(f"  SUCCEEDED en '{PROJECT_C2}': {len(jobs_c2)}")

if len(jobs_c2) == 0:
    print(f"  Buscando sin filtro de nombre...")
    jobs_c2 = hyp3.find_jobs(status_code="SUCCEEDED")
    print(f"  SUCCEEDED (todos): {len(jobs_c2)}")

# ============================================================
# PASO 4: Filtrar faltantes
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 4: Filtrado")
print(f"{'-' * 70}")

jobs_para_descargar = []
for job in jobs_c2:
    ya_existe = False
    try:
        files = job.files if hasattr(job, 'files') else []
        for fi in files:
            fname = fi.get('filename', '') if isinstance(fi, dict) else getattr(fi, 'filename', '')
            if fname and fname in zips_existentes:
                ya_existe = True
                break
    except:
        pass
    if not ya_existe:
        jobs_para_descargar.append(job)

print(f"  Ya descargados: {len(jobs_c2) - len(jobs_para_descargar)}")
print(f"  Por descargar: {len(jobs_para_descargar)}")

if len(jobs_para_descargar) == 0 and len(zips_existentes) < 1140:
    print(f"  Descargando todos — omitirá existentes")
    jobs_para_descargar = list(jobs_c2)

# ============================================================
# PASO 5: Descarga
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 5: Descarga")
print(f"{'-' * 70}")

descargados = 0
errores = 0
errores_detalle = []

if len(jobs_para_descargar) > 0:
    est_gb = len(jobs_para_descargar) * 0.33
    print(f"  Descargando {len(jobs_para_descargar)} productos (~{est_gb:.0f} GB)")
    print(f"  Destino: {PRODUCTS_DIR}\n")

    t_inicio = time.time()
    for i, job in enumerate(jobs_para_descargar, 1):
        try:
            job.download_files(PRODUCTS_DIR)
            descargados += 1
            if i % 10 == 0 or i == len(jobs_para_descargar):
                t_el = (time.time() - t_inicio) / 60
                t_rest = (t_el / i) * (len(jobs_para_descargar) - i)
                print(f"  [{i}/{len(jobs_para_descargar)}] "
                      f"{descargados} OK, {errores} err | "
                      f"{t_el:.1f} min | ~{t_rest:.0f} min rest.")
        except Exception as e:
            errores += 1
            errores_detalle.append({"job_id": getattr(job, 'job_id', 'N/A'), "error": str(e)})
            if errores <= 10:
                print(f"  ERROR [{i}]: {str(e)[:80]}")

    print(f"\n  Resultado: {descargados} OK, {errores} errores")
    print(f"  Tiempo: {(time.time() - t_inicio) / 60:.1f} min")
else:
    print(f"  ✓ Todo descargado.")

# ============================================================
# PASO 6: Verificación
# ============================================================
print(f"\n{'-' * 70}")
print("PASO 6: Verificación")
print(f"{'-' * 70}")

zips_final = [f for f in os.listdir(PRODUCTS_DIR) if f.endswith(".zip")]
size_final = sum(os.path.getsize(os.path.join(PRODUCTS_DIR, f)) for f in zips_final)

print(f"  Zips totales: {len(zips_final)}")
print(f"  Tamaño: {size_final / (1024**3):.1f} GB")
print(f"  Esperados: 1140")
print(f"  Completitud: {len(zips_final)/1140*100:.1f}%")

# ============================================================
# PASO 7: Inventario
# ============================================================
with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["archivo", "tamaño_MB"])
    for z in sorted(zips_final):
        size_mb = os.path.getsize(os.path.join(PRODUCTS_DIR, z)) / (1024**2)
        writer.writerow([z, f"{size_mb:.1f}"])
print(f"  ✓ {OUT_CSV}")

with open(OUT_REPORTE, "w", encoding="utf-8") as f:
    f.write(f"P09f — DESCARGA SBAS\nFecha: {datetime.now()}\n")
    f.write(f"C1: 768 | C2: {descargados} | Total: {len(zips_final)}\n")
    f.write(f"Tamaño: {size_final/(1024**3):.1f} GB\n")
print(f"  ✓ {OUT_REPORTE}")

print(f"\n{'=' * 70}")
print("P09f COMPLETADO")
print(f"{'=' * 70}")
print(f"  {len(zips_final)} zips ({size_final / (1024**3):.1f} GB)")
print(f"\n  ⚠ BORRAR CREDENCIALES DE ESTE SCRIPT")
print(f"  SIGUIENTE: P10 — Preparación para MintPy")
print(f"{'=' * 70}")