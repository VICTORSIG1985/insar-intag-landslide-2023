# =============================================================================
# P09e — DESCARGA DE PRODUCTOS SBAS COMPLETADOS (CUENTA 1)
# =============================================================================
# Proyecto: Detección de deslizamientos mediante InSAR — Zona de Intag
# Fase 3: Análisis SBAS de series temporales
#
# Propósito:
#   Descargar los 768 productos completados de cuenta 1.
#   (El diagnóstico de fallos y reenvío ya se ejecutó exitosamente)
#
# Nota: Los 9 pares fallidos ya fueron reenviados. Este script solo descarga.
# =============================================================================

import os
import sys
import time
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

OUTPUT_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"
DOWNLOAD_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3\hyp3_products"
HYP3_PROJECT_C1 = "INTAG_SBAS_Fase3_Path40"

# =============================================================================
# EJECUCIÓN
# =============================================================================

print("=" * 70)
print("P09e — DESCARGA DE PRODUCTOS SBAS (CUENTA 1)")
print("=" * 70)
print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"  Destino: {DOWNLOAD_DIR}")
print()

try:
    import hyp3_sdk
except ImportError:
    print("[ERROR] hyp3_sdk no instalado")
    sys.exit(1)

# --- Conectar ---
print("--- Conectando ---")
hyp3 = hyp3_sdk.HyP3()
print("  ✓ Conectado")

# --- Recuperar jobs completados ---
print("--- Recuperando jobs completados ---")
jobs = hyp3.find_jobs(name=HYP3_PROJECT_C1, status_code='SUCCEEDED')
job_list = list(jobs)
print(f"  Jobs SUCCEEDED: {len(job_list)}")

# --- Crear directorio ---
os.makedirs(DOWNLOAD_DIR, exist_ok=True)

# --- Verificar qué ya está descargado ---
ya_descargados = set()
for f in os.listdir(DOWNLOAD_DIR):
    if f.endswith('.zip'):
        ya_descargados.add(f)
print(f"  Ya descargados: {len(ya_descargados)} archivos .zip")

# --- Filtrar pendientes ---
# Cada job genera un .zip cuyo nombre contiene el job_id o los granules
# Usamos download_files() que automáticamente salta archivos existentes

por_descargar = []
for j in job_list:
    # Verificar si algún archivo de este job ya existe
    skip = False
    try:
        for fi in j.files:
            if fi['filename'] in ya_descargados:
                skip = True
                break
    except:
        pass
    if not skip:
        por_descargar.append(j)

print(f"  Pendientes de descarga: {len(por_descargar)}")

if len(por_descargar) == 0:
    print("  ✓ Todos los productos ya están descargados.")
    sys.exit(0)

# Estimar
est_gb = len(por_descargar) * 0.3
print(f"  Tamaño estimado: ~{est_gb:.0f} GB")
print()

# --- Descargar ---
print("-" * 70)
print(f"DESCARGANDO {len(por_descargar)} PRODUCTOS")
print("-" * 70)

t_start = time.time()
descargados = 0
errores = 0
errores_detalle = []

for i, job in enumerate(por_descargar, 1):
    try:
        job.download_files(DOWNLOAD_DIR)
        descargados += 1
    except Exception as e:
        errores += 1
        error_msg = f"{type(e).__name__}: {str(e)[:100]}"
        if errores <= 3:
            print(f"  ✗ Error job {job.job_id[:12]}...: {error_msg}")
            errores_detalle.append({'job_id': job.job_id, 'error': error_msg})

    # Progreso cada 25
    if i % 25 == 0 or i == len(por_descargar):
        elapsed = time.time() - t_start
        rate = descargados / elapsed * 3600 if elapsed > 0 and descargados > 0 else 0
        eta_h = (len(por_descargar) - i) / (i / elapsed * 3600) if elapsed > 0 else 0
        size_gb = sum(
            os.path.getsize(os.path.join(DOWNLOAD_DIR, f))
            for f in os.listdir(DOWNLOAD_DIR)
            if os.path.isfile(os.path.join(DOWNLOAD_DIR, f))
        ) / 1e9
        print(f"  [{i}/{len(por_descargar)}] {descargados} OK, {errores} err "
              f"| {elapsed/60:.1f} min | {size_gb:.1f} GB descargados")

# --- Resumen ---
elapsed_total = time.time() - t_start

# Espacio total
total_size = sum(
    os.path.getsize(os.path.join(DOWNLOAD_DIR, f))
    for f in os.listdir(DOWNLOAD_DIR)
    if os.path.isfile(os.path.join(DOWNLOAD_DIR, f))
) / 1e9
n_files = len([f for f in os.listdir(DOWNLOAD_DIR) if f.endswith('.zip')])

print()
print("=" * 70)
print("P09e COMPLETADO")
print("=" * 70)
print(f"  Descargados: {descargados}")
print(f"  Errores: {errores}")
print(f"  Duración: {elapsed_total/60:.1f} minutos")
print(f"  Espacio total: {total_size:.1f} GB ({n_files} archivos)")
print(f"  Directorio: {DOWNLOAD_DIR}")

if errores > 3:
    print(f"\n  ⚠ {errores} errores. Re-ejecutar el script descargará solo los faltantes.")

print()
print("  ESTADO:")
print(f"    Cuenta 1: {descargados + len(ya_descargados)} productos descargados")
print(f"    Cuenta 2: ejecutar P09d para monitorear, P09f para descargar")
print("=" * 70)