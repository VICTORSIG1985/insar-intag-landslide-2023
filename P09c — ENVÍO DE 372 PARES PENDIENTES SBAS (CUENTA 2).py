# =============================================================================
# P09c — ENVÍO DE 372 PARES PENDIENTES SBAS (CUENTA 2)
# =============================================================================
# Proyecto: Detección de deslizamientos mediante InSAR — Zona de Intag
# Fase 3: Análisis SBAS de series temporales
#
# Propósito:
#   Enviar los 372 pares interferométricos que no pudieron enviarse con la
#   cuenta principal (créditos agotados tras 754 jobs). Usa una segunda
#   cuenta NASA Earthdata para procesamiento en paralelo.
#
# Entrada:
#   P09b_pares_pendientes.csv (generado por P09b_diagnostico_envio_SBAS.py)
#
# Parámetros HyP3:
#   Idénticos a Fase 2 (P04/P04b) y a los 754 jobs ya enviados (cuenta 1)
#
# Salida:
#   P09c_jobs_cuenta2.json — registro completo de jobs enviados
#   P09c_reporte_envio.txt — resumen del envío
#
# Nomenclatura:
#   P09  = Diseño y envío del stack SBAS a HyP3 (754 jobs, cuenta 1)
#   P09b = Diagnóstico de envío parcial (754 enviados, 372 pendientes)
#   P09c = Envío de 372 pares pendientes (cuenta 2) [ESTE SCRIPT]
# =============================================================================

import os
import sys
import json
import time
import platform
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

OUTPUT_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"

# Archivo de pares pendientes (generado por P09b)
# Busca ambos nombres posibles por si el usuario usó la nomenclatura anterior
CSV_PENDIENTES_OPCIONES = [
    os.path.join(OUTPUT_DIR, "P09b_pares_pendientes.csv"),
    os.path.join(OUTPUT_DIR, "P09c_pares_pendientes.csv"),  # nombre anterior
]

# Credenciales cuenta 2 NASA Earthdata
EARTHDATA_USER = Colocar
EARTHDATA_PASS = Colocar

# Nombre del proyecto en HyP3 (mismo que cuenta 1 para trazabilidad)
HYP3_PROJECT_NAME = "Intag_SBAS_Fase3"

# Parámetros HyP3 — IDÉNTICOS a Fase 2 (P04/P04b) y cuenta 1 (P09)
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

# =============================================================================
# PASO 0: Configurar _netrc para cuenta 2
# =============================================================================
# hyp3_sdk usa la librería 'requests' que lee credenciales de _netrc.
# Configuramos _netrc antes de conectar, igual que funcionó con cuenta 1.
# =============================================================================

print("=" * 70)
print("P09c — ENVÍO DE PARES PENDIENTES SBAS (CUENTA 2)")
print("=" * 70)
print(f"  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"  Proyecto HyP3: {HYP3_PROJECT_NAME}")
print(f"  Cuenta: {EARTHDATA_USER}")
print(f"  Directorio: {OUTPUT_DIR}")
print()

# Configurar _netrc
print("--- Configurando credenciales ---")
if platform.system() == 'Windows':
    netrc_path = os.path.join(os.path.expanduser("~"), "_netrc")
else:
    netrc_path = os.path.join(os.path.expanduser("~"), ".netrc")

# Guardar backup de _netrc existente (cuenta 1)
netrc_backup = None
if os.path.exists(netrc_path):
    with open(netrc_path, 'r') as f:
        netrc_backup = f.read()
    print(f"  Backup de _netrc existente guardado")

# Escribir credenciales cuenta 2
netrc_content = f"""machine urs.earthdata.nasa.gov
login {EARTHDATA_USER}
password {EARTHDATA_PASS}
"""
with open(netrc_path, 'w') as f:
    f.write(netrc_content)
print(f"  ✓ _netrc configurado para: {EARTHDATA_USER}")
print(f"    Ruta: {netrc_path}")

# =============================================================================
# PASO 1: Cargar pares pendientes
# =============================================================================

print()
print("--- Cargando pares pendientes ---")

CSV_PENDIENTES = None
for path in CSV_PENDIENTES_OPCIONES:
    if os.path.exists(path):
        CSV_PENDIENTES = path
        break

if CSV_PENDIENTES is None:
    print(f"[ERROR] No se encontró archivo de pares pendientes:")
    for p in CSV_PENDIENTES_OPCIONES:
        print(f"  - {p}")
    print("  Ejecutar primero P09b (diagnóstico)")
    # Restaurar _netrc original
    if netrc_backup:
        with open(netrc_path, 'w') as f:
            f.write(netrc_backup)
    sys.exit(1)

df = pd.read_csv(CSV_PENDIENTES)
print(f"  Archivo: {os.path.basename(CSV_PENDIENTES)}")
print(f"  Pares a enviar: {len(df)}")

# Verificar columnas
for col in ['ref_scene', 'sec_scene']:
    if col not in df.columns:
        print(f"[ERROR] Columna '{col}' no encontrada. Columnas: {list(df.columns)}")
        if netrc_backup:
            with open(netrc_path, 'w') as f:
                f.write(netrc_backup)
        sys.exit(1)

# Detectar columna de baseline temporal
bt_col = None
for candidate in ['temporal_baseline', 'btemp_days', 'btemp']:
    if candidate in df.columns:
        bt_col = candidate
        break

# Resumen
if 'categoria' in df.columns:
    print("  Desglose:")
    for cat in sorted(df['categoria'].unique()):
        print(f"    {cat}: {(df['categoria'] == cat).sum()} pares")

if 'es_puente' in df.columns:
    n_puente = df['es_puente'].sum()
    if n_puente > 0:
        print(f"    PUENTE: {n_puente} pares")

print()
print("  Parámetros HyP3 (idénticos a Fase 2 y cuenta 1):")
for k, v in HYP3_PARAMS.items():
    print(f"    {k}: {v}")

# =============================================================================
# PASO 2: Conectar a HyP3 y verificar créditos
# =============================================================================

print()
print("-" * 70)
print("PASO 2: Conectar a HyP3 y verificar créditos")
print("-" * 70)

try:
    import hyp3_sdk
    print(f"  hyp3_sdk versión: {hyp3_sdk.__version__}")
except ImportError:
    print("[ERROR] hyp3_sdk no instalado. pip install hyp3_sdk")
    if netrc_backup:
        with open(netrc_path, 'w') as f:
            f.write(netrc_backup)
    sys.exit(1)

print("  Conectando...")
try:
    hyp3 = hyp3_sdk.HyP3()
    print(f"  ✓ Conectado exitosamente como: {EARTHDATA_USER}")
except Exception as e:
    print(f"  [ERROR] Conexión fallida: {type(e).__name__}: {e}")
    print()
    print("  POSIBLES CAUSAS:")
    print("  1. La cuenta necesita autorizar la aplicación ASF:")
    print("     → Ir a https://urs.earthdata.nasa.gov/profile")
    print("     → Sección 'Applications > Authorized Apps'")
    print("     → Buscar y autorizar 'Alaska Satellite Facility'")
    print("  2. Credenciales incorrectas")
    print("  3. Sin conexión a internet")
    # Restaurar _netrc
    if netrc_backup:
        with open(netrc_path, 'w') as f:
            f.write(netrc_backup)
        print(f"\n  _netrc restaurado a cuenta original")
    sys.exit(1)

# --- Verificar créditos disponibles ---
print()
print("  Verificando créditos...")
try:
    user_info = hyp3.my_info()
    creditos_restantes = user_info.get('remaining_credits', None)
    cuota_mensual = user_info.get('monthly_credits', None)

    if creditos_restantes is not None:
        print(f"  ✓ Créditos disponibles: {creditos_restantes:,.0f}")
        if cuota_mensual:
            print(f"    Cuota mensual: {cuota_mensual:,.0f}")
            print(f"    Usados: {cuota_mensual - creditos_restantes:,.0f}")

        # Cada job InSAR cuesta 10 créditos
        COSTO_POR_JOB = 10
        creditos_necesarios = len(df) * COSTO_POR_JOB
        print(f"\n  Créditos necesarios: {creditos_necesarios:,} ({len(df)} pares × {COSTO_POR_JOB} créditos)")

        if creditos_restantes < creditos_necesarios:
            pares_posibles = int(creditos_restantes // COSTO_POR_JOB)
            print(f"  ⚠ INSUFICIENTES: solo alcanza para {pares_posibles} de {len(df)} pares")
            print(f"    Se enviarán los {pares_posibles} primeros pares")
            # No salir, enviar lo que se pueda
        else:
            print(f"  ✓ SUFICIENTES: {creditos_restantes:,.0f} ≥ {creditos_necesarios:,}")
    else:
        print("  ⚠ No se pudo obtener info de créditos (continuando de todas formas)")
        # Mostrar toda la info disponible
        for k, v in user_info.items():
            print(f"    {k}: {v}")

except Exception as e:
    print(f"  ⚠ No se pudo verificar créditos: {type(e).__name__}: {e}")
    print("    Continuando con el envío...")

# =============================================================================
# PASO 3: Enviar jobs
# =============================================================================

print()
print("-" * 70)
print(f"PASO 3: ENVIANDO {len(df)} PARES A HyP3")
print("-" * 70)
print(f"  Estimado: ~{len(df) * 0.8 / 60:.0f}–{len(df) * 1.5 / 60:.0f} minutos")
print()

all_submitted = []
t_start = time.time()
total_pares = len(df)

for i, (_, par) in enumerate(df.iterrows(), 1):
    try:
        batch = hyp3.submit_insar_job(
            granule1=par['ref_scene'],
            granule2=par['sec_scene'],
            name=HYP3_PROJECT_NAME,
            **HYP3_PARAMS,
        )
        job = batch.jobs[0]

        registro = {
            'job_id': job.job_id,
            'ref_scene': par['ref_scene'],
            'sec_scene': par['sec_scene'],
            'status': job.status_code,
        }

        # Campos opcionales
        if 'ref_date' in par.index:
            registro['ref_date'] = str(par['ref_date'])
        if 'sec_date' in par.index:
            registro['sec_date'] = str(par['sec_date'])
        if bt_col and bt_col in par.index:
            registro['temporal_baseline'] = int(par[bt_col])
        if 'categoria' in par.index:
            registro['categoria'] = par['categoria']
        if 'es_puente' in par.index:
            registro['es_puente'] = bool(par['es_puente'])

        all_submitted.append(registro)

    except Exception as e:
        error_msg = str(e)
        print(f"  ⚠ ERROR par {i}: {error_msg[:120]}")

        registro_error = {
            'ref_scene': par['ref_scene'],
            'sec_scene': par['sec_scene'],
            'error': error_msg,
        }
        if 'ref_date' in par.index:
            registro_error['ref_date'] = str(par['ref_date'])
        if 'sec_date' in par.index:
            registro_error['sec_date'] = str(par['sec_date'])

        all_submitted.append(registro_error)

        # Si es error de créditos, parar inmediatamente
        if "credit" in error_msg.lower() or "remaining" in error_msg.lower():
            n_ok_parcial = len([j for j in all_submitted if 'job_id' in j])
            print(f"\n  ⛔ CRÉDITOS AGOTADOS en par {i}/{total_pares}")
            print(f"  Jobs enviados exitosamente: {n_ok_parcial}")
            break

    # Progreso cada 25 jobs
    if i % 25 == 0 or i == total_pares:
        elapsed = time.time() - t_start
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total_pares - i) / rate / 60 if rate > 0 else 0
        n_ok = len([j for j in all_submitted if 'job_id' in j])
        n_err = len([j for j in all_submitted if 'error' in j])
        print(f"  [{i}/{total_pares}] {n_ok} OK, {n_err} err "
              f"| {elapsed/60:.1f} min | ~{eta:.0f} min restantes")

    # Backup cada 100 jobs
    if i % 100 == 0:
        backup_path = os.path.join(OUTPUT_DIR, "P09c_backup_parcial.json")
        with open(backup_path, 'w') as f:
            json.dump(all_submitted, f, indent=2)

    # Pausa anti rate-limiting
    time.sleep(0.5)

# =============================================================================
# PASO 4: Guardar registro y restaurar _netrc
# =============================================================================

elapsed_total = time.time() - t_start
n_ok = len([j for j in all_submitted if 'job_id' in j])
n_err = len([j for j in all_submitted if 'error' in j])

# Guardar JSON de jobs
jobs_record = {
    'script': 'P09c_envio_pendientes_cuenta2.py',
    'proyecto': HYP3_PROJECT_NAME,
    'cuenta': EARTHDATA_USER,
    'fecha_envio': datetime.now().isoformat(),
    'duracion_minutos': round(elapsed_total / 60, 1),
    'total_enviados': n_ok,
    'total_errores': n_err,
    'parametros': HYP3_PARAMS,
    'nota': 'Complemento de 754 jobs enviados con cuenta principal (P09)',
    'jobs': all_submitted
}

jobs_path = os.path.join(OUTPUT_DIR, "P09c_jobs_cuenta2.json")
with open(jobs_path, 'w') as f:
    json.dump(jobs_record, f, indent=2, ensure_ascii=False)

# Generar reporte
reporte_lines = [
    "=" * 70,
    "P09c — REPORTE DE ENVÍO (CUENTA 2)",
    "=" * 70,
    f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    f"Cuenta: {EARTHDATA_USER}",
    f"Proyecto HyP3: {HYP3_PROJECT_NAME}",
    "",
    f"Pares procesados: {len(all_submitted)}",
    f"  ✓ Enviados: {n_ok}",
    f"  ✗ Errores: {n_err}",
    f"Duración: {elapsed_total/60:.1f} minutos",
    "",
    "ESTADO COMBINADO DE LA RED SBAS:",
    f"  Cuenta 1 (principal): 754 pares",
    f"  Cuenta 2 (esta):      {n_ok} pares",
    f"  Total en HyP3:        {754 + n_ok} pares",
    f"  Red completa (1,126): {'SÍ' if (754 + n_ok) >= 1126 else 'NO (' + str(1126 - 754 - n_ok) + ' faltantes)'}",
    "",
    "Parámetros (idénticos a Fase 2 y cuenta 1):",
]
for k, v in HYP3_PARAMS.items():
    reporte_lines.append(f"  {k}: {v}")
reporte_lines.extend(["", f"Registro: {jobs_path}", "=" * 70])

reporte_path = os.path.join(OUTPUT_DIR, "P09c_reporte_envio.txt")
with open(reporte_path, 'w', encoding='utf-8') as f:
    f.write("\n".join(reporte_lines))

# Limpiar backup parcial
backup_path = os.path.join(OUTPUT_DIR, "P09c_backup_parcial.json")
if os.path.exists(backup_path):
    os.remove(backup_path)

# RESTAURAR _netrc a cuenta original
print()
print("--- Restaurando _netrc ---")
if netrc_backup:
    with open(netrc_path, 'w') as f:
        f.write(netrc_backup)
    print(f"  ✓ _netrc restaurado a cuenta original")
else:
    # Si no había _netrc antes, eliminar el que creamos
    os.remove(netrc_path)
    print(f"  ✓ _netrc removido (no existía antes)")

# =============================================================================
# RESUMEN FINAL
# =============================================================================

print()
print("=" * 70)
print("P09c COMPLETADO")
print("=" * 70)
print()
print(f"  ✓ Jobs enviados: {n_ok}")
print(f"  ✗ Errores: {n_err}")
print(f"  Duración: {elapsed_total/60:.1f} minutos")
print()
print("  ESTADO COMBINADO:")
print(f"    Cuenta 1 (principal): 754 pares (2017-01 a 2021-10)")
print(f"    Cuenta 2 (esta):      {n_ok} pares (2018-05 a 2024-06)")
print(f"    Total en HyP3:        {754 + n_ok} de 1,126")
red_completa = (754 + n_ok) >= 1126
print(f"    Red completa:          {'SÍ ✓' if red_completa else 'NO'}")
print()
print(f"  Archivos generados:")
print(f"    P09c_jobs_cuenta2.json ({os.path.getsize(jobs_path)/1024:.1f} KB)")
print(f"    P09c_reporte_envio.txt")
print()
print("  Monitorear procesamiento:")
print("    https://hyp3-api.asf.alaska.edu/")
print("=" * 70)