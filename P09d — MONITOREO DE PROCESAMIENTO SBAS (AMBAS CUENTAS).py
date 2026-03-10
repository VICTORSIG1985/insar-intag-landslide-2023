# =============================================================================
# P09d — MONITOREO DE PROCESAMIENTO SBAS (AMBAS CUENTAS)
# =============================================================================
# Proyecto: Detección de deslizamientos mediante InSAR — Zona de Intag
# Fase 3: Análisis SBAS de series temporales
#
# Propósito:
#   Consultar el estado de los 1,126 jobs SBAS distribuidos en dos cuentas
#   NASA Earthdata y mostrar progreso combinado.
#
# Nomenclatura:
#   P09  = Diseño y envío del stack SBAS (754 jobs, cuenta 1)
#   P09b = Diagnóstico de envío parcial
#   P09c = Envío de 372 pares pendientes (cuenta 2)
#   P09d = Monitoreo de procesamiento [ESTE SCRIPT]
#
# Uso:
#   Ejecutar periódicamente para verificar avance del procesamiento.
#   No requiere que la PC esté encendida entre ejecuciones.
# =============================================================================

import os
import sys
import json
import platform
import time
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

OUTPUT_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1\sbas_fase3"

# Cuenta 1 (principal) — usa _netrc existente
# Cuenta 2
CUENTA2_USER = Colocar
CUENTA2_PASS = Colocar

# Nombres de proyecto por cuenta (P09 usó nombre diferente a P09c)
HYP3_PROJECT_C1 = "INTAG_SBAS_Fase3_Path40"   # Cuenta 1 (P09 original)
HYP3_PROJECT_C2 = "Intag_SBAS_Fase3"           # Cuenta 2 (P09c)

# =============================================================================
# FUNCIONES
# =============================================================================

def consultar_cuenta(hyp3, nombre_cuenta, n_esperados, project_name):
    """Consulta estado de jobs en una cuenta HyP3."""
    try:
        jobs = hyp3.find_jobs(name=project_name)
        job_list = list(jobs)
    except Exception as e:
        print(f"  [ERROR] {nombre_cuenta}: {e}")
        return None

    estados = {}
    for j in job_list:
        st = j.status_code
        estados[st] = estados.get(st, 0) + 1

    # Créditos
    creditos = None
    try:
        info = hyp3.my_info()
        creditos = info.get('remaining_credits', None)
    except:
        pass

    return {
        'nombre': nombre_cuenta,
        'total': len(job_list),
        'esperados': n_esperados,
        'estados': estados,
        'creditos': creditos,
        'jobs': job_list,
    }


def mostrar_cuenta(info):
    """Muestra resumen de una cuenta."""
    if info is None:
        return

    print(f"\n  {info['nombre']} ({info['total']}/{info['esperados']} jobs)")
    print(f"  {'─' * 50}")

    # Estados con iconos
    iconos = {
        'SUCCEEDED': '✓', 'RUNNING': '⚙', 'PENDING': '◌',
        'FAILED': '✗', 'EXPIRED': '⏰',
    }

    for st in ['SUCCEEDED', 'RUNNING', 'PENDING', 'FAILED', 'EXPIRED']:
        n = info['estados'].get(st, 0)
        if n > 0:
            icono = iconos.get(st, '?')
            pct = n / info['total'] * 100 if info['total'] > 0 else 0
            barra = '█' * int(pct / 5) + '░' * (20 - int(pct / 5))
            print(f"    {icono} {st:<12} {n:>5}  {barra} {pct:5.1f}%")

    if info['creditos'] is not None:
        print(f"    Créditos restantes: {info['creditos']:,.0f}")


# =============================================================================
# EJECUCIÓN
# =============================================================================

print("=" * 70)
print("P09d — MONITOREO DE PROCESAMIENTO SBAS")
print("=" * 70)
print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"  Proyecto cuenta 1: {HYP3_PROJECT_C1}")
print(f"  Proyecto cuenta 2: {HYP3_PROJECT_C2}")
print()

try:
    import hyp3_sdk
except ImportError:
    print("[ERROR] hyp3_sdk no instalado")
    sys.exit(1)

# --- Ruta _netrc ---
if platform.system() == 'Windows':
    netrc_path = os.path.join(os.path.expanduser("~"), "_netrc")
else:
    netrc_path = os.path.join(os.path.expanduser("~"), ".netrc")

# --- Cuenta 1 (principal, _netrc actual) ---
print("─" * 70)
print("Consultando cuenta 1 (principal)...")
try:
    hyp3_c1 = hyp3_sdk.HyP3()
    info_c1 = consultar_cuenta(hyp3_c1, "CUENTA 1 (principal)", 754, HYP3_PROJECT_C1)
except Exception as e:
    print(f"  [ERROR] Conexión cuenta 1: {e}")
    info_c1 = None

# --- Cuenta 2 (temporal via _netrc swap) ---
print("Consultando cuenta 2 (colocar)...")

# Backup _netrc
netrc_backup = None
if os.path.exists(netrc_path):
    with open(netrc_path, 'r') as f:
        netrc_backup = f.read()

# Escribir credenciales cuenta 2
netrc_c2 = f"""machine urs.earthdata.nasa.gov
login {CUENTA2_USER}
password {CUENTA2_PASS}
"""
with open(netrc_path, 'w') as f:
    f.write(netrc_c2)

try:
    hyp3_c2 = hyp3_sdk.HyP3()
    info_c2 = consultar_cuenta(hyp3_c2, "CUENTA 2 (vicsig1985)", 372, HYP3_PROJECT_C2)
except Exception as e:
    print(f"  [ERROR] Conexión cuenta 2: {e}")
    info_c2 = None

# Restaurar _netrc
if netrc_backup:
    with open(netrc_path, 'w') as f:
        f.write(netrc_backup)
elif os.path.exists(netrc_path):
    os.remove(netrc_path)

# =============================================================================
# RESULTADOS
# =============================================================================

print()
print("=" * 70)
print("ESTADO DEL PROCESAMIENTO")
print("=" * 70)

mostrar_cuenta(info_c1)
mostrar_cuenta(info_c2)

# --- Resumen combinado ---
print()
print("─" * 70)
print("RESUMEN COMBINADO")
print("─" * 70)

estados_total = {}
total_jobs = 0

for info in [info_c1, info_c2]:
    if info:
        total_jobs += info['total']
        for st, n in info['estados'].items():
            estados_total[st] = estados_total.get(st, 0) + n

succeeded = estados_total.get('SUCCEEDED', 0)
running = estados_total.get('RUNNING', 0)
pending = estados_total.get('PENDING', 0)
failed = estados_total.get('FAILED', 0)

print(f"\n  Total jobs:     {total_jobs} / 1,126")
print(f"  ✓ Completados:  {succeeded}")
print(f"  ⚙ Procesando:   {running}")
print(f"  ◌ En cola:       {pending}")
if failed > 0:
    print(f"  ✗ Fallidos:      {failed}")

if total_jobs > 0:
    pct_total = succeeded / total_jobs * 100
    barra = '█' * int(pct_total / 2) + '░' * (50 - int(pct_total / 2))
    print(f"\n  Progreso: [{barra}] {pct_total:.1f}%")

    if succeeded == total_jobs:
        print("\n  🎉 PROCESAMIENTO COMPLETO — Listo para P10 (descarga)")
    elif running + pending > 0:
        # Estimar tiempo restante (HyP3 ~20-40 min por job)
        restantes = running + pending
        eta_min_h = restantes * 20 / 60
        eta_max_h = restantes * 40 / 60
        print(f"\n  Estimación restante: {eta_min_h:.0f}–{eta_max_h:.0f} horas")
        print(f"  (basado en {restantes} jobs pendientes, ~20-40 min/job)")

# --- Guardar snapshot ---
snapshot = {
    'fecha': datetime.now().isoformat(),
    'proyecto_c1': HYP3_PROJECT_C1,
    'proyecto_c2': HYP3_PROJECT_C2,
    'total_jobs': total_jobs,
    'estados': estados_total,
    'cuenta1': {
        'total': info_c1['total'] if info_c1 else 0,
        'estados': info_c1['estados'] if info_c1 else {},
        'creditos': info_c1['creditos'] if info_c1 else None,
    },
    'cuenta2': {
        'total': info_c2['total'] if info_c2 else 0,
        'estados': info_c2['estados'] if info_c2 else {},
        'creditos': info_c2['creditos'] if info_c2 else None,
    },
}

snapshot_path = os.path.join(OUTPUT_DIR, "P09d_estado_procesamiento.json")
with open(snapshot_path, 'w') as f:
    json.dump(snapshot, f, indent=2, ensure_ascii=False)

print(f"\n  Snapshot guardado: {os.path.basename(snapshot_path)}")
print()
print("=" * 70)
print("  Ejecutar este script periódicamente para verificar avance.")
print("  Cuando 100% completado → ejecutar P10 (descarga).")
print("=" * 70)