# =============================================================================
# P05c — RENOMBRAR ARCHIVOS CLIPPED PARA COMPATIBILIDAD CON P06
# =============================================================================
# Proyecto: InSAR Deslizamientos Zona de Intag, Cotacachi, Ecuador
# Autor:    Víctor Pinto
# Fecha:    2026-02-27
#
# PROPOSITO:
#   Los archivos en hyp3_clipped/ tienen nombres largos de HyP3:
#     S1AA_20231124T233754_20231218T233753_VVP024_INT80_G_ueF_3ED6_corr.tif
#
#   P06 v2.0 espera nombres cortos por par_id:
#     P18_B1_corr.tif
#
#   Este script COPIA (no mueve) los archivos con nombres cortos,
#   preservando los originales para trazabilidad.
#
# ENTRADA:  hyp3_clipped/{par_id}/*.tif (nombres largos HyP3)
# SALIDA:   hyp3_clipped/{par_id}/{par_id}_corr.tif, etc. (nombres cortos)
#
# EJECUTAR EN: Spyder (Run file F5)
# PREREQUISITO: P05 v2.0 ejecutado (archivos recortados existentes)
# SIGUIENTE:    P06 v2.0 (extraccion de zonas de deformacion)
# =============================================================================

import os
import shutil
import glob
from datetime import datetime

print("=" * 70)
print("P05c — RENOMBRAR ARCHIVOS CLIPPED PARA P06")
print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# -- Configuracion ------------------------------------------------------------
SAVE_DIR    = r"D:\POSGRADOS\INTAG\data\sentinel1"
CLIPPED_DIR = os.path.join(SAVE_DIR, "hyp3_clipped")

# Pares individuales que P06 necesita procesar
# (los mosaicos M_ASC_CO_12d y M_ASC_CO_24d ya funcionan)
PARES = ["D1", "D2", "D3", "D4", "A1", "A2", "A3", "A4",
         "P18_A3", "P18_B1", "P18_B2", "P18_B3"]

# Sufijos que P06 necesita (los 3 criticos + DEM para verificacion)
SUFIJOS_CRITICOS = ["_corr.tif", "_los_disp.tif", "_vert_disp.tif", "_dem.tif"]

# Sufijos adicionales (se renombran tambien para completitud)
SUFIJOS_EXTRA = ["_unw_phase.tif", "_wrapped_phase.tif", "_inc_map.tif",
                 "_amp.tif", "_lv_theta.tif", "_lv_phi.tif"]

TODOS_SUFIJOS = SUFIJOS_CRITICOS + SUFIJOS_EXTRA

print(f"\n  Directorio: {CLIPPED_DIR}")
print(f"  Pares a procesar: {len(PARES)}")
print(f"  Sufijos criticos: {len(SUFIJOS_CRITICOS)}")

# -- Procesamiento ------------------------------------------------------------
print(f"\n{'─'*70}")
print(f"  {'PAR':<10} {'SUFIJO':<20} {'ACCION':<15} {'ARCHIVO ORIGINAL'}")
print(f"{'─'*70}")

total_copiados = 0
total_existentes = 0
total_faltantes = 0
pares_listos = []
pares_incompletos = []

for par_id in PARES:
    carpeta = os.path.join(CLIPPED_DIR, par_id)

    if not os.path.exists(carpeta):
        print(f"  {par_id:<10} {'---':<20} {'NO EXISTE':15} carpeta no encontrada")
        pares_incompletos.append(par_id)
        continue

    archivos = os.listdir(carpeta)
    criticos_ok = 0

    for sufijo in TODOS_SUFIJOS:
        nombre_corto = f"{par_id}{sufijo}"
        ruta_corta = os.path.join(carpeta, nombre_corto)

        # Si ya existe el nombre corto, no hacer nada
        if os.path.exists(ruta_corta):
            if sufijo in SUFIJOS_CRITICOS:
                criticos_ok += 1
            total_existentes += 1
            print(f"  {par_id:<10} {sufijo:<20} {'YA EXISTE':15}")
            continue

        # Buscar archivo con nombre largo que termine en este sufijo
        candidatos = [a for a in archivos if a.endswith(sufijo)]

        if len(candidatos) == 1:
            ruta_larga = os.path.join(carpeta, candidatos[0])
            shutil.copy2(ruta_larga, ruta_corta)
            total_copiados += 1
            if sufijo in SUFIJOS_CRITICOS:
                criticos_ok += 1
            print(f"  {par_id:<10} {sufijo:<20} {'COPIADO':15} {candidatos[0][:50]}")
        elif len(candidatos) == 0:
            total_faltantes += 1
            print(f"  {par_id:<10} {sufijo:<20} {'NO ENCONTRADO':15}")
        else:
            # Multiples candidatos: tomar el primero y advertir
            ruta_larga = os.path.join(carpeta, candidatos[0])
            shutil.copy2(ruta_larga, ruta_corta)
            total_copiados += 1
            if sufijo in SUFIJOS_CRITICOS:
                criticos_ok += 1
            print(f"  {par_id:<10} {sufijo:<20} {'COPIADO*':15} {candidatos[0][:50]}")
            print(f"  {'':10} {'':20} {'AVISO: ':15} {len(candidatos)} candidatos, usando primero")

    if criticos_ok >= 3:  # corr + los_disp + al menos 1 mas
        pares_listos.append(par_id)
    else:
        pares_incompletos.append(par_id)

# -- Resumen ------------------------------------------------------------------
print(f"\n{'='*70}")
print("RESUMEN")
print(f"{'='*70}")
print(f"  Archivos copiados (renombrados): {total_copiados}")
print(f"  Archivos ya existentes:          {total_existentes}")
print(f"  Archivos no encontrados:         {total_faltantes}")
print(f"\n  Pares listos para P06:     {len(pares_listos)}")
for p in pares_listos:
    print(f"    ok {p}")
print(f"\n  Pares incompletos:         {len(pares_incompletos)}")
for p in pares_incompletos:
    print(f"    !! {p}")

# -- Verificacion critica para la Prueba 1 ------------------------------------
print(f"\n{'─'*70}")
print("VERIFICACION CRITICA: Pares para Prueba Temporal")
print(f"{'─'*70}")

for par_id in ["P18_B1", "P18_B3"]:
    carpeta = os.path.join(CLIPPED_DIR, par_id)
    corr = os.path.join(carpeta, f"{par_id}_corr.tif")
    los  = os.path.join(carpeta, f"{par_id}_los_disp.tif")

    corr_ok = os.path.exists(corr)
    los_ok  = os.path.exists(los)

    status = "LISTO" if (corr_ok and los_ok) else "FALTA"
    print(f"  {par_id}: {status}")
    print(f"    _corr.tif:     {'ok' if corr_ok else 'FALTA'}")
    print(f"    _los_disp.tif: {'ok' if los_ok else 'FALTA'}")

print(f"\n{'='*70}")
if "P18_B1" in pares_listos and "P18_B3" in pares_listos:
    print("  RESULTADO: P18_B1 y P18_B3 listos.")
    print("  SIGUIENTE PASO: Ejecutar P06 v2.0 — procesara los 12 datasets")
else:
    print("  RESULTADO: Faltan archivos criticos.")
    print("  ACCION: Verificar que P05 v2.0 proceso estos pares correctamente")
print(f"{'='*70}")