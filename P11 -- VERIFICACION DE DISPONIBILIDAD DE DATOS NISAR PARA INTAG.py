# -*- coding: utf-8 -*-
# =============================================================================
# P11 -- VERIFICACION DE DISPONIBILIDAD DE DATOS NISAR PARA INTAG
# =============================================================================
# Proyecto : InSAR Deslizamientos Intag -- Fase 4 NISAR
# Autor    : Victor Pinto Paez
# Fecha    : 2026-03-08
# Version  : 1.0
#
# PROPOSITO:
#   Consultar el catalogo ASF para verificar si existen datos NISAR
#   L-band disponibles sobre la zona de Intag, Cotacachi, Ecuador.
#   Esto determina si la Fase 4 del Protocolo Maestro v2.0 puede
#   iniciarse inmediatamente.
#
# REFERENCIA:
#   - NISAR lanzado 30/07/2025, operacional desde 01/2026
#   - Primera liberacion de datos: 23/01/2026 (25 muestras)
#   - Segunda liberacion: 27/02/2026 (>100,000 archivos globales)
#   - Datos calibrados globales: previsto mayo/junio 2026
#   - Datos accesibles via ASF Vertex y Earthdata Search
# =============================================================================

import sys
from datetime import datetime

try:
    import asf_search as asf
    print(f"asf_search version: {asf.__version__}")
except ImportError:
    print("[ERROR] asf_search no instalado.")
    print("  Ejecutar: pip install asf_search --break-system-packages")
    print("  O en conda: pip install asf_search")
    sys.exit(1)

# ============================================================
# CONFIGURACION
# ============================================================
# AOI de Intag (EPSG:4326)
AOI_WEST  = -79.285032
AOI_EAST  = -78.328185
AOI_SOUTH =   0.193747
AOI_NORTH =   0.555718

# WKT del AOI
AOI_WKT = (f"POLYGON(({AOI_WEST} {AOI_SOUTH}, {AOI_EAST} {AOI_SOUTH}, "
           f"{AOI_EAST} {AOI_NORTH}, {AOI_WEST} {AOI_NORTH}, "
           f"{AOI_WEST} {AOI_SOUTH}))")

print("=" * 70)
print("P12 -- VERIFICACION DE DISPONIBILIDAD NISAR PARA INTAG")
print("=" * 70)
print(f"Fecha: {datetime.now()}")
print(f"AOI: W={AOI_WEST}, E={AOI_EAST}, S={AOI_SOUTH}, N={AOI_NORTH}")

# ============================================================
# BUSQUEDA 1: NISAR L-band (todos los productos)
# ============================================================
print(f"\n{'-' * 70}")
print("BUSQUEDA 1: Datos NISAR sobre Intag")
print(f"{'-' * 70}")

try:
    results = asf.search(
        platform=asf.PLATFORM.NISAR,
        intersectsWith=AOI_WKT,
    )
    n_nisar = len(results)
    print(f"  Resultados NISAR totales: {n_nisar}")
except Exception as e:
    print(f"  Error en busqueda NISAR: {e}")
    n_nisar = 0
    results = []

if n_nisar > 0:
    print(f"\n  DATOS NISAR ENCONTRADOS!")
    print(f"\n  Detalle de los primeros {min(n_nisar, 20)} resultados:")
    print(f"  {'Nombre':<60} {'Fecha':>12} {'Tipo':>15}")
    print(f"  {'-'*60} {'-'*12} {'-'*15}")
    for i, r in enumerate(results[:20]):
        props = r.properties
        nombre = props.get("sceneName", "N/A")[:58]
        fecha = props.get("startTime", "N/A")[:10]
        tipo = props.get("processingLevel", "N/A")
        print(f"  {nombre:<60} {fecha:>12} {tipo:>15}")

    # Analisis por tipo de producto
    print(f"\n  Resumen por tipo de producto:")
    tipos = {}
    for r in results:
        t = r.properties.get("processingLevel", "N/A")
        tipos[t] = tipos.get(t, 0) + 1
    for t, n in sorted(tipos.items()):
        print(f"    {t}: {n} archivos")

    # Analisis por geometria (asc/desc)
    print(f"\n  Resumen por direccion orbital:")
    dirs = {}
    for r in results:
        d = r.properties.get("flightDirection", "N/A")
        dirs[d] = dirs.get(d, 0) + 1
    for d, n in sorted(dirs.items()):
        print(f"    {d}: {n} archivos")

    # Rango de fechas
    fechas = []
    for r in results:
        f_str = r.properties.get("startTime", "")[:10]
        if f_str:
            fechas.append(f_str)
    if fechas:
        fechas.sort()
        print(f"\n  Rango temporal: {fechas[0]} a {fechas[-1]}")
        print(f"  Fechas unicas: {len(set(fechas))}")

    # Listar productos SLC (necesarios para interferometria)
    slc_results = [r for r in results
                   if "SLC" in r.properties.get("processingLevel", "").upper()
                   or "RSLC" in r.properties.get("processingLevel", "").upper()
                   or "L1" in r.properties.get("processingLevel", "").upper()]
    print(f"\n  Productos SLC/RSLC/L1 (para interferometria): {len(slc_results)}")

else:
    print(f"\n  NO se encontraron datos NISAR sobre Intag.")

# ============================================================
# BUSQUEDA 2: NISAR con AOI ampliado (Ecuador)
# ============================================================
print(f"\n{'-' * 70}")
print("BUSQUEDA 2: Datos NISAR sobre Ecuador (AOI ampliado)")
print(f"{'-' * 70}")

# AOI ampliado: Ecuador continental
ECU_WKT = "POLYGON((-81 -2, -75 -2, -75 2, -81 2, -81 -2))"

try:
    results_ecu = asf.search(
        platform=asf.PLATFORM.NISAR,
        intersectsWith=ECU_WKT,
        maxResults=100,
    )
    n_ecu = len(results_ecu)
    print(f"  Resultados NISAR sobre Ecuador: {n_ecu}")
except Exception as e:
    print(f"  Error en busqueda Ecuador: {e}")
    n_ecu = 0
    results_ecu = []

if n_ecu > 0:
    fechas_ecu = []
    for r in results_ecu:
        f_str = r.properties.get("startTime", "")[:10]
        if f_str:
            fechas_ecu.append(f_str)
    if fechas_ecu:
        fechas_ecu.sort()
        print(f"  Rango temporal: {fechas_ecu[0]} a {fechas_ecu[-1]}")
        print(f"  Fechas unicas: {len(set(fechas_ecu))}")

    print(f"\n  Primeros 10 resultados:")
    for i, r in enumerate(results_ecu[:10]):
        props = r.properties
        nombre = props.get("sceneName", "N/A")[:55]
        fecha = props.get("startTime", "N/A")[:10]
        tipo = props.get("processingLevel", "N/A")
        print(f"    {nombre:<55} {fecha:>12} {tipo:>15}")
else:
    print(f"  NO se encontraron datos NISAR sobre Ecuador.")

# ============================================================
# BUSQUEDA 3: ALOS-2 PALSAR-2 (alternativa L-band)
# ============================================================
print(f"\n{'-' * 70}")
print("BUSQUEDA 3: ALOS-2 PALSAR-2 (alternativa L-band)")
print(f"{'-' * 70}")

try:
    results_alos = asf.search(
        platform=asf.PLATFORM.ALOS,
        intersectsWith=AOI_WKT,
        maxResults=50,
    )
    n_alos = len(results_alos)
    print(f"  Resultados ALOS/ALOS-2 sobre Intag: {n_alos}")
except Exception as e:
    print(f"  Error: {e}")
    n_alos = 0
    results_alos = []

if n_alos > 0:
    fechas_alos = sorted(set(
        r.properties.get("startTime", "")[:10] for r in results_alos
    ))
    print(f"  Rango: {fechas_alos[0]} a {fechas_alos[-1]}")
    print(f"  Fechas unicas: {len(fechas_alos)}")

# ============================================================
# RESUMEN Y RECOMENDACION
# ============================================================
print(f"\n{'=' * 70}")
print("RESUMEN DE DISPONIBILIDAD")
print(f"{'=' * 70}")
print(f"  NISAR sobre Intag: {n_nisar} productos")
print(f"  NISAR sobre Ecuador: {n_ecu} productos")
print(f"  ALOS/ALOS-2 sobre Intag: {n_alos} productos")

print(f"\n{'=' * 70}")
print("RECOMENDACION")
print(f"{'=' * 70}")

if n_nisar > 0:
    print(f"  DATOS NISAR DISPONIBLES PARA INTAG.")
    print(f"  Se puede iniciar Fase 4 inmediatamente.")
    print(f"  NOTA: Datos pre-calibracion. Datos calibrados: mayo/junio 2026.")
    print(f"  Para analisis de coherencia C vs L, los datos pre-calibracion")
    print(f"  son suficientes (la coherencia no depende de la calibracion")
    print(f"  radiometrica absoluta).")
elif n_ecu > 0:
    print(f"  Hay datos NISAR sobre Ecuador pero NO sobre Intag.")
    print(f"  Esperar proxima liberacion de datos (mayo/junio 2026).")
    print(f"  Mientras tanto, avanzar con Fase 5 (Integracion geoespacial).")
else:
    print(f"  NO hay datos NISAR disponibles aun para Ecuador.")
    print(f"  Datos calibrados globales previstos: mayo/junio 2026.")
    print(f"  Avanzar con Fase 5 mientras se espera.")

print(f"{'=' * 70}")