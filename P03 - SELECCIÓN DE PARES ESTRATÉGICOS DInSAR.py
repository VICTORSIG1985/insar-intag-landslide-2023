# ============================================================
# P03 - SELECCIÓN DE PARES ESTRATÉGICOS DInSAR
# ============================================================
# Selecciona 8 pares interferométricos de 2 tracks verificados
# con cobertura real sobre las 6 parroquias.
# Criterios: cobertura, baseline temporal, propósito científico.
# Guardado: D:\POSGRADOS\INTAG\data\sentinel1\
# ============================================================
import json, os, csv
from datetime import datetime

SAVE_DIR = r"D:\POSGRADOS\INTAG\data\sentinel1"
EVENT = "2023-12-19"

print("=" * 70)
print("P03 - SELECCIÓN DE PARES ESTRATÉGICOS DInSAR")
print("=" * 70)

# ============================================================
# 1. DEFINICIÓN DE PARES SELECCIONADOS
# ============================================================
# Basado en:
#   - Verificación de cobertura espacial sobre 6 parroquias (GPKG)
#   - DESC Path:40 = 6/6 parroquias (100%)
#   - ASC Path:120 = 5/6 parroquias (García Moreno parcial 46.8%)
#   - ASC Path:18 = descartado (1/6 parroquias)
#   - DESC Path:142 = descartado (no cubre García Moreno)

pares = [
    # ── DESCENDENTE Path:40 Frame:590 (6/6 parroquias) ──
    {
        "id": "D1",
        "flightDirection": "DESCENDING",
        "pathNumber": 40,
        "frameNumber": 590,
        "cobertura": "6/6 parroquias (100%)",
        "pre_date": "2023-12-08",
        "pre_scene": "S1A_IW_SLC__1SDV_20231208T110106_20231208T110136_051562_063984_C1B0",
        "post_date": "2023-12-20",
        "post_scene": "S1A_IW_SLC__1SDV_20231220T110105_20231220T110134_051737_063FA1_0CB8",
        "baseline_dias": 12,
        "calidad": "★★★",
        "proposito": "Co-evento principal. Imagen post a 1 día del deslizamiento. Deformación inmediata.",
    },
    {
        "id": "D2",
        "flightDirection": "DESCENDING",
        "pathNumber": 40,
        "frameNumber": 590,
        "cobertura": "6/6 parroquias (100%)",
        "pre_date": "2023-11-26",
        "pre_scene": "S1A_IW_SLC__1SDV_20231126T110106_20231126T110136_051387_06337A_9395",
        "post_date": "2023-12-20",
        "post_scene": "S1A_IW_SLC__1SDV_20231220T110105_20231220T110134_051737_063FA1_0CB8",
        "baseline_dias": 24,
        "calidad": "★★★",
        "proposito": "Validación cruzada del co-evento con baseline mayor.",
    },
    {
        "id": "D3",
        "flightDirection": "DESCENDING",
        "pathNumber": 40,
        "frameNumber": 590,
        "cobertura": "6/6 parroquias (100%)",
        "pre_date": "2023-12-08",
        "pre_scene": "S1A_IW_SLC__1SDV_20231208T110106_20231208T110136_051562_063984_C1B0",
        "post_date": "2024-01-13",
        "post_scene": "S1A_IW_SLC__1SDV_20240113T110104_20240113T110133_052087_064B9D_5E75",
        "baseline_dias": 36,
        "calidad": "★★",
        "proposito": "Post-evento temprano. Detectar movimiento continuado tras el deslizamiento.",
    },
    {
        "id": "D4",
        "flightDirection": "DESCENDING",
        "pathNumber": 40,
        "frameNumber": 590,
        "cobertura": "6/6 parroquias (100%)",
        "pre_date": "2023-12-08",
        "pre_scene": "S1A_IW_SLC__1SDV_20231208T110106_20231208T110136_051562_063984_C1B0",
        "post_date": "2024-01-25",
        "post_scene": "S1A_IW_SLC__1SDV_20240125T110103_20240125T110133_052262_06517A_A2E4",
        "baseline_dias": 48,
        "calidad": "★★",
        "proposito": "Post-evento tardío. Evaluar estabilización o reactivación.",
    },
    # ── ASCENDENTE Path:120 Frame:1182 (5/6 parroquias) ──
    {
        "id": "A1",
        "flightDirection": "ASCENDING",
        "pathNumber": 120,
        "frameNumber": 1182,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_date": "2023-12-13",
        "pre_scene": "S1A_IW_SLC__1SDV_20231213T232947_20231213T233023_051642_063C64_14B1",
        "post_date": "2023-12-25",
        "post_scene": "S1A_IW_SLC__1SDV_20231225T232946_20231225T233022_051817_06426C_168C",
        "baseline_dias": 12,
        "calidad": "★★★",
        "proposito": "Co-evento desde geometría ascendente. Descomposición de deformación ASC+DESC.",
    },
    {
        "id": "A2",
        "flightDirection": "ASCENDING",
        "pathNumber": 120,
        "frameNumber": 1182,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_date": "2023-12-01",
        "pre_scene": "S1A_IW_SLC__1SDV_20231201T232947_20231201T233023_051467_06363D_8F6A",
        "post_date": "2023-12-25",
        "post_scene": "S1A_IW_SLC__1SDV_20231225T232946_20231225T233022_051817_06426C_168C",
        "baseline_dias": 24,
        "calidad": "★★★",
        "proposito": "Validación cruzada ASC con baseline mayor.",
    },
    {
        "id": "A3",
        "flightDirection": "ASCENDING",
        "pathNumber": 120,
        "frameNumber": 1182,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_date": "2023-12-13",
        "pre_scene": "S1A_IW_SLC__1SDV_20231213T232947_20231213T233023_051642_063C64_14B1",
        "post_date": "2024-01-06",
        "post_scene": "S1A_IW_SLC__1SDV_20240106T232945_20240106T233021_051992_06486B_CF9A",
        "baseline_dias": 24,
        "calidad": "★★★",
        "proposito": "Post-evento temprano desde geometría ascendente.",
    },
    {
        "id": "A4",
        "flightDirection": "ASCENDING",
        "pathNumber": 120,
        "frameNumber": 1182,
        "cobertura": "5/6 parroquias (García Moreno parcial 46.8%)",
        "pre_date": "2023-12-13",
        "pre_scene": "S1A_IW_SLC__1SDV_20231213T232947_20231213T233023_051642_063C64_14B1",
        "post_date": "2024-01-30",
        "post_scene": "S1A_IW_SLC__1SDV_20240130T232944_20240130T233020_052342_065447_619B",
        "baseline_dias": 48,
        "calidad": "★★",
        "proposito": "Post-evento tardío desde geometría ascendente.",
    },
]

# ============================================================
# 2. MOSTRAR RESUMEN
# ============================================================
print(f"\n--- 8 PARES SELECCIONADOS ---\n")

for p in pares:
    print(f"  {p['id']} {p['calidad']} {p['flightDirection']:10s} Path:{p['pathNumber']} "
          f"| {p['pre_date']} → {p['post_date']} ({p['baseline_dias']}d)")
    print(f"       Cobertura: {p['cobertura']}")
    print(f"       Propósito: {p['proposito']}\n")

# Estadísticas
desc_pares = [p for p in pares if p["flightDirection"] == "DESCENDING"]
asc_pares = [p for p in pares if p["flightDirection"] == "ASCENDING"]
tres = [p for p in pares if p["calidad"] == "★★★"]
dos = [p for p in pares if p["calidad"] == "★★"]

print(f"--- ESTADÍSTICAS ---")
print(f"  Total pares: {len(pares)}")
print(f"  Descendente (6/6 parroquias): {len(desc_pares)} pares")
print(f"  Ascendente (5/6 parroquias):  {len(asc_pares)} pares")
print(f"  Calidad ★★★ (≤24d): {len(tres)}")
print(f"  Calidad ★★  (25-48d): {len(dos)}")

# Escenas únicas a descargar
escenas_unicas = set()
for p in pares:
    escenas_unicas.add(p["pre_scene"])
    escenas_unicas.add(p["post_scene"])
print(f"  Escenas SLC únicas necesarias: {len(escenas_unicas)}")

# ============================================================
# 3. ANÁLISIS CIENTÍFICO
# ============================================================
print(f"\n--- CAPACIDAD ANALÍTICA ---")
print(f"  1. Deformación co-evento: D1 + A1 (12d, ambas geometrías)")
print(f"  2. Validación cruzada:    D1 vs D2, A1 vs A2")
print(f"  3. Serie post-evento:     D1→D3→D4 y A1→A3→A4")
print(f"  4. Descomposición 2D:     D1+A1 → vertical + horizontal")
print(f"  5. García Moreno:         Solo pares D1-D4 (cobertura 100%)")

# ============================================================
# 4. GUARDAR
# ============================================================
# JSON completo
output = {
    "generado": datetime.now().isoformat(),
    "evento": EVENT,
    "criterios_seleccion": {
        "cobertura_minima": "5/6 parroquias verificadas contra GPKG",
        "baseline_maxima_dias": 48,
        "tracks_descartados": [
            "ASC Path:18 Frame:1180 (1/6 parroquias)",
            "DESC Path:142 Frame:591 (no cubre García Moreno)",
        ],
        "tracks_seleccionados": [
            "DESC Path:40 Frame:590 (6/6 parroquias, 100%)",
            "ASC Path:120 Frame:1182 (5/6 parroquias, García Moreno parcial)",
        ],
    },
    "total_pares": len(pares),
    "escenas_unicas": len(escenas_unicas),
    "pares": pares,
}

json_path = os.path.join(SAVE_DIR, "pares_dinsar_optimos.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f"\n  ✓ Guardado: {json_path}")

# CSV
csv_path = os.path.join(SAVE_DIR, "pares_dinsar_resumen.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    campos = ["id", "flightDirection", "pathNumber", "frameNumber",
              "cobertura", "pre_date", "pre_scene", "post_date",
              "post_scene", "baseline_dias", "calidad", "proposito"]
    writer = csv.DictWriter(f, fieldnames=campos)
    writer.writeheader()
    writer.writerows(pares)
print(f"  ✓ Guardado: {csv_path}")

# Lista de escenas para descarga
escenas_path = os.path.join(SAVE_DIR, "escenas_para_descarga.txt")
with open(escenas_path, "w") as f:
    for s in sorted(escenas_unicas):
        f.write(s + "\n")
print(f"  ✓ Guardado: {escenas_path}")

print(f"\n" + "=" * 70)
print("Siguiente paso: P04 - Enviar 8 pares a HyP3 para procesamiento")
print("=" * 70)