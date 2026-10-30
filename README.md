# Sgt-GENETICO — Dirección A (equipo T16, Investigathón 2026)

Búsqueda numérica sobre la familia $S_7$ de O'Keeffe en el modelo de Kuramoto homogéneo: evaluador SDP con certificación espectral en 50 dígitos, algoritmo genético y ajuste local de masas. **No obtiene una cota nueva**; documenta un resultado negativo. Detalles en el informe técnico de la entrega.

> Todo el código de este repositorio fue escrito por un agente de IA (Google Antigravity). Las cifras del informe se reprodujeron ejecutando los scripts de abajo.

## Reproducir rápido

Requiere Python 3.12 y [`uv`](https://docs.astral.sh/uv/). Las versiones están fijadas en `uv.lock` (con otras versiones de NumPy/SciPy algunos puntos del solver fallan).

```bash
make env                              # crea el entorno e instala dependencias
uv run python scripts/final_thesis_front.py   # ~1 min · Tabla principal del informe
```

`final_thesis_front.py` imprime μ\* certificado de la plantilla $S_7$, de la semilla R3 y del campeón del GA, a tres márgenes λ_M. Reproduce la tabla central: a λ_M = 10⁻¹⁰, $S_7$ = 0,692009, R3 = 0,692557, campeón GA = 0,692453. Ninguno alcanza μ₀ = 0,6929077.

## Cada cifra del informe → su script

Todos se corren con `uv run python scripts/<archivo>`.

| Resultado | Script |
| --- | --- |
| Tabla principal (μ\* de S7, R3, GA a 3 márgenes) | `final_thesis_front.py` |
| Plantilla base vs O'Keeffe (barrido en h) | `run_s7_h_sweep.py` |
| Ajuste local de masas desde R3 → 0,69271 | `optimize_true_symmetry.py` |
| Óptimo con 6 grados de libertad (10 arranques) | `final_1_true_6d.py` |
| Ganancia vs h | `final_2_scaling_law.py` |
| Ganancia vs λ_M | `final_3_pareto_reopt.py` |
| Calibración del margen realizado | `audit_reopt_106.py` |
| Clarabel vs SCS (sensibilidad numérica) | `def_4_pareto_cross_solver.py` |
| Candado de escala (35,9 % NOT_SCALABLE) | `final_4_candado_stress.py` |

Corrida del GA de 12 h: registro en `runs/night12h_seeded_mu06923_s42/`. El ajuste local tiene dispersión de hasta 3·10⁻⁵ entre ejecuciones; los puntos fijos reproducen a 10⁻⁶.

## Verificación y tests

```bash
make reproduce-eval   # reporte JSON auditable del evaluador (results/reproduce_eval_report.json)
make test             # tests unitarios, de integración y de arquitectura
make lint             # ruff + mypy --strict
```

## Estructura

- `src/kuramoto_templates/` — paquete: `domain/` (matemática pura), `evaluators/` (SDP), `cases/` (S7, C4), `search/` (GA).
- `scripts/` — scripts que reproducen las cifras del informe (tabla de arriba).
- `exploraciones/` — Direcciones B y C, independientes del paquete. Abrir con `uv run --with jupyter --with matplotlib jupyter notebook exploraciones/`.
  - `codigo_direccion_B.ipynb` — cotas cuantitativas cerca de μ = 3/4.
  - `simulacion_direccion_B.html` — página interactiva de los cuatro grupos (abrir en el navegador).
  - `codigo_direccion_C.ipynb` — clústeres de cinco grupos no equidistantes y μ crítico.
- `runs/`, `results/` — registros y reportes de las corridas.

## Clasificación de resultados

Cada afirmación se etiqueta como **DEMOSTRADO** (probado o verificado con intervalos aquí), **EVIDENCIA NUMÉRICA** (reproducido con semillas y tolerancias explícitas, sin demostración) o **CITADO** (tomado de la literatura, no verificado). El preprint de O'Keeffe (arXiv:2610.21271v1) **no tiene referato por pares**. No se inventan datos ausentes: lo que no tiene respaldo se marca como `TODO` o `XFAIL`.

## Uso de IA

Código escrito por Google Antigravity. Análisis, documentación y contrastación de cifras con asistencia de Claude (Anthropic) y Gemini (Google). Historiales y material adicional en la entrega.
