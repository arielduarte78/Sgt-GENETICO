# kuramoto_templates

Marco de investigación científica auditable para el estudio y búsqueda de plantillas de clusters $(x, \phi, P)$ en el modelo de Kuramoto homogéneo sobre grafos no dirigidos.

## 1. Propósito Científico

El objetivo de este paquete es modelar y analizar plantillas de clusters según la formulación de O'Keeffe (2026, Sección III) y Canale (2022). Una plantilla $T = (x, \phi, P)$ para $m$ clusters está definida por:
- **Fracciones de osciladores (tamaños):** $x_a > 0$ con $\sum_{a=1}^m x_a = 1$.
- **Fases de cluster:** $\phi_a \in (-\pi, \pi]$.
- **Densidades de enlace:** $P_{ab} = P_{ba} \in [0, 1]$.

Bajo $(x, \phi)$ fijos, las condiciones del modelo son afines en $P$:
1. **Conectividad:** $\mu(T) = \min_a d_a$, donde $d_a = \sum_b P_{ab} x_b$.
2. **Equilibrio de fases:** $F_a = \sum_b P_{ab} x_b \sin(\phi_a - \phi_b) = 0$ para todo $a$.
3. **Estabilidad intra-cluster:** $c_a = \sum_b P_{ab} x_b \cos(\phi_a - \phi_b) > 0$ para todo $a$.
4. **Estabilidad inter-cluster:** $Q^T L_w(P) Q \succ 0$, donde $w_{ab} = x_a x_b P_{ab} \cos(\phi_a - \phi_b)$ y $Q \in \mathbb{R}^{m \times (m-1)}$ es una base ortonormal de $\mathbf{1}^\perp$ (deflación de la rotación global neutra).

---

## 2. Reglas de Honestidad Epistémica (Obligatorias)

1. **Etiquetado de Valores Citados:** Todo valor numérico proveniente de un preprint o artículo científico se etiqueta explícitamente como `CITADO (no verificado)`, citando autor, año y sección. En particular, el preprint de K. P. O'Keeffe (*Beyond 11/16: dense Kuramoto networks that resist synchrony*, arXiv:2610.21271v1, 15 de octubre de 2026) **NO cuenta con referato por pares** y declara haber sido producido mayoritariamente con agentes de IA.
2. **Taxonomía de Estados:** Toda afirmación o resultado en la documentación y logs se clasifica bajo una de tres categorías estrictas:
   - **DEMOSTRADO:** Proposiciones matemáticas probadas formalmente o verificadas por métodos de intervalo rigurosos en este repositorio.
   - **EVIDENCIA NUMÉRICA:** Resultados de optimización numérica (SDP, GA, Newton en mpmath) reproducidos localmente con tolerancias y semillas explícitas, que no constituyen una demostración analítica ni intervalo-certificada.
   - **CITADO:** Afirmaciones, cotas o valores numéricos tomados de la literatura externa aún no reproducidos ni verificados de manera independiente en nuestro entorno.
3. **Ausencia de Fabricación de Datos:** El *Supplementary Material* de O'Keeffe no se encuentra disponible en este repositorio. Queda estrictamente prohibido inventar o extrapolar números ausentes (como matrices racionales de 42 dígitos o correcciones exactas $\Delta(h)$). Todo caso de prueba o parámetro sin respaldo documental directo se implementa como `stub` con un `TODO` explícito.

---

## 3. Arquitectura y Restricciones de Dependencia

El código sigue una estricta jerarquía unidireccional de capas verificada automáticamente mediante tests de AST (`tests/unit/test_architecture.py`):

$$\text{domain} \longrightarrow \text{evaluators} \longrightarrow \text{search}$$

- `domain/` (**Dominio Matemático Puro**):
  - No depende de ningún otro módulo del paquete ni de librerías de optimización (`cvxpy`, solvers externos, generadores estocásticos).
  - Contiene las definiciones algebraicas de `Template`, firmas de residuos y formas cuadráticas, y transformaciones de gauge.
- `evaluators/` (**Evaluación de Aptitud SDP**):
  - Depende únicamente de `domain/`.
  - Define la interfaz abstracta `Evaluator` y el dataclass de resultados `EvalResult`.
  - No depende de `search/`.
- `search/` (**Motor del Algoritmo Genético**):
  - Depende de `evaluators/` y `domain/`.
  - No es importado por `evaluators/` ni por `domain/`.
- `cases/` (**Casos de Prueba con Procedencia**):
  - Depende de `domain/`.
- `io/` (**Registro Auditable**):
  - Generación de manifiestos y logs estructurados (JSONL) con metadatos de reproducibilidad.

---

## 4. Guía de Reproducción

El entorno y las herramientas se gestionan mediante `uv` y `make`:

### Creación del entorno e instalación
```bash
make env
```
Crea el entorno virtual con Python 3.12 y descarga las dependencias pineadas exactas.

### Ejecución de tests
```bash
make test
```
Ejecuta la suite de pruebas unitarias, de integración y de arquitectura con `pytest`.

### Análisis estático y tipado estricto
```bash
make lint
```
Ejecuta el linter `ruff` y el verificador estático de tipos `mypy --strict`.

### Flujo de reproducción
```bash
make reproduce
```
Ejecuta las validaciones de baseline y verificación del estado del pipeline.