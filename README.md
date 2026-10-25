# kuramoto_templates

Marco de investigación científica auditable para el estudio y búsqueda de plantillas de clusters $(x, \phi, P)$ en el modelo de Kuramoto homogéneo sobre grafos no dirigidos.

## 1. Propósito Científico

El objetivo de este paquete es modelar y analizar plantillas de clusters según la formulación de O'Keeffe (2026, Sección III) y Canale (2022). Una plantilla $T = (x, \phi, P)$ para $m$ clusters está definida por:
- **Fracciones de osciladores (tamaños):** $x_a > 0$ con $\sum_{a=1}^m x_a = 1$.
- **Fases de cluster:** $\phi_a \in (-\pi, \pi]$.
- **Densidades de enlace:** $P_{ab} = P_{ba} \in [0, 1]$.

Bajo $(x, \phi)$ fijos, las condiciones del modelo son afines en $P$:
1. **Conectividad:** $\mu(T) = \min_a d_a$, donde $d_a = x_a + \sum_{b \neq a} x_b P_{ab}$ (con $P_{aa} = 1$).
2. **Equilibrio de fases:** $F_a = \sum_{b \neq a} x_b \sin(\phi_a - \phi_b) P_{ab} = 0$ para todo $a$.
3. **Estabilidad intra-cluster:** $c_a = x_a + \sum_{b \neq a} x_b \cos(\phi_a - \phi_b) P_{ab} > 0$ para todo $a$.
4. **Estabilidad inter-cluster:** $M(p) = Q^T L_w(p) Q \succ 0$, donde $L_w(p) = \sum_{a < b} x_a x_b \cos(\phi_a - \phi_b) P_{ab} (e_a - e_b)(e_a - e_b)^T$ y $Q \in \mathbb{R}^{m \times (m-1)}$ es una base ortonormal de $\mathbf{1}^\perp$ (deflación de la rotación global neutra).

---

## 2. Reglas de Honestidad Epistémica (Obligatorias)

1. **Etiquetado de Valores Citados:** Todo valor numérico proveniente de un preprint o artículo científico se etiqueta explícitamente como `CITADO (no verificado)`, citando autor, año y sección. En particular, el preprint de K. P. O'Keeffe (*Beyond 11/16: dense Kuramoto networks that resist synchrony*, arXiv:2610.21271v1, 15 de octubre de 2026) **NO cuenta con referato por pares** y declara haber sido producido mayoritariamente con agentes de IA.
2. **Taxonomía de Estados:** Toda afirmación o resultado en la documentación y logs se clasifica bajo una de tres categorías estrictas:
   - **DEMOSTRADO:** Proposiciones matemáticas probadas formalmente o verificadas por métodos de intervalo rigurosos en este repositorio.
   - **EVIDENCIA NUMÉRICA:** Resultados de optimización numérica (SDP, GA, Newton en mpmath) reproducidos localmente con tolerancias y semillas explícitas, que no constituyen una demostración analítica ni intervalo-certificada.
   - **CITADO:** Afirmaciones, cotas o valores numéricos tomados de la literatura externa aún no reproducidos ni verificados de manera independiente en nuestro entorno.
3. **Ausencia de Fabricación de Datos:** El *Supplementary Material* de O'Keeffe no se encuentra disponible en este repositorio. Queda estrictamente prohibido inventar o extrapolar números ausentes (como matrices racionales de 42 dígitos o correcciones exactas $\Delta(h)$). Todo caso de prueba o parámetro sin respaldo documental directo se implementa como `stub` con un `TODO` explícito o se reporta como `XFAIL`.

---

## 3. Arquitectura y Restricciones de Dependencia

El código sigue una estricta jerarquía unidireccional de capas verificada automáticamente mediante tests de AST (`tests/unit/test_architecture.py`):

$$\text{domain} \longrightarrow \text{evaluators} \longrightarrow \text{search}$$

- `domain/` (**Dominio Matemático Puro**):
  - No depende de ningún otro módulo del paquete ni de librerías de optimización (`cvxpy`, solvers externos, generadores estocásticos).
  - Contiene las definiciones algebraicas de `Template`, matrices y operadores puros en numpy ($A_F, A_c, A_d, L_w, M$), y la verificación independiente `verify_solution`.
- `evaluators/` (**Evaluación de Aptitud SDP**):
  - Depende únicamente de `domain/`.
  - Implementa la clase `SDPEvaluator(Evaluator)` y el dataclass inmutable `EvalConfig`.
  - Registra diagnósticos completos en `EvalResult`: estatus del solver, versiones de cvxpy y Clarabel/SCS, tiempos de ejecución y banderas de inexactitud numérica.
- `cases/` (**Casos de Prueba con Procedencia**):
  - Depende de `domain/`.
  - Contiene las definiciones de $C_4$, $S_7$ y Canale (stub).
- `search/` (**Motor del Algoritmo Genético**):
  - Depende de `evaluators/` y `domain/`.
- `io/` (**Registro Auditable**):
  - Generación de manifiestos y logs estructurados (JSONL) con metadatos de reproducibilidad.

---

## 4. Evaluador SDP en Dos Fases y Piso Numérico

La evaluación de aptitud de un candidato $(x, \phi)$ se realiza mediante optimización semidefinida convexa (SDP) en dos fases, controladas por el dataclass `EvalConfig`:

### Fase I (Elástica, siempre factible):
$$\min_{\tau \ge 0, \, \delta \ge 0, \, p \in [0, 1]^K} \Phi = \frac{\tau}{\varepsilon \cdot \theta_{\text{tol}}} + \frac{\delta}{\varepsilon}$$
$$\text{s.a.} \quad |F_a(p)| \le \tau, \quad M(p) \succeq (\varepsilon - \delta) I_{m-1}, \quad c_a(p) \ge \varepsilon_c$$
- **Propiedad teórica:** $p = 0$ es factible con $\tau = 0$ y $\delta = \varepsilon$, garantizando $0 \le \Phi \le 1$.

### Fase II (Solo si $\Phi \le 10^{-6}$):
Elimina el residuo de equilibrio parametrizando el espacio nulo $A_F p = 0 \implies p = N z$ con $N = \text{null\_space}(A_F)$:
$$\max_{t, \, z} t \quad \text{s.a.} \quad d_a(N z) \ge t, \quad 0 \le N z \le 1, \quad M(N z) \succeq \varepsilon I_{m-1}, \quad c_a(N z) \ge \varepsilon_c$$
Devuelve la conectividad óptima $\mu^* = t$.

### Fitness Escalar $f$:
$$f = \begin{cases} \mu^* & \text{si Fase I ok, Fase II óptima y validada independientemente} \\ -\log_{10}\left(1 + \frac{\max(\Phi, 10^{-6})}{10^{-6}}\right) & \text{en otro caso (penalización suave en } [-6, 0)\text{)} \end{cases}$$

### Verificación Independiente en numpy Puro:
Todo resultado reportado como factible pasa obligatoriamente por `verify_solution`, que recalcula con numpy puro (sin cvxpy) el residuo $\|F\|_\infty$, los autovalores de $M$ mediante `scipy.linalg.eigvalsh`, $c_a$ y $\mu = \min_a d_a$.

### Piso Numérico de Resolución Espectral:
Márgenes espectrales $\varepsilon < \sim 10^{-7}$ no son resolubles con tolerancias de solver de $10^{-8} \dots 10^{-10}$ en aritmética IEEE 754 de doble precisión:
- Para $\varepsilon \ge 10^{-5}$, Clarabel converge con estatus `OPTIMAL`.
- Para $\varepsilon \approx 10^{-6} \dots 10^{-7}$, el solver entra en zona de inexactitud y reporta `OPTIMAL_INACCURATE` (registrado explícitamente en `EvalResult.solver_info["optimal_inaccurate"] = 1`).
- Para $\varepsilon \le 10^{-8}$, el cono semidefinido colapsa frente al ruido de redondeo, produciendo fallos de convergencia o soluciones no acotadas.
- La curva completa de degradación está documentada y guardada como artefacto en `results/numerical_floor_s7.json` y `results/numerical_floor_s7.md`.
- **Regla de interpretación física:** En la Fase I ($\max \lambda - \rho s$), $s^* \approx 0$ **NO implica estabilidad**; solo $\lambda^* \ge \text{piso}$ la implica. El slack $s^* \approx 0$ puede coexistir con $\lambda^* < 0$ cuando el solver compensa grado a costa de estabilidad espectral.

### 4.1 Heurísticas de Optimización: Barrera de Régimen y Piso Numérico

En el evaluador escalado (`ScaledSDPEvaluator`), el paisaje de búsqueda incorpora dos constantes heurísticas fundamentales para guiar la optimización genética:

#### 1. Barrera de Régimen ($-10.0$) para Estados `NOT_SCALABLE` [HEURÍSTICA]:
- **Justificación técnica:**
  - En la formulación asintótica escalada, la aptitud se define sobre el margen normalizado $\hat{\lambda} = \lambda / h_{\text{eff}}^2$.
  - Los candidatos degenerados o colapsados ($h_{\text{eff}} < 10^{-4}$, $m < 5$ o particiones sin hojas resolubles) no pueden evaluarse en el cono acoplado $6 \times 6$ y deben ser procesados por el evaluador plano.
  - En el evaluador plano, la cota espectral residual $\lambda_{\text{plain}} \sim -10^{-10}$ produce una penalización aparente $f_{\text{plain}} \approx -10^{-5}$. En contraste, un candidato escalable prometedor en la frontera crítica ($\mu_t > 0.691$, $h_{\text{eff}} \sim 2 \cdot 10^{-3}$) con $\hat{\lambda} \sim -10^{-3}$ recibe $f \sim -10^{-3}$.
  - Sin una barrera divisoria, $-10^{-5} > -10^{-3}$, lo que desata un fenómeno severo de *Reward Hacking*: la selección por torneo extingue a la población escalable en favor de configuraciones triviales no escalables como $C_4$ ($m=4, h=0$).
  - La barrera penaliza todo estado no escalable con $f = -10.0 + \dots$, asegurando una separación topológica estricta: cualquier individuo escalable ($f \in [-1, +\infty)$) domina a cualquier candidato colapsado ($f \le -10.0$).
- **Limitaciones:**
  - **Discontinuidad artificial:** Introduce un escalón pronunciado en el paisaje de aptitud que no refleja una propiedad física continua del modelo de Kuramoto.
  - **Pérdida de gradiente:** Por debajo de $-10.0$, el algoritmo genético pierde sensibilidad direccional para guiar la reconstitución de agrupamientos válidos, dependiendo puramente de operadores de mutación y cruce estocásticos.
  - **Arbitrariedad de escala:** El valor $-10.0$ es una convención arbitraria; cualquier valor inferior a $-2.0$ satisface el requerimiento de separación.

#### 2. Piso Espectral $\hat{\lambda}_{\text{floor\_hat}} = 10^{-5}$ [HEURÍSTICA]:
- **Justificación técnica:**
  - Distingue la estabilidad asintótica física genuina de las fluctuaciones numéricas de punto interior de los solvers primarios y secundarios.
  - En la calibración empírica sobre 200 individuos independientes contrastando `Clarabel(1e-8)`, `Clarabel(1e-9)` y `SCS(1e-8)`:
    1. El desacuerdo interno de punto interior de Clarabel no supera $7.73 \times 10^{-8}$.
    2. El desacuerdo máximo con SCS cuando ambos alcanzan `OPTIMAL` es de $5.04 \times 10^{-6}$.
    3. Las divergencias mayores a $10^{-4}$ se deben exclusivamente a fallos de convergencia de SCS (`OPTIMAL_INACCURATE`).
  - Fijar $\hat{\lambda}_{\text{floor\_hat}} = 10^{-5}$ provee un margen de seguridad de $>100\times$ respecto al ruido interno de Clarabel y $>2\times$ respecto a SCS convergido, garantizando además que plantillas físicas como $S_7$ ($\hat{\lambda} \approx 3.2 \times 10^{-4}$ para $h=0.002$) sean aceptadas inequívocamente.
- **Limitaciones:**
  - **Falsos negativos marginales:** Descarta candidatos teóricamente estables cuya curvatura espectral normalizada se encuentre en el intervalo marginal $(0, 10^{-5})$, si bien tales configuraciones carecen de robustez física frente a perturbaciones.
  - **Dependencia de la aritmética:** El umbral $10^{-5}$ está calibrado para aritmética IEEE 754 float64; en aritmética de precisión arbitraria o variable (p. ej. `mpmath`), un piso analítico inferior podría ser resoluble.
  - **Condicionamiento en $h$ pequeño:** Para $h_{\text{eff}} \ll 10^{-3}$, el margen absoluto en la matriz original es $\lambda_{\min}(M) \sim \hat{\lambda} \cdot h_{\text{eff}}^2 \approx 10^{-11}$, aproximándose al umbral de sensibilidad de las descomposiciones espectrales en doble precisión.

---

## 5. Guía de Reproducción

El entorno y los experimentos se gestionan mediante `uv` y `make`:

### Creación del entorno e instalación
```bash
make env
```

### Ejecución de tests
```bash
make test
```
Ejecuta la suite completa de pruebas unitarias, de integración, regresión citada y arquitectura.

### Análisis estático y tipado estricto
```bash
make lint
```
Ejecuta `ruff check` y `mypy --strict`.

### Reporte auditable del evaluador
```bash
make reproduce-eval
```
Ejecuta el protocolo completo de evaluación, verifica identidades de dominio, contrastes Clarabel vs SCS, regresión contra literatura citada y genera el reporte JSON en `results/reproduce_eval_report.json`.