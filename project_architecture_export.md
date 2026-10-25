# Dossier Arquitectónico y Físico-Matemático: Sistema Sgt-GENETICO
**Framework:** `kuramoto_templates`  
**Rol:** Ingeniero de Software Principal / Arquitecto de Sistemas  
**Destinatario:** Auditor Externo de Inteligencia Artificial  
**Estado:** Documento de Especificación Técnica Integral (Matemática, Código y Arquitectura)

---

## 1. Arquitectura de 3 Capas y Jerarquía Estricta

El sistema `Sgt-GENETICO` está diseñado bajo un principio estricto de separación de responsabilidades y jerarquía unidireccional de dependencias:

$$\mathbf{domain} \longrightarrow \mathbf{evaluators} \longrightarrow \mathbf{search}$$

Esta jerarquía no es meramente conceptual: se encuentra formalizada y garantizada mecánicamente mediante tests de Árbol de Sintaxis Abstracta (AST) en la suite de pruebas unitarias (`tests/unit/test_architecture.py`). Cualquier importación inversa (por ejemplo, que `domain` importe algo de `evaluators`, o `evaluators` de `search`) rompe la compilación y los tests de CI.

```
+-------------------------------------------------------------------------------+
|                                    search/                                    |
|   - Algoritmo Genético de Islas con migración en anillo (GAEngine)            |
|   - Genotipo híbrido discreto-continuo (Genome) y reducción de Gauge          |
|   - Operadores de cruce, mutación gaussiana, salto angular, Split & Merge     |
|   - Selección por Torneo k y Deterministic Crowding en espacio fenotípico     |
+---------------------------------------+---------------------------------------+
                                        | (invoca evaluación de aptitud)
                                        v
+-------------------------------------------------------------------------------+
|                                  evaluators/                                  |
|   - ScaledSDPEvaluator (Lupa de Sylvester, reescalado asintótico O(1))        |
|   - SDPEvaluator plano (Formulación estándar de O'Keeffe)                     |
|   - Protocolo en Dos Fases (Fase I: Slack elástico; Fase II: Maximización)    |
|   - Antidoping Robusto (Clarabel 1e-9, NumPy float64, mpmath 50d, SCS)        |
+---------------------------------------+---------------------------------------+
                                        | (utiliza definiciones algebraicas puras)
                                        v
+-------------------------------------------------------------------------------+
|                                    domain/                                    |
|   - Plantilla inmutable Template (x, phi, P)                                  |
|   - Operadores algebraicos lineales puros en NumPy: A_F, A_c, A_d             |
|   - Formas bilineales cuadráticas: L_w(p) y M(p) = Q^T L_w(p) Q               |
|   - Verificador físico independiente de condiciones de Kuramoto (sin cvxpy)  |
+-------------------------------------------------------------------------------+
```

### 1.1. Capa `domain/` (Dominio Matemático Puro)
- **Aislamiento absoluto:** No depende de ningún otro módulo del paquete ni de frameworks de optimización matemática (`cvxpy`), solvers numéricos externos (`Clarabel`, `SCS`) ni generadores estocásticos.
- **Estructura inmutable:** Define la estructura de datos `Template(x, phi, P)` con validación de invariantes fundamentales en su `__post_init__` (positividad de tamaños $x_a > 0$, suma $\sum x_a = 1$, fases $\phi_a \in (-\pi, \pi]$, simetría de densidades $P_{ab} = P_{ba} \in [0, 1]$ con diagonal unidad $P_{aa} = 1$).
- **Operadores matriciales en NumPy:** Construye de forma analítica y determinista las matrices de incidencia y acoplamiento:
  - $A_F$: Matriz de equilibrio de fases ($m \times K$).
  - $A_c$: Matriz de estabilidad intra-cluster ($m \times K$).
  - $A_d$: Matriz de conectividad/grados de cluster ($m \times K$).
  - $L_w(p)$: Matriz Laplaciana ponderada de Kuramoto ($m \times m$).
  - $M(p) = Q^T L_w(p) Q$: Proyección de $L_w(p)$ sobre el subespacio ortogonal al vector de unos $\mathbf{1}_m^\perp$ mediante una base ortonormal $Q \in \mathbb{R}^{m \times (m-1)}$, deflando la simetría de rotación global.
- **Verificación pura:** Contiene `verify_solution`, función que audita cualquier triplete $(x, \phi, P)$ evaluando de forma aislada $\|F\|_\infty$, $\min_a c_a$, $\min_a d_a$ y $\lambda_{\min}(M)$ mediante `scipy.linalg.eigvalsh`.

### 1.2. Capa `evaluators/` (Optimización Semidefinida y Lupa de Sylvester)
- **Responsabilidad:** Resolver el problema inverso de grafos: dados los tamaños de cluster $x$ y las fases $\phi$, determinar si existe una matriz de densidades de enlace $P \in [0, 1]^{m \times m}$ que satisfaga el equilibrio dinámico, la estabilidad espectral y una conectividad mínima fijada $\mu_t$, optimizando dicha estructura.
- **`SDPEvaluator`:** Implementa la formulación semidefinida convexa (SDP) estándar de dos fases sobre la proyección física $M(p)$.
- **`ScaledSDPEvaluator`:** Resuelve la patología de escala cuando clusters macroscópicos se dividen en hojas separadas por ángulos infinitesimales $h \ll 1$. Introduce una base reducida no ortonormal $T$ y el operador reescalado $\hat{G}(p) = T^T L_w(p) T / h_{\text{eff}}^2$, gobernado por la **Ley de Inercia de Sylvester**.
- **Protocolo en Dos Fases:**
  - *Fase I (Elástica sobre espacio nulo):* Garantiza equilibrio exacto $\|A_F p\|_\infty \le 10^{-14}$ parametrizando $p = N z$ con $N = \text{null\_space}(A_F)$. Maximiza $\lambda - \rho s$ (o $\hat{\lambda} - \rho s$), donde $s$ mide el déficit de conectividad respecto a $\mu_t$.
  - *Fase II (Maximización estricta):* Si $s^* \le \text{tol}$, congela el margen espectral sobre el piso numérico $\hat{\lambda} \ge \hat{\lambda}_{\text{floor\_hat}}$ y maximiza directamente el grado mínimo $\mu^* = \min_a d_a$.
- **Antidoping Robusto:** Protocolo de cuádruple verificación cruzada para descartar artefactos numéricos del solver antes de admitir cualquier individuo en el *Hall of Fame*: Clarabel con tolerancia $10^{-9}$, NumPy float64 puro, `mpmath` a 50 dígitos de precisión decimal (`mp.eigsy`) y SCS con 20,000 iteraciones.

### 1.3. Capa `search/` (Motor Evolutivo por Islas)
- **Genotipo compacto:** Representación de genoma $g = (m, K_{\text{ref}}, k, \text{sgn}, \ell, u)$ que desacopla la topología combinatoria (anclas $k$, signos $\text{sgn}$) de la separación angular continua en escala logarítmica ($\ell = \log_{10}|\Delta\phi| \in [-6, -0.5]$) y los logits de masa $u$.
- **Reducción de Gauge Canónico:** Transforma cualquier genoma a su representante canónico único en su órbita de simetría (rotación global para anclar el cluster de mayor masa en $k=0$ y reflexión especular $\phi \to -\phi$ eligiendo el mínimo lexicográfico).
- **Topología de Islas y Migración en Anillo:** Divide la población en múltiples islas independientes que evolucionan en paralelo. Periódicamente se ejecuta una migración en anillo unidireccional de los mejores individuos para transferir innovaciones estructurales sin homogenizar prematuramente el acervo genético.
- **Reemplazo por Deterministic Crowding:** Los descendientes compiten directamente contra sus progenitores más cercanos en distancia fenotípica circular. Esto mantiene múltiples nichos ecológicos activos y evita que un atractor local colapse la diversidad poblacional.

---

## 2. El Modelo Físico-Matemático (Las Ecuaciones)

### 2.1. Dinámica de Kuramoto en Grafos y Reducción por Plantillas

El modelo clásico de Kuramoto describe $N$ osciladores de fase acoplados sobre un grafo no dirigido simple con matriz de adyacencia binaria $A_{ij} \in \{0, 1\}$ y frecuencias naturales idénticas (asumidas nulas bajo el sistema de referencia rotante):

$$\dot{\theta}_i = \sum_{j=1}^N A_{ij} \sin(\theta_j - \theta_i), \quad i = 1, \dots, N$$

Un estado de sincronización de clusters consiste en una partición de los $N$ osciladores en $m$ grupos disjuntos $C_1, \dots, C_m$, donde los osciladores en el cluster $a$ tienen una masa normalizada $x_a = |C_a|/N > 0$ con $\sum_{a=1}^m x_a = 1$, y comparten una fase común $\theta_i = \phi_a$ para todo $i \in C_a$.

La interacción entre clusters queda gobernada por la matriz de densidades de enlace $P \in [0, 1]^{m \times m}$, donde $P_{ab} = P_{ba}$ denota la fracción de enlaces presentes entre el cluster $a$ y el cluster $b$ respecto al total de posibles enlaces $N^2 x_a x_b$. Para enlaces intra-cluster, se asume densidad completa $P_{aa} = 1$.

En el límite termodinámico ($N \to \infty$), el comportamiento microscópico colapsa de forma exacta en el triplete macroscópico $T = (x, \phi, P)$.

### 2.2. Condiciones de Equilibrio y Estabilidad Dinámica

Para que una plantilla $T = (x, \phi, P)$ represente un estado de equilibrio dinámico linealmente estable del grafo subyacente, debe satisfacer simultáneamente 4 condiciones analíticas:

#### 1. Conectividad Normalizada (Grados de Cluster):
El grado asintótico de los osciladores del cluster $a$ viene dado por la combinación afín:
$$d_a(p) = x_a P_{aa} + \sum_{b \neq a} x_b P_{ab} = x_a + \sum_{b \neq a} x_b P_{ab}$$
La conectividad global del grafo se define como el grado mínimo:
$$\mu(T) = \min_{a \in \{1, \dots, m\}} d_a(p)$$
El objetivo del marco es encontrar configuraciones estables que resistan la sincronización global alcanzando $\mu(T) > \frac{11}{16} = 0.6875$.

#### 2. Ecuación de Equilibrio de Fases:
Cada cluster debe tener torque neto nulo:
$$F_a(p) = \sum_{b \neq a} x_b \sin(\phi_a - \phi_b) P_{ab} = 0, \quad \forall a \in \{1, \dots, m\}$$
En forma vectorial compacta:
$$A_F p = \mathbf{0}$$
donde $p \in [0, 1]^K$ es el vector que apila las $K = \frac{m(m-1)}{2}$ entradas fuera de la diagonal de $P$, y $A_F \in \mathbb{R}^{m \times K}$ tiene entradas no nulas para el par $k = (a, b)$ con $a < b$:
$$(A_F)_{a, k} = x_b \sin(\phi_a - \phi_b), \quad (A_F)_{b, k} = -x_a \sin(\phi_a - \phi_b)$$

#### 3. Estabilidad Intra-Cluster (Modos Internos):
Para que osciladores con fases idénticas dentro del mismo cluster no se dispersen frente a perturbaciones microscópicas, el acoplamiento restaurador medio debe ser estrictamente positivo:
$$c_a(p) = x_a + \sum_{b \neq a} x_b \cos(\phi_a - \phi_b) P_{ab} \ge \varepsilon_c > 0, \quad \forall a \in \{1, \dots, m\}$$
En forma matricial:
$$x + A_c p \ge \varepsilon_c \mathbf{1}_m$$
donde $(A_c)_{a, k} = x_b \cos(\phi_a - \phi_b)$ y $(A_c)_{b, k} = x_a \cos(\phi_a - \phi_b)$.

#### 4. Estabilidad Inter-Cluster (Modos Colectivos):
La estabilidad del movimiento relativo entre clusters está gobernada por la matriz Hessiana/Laplaciana ponderada $L_w(p) \in \mathbb{R}^{m \times m}$:
$$L_w(p) = \sum_{a < b} x_a x_b \cos(\phi_a - \phi_b) P_{ab} (e_a - e_b)(e_a - e_b)^T$$
Dado que $\sum_b (L_w)_{ab} = 0$, el vector $\mathbf{1}_m$ siempre pertenece al núcleo de $L_w$ con autovalor cero, correspondiente a la simetría de norma global (rotación simultánea de todas las fases $\theta_i \to \theta_i + \alpha$).

Para desacoplar este modo neutro, se proyecta $L_w$ sobre $\mathbf{1}_m^\perp$ mediante una base ortonormal $Q \in \mathbb{R}^{m \times (m-1)}$ tal que $Q^T Q = I_{m-1}$ y $Q^T \mathbf{1}_m = \mathbf{0}$. La matriz proyectada es:
$$M(p) = Q^T L_w(p) Q \in \mathbb{R}^{(m-1) \times (m-1)}$$
La condición de estabilidad lineal exige que $M(p)$ sea estrictamente definida positiva:
$$M(p) \succeq \lambda_{\text{floor}} I_{m-1}, \quad \text{con } \lambda_{\text{floor}} > 0$$

---

### 2.3. La Ceguera Numérica del Solver y la Lupa de Sylvester

#### La Patología de Escala Asintótica ($h \to 0$):
En las construcciones teóricas de alta densidad (como la familia $S_7$ de O'Keeffe), los estados estables se obtienen como deformaciones perturbativas de configuraciones macroscópicas altamente simétricas (como $C_4$, el estado de 4 clusters ortogonales en $0, \pi/2, \pi, 3\pi/2$).

Para romper la barrera de conectividad $\mu = 11/16$, los clusters macroscópicos se dividen en sub-clusters ("hojas" o *sheets*) separados por un ángulo de splitting infinitesimal $2h$:
$$\phi_{a^\pm} = \phi_a \pm h, \quad \text{con } h \sim 10^{-3} \dots 10^{-4}$$

Cuando $h \ll 1$, las diferencias angulares internas satisfacen:
$$\cos(\phi_{a^+} - \phi_{a^-}) = \cos(2h) = 1 - 2h^2 + \mathcal{O}(h^4)$$
$$\sin(\phi_{a^+} - \phi_{a^-}) = \sin(2h) = 2h + \mathcal{O}(h^3)$$

Al calcular el Laplaciano proyectado físico $M(p)$, la matriz exhibe dos escalas espectrales radicalmente dispares:
1. **Modos macroscópicos inter-hojas:** $\lambda \sim \mathcal{O}(1)$.
2. **Modos microscópicos de splitting intra-hoja:** $\lambda \sim \mathcal{O}(h^2)$.

Para $h = 2 \cdot 10^{-3}$, el autovalor mínimo físico es $\lambda_{\min}(M) \sim h^2 \approx 4 \cdot 10^{-6}$. Para $h = 5 \cdot 10^{-4}$, $\lambda_{\min}(M) \sim 2.5 \cdot 10^{-7}$.

#### El Colapso Numérico en Punto Interior:
Los solvers de optimización cónica de precisión estándar en float64 (Clarabel, SCS) operan con tolerancias numéricas de convergencia entre $10^{-8}$ y $10^{-9}$. Cuando el solver evalúa la restricción de cono semidefinido:
$$M(p) - \lambda I \succeq 0$$
para $\lambda \sim 10^{-7}$, la distancia al cono cae por debajo de la resolución de los residuos de primalidad y dualidad. El solver entra en **ceguera numérica**:
- Se vuelve incapaz de discernir entre curvatura física estable y ruido de punto interior.
- Emite estatus `OPTIMAL_INACCURATE` o `SOLVER_ERROR`.
- Peor aún: el gradiente de optimización se aplana, provocando *Reward Hacking*, donde el algoritmo genético extingue poblaciones prometedoras al confundir perturbaciones microscópicas legítimas con estados singulares.

#### La Lupa de Sylvester (Sylvester's Magnifier):
Para curar la ceguera numérica sin incurrir en el costo computacional prohibitivo de los solvers de precisión arbitraria, implementamos un reescalado analítico basado en el **Teorema de Inercia de Sylvester**.

> **Teorema (Ley de Inercia de Sylvester):**  
> Sean $A \in \mathbb{R}^{n \times n}$ una matriz simétrica real y $S \in \mathbb{R}^{n \times n}$ una matriz no singular ($\det(S) \neq 0$). Entonces, la matriz congruente:
> $$B = S^T A S$$
> posee exactamente la misma inercia que $A$, es decir, el mismo número de autovalores estrictamente positivos, estrictamente negativos y nulos. En particular:
> $$A \succ 0 \iff S^T A S \succ 0$$

#### Construcción de la Escala Efectiva $h_{\text{eff}}$:
Los clusters se particionan en hojas $C_1, \dots, C_{n_s}$ según su índice de ancla discreta $k_a$ del genoma. Para cada hoja $j$, se calcula la media circular de las fases:
$$\bar{\phi}_j = \text{atan2}\left(\sum_{i \in C_j} x_i \sin\phi_i, \;\; \sum_{i \in C_j} x_i \cos\phi_i\right)$$
La escala efectiva $h_{\text{eff}}$ se define rigurosamente como la raíz cuadrática ponderada por masa (RMS) de los desvíos angulares respecto a la media de su respectiva hoja:
$$h_{\text{eff}} = \sqrt{\sum_{j=1}^{n_s} \sum_{i \in C_j} x_i \left( \text{wrap}_{(-\pi, \pi]}(\phi_i - \bar{\phi}_j) \right)^2}$$
- **Rango admisible:** Se exige $h_{\text{eff}} \in [10^{-4}, 0.1]$. Si $h_{\text{eff}} < 10^{-4}$, el individuo se cataloga como `NOT_SCALABLE` y se le aplica la barrera de régimen.

#### Construcción de la Base Reducida $T \in \mathbb{R}^{m \times (m-1)}$:
La matriz $T$ se construye ensamblando dos bloques de columnas que desacoplan el espacio tangente:
1. **Columnas Macro ($\eta$):** Para cada una de las $n_s - 1$ hojas que no son la hoja de referencia $s_{\text{ref}}$ (la hoja de mayor masa), se define un vector indicador:
   $$\eta_j = \sum_{a \in C_j} e_a \in \mathbb{R}^m$$
2. **Columnas Micro ($\zeta$):** Para cada hoja que contiene $n_j \ge 2$ clusters, se calcula una base ortonormal del espacio nulo de $\mathbf{1}_{n_j}^T$ mediante factorización QR: $Z_j \in \mathbb{R}^{n_j \times (n_j - 1)}$, con $Z_j^T \mathbf{1}_{n_j} = \mathbf{0}$ y $Z_j^T Z_j = I$. Cada una de estas direcciones describe el desbalance relativo interno entre osciladores de la misma hoja, y **se multiplica analíticamente por la escala $h_{\text{eff}}$**:
   $$\zeta_{j, \ell} = h_{\text{eff}} \cdot \tilde{Z}_{j, \ell}$$
   donde $\tilde{Z}_{j, \ell} \in \mathbb{R}^m$ contiene la columna $\ell$ de $Z_j$ en las posiciones de $C_j$ y ceros en el resto.

La matriz resultante $T = [\eta \mid \zeta] \in \mathbb{R}^{m \times (m-1)}$ satisface estrictamente:
$$\text{rank}([T, \mathbf{1}_m]) = m \implies \text{span}(T) \oplus \text{span}(\mathbf{1}_m) = \mathbb{R}^m$$

#### El Operador Escalado $\hat{G}(p)$ y la Normalización $\hat{\lambda}$:
Definimos el operador matricial reescalado:
$$\hat{G}(p) = \frac{T^T L_w(p) T}{h_{\text{eff}}^2}$$
- **Mecanismo de magnificación:** Las direcciones microscópicas $\zeta$ aportan términos cuadráticos proporcionales a $h_{\text{eff}}^2$ en el numerador $T^T L_w T$. Al dividir exactamente por $h_{\text{eff}}^2$, la curvatura espectral intra-hoja se amplifica a orden $\mathcal{O}(1)$.
- Por la Ley de Inercia de Sylvester, como $\text{span}(T) = \mathbf{1}_m^\perp$, se cumple de forma idéntica:
  $$M(p) \succ 0 \iff \hat{G}(p) \succ 0$$
- El margen de estabilidad se mide sobre $\hat{\lambda} = \lambda_{\min}(\hat{G})$. En lugar de lidiar con valores numéricos microscópicos de $10^{-7}$, el solver opera confortablemente con $\hat{\lambda} \sim 10^{-3} \dots 10^{-1}$.
- El piso espectral se fija en $\hat{\lambda}_{\text{floor\_hat}} = 10^{-5}$, lo que proporciona un margen de seguridad de más de $100\times$ sobre el ruido interno de Clarabel y valida inequívocamente plantillas físicas como $S_7$ ($\hat{\lambda} \approx 3.2 \cdot 10^{-4}$ para $h=0.002$).

---

## 3. El Ciclo Evolutivo y la Función de Fitness

### 3.1. Genotipo, Decodificación Fenotípica y Canonicalización

El espacio de búsqueda continuo del modelo de Kuramoto está plagado de singularidades combinatorias y de simetrías continuas. Para estructurarlo eficientemente, el genoma parametriza la plantilla mediante variables discretas y continuas acopladas:

$$g = (m, K_{\text{ref}}, k, \text{sgn}, \ell, u)$$

| Variable | Tipo / Dominio | Rol en el Modelo Físico |
|---|---|---|
| $m$ | Entero $\in [2, 10]$ | Número de clusters del fenotipo. |
| $K_{\text{ref}}$ | Entero $\in \{4, 8, 12, 16\}$ | Número de anclas del retículo de referencia discreto en el círculo. |
| $k$ | Tupla de enteros $\in \{0, \dots, K_{\text{ref}}-1\}^m$ | Índice del ancla base para cada cluster. |
| $\text{sgn}$ | Tupla $\in \{-1, 0, 1\}^m$ | Dirección del desvío respecto al ancla (-1: antihorario, 0: exacto, +1: horario). |
| $\ell$ | Tupla flotante $\in [-6.0, -0.5]^m$ | Escala logarítmica del desvío angular: $|\Delta\phi| = 10^\ell \in [10^{-6}, 0.316]$. |
| $u$ | Tupla flotante $\in \mathbb{R}^m$ | Logits continuos sin restricciones que definen las masas $x_a$. |

#### Decodificación Fenotípica:
1. **Tamaños de cluster $x \in \Delta^{m-1}$ (Símplex con piso estricto):**
   $$x_a = x_{\min} + (1 - m \cdot x_{\min}) \cdot \frac{e^{u_a - \max(u)}}{\sum_{b=1}^m e^{u_b - \max(u)}}, \quad x_{\min} = 10^{-3}$$
   Garantiza algebraicamente que $x_a \ge 10^{-3} > 0$ y $\sum_{a=1}^m x_a = 1.0$ de manera exacta.
2. **Fases de cluster $\phi \in (-\pi, \pi]$:**
   $$\phi_a = \text{wrap}_{(-\pi, \pi]}\left( \frac{2\pi (k_a \bmod K_{\text{ref}})}{K_{\text{ref}}} + \text{sgn}_a \cdot 10^{\ell_a} \right)$$

#### Reducción Canónica de Gauge (`canonicalize`):
El modelo de Kuramoto posee dos simetrías continuas/discretas fundamentales:
1. Invarianza bajo rotación de fase global: $\theta_i \to \theta_i + \alpha$.
2. Invarianza bajo reflexión quiral: $\theta_i \to -\theta_i$.

Para evitar que el algoritmo genético explore copias idénticas trasladadas o reflejadas, todo individuo generado se reduce a su forma canónica normalizada:
1. Decodifica $(x, \phi)$ y ordena los clusters monótonamente por fase en $(-\pi, \pi]$.
2. Identifica el cluster de mayor tamaño $x_{\max}$ y calcula la traslación de anclas requerida para que su ancla quede en $k=0$.
3. Construye el genoma reflejado conjugado ($k \to -k \bmod K_{\text{ref}}$, $\text{sgn} \to -\text{sgn}$).
4. Compara lexicográficamente la representación hacia adelante y la reflejada, seleccionando de forma determinista la menor de ambas.

---

### 3.2. Formulación SDP en Dos Fases

Dado el fenotipo $(x, \phi)$, el evaluador determina $P$ resolviendo optimización semidefinida convexa.

#### Fase I: Problema Elástico sobre el Espacio Nulo de Equilibrio
Para evitar que el solver introduzca residuos de desequilibrio dinámico, se calcula la base ortonormal del núcleo de $A_F$ mediante descomposición en valores singulares:
$$A_F N = \mathbf{0}, \quad N \in \mathbb{R}^{K \times d}, \quad d = \dim(\ker A_F)$$
Se aplica refinamiento iterativo por mínimos cuadrados garantizando $\|A_F N\|_\infty \le 10^{-14}$. Las densidades de enlace quedan exactamente parametrizadas por $p = N z$, con $z \in \mathbb{R}^d$.

El problema de optimización en cvxpy formula:
$$\max_{z \in \mathbb{R}^d, \; \hat{\lambda} \in \mathbb{R}, \; s \ge 0} \left( \hat{\lambda} - \rho \cdot s \right)$$
sujeto a:
$$\begin{aligned}
0 \le N z &\le 1 \\
\sum_{k=1}^K (N z)_k G_k &\succeq \hat{\lambda} I_{m-1} \\
x + A_c (N z) &\ge \varepsilon_c \mathbf{1}_m \\
x + A_d (N z) &\ge (\mu_t - s) \mathbf{1}_m \\
s &\ge 0
\end{aligned}$$
donde $\rho = 1.0$, $\mu_t = 0.6929077$, $\varepsilon_c = 10^{-4}$ y $G_k$ son las matrices precalculadas del operador reescalado $\hat{G}$.

- **Significado del slack $s$:** Mide el déficit de grado respecto al objetivo $\mu_t$. Si $s^* \le \text{s\_tol} = 10^{-6}$, el candidato alcanza o supera la conectividad crítica $\mu_t$.

#### Fase II: Maximización de Conectividad Pura
Si la Fase I converge con $s^* \le 10^{-6}$ y $\hat{\lambda}^* \ge \hat{\lambda}_{\text{floor\_hat}}$, se fija el piso espectral y se maximiza el umbral de grado $t$:
$$\max_{z \in \mathbb{R}^d, \; t \in \mathbb{R}} t$$
sujeto a:
$$\begin{aligned}
x + A_d (N z) &\ge t \mathbf{1}_m \\
0 \le N z &\le 1 \\
\sum_{k=1}^K (N z)_k G_k &\succeq \hat{\lambda}_{\text{floor\_hat}} I_{m-1} \\
x + A_c (N z) &\ge \varepsilon_c \mathbf{1}_m
\end{aligned}$$
La solución entrega la conectividad óptima certificada $\mu^* = t^*$.

---

### 3.3. Función de Fitness Continua y Jerarquía Estricta

La asignación de aptitud $f$ y aptitud de selección $f_{\text{sel}}$ descarta por completo funciones con pesos artificiales (como penalizaciones logarítmicas o constantes arbitrarias $10^{10}$) que deforman el gradiente. Se estructura en una **jerarquía estricta de 3 niveles físicos** más una **barrera de régimen**:

```
Paisaje de Aptitud f:
====================================================================================
 f > 0.0          NIVEL 1: FACTIBLE ESTABLE           f = +lam_hat
                  (lam_hat >= 1e-5 Y s* <= 1e-6)     (Compite por maximizar estabilidad)
------------------------------------------------------------------------------------
 0.0 >= f >= -1.0 NIVEL 2: ESTABLE INCOMPLETO         f = -s*
                  (lam_hat >= 1e-5 Y s* > 1e-6)      (Déficit de conectividad, s* in (0, 1])
------------------------------------------------------------------------------------
 f < -1.0         NIVEL 3: INESTABLE FÍSICO           f = -1.0 + clip(lam_hat, -1, 0)
                  (lam_hat < 1e-5)                   (SIEMPRE peor que cualquier estable)
====================================================================================
 f <= -10.0       BARRERA DE RÉGIMEN: NOT_SCALABLE    f = -10.0 - s_plain
                  (h_eff < 1e-4 o rango colapsado)   (Protege contra Reward Hacking de C4)
====================================================================================
```

#### Regla de Dominancia Absoluta:
1. **Nivel 1 (Factible y Estable):** $\hat{\lambda} \ge \hat{\lambda}_{\text{floor\_hat}}$ y $s^* \le \text{tol}$.  
   $$f = +\hat{\lambda} > 0$$
   El individuo satisface todas las restricciones físicas y supera el grado $\mu_t$.
2. **Nivel 2 (Físicamente Estable, falta conectividad):** $\hat{\lambda} \ge \hat{\lambda}_{\text{floor\_hat}}$ y $s^* > \text{tol}$.  
   $$f = -s^* \in [-1.0, 0.0)$$
   El individuo es un estado de Kuramoto físicamente estable, pero su grado no alcanza $\mu_t$. Su distancia al objetivo es exactamente $s^*$.
3. **Nivel 3 (Físicamente Inestable):** $\hat{\lambda} < \hat{\lambda}_{\text{floor\_hat}}$.  
   $$f = -1.0 + \text{clip}(\hat{\lambda}, -1.0, 0.0) < -1.0$$
   **Teorema de Selección:** Dado que para todo individuo del Nivel 2 se cumple $s^* \le 1.0 \implies f \ge -1.0$, **ningún individuo físicamente inestable puede superar jamás a un individuo estable**, independientemente de cuán alta sea su conectividad o si tiene $s^* = 0$. Un individuo con $s^* = 0$ pero inestable ($\hat{\lambda} < 10^{-5}$) recibe $f \le -1.0$, siendo derrotado en torneo por un individuo estable con déficit $s^* = 0.5$ ($f = -0.5$).
4. **Barrera de Régimen ($-10.0$) para `NOT_SCALABLE`:**  
   Si $h_{\text{eff}} < 10^{-4}$ o $m < 5$, el individuo no puede evaluarse en el cono escalado. En el evaluador plano, la cota residual daría $f_{\text{plain}} \approx -10^{-5}$. Un candidato escalable legítimo en la frontera crítica tendría $f \sim -10^{-3}$. Sin la barrera, $-10^{-5} > -10^{-3}$, lo que causaría la extinción de las poblaciones deformadas frente a plantillas rígidas triviales como $C_4$. La barrera impone $f \le -10.0$, garantizando separación topológica.

---

## 4. Extractos de Código Crítico

A continuación se transcriben textualmente los bloques algorítmicos centrales del sistema, omitiendo comentarios periféricos y boilerplate.

### 4.1. Decodificación del Genoma a Fenotipo (`search/genome.py`)

```python
def decode_x(u: tuple[float, ...], xmin: float = X_MIN) -> tuple[float, ...]:
    m = len(u)
    if m < 2 or m > 10:
        raise ValueError(f"m debe estar en [2, 10], recibido m={m}.")
    if xmin <= 0.0 or m * xmin >= 1.0:
        raise ValueError(f"xmin={xmin} inválido para m={m}.")

    max_u = max(u)
    exp_u = [math.exp(v - max_u) for v in u]
    sum_exp = sum(exp_u)
    softmax_u = [v / sum_exp for v in exp_u]

    slack = 1.0 - m * xmin
    raw_x = [xmin + slack * sm for sm in softmax_u]
    total = sum(raw_x)
    return tuple(float(v / total) for v in raw_x)


def decode_phi(
    K_ref: int,
    k: tuple[int, ...],
    sgn: tuple[int, ...],
    ell: tuple[float, ...],
) -> tuple[float, ...]:
    m = len(k)
    two_pi = 2.0 * math.pi
    phases: list[float] = []
    for a in range(m):
        base_angle = two_pi * (k[a] % K_ref) / float(K_ref)
        offset = float(sgn[a]) * (10.0 ** ell[a])
        phases.append(wrap_minus_pi_to_pi(base_angle + offset))
    return tuple(phases)


def to_phenotype(genome: Genome) -> tuple[tuple[float, ...], tuple[float, ...]]:
    x = decode_x(genome.u, xmin=X_MIN)
    phi = decode_phi(genome.K_ref, genome.k, genome.sgn, genome.ell)
    return x, phi
```

---

### 4.2. Construcción de la Escala y Matriz Escalada (`evaluators/scaled.py`)

```python
def compute_h_eff(
    x: tuple[float, ...] | np.ndarray,
    phi: tuple[float, ...] | np.ndarray,
    sheets: dict[int, list[int]],
) -> float:
    x_arr = np.asarray(x, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    two_pi = 2.0 * math.pi

    sum_sq_devs = 0.0
    for idxs in sheets.values():
        if len(idxs) == 0:
            continue
        sin_sum = sum(x_arr[i] * math.sin(phi_arr[i]) for i in idxs)
        cos_sum = sum(x_arr[i] * math.cos(phi_arr[i]) for i in idxs)
        mean_phi = math.atan2(sin_sum, cos_sum)

        for i in idxs:
            delta = (phi_arr[i] - mean_phi + math.pi) % two_pi - math.pi
            sum_sq_devs += float(x_arr[i]) * (delta**2)

    return float(math.sqrt(max(0.0, sum_sq_devs)))


def build_scaled_basis_t(
    x: tuple[float, ...] | np.ndarray,
    phi: tuple[float, ...] | np.ndarray,
    sheets: dict[int, list[int]],
    ref_sheet: int,
    h_eff: float,
) -> np.ndarray | None:
    m = len(x)
    eta_cols: list[np.ndarray] = []

    # 1. Columnas eta: indicadores de hojas no referencia
    for s, idxs in sorted(sheets.items(), key=lambda it: it[0]):
        if s == ref_sheet:
            continue
        col = np.zeros(m, dtype=float)
        col[idxs] = 1.0
        eta_cols.append(col)

    # 2. Columnas zeta: base ortonormal de suma cero en hojas con n_j >= 2, escaladas por h_eff
    zeta_cols: list[np.ndarray] = []
    for _s, idxs in sorted(sheets.items(), key=lambda it: it[0]):
        n_j = len(idxs)
        if n_j >= 2:
            Z = scipy.linalg.null_space(np.ones((1, n_j), dtype=float))
            for k in range(Z.shape[1]):
                col_z = Z[:, k]
                col = np.zeros(m, dtype=float)
                col[idxs] = col_z * h_eff
                zeta_cols.append(col)

    all_cols = eta_cols + zeta_cols
    if len(all_cols) != m - 1:
        return None

    T = np.column_stack(all_cols)
    aug = np.column_stack([T, np.ones(m, dtype=float)])
    if np.linalg.matrix_rank(aug) < m:
        return None

    return T


def build_g_hat_basis(
    x: np.ndarray,
    phi: np.ndarray,
    T: np.ndarray,
    h_eff: float,
) -> list[np.ndarray]:
    m = len(x)
    pairs = get_pair_indices(m)
    h_eff_sq = h_eff**2

    G_basis: list[np.ndarray] = []
    for a, b in pairs:
        d_ab = np.zeros(m, dtype=float)
        d_ab[a] = 1.0
        d_ab[b] = -1.0
        Td = T.T @ d_ab
        weight = (x[a] * x[b] * math.cos(phi[a] - phi[b])) / h_eff_sq
        Gk = weight * np.outer(Td, Td)
        Gk = 0.5 * (Gk + Gk.T)
        G_basis.append(Gk)

    return G_basis
```

---

### 4.3. Formulación y Resolución de Fase I y Fase II (`evaluators/scaled.py`)

```python
# --- FASE I: Maximización de margen y relajación elástica ---
z = cp.Variable(d)
lam_var = cp.Variable()
s_var = cp.Variable(nonneg=True)
p_z = N @ z

G_z_expr = sum(p_z[k] * G_basis[k] for k in range(K))
obj = cp.Maximize(lam_var - self.config.rho * s_var)

constraints = [
    p_z >= 0.0,
    p_z <= 1.0,
    G_z_expr >> lam_var * np.eye(m - 1),
    x_arr + A_c @ p_z >= self.config.eps_c,
    x_arr + A_d @ p_z >= self.config.mu_t - s_var,
    s_var >= 0.0,
]

prob = cp.Problem(obj, constraints)
prob.solve(**solver_kwargs)

# --- FASE II: Maximización de conectividad mínima mu* ---
z = cp.Variable(d)
t_var = cp.Variable()
p_z = N @ z

floor_margin = self.lam_floor_hat if eps is None else float(eps)
G_z_expr = sum(p_z[k] * G_basis[k] for k in range(K))

constraints = [
    x_arr + A_d @ p_z >= t_var,
    p_z >= 0.0,
    p_z <= 1.0,
    G_z_expr >> floor_margin * np.eye(m - 1),
    x_arr + A_c @ p_z >= self.config.eps_c,
]

prob = cp.Problem(cp.Maximize(t_var), constraints)
prob.solve(**solver_kwargs)
```

---

### 4.4. Función de Fitness Continua y Candados Físicos (`evaluators/fitness.py` y `evaluators/scaled.py`)

```python
# --- evaluators/fitness.py ---
def compute_continuous_fitness(
    lam: float,
    s: float,
    tol: float = 1e-6,
    lam_floor: float = 1e-7,
) -> float:
    lam_val = float(lam)
    s_val = float(s)

    if lam_val >= lam_floor:
        if s_val <= tol:
            return lam_val
        return -s_val
    else:
        # lam < lam_floor: clip a >= -1, SIEMPRE peor que cualquier estable
        lam_clipped = max(lam_val, -1.0)
        if lam_clipped > 0.0:
            lam_clipped = lam_clipped - float(lam_floor)
        return -1.0 + lam_clipped


# --- Candado de escala h_eff_min y barrera de régimen en evaluators/scaled.py ---
if h_eff < self.h_eff_min or h_eff > 0.1:
    res_plain = self.plain_evaluator.evaluate_phase1(x, phi)
    s_plain = float(res_plain.s) if res_plain.s is not None else 1.0
    lam_plain = float(res_plain.lam) if res_plain.lam is not None else float("-inf")
    is_stable_macro = (
        res_plain.status in ("OPTIMAL", "OPTIMAL_INACCURATE")
        and lam_plain >= self.config.lam_floor
        and s_plain <= self.config.s_tol
    )
    if is_stable_macro:
        f_val = lam_plain
    else:
        # Barrera de régimen (-10.0) para proteger la optimización asintótica
        f_val = min(lam_plain / 1e-8 - self.lam_floor_hat, -10.0) - s_plain

    return EvalResult(
        f=f_val,
        status="NOT_SCALABLE",
        solver_info={"h_eff": h_eff, "reason": "h_eff outside [1e-4, 0.1]"},
        ...
    )
```

---

## 5. Resumen de Constantes Críticas y Parámetros del Sistema

| Parámetro | Valor por Defecto | Justificación Teórica y Numérica |
|---|---|---|
| `mu_t` | $0.6929077$ | Conectividad objetivo ($\mu_0$ de O'Keeffe para superar $11/16 = 0.6875$). |
| `h_eff_min` | $10^{-4}$ | Límite inferior de escala; evita singularidad en la división $T^T L_w T / h_{\text{eff}}^2$. |
| `h_eff_max` | $0.1$ | Límite superior; desvíos mayores no corresponden al régimen asintótico de splitting. |
| `lam_floor_hat` | $10^{-5}$ | Piso espectral en operador escalado. Desacuerdo interno de Clarabel $\le 7.7 \cdot 10^{-8}$. |
| `lam_floor` | $10^{-7}$ | Piso espectral físico en evaluador plano (límite de resolución IEEE 754). |
| `eps_c` | $10^{-4}$ | Cota de estabilidad intra-cluster ($c_a \ge \varepsilon_c$). |
| `s_tol` | $10^{-6}$ | Tolerancia en el déficit de conectividad para declarar factibilidad de Fase I. |
| `theta_tol` | $10^{-4}$ | Tolerancia máxima permitida en residuo de equilibrio en `verify_solution`. |
| `solver_tol` | $10^{-8}$ | Tolerancia de convergencia de conos primal-dual para Clarabel y SCS. |
| `fitness_noise` | $10^{-9}$ | Umbral de desempate en torneo para evitar decisiones espurias por ruido de redondeo. |
| Barrera de Régimen | $-10.0$ | Desplaza estados `NOT_SCALABLE` por debajo de cualquier candidato en búsqueda activa. |

---

*Fin del Dossier de Arquitectura Técnica - Sistema Sgt-GENETICO.*
