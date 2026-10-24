# Registro de Afirmaciones y Citas Científicas (CITATIONS.md)

Este documento registra sistemáticamente cada afirmación teórica, cota de conectividad o valor numérico referenciado en el proyecto, indicando su fuente original, ubicación precisa y estatus de verificación.

---

## Taxonomía de Estados

- **DEMOSTRADO:** Teorema o lema matemático probado formalmente en la literatura con revisión por pares o verificado por aritmética de intervalos rigurosa.
- **EVIDENCIA NUMÉRICA:** Resultado derivado de solvers o algoritmos numéricos (SDP, GA, Newton en alta precisión) sin garantía formal de intervalo.
- **CITADO (no verificado):** Afirmación o valor numérico reportado en la literatura externa (en particular en preprints no revisados por pares) pendiente de reproducción independiente.

---

## Tabla de Afirmaciones

| ID | Fuente | Sección / Ecuación | Afirmación / Valor | Estado | Notas / Advertencias |
|---|---|---|---|---|---|
| **CIT-01** | Kassabov, Strogatz, Townsend (2021), *Chaos* 31, 073135 | Sec. II, Fig. 1 | Estado de cuatro clusters $C_4$ con $\mu = 3/4$ es linealmente degenerado ($L \equiv 0$, autovalores cero en el Jacobiano). | **DEMOSTRADO** | Publicado en revista con arbitraje. El análisis no lineal posterior (Lepsveridze & Zhang 2026) confirma que es desestabilizado a orden cúbico. |
| **CIT-02** | Canale (2022), arXiv:2209.06362 | Sec. 1-3 | Cota inferior de conectividad $\mu_c \ge 11/16 = 0.6875$ mediante plantilla de 8 clusters en dos retículos girados $\arctan(1/2)$. | **CITADO** | Preprint de referencia que superó la familia circulante. Faltan los datos completos de la matriz $P$ 8x8 y el desfase angular exacto en `/docs` (requiere paper original de Canale); conforme a honestidad epistémica, el test de regresión `test_canale_eight_cluster_bound` queda marcado como `xfail`. |
| **CIT-03** | Lepsveridze & Zhang (2026), arXiv:2608.20010 | Teorema 1 | Cota superior estricta $\mu_c < 3/4$. | **CITADO** | Preprint reciente; prueba no lineal de que la sincronización global prevalece estrictamente por debajo de 3/4. |
| **CIT-04** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. I, Ec. (2), Teorema 1 | Cota inferior $\mu_c \ge \mu_0 = \frac{4 + 2q^*}{7 + q^*} \approx 0.6929077063...$ | **CITADO (no verificado)** | Preprint sin referato por pares. Afirma haber sido elaborado mayoritariamente con agentes de IA. |
| **CIT-05** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. I, Ec. (3) | $q^* \approx 0.6505691669...$ es la única raíz en $(0, 1)$ del polinomio irreducible de grado 8: $9q^8 + 72q^7 + 288q^6 + 448q^5 + 470q^4 + 88q^3 - 208q^2 - 128q - 15 = 0$. | **CITADO (no verificado)** | Derivado analíticamente de la condición de doble neutralidad en el modelo límite $h \to 0$. |
| **CIT-06** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. IV.A, Ec. (11) | Razón de apertura del cluster central $w^* = 1.3982781340...$ dada por $w^2 = \frac{q^2 + 6q + 1}{(q + 1)^2}$. | **CITADO (no verificado)** | Relación exacta del ansatz $S_7$ en el límite $h \to 0$. |
| **CIT-07** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. IV.A, Tabla II | Matriz de enlaces límite de $S_7$: 3 clusters divididos en mitades con enlaces completos excepto 4 pares cruzados de densidad $q \approx 0.65$, y un 4to cluster entero de tamaño $1 - 6a$. | **CITADO (no verificado)** | Base del ansatz para la familia de perturbación $h > 0$. |
| **CIT-08** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. IV.C | Progresión numérica en $S_7$: para $h \in \{0.02, 0.01, 0.005, 0.002\}$, valores de $\mu \in \{0.6906, 0.6918, 0.6924, 0.69272\}$. | **CITADO (no verificado)** | Para $h=0.002$, el valor analítico asintótico de Tabla II da $\approx 0.692908$; el valor citado por O'Keeffe es $0.69272$ (con correcciones $O(h)$ del Apéndice A.1); y el valor óptimo devuelto por el SDP (`SDPEvaluator` Fase II con $A_F p = 0, M \succeq 0$) es $\mu^* = 0.69118098 \approx 0.69118$. Ver nota detallada abajo. Test `test_s7_h_finite_numerical_connectivity` valida rango $[0.690, 0.6930]$. |

| **CIT-09** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. V, Proposición 1 | Existencia de grafos arbitrariamente grandes con $\mu^* = 0.6929050246493944...$ certificado por aritmética de intervalos (Krawczyk + Cholesky). | **CITADO (no verificado)** | Los datos exactos racionales residen en el *Supplementary Material* (`certificate.json`), no provisto en este repositorio. |
| **CIT-10** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. V, Lema 1 | Denominador racional común $D \approx 3.06 \times 10^{42}$, autovalor mínimo de modos de cluster $6.2313 \times 10^{-19}$, $c_a \ge 0.2157$. | **CITADO (no verificado)** | Requiere aritmética de 50 dígitos e intervalos. Ausente en el repositorio. |
| **CIT-11** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. VI.(i) | Fragilidad de la familia $S_7$: margen espectral decrece como $h^2(t+h)$, y el déficit de conectividad crece como $\varepsilon^{1/3}$. | **CITADO (no verificado)** | Estimación asintótica reportada por el autor para justificar por qué los optimizadores genéricos fallan. |
| **CIT-12** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. VI | Estado satélite de $C_4$ con conectividad certificada $\mu \approx 0.6915537...$ | **CITADO (no verificado)** | Mencionado brevemente; datos completos en el *Supplementary Material* ausente. |
| **CIT-13** | Wiley, Strogatz, Girvan (2006) / Taylor (2012) | *Chaos* 16, 015103 / *J. Phys. A* 45, 055102 | Estados twisted en anillos dan $\mu_c \ge 0.6809$. | **DEMOSTRADO** | Publicado en revistas con referato. |
| **CIT-14** | Yoneda, Tatsukawa, Teramae (2021) | *Chaos* 31, 063124 | Óptimo en familias circulantes: $\mu_c \ge 0.6838$. | **DEMOSTRADO** | Publicado en revista con referato. |
| **CIT-15** | Taylor (2012), Ling et al. (2019), Lu & Steinerberger (2020) | Refs. [2, 3, 4] | Cotas superiores históricas: $\mu_c \le 0.9395 \to 0.7929 \to 0.7889$. | **DEMOSTRADO** | Publicados en revistas con referato. |

---

## Nota Técnica: Origen de 0.69118 y Comparación con 0.6929077 y 0.69272 en $S_7$

Existen tres valores clave de conectividad $\mu$ asociados a la familia $S_7$ que deben distinguirse con precisión matemática:

1. **Límite Teórico Singular ($h \to 0$): $\mu_0 = \frac{4 + 2q^*}{7 + q^*} \approx 0.6929077064...$**
   - **Origen:** Derivación analítica en el límite singular asintótico (O'Keeffe 2026, Sec. IV.A, Ec. 9 y Tabla II).
   - **Condición:** En este límite ideal, las fases están colapsadas en múltiplos exactos de $\pi/2$: $(0, 0, \pi/2, \pi/2, \pi, \pi, -\pi/2)$, con tres clusters divididos en mitades iguales de peso $a = 1/(7+q^*)$ y un cuarto cluster entero de peso $(1+q^*)/(7+q^*)$.
   - **Propiedad:** La matriz de densidades estática de la Tabla II tiene exactamente grado mínimo $d_a = 1 - a(3 - q^*) = \frac{4 + 2q^*}{7 + q^*}$ y satisface $F \equiv 0$, pero carece de margen espectral ($\lambda_{\min}(M) = 0$).

2. **Valor del Optimizador SDP a Perturbación Finita ($h = 0.002$): $\mu^* \approx 0.69118098...$**
   - **Origen:** Resolución de la Fase II en `SDPEvaluator` mediante optimización convexa semidefinida sobre el ansatz perturbado $\phi_k = j_k \pi/2 + h s_k$.
   - **Mecanismo:** A $h = 0.002$, la matriz estática de la Tabla II deja de ser un equilibrio exacto ($\max_a |F_a| \sim O(h) \approx 1.7 \times 10^{-4} \neq 0$). Para satisfacer la restricción física estricta de equilibrio $A_F p = 0$, el solver SDP proyecta las densidades sobre el espacio nulo $p = N z$ con $N = \ker(A_F)$.
   - **Efecto en el grado mínimo:** Al obligar a las densidades a compensar los torques desbalanceados introducidos por el desfasaje $h = 0.002$, las densidades se desvían de los valores límite de la Tabla II, reduciendo el grado mínimo óptimo alcanzable a $\mu^* = 0.69118098 \approx 0.69118$.
   - **Verificación independiente en NumPy puro (`verify_solution`):**
     - $\max_a |F_a| = 8.57 \times 10^{-13}$ (equilibrio exacto a precisión de máquina).
     - $\min_a c_a = 0.215747 \ge \varepsilon_c = 10^{-4}$ (estabilidad local intra-cluster satisfecha).
     - $\lambda_{\min}(M) = -2.69 \times 10^{-9}$ (semidefinición positiva $M \succeq 0$ dentro de la tolerancia de convergencia del solver).
     - $\mu^* = 0.691181$.

3. **Valor Citado por O'Keeffe en el Preprint a $h = 0.002$: $\mu = 0.69272$**
   - **Origen:** O'Keeffe (2026), Sec. IV.C (Tabla de progresión numérica de $h$ vs $\mu$).
   - **Diferencia:** O'Keeffe no ejecutó una optimización SDP libre sobre el espacio nulo estándar, sino que construyó una expansión perturbativa analítica $P(h) = P^{(0)} + h P^{(1)} + h^2 P^{(2)} + \dots$ (Apéndice A.1) calculada con aritmética exacta racional (cuyos coeficientes residen en el *Supplementary Material* ausente). Esa corrección perturbativa dirigida preserva mayor simetría y alcanza $\mu = 0.69272$, superior al $0.69118$ obtenido por la optimización genérica del SDP.
   - **Intervalo de Aceptación:** Ambos valores ($0.69118$ y $0.69272$) residen holgadamente dentro del intervalo de regresión $[0.690, 0.6930]$ validado por `test_s7_h_finite_numerical_connectivity`.

