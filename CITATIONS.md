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
| **CIT-08** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. IV.C | Progresión numérica en $S_7$: para $h \in \{0.02, 0.01, 0.005, 0.002\}$, valores de $\mu \in \{0.6906, 0.6918, 0.6924, 0.69272\}$. | **CITADO (no verificado)** | Para $h=0.002$, el SDP da $\mu^* = 0.69118$ mientras que el preprint cita $0.69272$. El SDP da el óptimo para (x,phi) fijos del ansatz; el miembro citado probablemente usa otros (x,phi); no verificable sin el Supplementary Material. Test `test_s7_h_finite_numerical_connectivity` valida rango $[0.690, 0.6930]$. |

| **CIT-09** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. V, Proposición 1 | Existencia de grafos arbitrariamente grandes con $\mu^* = 0.6929050246493944...$ certificado por aritmética de intervalos (Krawczyk + Cholesky). | **CITADO (no verificado)** | Los datos exactos racionales residen en el *Supplementary Material* (`certificate.json`), no provisto en este repositorio. |
| **CIT-10** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. V, Lema 1 | Denominador racional común $D \approx 3.06 \times 10^{42}$, autovalor mínimo de modos de cluster $6.2313 \times 10^{-19}$, $c_a \ge 0.2157$. | **CITADO (no verificado)** | Requiere aritmética de 50 dígitos e intervalos. Ausente en el repositorio. |
| **CIT-11** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. VI.(i) | Fragilidad de la familia $S_7$: margen espectral decrece como $h^2(t+h)$, y el déficit de conectividad crece como $\varepsilon^{1/3}$. | **CITADO (no verificado)** | Estimación asintótica reportada por el autor para justificar por qué los optimizadores genéricos fallan. |
| **CIT-12** | O'Keeffe (2026), arXiv:2610.21271v1 | Sec. VI | Estado satélite de $C_4$ con conectividad certificada $\mu \approx 0.6915537...$ | **CITADO (no verificado)** | Mencionado brevemente; datos completos en el *Supplementary Material* ausente. |
| **CIT-13** | Wiley, Strogatz, Girvan (2006) / Taylor (2012) | *Chaos* 16, 015103 / *J. Phys. A* 45, 055102 | Estados twisted en anillos dan $\mu_c \ge 0.6809$. | **DEMOSTRADO** | Publicado en revistas con referato. |
| **CIT-14** | Yoneda, Tatsukawa, Teramae (2021) | *Chaos* 31, 063124 | Óptimo en familias circulantes: $\mu_c \ge 0.6838$. | **DEMOSTRADO** | Publicado en revista con referato. |
| **CIT-15** | Taylor (2012), Ling et al. (2019), Lu & Steinerberger (2020) | Refs. [2, 3, 4] | Cotas superiores históricas: $\mu_c \le 0.9395 \to 0.7929 \to 0.7889$. | **DEMOSTRADO** | Publicados en revistas con referato. |

---

---

## Nota Técnica: Comparación de Conectividad en $S_7$ ($0.69118$ vs $0.69272$)

El SDP da el óptimo para (x,phi) fijos del ansatz; el miembro citado probablemente usa otros (x,phi); no verificable sin el Supplementary Material.

## Nota Técnica: Estabilidad Espectral y Slack de Grado

En la formulación de Fase I ($\max \lambda - \rho s$ s.a. $d \ge \mu_t - s$, $M \succeq \lambda I$, $s \ge 0$):
$s^* \approx 0$ **NO implica estabilidad**; solo $\lambda^* \ge \text{piso}$ la implica. El slack $s^* \approx 0$ puede coexistir con $\lambda^* < 0$ cuando el solver compensa grado a expensas de estabilidad espectral.

## Nota Técnica: Conectividad Apéndice A en $S_7$ ($h=0.02$)

Para $h = 0.02, t = 0$, la matriz $P$ analítica del Apéndice A da una conectividad $\mu(\text{Apéndice A}, h=0.02) = 0.6735$, frente a $0.6906$ citado en el preprint.



