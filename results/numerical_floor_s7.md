# Curva de Degradación del Piso Numérico en S7 (h=0.002)

Análisis empírico del comportamiento del solver Clarabel con tolerancias 1e-10 al variar el margen espectral $\varepsilon$.

| $\varepsilon$ | $\Phi$ | Estatus del Solver | $\delta$ | Tiempo (s) |
|---|---|---|---|---|
| 1.0e-03 | 0.999730 | `OPTIMAL` | 1.00e-03 | 0.0211 |
| 1.0e-04 | 1.000000 | `SOLVER_ERROR` | -1.00e+00 | 0.0187 |
| 1.0e-05 | 0.972971 | `OPTIMAL` | 9.73e-06 | 0.0764 |
| 1.0e-06 | 0.644094 | `OPTIMAL_INACCURATE` | 6.44e-07 | 0.0237 |
| 1.0e-07 | 0.000000 | `OPTIMAL_INACCURATE` | 0.00e+00 | 0.0276 |
| 1.0e-08 | 0.000000 | `OPTIMAL_INACCURATE` | 0.00e+00 | 0.0242 |
| 1.0e-09 | 1.000000 | `UNBOUNDED` | -1.00e+00 | 0.0164 |
| 1.0e-10 | 1.000000 | `UNBOUNDED` | -1.00e+00 | 0.0121 |


**Conclusión**: Para $\varepsilon < 10^{-7}$, el problema semidefinido entra en la zona de ruido numérico de las tolerancias $10^{-8}\dots 10^{-10}$, provocando estados `OPTIMAL_INACCURATE` o fallos de convergencia.
