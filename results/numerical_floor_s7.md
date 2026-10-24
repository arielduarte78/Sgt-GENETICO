# Curva de Degradación del Piso Numérico en S7 (h=0.002)

Análisis empírico del comportamiento del solver Clarabel con tolerancias 1e-10 al variar el margen espectral $\varepsilon$.

| $\varepsilon$ | $\mu^*$ | Estatus del Solver | Tiempo (s) |
|---|---|---|---|
| 1.0e-03 | None | `INFEASIBLE` | 0.0194 |
| 1.0e-04 | None | `INFEASIBLE` | 0.0144 |
| 1.0e-05 | None | `INFEASIBLE` | 0.0203 |
| 1.0e-06 | None | `INFEASIBLE` | 0.0155 |
| 1.0e-07 | None | `VERIFICATION_FAILED` | 0.0274 |
| 1.0e-08 | 0.688776 | `OPTIMAL_INACCURATE` | 0.0156 |
| 1.0e-09 | 0.691147 | `OPTIMAL_INACCURATE` | 0.0173 |
| 1.0e-10 | None | `VERIFICATION_FAILED` | 0.0158 |


**Conclusión**: Para $\varepsilon \ge 10^{-6}$, el requerimiento espectral es infactible. Para $\varepsilon < 10^{-7}$, el problema semidefinido entra en la zona de ruido numérico de las tolerancias $10^{-8}\dots 10^{-10}$, provocando estados `OPTIMAL_INACCURATE` o fallos de verificación.
