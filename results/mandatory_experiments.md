# Evidencia Numérica: Experimentos Obligatorios por $\mu_t$ y por $m$

> [!IMPORTANT]
> **Evidencia numérica**: Ningún resultado aquí tabulado se declara cota inferior.
> Todas las métricas representan valores observados por optimización numérica con tolerancia $10^{-8}$.

## Resumen Global por $\mu_t$

| $\mu_t$ | Mejor $f$ | Mejor $m$ | $\lambda^*$ crudo | $s^*$ | Evals | Cache Hits |
|---|---|---|---|---|---|---|
| 0.69 | 4.589500e-09 | 7 | 4.589500e-09 | 2.725566e-11 | 81 | 29 |
| 0.6929077 | -1.193758e-10 | 4 | -1.193758e-10 | 3.207884e-12 | 67 | 43 |
| 0.6935 | -1.269302e-10 | 4 | -1.269302e-10 | 4.046245e-12 | 67 | 43 |
| 0.695 | -1.475711e-10 | 4 | -1.475711e-10 | 6.389292e-12 | 67 | 43 |

## Mejor Métrica por cada $m \in [2, 10]$ (rotulado evidencia numérica)

| $\mu_t$ | m=2 | m=3 | m=4 | m=5 | m=6 | m=7 | m=8 | m=9 | m=10 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.69 | -s*=-3.42e-01 | -s*=-2.72e-01 | lam*=-8.68e-11 | -s*=-1.02e-01 | lam*=-1.53e-07 | lam*=4.59e-09 | lam*=-1.38e-06 | - | - |
| 0.6929077 | -s*=-3.45e-01 | -s*=-2.75e-01 | lam*=-1.19e-10 | lam*=-1.67e-06 | lam*=-1.02e-04 | lam*=-6.81e-05 | -s*=-2.14e-01 | - | - |
| 0.6935 | -s*=-3.46e-01 | -s*=-2.76e-01 | lam*=-1.27e-10 | lam*=-1.68e-06 | lam*=-2.32e-04 | lam*=-2.15e-04 | -s*=-2.14e-01 | - | - |
| 0.695 | -s*=-3.47e-01 | -s*=-2.77e-01 | lam*=-1.48e-10 | lam*=-1.69e-06 | lam*=-8.20e-04 | lam*=-7.94e-04 | -s*=-2.16e-01 | - | - |
