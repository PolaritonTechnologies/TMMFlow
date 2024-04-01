# TMM

Simulation and Optimisation software based on the Transfer Matrix Method

## Optimization Methods

Optimization benchmark (test_optimisation_DBRCavity.json on Julian's laptop
varying all layer thicknesses plus the organic's position )

| Method                 | Time to Optimize | Comments                                                               |
| ---------------------- | ---------------- | ---------------------------------------------------------------------- |
| minimize (Nelder-Mead) | 60 s             | Using 1e-9 tolerance                                                   |
| dual_annealing         | 1802 s           | Using 1e-5 tolerance                                                   |
| basinhopping           | 126 s            | Does not support bounds (recovers results of minimize for the example) |
| direct                 |                  | Did not yield sensible results after 2000 iterations                   |
| differential_evolution |                  |                                                                        |
| shgo                   |                  |                                                                        |

dual_annealing, optimized weights:
[ 43.44974405 32.76934884 31.45234875 190.87887203 65.24474902
140.11114023 19.47791629 132.36720361 324.52413575 81.47214849
69.120849 239.70597544 62.72602507 282.03428307 447.86895182
459.8029801 143.59800317 484.79319159 178.5618158 24.86932182
221.55959766 15.84151904 2.68046251]

![Dual annealing optimized image](/doc/dual-annealing-optimized.png)
