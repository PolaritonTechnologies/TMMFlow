# TMM

Simulation and Optimisation software based on the Transfer Matrix Method

## Optimization Methods

Optimization benchmark (test_optimisation_DBRCavity.json on Julian's laptop
varying all layer thicknesses plus the organic's position )

| Method                 | Time to Optimize | Comments             |
| ---------------------- | ---------------- | -------------------- |
| minimize (Nelder-Mead) | 60 s             | Using 1e-9 tolerance |
| dual_annealing         |                  |                      |
| shgo                   |                  |                      |
| direct                 |                  |                      |
| differential_evolution |                  |                      |
| basinhopping           |                  |                      |
