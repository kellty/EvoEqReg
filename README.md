# Implementation of Spatio-temporal Regression with Evolution Equation Constraint

## Data
Fluorescence recovery after photobleaching (FRAP) data from [https://zenodo.org/records/3874218](https://zenodo.org/records/3874218) was used for this work. An *.RData* file in "***real_prep.zip***" was obtained after pre-processing by "***real_prep.R***".

## Code
- "***simu_begin.R***": Generate the simulation setting.
- "***simu_rand.R***": Perform regression on simulated data under different hyperparameters.
- "***real_cut.R***" ＆ "***real_pen.R***": Perform regression on FRAP data using different methods, with tuning parameters selected to be optimal.
- "***real_sum.R***": Collect regression results for FRAP data and draw pictures.

## Workflow
- Simulation: Run "***simu_begin.R***" and *then* "***simu_rand.R***".
- Real data example: Make sure that "***real_prep.RData***" is in the working directory. Run "***real_cut.R***" ＆ "***real_pen.R***" and *then* "***real_sum.R***".
