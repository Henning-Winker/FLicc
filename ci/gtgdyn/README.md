# gtg.dyn TMB trial

Standalone trial of the dynamic length-indexed GTG model (`gtgdyn.cpp`); FLicc itself is unchanged.

- `gtgdyn.cpp`: TMB template (pop 0 = equilibrium GTG, 1 = gtg.dyn).
- `model_q.py`: Python reference implementation (quadrature form) used to produce `data/ref_*.csv`.
- `run_trial.R`: reference test (TMB vs Python, must agree to 1e-8), then fits to the 8 simulated
  rebuild datasets in `data/` (penalised as in the Python trial, plus random-effects variants), with timings.
- Workflow `.github/workflows/gtgdyn-trial.yml` publishes `summary.md`, logs and csv files to the
  branch `ci-out/gtgdyn-trial`.

Design and Python results: project doc `claude/gtgdyn_design.md`.
