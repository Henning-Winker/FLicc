# FLicc CI: dev-gtgdyn @ 84d1311

Run: 2026-10-12 11:53:10 UTC

| step | result |
|---|---|
| build | passed |
| test:test_dev.R | passed |
| test:test_frw.R | passed |
| test:test_gtgdyn.R | passed |
| compare_main | passed |
| render:docs/FLicc_dev_tests.Rmd | passed |
| render:docs/gtgdyn_alfonsino_OM_test.Rmd | passed |

## Regression against main (alfonsino, default settings)

| case | same-params diff (obj, SPR, LFD) | fitted rel. obj diff | SPR (last yr) main / branch | conv | result |
|---|---|---|---|---|---|
| gtg_mn | 2.3e-15, 0, 0 | 7.7e-13 | 0.2527 / 0.2527 | 0 / 0 | unchanged |
| gtg_dm | 5.9e-15, 0, 0 | 3.0e-12 | 0.2722 / 0.2722 | 0 / 0 | unchanged |
| gtg_nb | 9.8e-13, 0, 0 | 1.0e-12 | 0.2562 / 0.2562 | 1 / 0 | unchanged |
| gamma_nb | 1e-12, 7.9e-13, 4.7e-14 | 5.3e-12 | 0.2478 / 0.2478 | 1 / 0 | unchanged |

Logs: build.log, tests/*.log, compare_main.log, render_*.log; rendered files in html/.
