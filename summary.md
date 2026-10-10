# FLicc CI: dev-frw @ a43935f

Run: 2026-10-10 15:54:19 UTC

| step | result |
|---|---|
| build | passed |
| test:test_dev.R | passed |
| test:test_frw.R | passed |
| compare_main | passed |

## Regression against main (alfonsino, default settings)

| case | same-params diff (obj, SPR, LFD) | fitted rel. obj diff | SPR (last yr) main / branch | conv | result |
|---|---|---|---|---|---|
| gtg_mn | 0, 0, 0 | 0.0e+00 | 0.2527 / 0.2527 | 0 / 0 | unchanged |
| gtg_dm | 0, 0, 0 | 0.0e+00 | 0.2722 / 0.2722 | 0 / 0 | unchanged |
| gtg_nb | 0, 0, 0 | 0.0e+00 | 0.2562 / 0.2562 | 1 / 1 | unchanged |
| gamma_nb | 0, 0, 0 | 0.0e+00 | 0.2478 / 0.2478 | 1 / 1 | unchanged |

Logs: build.log, tests/*.log, compare_main.log, render_*.log; rendered files in html/.
