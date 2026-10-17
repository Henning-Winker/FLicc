# FLicc CI: dev-merge @ 2304b2f

Run: 2026-10-17 11:20:02 UTC

| step | result |
|---|---|
| build | passed |
| test:test_dev.R | passed |
| test:test_catch_weight.R | passed |
| test:test_domesens.R | passed |
| test:test_frw.R | passed |
| test:test_quickstart.R | passed |
| test:test_resid.R | passed |
| test:test_selF.R | passed |
| compare_main | changed (see table) |
| render:docs/FLicc_dev_tests.Rmd | passed |
| render:inst/guide/FLicc_user_guide.Rmd | passed |

## Regression against main (alfonsino, default settings)

| case | same-params diff (obj, SPR, LFD) | fitted rel. obj diff | SPR (last yr) main / branch | conv | result |
|---|---|---|---|---|---|
| gtg_mn | 0.11, 0.033, 0.0066 | 2.0e-04 | 0.2527 / 0.1015 | 0 / 0 | CHANGED |
| gtg_dm | 0.1, 0.033, 0.0066 | 2.5e-04 | 0.2722 / 0.1142 | 0 / 0 | CHANGED |
| gtg_nb | 1.2, 0.033, 0.0066 | 4.3e-03 | 0.2562 / 0.1050 | 1 / 1 | CHANGED |
| gamma_nb | 1.3, 0.033, 0.0066 | 3.6e-03 | 0.2478 / 0.0999 | 1 / 0 | CHANGED |

Logs: build.log, tests/*.log, compare_main.log, render_*.log; rendered files in html/.
