# FLicc CI: dev-gtgdyn @ 84d1311

Run: 2026-11-24 08:19:01 UTC

| step | result |
|---|---|
| build | passed |
| test:test_dev.R | passed |
| test:test_frw.R | passed |
| test:test_gtgdyn.R | passed |
| compare_main | changed (see table) |
| render:docs/gtgdyn_alfonsino_OM_test.Rmd | passed |

## Regression against main (alfonsino, default settings)

| case | same-params diff (obj, SPR, LFD) | fitted rel. obj diff | SPR (last yr) main / branch | conv | result |
|---|---|---|---|---|---|
| gtg_mn | 0.11, 0.013, NA | 3.3e-02 | 0.2327 / 0.2527 | 0 / 0 | CHANGED |
| gtg_dm | 0.12, 0.013, NA | 3.4e-02 | 0.2264 / 0.2722 | 0 / 0 | CHANGED |
| gtg_nb | 0.49, 0.013, NA | 4.0e-02 | 0.2195 / 0.2562 | 0 / 0 | CHANGED |
| gamma_nb | 0.49, 0.013, NA | 3.9e-02 | 0.2187 / 0.2478 | 0 / 0 | CHANGED |

Logs: build.log, tests/*.log, compare_main.log, render_*.log; rendered files in html/.
