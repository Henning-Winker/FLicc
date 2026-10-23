# FLicc CI: dev-docs @ 2997f35

Run: 2026-10-23 11:08:01 UTC

| step | result |
|---|---|
| build | passed |
| test:test_dev.R | passed |
| test:test_binormal.R | passed |
| test:test_catch_weight.R | passed |
| test:test_combine_lfd.R | passed |
| test:test_dm_theta.R | passed |
| test:test_domesens.R | passed |
| test:test_ensemble.R | passed |
| test:test_frw.R | passed |
| test:test_growth.R | passed |
| test:test_lbb_data.R | passed |
| test:test_lbb_example.R | passed |
| test:test_lbb_illex.R | passed |
| test:test_lbb.R | passed |
| test:test_lbi.R | passed |
| test:test_make_lbspr_data.R | passed |
| test:test_maxsd_default.R | passed |
| test:test_mc.R | passed |
| test:test_newton_dm.R | passed |
| test:test_quickstart.R | passed |
| test:test_resid.R | passed |
| test:test_selF.R | passed |
| test:test_spr_gtg.R | passed |
| compare_main | passed |
| render:inst/guide/preview.Rmd | passed |

## Regression against main (alfonsino, default settings)

| case | same-params diff (obj, SPR, LFD) | fitted rel. obj diff | SPR (last yr) main / branch | conv | result |
|---|---|---|---|---|---|
| gtg_mn | 0, 0, 0 | 0.0e+00 | 0.2327 / 0.2327 | 0 / 0 | unchanged |
| gtg_dm | 0, 0, 0 | 0.0e+00 | 0.2264 / 0.2264 | 0 / 0 | unchanged |
| gtg_nb | 0, 0, 0 | 0.0e+00 | 0.2195 / 0.2195 | 0 / 0 | unchanged |
| gamma_nb | 0, 0, 0 | 0.0e+00 | 0.2187 / 0.2187 | 0 / 0 | unchanged |

Logs: build.log, tests/*.log, compare_main.log, render_*.log; rendered files in html/.
