# gtg.dyn TMB trial

Stage: penalised

Reference test (TMB vs Python, same inputs): dyn exact 2.7e-12, dyn grid 3.1e-12, eq 2.5e-12 -> **PASS**

## Median fit time (s)

| model | TMB | Python |
|---|---|---|
| dyn_dome | 1.6 | 216 |
| dyn_dome_Rdev | 6.4 | 620 |
| dyn_dome_Rrw | 13.8 | 695 |
| eq_dome | 0.3 | 5 |

## Relative error in SPR: median (RMSE)

| source | scenario | model | decline | rebuild | recent | final |
|---|---|---|---|---|---|---|
| Python | noRdev | dyn_dome | +1% (6%) | -2% (8%) | -6% (11%) | -17% (18%) |
| TMB | noRdev | dyn_dome | -0% (6%) | -3% (8%) | -7% (11%) | -17% (18%) |
| Python | noRdev | dyn_dome_Rdev | +2% (7%) | -6% (10%) | -2% (15%) | -6% (14%) |
| TMB | noRdev | dyn_dome_Rdev | +0% (7%) | -6% (9%) | -5% (13%) | -10% (15%) |
| Python | noRdev | dyn_dome_Rrw | +1% (7%) | -7% (10%) | -4% (15%) | -9% (16%) |
| TMB | noRdev | dyn_dome_Rrw | -0% (8%) | -8% (10%) | -5% (15%) | -8% (17%) |
| Python | noRdev | dyn_logistic | +1% (4%) | -7% (10%) | -2% (12%) | -27% (27%) |
| Python | noRdev | eq_dome | -9% (12%) | -22% (22%) | -3% (8%) | -8% (8%) |
| TMB | noRdev | eq_dome | -8% (11%) | -19% (22%) | +1% (28%) | -8% (36%) |
| Python | noRdev | eq_logistic | -9% (8%) | -29% (27%) | -9% (9%) | -15% (15%) |
| Python | Rdev | dyn_dome | +3% (20%) | +13% (23%) | -2% (17%) | -17% (24%) |
| TMB | Rdev | dyn_dome | +5% (20%) | +11% (22%) | -7% (16%) | -20% (24%) |
| Python | Rdev | dyn_dome_Rdev | +0% (9%) | +6% (14%) | -0% (11%) | -13% (15%) |
| TMB | Rdev | dyn_dome_Rdev | +1% (10%) | +2% (11%) | -3% (11%) | -19% (18%) |
| Python | Rdev | dyn_dome_Rrw | +1% (9%) | +5% (15%) | +3% (14%) | -6% (13%) |
| TMB | Rdev | dyn_dome_Rrw | +1% (10%) | +2% (10%) | +1% (10%) | -15% (15%) |
| Python | Rdev | dyn_logistic | -7% (13%) | +20% (24%) | -2% (16%) | -32% (32%) |
| Python | Rdev | eq_dome | -8% (11%) | -20% (23%) | +9% (16%) | -1% (6%) |
| TMB | Rdev | eq_dome | -7% (11%) | -18% (24%) | +20% (29%) | -1% (28%) |
| Python | Rdev | eq_logistic | -14% (15%) | -27% (26%) | +2% (4%) | -13% (13%) |

## Max |SPR TMB - SPR Python| for the same penalised model

| model | scenario | max diff |
|---|---|---|
| dyn_dome | noRdev | 0.013 |
| dyn_dome_Rdev | noRdev | 0.033 |
| dyn_dome_Rrw | noRdev | 0.011 |
| eq_dome | noRdev | 0.312 |
| dyn_dome | Rdev | 0.045 |
| dyn_dome_Rdev | Rdev | 0.065 |
| dyn_dome_Rrw | Rdev | 0.114 |
| eq_dome | Rdev | 0.208 |

Convergence: 25 of 32 fits
