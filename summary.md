# gtg.dyn with catch in the likelihood

gtg.dyn (19 GTG, fast numerics), penalised R deviations (sigmaR 0.5) unless RE; catch weight 1995-2024, LFD 2005-2024.
Relative error: median (RMSE) over 4 datasets per scenario. 'final' = last 3 years 2022-2024.

| scenario | model | median s | SPR rebuild | SPR recent | SPR final 3 yr | F final 3 yr |
|---|---|---|---|---|---|---|
| noRdev | eq | 0.2 | -22% (22%) | -3% (8%) | -5% (8%) | +29% (31%) |
| noRdev | no_catch | 1.5 | -7% (10%) | -4% (14%) | -3% (16%) | +4% (40%) |
| noRdev | c_exact | 4.7 | +1% (5%) | -1% (9%) | -3% (11%) | +7% (24%) |
| noRdev | c_cv10 | 3.3 | -1% (6%) | -1% (11%) | -5% (11%) | +12% (25%) |
| noRdev | c_under | 3.4 | -1% (6%) | -1% (11%) | -5% (11%) | +12% (25%) |
| noRdev | c_cv10_RE | 30.4 | +5% (8%) | +2% (7%) | -2% (4%) | +5% (9%) |
| Rdev | eq | 0.2 | -20% (23%) | +9% (16%) | +7% (13%) | +5% (18%) |
| Rdev | no_catch | 1.4 | +6% (14%) | -1% (11%) | -9% (13%) | +20% (33%) |
| Rdev | c_exact | 4.2 | +5% (9%) | +6% (10%) | +6% (10%) | -8% (14%) |
| Rdev | c_cv10 | 3.2 | +4% (12%) | +1% (10%) | +2% (10%) | -0% (15%) |
| Rdev | c_under | 3.2 | +4% (12%) | +1% (10%) | +2% (10%) | -0% (15%) |
| Rdev | c_cv10_RE | 24.4 | +5% (12%) | +1% (11%) | +2% (11%) | -1% (15%) |
