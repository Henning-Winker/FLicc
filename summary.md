# gtg.dyn with catch in the likelihood

gtg.dyn (19 GTG, fast numerics), penalised R deviations (sigmaR 0.5) unless RE; catch weight 1995-2024, LFD 2005-2024.
Relative error: median (RMSE) over 4 datasets per scenario. 'final' = last 3 years 2022-2024.

| scenario | model | median s | SPR rebuild | SPR recent | SPR final 3 yr | F final 3 yr |
|---|---|---|---|---|---|---|
| noRdev | eq | 0.7 | -22% (22%) | -3% (8%) | -5% (8%) | +29% (31%) |
| noRdev | no_catch | 3.8 | -7% (10%) | -4% (14%) | -3% (16%) | +4% (40%) |
| noRdev | c_exact | 17.3 | +1% (5%) | -1% (9%) | -3% (11%) | +7% (24%) |
| noRdev | c_cv10 | 13.8 | -3% (6%) | -6% (10%) | -7% (11%) | +15% (26%) |
| noRdev | c_under | 15.9 | -1% (6%) | -1% (11%) | -5% (11%) | +12% (25%) |
| noRdev | c_cv10_RE | 168.1 | +2% (6%) | -1% (7%) | -3% (5%) | +7% (11%) |
| Rdev | eq | 1.0 | -20% (23%) | +9% (16%) | +7% (13%) | +5% (18%) |
| Rdev | no_catch | 3.0 | +6% (14%) | -1% (11%) | -9% (13%) | +20% (33%) |
| Rdev | c_exact | 11.2 | +5% (9%) | +6% (10%) | +6% (10%) | -8% (14%) |
| Rdev | c_cv10 | 7.6 | +4% (12%) | +1% (10%) | +2% (10%) | -0% (15%) |
| Rdev | c_under | 7.0 | +4% (12%) | +1% (10%) | +2% (10%) | -0% (15%) |
| Rdev | c_cv10_RE | 76.7 | +5% (12%) | +1% (11%) | +2% (11%) | -1% (15%) |
