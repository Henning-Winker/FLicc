# FLicc-SPM prototype (branch dev-spm, from dev-sex)

One TMB objective: FLicc length likelihood (gtg, one or two sexes; F as fixed
effects, a nuisance that tunes SPR; `F_re = TRUE` still possible) + a Pella-Tomlinson state-space surplus production model in spawning
biomass fitted to catch and CPUE + coupling. Switched on by `settings$spm`
(`spm_settings_flicc()`); without it FLicc is unchanged (`spm_use = 0`, SPM
parameters mapped off).

## What is in the tape
- **Per-recruit block** (AD-safe copy of the brp block): SBPR (female), YPR and
  EBPR (both sexes) at annual apical F with the per-recruit selectivity.
- **Reference F**: Fmsy (`shape = "fmsy"`, maximum of R(F)·YPR(F) under
  Beverton-Holt steepness) or Fspr (`shape = "spr"`), by taped damped Newton
  iterations on log F from the best node of the F grid. TMB tapes once, so a
  root solved on doubles would be frozen at the taping values; the iterations
  are therefore on the tape (`newton.it`, default 12; 6 is enough from the grid).
- **Shape and productivity**: BmsyK = SB(Fref)/SB0, m from BmsyK = m^(-1/(m-1))
  (taped Newton), Hmsy = YPR/EBPR at Fref, r = m·Hmsy (JABBA3 convention).
  r enters as a penalty log r ~ N(log r(theta), r.sd) (default 0.2); m fixed at
  m(theta) unless `m.sd > 0`.
- **EB/SB curves**: equilibrium SB/SB0 and EB/SB on the F grid (`xgrid`,
  multiples of M at Linf), Nadaraya-Watson kernel smoother in log SB/SB0
  (Gaussian, bandwidth 0.25). Must be smooth: SB/SB0 is a random effect, and
  with linear interpolation the Laplace log-determinant jumps when a state
  crosses a node (negative Hessian eigenvalues, false convergence). Column 0
  for removals (fishery), one per index (gear selectivity, joint fishery, or SB).
- **SPM**: translation of the subset of `JABBA3::jb3_nll()` (Selex state):
  constant process error, catch without error, production lag, soft bounds on
  P, JABBA priors on bkinit, sdb, sdi, and JABBA3's one-sided K prior
  ("half8": flat below 5.7 x max catch, half-normal on log K above). Without
  it the fit can run up the K-r ridge (seen in CI: K 60 x larger, r and sdb
  compensating). States `spm_logP` are random effects.
- **Links**: SPR link log D_spr,y ~ N(log P_t, sigma.link) (default on), with
  D_spr = (4h·SPR − (1−h))/(5h−1); optional F link log H_t ~ N(log H_eq(F_y), sigma.Flink).

## Speed
The per-recruit block is large (grid + Newton). If its selectivity depends on
the random F (`sel.weights = "F"`, the F-weighted joint selectivity of the
terminal years), every inner Laplace step goes through it: 0.5 s per objective
evaluation. With `sel.weights = "catch"` (gear selectivities weighted by catch
shares, data) it depends on fixed effects only: 0.03 s per evaluation, a full
fit with SEs about 30 s. The catch weighting is an approximation of the
joint selectivity (F shares are not catch shares). To do: exact fixed-effect
F shares.

## First results, with the K prior (alfonsino OM, `ci/sim/om_sim_spm.R`, CPUE = Gillnet EB, CV 0.2, sigR 0.3)
| scenario | terminal SB/SB0 true | FLicc-SPM | FLicc alone (D_spr) | F/Fmsy true / SPM | conv |
|---|---|---|---|---|---|
| fudc 2001 | 0.130 | 0.127 | 0.137 | 1.33 / 1.72 | ok |
| rebuild 11 | 0.164 | 0.117 | 0.079 | 0.78 / 1.31 | ok |
| rebuild 12 | 0.282 | 0.180 | 0.127 | 0.77 / 1.40 | fail |
| decline 11 | 0.105 | 0.136 | 0.295 | 2.78 / 2.03 | ok |
| decline 12 | 0.185 | 0.188 | 0.252 | 2.74 / 2.42 | ok |
| flat 11 | 0.131 | 0.091 | 0.089 | 1.67 / 2.70 | fail |
| flat 12 | 0.226 | 0.177 | 0.145 | 1.65 / 2.39 | ok |

The joint model removes most of FLicc's lag under rising F (decline). Under
falling F (rebuild) it improves on FLicc alone but stays pessimistic: the SPR
link with sigma.link = 0.2 still pulls towards the lagging D_spr (lag errors of
20-45% in these periods). F/Fmsy is biased high in several runs. Two of seven
fits do not converge.

## Next
1. sigma.link: estimate, or larger / smoothed; compare with the F link.
2. Convergence: starting values (K grid), bounds on r, profile the failures.
3. Exact F shares for the per-recruit selectivity without random-F dependence.
4. Gradient check of the taped Newton steps against finite differences.
5. Regression test of the SPM block against `JABBA3:::jb3_nll` at fixed parameters.
6. Output in the `jabba3` layout for JABBA3's plots and diagnostics.
7. Simulation grid (scenarios x iterations), sexes (FLOSS `om_sex()`).

## Status 16 Nov (work in progress, stopped here)

Changes after the first push (commit "WIP"):
- **F not random by default** (HW: in FLicc F is a nuisance that tunes SPR; the
  random effects belong to the SPM states, as in JABBA3).
- **Per-recruit selectivity from catch-share weights by default**
  (`sel.weights = "catch"`). With F fixed, the F-weighted joint selectivity
  (`"F"`) gives a flat direction: the whole F level of one gear moves with its
  selectivity, K and q (4 of 7 fits failed). The production layer should not
  depend on the nuisance F anyway.
- **AR1 link errors** (`link.rho = c(spr = , F = )`, value fixes, NA
  estimates), rho^gap over gaps between length-data years. Estimated rho on
  the SPR link runs to 1 (trend-only link, process error collapses): needs a
  prior or a fixed value.
- Smooth max (p-norm) for the apical scaling of the per-recruit selectivity.
- Fixed: `pow(rho, gap)` has no derivative at rho = 0 in TMB (products instead).

Partial comparison (F fixed, catch weights; 15 of 28 runs before stopping):
- SPR link, independent errors: 5 of 7 converged; terminal SB/SB0 close in
  decline and fudc runs, pessimistic under rebuild and flat F.
- SPR link, rho = 0.7: little change in depletion; 4 of 7 converged; one
  nonsense F/Fmsy (0.03).
- F link, rho = 0.7: one run (fudc 2001) converged, F/Fmsy 1.51 vs true 1.33
  (SPR link 1.75-1.78), SB/SB0 0.157 vs 0.130.

Next when resumed: finish the F-link (AR1) runs and the combined links; prior
on rho; convergence failures; F/Fmsy bias. Scripts used: `ci/sim/om_sim_spm.R`
and an evaluation loop (seeds 2001/11/12; fudc, rebuild, decline, flat).

## Status 18 Nov

- Rebased on `dev-sex` (74c62b6).
- **Smooth EB/SB curve** (kernel smoother instead of linear interpolation):
  the convergence failures came from the kinks. With it all fits converge
  with the SPR link and with both links (30-50 s per fit). `test_spm.R`
  passes again (terminal SB/SB0 0.133 vs 0.130 true).
- **Link comparison** (F fixed, catch-weighted per-recruit selectivity, 7
  scenarios; `ci/spm/spm_case.R`, results in `ci/spm/results_links_smooth.txt`;
  errors are log(est/true) in the terminal year):

| link | converged | SB/SB0 median error | median abs | F/Fmsy median error | median abs |
|---|---|---|---|---|---|
| SPR | 7/7 | -0.10 | 0.23 | +0.36 | 0.36 |
| F (rho 0.7) | 2/7 | +0.31 | 0.45 | -0.14 | 0.54 |
| SPR + F (rho 0.7 both) | 7/7 | -0.08 | 0.16 | +0.28 | 0.37 |

  The F link alone is not usable: FLicc's equilibrium F does not carry the
  SPM. Added to the SPR link it improves depletion (flat F especially).
  Remaining problems: F/Fmsy biased high under flat and rebuilding F
  (+0.3 to +0.6) and low under rising F; depletion pessimistic under
  rebuilding F (seed 12: -0.44 in all variants).

Next: the F/Fmsy bias (per-recruit BmsyK 0.296 vs 0.275 true; H = C/EB vs
the OM's Fbar), prior on rho, per-recruit gradient check, regression test
against `jb3_nll`.

## Status 18 Nov (evening): where the F/Fmsy offset comes from

Trajectory plot: `ci/spm/plot_status.R`. Biomass (SB/SB0, B/Bmsy) tracks the
OM closely over the whole history; F/Fmsy is about 30% high where biomass is
right.

- **Not a harvest-rate/F conversion.** The SPM's harvest rate C/EB is within
  about 10% of the OM's (same definition, `ci/spm/hmsy_check.R`); the
  reference level is off: Hmsy 0.117 vs 0.136 true.
- **Not the links.** Without any link Hmsy is 0.122; the SPR link is what
  makes biomass precise (no link: terminal SB/SB0 0.335 vs 0.130 true in the
  up-down run; `results_nolink.txt`). Dropping the F link does not change F/Fmsy.
- **FLicc alone vs FLicc-SPM** (`fspr40_compare.R`, length-data years): FLicc
  F/Fspr40 vs the OM's F/Fspr40 median absolute log error 0.75 (equilibrium
  lag: about +0.95 under falling F, -0.85 under rising F); FLicc-SPM F/Fmsy
  vs the OM's F/Fmsy 0.29.
- **Cause: the estimated dome** (`prtrue_check.R`). The per-recruit block at
  the true selectivity gives Hmsy 0.132 vs 0.136 true and r 0.083 vs 0.080
  (the gtg approximation of the age-based OM is good; catch-share weights
  barely matter). At the estimated selectivity Hmsy is 0.109 (-19%): the trawl
  right-hand SD is estimated at 17.6 vs 11.2 true, so too many large fish count
  as exploitable, EB is too high and Hmsy = yield / EB too low.
- **Biomass timing** (`timing = "start"`, new default): EB and SB per recruit
  at the start of the year (mean numbers x Z / (1 - exp(-Z))), as the discrete
  SPM's harvest rate C_t / EB_t. Correct by definition; lowers Hmsy further
  while the dome is too wide. Results (`results_links_start_timing.txt`, median
  log error, median absolute in brackets): SPR + F links B/Bmsy +0.00 (0.17),
  F/Fmsy +0.17 (0.28), 6/7 converged; SPR link alone B/Bmsy -0.32 (0.32),
  F/Fmsy +0.42 (0.42). `timing = "mean"` keeps the previous behaviour.

Next: the right-hand limb. Options: a prior on the dome's right SD; the
production model informing it (catch/CPUE already in the same objective, but
not enough here); an index with logistic selectivity in the length data.
