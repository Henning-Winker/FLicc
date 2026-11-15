# FLicc-SPM prototype (branch dev-spm, from dev-sex)

One TMB objective: FLicc length likelihood (gtg, one or two sexes, F as random
effects) + a Pella-Tomlinson state-space surplus production model in spawning
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
  multiples of M at Linf), linear interpolation in SB/SB0 (CondExp). Column 0
  for removals (fishery), one per index (gear selectivity, joint fishery, or SB).
- **SPM**: translation of the subset of `JABBA3::jb3_nll()` (Selex state):
  constant process error, catch without error, production lag, soft bounds on
  P, JABBA priors on bkinit, sdb, sdi. States `spm_logP` are random effects.
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

## First results (alfonsino OM, `ci/sim/om_sim_spm.R`, CPUE = Gillnet EB, CV 0.2, sigR 0.3)
| scenario | terminal SB/SB0 true | FLicc-SPM | FLicc alone (D_spr) | F/Fmsy true / SPM | conv |
|---|---|---|---|---|---|
| fudc 2001 | 0.130 | 0.127 | 0.137 | 1.33 / 1.72 | ok |
| rebuild 11 | 0.164 | 0.117 | 0.079 | 0.78 / 1.31 | ok |
| rebuild 12 | 0.282 | 0.180 | 0.127 | 0.77 / 1.40 | fail |
| decline 11 | 0.105 | 0.136 | 0.295 | 2.78 / 2.03 | ok |
| decline 12 | 0.185 | 0.188 | 0.252 | 2.74 / 2.42 | ok |
| flat 11 | 0.131 | 0.082 | 0.089 | 1.67 / 0.05 | fail |
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
