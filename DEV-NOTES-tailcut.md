# dev-tailcut: changes and test checklist

Branch off `main` at `b3f7879` ("push rfb rule"). Nothing on this branch has
been compiled or run yet: it was written without access to R.

## Changes

### Tail cut in the likelihood (`settings$tail_cut`)

* `fiticc(..., settings = list(tail_cut = 0.9))` excludes length bins with
  lower bound `>= tail_cut * linf` (input Linf, fixed during estimation) from
  the composition likelihood. Observed and predicted compositions are
  truncated and the predicted proportions renormalised over the retained
  bins; the observed total over retained bins is the ESS / precision in the
  nb, multinomial and Dirichlet-multinomial likelihoods. Without a cut the
  likelihood is unchanged.
* `src/FLicc.cpp`: new `DATA_INTEGER(cut_bin)` (0-based first excluded bin,
  `-1` = off).
* `R/flicc_tmb.R`: `settings$tail_cut` (default `NULL`) -> `cut_bin`;
  documented in `?fiticc`. `Lplus` at or above the cut is ignored with a
  warning. Catch likelihood and reported predictions are not truncated.
* `R/fliccit.R`: `cut_bin = -1` for fits created before this change.
* `R/LBIspr.R`: `LBIspr()` and `LBImean()` gain `tail_cut` (default: the
  fit's own cut; `FALSE` = all bins; numeric = fraction of input Linf), and
  return `attr(, "Lcut")`.
* `NAMESPACE`: `export(LBImean)` (was missing).
* Removed the tracked, now stale `src/FLicc.o` and `src/FLicc.dll`.

### MSE functions (`R/flicc_mse.R`)

* Removed the duplicated `lbi.hcr()` definition (two identical copies).
* `rfb.flicc.hcr()`: default `index = "SPR"`.
* `lbi.hcr()` and `rfb.flicc.hcr()`: new `trend = c("fit", "track")`.
  `"track"` computes the n1-over-n2 trend from the terminal-year values of
  previous cycles (tracked as `index.hcr`) plus the current one, instead of
  the current fit's series (end-of-window noise). Default `"fit"` keeps the
  previous behaviour. Both HCRs now track `index.hcr` every cycle.
* `flicc.is()`: new `Cmin`, a floor on the TAC the multiplier is applied to,
  so a multiplicative rule cannot get stuck near zero after a crash.

### Utilities (`R/utils.R`)

* `lfdess()`: empty gear-years give 0 instead of NaN; `ess.g` matched to
  gears by name (or position); documentation fixed (gears are list
  elements, not units).

## Test checklist

```r
devtools::install()          # recompiles src/FLicc.cpp (no stale .o/.dll)
library(FLicc)
data(alfonsino)
stklen <- stocklen(lfd_alfonsino, lhpar_alfonsino)

## 1. no cut: identical to main
f0 <- fiticc(lfd_alfonsino, stklen, sel_fun = c("dsnormal", "logistic"),
             catch_by_gear = c(0.7, 0.3))
## compare with the same fit on main: logLik(f0), f0$report$spr, f0$opt$par

## 2. with a cut
f9 <- fiticc(lfd_alfonsino, stklen, sel_fun = c("dsnormal", "logistic"),
             catch_by_gear = c(0.7, 0.3), settings = list(tail_cut = 0.9))
f9$tmb_data$cut_bin                      # first excluded bin (0-based)
c(f0$report$spr[, ac(dims(f0$report$spr)$maxyear)],
  f9$report$spr[, ac(dims(f9$report$spr)$maxyear)])

## 3. indicators use the same bins
l9 <- LBIspr(f9); attr(l9, "Lcut"); attr(l9, "Lref")   # Lref < Lcut
m9 <- LBImean(f9); attr(m9, "Lcut")
attr(LBIspr(f9, tail_cut = FALSE), "Lcut")            # NA

## 4. all observation models run with a cut
for (om in c("nb", "mn", "dm"))
  fiticc(lfd_alfonsino, stklen, sel_fun = c("dsnormal", "logistic"),
         catch_by_gear = c(0.7, 0.3),
         settings = list(tail_cut = 0.9, obs_model = om))

## 5. reference points on a cut fit and on an old saved fit
fspr_flicc(f9, spr = 40); brp_tmb_flicc(f9, spr = 40)

## 6. lfdess with an empty gear-year
x <- lfd_alfonsino; x[[1]][, 1] <- 0
any(is.nan(unlist(lapply(lfdess(x, 100), c))))        # FALSE

## 7. MP smoke test (FLicc_lfd_mp_test.Rmd), short horizon
##    est: flicc.sa(settings = list(..., tail_cut = 0.9))
##    hcr: lbi.hcr(gear = "SPR", trend = "track") or
##         rfb.flicc.hcr(index = "SPR", trend = "track", f = "Zrel")
##    isys: flicc.is(initac = initac, Cmin = 0.1 * initac)
##    check tracking: index.hcr, trend.hcr / r.hcr, Lref.<gear>, Lmeanref.<gear>
```
