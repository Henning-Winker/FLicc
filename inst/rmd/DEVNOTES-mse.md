
# Developer notes: FLicc MSE modules

Companion to `inst/rmd/mse-01-flicc-modules.Rmd`. Covers the documentation
set-up, issues found in the current code (`R/flicc_mse.R`, `R/LBIspr.R`,
`R/flicc_tmb.R`) while writing the guideline, and proposed changes for speed.

## 1. Documentation set-up (as in FLRef and cpueTargeting)

```
FLicc/
  inst/rmd/mse-01-flicc-modules.Rmd   # sources
  inst/rmd/mse-02-alfonsino-mse.Rmd   # (next) OM from FLicc_MSE_lite_OM_v3 + sampling.lfd.oem + flicc.sa
  inst/rmd/mse-03-robustness.Rmd      # (next) misspecified M, Linf, CVL, selectivity
  inst/rmd/DEVNOTES-mse.md            # this file
  inst/scripts/render_docs.R          # renders inst/rmd -> docs/
  docs/index.html                     # landing page (hand-written)
  docs/.nojekyll
  rdata/                              # cached MSE runs (not pushed)
```

* `.gitignore`: `rdata/`, `inst/rmd/*.html`, `inst/rmd/rdata/` (added).
* `.Rbuildignore`: `^docs$`, `^rdata$`, `^inst/rmd$`, `^inst/scripts$` (added).
* GitHub: Settings > Pages > Deploy from a branch > `main` / `/docs`
  (after merging `docs-dev`; for previews, select `docs-dev` / `/docs`).
* Render from the package root after installing the current `R/`:
  `source("inst/scripts/render_docs.R")`.
* The existing `docs/*.Rmd` (intro vignette, structural uncertainty, 1.0.9
  updates) could move to `inst/rmd/` later so that all sources live in one
  place; they are linked from `docs/index.html` as they are.
* `docs/mse-01-flicc-modules.html` in the first commit is a **preview**
  converted with pandoc without running R (no output, no figures). Replace it
  by running `render_docs.R`.

## 2. Issues found in the current code

1. **`lbi.hcr()` is defined twice** in `R/flicc_mse.R` (lines ~903 and ~983),
   with identical code and roxygen blocks. The second wins at collation and
   roxygen writes the Rd twice. Delete one copy.

2. **Opposite `decision.hcr` conventions.** `lbi.hcr()` sets
   `decision.hcr = 3` (limits applied in `flicc.is()`) when SPR is *below*
   `spr_trigger`; `rfb.flicc.hcr()` sets it when the safeguard is *inactive*
   (b = 1, SPR above the trigger), as in the ICES rfb rule. Both are
   defensible but should be documented side by side in `?flicc.is`, or made
   explicit with an argument (e.g. `limits = c("always", "above", "below")`).

3. **`prior_sigmaF` with a zero element switches the random walk off.**
   `fiticc_core()` uses `if (any(settings$prior_sigmaF == FALSE))`, so a prior
   mean of `log(1) = 0` (or a zero sd) silently sets
   `prior_sigmaF = c(NA, NA, 0)`. Use
   `if (isFALSE(settings$prior_sigmaF) || identical(settings$prior_sigmaF, 0))`
   or check `length(settings$prior_sigmaF) == 1`.

4. **With the random walk off, `log_sigmaF` is a free, unidentified
   parameter.** When `prior_sigmaF_use = 0`, `log_sigmaF` does not enter the
   objective, but it is only mapped out (`map$log_sigmaF <- factor(NA)`) for
   a single year. Its gradient is zero and the Hessian is singular, which can
   flag otherwise good fits as not converged in `need_refit_flicc()` /
   `flicc_convergence()`. Map it whenever `prior_sigmaF[3] == 0`.

5. **Stale statement in `docs/FLicc_updates_1.0.9.Rmd`**: "`flicc.sa()` sets
   `spr_ref` to its `spr` target automatically". The current `flicc.sa()`
   passes `settings` unchanged; it only works because both default to 40.
   Either add `settings$spr_ref <- unique(c(spr, settings$spr_ref))` in
   `flicc.sa()` (preferred) or correct the text.

6. **`thresh` defaults differ**: `LBIspr()`, `flicc.sa()`, `hcr_sprlbi()` and
   `plot_LBIp()` use 0.75, `plot_LBIspr()` uses 0.8. Align to 0.75.

7. **Roxygen names**: the block above `plot_hcrspr()` is titled
   `plot_hcrflicc` and its example calls `plot_hcrflicc()`; the
   `hcr_sprlbi()` examples call `hcr_sprr()`. The `spr_rule()` title says
   `hcr_spr (flicc)`.

8. **Seasonal OEM output is not pooled.** `flicc.sa()` applies `unitSums()`
   but not `seasonSums()` to `attr(stk, "lfd")`. With
   `sampling.lfd.oem(season = TRUE)` the fit would receive seasonal
   compositions. Add `seasonSums()` (or stop with a message).

9. **`need_refit_flicc(fit, grad_tol)` ignores `grad_tol`.** It checks the
   optimiser code, a finite objective and `pdHess`, but not the gradient,
   although `fiticc()` passes `grad_tol` to it. Either check
   `fit$max_gradient < grad_tol` or drop the argument.

10. **`browser()` comments** remain in `flicc.sa()`, `lbi.hcr()`,
   `flicc.is()`; harmless but should be removed before release.

## 3. Proposed changes for speed

1. **Fix `sigmaF` through the map.** Add a setting, e.g.
   `settings$sigmaF_fixed = TRUE` (or accept `prior_sigmaF = c(mu, 0, 1)`),
   that maps `log_sigmaF` (`map$log_sigmaF <- factor(NA)`) at its initial
   value, which is already `prior_sigmaF_mean`. The prior term must then be
   skipped (sd = 0 gives `-Inf` in `dnorm()`): set `prior_sigmaF_sd = 1` in
   `tmb_data` when fixed, since the term is constant. The random-walk
   penalty itself stays. This replaces the tight-prior workaround in the
   guideline.

2. **Per-iteration `sigmaF` in `flicc.sa()`.** Accept `sigmaF` as a numeric
   vector of length `it` (from a conditioning step before `mp()`), and set
   `settings$prior_sigmaF[1] <- log(sigmaF[i])` inside `fit_one(i)`.

3. **Warm starts.** With a fixed `nyrs` window the number of estimated
   parameters is constant across cycles, so `fit$opt$par` of the previous
   cycle is a valid `start` for the next one (the log-F random walk shifts by
   `frq` years; shifting `log_Fk` rows or re-initialising them only is a
   refinement). Store `opt$par` by iteration in an environment or in the
   tracking object and pass it to `fiticc(start = )`. Fall back to the
   default start when the previous fit was unusable.

4. **Conditioning helper.** A function
   `flicc_condition(oem, lhpar, sel_fun, settings, years)` that fits the
   historical OEM samples once per iteration and returns `sigmaF`, the
   starting parameters and the convergence table, to be used in (2) and (3).

5. **Profile the remaining cost.** With Fspr in TMB, the time of a cycle
   should be dominated by `nlminb()`. Worth timing `stocklen()` per
   iteration (built from scratch in every cycle although `lhpar` is static;
   only the length bins of `lfd` would change) and `sdreport()`. The
   indicators do not need standard errors, but `need_refit_flicc()` uses
   `rep$pdHess`; a cheaper option inside `flicc.sa()` is to replace
   `sdreport()` by a Cholesky check of `obj$he()` at the optimum (or
   `optimHess()`), keeping the convergence screen without the delta-method
   step.
