# Dirichlet-multinomial with estimated precision: sensitivity to the input sample size
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
se_last <- function(f) tail(spr_ci_flicc(f)$se_log, 1)
ess_values <- c(50, 100, 200, 500, 1000, NA)          # NA = raw counts
res <- list()
for (dt in c("fixed", "linear", "saturating")) for (ess in ess_values) {
  d <- if (is.na(ess)) lfd else lfdess(lfd, ess.g = ess)
  t0 <- system.time(f <- tryCatch(fiticc(d, stklen, sel_fun = sf, catch_by_gear = cbg,
                                         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = dt)),
                                  error = function(e) e))[["elapsed"]]
  if (inherits(f, "error")) { cat(dt, ess, "ERROR", conditionMessage(f), "\n"); next }
  ne <- neff_flicc(f)
  ld <- f$obj$env$parList()$log_dm
  res[[length(res) + 1]] <- data.frame(
    dm = dt, ESS = ifelse(is.na(ess), "raw", ess), conv = f$opt$convergence,
    maxgrad = signif(f$max_gradient, 2), SPR_last = round(sprcur_flicc(f), 3),
    se_log = round(se_last(f), 3),
    par_Trawl = round(exp(ld[1]), 3), par_Gillnet = round(exp(ld[2]), 3),
    neff_Trawl = round(mean(ne$neff[ne$gear == "Trawl"]), 1),
    neff_Gillnet = round(mean(ne$neff[ne$gear == "Gillnet"]), 1), sec = round(t0, 1))
  if (dt == "saturating" && identical(ess, 200)) fsat <- f
}
tab <- do.call(rbind, res); print(tab, row.names = FALSE)
# spread of last-year SPR across the input sample sizes, by option
print(aggregate(SPR_last ~ dm, data = tab, FUN = function(v) round(diff(range(v)), 3)))
# residual diagnostics with the estimated precision
r <- resid_flicc(fsat, "osa")
cat("saturating ESS 200: OSA mean", round(mean(r$resid, na.rm = TRUE), 2), "sd", round(sd(r$resid, na.rm = TRUE), 2), "\n")
print(lmean_resid_flicc(fsat)$runs)
print(ppc_flicc(fsat, nsim = 300)$coverage)
# fixed default unchanged
# default (no obs_model given) = Dirichlet-multinomial with estimated (saturating) precision
f0 <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = list(CVL = 0.1))
f1 <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
stopifnot(identical(f0$obs_model, "dm"), f0$tmb_data$dm_type == 2L)
stopifnot(isTRUE(all.equal(f0$opt$objective, f1$opt$objective)))
