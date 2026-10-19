# mc_flicc: convergence flags, sel_fixed carried over, NA for failed draws, equal weights
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
gd <- names(lfd)[sf == "dsnormal"]          # the double-normal gear
sfx <- setNames(list(list(rsd = 30)), gd)
fit <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = list(CVL = 0.1),
              sel_fixed = sfx)
stopifnot(identical(fit$sel_fixed[[gd]]$rsd, 30))
t0 <- Sys.time()
mc <- mc_flicc(fit, nsim = 12, cv = list(linf = 0.1, k = 0.1, M = 0.2, L50 = 0.1, CVL = 0.2), seed = 1)
cat("time:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
print(table(mc$conv)); print(mc$ok)
stopifnot(is.logical(mc$ok), !anyNA(mc$ok), !anyNA(mc$conv))
# rsd stayed fixed in every successful refit
print(mc$selpars[[gd]])
rsd <- as.numeric(mc$selpars[[gd]]["rsd", ])[mc$ok]
print(rsd)
stopifnot(all(abs(rsd - 30) < 1e-6))
# failed / non-converged draws carry NA, not copies of the reference fit
stopifnot(all(is.na(mc$logLik[!mc$ok])), all(is.finite(mc$logLik[mc$ok])))
w_aic <- mc_aic_weights(mc); w_eq <- mc_aic_weights(mc, weights = "equal")
stopifnot(abs(sum(w_aic$weight) - 1) < 1e-8, abs(sum(w_eq$weight) - 1) < 1e-8,
          all(w_eq$weight[!mc$ok] == 0), length(unique(round(w_eq$weight[mc$ok], 10))) == 1)
print(mc_spr_tab(mc, dAIC_cut = 1e6, weights = "equal"))
print(mc_spr_tab(mc, dAIC_cut = 1e6, weights = "aic"))
stopifnot(nrow(mc$draws) == 12)
