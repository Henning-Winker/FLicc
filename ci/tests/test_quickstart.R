# Quick-start fit variants: which converge cleanly on the alfonsino data?
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)[, c("len", "year", "data", "qname")]),
                  unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
show <- function(nm, f) {
  cc <- flicc_convergence(f)
  cat(sprintf("%-28s maxgrad %.4g  conv %d  passed %s  SPR %.3f\n", nm, f$max_gradient, f$opt$convergence,
              paste(cc$passed, collapse = ","), sprcur_flicc(f)))
}
show("raw mn default", fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg))
show("raw mn n_restart 3", fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, n_restart = 3))
show("ess200 mn", fiticc(lfdess(lfd, ess.g = 200), stklen, sel_fun = sf, catch_by_gear = cbg))
show("ess200 mn n_restart 3", fiticc(lfdess(lfd, ess.g = 200), stklen, sel_fun = sf, catch_by_gear = cbg, n_restart = 3))
show("ess200 dm", fiticc(lfdess(lfd, ess.g = 200), stklen, sel_fun = sf, catch_by_gear = cbg,
                         settings = list(obs_model = "dm")))
show("raw dm", fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = list(obs_model = "dm")))
