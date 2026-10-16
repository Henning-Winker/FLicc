# dome_sens_flicc diagnostics on the alfonsino example
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
lfd200 <- lfdess(lfd, ess.g = 200)
sel_fun <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
fit200 <- fiticc(lfd200, stklen, sel_fun = sel_fun, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = "dm"))
str(fit200$settings)
print(fit200$tmb_data$gear_names); print(fit200$tmb_data$sel_type)
for (se in c(TRUE, FALSE)) {
  ds <- dome_sens_flicc(fit200, gear = "Gillnet", se = se)
  cat("se =", se, "\n"); print(ds$table); print(ds$errors)
}
r <- tryCatch(fiticc(fit200$report$obslen, fit200$stklen, sel_fun = c("dsnormal", "dsnormal"), catch_by_gear = fit200$report$catch_by_gear,
                     settings = fit200$settings, sel_fixed = list(Gillnet = list(rsd = 22))), error = function(e) e)
print(if (inherits(r, "error")) conditionMessage(r) else r$opt$convergence)
if (inherits(r, "error")) print(sys.calls())
# priors
fp <- fiticc(lfd200, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
             settings = list(CVL = 0.1, obs_model = "dm", linf.sd = 0.05, Mk.sd = 0.1, CVL.sd = 0.25))
print(fp$par[fp$par$par %in% c("log_Linf", "log_Mk", "log_Galpha", "log_sigmaF"), ]); print(fp$rep$pdHess)
print(fp$opt$convergence); print(fp$max_gradient)
for (s in list(list(linf.sd = 0.05), list(Mk.sd = 0.1), list(CVL.sd = 0.25), list(linf.sd = 0.05, Mk.sd = 0.1))) {
  f <- fiticc(lfd200, stklen, sel_fun = sel_fun, catch_by_gear = cbg, settings = c(list(CVL = 0.1, obs_model = "dm"), s))
  cat(names(s), "pdHess", f$rep$pdHess, "\n"); print(f$par[f$par$par %in% c("log_Linf", "log_Mk", "log_Galpha"), ])
}
