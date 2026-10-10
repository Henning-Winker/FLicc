# F random walk as random effects (settings$F_re)
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
fitf <- function(...) fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
                             settings = c(list(CVL = 0.1), list(...)))
summ <- function(f) c(nll = round(f$opt$objective, 2), conv = f$opt$convergence,
                      maxgrad = signif(f$opt$max_gradient, 2), pdHess = f$rep$pdHess,
                      sigmaF = round(exp(f$par$mpd[f$par$par == "log_sigmaF"][1]), 3),
                      sigmaF_se_log = round(f$par$se[f$par$par == "log_sigmaF"][1], 3),
                      SPR_last = round(tail(c(f$report$spr), 1), 4))

## 1. penalised (default) vs random effects, three observation models
res <- list()
for (om in c("mn", "dm", "nb")) {
  t0 <- system.time(f0 <- fitf(obs_model = om, prior_sigmaF = c(log(0.5), 0.3, 1)))[["elapsed"]]
  t1 <- system.time(f1 <- fitf(obs_model = om, prior_sigmaF = c(log(0.5), 0.3, 1), F_re = TRUE))[["elapsed"]]
  t2 <- system.time(f2 <- fitf(obs_model = om, prior_sigmaF = c(log(0.5), NA, 1), F_re = TRUE))[["elapsed"]]
  res[[om]] <- list(pen = f0, re = f1, re_noprior = f2)
  print(rbind(penalised = c(summ(f0), secs = t0), RE_prior = c(summ(f1), secs = t1),
              RE_noprior = c(summ(f2), secs = t2)))
  stopifnot(isTRUE(f1$F_re), !isTRUE(f0$F_re), length(f1$opt$par) < length(f0$opt$par))
}

## SPR by year: penalised vs random effects (dm)
print(round(rbind(penalised = c(res$dm$pen$report$spr), RE = c(res$dm$re$report$spr),
                  RE_noprior = c(res$dm$re_noprior$report$spr)), 3))

## SE of log SPR (last year) includes the F process with F_re
se_last <- function(f) { s <- summary(f$rep, "report"); tail(s[rownames(s) == "log_spr_y", 2], 1) }
round(c(penalised = se_last(res$dm$pen), RE = se_last(res$dm$re)), 4)

## 2. report is at the random-effect mode
f1 <- res$dm$re
r <- f1$obj$report(f1$obj$env$last.par.best)
stopifnot(max(abs(r$spr_y - c(f1$report$spr))) < 1e-10)

## 3. downstream functions on an F_re fit
b <- brp_tmb_flicc(f1, spr = c(40, 20))
stopifnot(abs(fspr_flicc(f1, spr = 40) - b$Fspr[1]) < 1e-8)
options(FLicc.brp_tmb = FALSE); f40_old <- fspr_flicc(f1, spr = 40); options(FLicc.brp_tmb = TRUE)
c(Fspr40_tmb = b$Fspr[1], Fspr40_R = f40_old)
stopifnot(abs(b$Fspr[1] / f40_old - 1) < 1e-3)
invisible(LBIspr(f1)); invisible(flicc2FLStockR(f1, rel = TRUE, s = 0.75))
invisible(eqstklen(f1, s = 0.75)); invisible(selpars_flicc(f1))
print(flicc_convergence(f1))
print(LLflicc(f1)[[1]])

## 4. restart from a previous F_re fit
f1b <- fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
              settings = list(CVL = 0.1, obs_model = "dm", prior_sigmaF = c(log(0.5), 0.3, 1), F_re = TRUE),
              start = f1$opt$par)
stopifnot(abs(f1b$opt$objective - f1$opt$objective) < 1e-4)

## 5. gamma model with F_re
fg <- fitf(pop_model = "gamma", obs_model = "nb", prior_sigmaF = c(log(0.5), NA, 1), F_re = TRUE)
print(summ(fg))

## 6. single year: F_re ignored with a warning
w <- tryCatch(fiticc(window(lfd, start = 2024, end = 2024), window(stklen, start = 2024, end = 2024),
                     sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
                     settings = list(F_re = TRUE)), warning = function(w) conditionMessage(w))
print(w)

## 7. two-step: sigmaF estimated once with F_re, then fixed (fast)
sig <- exp(f1$par$mpd[f1$par$par == "log_sigmaF"])
f3 <- fitf(obs_model = "dm", prior_sigmaF = c(log(0.5), 0.3, 1), sigmaF = sig)
stopifnot(!"log_sigmaF" %in% names(f3$opt$par))
d3 <- max(abs(c(f3$report$spr) - c(f1$report$spr)))
cat(sprintf("sigmaF fixed at %.3f: max |SPR - SPR(F_re)| = %.4f\n", sig, d3))
stopifnot(d3 < 0.01)

## 8. se = FALSE skips sdreport; downstream still works
fse <- fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
              settings = list(CVL = 0.1, obs_model = "dm", sigmaF = sig), se = FALSE)
stopifnot(is.null(fse$rep), max(abs(c(fse$report$spr) - c(f3$report$spr))) < 1e-8)
invisible(LBIspr(fse)); invisible(fspr_flicc(fse)); invisible(flicc2FLStockR(fse))
print(tryCatch(need_refit_flicc(fse), error = function(e) conditionMessage(e)))

## 9. timing (dm)
tm <- function(...) system.time(fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"),
  catch_by_gear = c(0.7, 0.3), ...))[["elapsed"]]
base <- list(CVL = 0.1, obs_model = "dm", prior_sigmaF = c(log(0.5), 0.3, 1))
tt <- c(penalised          = tm(settings = base),
        penalised_se_FALSE = tm(settings = base, se = FALSE),
        F_re               = tm(settings = c(base, list(F_re = TRUE))),
        F_re_se_FALSE      = tm(settings = c(base, list(F_re = TRUE)), se = FALSE),
        sigmaF_fixed       = tm(settings = c(base, list(sigmaF = sig))),
        sigmaF_fixed_se_FALSE = tm(settings = c(base, list(sigmaF = sig)), se = FALSE))
print(round(tt, 2))
cat("\nF_re tests passed.\n")
