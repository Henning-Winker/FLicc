# Speed options for F_re (diagnostic, temporary)
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
fitf <- function(..., n_restart = 3) fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"),
  catch_by_gear = c(0.7, 0.3), n_restart = n_restart,
  settings = c(list(CVL = 0.1, obs_model = "dm"), list(...)))
invisible(fitf())
tm <- function(e) system.time(e)[["elapsed"]]
lastspr <- function(f) round(tail(c(f$report$spr), 1), 4)

t0 <- tm(f0 <- fitf(prior_sigmaF = c(log(0.5), 0.3, 1)))
t1 <- tm(f1 <- fitf(prior_sigmaF = c(log(0.5), 0.3, 1), F_re = TRUE))
sig <- exp(f1$par$mpd[f1$par$par == "log_sigmaF"])
cat(sprintf("\npenalised %.2f s   F_re %.2f s   sigmaF(F_re) = %.3f\n", t0, t1, sig))

## A. inner optimiser settings, same object, restart nlminb from the initial values
o <- f1$obj; p0 <- o$par
run_inner <- function(label, ...) {
  if (length(list(...))) TMB::newtonOption(o, ...)
  o$env$last.par <- o$env$par; o$env$last.par.best <- o$env$par
  t <- tm(op <- nlminb(p0, o$fn, o$gr, control = list(iter.max = 1000, eval.max = 2000, rel.tol = 1e-9, x.tol = 1e-7)))
  cat(sprintf("%-38s nlminb %.2f s  obj %.4f  evals %d/%d\n", label, t, op$objective,
              op$evaluations[1], op$evaluations[2]))
}
cat("\nA. inner Newton options (one nlminb from start, no restarts / sdreport):\n")
run_inner("default")
run_inner("smartsearch = FALSE", smartsearch = FALSE)
run_inner("smartsearch = FALSE, maxit = 20", smartsearch = FALSE, maxit = 20)
run_inner("tol10 = 1e-3 (looser)", smartsearch = TRUE, maxit = 1000, tol10 = 1e-3)

## B. cost of sdreport
cat(sprintf("\nB. sdreport: F_re %.2f s, penalised %.2f s\n", tm(TMB::sdreport(f1$obj)), tm(TMB::sdreport(f0$obj))))

## C. two-step: estimate sigmaF once with F_re, then penalised fits with sigmaF fixed
t3 <- tm(f3 <- fitf(prior_sigmaF = c(log(sig), 1e-4, 1)))
cat(sprintf("\nC. penalised with sigmaF fixed at %.3f: %.2f s\n", sig, t3))
print(round(rbind(penalised_prior = c(f0$report$spr), F_re = c(f1$report$spr),
                  penalised_sigmaF_fixed = c(f3$report$spr)), 3))
cat(sprintf("max |SPR(F_re) - SPR(two-step)| = %.4f\n", max(abs(c(f1$report$spr) - c(f3$report$spr)))))
