# Where does the time go with F_re = TRUE? (diagnostic, temporary)
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
fitf <- function(..., n_restart = 3) fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"),
  catch_by_gear = c(0.7, 0.3), n_restart = n_restart,
  settings = c(list(CVL = 0.1, obs_model = "dm", prior_sigmaF = c(log(0.5), 0.3, 1)), list(...)))
invisible(fitf())   # warm-up (DLL load)

prof <- function(expr_fun, label) {
  f <- tempfile(); Rprof(f, interval = 0.005)
  t <- system.time(x <- expr_fun())[["elapsed"]]
  Rprof(NULL)
  s <- summaryRprof(f)$by.total
  keep <- c("fiticc", "fiticc_core", "TMB::MakeADFun", "MakeADFun", "nlminb", "TMB::sdreport", "sdreport",
            "obj$fn", "obj$gr", "fn", "gr", "as_FLQuants", "optimHess", "need_refit_flicc", "retape",
            "MakeADFunObject", "eval.parent", "f", "h", "ff", "newton", "LBIspr", "fspr_flicc")
  cat(sprintf("\n==== %s: %.2f s ====\n", label, t))
  print(round(s[rownames(s) %in% paste0('"', keep, '"'), c("total.time", "total.pct")], 2))
  invisible(x)
}
f0 <- prof(function() fitf(), "penalised")
f1 <- prof(function() fitf(F_re = TRUE), "F_re")
f2 <- prof(function() fitf(F_re = TRUE, n_restart = 0), "F_re, n_restart = 0")

# component timings for F_re by hand
tm <- function(e) system.time(e)[["elapsed"]]
o <- f1$obj
cat("\nsingle evaluations on the F_re object:\n")
print(c(fn = tm(o$fn(f1$opt$par)), gr = tm(o$gr(f1$opt$par)),
        sdreport = tm(TMB::sdreport(o)), report = tm(o$report(o$env$last.par.best))))
cat("penalised object:\n")
o0 <- f0$obj
print(c(fn = tm(o0$fn(f0$opt$par)), gr = tm(o0$gr(f0$opt$par)), sdreport = tm(TMB::sdreport(o0))))
cat("\nnlminb iterations / evaluations:\n")
print(rbind(penalised = unlist(f0$opt$evaluations), F_re = unlist(f1$opt$evaluations)))
print(f1$restart_log)
