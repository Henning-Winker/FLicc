suppressPackageStartupMessages(library(FLicc, lib.loc = "lib_new"))
cat("FLicc", as.character(packageVersion("FLicc", lib.loc = "lib_new")), "\n")
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
fail <- FALSE
rel <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)) / pmax(abs(as.numeric(b)), 1e-12), na.rm = TRUE)
for (pm in c("gtg", "gamma")) for (mm in c("constant", "inverse")) {
  stklen <- stocklen(lfd, lhpar, m_model = mm)
  set <- list(prior_sigmaF = c(log(0.5), 0.3, 1), CVL = 0.1, pop_model = pm,
              obs_model = if (pm == "gtg") "mn" else "nb")
  fit <- fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3), settings = set)
  cat("\n==========", pm, "/ M", mm, "==========\n")
  run <- function(fast) {
    options(FLicc.brp_tmb = fast)
    t <- list(); v <- list()
    t$fspr40   <- system.time(v$fspr40   <- fspr_flicc(fit, spr = 40))["elapsed"]
    t$fspr20FM <- system.time(v$fspr20FM <- fspr_flicc(fit, spr = 20, input = "FM"))["elapsed"]
    t$spr      <- system.time(v$spr      <- spr_flicc(fit, F = 0.1))["elapsed"]
    t$pr       <- system.time(pr <- pr_flicc(fit, F = 0.1))["elapsed"]
    v$ypr <- pr$YPR; v$sbpr <- pr$SBPR; v$prN <- c(pr$N); v$prC <- c(pr$Cn)
    t$nf       <- system.time(v$nf <- c(nf_flicc(fit, F = 0.05)))["elapsed"]
    t$prbrp    <- system.time(pb <- prbrp_flicc(fit))["elapsed"]
    v$pbY <- c(pb$YPR); v$pbS <- c(pb$SBPR); v$pbN <- c(pb$NPR); v$spr0 <- pb$spr0
    t$eqstklen <- system.time(eq <- eqstklen(fit, s = 0.75))["elapsed"]
    v$eqref <- c(attr(eq, "refpts")); v$eqcatch <- c(catch(eq))
    t$LBIspr   <- system.time(lb <- LBIspr(fit))["elapsed"]
    v$lbi <- unlist(lapply(lb, function(x) c(index(x))))
    t$stkr     <- system.time(st <- flicc2FLStockR(fit, rel = TRUE))["elapsed"]
    v$stkr <- c(unlist(as.list(st@refpts)), c(st@stock), c(st@harvest))
    list(t = unlist(t), v = v)
  }
  old <- run(FALSE); new <- run(TRUE)
  cat(sprintf("%-10s %9s %9s %8s   %s\n", "function", "old(s)", "new(s)", "speedup", "max rel diff"))
  map <- list(fspr40 = "fspr40", fspr20FM = "fspr20FM", spr = "spr", pr = c("ypr","sbpr","prN","prC"),
              nf = "nf", prbrp = c("pbY","pbS","pbN","spr0"), eqstklen = c("eqref","eqcatch"),
              LBIspr = "lbi", stkr = "stkr")
  for (f in names(map)) {
    d <- max(sapply(map[[f]], function(k) rel(new$v[[k]], old$v[[k]])))
    cat(sprintf("%-10s %9.3f %9.3f %7.0fx   %.2e\n", f, old$t[f], new$t[f], old$t[f] / max(new$t[f], 1e-3), d))
    if (!is.finite(d) || d > 1e-3) { fail <- TRUE; cat("   ^ FAIL\n") }
  }
  b <- brp_tmb_flicc(fit, spr = c(40, 30, 20, 10))
  cat("brp_tmb_flicc Fspr40/30/20/10:", round(b$Fspr, 5), "\n")
}
options(FLicc.brp_tmb = TRUE)
if (fail) quit(status = 1)
