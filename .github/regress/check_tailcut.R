# Test dev-tailcut against main (lib_base = main, lib_new = dev-tailcut)
args <- commandArgs(trailingOnly = TRUE); lib <- args[1]; out <- args[2]
suppressPackageStartupMessages(library(FLicc, lib.loc = lib))
cat("\n####", lib, "FLicc", as.character(packageVersion("FLicc", lib.loc = lib)), "\n")
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
fitf <- function(...) fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
                             settings = c(list(prior_sigmaF = c(log(0.5), 0.3, 1), CVL = 0.1), list(...)))
last <- function(f) round(tail(c(f$report$spr), 1), 4)
res <- list()
for (om in c("mn", "dm", "nb")) {
  f <- fitf(obs_model = om)
  res[[paste0("nocut_", om)]] <- list(obj = f$opt$objective, spr = c(f$report$spr), par = f$opt$par,
                                      obj0 = f$obj$fn(f$obj$par), fspr = fspr_flicc(f))
  cat(sprintf("no cut %-3s nll %.4f  SPR %.4f  Fspr40 %.5f\n", om, f$opt$objective, last(f), fspr_flicc(f)))
}
if (packageVersion("FLicc", lib.loc = lib) >= "0" && "cut_bin" %in% names(res) == FALSE &&
    !is.null(formals(FLicc:::data_tmb_flicc))) {}
new <- grepl("new", lib)
if (new) {
  cat("\n--- tail_cut tests ---\n")
  for (x in c(1.0, 0.9)) for (om in c("mn", "dm", "nb")) {
    f <- tryCatch(fitf(obs_model = om, tail_cut = x), error = function(e) e)
    if (inherits(f, "error")) { cat("ERROR tail_cut", x, om, ":", conditionMessage(f), "\n"); next }
    cb <- f$tmb_data$cut_bin
    cat(sprintf("tail_cut %.1f %-3s cut_bin %d (Lcut %s cm)  nll %.3f  conv %d  SPR %.4f  Fspr40 %.5f\n",
                x, om, cb, f$tmb_data$LLB[cb + 1], f$opt$objective, f$opt$convergence, last(f),
                tryCatch(fspr_flicc(f), error = function(e) NA)))
    if (x == 0.9 && om == "dm") f9 <- f
  }
  cat("\n--- indicators ---\n")
  l9 <- tryCatch(LBIspr(f9), error = function(e) e); print(if (inherits(l9, "error")) conditionMessage(l9) else
    c(Lcut = attr(l9, "Lcut"), Lref = attr(l9, "Lref")))
  l0 <- tryCatch(LBIspr(f9, tail_cut = FALSE), error = function(e) e); print(if (inherits(l0, "error")) conditionMessage(l0) else
    c(Lcut_nocut = attr(l0, "Lcut")))
  m9 <- tryCatch(LBImean(f9), error = function(e) e); print(if (inherits(m9, "error")) conditionMessage(m9) else
    c(Lcut_mean = attr(m9, "Lcut")))
  if (!inherits(l9, "error")) print(round(c(index(l9[[1]])), 3))
  cat("\n--- reference points on cut fit ---\n")
  print(c(fspr = fspr_flicc(f9, spr = 40), brp = brp_tmb_flicc(f9, spr = 40)$Fspr))
  cat("\n--- lfdess with empty gear-year ---\n")
  x <- lfd; x[[1]][, 1] <- 0
  print(any(is.nan(unlist(lapply(lfdess(x, 100), c)))))
  st <- tryCatch(flicc2FLStockR(f9, rel = TRUE, s = 0.75), error = function(e) e)
  cat("flicc2FLStockR on cut fit:", if (inherits(st, "error")) conditionMessage(st) else "ok", "\n")
}
saveRDS(res, out)
