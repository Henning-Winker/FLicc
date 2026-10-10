# Regression of the branch against main on the alfonsino example.
# Fits a fixed set of default cases with both builds (separate processes) and
# compares (a) objective, SPR and predicted LFD at identical parameter values
# and (b) fitted results. Differences are reported, not treated as failures,
# unless a case errors or default results change beyond optimiser tolerance.
# Usage: Rscript ci/compare_main.R <repo> <lib_base> <lib_new> <outdir>
args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]; lib_b <- normalizePath(args[2]); lib_n <- normalizePath(args[3]); out <- args[4]

runner <- file.path(tempdir(), "cases.R")
writeLines(c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "suppressPackageStartupMessages(library(FLicc, lib.loc = args[1]))",
  "data('alfonsino')",
  "lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = 'cm', midL = FALSE)",
  "lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)",
  "stklen <- stocklen(lfd, lhpar, m_model = 'constant')",
  "base <- list(prior_sigmaF = c(log(0.5), 0.3, 1), CVL = 0.1)",
  "cases <- list(gtg_mn = list(obs_model = 'mn'), gtg_dm = list(obs_model = 'dm'),",
  "              gtg_nb = list(obs_model = 'nb'),",
  "              gamma_nb = list(pop_model = 'gamma', obs_model = 'nb'))",
  "res <- lapply(cases, function(cs) tryCatch({",
  "  f <- fiticc(lfd, stklen, sel_fun = c('dsnormal', 'logistic'), catch_by_gear = c(0.7, 0.3),",
  "              settings = c(base, cs))",
  "  r0 <- f$obj$report(f$obj$par)",
  "  list(obj = f$opt$objective, conv = f$opt$convergence, spr = c(f$report$spr),",
  "       fspr = fspr_flicc(f), obj0 = f$obj$fn(f$obj$par), spr0 = as.numeric(r0$spr_y),",
  "       plen0 = as.numeric(r0$plen))",
  "}, error = function(e) list(error = conditionMessage(e))))",
  "saveRDS(res, args[2])"), runner)

rb <- file.path(out, "cmp_main.rds"); rn <- file.path(out, "cmp_branch.rds")
log <- file.path(out, "compare_main.log")
system2("Rscript", c(runner, lib_b, rb), stdout = log, stderr = log)
system2("Rscript", c(runner, lib_n, rn), stdout = log, stderr = log)
if (!file.exists(rb) || !file.exists(rn)) {
  cat("compare_main=1\n", file = file.path(out, "status.txt"), append = TRUE); quit(status = 0)
}
b <- readRDS(rb); n <- readRDS(rn)
lines <- c("| case | same-params diff (obj, SPR, LFD) | fitted rel. obj diff | SPR (last yr) main / branch | conv | result |",
           "|---|---|---|---|---|---|")
fail <- FALSE
for (k in names(b)) {
  if (!is.null(b[[k]]$error) || !is.null(n[[k]]$error)) {
    lines <- c(lines, sprintf("| %s | | | | | ERROR: %s %s |", k, b[[k]]$error, n[[k]]$error)); fail <- TRUE; next
  }
  same <- c(abs(b[[k]]$obj0 - n[[k]]$obj0) / abs(b[[k]]$obj0),
            if (length(b[[k]]$spr0) == length(n[[k]]$spr0)) max(abs(b[[k]]$spr0 - n[[k]]$spr0)) else NA,
            if (length(b[[k]]$plen0) == length(n[[k]]$plen0)) max(abs(b[[k]]$plen0 - n[[k]]$plen0)) else NA)
  fitd <- abs(b[[k]]$obj - n[[k]]$obj) / abs(b[[k]]$obj)
  res <- if (all(is.finite(same)) && all(same < 1e-10) && fitd < 1e-8) "unchanged" else "CHANGED"
  if (res == "CHANGED") fail <- TRUE
  lines <- c(lines, sprintf("| %s | %s | %.1e | %.4f / %.4f | %d / %d | %s |", k,
                            paste(signif(same, 2), collapse = ", "), fitd,
                            tail(b[[k]]$spr, 1), tail(n[[k]]$spr, 1), b[[k]]$conv, n[[k]]$conv, res))
}
writeLines(lines, file.path(out, "compare_main.md"))
cat(lines, sep = "\n")
# a change in default results is reported as 'changed' (exit code 2), not a failure
cat(sprintf("compare_main=%d\n", if (fail) 2 else 0), file = file.path(out, "status.txt"), append = TRUE)
