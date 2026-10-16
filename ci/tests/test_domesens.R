# dome_sens_flicc and life-history priors on the alfonsino example
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
lfd200 <- lfdess(lfd, ess.g = 200)
sel_fun <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
for (g in names(lfd)) { q <- lfd[[g]]; L <- as.numeric(dimnames(q)$len); n <- apply(q@.Data, 1, sum)
  cat(g, "largest bin with fish:", max(L[n > 0]), " n in bins >= 60:", sum(n[L >= 60]), "\n") }
fs <- function(...) fiticc(lfd200, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                           settings = modifyList(list(CVL = 0.1, obs_model = "dm"), list(...)))
fit200 <- fs()
print(fit200$par[grepl("Sm", fit200$par$par), ])
cat("base pdHess", fit200$rep$pdHess, "maxgrad", fit200$max_gradient, "SPR", tail(c(fit200$report$spr), 1), "\n")

# --- where does sdreport fail? ---
nonfin <- function(f) {
  obj <- f$obj; p0 <- f$opt$par
  g0 <- obj$gr(p0); cat("fn", obj$fn(p0), " non-finite gr at opt:", which(!is.finite(g0)), "\n")
  bad <- integer()
  for (i in seq_along(p0)) for (h in c(-1e-3, 1e-3)) {
    p <- p0; p[i] <- p[i] + h; gi <- obj$gr(p)
    if (any(!is.finite(gi)) || !is.finite(obj$fn(p))) bad <- c(bad, i)
  }
  bad <- unique(bad)
  if (length(bad)) print(data.frame(i = bad, name = names(p0)[bad], value = p0[bad]))
  else cat("all perturbed gradients finite\n")
  print(round(cbind(par = p0, gr = as.numeric(g0)), 4)[grepl("Sm|sigma|phi", names(p0)), ])
}
sf_base <- c("logistic", "dsnormal", "normal")[fit200$tmb_data$sel_type]
r0 <- tryCatch(fiticc(fit200$report$obslen, fit200$stklen, sel_fun = sf_base, catch_by_gear = fit200$report$catch_by_gear,
                      settings = fit200$settings), error = function(e) conditionMessage(e))
cat("refit via report route, no change:", if (is.character(r0)) r0 else "ok", "\n")
r1 <- tryCatch(fiticc(lfd200, stklen, sel_fun = sf_base, catch_by_gear = cbg, settings = fit200$settings),
               error = function(e) conditionMessage(e))
cat("refit with stored settings, original data:", if (is.character(r1)) r1 else "ok", "\n")
r2 <- tryCatch(fiticc(fit200$report$obslen, stklen, sel_fun = sf_base, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = "dm")),
               error = function(e) conditionMessage(e))
cat("refit with report obslen, fresh settings:", if (is.character(r2)) r2 else "ok", "\n")
print(all.equal(c(fit200$report$obslen[[1]]@.Data), c(lfd200[[1]]@.Data)))
for (r in c(1000, 22.3)) {
  f <- fiticc(lfd200, stklen, sel_fun = c("dsnormal", "dsnormal"), catch_by_gear = cbg,
              settings = list(CVL = 0.1, obs_model = "dm"), sel_fixed = list(Gillnet = list(rsd = r)), se = FALSE)
  cat("\n-- Gillnet dsnormal rsd", r, "conv", f$opt$convergence, "\n"); nonfin(f)
  cat("opt$par:\n"); print(f$opt$par); cat("non-finite pars:", which(!is.finite(f$opt$par)), "\n")
  pb <- f$obj$env$last.par.best
  cat("last.par.best length", length(pb), "opt length", length(f$opt$par), "same", isTRUE(all.equal(unname(pb), unname(f$opt$par))), "\n")
  print(names(f$obj$env$last.par.best)[!names(f$obj$env$last.par.best) %in% names(f$opt$par)])
  h <- tryCatch(optimHess(f$opt$par, f$obj$fn, f$obj$gr), error = function(e) conditionMessage(e))
  cat("optimHess at opt$par:", if (is.character(h)) h else paste("ok, min eigen", signif(min(eigen(h)$values), 3)), "\n")
  # replicate optimHess and stop at the first non-finite gradient
  dp <- f$opt$par; hit <- FALSE
  for (i in seq_along(dp)) {
    for (h in c(1e-3, -2e-3)) {
      dp[i] <- dp[i] + h; gi <- as.numeric(f$obj$gr(dp))
      if (any(!is.finite(gi))) {
        cat("non-finite gradient at i =", i, names(dp)[i], "step", h, "; fn =", f$obj$fn(dp),
            "; bad components:", which(!is.finite(gi)), "\n")
        rp <- f$obj$report(dp)
        for (nm in names(rp)) if (any(!is.finite(unlist(rp[[nm]])))) cat("  non-finite in report:", nm, "\n")
        print(range(rp$Sel)); hit <- TRUE; break
      }
    }
    if (hit) break
    dp[i] <- dp[i] + 1e-3
  }
  if (!hit) cat("replicated optimHess: no non-finite gradient\n")
  sdr <- withCallingHandlers(tryCatch(TMB::sdreport(f$obj), error = function(e) { cat("sdreport error:", conditionMessage(e), "\n"); NULL }),
                             error = function(e) print(tail(sys.calls(), 6)))
  if (!is.null(sdr)) cat("sdreport ok, pdHess", sdr$pdHess, "\n")
  f2 <- tryCatch(fiticc(lfd200, stklen, sel_fun = c("dsnormal", "dsnormal"), catch_by_gear = cbg,
              settings = list(CVL = 0.1, obs_model = "dm"), sel_fixed = list(Gillnet = list(rsd = r))), error = function(e) conditionMessage(e))
  cat("with se:", if (is.character(f2)) f2 else "ok", "\n")
}
ds <- dome_sens_flicc(fit200, gear = "Gillnet")
print(ds$table); print(ds$errors)
for (nm in names(ds$fits)) { f <- ds$fits[[nm]]; cat(nm, "pdHess", f$rep$pdHess, "maxgrad", signif(f$max_gradient, 2), "\n") }
fit_cut <- fs(tail_cut = 47)
dsc <- dome_sens_flicc(fit_cut, gear = "Gillnet")
print(dsc$table); print(dsc$errors)
# Trawl dsnormal with rsd fixed, too
dst <- dome_sens_flicc(fit200, gear = "Trawl", se = FALSE)
print(dst$table)

# priors
pr <- function(f) { p <- f$par[f$par$par %in% c("log_Linf", "log_Mk", "log_Galpha"), ]
  p$value <- ifelse(p$par == "log_Galpha", 1 / sqrt(exp(p$mpd)), exp(p$mpd)); p }
for (s in list(list(linf.sd = 0.05), list(Mk.sd = 0.1), list(CVL.sd = 0.25),
               list(linf.sd = 0.05, Mk.sd = 0.1, CVL.sd = 0.25),
               list(linf.sd = 0.05, tail_cut = 47), list(linf.sd = 0.05, Mk.sd = 0.1, CVL.sd = 0.25, tail_cut = 47),
               list(linf.sd = 0.05, maxsd = 3, ngtg = 19), list(linf.sd = 0.05, rob_eps = 1e-4))) {
  f <- do.call(fs, s)
  cat("\n==", paste(names(s), unlist(s), sep = "=", collapse = ", "), " conv", f$opt$convergence,
      "pdHess", f$rep$pdHess, "maxgrad", signif(f$max_gradient, 3), "SPR", round(tail(c(f$report$spr), 1), 3), "\n")
  print(pr(f))
}
# profile of the objective in Linf (all else re-estimated)
prof <- sapply(c(54, 54.5, 55, 55.5, 56, 57), function(li) {
  lh2 <- lhpar; lh2["linf"] <- li; st <- stocklen(lfd, lh2, m_model = "constant")
  f <- fiticc(lfd200, st, sel_fun = sel_fun, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = "dm"), se = FALSE)
  c(Linf = li, nll = f$opt$objective)
})
print(round(t(prof), 2))
