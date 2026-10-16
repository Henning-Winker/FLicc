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

ds <- dome_sens_flicc(fit200, gear = "Gillnet")
print(ds$table); print(ds$errors)
for (nm in names(ds$fits)) { f <- ds$fits[[nm]]; cat(nm, "pdHess", f$rep$pdHess, "maxgrad", signif(f$max_gradient, 2), "\n") }
stopifnot(nrow(ds$table) == 5, length(ds$errors) == 0)
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
  st <- stklen; st@lhpar["linf"] <- li
  f <- fiticc(lfd200, st, sel_fun = sel_fun, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = "dm"), se = FALSE)
  c(Linf = li, nll = f$opt$objective)
})
print(round(t(prof), 2))
