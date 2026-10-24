# Residual diagnostics and posterior predictive checks
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
lfd200 <- lfdess(lfd, ess.g = 200)
od <- if (dir.exists("out/tests")) "out/tests" else "."

## 1. real data, dm and mn
for (om in c("dm", "mn")) {
  f <- fiticc(lfd200, stklen, sel_fun = sf, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = om))
  rp <- resid_flicc(f, "pearson"); ro <- resid_flicc(f, "osa")
  cat("\n==", om, " Pearson sd", round(sd(rp$resid), 2), " OSA mean", round(mean(ro$resid, na.rm = TRUE), 2),
      "sd", round(sd(ro$resid, na.rm = TRUE), 2), " n", sum(is.finite(ro$resid)), "\n")
  lm <- lmean_resid_flicc(f); print(lm$runs)
  t0 <- system.time(pc <- ppc_flicc(f, nsim = 300))[["elapsed"]]
  cat("ppc seconds", t0, "\n"); print(pc$coverage)
  ggplot2::ggsave(file.path(od, paste0("resid_bubble_", om, ".png")), plot_resid_bubble(f, cohorts = TRUE),
                  width = 9, height = 5, dpi = 100)
  ggplot2::ggsave(file.path(od, paste0("resid_bubble_osa_", om, ".png")), plot_resid_bubble(ro, fit = f, cohorts = TRUE),
                  width = 9, height = 5, dpi = 100)
  ggplot2::ggsave(file.path(od, paste0("resid_qq_", om, ".png")), plot_resid_qq(ro), width = 7, height = 3.5, dpi = 100)
  ggplot2::ggsave(file.path(od, paste0("lmean_", om, ".png")), lm$plot, width = 8, height = 3.5, dpi = 100)
  ggplot2::ggsave(file.path(od, paste0("ppc_", om, ".png")), pc$plot, width = 8, height = 7, dpi = 100)
  if (om == "dm") fdm <- f
}
pc2 <- ppc_flicc(fdm, nsim = 200, par_draws = TRUE); print(pc2$coverage)

## 2. self-check: data simulated from the fitted model give N(0,1) OSA residuals
selfcheck <- function(f, om) {
  cd <- FLicc:::comp_data_flicc(f)
  set.seed(42)
  sim <- f$tmb_data$obs
  td <- f$tmb_data
  rep0 <- f$obj$report(f$obj$env$last.par.best)
  rpp <- rep0$plen
  for (g in seq_len(dim(sim)[3])) for (y in seq_len(dim(sim)[2])) {
    n <- round(sum(sim[, y, g])); if (n <= 0) next
    p <- rpp[, y, g] / sum(rpp[, y, g])
    a0 <- if (!is.null(rep0$dm_alpha0) && rep0$dm_alpha0[y, g] > 0) rep0$dm_alpha0[y, g] else n   # fitted DM precision
    sim[, y, g] <- if (om == "dm") { a <- rgamma(length(p), shape = a0 * p); rmultinom(1, n, a / sum(a))[, 1] }
                   else rmultinom(1, n, p)[, 1]
  }
  f2 <- f; f2$tmb_data$obs <- sim
  r <- resid_flicc(f2, "osa")$resid
  r <- r[is.finite(r)]
  c(mean = mean(r), sd = sd(r), ks_p = suppressWarnings(ks.test(r, "pnorm")$p.value))
}
sc_dm <- selfcheck(fdm, "dm"); print(round(sc_dm, 3))
stopifnot(abs(sc_dm["mean"]) < 0.1, abs(sc_dm["sd"] - 1) < 0.1)

# residuals by length (single year and several years)
od <- if (dir.exists("out/tests")) "out/tests" else "."
p1 <- plot_resid_len(fdm)
p2 <- plot_resid_len(fdm, years = "each", cap = 3)       # deprecated argument still works
p3 <- plot_resid_len(fdm, type = "osa", pool_years = TRUE)
p4 <- plot_resid_len(fdm, by_gear = FALSE, year = tail(fdm$tmb_data$year_names, 2))
p5 <- plot_resid_len(fdm, gear = 1, year = 1:2)
for (pp in list(p3, p4, p5)) stopifnot(inherits(pp, "ggplot")); invisible(lapply(list(p3, p4, p5), ggplot2::ggplot_build))
ggplot2::ggsave(file.path(od, "resid_len.png"), p1, width = 8, height = 4, dpi = 100)
ggplot2::ggsave(file.path(od, "resid_len_each.png"), p2, width = 8, height = 4, dpi = 100)
print(nrow(attr(p2, "capped")))

# lfd_trim: common length range, counts removed, fit runs
lt <- lfd_trim(lfd, Lmin = 20, Lmax = 60)
rm_tab <- attr(lt, "removed"); print(rm_tab)
L1 <- lapply(lt, function(x) as.numeric(dimnames(x)$len))
stopifnot(all(vapply(L1, function(l) min(l) >= 20 && max(l) < 60, logical(1))),
          all(abs(vapply(names(lfd), function(g) sum(lfd[[g]]) - sum(lt[[g]]), numeric(1)) -
                  (rm_tab$below_Lmin + rm_tab$above_Lmax)) < 1e-6))
ft <- fiticc(lt, stocklen(lt, lhpar, m_model = "constant"), sel_fun = sf, catch_by_gear = cbg,
             settings = list(CVL = 0.1))
print(ft$diagnostics)

## distribution view of the predictive check
pl <- plot_ppc_lfd(fdm, nsim = 100, pool_years = TRUE)
d <- attr(pl, "data")
stopifnot(inherits(pl, "ggplot"), all(c("obs", "lo", "med", "hi", "outside") %in% names(d)),
          all(abs(tapply(d$obs, d$gear, sum) - 1) < 1e-8), all(d$lo <= d$hi))
pc3 <- ppc_flicc(fdm, nsim = 100)
stopifnot(all(pc3$coverage_bins$share_bins_outside_95 < 0.2), inherits(pc3$plot_lfd, "ggplot"))
d1 <- attr(plot_ppc_lfd(pc3, year = 1, gear = "Trawl"), "data")
stopifnot(nlevels(d1$year) == 1, all(d1$gear == "Trawl"))
ggplot2::ggsave(file.path(od, "ppc_lfd_pooled.png"), pl, width = 8, height = 3.5, dpi = 100)
cat("plot_ppc_lfd ok\n")
