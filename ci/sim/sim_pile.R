# Simulation: pile-up correction (LBB Pile 0/1), selectivity shape and tail cut
# ---------------------------------------------------------------------------
# Truth: independent annual age-structured OM for alfonsino (ci/sim/om_sim.R,
# from the dev-gtgdyn OM work): length-at-age CV 0.1, recruitment as an annual
# pulse (sampled as a mid-year snapshot of all cohorts) or spread over the year,
# one gillnet fishery (logistic, or dome fitted as logistic), equilibrium at a
# given F/Fmsy, two-stage length sampling (100 fish for ages, 500 lengths).
# Fits per data set (3 years, each fitted separately as LBB):
#   LBB-style (one growth type, Linf and M/K at true values with LBB's narrow priors)
#     P0     Pile = 0 (no pile-up correction, LBB default)
#     P1     Pile = 1 (numbers integrated over the class)
#     P0cut  Pile = 0, tail cut at 0.85 Linf
#     P1cut  Pile = 1, tail cut at 0.85 Linf
#   GTG    FLicc default (13 growth-type groups, CVL 0.1, estimated DM), years jointly
# Metrics: relative error log(est/true) of apical F/M, EPR (= LBB B/B0) and SPR.
library(FLicc)
here <- if (file.exists("new/ci/sim/om_sim.R")) "new/ci/sim" else "ci/sim"
source(file.path(here, "om_sim.R"))
od <- if (dir.exists("out/tests")) "out/tests" else "."
set.seed(20261019)
nit <- as.integer(Sys.getenv("SIM_NIT", "30"))
ncores <- max(1L, parallel::detectCores() - 1L)
cat("iterations", nit, "cores", ncores, "\n")

lh <- om_lh
MK_true <- lh$M / lh$k
Lm50 <- lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0)))
lhpar <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = Lm50, a = lh$a, b = lh$b)

make_st <- function(rec, sel) {
  st <- om_structure(persist = FALSE, continuous = rec == "continuous", gill = if (sel == "logistic") "logistic" else "dome",
                     gill_sdr = 10)
  st$Fsel[, "Trawl"] <- 0                       # one gear: gillnet only
  st
}
truth <- function(st, Fbar) {
  N <- om_eq(st, Fbar); N0 <- om_eq(st, 0)
  vul <- st$Fsel[, "Gillnet"]                   # proportional to gear selectivity of each cell
  c(FM = Fbar * st$Fap[["Gillnet"]] / st$lh$M,
    EPR = sum(N * st$cells$wt * vul) / sum(N0 * st$cells$wt * vul),
    SPR = om_spr(st, Fbar))
}
sample_lfd <- function(st, Fbar, ny = 3, n_age = 100, n_len = 500) {
  N <- om_eq(st, Fbar); Fc <- Fbar * st$Fsel[, "Gillnet"]; Z <- st$lh$M + Fbar * rowSums(st$Fsel)
  p <- N * Fc / Z * (1 - exp(-Z)); p <- p / sum(p)
  m <- sapply(seq_len(ny), function(y) {
    na <- as.numeric(rmultinom(1, n_age, p)) / n_age
    pl <- as.numeric(na %*% st$cond$Gillnet)
    as.numeric(rmultinom(1, n_len, pl / sum(pl)))
  })
  FLQuants(Gillnet = FLQuant(m, dimnames = list(len = st$lb, year = 2001:(2000 + ny)), units = "cm"))
}
fit_all <- function(lfd) {
  est <- function(tab) c(FM = stats::median(tab$FM), EPR = stats::median(tab$EPR), SPR = stats::median(tab$SPR))
  lbb <- function(pile, cut) {
    set <- if (cut) list(tail_cut = 0.85 * lh$linf) else list()
    r <- suppressMessages(lbb_flicc(lfd, Linf = lh$linf, MK = MK_true, k = lh$k, Lm50 = Lm50, a = lh$a, b = lh$b,
                                    pile = pile, settings = set, verbose = FALSE, se = FALSE))
    est(r$table)
  }
  out <- list(P0 = lbb(0, FALSE), P1 = lbb(1, FALSE), P0cut = lbb(0, TRUE), P1cut = lbb(1, TRUE))
  stk <- stocklen(lfd, lhpar, m_model = "constant")
  g <- suppressMessages(fiticc(lfd, stk, sel_fun = "logistic", catch_by_gear = 1, se = FALSE, verbose = FALSE,
                               settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating")))
  out$GTG <- c(FM = stats::median(as.numeric(g$report$Fap) / lh$M),
               EPR = stats::median(as.numeric(g$report$epr)), SPR = stats::median(as.numeric(g$report$spr)))
  out
}

cells <- expand.grid(rec = c("pulse", "continuous"), sel = c("logistic", "dome"), FFmsy = c(0.5, 1, 2),
                     stringsAsFactors = FALSE)
res <- list()
for (ci in seq_len(nrow(cells))) {
  cc <- cells[ci, ]
  st <- make_st(cc$rec, cc$sel)
  Fbar <- cc$FFmsy * om_fmsy(st)
  tr <- truth(st, Fbar)
  t0 <- Sys.time()
  one <- function(it) {
    lfd <- sample_lfd(st, Fbar)
    f <- tryCatch(fit_all(lfd), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    do.call(rbind, lapply(names(f), function(v) data.frame(rec = cc$rec, sel = cc$sel, FFmsy = cc$FFmsy, it = it,
      fit = v, FM = f[[v]][["FM"]], EPR = f[[v]][["EPR"]], SPR = f[[v]][["SPR"]],
      FM_true = tr[["FM"]], EPR_true = tr[["EPR"]], SPR_true = tr[["SPR"]])))
  }
  r <- parallel::mclapply(seq_len(nit), one, mc.cores = ncores)
  res[[ci]] <- do.call(rbind, r)
  cat(sprintf("cell %d/%d %s %s F/Fmsy %.1f: %d ok, %.0f s\n", ci, nrow(cells), cc$rec, cc$sel, cc$FFmsy,
              length(unique(res[[ci]]$it)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- do.call(rbind, res)
res$eFM <- log(res$FM / res$FM_true); res$eEPR <- log(res$EPR / res$EPR_true); res$eSPR <- log(res$SPR / res$SPR_true)
utils::write.csv(res, file.path(od, "sim_pile.csv"), row.names = FALSE)
q <- function(v) sprintf("%+.2f [%+.2f,%+.2f]", stats::median(v, na.rm = TRUE),
                         stats::quantile(v, 0.25, na.rm = TRUE), stats::quantile(v, 0.75, na.rm = TRUE))
sm <- aggregate(cbind(eFM, eEPR, eSPR) ~ rec + sel + FFmsy + fit, res, q)
sm <- sm[order(sm$sel, sm$rec, sm$FFmsy, match(sm$fit, c("P0", "P1", "P0cut", "P1cut", "GTG"))), ]
options(width = 200)
cat("\nMedian [IQR] of log(estimate / truth)\n")
print(sm, row.names = FALSE)
tr <- unique(res[, c("rec", "sel", "FFmsy", "FM_true", "EPR_true", "SPR_true")]); print(tr, row.names = FALSE)
library(ggplot2)
d <- rbind(data.frame(res[, c("rec", "sel", "FFmsy", "fit")], what = "F/M", err = res$eFM),
           data.frame(res[, c("rec", "sel", "FFmsy", "fit")], what = "EPR (B/B0)", err = res$eEPR))
d$fit <- factor(d$fit, levels = c("P0", "P1", "P0cut", "P1cut", "GTG"))
p <- ggplot(d, aes(fit, err, fill = fit)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.4) + facet_grid(what ~ sel + rec + FFmsy, labeller = label_both) +
  coord_cartesian(ylim = c(-2, 2)) + labs(x = NULL, y = "log(estimate / truth)") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
                     strip.text = element_text(size = 6))
ggsave(file.path(od, "sim_pile.png"), p, width = 14, height = 6, dpi = 100)
