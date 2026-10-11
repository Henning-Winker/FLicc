# Catch in the likelihood (sourced by run_trial.R when mode.txt = catch).
# gtg.dyn with penalised R deviations (sigmaR 0.5) and catch weight for all model
# years (1995-2024; length data 2005-2024). Catch scenarios:
#   c_exact : true catch, CV 0.05 in the likelihood
#   c_cv10  : lognormal observation error CV 0.1, CV 0.1 in the likelihood
#   c_under : 20% constant under-reporting + CV 0.1 error
# plus R deviations as random effects (sigmaR estimated) with c_cv10.
cfgc <- make_cfg(ngtg = 19, maxsd = 3, nq = 1, db = 0.5, h = 4)   # 'fast' numerics, base GTG set
rows <- list()
for (scen in c("noRdev", "Rdev")) for (seed in 1:4) {
  obs <<- as.matrix(read.csv(file.path(dd, sprintf("obs_%s_%d.csv", scen, seed)), header = FALSE))
  tru <- read.csv(file.path(dd, sprintf("truth_%s_%d.csv", scen, seed)))
  ct  <- read.csv(file.path(dd, sprintf("catch_%s_%d.csv", scen, seed)))
  set.seed(100 * seed + nchar(scen))
  eps <- rnorm(nrow(ct), 0, 0.1) - 0.1^2 / 2
  cobs <- list(c_exact = ct$catch_w, c_cv10 = ct$catch_w * exp(eps), c_under = 0.8 * ct$catch_w * exp(eps))
  ccv  <- c(c_exact = 0.05, c_cv10 = 0.1, c_under = 0.1)
  feq  <- fit1(obs, list(pop = 0, cfg = cfgc))
  f0   <- fit1(obs, list(pop = 1, rmode = 1, sigR = 0.5, cfg = cfgc), start = mk_start(feq$pl))
  fits <- list(eq = feq, no_catch = f0)
  for (cn in names(cobs))
    fits[[cn]] <- fit1(obs, list(pop = 1, rmode = 1, sigR = 0.5, cfg = cfgc, cobs = cobs[[cn]], catch_cv = ccv[[cn]]),
                       start = mk_start(f0$pl))
  fits$c_cv10_RE <- fit1(obs, list(pop = 1, rmode = 1, sigR = 0.5, est_sigR = TRUE, re = "logR", cfg = cfgc,
                                   cobs = cobs$c_cv10, catch_cv = 0.1), start = mk_start(fits$c_cv10$pl), limit = 600)
  kd <- 10 + seq_len(nrow(tru))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    say("%-6s seed %d %-10s conv %-5s %5.1f s  R0 %s  sigR %s  F final est/true %s / %.3f", scen, seed, nm, f$conv, f$secs,
        if (is.null(f$R0)) "-" else sprintf("%.2f", f$R0), if (is.numeric(f$sig)) sprintf("%.2f", f$sig[2]) else "-",
        if (is.null(f$F)) "NA" else sprintf("%.3f", tail(f$F, 1)), tail(tru$F, 1))
    if (!is.null(f$spr)) {
      Fd <- if (nm == "eq") f$F else f$F[kd]
      rows[[length(rows) + 1]] <- data.frame(scen = scen, seed = seed, model = nm, year = tru$year, spr = f$spr,
                                             spr_true = tru$spr, F = Fd, F_true = tru$F, secs = f$secs, conv = f$conv)
    }
  }
  write.csv(do.call(rbind, rows), file.path(out, "catch_fits.csv"), row.names = FALSE)
}
d <- do.call(rbind, rows)
d$e <- (d$spr - d$spr_true) / d$spr_true; d$eF <- (d$F - d$F_true) / d$F_true
per <- function(x, v, a, b) { s <- x$year >= a & x$year <= b; sprintf("%+.0f%% (%.0f%%)", 100 * median(x[[v]][s]), 100 * sqrt(mean(x[[v]][s]^2))) }
mods <- c("eq", "no_catch", "c_exact", "c_cv10", "c_under", "c_cv10_RE")
tab <- do.call(rbind, lapply(split(d, list(d$scen, d$model), drop = TRUE), function(x) data.frame(
  scen = x$scen[1], model = x$model[1], secs = median(unique(x[, c("seed", "secs")])$secs),
  rebuild = per(x, "e", 2011, 2018), recent = per(x, "e", 2019, 2023), final = per(x, "e", 2022, 2024),
  F_final = per(x, "eF", 2022, 2024))))
tab <- tab[order(tab$scen, match(tab$model, mods)), ]
md <- c("# gtg.dyn with catch in the likelihood", "",
        "gtg.dyn (19 GTG, fast numerics), penalised R deviations (sigmaR 0.5) unless RE; catch weight 1995-2024, LFD 2005-2024.",
        "Relative error: median (RMSE) over 4 datasets per scenario. 'final' = last 3 years 2022-2024.", "",
        "| scenario | model | median s | SPR rebuild | SPR recent | SPR final 3 yr | F final 3 yr |", "|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %.1f | %s | %s | %s | %s |", tab$scen, tab$model, tab$secs, tab$rebuild, tab$recent, tab$final, tab$F_final))
writeLines(md, file.path(out, "summary.md")); cat(md, sep = "\n")
