# Speed / accuracy trade-off for gtg.dyn (sourced by run_trial.R when mode.txt = speed).
# Model: gtg.dyn, dome selectivity, penalised recruitment deviations (sigmaR fixed 0.5),
# penalised RW on log F (sigmaF 0.3) -- the MSE candidate. Each setting is started from
# the equilibrium fit with the same GTG set.
settings <- list(
  base          = list(),
  gtg13_sd2     = list(ngtg = 13, maxsd = 2),
  nq1           = list(nq = 1),
  db0.5         = list(db = 0.5),
  h4            = list(h = 4),
  npre5         = list(n_pre = 5),
  fast          = list(ngtg = 13, maxsd = 2, nq = 1, db = 0.5, h = 4),
  fastest       = list(ngtg = 13, maxsd = 2, nq = 1, db = 0.5, h = 4, n_pre = 5)
)
cfgs <- lapply(settings, function(a) do.call(make_cfg, a))
for (nm in names(cfgs)) say("setting %-10s GTG %2d  nodes %5d  tq %s  db %.2f  n_pre %d", nm, nrow(cfgs[[nm]]$gtg),
                            length(cfgs[[nm]]$grid$node_g), paste(cfgs[[nm]]$tq, collapse = ","), cfgs[[nm]]$db, cfgs[[nm]]$n_pre)
rows <- list()
for (scen in c("noRdev", "Rdev")) for (seed in 1:4) {
  obs <<- as.matrix(read.csv(file.path(dd, sprintf("obs_%s_%d.csv", scen, seed)), header = FALSE))
  tru <- read.csv(file.path(dd, sprintf("truth_%s_%d.csv", scen, seed)))
  for (nm in names(cfgs)) {
    cf <- cfgs[[nm]]
    feq <- fit1(obs, list(pop = 0, cfg = cf))
    f <- fit1(obs, list(pop = 1, rmode = 1, sigR = 0.5, cfg = cf), start = mk_start(feq$pl))
    say("%-6s seed %d %-10s conv %-5s eq %4.1f s  dyn %5.1f s  it %s", scen, seed, nm, f$conv, feq$secs, f$secs, format(f$iter %||% NA))
    if (!is.null(f$spr))
      rows[[length(rows) + 1]] <- data.frame(scen = scen, seed = seed, setting = nm, year = tru$year, spr = f$spr,
                                             spr_true = tru$spr, secs = f$secs, conv = f$conv)
  }
  write.csv(do.call(rbind, rows), file.path(out, "speed_fits.csv"), row.names = FALSE)
}
d <- do.call(rbind, rows)
base <- d[d$setting == "base", c("scen", "seed", "year", "spr")]; names(base)[4] <- "spr_base"
d <- merge(d, base)
d$e <- (d$spr - d$spr_true) / d$spr_true
per <- function(x, a, b) { s <- x$year >= a & x$year <= b; sprintf("%+.0f%% (%.0f%%)", 100 * median(x$e[s]), 100 * sqrt(mean(x$e[s]^2))) }
tab <- do.call(rbind, lapply(split(d, d$setting), function(x) data.frame(
  setting = x$setting[1], secs = median(unique(x[, c("scen", "seed", "secs")])$secs),
  conv = sprintf("%d/%d", sum(unique(x[, c("scen", "seed", "conv")])$conv), nrow(unique(x[, c("scen", "seed")]))),
  maxdiff_base = max(abs(x$spr - x$spr_base)), rebuild = per(x, 2011, 2018), recent = per(x, 2019, 2023), final = per(x, 2024, 2024))))
tab <- tab[match(names(settings), tab$setting), ]
md <- c("# gtg.dyn speed settings", "",
        "Model: gtg.dyn, dome selectivity, penalised R deviations (sigmaR 0.5 fixed), penalised RW log F.",
        "8 datasets (both scenarios). SPR error vs truth: median (RMSE).", "",
        "| setting | median s | converged | max abs SPR diff vs base | rebuild | recent | final |", "|---|---|---|---|---|---|---|",
        sprintf("| %s | %.1f | %s | %.3f | %s | %s | %s |", tab$setting, tab$secs, tab$conv, tab$maxdiff_base, tab$rebuild, tab$recent, tab$final))
writeLines(md, file.path(out, "summary.md")); cat(md, sep = "\n")
