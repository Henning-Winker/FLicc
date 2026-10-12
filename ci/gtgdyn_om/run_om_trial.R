# gtg vs gtg.dyn on OM-like simulations, switching OM-vs-model mismatches on in turn.
# Usage: Rscript run_om_trial.R <dir of this script> <output dir> [n iterations]
args <- commandArgs(TRUE)
here <- if (length(args) > 0) args[1] else "ci/gtgdyn_om"
out  <- if (length(args) > 1) args[2] else "out_om"
nit  <- if (length(args) > 2) as.integer(args[3]) else 8
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source(file.path(here, "om_sim.R"))
log_con <- file(file.path(out, "trial.log"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = log_con); flush(log_con) }

# scenarios: from gtg.dyn's own assumptions (base) to the OM (all mismatches on)
scen <- list(
  base       = list(persist = TRUE,  continuous = TRUE,  sigR = 0),
  pulse      = list(persist = TRUE,  continuous = FALSE, sigR = 0),
  nonpersist = list(persist = FALSE, continuous = TRUE,  sigR = 0),
  rdev       = list(persist = TRUE,  continuous = TRUE,  sigR = 0.5),
  om         = list(persist = FALSE, continuous = FALSE, sigR = 0.5)
)
base_set <- list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3)
models <- list(
  gtg        = base_set,
  dyn        = modifyList(base_set, list(pop_model = "gtg.dyn", sigmaR = 0.5, R_re = FALSE)),
  dyn_Rrw    = modifyList(base_set, list(pop_model = "gtg.dyn", dyn_rec = "rw", sigmaR = 0.3, R_re = FALSE)),
  # Beverton-Holt recruitment (steepness 0.75, as the OM) on the dynamic SSB
  dyn_BH     = modifyList(base_set, list(pop_model = "gtg.dyn", dyn_rec = "bh", dyn_s = 0.75,
                                         sigmaR = 0.5, R_re = FALSE)),
  dyn_BH_catch10  = modifyList(base_set, list(pop_model = "gtg.dyn", dyn_rec = "bh", dyn_s = 0.75,
                                              sigmaR = 0.5, R_re = FALSE, catch_cv = 0.1)),
  dyn_BH_catchAll = modifyList(base_set, list(pop_model = "gtg.dyn", dyn_rec = "bh", dyn_s = 0.75,
                                              sigmaR = 0.5, R_re = FALSE, catch_cv = 0.1, dyn_npre = "catch"))
)
lh <- om_lh
lhpar <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))),
               a = lh$a, b = lh$b)

rows <- list(); tim <- list()
for (sn in names(scen)) {
  sc <- scen[[sn]]
  st <- om_structure(persist = sc$persist, continuous = sc$continuous)
  for (it in seq_len(nit)) {
    set.seed(1000 + it)                                   # same F history / noise across scenarios
    sim <- om_simulate(st, sigR = sc$sigR)
    lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
    stk <- stocklen(lfd, lhpar, m_model = "constant")
    kd <- match(sim$data_years, sim$years)
    cbg <- colSums(sim$catch_w[kd, , drop = FALSE]); cbg <- cbg / sum(cbg)
    cobs <- sim$catch_w * exp(matrix(rnorm(length(sim$catch_w), 0, 0.1) - 0.005, nrow(sim$catch_w)))
    cq <- FLQuants(lapply(colnames(cobs), function(g) FLQuant(cobs[, g], dimnames = list(year = sim$years))))
    names(cq) <- colnames(cobs)
    for (mn in names(models)) {
      t0 <- Sys.time()
      cin <- if (is.null(models[[mn]]$catch_cv)) cbg else cq
      f <- tryCatch(fiticc(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cin,
                           settings = models[[mn]], se = FALSE),
                    error = function(e) { say("ERROR %s %d %s: %s", sn, it, mn, conditionMessage(e)); NULL })
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(f)) next
      spr <- c(f$report$spr)
      say("%-10s it %d %-10s conv %d %5.1f s  SPR first/last est %.3f/%.3f true %.3f/%.3f", sn, it, mn,
          f$opt$convergence, secs, spr[1], tail(spr, 1), sim$spr[kd[1]], sim$spr[tail(kd, 1)])
      rows[[length(rows) + 1]] <- data.frame(scen = sn, iter = it, model = mn, year = sim$data_years,
                                             spr = spr, spr_true = sim$spr[kd], secs = secs,
                                             conv = f$opt$convergence)
    }
    write.csv(do.call(rbind, rows), file.path(out, "om_fits.csv"), row.names = FALSE)
  }
  # F history of the last iteration (for the record)
  if (sn == "om") write.csv(data.frame(year = sim$years, Fbar = sim$Fbar, Fmsy = sim$Fmsy, spr = sim$spr,
                                       ssb = sim$ssb, dev = sim$dev), file.path(out, "om_example.csv"),
                            row.names = FALSE)
}

d <- do.call(rbind, rows)
d$e <- (d$spr - d$spr_true) / d$spr_true
yr <- sort(unique(d$year))
per <- list(all = range(yr), first5 = c(yr[1], yr[5]), last3 = c(max(yr) - 2, max(yr)))
fmt <- function(e) sprintf("%+.0f%% (%.0f%%)", 100 * median(e), 100 * sqrt(mean(e^2)))
tab <- do.call(rbind, lapply(split(d, list(d$scen, d$model), drop = TRUE), function(x) {
  # drift: median over iterations of the change in error from the first to the last 3 years
  dr <- sapply(split(x, x$iter), function(z) mean(z$e[z$year >= max(yr) - 2]) - mean(z$e[z$year <= yr[3]]))
  data.frame(scen = x$scen[1], model = x$model[1], n = length(unique(x$iter)),
             secs = round(median(x$secs), 1), conv = sum(tapply(x$conv, x$iter, `[`, 1) == 0),
             all = fmt(x$e[x$year >= per$all[1]]), first5 = fmt(x$e[x$year <= per$first5[2]]),
             last3 = fmt(x$e[x$year >= per$last3[1]]), drift = sprintf("%+.0f%%", 100 * median(dr)))
}))
tab <- tab[order(match(tab$scen, names(scen)), match(tab$model, names(models))), ]
md <- c("# gtg vs gtg.dyn on OM-like simulations", "",
        "Scenarios switch OM-vs-model mismatches on in turn: `base` = gtg.dyn's own assumptions",
        "(persistent growth groups, recruitment spread over the year, no recruitment deviations);",
        "`pulse` = annual pulse recruitment; `nonpersist` = length-at-age with fixed CV 0.1 (no persistent",
        "groups); `rdev` = recruitment deviations (sigma 0.5, rho 0.3); `om` = all mismatches (as the OM).",
        "F history shaped like fudc (peak 3 x Fmsy ~2003, 1.5 x Fmsy from ~2014), length data 2012-2024,",
        "two gears (Trawl dome, Gillnet logistic), two-stage sampling (n_age 50/30, n_len 150/100).", "",
        "Relative error in SPR at the year's F: median (RMSE). drift = change in error from the first",
        "3 to the last 3 data years (median over iterations).", "",
        "| scenario | model | n | median s | converged | all years | first 5 | last 3 | drift |",
        "|---|---|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %d | %.1f | %d | %s | %s | %s | %s |", tab$scen, tab$model, tab$n, tab$secs,
                tab$conv, tab$all, tab$first5, tab$last3, tab$drift))
writeLines(md, file.path(out, "summary.md")); cat(md, sep = "\n")
close(log_con)
