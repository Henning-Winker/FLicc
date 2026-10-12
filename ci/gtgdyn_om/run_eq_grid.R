# Equilibrium FLicc (gtg) on OM-like simulations: lag of SPR under changing F
# (tail cut x random walk on F) and sensitivity to dome-shaped selectivity,
# input Linf and CVL (x tail cut). Tail cuts are absolute lengths (fractions of
# the true Linf), so all fits of a dataset use the same bins.
# Usage: Rscript run_eq_grid.R <dir> <out> [n iterations]
args <- commandArgs(TRUE)
here <- args[1]; out <- args[2]; nit <- if (length(args) > 2) as.integer(args[3]) else 10
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source(file.path(here, "om_sim.R"))
log_con <- file(file.path(out, "grid.log"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = log_con); flush(log_con) }
lh <- om_lh
cut_frac <- c(none = NA, c0.90 = 0.90, c0.85 = 0.85, c0.80 = 0.80, c0.75 = 0.75)
cut_len  <- round(cut_frac * lh$linf, 1)
fmods <- list(sF_est = list(), sF0.3 = list(sigmaF = 0.3), sF0.5 = list(sigmaF = 0.5),
              sF1 = list(sigmaF = 1), noRW = list(prior_sigmaF = FALSE))
base <- list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3)

fit1 <- function(lfd, stk, cbg, set) {
  t0 <- Sys.time()
  f <- tryCatch(suppressWarnings(fiticc(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                                        settings = set, se = FALSE)), error = function(e) NULL)
  list(spr = if (is.null(f)) NA else c(f$report$spr), conv = if (is.null(f)) 9 else f$opt$convergence,
       secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
mk_data <- function(sim, linf_in) {
  lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
  lhp <- FLPar(linf = linf_in, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))),
               a = lh$a, b = lh$b)
  kd <- match(sim$data_years, sim$years)
  cbg <- colSums(sim$catch_w[kd, , drop = FALSE]); cbg <- cbg / sum(cbg)
  list(lfd = lfd, stk = stocklen(lfd, lhp, m_model = "constant"), cbg = cbg, kd = kd)
}

rows <- list()
add <- function(block, scen, it, sim, d, model, cut, r, extra = NULL)
  rows[[length(rows) + 1]] <<- data.frame(block = block, scen = scen, iter = it, model = model, cut = cut,
                                         year = sim$data_years, spr = r$spr, spr_true = sim$spr[d$kd],
                                         conv = r$conv, secs = r$secs)

## ---- block A: lag under changing F ------------------------------------------
for (shape in c("rebuild", "decline", "flat")) {
  st <- om_structure()
  for (it in seq_len(nit)) {
    set.seed(2000 + it)
    sim <- om_simulate(st, sigR = 0.5, fshape = shape)
    d <- mk_data(sim, lh$linf)
    for (cn in names(cut_frac)) for (fm in names(fmods)) {
      set <- modifyList(base, fmods[[fm]])
      if (!is.na(cut_len[cn])) set$tail_cut <- cut_len[[cn]]
      add("lag", shape, it, sim, d, fm, cn, fit1(d$lfd, d$stk, d$cbg, set))
    }
    say("lag %-8s it %2d done (%d fits)", shape, it, length(rows))
    write.csv(do.call(rbind, rows), file.path(out, "grid_fits.csv"), row.names = FALSE)
  }
}

## ---- block B: dome, Linf and CVL sensitivity x tail cut ---------------------
for (shape in c("flat", "rebuild")) for (gs in c("logistic", "dome20", "dome10")) {
  st <- om_structure(gill = if (gs == "logistic") "logistic" else "dome",
                     gill_sdr = if (gs == "dome10") 10 else 20)
  for (it in seq_len(nit)) {
    set.seed(3000 + it)
    sim <- om_simulate(st, sigR = 0.5, fshape = shape)
    for (lm in c(0.9, 1, 1.1)) {
      d <- mk_data(sim, lh$linf * lm)
      for (cv in c(0.07, 0.1, 0.13)) for (cn in names(cut_frac)[1:4]) {
        set <- modifyList(base, list(sigmaF = 0.3, CVL = cv))
        if (!is.na(cut_len[cn])) set$tail_cut <- cut_len[[cn]]
        add("sens", paste(shape, gs, sep = "_"), it, sim, d, sprintf("Linf%.1f_CVL%.2f", lm, cv), cn,
            fit1(d$lfd, d$stk, d$cbg, set))
      }
    }
    say("sens %-8s %-8s it %2d done (%d fits)", shape, gs, it, length(rows))
    write.csv(do.call(rbind, rows), file.path(out, "grid_fits.csv"), row.names = FALSE)
  }
}

## ---- summaries ----------------------------------------------------------------
all <- do.call(rbind, lapply(rows, function(x) x[, c("block", "scen", "iter", "model", "cut", "year", "spr", "spr_true", "conv", "secs")]))
all$e <- (all$spr - all$spr_true) / all$spr_true
yr <- sort(unique(all$year))
pc <- function(x) sprintf("%+.0f%%", 100 * x)
lagtab <- do.call(rbind, lapply(split(all[all$block == "lag", ], list(all$scen[all$block == "lag"], all$model[all$block == "lag"],
                                                                     all$cut[all$block == "lag"]), drop = TRUE), function(x) {
  ok <- x[x$conv == 0 & is.finite(x$spr), ]
  per_it <- lapply(split(ok, ok$iter), function(z) {
    se <- coef(lm(spr ~ year, z))[2]; stt <- coef(lm(spr_true ~ year, z))[2]
    c(slope_ratio = se / stt, noise = sd(diff(z$spr - z$spr_true)))
  })
  pm <- do.call(rbind, per_it)
  data.frame(scen = x$scen[1], model = x$model[1], cut = x$cut[1], conv = sprintf("%d/%d", length(unique(ok$iter)), length(unique(x$iter))),
             bias = pc(median(ok$e)), rmse = pc(sqrt(mean(ok$e^2))), bias_last3 = pc(median(ok$e[ok$year >= max(yr) - 2])),
             slope_ratio = round(median(pm[, 1]), 2), noise = round(median(pm[, 2]), 3))
}))
lagtab <- lagtab[order(lagtab$scen, match(lagtab$model, names(fmods)), match(lagtab$cut, names(cut_frac))), ]
senstab <- do.call(rbind, lapply(split(all[all$block == "sens", ], list(all$scen[all$block == "sens"], all$model[all$block == "sens"],
                                                                       all$cut[all$block == "sens"]), drop = TRUE), function(x) {
  ok <- x[x$conv == 0 & is.finite(x$spr), ]
  data.frame(scen = x$scen[1], assumed = x$model[1], cut = x$cut[1], conv = sprintf("%d/%d", length(unique(ok$iter)), length(unique(x$iter))),
             bias = pc(median(ok$e)), rmse = pc(sqrt(mean(ok$e^2))))
}))
senstab <- senstab[order(senstab$scen, senstab$assumed, match(senstab$cut, names(cut_frac))), ]
# compact sensitivity view: median bias by assumed Linf / CVL, per scenario and cut
sens_wide <- reshape(senstab[, c("scen", "assumed", "cut", "bias")], idvar = c("scen", "assumed"), timevar = "cut", direction = "wide")
write.csv(lagtab, file.path(out, "lag_table.csv"), row.names = FALSE)
write.csv(senstab, file.path(out, "sens_table.csv"), row.names = FALSE)
mdt <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", strrep("---|", ncol(df))),
                      apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
md <- c("# Equilibrium FLicc on OM-like simulations: lag and sensitivity", "",
        sprintf("Alfonsino-like OM (fixed-CV length-at-age, annual pulse recruitment, BH h = 0.75, sigmaR 0.5), two gears,
two-stage sampling (n_age 50/30, n_len 150/100), length data 2012-2024, %d iterations. Tail cuts are absolute lengths
(%s cm = 0.90/0.85/0.80/0.75 x true Linf). Relative error in SPR at the year's F.", nit, paste(cut_len[-1], collapse = "/")), "",
        "## A. Lag under changing F (correct biology, logistic Gillnet)", "",
        "slope_ratio = estimated / true linear trend in SPR over the data years (1 = no lag; rebuild and decline only).",
        "noise = SD of the year-to-year change in the error.", "", mdt(lagtab), "",
        "## B. Sensitivity: Gillnet dome (fitted as logistic), input Linf and CVL, x tail cut (sigmaF 0.3)", "",
        "Median bias by assumed Linf multiplier and CVL (rows) and tail cut (columns).", "", mdt(sens_wide))
writeLines(md, file.path(out, "summary.md")); close(log_con)
