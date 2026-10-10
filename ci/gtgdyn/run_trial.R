# gtg.dyn TMB trial
#   1. reference test: TMB predictions vs the Python prototype for identical inputs
#   2. fits to the 8 simulated rebuild datasets (same data as the Python trial),
#      penalised (as in Python) and with random effects (Laplace), with timings
# Usage: Rscript run_trial.R <dir of this script> <output dir>
args <- commandArgs(TRUE)
here <- if (length(args) > 0) args[1] else "ci/gtgdyn"
out  <- if (length(args) > 1) args[2] else "out_gtgdyn"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(TMB))
`%||%` <- function(a, b) if (is.null(a)) b else a
dd <- file.path(here, "data")
log_con <- file(file.path(out, "trial.log"), open = "wt")
say <- function(...) { msg <- sprintf(...); cat(msg, "\n"); cat(msg, "\n", file = log_con); flush(log_con) }

# ---- compile -----------------------------------------------------------------
cpp <- file.path(out, "gtgdyn.cpp")
file.copy(file.path(here, "gtgdyn.cpp"), cpp, overwrite = TRUE)
t0 <- Sys.time()
fw <- Sys.getenv("TMB_FRAMEWORK", "TMBad")
TMB::compile(cpp, flags = "-O2", framework = fw)
dyn.load(TMB::dynlib(file.path(out, "gtgdyn")))
say("framework %s, compiled in %.0f s", fw, as.numeric(difftime(Sys.time(), t0, units = "secs")))

# ---- biology / grid (exported from the Python prototype) ----------------------
bio <- read.csv(file.path(dd, "bio.csv")); LB <- read.csv(file.path(dd, "lb.csv"))$lb
gtg <- read.csv(file.path(dd, "gtg.csv"))
K <- 0.08; M <- 0.162; nb <- nrow(bio); ng <- nrow(gtg)

make_grid <- function(Tcap = 40, h = 2) {
  A <- dT <- NULL
  A <- t(sapply(gtg$linf, function(li) { r <- 1 - LB / li; ifelse(r > 0, -log(pmax(r, 1e-300)) / K, Inf) }))
  reach <- matrix(0L, ng, nb); dT <- matrix(0, ng, nb)
  for (g in 1:ng) for (j in 1:nb) {
    if (is.finite(A[g, j])) {
      if (is.finite(A[g, j + 1])) { reach[g, j] <- 1L; dT[g, j] <- A[g, j + 1] - A[g, j] }
      else { reach[g, j] <- 2L; dT[g, j] <- Tcap }
    }
  }
  x <- 0.5 + c(-1, 1) * 0.5 / sqrt(3); w <- c(0.5, 0.5)       # 2-point Gauss-Legendre on [0,1]
  ng_ <- nj_ <- tau_ <- w_ <- NULL
  for (g in 1:ng) for (j in 1:nb) if (reach[g, j] > 0) {
    D <- dT[g, j]; np <- ceiling(D / h); hp <- D / np
    tau <- as.vector(outer(x, 0:(np - 1), "+")) * hp
    ng_ <- c(ng_, rep(g - 1L, length(tau))); nj_ <- c(nj_, rep(j - 1L, length(tau)))
    tau_ <- c(tau_, tau); w_ <- c(w_, rep(w, np) * hp)
  }
  Afin <- A; Afin[!is.finite(Afin)] <- 1e10
  list(A = Afin, dT = dT, reach = reach, node_g = as.integer(ng_), node_j = as.integer(nj_),
       node_tau = tau_, node_w = w_)
}
grid <- make_grid()
say("GTG %d, bins %d, quadrature nodes %d", ng, nb, length(grid$node_g))
say("dims: A %s, dT %s, reach %s, lb %d, node_j range %s", paste(dim(grid$A), collapse = "x"), paste(dim(grid$dT), collapse = "x"), paste(dim(grid$reach), collapse = "x"), length(LB), paste(range(grid$node_j), collapse = "-"))

b0 <- -max(grid$A[, 1:nb][grid$reach > 0]) - 1

make_obj <- function(obs, pop, n_pre = 10, db = 0.25, sel_type = 1, rmode = 0, sigF = 0.3, sigR = 0.5,
                     est_sigF = FALSE, est_sigR = FALSE, re = character(0),
                     prior_sigF = c(log(0.3), 0.5, 0), prior_sigR = c(log(0.5), 0.5, 0),
                     par = NULL, map_all = FALSE) {
  ny <- nrow(obs); npre <- if (pop == 1) n_pre else 0; nF <- ny + npre
  mode <- bio$lmid[which.max(colSums(obs))]
  storage.mode(obs) <- "double"
  data <- c(list(pop = as.integer(pop), obs = obs, kd = as.integer(npre + 0:(ny - 1)),
                 tq = c(0.25, 0.75), lmid = bio$lmid, mw = bio$mw, wg = gtg$w, M = M,
                 sel_type = as.integer(sel_type), rmode = as.integer(rmode),
                 prior_sigF = prior_sigF, prior_sigR = prior_sigR,
                 db = db, b0 = b0, nbg = as.integer(ceiling((nF + 1 - b0) / max(db, 0.01)) + 2)), grid)
  if (is.null(par))
    par <- list(logF = rep(log(0.1), nF),
                theta = if (sel_type == 1) c(mode, log(5), log(20)) else c(mode, log(5)),
                logR = rep(0, nF), log_sigF = log(sigF), log_sigR = log(sigR))
  map <- list()
  if (!est_sigF) map$log_sigF <- factor(NA)
  if (!est_sigR || rmode == 0) map$log_sigR <- factor(NA)
  if (rmode == 0) map$logR <- factor(rep(NA, nF))
  if (rmode == 2) map$logR <- factor(c(NA, 1:(nF - 1)))   # level of R not identified by proportions: anchor
  if (map_all) map <- lapply(par[names(par) != "logF"], function(p) factor(rep(NA, length(p))))  # TMB needs >= 1 free parameter
  MakeADFun(data, par, map = map, random = if (length(re)) re else NULL,
            DLL = "gtgdyn", silent = TRUE)
}

# ---- 1. reference test ---------------------------------------------------------
rp <- read.csv(file.path(dd, "ref_par.csv")); th <- read.csv(file.path(dd, "ref_sel.csv"))$theta
obs1 <- as.matrix(read.csv(file.path(dd, "obs_Rdev_1.csv"), header = FALSE))
ref_dyn <- as.matrix(read.csv(file.path(dd, "ref_pred_dyn.csv"), header = FALSE))
ref_eq  <- as.matrix(read.csv(file.path(dd, "ref_pred_eq.csv"), header = FALSE))
par <- list(logF = rp$logF, theta = th, logR = rp$logR, log_sigF = log(0.3), log_sigR = log(0.5))
o <- make_obj(obs1, pop = 1, rmode = 1, par = par, map_all = TRUE, db = 0)
o$fn(o$par); d_dyn <- max(abs(o$report()$pred - ref_dyn))
ref_dyng <- as.matrix(read.csv(file.path(dd, "ref_pred_dyng.csv"), header = FALSE))
o <- make_obj(obs1, pop = 1, rmode = 1, par = par, map_all = TRUE, db = 0.25)
o$fn(o$par); d_dyng <- max(abs(o$report()$pred - ref_dyng))
say("grid vs exact (TMB): %.1e", max(abs(o$report()$pred - ref_dyn)))
k <- 10 + 1:20
par_eq <- list(logF = rp$logF[k], theta = th, logR = rp$logR[k], log_sigF = log(0.3), log_sigR = log(0.5))
o <- make_obj(obs1, pop = 0, rmode = 0, par = par_eq, map_all = TRUE)
o$fn(o$par); d_eq <- max(abs(o$report()$pred - ref_eq))
ref_ok <- d_dyn < 1e-8 && d_eq < 1e-8 && d_dyng < 1e-8
say("REFERENCE TEST: max |TMB - Python|  dyn exact %.2e   dyn grid %.2e   eq %.2e   -> %s", d_dyn, d_dyng, d_eq, if (ref_ok) "PASS" else "FAIL")

# ---- 2. fits -------------------------------------------------------------------
lohi <- function(o, nF, rmode) {
  nm <- names(o$par)
  lo <- ifelse(nm == "logF", log(1e-4), ifelse(nm == "logR", -3, ifelse(nm == "log_sigR", log(0.01), -Inf)))
  hi <- ifelse(nm == "logF", log(3), ifelse(nm == "logR", 3, Inf))
  th <- which(nm == "theta"); lo[th] <- c(10, log(0.5), log(0.5)); hi[th] <- c(60, log(40), log(200))
  list(lo = lo, hi = hi)
}
fit1 <- function(obs, spec, start = NULL, limit = 900) {
  t0 <- Sys.time()
  setTimeLimit(elapsed = limit, transient = TRUE); on.exit(setTimeLimit(elapsed = Inf))
  if (!is.null(start)) spec$par <- start(spec)
  o <- do.call(make_obj, c(list(obs = obs), spec))
  b <- lohi(o)
  opt <- try(nlminb(o$par, o$fn, o$gr, lower = b$lo, upper = b$hi,
                    control = list(eval.max = 2000, iter.max = 1000)), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(opt, "try-error")) return(list(conv = FALSE, secs = secs, msg = as.character(opt)))
  if (length(spec$re)) o$fn(opt$par)
  rep <- o$report(o$env$last.par.best)
  pl <- o$env$parList(par = o$env$last.par.best)
  list(pl = pl, iter = opt$iterations, conv = opt$convergence == 0, nll = opt$objective, secs = secs, spr = rep$spr,
       ssb_rel = rep$ssb_rel, sel = rep$sel, F = rep$F, msg = opt$message,
       sig = exp(o$env$last.par.best[c("log_sigF", "log_sigR")]))
}
specs <- list(
  eq_dome       = list(pop = 0),
  dyn_dome      = list(pop = 1),
  dyn_dome_Rdev = list(pop = 1, rmode = 1, sigR = 0.5),
  dyn_dome_Rrw  = list(pop = 1, rmode = 2, sigR = 0.2),
  dyn_Rdev_RE   = list(pop = 1, rmode = 1, sigR = 0.5, est_sigR = TRUE, re = "logR"),
  dyn_FR_RE     = list(pop = 1, rmode = 1, sigR = 0.5, est_sigR = TRUE, est_sigF = TRUE,
                       prior_sigF = c(log(0.3), 0.5, 1), re = c("logF", "logR"))
)
# starting values: dyn from the equilibrium fit (pre-data years at the first F),
# recruitment variants from dyn_dome, random-effects variants from their penalised twin
start_from <- list(dyn_dome = "eq_dome", dyn_dome_Rdev = "dyn_dome", dyn_dome_Rrw = "dyn_dome",
                   dyn_Rdev_RE = "dyn_dome_Rdev", dyn_FR_RE = "dyn_dome_Rdev")
mk_start <- function(prev) function(spec) {
  lf <- prev$logF; nF_new <- nrow(obs) + if (spec$pop == 1) 10 else 0
  if (length(lf) < nF_new) lf <- c(rep(lf[1], nF_new - length(lf)), lf)
  list(logF = lf, theta = prev$theta, logR = if (length(prev$logR) == nF_new) prev$logR else rep(0, nF_new),
       log_sigF = log(if (is.null(spec$sigF)) 0.3 else spec$sigF), log_sigR = log(if (is.null(spec$sigR)) 0.5 else spec$sigR))
}
res <- list(); tim <- list(); done_all <- list()
run_pass <- function(models) {
  for (scen in c("noRdev", "Rdev")) for (seed in 1:4) {
    obs <<- as.matrix(read.csv(file.path(dd, sprintf("obs_%s_%d.csv", scen, seed)), header = FALSE))
    tru <- read.csv(file.path(dd, sprintf("truth_%s_%d.csv", scen, seed)))
    key <- paste(scen, seed); done <- done_all[[key]] %||% list()
    for (nm in models) {
      prev <- done[[start_from[[nm]] %||% ""]]
      f <- fit1(obs, specs[[nm]], start = if (!is.null(prev$pl)) mk_start(prev$pl),
                limit = if (length(specs[[nm]]$re)) 600 else 300)
      done[[nm]] <- f
      say("%-6s seed %d %-14s conv %-5s %6.1f s  it %s  nll %s  sig %s  %s", scen, seed, nm, f$conv, f$secs, format(f$iter %||% NA),
          if (is.numeric(f$nll)) format(round(f$nll, 1)) else "NA",
          if (is.numeric(f$sig)) paste(round(f$sig, 3), collapse = "/") else "NA", substr(f$msg %||% "", 1, 60))
      tim[[length(tim) + 1]] <<- data.frame(scen = scen, seed = seed, model = nm, secs = f$secs, conv = f$conv)
      if (!is.null(f$spr))
        res[[length(res) + 1]] <<- data.frame(scen = scen, seed = seed, model = nm, year = tru$year,
                                              spr = f$spr, ssb_rel = f$ssb_rel, spr_true = tru$spr, bb0 = tru$bb0)
    }
    done_all[[key]] <<- done
  }
}
write_summary <- function(stage) {
res <- do.call(rbind, res); tim <- do.call(rbind, tim)
write.csv(res, file.path(out, "tmb_fits.csv"), row.names = FALSE)
write.csv(tim, file.path(out, "tmb_times.csv"), row.names = FALSE)

# ---- metrics: TMB vs Python ----------------------------------------------------
py <- read.csv(file.path(dd, "py_fits.csv"))
tr_all <- unique(res[, c("scen", "seed", "year", "spr_true")])
py <- merge(py, tr_all)
periods <- list(decline = c(2005, 2010), rebuild = c(2011, 2018), recent = c(2019, 2023), final = c(2024, 2024))
metr <- function(d, src) {
  do.call(rbind, lapply(split(d, list(d$scen, d$model), drop = TRUE), function(x) {
    e <- (x$spr - x$spr_true) / x$spr_true
    data.frame(src = src, scen = x$scen[1], model = x$model[1],
               t(sapply(periods, function(p) { s <- x$year >= p[1] & x$year <= p[2]
                 sprintf("%+.0f%% (%.0f%%)", 100 * median(e[s]), 100 * sqrt(mean(e[s]^2))) })))
  }))
}
tab <- rbind(metr(res, "TMB"), metr(py, "Python"))
tab <- tab[order(tab$scen, tab$model, tab$src), ]
write.csv(tab, file.path(out, "metrics.csv"), row.names = FALSE)
# agreement of TMB and Python point estimates for identical penalised models
cmp <- merge(res[, c("scen", "seed", "model", "year", "spr")], py[, c("scen", "seed", "model", "year", "spr")],
             by = c("scen", "seed", "model", "year"), suffixes = c("_tmb", "_py"))
agree <- aggregate(cbind(maxdiff = abs(spr_tmb - spr_py)) ~ model + scen, cmp, max)
pysec <- aggregate(secs ~ model, unique(read.csv(file.path(dd, "py_fits.csv"))[, c("scen", "seed", "model", "secs")]), median)
tmbsec <- aggregate(secs ~ model, tim, median)

md <- c("# gtg.dyn TMB trial", "", paste("Stage:", stage), "",
        sprintf("Reference test (TMB vs Python, same inputs): dyn exact %.1e, dyn grid %.1e, eq %.1e -> **%s**", d_dyn, d_dyng, d_eq,
                if (ref_ok) "PASS" else "FAIL"), "",
        "## Median fit time (s)", "", "| model | TMB | Python |", "|---|---|---|",
        sprintf("| %s | %.1f | %s |", tmbsec$model, tmbsec$secs,
                ifelse(is.na(match(tmbsec$model, pysec$model)), "-",
                       sprintf("%.0f", pysec$secs[match(tmbsec$model, pysec$model)]))), "",
        "## Relative error in SPR: median (RMSE)", "",
        "| source | scenario | model | decline | rebuild | recent | final |", "|---|---|---|---|---|---|---|",
        apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "",
        "## Max |SPR TMB - SPR Python| for the same penalised model", "",
        "| model | scenario | max diff |", "|---|---|---|",
        sprintf("| %s | %s | %.3f |", agree$model, agree$scen, agree$maxdiff), "",
        sprintf("Convergence: %d of %d fits", sum(tim$conv), nrow(tim)))
writeLines(md, file.path(out, "summary.md"))
cat(md, sep = "\n")
}
run_pass(c("eq_dome", "dyn_dome", "dyn_dome_Rdev", "dyn_dome_Rrw"))
write_summary("penalised")
run_pass(c("dyn_Rdev_RE", "dyn_FR_RE"))
write_summary("all")
close(log_con)
if (!ref_ok) quit(status = 1)
