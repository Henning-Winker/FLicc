# gtg.dyn evaluated at the simulation's true parameters (F by gear, selectivity,
# no recruitment deviations, Beverton-Holt s = 0.75), R0 free: does the
# predicted catch follow the observed catch, and how does the objective compare
# with the fitted one? Usage: Rscript check_truth.R <dir> <out>
args <- commandArgs(TRUE); here <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source(file.path(here, "om_sim.R"))
sink(file.path(out, "check_truth.log"), split = TRUE)
lh <- om_lh
lhpar <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))), a = lh$a, b = lh$b)
for (sn in c("base", "om")) {
  sc <- list(base = list(TRUE, TRUE, 0), om = list(FALSE, FALSE, 0.5))[[sn]]
  st <- om_structure(persist = sc[[1]], continuous = sc[[2]])
  set.seed(1001); sim <- om_simulate(st, sigR = sc[[3]])
  lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
  stk <- stocklen(lfd, lhpar, m_model = "constant")
  cq <- FLQuants(lapply(colnames(sim$catch_w), function(g) FLQuant(sim$catch_w[, g], dimnames = list(year = sim$years))))
  names(cq) <- colnames(sim$catch_w)
  set <- list(pop_model = "gtg.dyn", obs_model = "dm", CVL = 0.1, maxsd = 3, dyn_rec = "bh", dyn_s = 0.75,
              sigmaR = 0.5, R_re = FALSE, catch_cv = 0.1, dyn_npre = "catch")
  f <- fiticc(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cq, settings = set, se = FALSE)
  td <- f$tmb_data; npre <- td$dyn_npre
  pl <- f$obj$env$parList(f$opt$par)
  # true parameters
  Ftrue <- outer(sim$Fbar, st$Fap)                       # years x gears, annual apical F
  pt <- pl
  pt$log_Fk[]   <- log(Ftrue[(npre + 1):(npre + td$nyear), ] / lh$k)
  pt$log_Fpre[] <- log(Ftrue[seq_len(npre), , drop = FALSE] / lh$k)
  pt$Sm <- c(log(23.35 / lh$linf), log(3.19 / lh$linf), log(11.22 / lh$linf),
             log(29.7 / lh$linf), log((34.1 - 29.7) / lh$linf))
  pt$log_Rdev[] <- 0
  dd <- td[setdiff(names(td), c("gear_names", "year_names", "Linf_init", "Mk_init", "CVL_init", "Galpha_init"))]
  map <- lapply(pt, function(x) factor(rep(NA, length(x)))); map$log_R0 <- NULL
  o <- TMB::MakeADFun(dd, pt, map = map, DLL = "FLicc", silent = TRUE)
  ot <- nlminb(o$par, o$fn, o$gr)
  rt <- o$report(); rf <- f$obj$report(f$obj$env$last.par.best)
  cat("\n==", sn, ": objective at truth (R0 fitted)", round(ot$objective, 2), " at estimate", round(f$opt$objective, 2), "\n")
  yrs <- sim$years[(length(sim$years) - nrow(rt$catch_pred) + 1):length(sim$years)]
  for (g in 1:2) {
    cat("gear", colnames(sim$catch_w)[g], ": observed / predicted catch at truth (every 4th year)\n")
    ii <- seq(1, length(yrs), by = 4)
    print(round(rbind(year = yrs[ii], ratio_truth = (sim$catch_w[match(yrs, sim$years), g] / rt$catch_pred[, g])[ii],
                      ratio_fit = (sim$catch_w[match(yrs, sim$years), g] / rf$catch_pred[, g])[ii]), 3))
  }
  cat("SPR truth / at true pars / fitted (data years):\n")
  print(round(rbind(truth = sim$spr[match(sim$data_years, sim$years)], at_true = rt$spr_y, fitted = rf$spr_y), 3))
  cat("SSB/SSB0: sim / model at true pars (data years):\n")
  print(round(rbind(sim = sim$ssb[match(sim$data_years, sim$years)] / om_sbpr(st, 0), model = rt$dep_y), 3))
  cat("estimated vs true apical F, last 5 years:\n")
  print(round(rbind(true = Ftrue[nrow(Ftrue) - 4:0, 1], est = exp(tail(pl$log_Fk[, 1], 5)) * lh$k), 3))
}
sink()
