# FLicc-SPM prototype: joint fit of the gtg length model (F fixed effects,
# a nuisance that tunes SPR) and a state-space surplus production model
# (states random) to catch and CPUE, on one
# OM-like simulation (alfonsino, F history up-down-constant). Checks that the
# fit converges, that the reference F from the taped Newton iterations is the
# maximum of the equilibrium yield on the grid, and that the terminal SB/SB0 is
# close to the truth.
suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
cand <- c(if (!is.null(of)) file.path(dirname(of), "..", "sim"), "ci/sim", "new/ci/sim")
here <- cand[file.exists(file.path(cand, "om_sim.R"))][1]
source(file.path(here, "om_sim.R")); source(file.path(here, "om_sim_spm.R"))
lh <- om_lh; st <- om_structure(); set.seed(2001)
sim <- om_simulate_spm(st, sigR = 0.3, fshape = "fudc")
lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
lhp <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))), a = lh$a, b = lh$b)
stk <- stocklen(lfd, lhp, m_model = "constant")
kd <- match(sim$data_years, sim$years); cbg <- colSums(sim$catch_w[kd, , drop = FALSE]); cbg <- cbg / sum(cbg)
spm <- list(catch = sim$catch, index = sim$index, index.gear = "Gillnet", h = lh$h, sel.weights = "catch",
            xgrid = c(0, exp(seq(log(0.02), log(40), length.out = 20))), newton.it = 6)
t0 <- Sys.time()
fit <- FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                           settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3, spm = spm))
cat("FLicc-SPM fit:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
s <- spm_flicc(fit)
print(round(s$refpts, 4))
rp <- fit$report
dep <- s$ts$est[s$ts$quantity == "SB/SB0"]
cat(sprintf("terminal SB/SB0 true %.3f est %.3f | BmsyK true %.3f per-recruit %.3f\n",
            tail(sim$depletion, 1), tail(dep, 1), sim$sbmsy / sim$sb0, rp$BmsyK_pr))
stopifnot(
  fit$diagnostics$conv_tier == "ok",
  all(is.finite(dep)),
  # taped Newton: equilibrium yield at Fref at least the best grid value
  rp$Yref >= max(rp$Yg) * (1 - 1e-8),
  rp$Fref > 0, rp$m_pr > 0, rp$r_pr > 0,
  abs(log(tail(dep, 1) / tail(sim$depletion, 1))) < 0.25,
  abs(log(rp$BmsyK_pr / (sim$sbmsy / sim$sb0))) < 0.2
)
# spm off: SPM parameters mapped off, objective identical to FLicc without spm
f0 <- FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                          settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3), se = FALSE)
stopifnot(!any(grepl("^spm_", names(f0$opt$par))), is.null(f0$report$P))
cat("test_spm passed\n")
