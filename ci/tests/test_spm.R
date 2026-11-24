# FLicc-SPM prototype: joint fit of the gtg length model (F fixed effects,
# a nuisance that tunes SPR) and a state-space surplus production model
# (states random) to catch and CPUE, on one
# OM-like simulation (alfonsino, F history up-down-constant). Checks that the
# fit converges, that the reference F from the taped Newton iterations is the
# maximum of the equilibrium yield on the grid, and that the terminal SB/SB0 is
# close to the truth. Also: identical sexes reproduce the one-sex fit, and growth
# from an age-length key (k estimated) works with the SPM on.
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
# penalised likelihood (re = FALSE): states as fixed effects, process SD fixed; close to the Laplace fit
spm_p <- spm; spm_p$re <- FALSE; spm_p$sdb <- c(exp(fit$opt$par[["spm_logsdb"]]), 0.5, 1)
tp <- system.time(fp <- FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                  settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3, spm = spm_p)))[3]
rpp <- fp$report
cat(sprintf("penalised: %.1f s, conv %d, SB/SB0 end %.3f (Laplace %.3f, true %.3f), SPR end %.3f (Laplace %.3f)\n", tp,
            fp$opt$convergence, tail(rpp$P, 1), tail(rp$P, 1), tail(sim$depletion, 1), tail(rpp$spr_y, 1), tail(rp$spr_y, 1)))
stopifnot(is.null(fp$obj$env$random) || !any(names(fp$obj$env$par)[fp$obj$env$random] == "spm_logP"),
          !any(c("spm_logsdb", "spm_logsdi") %in% names(fp$opt$par)), fp$opt$convergence == 0,
          abs(log(tail(rpp$P, 1) / tail(rp$P, 1))) < 0.15, abs(log(tail(rpp$spr_y, 1) / tail(rp$spr_y, 1))) < 0.1)
# sex structure: identical sexes reproduce the one-sex FLicc-SPM fit
stk2 <- stocklen_sex(lfd, list(F = lhp, M = lhp), m_model = "constant")
fs <- fiticc(lfd, stk2, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                          settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3, spm = spm), verbose = FALSE)
cat(sprintf("sex identity: nll one-sex %.6f two-sex %.6f\n", fit$opt$objective, fs$opt$objective))
stopifnot(fs$tmb_data$nsex == 2, abs(fit$opt$objective - fs$opt$objective) < 1e-4,
          abs(tail(fs$obj$report()$P, 1) / tail(rp$P, 1) - 1) < 1e-3)
# ALK: growth (k) estimated with the SPM on; the SPM uses the estimated k
set.seed(5); yk <- tail(sim$data_years, 1); ages <- rep(1:12, each = 40)
la <- lh$linf * (1 - exp(-lh$k * (ages + 0.5 - lh$t0))) * exp(rnorm(length(ages), 0, 0.1))
key <- alk_flicc(data.frame(age = ages, len = floor(la)), gears = "Gillnet", year = yk)
fa <- fiticc(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg, alk = key,
                          settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3, spm = spm,
                                          alk_t0 = lh$t0), verbose = FALSE)
ka <- exp(fa$opt$par[["log_k"]]); ra <- fa$obj$report()
cat(sprintf("ALK + SPM: conv %s, k est %.3f (true %.3f), r_pr %.3f\n", fa$diagnostics$conv_tier, ka, lh$k, ra$r_pr))
stopifnot(fa$tmb_data$alk_on == 1, "log_k" %in% names(fa$opt$par), abs(log(ka / lh$k)) < 0.3,
          fa$diagnostics$conv_tier == "ok", is.finite(ra$r_pr), all(is.finite(ra$P)))
# removed / not estimable options fail with a clear message
for (bad in list(list(link = c("spr", "F")), list(link.rho = NA))) {
  e <- tryCatch(FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3, spm = modifyList(spm, bad)),
                se = FALSE), error = function(e) conditionMessage(e))
  stopifnot(is.character(e), grepl("memory/spm-flink|cannot be estimated", e))
}
cat("test_spm passed\n")
