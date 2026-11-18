suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source("./ci/sim/om_sim.R"); source("./ci/sim/om_sim_spm.R")
a <- commandArgs(TRUE); seed <- as.integer(a[1]); shape <- a[2]; cfg <- a[3]
lh <- om_lh; st <- om_structure(); set.seed(seed)
sim <- om_simulate_spm(st, sigR = 0.3, fshape = shape)
lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
lhp <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))), a = lh$a, b = lh$b)
stk <- stocklen(lfd, lhp, m_model = "constant")
kd <- match(sim$data_years, sim$years); cbg <- colSums(sim$catch_w[kd, , drop = FALSE]); cbg <- cbg / sum(cbg)
spm <- list(catch = sim$catch, index = sim$index, index.gear = "Gillnet", h = lh$h,
            xgrid = c(0, exp(seq(log(0.02), log(40), length.out = 20))), newton.it = 6)
spm <- modifyList(spm, switch(cfg,
  spr     = list(),
  sprSdb  = list(sdb = c(0.1, 0.25, 1)),
  FAR     = list(link = "F", link.rho = c(0, 0.7)),
  FARsdb  = list(link = "F", link.rho = c(0, 0.7), sdb = c(0.1, 0.25, 1)),
  both    = list(link = c("spr", "F"), link.rho = c(0.7, 0.7)),
  bothSdb = list(link = c("spr", "F"), link.rho = c(0.7, 0.7), sdb = c(0.1, 0.25, 1))))
t0 <- Sys.time()
f <- tryCatch(FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
       settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3, spm = spm)), error = function(e) e)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (inherits(f, "error")) { cat(sprintf("RES %-7s %-7s %4d ERROR %s\n", cfg, shape, seed, conditionMessage(f))); quit() }
s <- spm_flicc(f); g <- function(q) tail(s$ts$est[s$ts$quantity == q], 1)
cat(sprintf("RES %-7s %-7s %4d %4.0fs %-4s K %6.2f r %.3f sdb %.3f | dep %.3f/%.3f | B/Bmsy %.2f/%.2f | F/Fmsy %.2f/%.2f\n",
  cfg, shape, seed, secs, f$diagnostics$conv_tier, s$refpts["K"], s$refpts["r"], s$refpts["sdb"],
  tail(sim$depletion, 1), g("SB/SB0"), tail(sim$bbmsy, 1), g("B/Bmsy"), tail(sim$ffmsy, 1), g("F/Fmsy")))
saveRDS(list(sim = sim, f = f), sprintf("/tmp/case_%s_%s_%d.rds", cfg, shape, seed))
