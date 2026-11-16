suppressPackageStartupMessages({ library(FLCore); library(FLicc, ) })
source("./ci/sim/om_sim.R"); source("./ci/sim/om_sim_spm.R")
args <- commandArgs(TRUE); seed <- as.integer(args[1]); shape <- args[2]; cfg <- args[3]
lh <- om_lh; st <- om_structure(); set.seed(seed)
sim <- om_simulate_spm(st, sigR = 0.3, fshape = shape)
lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
lhp <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))), a = lh$a, b = lh$b)
stk <- stocklen(lfd, lhp, m_model = "constant")
kd <- match(sim$data_years, sim$years); cbg <- colSums(sim$catch_w[kd, , drop = FALSE]); cbg <- cbg / sum(cbg)
base <- list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3)
f0 <- FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg, settings = base, se = FALSE)
lk <- switch(cfg, spr0 = list(link = "spr", link.rho = c(0, 0)), sprAR = list(link = "spr", link.rho = c(0.7, 0)),
  FAR = list(link = "F", link.rho = c(0, 0.7)), bothAR = list(link = c("spr", "F"), link.rho = c(0.7, 0.7)))
spm <- list(catch = sim$catch, index = sim$index, index.gear = "Gillnet", h = lh$h, link = lk$link, link.rho = lk$link.rho, sel.weights = if (length(args) > 3) args[4] else "F",
            xgrid = c(0, exp(seq(log(0.02), log(40), length.out = 20))), newton.it = 6)
t0 <- Sys.time()
f1 <- FLicc:::fiticc_core(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                          settings = modifyList(base, list(spm = spm)), se = TRUE)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
s <- spm_flicc(f1)
g <- function(q) s$ts$est[s$ts$quantity == q]
d <- data.frame(year = sim$years, dep = sim$depletion, dep_spm = g("SB/SB0"), ffmsy = sim$ffmsy, ffmsy_spm = g("F/Fmsy"),
                bbmsy = sim$bbmsy, bbmsy_spm = g("B/Bmsy"))
d$dspr_flicc <- NA; d$dspr_flicc[kd] <- pmax((4 * lh$h * c(f0$report$spr_y) - (1 - lh$h)) / (5 * lh$h - 1), 0)
cat(sprintf("seed %d %s: %.1f s conv %d tier %s maxgrad %.2g | BmsyK pr %.3f true %.3f | r %.3f m %.2f sdb %.3f\n",
            seed, shape, secs, f1$opt$convergence, f1$diagnostics$conv_tier, f1$max_gradient,
            s$refpts["BmsyK_pr"], sim$sbmsy / sim$sb0, s$refpts["r"], s$refpts["m"], s$refpts["sdb"]))
print(round(d[d$year %in% c(1975, 1990, 2000, 2005, 2010, kd + 1970), ], 3), row.names = FALSE)
saveRDS(list(sim = sim, f0 = f0, f1 = f1, d = d), sprintf("/tmp/spm_eval_%d_%s_%s.rds", seed, shape, cfg))
ty <- nrow(d)
cat(sprintf("SUMMARY %s %s seed %d |", cfg, shape, seed)); cat(sprintf(" time %.0fs rho %.2f/%.2f |", secs, s$refpts["rho_spr"], s$refpts["rho_F"])); cat(sprintf("%s seed %d | terminal dep true %.3f spm %.3f flicc-Dspr %.3f | B/Bmsy true %.2f spm %.2f | F/Fmsy true %.2f spm %.2f | conv %s\n",
  shape, seed, d$dep[ty], d$dep_spm[ty], d$dspr_flicc[ty], d$bbmsy[ty], d$bbmsy_spm[ty], d$ffmsy[ty], d$ffmsy_spm[ty], f1$diagnostics$conv_tier))
