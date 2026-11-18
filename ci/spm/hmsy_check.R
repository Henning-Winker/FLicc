suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source("./ci/sim/om_sim.R"); source("./ci/sim/om_sim_spm.R")
lh <- om_lh; st <- om_structure()
sc <- data.frame(seed = c(2001, 11, 12, 11, 12, 11, 12), shape = c("fudc", "rebuild", "rebuild", "decline", "decline", "flat", "flat"))
out <- list()
for (i in seq_len(nrow(sc))) {
  set.seed(sc$seed[i]); sim <- om_simulate_spm(st, sigR = 0.3, fshape = sc$shape[i])
  x <- readRDS(sprintf("ci/spm/fits/case_spr_%s_%d.rds", sc$shape[i], sc$seed[i]))
  stopifnot(all.equal(sim$catch, x$sim$catch))             # same draws as the saved fit
  rp <- x$f$report
  out[[i]] <- data.frame(scen = paste(sc$shape[i], sc$seed[i]), year = sim$years,
    ffmsy = sim$ffmsy, hh_start = sim$hhmsy_start, hh_mean = sim$hhmsy_mean, hh_spm = rp$HHmsy,
    H_start = sim$H_start, H_spm = rp$H)
  if (i == 1) cat(sprintf("Hmsy true: start-of-year %.4f, mean-year %.4f | SPM Hmsy %.4f (per-recruit %.4f)\n",
                          sim$Hmsy["start"], sim$Hmsy["mean"], rp$Hmsy_spm, rp$Hmsy_pr))
}
D <- do.call(rbind, out); saveRDS(D, "/tmp/hcheck.rds")
agg <- do.call(rbind, lapply(split(D, D$scen), function(d) data.frame(scen = d$scen[1],
  `Fbar/Fmsy` = median(log(d$hh_spm / d$ffmsy)), `H/Hmsy start` = median(log(d$hh_spm / d$hh_start)),
  `H/Hmsy mean` = median(log(d$hh_spm / d$hh_mean)), `H level` = median(log(d$H_spm / d$H_start)), check.names = FALSE)))
cat("median log(SPM / truth) over all years:\n"); print(agg, digits = 2, row.names = FALSE)
