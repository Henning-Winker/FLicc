suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source("./ci/sim/om_sim.R"); source("./ci/sim/om_sim_spm.R")
st <- om_structure()
fspr40 <- uniroot(function(F) om_spr(st, F) - 0.4, c(1e-4, 2))$root
cat(sprintf("OM: Fmsy %.3f, Fspr40 %.3f (Fspr40/Fmsy %.2f); SPR at Fmsy %.2f\n", om_fmsy(st), fspr40, fspr40 / om_fmsy(st), om_spr(st, om_fmsy(st))))
sc <- data.frame(seed = c(2001, 11, 12, 11, 12, 11, 12), shape = c("fudc", "rebuild", "rebuild", "decline", "decline", "flat", "flat"))
R <- list(); T <- list()
for (i in seq_len(nrow(sc))) {
  a <- readRDS(sprintf("ci/spm/fits/spm_eval_%d_%s_spr0.rds", sc$seed[i], sc$shape[i]))
  b <- readRDS(sprintf("ci/spm/fits/case_both_%s_%d.rds", sc$shape[i], sc$seed[i]))
  sim <- b$sim; kd <- match(sim$data_years, sim$years)
  stopifnot(all.equal(a$sim$catch, sim$catch))
  lh <- om_lh
  lfd <- FLQuants(lapply(sim$lfd, function(x) FLQuant(x, dimnames = list(len = sim$lb, year = sim$data_years))))
  lhp <- FLPar(linf = lh$linf, k = lh$k, M = lh$M, L50 = lh$linf * (1 - exp(-lh$k * (lh$tm - lh$t0))), a = lh$a, b = lh$b)
  stk <- stocklen(lfd, lhp, m_model = "constant")
  cbg <- colSums(sim$catch_w[kd, , drop = FALSE]); cbg <- cbg / sum(cbg)
  f0 <- suppressMessages(fiticc(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                                settings = list(pop_model = "gtg", obs_model = "dm", CVL = 0.1, maxsd = 3)))
  ff_flicc <- as.numeric(harvest(flicc2FLStockR(f0, rel = TRUE)))
  ff_om_spr <- sim$Fbar[kd] / fspr40
  s <- spm_flicc(b$f); ff_spm <- s$ts$est[s$ts$quantity == "F/Fmsy"][kd]
  ff_om_msy <- sim$ffmsy[kd]
  e1 <- log(ff_flicc / ff_om_spr); e2 <- log(ff_spm / ff_om_msy)
  lab <- sprintf("%s %d", sc$shape[i], sc$seed[i])
  R[[i]] <- data.frame(scen = lab, `FLicc F/Fspr40 median err` = median(e1), `terminal` = tail(e1, 1),
                       `SPM F/Fmsy median err` = median(e2), `terminal ` = tail(e2, 1), check.names = FALSE)
  T[[i]] <- data.frame(scen = lab, year = sim$data_years, om_fspr = ff_om_spr, flicc = ff_flicc, om_fmsy = ff_om_msy, spm = ff_spm)
}
print(do.call(rbind, R), digits = 2, row.names = FALSE)
X <- do.call(rbind, R)
cat(sprintf("\nmedian |err| over scenarios (median over years): FLicc F/Fspr40 %.2f | SPM F/Fmsy %.2f\n",
            median(abs(X[[2]])), median(abs(X[[4]]))))
cat(sprintf("terminal year |err|: FLicc %.2f | SPM %.2f\n", median(abs(X[[3]])), median(abs(X[[5]]))))
saveRDS(do.call(rbind, T), "/tmp/fspr_traj.rds")
