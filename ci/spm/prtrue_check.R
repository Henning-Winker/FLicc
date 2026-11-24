suppressPackageStartupMessages({ library(FLCore); library(FLicc) })
source("./ci/sim/om_sim.R"); source("./ci/sim/om_sim_spm.R")
lh <- om_lh; st <- om_structure(); set.seed(2001); sim <- om_simulate_spm(st, sigR = 0.3, fshape = "fudc")
x <- readRDS("ci/spm/fits/case_spr_fudc_2001.rds"); f <- x$f; obj <- f$obj
par <- obj$env$last.par.best; nm <- names(par)
rp_at <- function(p) { r <- obj$report(p); c(Fref = r$Fref, BmsyK = r$BmsyK_pr, Hmsy = r$Hmsy_pr, r = r$r_pr, m = r$m_pr) }
est <- rp_at(par)
L <- lh$linf
tru <- par
tru[nm == "log_Linf"] <- log(L)
tru[nm == "log_Mk"] <- log(lh$M / lh$k)
tru[nm == "log_Galpha"] <- -2 * log(0.1)
sm <- which(nm == "Sm")
tru[sm] <- c(log(23.35 / L), log(3.19 / L), log(11.22 / L), log(29.7 / L), log((34.1 - 29.7) / L))
cat("estimated vs true (Linf, M/k, CVL, Sm):\n")
print(rbind(est = par[c(which(nm %in% c("log_Linf", "log_Mk", "log_Galpha")), sm)], true = tru[c(which(nm %in% c("log_Linf", "log_Mk", "log_Galpha")), sm)]), digits = 3)
r_true <- rp_at(tru)
# true F weights of the gear selectivities (apical F shares) instead of catch shares
td <- obj$env$data
w_true <- st$Fap / sum(st$Fap)
cat("selectivity weights: catch shares", round(td$spm_selw, 3), " | true apical F shares", round(w_true, 3), "\n")
obj$env$data$spm_selw <- w_true
r_true_w <- rp_at(tru)
obj$env$data$spm_selw <- td$spm_selw
# OM truth on the same definitions
Fmsy <- om_fmsy(st); h <- lh$h; sb0 <- om_sbpr(st, 0); sbm <- om_sbpr(st, Fmsy)
Rm <- (4 * h * sbm - sb0 * (1 - h)) / ((5 * h - 1) * sbm)
sj <- rowSums(st$Fsel); Fap_j <- max(sj)
om <- c(Fref = Fmsy * Fap_j, BmsyK = Rm * sbm / sb0, Hmsy = unname(sim$Hmsy["start"]), r = NA, m = NA)
mfun <- function(b) uniroot(function(m) m^(-1 / (m - 1)) - b, c(0.05, 0.999))$root
om["m"] <- mfun(om["BmsyK"]); om["r"] <- om["m"] * om["Hmsy"]
out <- rbind(`OM truth` = om, `per-recruit, estimated pars` = est, `per-recruit, true pars` = r_true,
             `per-recruit, true pars + true F weights` = r_true_w)
print(round(out, 4))
