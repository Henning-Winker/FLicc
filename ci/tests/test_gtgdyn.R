# pop_model = "gtg.dyn" (dynamic length-indexed GTG)
#  1. reference: FLicc gtg.dyn reproduces the standalone prototype (ci/gtgdyn,
#     itself checked against the Python prototype) for identical inputs
#  2. at constant F and recruitment, gtg.dyn reproduces the equilibrium gtg
#  3. fits on the alfonsino data (two gears) and downstream functions
library(FLicc)
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "ci/tests")
ddir <- normalizePath(file.path(here, "..", "gtgdyn", "data"))

# parameter list of an object at given values, all mapped except log_phi
# (TMB needs one free parameter); returns obj$report()
rep_at <- function(td, pl) {
  td <- td[setdiff(names(td), c("gear_names", "year_names", "Linf_init", "Mk_init", "CVL_init", "Galpha_init"))]
  map <- lapply(pl, function(x) factor(rep(NA, length(x)))); map$log_phi <- NULL
  obj <- TMB::MakeADFun(td, pl, map = map, DLL = "FLicc", silent = TRUE)
  obj$report()
}
add_rw <- function(td) {
  td$prior_sigmaF_mean <- log(0.3); td$prior_sigmaF_sd <- 0; td$prior_sigmaF_use <- 1L
  td$FM_min <- 0.05; td$FM_max <- 4; td$FMpen_sd <- 0
  td
}

## ---- 1. reference test against the standalone prototype --------------------
obs <- as.matrix(read.csv(file.path(ddir, "obs_Rdev_1.csv"), header = FALSE))
rp  <- read.csv(file.path(ddir, "ref_par.csv"))
th  <- read.csv(file.path(ddir, "ref_sel.csv"))$theta
ref <- as.matrix(read.csv(file.path(ddir, "ref_pred_dyng.csv"), header = FALSE))
lens <- seq(0, 66, by = 2); yrs <- 2005:2024
lfd1 <- FLQuants(G1 = FLQuant(t(obs), dimnames = list(len = lens, year = yrs)))
k <- 0.08; linf <- 55.662
lhp <- FLPar(linf = linf, k = k, M = 0.162, L50 = linf * (1 - exp(-k * 6)), a = 3.24e-05, b = 2.932)
stk1 <- stocklen(lfd1, lhp, m_model = "constant")
set_ref <- list(CVL = 0.1, pop_model = "gtg.dyn", obs_model = "mn", ngtg = 19, maxsd = 3,
                dyn_npre = 10, dyn_nq = 2, dyn_db = 0.25, dyn_h = 2, dyn_Tcap = 40, dyn_rec = "iid", GL = 30)
cbg1 <- FLQuants(G1 = quantSums(lfd1[[1]]) * 0 + 1)
td <- add_rw(FLicc:::data_tmb_flicc(lfd1, stk1, sel_fun = "dsnormal", catch_by_gear = cbg1,
                                    settings = set_ref))
pl <- FLicc:::init_tmb_flicc(td)
pl$Sm <- c(log(th[1] / linf), th[2] - log(linf), th[3] - log(linf))
pl$log_Fk[] <- rp$logF[11:30] - log(k)          # FLicc F is F / k
pl$log_Fpre[] <- rp$logF[1:10] - log(k)
pl$log_Rdev <- rp$logR
r1 <- rep_at(td, pl)
d_ref <- max(abs(t(r1$plen[, , 1]) - ref))
cat(sprintf("REFERENCE: max |FLicc gtg.dyn - prototype| = %.2e\n", d_ref))
stopifnot(d_ref < 1e-8)

## ---- 2. constant F: gtg.dyn = gtg ------------------------------------------
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
cbg <- c(0.7, 0.3)
fitf <- function(...) fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                             settings = c(list(CVL = 0.1, obs_model = "dm"), list(...)))
t_gtg <- system.time(fg <- fitf(pop_model = "gtg"))[["elapsed"]]
cbgq <- FLQuants(Map(function(x, y) { q <- quantSums(x); q[] <- y; q }, lfd, cbg))
for (h in c(1, 4)) {
  tdd <- add_rw(FLicc:::data_tmb_flicc(lfd, stklen, c("dsnormal", "logistic"), cbgq,
                                       settings = list(CVL = 0.1, GL = 30, pop_model = "gtg.dyn", dyn_h = h)))
  tdg <- add_rw(FLicc:::data_tmb_flicc(lfd, stklen, c("dsnormal", "logistic"), cbgq,
                                       settings = list(CVL = 0.1, GL = 30, pop_model = "gtg")))
  p0 <- fg$obj$env$parList(fg$opt$par)
  Fc <- matrix(colMeans(p0$log_Fk), nrow(p0$log_Fk), ncol(p0$log_Fk), byrow = TRUE)
  pg <- FLicc:::init_tmb_flicc(tdg); pg[c("Sm", "log_Mk", "log_phi")] <- p0[c("Sm", "log_Mk", "log_phi")]
  pg$log_Fk[] <- Fc
  pd <- FLicc:::init_tmb_flicc(tdd); pd[c("Sm", "log_Mk", "log_phi")] <- p0[c("Sm", "log_Mk", "log_phi")]
  pd$log_Fk[] <- Fc; pd$log_Fpre[] <- rep(Fc[1, ], each = nrow(pd$log_Fpre))
  rg <- rep_at(tdg, pg); rd <- rep_at(tdd, pd)
  d_eq <- max(abs(rg$plen - rd$plen))
  cat(sprintf("CONSTANT F (dyn_h = %g): max |gtg.dyn - gtg| in predicted proportions = %.1e (peak %.3f); SPR identical: %s\n",
              h, d_eq, max(rg$plen), isTRUE(all.equal(rg$spr_y, rd$spr_y))))
  if (h == 1) stopifnot(d_eq < 2e-4)
  stopifnot(isTRUE(all.equal(rg$spr_y, rd$spr_y)))
}

## ---- 3. fits on the alfonsino data -------------------------------------------
summ <- function(f, secs) c(secs = round(secs, 1), nll = round(f$opt$objective, 2), conv = f$opt$convergence,
                            maxgrad = signif(f$opt$max_gradient, 2),
                            sigmaR = if (!is.null(f$report$sigmaR)) round(f$report$sigmaR, 3) else NA,
                            SPR_last = round(tail(c(f$report$spr), 1), 3))
runs <- list(
  gtg              = list(pop_model = "gtg"),
  dyn_default      = list(pop_model = "gtg.dyn"),
  dyn_sigR_fixed   = list(pop_model = "gtg.dyn", sigmaR = 0.5, R_re = FALSE),
  dyn_constR       = list(pop_model = "gtg.dyn", dyn_rec = "none"),
  dyn_Rrw          = list(pop_model = "gtg.dyn", dyn_rec = "rw", sigmaR = 0.2, R_re = FALSE),
  dyn_noSE         = list(pop_model = "gtg.dyn", sigmaR = 0.5, R_re = FALSE)
)
fits <- list(); tab <- NULL
for (nm in names(runs)) {
  t0 <- system.time(f <- fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = cbg,
                                settings = c(list(CVL = 0.1, obs_model = "dm"), runs[[nm]]),
                                se = nm != "dyn_noSE"))[["elapsed"]]
  fits[[nm]] <- f
  tab <- rbind(tab, c(model = nm, summ(f, t0)))
}
print(as.data.frame(tab), row.names = FALSE)
stopifnot(all(as.numeric(tab[, "conv"]) == 0))

fd <- fits$dyn_default
cat("\nSPR by year:\n")
print(round(rbind(gtg = c(fits$gtg$report$spr), gtg.dyn = c(fd$report$spr),
                  dyn_sigR_fixed = c(fits$dyn_sigR_fixed$report$spr),
                  dyn_constR = c(fits$dyn_constR$report$spr)), 3))
cat("\ngtg.dyn dynamic SSB / unfished (per mean recruit):\n"); print(round(c(fd$report$ssb_dyn), 3))
cat("\ngtg.dyn recruitment multipliers:\n"); print(round(c(fd$report$rec_dev), 2))
cat("\nsigmaR:", fd$report$sigmaR, "\n")
stopifnot(!is.null(fd$report$Ndyn), !is.null(fd$report$ssb_dyn), !is.null(fd$report$rec_dev),
          all(dim(fd$report$Ndyn)[1:2] == dim(fd$report$N)[1:2]))

## downstream functions on a gtg.dyn fit
st <- flicc2FLStockR(fd, spr.btgt = 0.4, rel = TRUE)
fs <- fspr_flicc(fd)
print(fs)
lbi <- LBIspr(fd)
print(lbi)
stopifnot(is.finite(c(fs)[1]))
cat("\ngtg.dyn tests passed\n")
