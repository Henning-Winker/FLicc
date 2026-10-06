#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Tests for FLicc dev-lbsrp additions (versions 1.0.6 - 1.0.9)
#  1. GTG uses estimated Linf, M/K and CVL
#  2. GTG tails: maxsd / ngtg (defaults as in LBSPR: maxsd = 2, ngtg = 13)
#  3. Optional warning for fish beyond the largest GTG Linf (tail_warning)
#  4. Robust likelihood options: rob_eps and Lplus
#  5. Steepness: spr2brel(), sprbrel_flicc(), fbrel_flicc(),
#     flicc2FLStockR(s = )
#  6. Speed: Fspr reported by the fit (spr_ref) and TMB per-recruit
#     calculations (brp_tmb_flicc)
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

library(FLicc)

data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)),
                  unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859,
               a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")

sel_fun <- c("dsnormal", "logistic")
cbg     <- c(0.7, 0.3)
base    <- list(prior_sigmaF = c(log(0.5), 0.3, 1), CVL = 0.1)
fitf <- function(...) fiticc(lfd, stklen, sel_fun = sel_fun,
                             catch_by_gear = cbg, settings = c(base, list(...)))
lastspr <- function(fit) round(tail(c(fit$report$spr), 1), 4)

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# 1. GTG uses estimated Linf, M/K and CVL
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
fit.fix <- fitf()                                            # all fixed
fit.est <- fitf(linf.sd = 0.1, Mk.sd = 0.2, CVL.sd = 0.2)     # estimated

rbind(fixed     = c(Linf = c(fit.fix$report$lhpar["linf"]), Mk = c(fit.fix$report$lhpar["Mk"]),
                    CVL = c(fit.fix$report$pars["CVL"]), SPR = lastspr(fit.fix)),
      estimated = c(Linf = c(fit.est$report$lhpar["linf"]), Mk = c(fit.est$report$lhpar["Mk"]),
                    CVL = c(fit.est$report$pars["CVL"]), SPR = lastspr(fit.est)))

# GTG Linf groups follow the estimates: Linf * (1 + CVL * z)
z <- seq(-2, 2, length.out = 13)                             # default maxsd = 2, ngtg = 13
stopifnot(isTRUE(all.equal(fit.est$report$gtgLinfs,
  c(fit.est$report$lhpar["linf"]) * (1 + c(fit.est$report$pars["CVL"]) * z))))
# estimated parameters actually change the population model
stopifnot(abs(lastspr(fit.est) - lastspr(fit.fix)) > 1e-4)

# fit$report$N has 'len' as its first dimension
stopifnot(names(dimnames(fit.fix$report$N))[1] == "len")

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# 2. GTG tails: maxsd changes the actual Linf spread, ngtg only resolution
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
fit.wide <- fitf(maxsd = 3, ngtg = 19)                       # CVL close to actual spread
fit.ng30 <- fitf(maxsd = 3, ngtg = 30)
c(default_maxsd2_ngtg13 = lastspr(fit.fix), maxsd3_ngtg19 = lastspr(fit.wide),
  maxsd3_ngtg30 = lastspr(fit.ng30))

# actual CV of the GTG Linf distribution versus CVL = 0.1
cv_gtg <- function(fit) {
  L <- fit$report$gtgLinfs; w <- fit$report$recP
  sqrt(sum(w * (L - sum(w * L))^2)) / sum(w * L)
}
round(c(maxsd2 = cv_gtg(fit.fix), maxsd3 = cv_gtg(fit.wide)), 3)     # ~0.091 vs ~0.099

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# 3. Warning when observed fish lie beyond the largest GTG Linf
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
w <- NULL
fit.nowarn <- fitf(maxsd = 0.5, ngtg = 7)                    # off by default: silent
fit.narrow <- withCallingHandlers(fitf(maxsd = 0.5, ngtg = 7, tail_warning = TRUE),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
print(w[grepl("largest GTG Linf", w)])
stopifnot(any(grepl("largest GTG Linf", w)))

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# 4. Robust likelihood options
#    rob_eps: p = (1 - rob_eps) * p + rob_eps / nbins
#    Lplus:   pool bins with lower bound >= Lplus into a plus group
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
fit.r0    <- fitf(obs_model = "dm", rob_eps = 0)
fit.r001  <- fitf(obs_model = "dm", rob_eps = 0.001)
fit.r01   <- fitf(obs_model = "dm", rob_eps = 0.01)
fit.plus  <- fitf(obs_model = "dm", Lplus = 50)
fit.both  <- fitf(obs_model = "dm", Lplus = 50, rob_eps = 0.001)
fit.dm    <- fitf(obs_model = "dm")

rob <- list(dm = fit.dm, rob0 = fit.r0, rob0.001 = fit.r001, rob0.01 = fit.r01,
            plus50 = fit.plus, plus50_rob0.001 = fit.both)
data.frame(nll = sapply(rob, function(f) round(f$opt$objective, 2)),
           SPR = sapply(rob, lastspr),
           Fspr40 = sapply(rob, function(f) round(fspr_flicc(f), 4)))
# rob_eps = 0 is the default likelihood
stopifnot(abs(fit.r0$opt$objective - fit.dm$opt$objective) < 1e-6)
plot_spr(rob)
plot_len(fit.both)          # predicted LFD stays unpooled in plots

# all observation models accept both options
fit.nb <- fitf(obs_model = "nb", Lplus = 50, rob_eps = 0.001)
fit.mn <- fitf(obs_model = "mn", Lplus = 50, rob_eps = 0.001)
c(nb = lastspr(fit.nb), mn = lastspr(fit.mn))

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# 5. Steepness: biomass-scale status and reference points
#    B/B0 = (4 s SPR - (1 - s)) / (5 s - 1) under Beverton-Holt
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
s <- 0.75
sprbrel_flicc(c(1, 0.5, 0.25), s = s)            # SPR (%) giving B/B40 = 1, 0.5, 0.25
spr2brel(c(0.4, 0.2, 0.16, 0.1), s = s)          # B/B40 at these SPR values
stopifnot(isTRUE(all.equal(spr2brel(sprbrel_flicc(0.25, s = s) / 100, s = s), 0.25)))

fbrel_flicc(fit.fix, brel = 0.25, s = s)         # F giving 0.25 * B40 (Blim proxy)
fspr_flicc(fit.fix, spr = 10)                    # compare: Fspr10 is much higher

spr2brel(fit.fix$report$spr, s = s)              # status as B / B40 by year

stk.spr <- flicc2FLStockR(fit.fix, rel = TRUE)           # SPR / SPR40 (as before)
stk.b   <- flicc2FLStockR(fit.fix, rel = TRUE, s = s)    # B / B40 under s
rbind(SPR_scale = tail(c(stk.spr@stock), 3), B_scale = tail(c(stk.b@stock), 3))
stk.b@refpts
plot_LBAdvice(stk.b)
stopifnot(isTRUE(all.equal(c(stk.b@stock), c(spr2brel(c(fit.fix$report$spr), s = s)))))

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# 6. Speed: Fspr from the fit and per-recruit calculations in TMB
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
fit.sp <- fitf(spr_ref = c(40, 20, 10))
fit.sp$report$Fspr                               # reported by the fit

system.time(f40 <- fspr_flicc(fit.sp, spr = 40))     # stored: instant
system.time(f30 <- fspr_flicc(fit.sp, spr = 30))     # other target: one TMB call
b <- brp_tmb_flicc(fit.sp, Fseq = c(0, 0.05, 0.1), spr = c(40, 30))
b$Fspr; b$SBPR / b$SBPR0                            # SPR at Fseq

# compare with the previous R implementation
timeit <- function() c(
  fspr40  = system.time(a1 <<- fspr_flicc(fit.sp, spr = 40))[["elapsed"]],
  prbrp   = system.time(a2 <<- prbrp_flicc(fit.sp))[["elapsed"]],
  eqstk   = system.time(a3 <<- eqstklen(fit.sp, s = s))[["elapsed"]],
  LBIspr  = system.time(a4 <<- LBIspr(fit.sp))[["elapsed"]],
  stkR    = system.time(a5 <<- flicc2FLStockR(fit.sp))[["elapsed"]])
options(FLicc.brp_tmb = FALSE); t.old <- timeit(); old <- list(a1, c(a2$SBPR), c(a3@refpts))
options(FLicc.brp_tmb = TRUE);  t.new <- timeit(); new <- list(a1, c(a2$SBPR), c(a3@refpts))
rbind(old = t.old, new = t.new, speedup = round(t.old / pmax(t.new, 1e-3)))

# results agree (Fspr within the old uniroot tolerance)
reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12), na.rm = TRUE)
stopifnot(reldiff(new[[1]], old[[1]]) < 1e-3,
          reldiff(new[[2]], old[[2]]) < 1e-8,
          reldiff(new[[3]], old[[3]]) < 1e-3)

cat("\nAll dev-lbsrp tests passed.\n")
