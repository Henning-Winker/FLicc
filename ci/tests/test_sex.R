# Sex-structured FLicc: identity with the one-sex model, sexed composition
# types (SS3 gender codes 0-3) and recovery of a dimorphic stock
library(FLicc)
data(alfonsino)
lfd <- lfdess(FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE), ess.g = 200)
lh <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stk <- stocklen(lfd, lh, m_model = "constant")
stk2 <- stocklen_sex(lfd, list(F = lh, M = lh), m_model = "constant")
sf <- c("dsnormal", "logistic"); cb <- c(0.7, 0.3)
sexq <- function(x, f = 0.5, units = c("F", "M")) {
  d <- dimnames(x); d$unit <- units; q <- FLQuant(0, dimnames = d)
  if ("F" %in% units) q[, , "F"] <- x * f
  if ("M" %in% units) q[, , "M"] <- x * (1 - f)
  if ("U" %in% units) q[, , "U"] <- x * (1 - f)
  q
}

## 0. sexed length data from a long table with a sex column (F/M/U, female/male, SS3 codes)
dl <- as.data.frame(lfd_alfonsino)
dl2 <- rbind(transform(dl[dl$qname == "Trawl", ], sex = "U"),
             transform(dl[dl$qname == "Gillnet", ], data = data * 0.6, sex = "female"),
             transform(dl[dl$qname == "Gillnet", ], data = data * 0.4, sex = 2))
ls <- FLQuantLen(lfd_long_to_wide(dl2), unit = "cm", midL = FALSE)
stopifnot(identical(names(ls), c("Trawl", "Gillnet")), identical(dimnames(ls$Gillnet)$unit, c("F", "M")),
          abs(sum(ls$Gillnet[, , "F"]) / sum(ls$Gillnet) - 0.6) < 1e-8)
lw <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
stopifnot(dimnames(lw[[1]])$unit == "unique")

## 1. identical sexes reproduce the one-sex fit
for (om in c("dm", "mn")) {
  f1 <- fiticc(lfd, stk, sel_fun = sf, catch_by_gear = cb, settings = list(obs_model = om), verbose = FALSE)
  f2 <- fiticc(lfd, stk2, sel_fun = sf, catch_by_gear = cb, settings = list(obs_model = om), verbose = FALSE)
  stopifnot(f2$tmb_data$nsex == 2, abs(f1$opt$objective - f2$opt$objective) < 1e-6,
            max(abs(f1$obj$report()$spr_y - f2$obj$report()$spr_y)) < 1e-6)
}
# sexed halves (joint F + M), multinomial: same estimates (likelihood shifted by a constant)
lfds <- FLQuants(lapply(lfd, sexq))
f3 <- fiticc(lfds, stk2, sel_fun = sf, catch_by_gear = cb, settings = list(obs_model = "mn"), verbose = FALSE)
stopifnot(all(f3$tmb_data$sex_type == 3L), max(abs(f1$opt$par - f3$opt$par)) < 1e-5)
# per-recruit: same Fspr and YPR per total recruit; female SBPR = sex_ratio x one-sex
b1 <- FLicc:::brp_tmb_flicc(f1, Fseq = c(0.05, 0.2), spr = c(30, 40))
b2 <- FLicc:::brp_tmb_flicc(f2, Fseq = c(0.05, 0.2), spr = c(30, 40))
stopifnot(max(abs(b1$Fspr - b2$Fspr)) < 1e-8, max(abs(b2$YPR / b1$YPR - 1)) < 1e-8,
          max(abs(b2$SBPR / b1$SBPR - 0.5)) < 1e-8)
# equilibrium objects: catch both sexes, SSB female (identical sexes: catch equal, SSB x 0.5)
e1 <- eqstklen(f1, s = 0.75); e2 <- eqstklen(f2, s = 0.75)
stopifnot(max(abs(an(catch(e2)) / an(catch(e1)) - 1), na.rm = TRUE) < 1e-5,
          max(abs(an(ssbl(e2)) / an(ssbl(e1)) - 0.5), na.rm = TRUE) < 1e-6,
          max(abs(attr(e2, "refpts")[c("Fcur", "Fspr", "Ycur", "SSBcur", "SPRcur")] -
                  attr(e1, "refpts")[c("Fcur", "Fspr", "Ycur", "SSBcur", "SPRcur")])) < 1e-5)
# contour engine with two sexes equals the TMB per-recruit block
eng <- FLicc:::pr_engine_flicc(f2); bb <- FLicc:::brp_tmb_flicc(f2, Fseq = 0.1)
ee <- eng(0.1 * as.numeric(f2$obj$report(f2$obj$env$last.par.best)$sA))
stopifnot(abs(ee["SPR"] / (bb$SBPR / bb$SBPR0) - 1) < 1e-6, abs(ee["YPR"] / bb$YPR - 1) < 1e-6)
# EPR, numbers by sex and the fitted-years stock object (identical sexes)
stopifnot(max(abs(f1$obj$report()$epr_y - f2$obj$report()$epr_y)) < 1e-6,
          identical(dimnames(f2$report$N_sex)$unit, c("F", "M")))
k1 <- flicc_stklen(f1); k2 <- flicc_stklen(f2)
stopifnot(max(abs(an(catch(k2)) / an(catch(k1)) - 1)) < 1e-6)
stopifnot(inherits(plot_m(f2), "ggplot"), inherits(plot_mat(stk2), "ggplot"), inherits(plot_lw(f2), "ggplot"))
stopifnot(all(c("Linf_M", "Mk_M") %in% names(fits_summary(list(a = f2)))))
cat("identity ok\n")

## 2. composition types: females only, males only, combined + sexed in one gear-year
mixed <- FLQuants(Trawl = sexq(lfd$Trawl, 0.6, c("F", "U")),      # partly sexed: females + unsexed
                  Gillnet = sexq(lfd$Gillnet, 0, c("M")))          # males only
fm <- fiticc(mixed, stk2, sel_fun = sf, catch_by_gear = cb, settings = list(obs_model = "dm"), verbose = FALSE)
stopifnot(all(fm$tmb_data$sex_type[, 1] == 1L), all(fm$tmb_data$sex_type[, 2] == 2L),
          sum(fm$tmb_data$obs[, , 1]) > 0, sum(fm$tmb_data$obs[, , 2]) == 0,
          is.finite(fm$opt$objective))
cat("composition types ok\n")
err <- try(fiticc(lfds, stk, sel_fun = sf, catch_by_gear = cb, verbose = FALSE), silent = TRUE)
stopifnot(inherits(err, "try-error"), grepl("stocklen_sex", err))

## 3. recovery of a dimorphic stock (simulated from the model)
set.seed(42)
lens <- seq(2, 80, by = 2); yrs <- 2016:2025
mkq <- function(u = "unique") FLQuant(1, dimnames = list(len = lens, year = yrs, unit = u))
lfd0 <- FLQuants(trawl = mkq())
lhF <- FLPar(linf = 60, k = 0.15, M = 0.3, L50 = 35, a = 1e-5, b = 3)
lhM <- FLPar(linf = 45, k = 0.25, M = 0.3, L50 = 30, a = 1e-5, b = 3)
st2 <- stocklen_sex(lfd0, list(F = lhF, M = lhM), m_model = "constant")
td <- FLicc:::data_tmb_flicc(lfd0, st2$F, "logistic", FLQuants(trawl = quantSums(lfd0[[1]])),
                             settings = list(CVL = 0.1, GL = 30, pop_model = "gtg", obs_model = "mn",
                                             ngtg = 13, maxsd = 3, sex_ratio = 0.5), stklen_m = st2$M)
td$prior_sigmaF_mean <- log(0.3); td$prior_sigmaF_sd <- 0; td$prior_sigmaF_use <- 0L
td$FM_min <- 0.05; td$FM_max <- 4; td$FMpen_sd <- 0
pars <- FLicc:::init_tmb_flicc(td)
pars$log_Fk[] <- log(seq(1.5, 3, length.out = length(yrs))); pars$Sm <- c(log(25 / 60), log(4 / 60))
drop <- c("gear_names", "year_names", "Linf_init", "Linf_m_init", "Mk_m_init", "CVL_m_init", "Mk_init", "CVL_init", "Galpha_init")
rp <- TMB::MakeADFun(td[setdiff(names(td), drop)], pars, DLL = "FLicc", silent = TRUE)$report()
ps <- rp$plen_sex
sexed <- FLQuants(trawl = mkq(c("F", "M"))); uns <- FLQuants(trawl = mkq("U"))
for (iy in seq_along(yrs)) {
  pj <- c(ps[, iy, 1, 1], ps[, iy, 1, 2]); x <- rmultinom(1, 500, pj)[, 1]
  sexed$trawl[, iy, "F"] <- x[seq_along(lens)]; sexed$trawl[, iy, "M"] <- x[length(lens) + seq_along(lens)]
  uns$trawl[, iy, "U"] <- rmultinom(1, 500, ps[, iy, 1, 1] + ps[, iy, 1, 2])[, 1]
}
st <- list(obs_model = "mn")
stw <- stocklen_sex(lfd0, list(F = lhF, M = FLPar(linf = 52, k = 0.25, M = 0.3, L50 = 30, a = 1e-5, b = 3)),
                    m_model = "constant")
fa <- fiticc(sexed, stw, sel_fun = "logistic", catch_by_gear = 1, settings = modifyList(st, list(linf_m.sd = 0.5)), verbose = FALSE)
fb <- fiticc(uns, st2, sel_fun = "logistic", catch_by_gear = 1, settings = st, verbose = FALSE)
fc <- fiticc(FLQuants(trawl = uns$trawl[, , "U"]), stocklen(lfd0, lhF, m_model = "constant"),
             sel_fun = "logistic", catch_by_gear = 1, settings = st, verbose = FALSE)
bias <- function(f) median(f$obj$report()$spr_y / rp$spr_y - 1)
cat(sprintf("male Linf %.1f (true 45, start 52); SPR bias: sexed %+.2f, unsexed 2-sex %+.2f, 1-sex %+.2f\n",
            fa$obj$report()$Linf_m, bias(fa), bias(fb), bias(fc)))
stopifnot(abs(fa$obj$report()$Linf_m - 45) < 2.5, abs(bias(fa)) < 0.1, abs(bias(fb)) < 0.1, bias(fc) < -0.15)

## 4. survey gears: selectivity only (no F, no catch share), sexed survey + unsexed fishery
lfd0 <- FLQuants(trawl = mkq(), survey = mkq())
st2 <- stocklen_sex(lfd0, list(F = lhF, M = lhM), m_model = "constant")
set0 <- list(CVL = 0.1, GL = 30, pop_model = "gtg", obs_model = "mn", ngtg = 13, maxsd = 3,
             sex_ratio = 0.5, surveys = "survey")
td <- FLicc:::data_tmb_flicc(lfd0, st2$F, c("logistic", "logistic"),
                             FLQuants(trawl = quantSums(lfd0[[1]]), survey = quantSums(lfd0[[1]]) * 0),
                             settings = set0, stklen_m = st2$M)
stopifnot(identical(td$gear_type, c(0L, 1L)))
td$prior_sigmaF_mean <- log(0.3); td$prior_sigmaF_sd <- 0; td$prior_sigmaF_use <- 0L
td$FM_min <- 0.05; td$FM_max <- 4; td$FMpen_sd <- 0
pars <- FLicc:::init_tmb_flicc(td)
pars$log_Fk[, 1] <- log(seq(1.5, 3, length.out = length(yrs)))
pars$Sm <- c(log(30 / 60), log(4 / 60), log(12 / 60), log(3 / 60))
rp <- TMB::MakeADFun(td[setdiff(names(td), drop)], pars, DLL = "FLicc", silent = TRUE)$report()
stopifnot(all(rp$catch_share[, 2] == 0))
set.seed(7)
lfdv <- FLQuants(trawl = mkq("U"), survey = mkq(c("F", "M")))
for (iy in seq_along(yrs)) {
  lfdv$trawl[, iy, "U"] <- rmultinom(1, 500, rp$plen[, iy, 1])[, 1]
  x <- rmultinom(1, 300, c(rp$plen_sex[, iy, 2, 1], rp$plen_sex[, iy, 2, 2]))[, 1]
  lfdv$survey[, iy, "F"] <- x[seq_along(lens)]; lfdv$survey[, iy, "M"] <- x[length(lens) + seq_along(lens)]
}
fs <- fiticc(lfdv, stw, sel_fun = c("logistic", "logistic"), catch_by_gear = 1,
             settings = list(obs_model = "mn", surveys = "survey", linf_m.sd = 0.5), verbose = FALSE)
rs <- fs$obj$report()
cat(sprintf("survey: male Linf %.1f, SPR bias %+.2f\n", rs$Linf_m, median(rs$spr_y / rp$spr_y - 1)))
stopifnot(fs$opt$convergence == 0, abs(rs$Linf_m - 45) < 2.5, abs(median(rs$spr_y / rp$spr_y - 1)) < 0.1,
          all(rs$catch_share[, 2] == 0), all(rs$Fk[, 2] < 1e-6))
# life history by sex: females fixed, male Linf and M/k estimated from a wrong start
stw2 <- stocklen_sex(lfd0, list(F = lhF, M = FLPar(linf = 52, k = 0.25, M = 0.4, L50 = 30, a = 1e-5, b = 3)),
                     m_model = "constant")
fl <- fiticc(lfdv, stw2, sel_fun = c("logistic", "logistic"), catch_by_gear = 1,
             settings = list(obs_model = "mn", surveys = "survey", linf_m.sd = 0.5, Mk_m.sd = 0.5), verbose = FALSE)
rl <- fl$obj$report()
cat(sprintf("male Linf %.1f, male M/k %.2f (true 45, 1.2); female fixed Linf %.1f\n", rl$Linf_m, rl$Mk_m, rl$Linf))
stopifnot(abs(rl$Linf_m - 45) < 3, abs(rl$Mk_m - 1.2) < 0.2, abs(rl$Linf - 60) < 1e-8, abs(rl$Mk - 2) < 1e-8)
# diagnostics by sex: residuals, mean length, predictive check and plot_len include the sexed parts
r <- resid_flicc(fs, "osa")
stopifnot(all(c("U", "F", "M") %in% r$sex), all(r$sex[r$gear == "survey"] %in% c("F", "M")))
stopifnot(inherits(plot_resid_bubble(r), "ggplot"), inherits(plot_resid_len(r, pool_years = TRUE), "ggplot"))
lm <- lmean_resid_flicc(fs); stopifnot("survey F+M" %in% as.character(lm$table$gear))
pc <- ppc_flicc(fs, nsim = 50)
stopifnot(all(c("trawl", "survey F", "survey M") %in% pc$gear_names), inherits(plot_ppc_lfd(pc, pool_years = TRUE), "ggplot"))
pl <- plot_len(fs, by_gear = TRUE, sex = "panels")
stopifnot(all(c("survey F", "survey M") %in% levels(pl$data$gear)), !is.null(fs$report$obslen_sex))
# one-sex fits keep their residual columns (no sex column)
stopifnot(is.null(resid_flicc(f1)$sex))
# dimorphic: effective weight and maturity reproduce yield and female SBPR
bd <- FLicc:::brp_tmb_flicc(fs, Fseq = c(0.1, 0.3))
stopifnot(max(abs(colSums(bd$C * bd$sex$w_eff) / bd$YPR - 1)) < 1e-8,
          max(abs(colSums(bd$N * bd$sex$mat_eff * bd$sex$w_eff) / bd$SBPR - 1)) < 1e-8)
ed <- FLicc:::pr_engine_flicc(fs)(0.3 * as.numeric(fs$obj$report(fs$obj$env$last.par.best)$sA))
stopifnot(abs(ed["YPR"] / bd$YPR[2] - 1) < 1e-6)
stopifnot(identical(FLicc:::gear_state_flicc(fs)$gears, "trawl"))   # surveys not in contours
# Monte Carlo with two sexes; growth plots by sex
mcs <- mc_flicc(fs, nsim = 3, verbose = FALSE)
stopifnot(sum(mcs$conv == "ok") >= 2)
stopifnot(inherits(plot_growth(fs, sex = "M"), "ggplot"))
# indicators by sex, catch curves and pooling with sexed samples
lb <- lbi_flicc(fs)
stopifnot(all(c("trawl", "survey F", "survey M") %in% levels(lb$values$gear)),
          abs(lb$linf[["survey M"]] - as.numeric(fs$report$lhpar_m["linf"])) < 1e-6)
stopifnot("pooled" %in% levels(lbi_flicc(fs, sex = "total", pool = TRUE)$values$gear))
stopifnot(inherits(combine_lfd(lfdv, c(1, 1), lhpar = lhF, ess.n = NULL), "FLQuants"))
stopifnot(inherits(LBIspr(fs, spr = 40, thresh = 0.75), "FLIndices"))
p <- plot_sexratio(fs)
stopifnot(inherits(p, "ggplot"), all(attr(p, "data")$gear == "survey"))
# MSE estimator: flicc.sa() with sexed OEM samples and life history by sex
# (mse::track<- replaced by a stand-in so the test runs without mse)
`track<-` <- function(object, quant, year, value) object
lfd2 <- FLQuants(lapply(lfdv, function(q) propagate(q, 2)))
stk_oem <- FLQuant(0, dimnames = list(year = yrs)); attr(stk_oem, "lfd") <- lfd2
attr(stk_oem, "catch_gear") <- FLQuants(trawl = propagate(FLQuant(1, dimnames = list(year = yrs)), 2),
                                        survey = propagate(FLQuant(0, dimnames = list(year = yrs)), 2))
sa <- flicc.sa(stk_oem, args = list(it = 2, ay = 2026, dy = 2025), tracking = NULL,
               lhpar = list(F = lhF, M = lhM), sel_fun = c("logistic", "logistic"),
               settings = list(CVL = 0.1, GL = 30, pop_model = "gtg", obs_model = "mn", ngtg = 13, maxsd = 3,
                               surveys = "survey", sex_ratio = 0.5), lbi = "spr", ess.g = NULL)
stopifnot(all(is.finite(an(sa$ind$SPR))), abs(median(an(sa$ind$SPR[, ac(2025)])) / rp$spr_y[length(yrs)] - 1) < 0.15)
cat("sex tests passed\n")


# mirrored compositions by sex (females above, males below the line)
pm <- plot_len(fs)
stopifnot(inherits(pm, "ggplot"), all(c("F", "M") %in% as.character(pm$data$sex)))
jm <- pm$data[pm$data$sex %in% c("F", "M"), ]
sj <- tapply(jm$observed, paste(jm$panel, jm$year), sum)
stopifnot(all(abs(sj - 1) < 1e-8))            # joint composition sums to 1 over both sexes
stopifnot(all(pm$data$obs_s[pm$data$sex == "M"] <= 0))
stopifnot(inherits(plot_len(fs, sex = "panels", by_gear = TRUE), "ggplot"))
data(sexsim); pl2 <- plot_lfd(lfd_sexsim, type = "prop")
stopifnot(inherits(pl2, "ggplot"), "U" %in% as.character(pl2$data$sex))
data(alfonsino); stopifnot(inherits(plot_lfd(lfd_alfonsino, type = "prop"), "ggplot"))
cat("mirror plots OK\n")

# life history by sex from a one-sex FLPar plus the male values that differ
lhF <- FLPar(lhpar_sexsim@.Data[, "F", 1], params = dimnames(lhpar_sexsim)$params)
stopifnot(isTRUE(all.equal(c(lhpar_sex(lhF, linf = 45, k = 0.25, L50 = 28)@.Data), c(lhpar_sexsim@.Data))))
s1 <- stocklen_sex(lfd_sexsim, lhF, male = list(linf = 45, k = 0.25, L50 = 28), m_model = "constant")
s2 <- stocklen_sex(lfd_sexsim, lhpar_sexsim, m_model = "constant")
stopifnot(isTRUE(all.equal(s1$M, s2$M)), isTRUE(all.equal(s1$F, s2$F)))
stopifnot(dim(lhpar_sex(propagate(lhF, 3), linf = c(40, 45, 50)))[3] == 3)
cat("lhpar_sex OK\n")
