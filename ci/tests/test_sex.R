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
drop <- c("gear_names", "year_names", "Linf_init", "Linf_m_init", "Mk_init", "CVL_init", "Galpha_init")
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
cat("sex tests passed\n")
