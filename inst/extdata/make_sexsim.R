# Simulated sex-dimorphic stock for the guide and tests (data(sexsim)).
# Females grow larger and slower than males; one trawl fishery with unsexed
# length samples and one survey with sexed samples (females and males in the
# same sample, from 2016). Simulated from the FLicc sex-structured model itself.
library(FLicc)
set.seed(2026)
lens <- seq(2, 80, by = 2); yrs <- 2011:2025
mkq <- function(u = "unique") FLQuant(1, dimnames = list(len = lens, year = yrs, unit = u))
lfd0 <- FLQuants(Trawl = mkq(), Survey = mkq())
lhF <- FLPar(linf = 60, k = 0.15, M = 0.3, L50 = 35, a = 1e-5, b = 3)
lhM <- FLPar(linf = 45, k = 0.25, M = 0.3, L50 = 28, a = 1e-5, b = 3)
st <- stocklen_sex(lfd0, list(F = lhF, M = lhM), m_model = "constant")
td <- FLicc:::data_tmb_flicc(lfd0, st$F, c("dsnormal", "logistic"),
                             FLQuants(Trawl = quantSums(lfd0[[1]]), Survey = quantSums(lfd0[[1]]) * 0),
                             settings = list(CVL = 0.1, GL = 30, pop_model = "gtg", obs_model = "mn",
                                             ngtg = 13, maxsd = 3, sex_ratio = 0.5, surveys = "Survey"),
                             stklen_m = st$M)
td$prior_sigmaF_mean <- log(0.3); td$prior_sigmaF_sd <- 0; td$prior_sigmaF_use <- 0L
td$FM_min <- 0.05; td$FM_max <- 4; td$FMpen_sd <- 0
pars <- FLicc:::init_tmb_flicc(td)
Fk <- 1 + 1.5 * plogis(seq(-3, 3, length.out = length(yrs)))            # F/k_female rising
pars$log_Fk[, 1] <- log(Fk)
# trawl: dome (mode 34, left SD 6, right SD 20); survey: logistic (SL50 12, SL95 16)
pars$Sm <- c(log(34 / 60), log(6 / 60), log(20 / 60), log(12 / 60), log(4 / 60))
drop <- c("gear_names", "year_names", "Linf_init", "Linf_m_init", "Mk_m_init", "CVL_m_init",
          "Mk_init", "CVL_init", "Galpha_init")
rp <- TMB::MakeADFun(td[setdiff(names(td), drop)], pars, DLL = "FLicc", silent = TRUE)$report()
tr <- mkq("U"); sv <- mkq(c("F", "M")); sv[] <- 0
for (iy in seq_along(yrs)) {
  tr[, iy, "U"] <- rmultinom(1, 400, rp$plen[, iy, 1])[, 1]
  if (yrs[iy] >= 2016) {
    x <- rmultinom(1, 250, c(rp$plen_sex[, iy, 2, 1], rp$plen_sex[, iy, 2, 2]))[, 1]
    sv[, iy, "F"] <- x[seq_along(lens)]; sv[, iy, "M"] <- x[length(lens) + seq_along(lens)]
  }
}
lfd_sexsim <- FLQuants(Trawl = tr, Survey = sv)
lhpar_sexsim <- FLPar(array(c(as.numeric(lhF), as.numeric(lhM)), dim = c(6, 2, 1),
                            dimnames = list(params = dimnames(lhF)$params, unit = c("F", "M"), iter = 1)))
truth_sexsim <- data.frame(year = yrs, spr = rp$spr_y, F_k = Fk, F = Fk * 0.15)
save(lfd_sexsim, lhpar_sexsim, truth_sexsim, file = "data/sexsim.rda", compress = "xz")
