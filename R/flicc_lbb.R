#' LBB priors from length-frequency data
#'
#' Reproduces the start values and priors of LBB (Froese et al. 2018,
#' `LBB_33a.R`) from the length data: the length frequencies of all years are
#' standardised per year, weighted by the square root of the annual sample size
#' and summed; the selection ogive is approximated from the lengths at 10% and
#' 90% of the peak; Linf and Z/K are estimated by `nls()` from the fully
#' selected part of the aggregated frequencies (Froese et al. 2018, Eq. 4),
#' starting from the median of the annual largest lengths. With `Linf.user`
#' only Z/K is estimated.
#'
#' Prior SDs as in LBB: Linf with a CV of at most 1%; M/K = 1.5 with SD 0.15
#' (0.075 when `MK.user` is given).
#'
#' @param lfd An `FLQuants` (gears are pooled by summing) or an `FLQuant` of
#'   length frequencies; length dimnames are the lower bounds of the bins.
#' @param Linf.user,MK.user,Lc.user,Lstart.user Optional user values as in the
#'   LBB ID file (same length units as `lfd`).
#' @return A list with Linf, Linf.sd, ZK, ZK.sd, MK, MK.sd, FK, Lc, alpha
#'   (logistic slope per length unit), Lstart, Lmax, Lmax.med and the
#'   aggregated frequencies `LF.all`.
#' @references Froese, R., Winker, H., Coro, G., Demirel, N., Tsikliras, A.C.,
#'   Dimarchopoulou, D., Scarcella, G., Probst, W.N., Dureuil, M. and Pauly, D.
#'   (2018). A new approach for estimating stock status from length frequency
#'   data. ICES Journal of Marine Science 75(6): 2004-2015.
#' @export
lbb_priors_flicc <- function(lfd, Linf.user = NULL, MK.user = NULL, Lc.user = NULL,
                             Lstart.user = NULL) {
  x <- lbb_pool_flicc(lfd)
  L <- as.numeric(dimnames(x)$len)
  m <- matrix(x@.Data[, , 1, 1, 1, 1], nrow = length(L))
  m[is.na(m)] <- 0
  m <- m[, colSums(m) > 0, drop = FALSE]
  if (ncol(m) == 0) stop("lfd has no observations.")
  Lmax <- max(L[rowSums(m) > 0])
  Lmax.med <- stats::median(apply(m, 2, function(v) max(L[v > 0])))
  # aggregate: per year standardised and weighted by sqrt(n) (LBB function AG)
  agg <- rowSums(sweep(m, 2, sqrt(colSums(m)) / colSums(m), "*"))
  agg <- agg / max(agg)
  i1 <- which(agg > 0)[1]; i2 <- max(which(agg > 0))
  LF <- data.frame(Length = L[i1:i2], Freq = agg[i1:i2])

  L10 <- LF$Length[which(LF$Freq > 0.1)[1]]
  L90 <- LF$Length[which(LF$Freq > 0.9)[1]]
  Lc.st <- if (is.null(Lc.user)) (L10 + L90) / 2 else Lc.user
  alpha.st <- -log(1 / LF$Freq[which(LF$Freq > 0.1)[1]]) / (L10 - Lc.st)
  if (!is.finite(alpha.st) || alpha.st <= 0) alpha.st <- 4 / max(L90 - L10, diff(L)[1])
  Linf.st <- if (is.null(Linf.user)) Lmax.med else Linf.user
  sel <- LF$Length >= Lc.st
  Lmean.st <- sum(LF$Length[sel] * LF$Freq[sel]) / sum(LF$Freq[sel])
  MK.st <- if (is.null(MK.user)) 1.5 else MK.user
  ZK.st <- (Linf.st - Lmean.st) / (Lmean.st - Lc.st)
  if (!is.finite(ZK.st) || ZK.st <= 0) ZK.st <- 2 * MK.st
  FK.st <- if (ZK.st - MK.st > 0) ZK.st - MK.st else 0.3

  if (!is.null(Lstart.user)) {
    Lstart <- Lstart.user
  } else {
    Lstart <- (alpha.st * Lc.st - log(1 / 0.95 - 1)) / alpha.st
    Lstart.i <- which(LF$Length >= Lstart)[1]
    if (is.na(Lstart.i)) Lstart.i <- which.max(LF$Freq) + 1
    peak.i <- which.max(LF$Freq)
    if (Lstart.i < (peak.i + 1)) Lstart <- LF$Length[min(peak.i + 1, nrow(LF))]
    Lstart.i <- which(LF$Length >= Lstart)[1]
    if ((nrow(LF) - Lstart.i) < 4) Lstart <- LF$Length[max(Lstart.i - 1, 1)]
  }
  keep <- LF$Length >= Lstart & LF$Length < Linf.st
  L.L <- LF$Length[keep]; L.F <- LF$Freq[keep]
  if (length(L.L) < 4) stop("Too few fully selected length classes for the LBB Linf/ZK guess: set Lstart.user.")
  L.F <- L.F / sum(L.F)

  Linf <- Linf.st; Linf.sd <- 0.01 * Linf.st; ZK <- ZK.st; ZK.sd <- NA_real_
  if (is.null(Linf.user)) {
    mod <- tryCatch(stats::nls(L.F ~ ((Linf - L.L) / (Linf - Lstart))^ZK / sum(((Linf - L.L) / (Linf - Lstart))^ZK),
                               start = list(ZK = ZK.st, Linf = Linf.st),
                               lower = c(0.5 * ZK.st, 0.999 * Linf.st), upper = c(1.5 * ZK.st, 1.2 * Linf.st),
                               algorithm = "port"), error = function(e) NULL)
    if (is.null(mod)) {
      message("LBB guess: nls for Linf and Z/K failed; using the median of the annual largest lengths (", Linf.st, ") as Linf.")
    } else {
      cf <- summary(mod)$coefficients
      ZK <- cf["ZK", 1]; ZK.sd <- cf["ZK", 2]; Linf <- cf["Linf", 1]; Linf.sd <- cf["Linf", 2]
    }
  } else {
    mod <- tryCatch(stats::nls(L.F ~ exp(ZK * (log(1 - L.L / Linf.user) - log(1 - L.L[1] / Linf.user))) /
                                 sum(exp(ZK * (log(1 - L.L / Linf.user) - log(1 - L.L[1] / Linf.user)))),
                               start = list(ZK = ZK.st), lower = 0.7 * ZK.st, upper = 1.3 * ZK.st,
                               algorithm = "port"), error = function(e) NULL)
    if (!is.null(mod)) { cf <- summary(mod)$coefficients; ZK <- cf["ZK", 1]; ZK.sd <- cf["ZK", 2] }
  }
  Linf.sd <- if (is.finite(Linf.sd) && Linf.sd / Linf < 0.01) Linf.sd else 0.01 * Linf
  list(Linf = Linf, Linf.sd = Linf.sd, ZK = ZK, ZK.sd = ZK.sd,
       MK = MK.st, MK.sd = if (is.null(MK.user)) 0.15 else 0.075,
       FK = if (ZK - MK.st > 0) ZK - MK.st else 0.3,
       Lc = if (is.null(Lc.user)) 1.02 * Lc.st else Lc.user, alpha = alpha.st,
       Lstart = Lstart, Lmax = Lmax, Lmax.med = Lmax.med, LF.all = LF,
       user = list(Linf = !is.null(Linf.user), MK = !is.null(MK.user)))
}

#' Pool gears of a length-frequency object (internal)
#' @noRd
lbb_pool_flicc <- function(lfd) {
  if (inherits(lfd, "FLQuant")) return(lfd)
  if (!inherits(lfd, "FLQuants")) stop("lfd must be an FLQuants or FLQuant.")
  if (length(lfd) == 1) return(lfd[[1]])
  out <- lfd[[1]]
  for (i in 2:length(lfd)) out <- out + lfd[[i]]
  out
}

#' Fit LBB as a special case of FLicc
#'
#' Fits the LBB model (Froese et al. 2018) within FLicc: one growth type
#' (`pop_model = "gtg"`, `ngtg = 1`), one gear (gears in `lfd` are pooled, as
#' LBB pools the selected gears), logistic (trawl-like) or normal (gillnet-like)
#' selectivity, Linf and M/K estimated under LBB's priors from
#' [lbb_priors_flicc()] (or user values), lengths at or above the Linf prior
#' excluded, and each year fitted separately (`by_year = TRUE`, as LBB) or all
#' years jointly with a random walk on F (`by_year = FALSE`).
#'
#' Differences from LBB, by design: penalised maximum likelihood with standard
#' errors instead of MCMC; the composition likelihood is FLicc's (default
#' Dirichlet-multinomial with estimated precision, instead of LBB's Dirichlet
#' with a fixed concentration); selectivity parameters and F/K are estimated
#' without priors (LBB uses informative priors on Lc and alpha and a wide prior
#' on F/K); numbers in a length class are integrated over the class (LBB
#' `Pile = 1`). Constraints are explicit settings that can be changed.
#'
#' @param lfd `FLQuants` (or `FLQuant`) of length frequencies; gears are pooled.
#' @param Linf,MK Optional user values for Linf and M/K (LBB's `Linf.user`,
#'   `MK.user`); with `MK` the prior SD is 0.075 instead of 0.15.
#' @param Lc,Lstart Optional user values for the LBB guess (`Lc.user`,
#'   `Lstart.user`).
#' @param Lcut Optional length below which data are ignored (`Lcut.user`).
#' @param Lm50 Length at 50% maturity (used for SPR and the share of mature
#'   fish; default 2/3 of the Linf prior, with a message).
#' @param a,b Length-weight parameters (default `a = 0.01`, `b = 3`, as LBB's
#'   cubic weight).
#' @param k Von Bertalanffy k (default 1, so F is reported as F/K).
#' @param sel_fun `"logistic"` (LBB default, trawl-like) or `"normal"`
#'   (LBB `GausSel`).
#' @param years Optional years to use (LBB `Years.user`).
#' @param ess Optional effective sample size per year. `NULL` (default): raw
#'   counts with the Dirichlet-multinomial precision estimated from the data
#'   (FLicc's recommended likelihood). A number: counts rescaled to this ESS
#'   with [lfdess()] and the DM precision fixed (`dm_theta = "fixed"`), for
#'   comparison with LBB. LBB's Dirichlet has a concentration of 10 x 1000 =
#'   10^4 on the proportions, i.e. a variance close to a multinomial with
#'   n = 10^4; a fixed DM at ESS n has n_eff about n/2, so `ess = 20000`
#'   matches LBB's weighting (`ess = 1000` weights the data 20 times less).
#' @param by_year `TRUE` (default, as LBB): each year fitted separately;
#'   `FALSE`: all years jointly with a random walk on F.
#' @param settings Further [fiticc()] settings; they override the LBB defaults
#'   (e.g. `obs_model`, `linf.sd`, `Mk.sd`).
#' @param ... Passed to [fiticc()].
#' @return For `by_year = TRUE` a list of class `"flicc_lbb"` with `fits` (one
#'   fit per year), `priors`, `table` ([lbbref_flicc()]) and the inputs; for
#'   `by_year = FALSE` the same with one joint fit in `fits$joint`.
#' @seealso [lbb_priors_flicc()], [lbbref_flicc()], [plot_lbb()]
#' @export
lbb_flicc <- function(lfd, Linf = NULL, MK = NULL, Lc = NULL, Lstart = NULL, Lcut = NULL,
                      Lm50 = NULL, a = 0.01, b = 3, k = 1,
                      sel_fun = c("logistic", "normal"), years = NULL, ess = NULL,
                      by_year = TRUE, settings = list(), ...) {
  sel_fun <- match.arg(sel_fun)
  x <- lbb_pool_flicc(lfd)
  x[is.na(x)] <- 0
  if (!is.null(years)) x <- x[, as.character(intersect(years, as.numeric(dimnames(x)$year)))]
  if (!is.null(Lcut)) x <- x[as.numeric(dimnames(x)$len) >= Lcut, ]
  pr <- lbb_priors_flicc(x, Linf.user = Linf, MK.user = MK, Lc.user = Lc, Lstart.user = Lstart)
  if (is.null(Lm50)) {
    Lm50 <- round(2 / 3 * pr$Linf, 2)
    message("lbb_flicc: Lm50 not given; SPR uses maturity at 2/3 of the Linf prior (", Lm50, ").")
  }
  lhpar <- FLCore::FLPar(linf = pr$Linf, k = k, M = pr$MK * k, L50 = Lm50, a = a, b = b)
  set0 <- list(pop_model = "gtg", ngtg = 1, CVL = 0.1, maxsd = 0,
               linf.sd = pr$Linf.sd / pr$Linf, Mk.sd = pr$MK.sd / pr$MK,
               tail_cut = pr$Linf, obs_model = "dm",
               dm_theta = if (is.null(ess)) "saturating" else "fixed")
  set <- utils::modifyList(set0, settings)
  gname <- if (inherits(lfd, "FLQuants") && length(lfd) == 1) names(lfd) else "pooled"
  yrs <- dimnames(x)$year
  has <- vapply(yrs, function(y) sum(x[, y], na.rm = TRUE) > 0, logical(1))
  yrs <- yrs[has]
  fit_one <- function(xx) {
    lq <- FLCore::FLQuants(stats::setNames(list(xx), gname))
    if (!is.null(ess)) lq <- lfdess(lq, ess.g = ess)
    stk <- stocklen(lq, lhpar, m_model = "constant")
    fiticc(lq, stk, sel_fun = sel_fun, catch_by_gear = 1, settings = set, ...)
  }
  fits <- if (isTRUE(by_year)) {
    stats::setNames(lapply(yrs, function(y) tryCatch(fit_one(x[, y]), error = function(e) {
      warning("LBB fit for year ", y, " failed: ", conditionMessage(e), call. = FALSE); NULL })), yrs)
  } else list(joint = fit_one(x[, yrs]))
  out <- structure(list(fits = fits, priors = pr, by_year = isTRUE(by_year), sel_fun = sel_fun,
                        settings = set, lhpar = lhpar, lfd = x, Lm50 = Lm50, ess = ess),
                   class = "flicc_lbb")
  out$table <- lbbref_flicc(out)
  out
}

#' LBB per-recruit B/B0 and relative yield (Beverton-Holt closed form)
#'
#' Port of the `BH()` function of `LBB_33a.R` (Froese et al. 2018, Eqs 12-19):
#' exploited biomass per recruit relative to the unfished state (B/B0, the
#' biomass vulnerable to the gear) and relative yield per recruit Y/R', from
#' Linf, M/K, F/K and the selectivity, with cubic weight; Y/R' is reduced in
#' proportion when B/B0 < 0.25.
#' @param lens Length classes (lower bounds) present in the data.
#' @param Linf,MK,FK Asymptotic length, M/K and F/K.
#' @param gaus `FALSE` (logistic) or `TRUE` (normal selectivity).
#' @param selpar1,selpar2 Lc/Linf and alpha (relative to Linf) for logistic;
#'   mean/Linf and SD/Linf for normal.
#' @return `c(BB0 = , YR = )`.
#' @export
lbb_bh_flicc <- function(lens, Linf, MK, FK, gaus = FALSE, selpar1, selpar2) {
  lens <- sort(unique(lens))
  Lx <- if (!gaus) lens[lens >= Linf * (selpar1 - 4.59 / selpar2)][1] else
    lens[lens >= Linf * (selpar1 - 3 * selpar2)][1]
  if (is.na(Lx)) Lx <- lens[1]
  cw <- stats::median(diff(lens))
  FM <- FK / MK
  L.bh <- seq(from = Lx, to = Linf, by = cw)
  r <- L.bh / Linf
  n <- length(r)
  if (n < 2) return(c(BB0 = NA_real_, YR = NA_real_))
  s1 <- function(v) if (!gaus) 1 / (1 + exp(-selpar2 * (v - selpar1))) else exp(-((v - selpar1)^2 / (2 * selpar2^2)))
  SL <- vapply(seq_len(n), function(o) if (o < n) mean(c(s1(r[o]), s1(r[o + 1]))) else s1(r[o]), numeric(1))
  yfun <- function(rr, sl) FM * sl / (1 + FM * sl) * (1 - rr)^MK *
    (1 - 3 * (1 - rr) / (1 + 1 / (MK + FK * sl)) + 3 * (1 - rr)^2 / (1 + 2 / (MK + FK * sl)) -
       (1 - rr)^3 / (1 + 3 / (MK + FK * sl)))
  bfun <- function(rr) (1 - rr)^MK * (1 - 3 * (1 - rr) / (1 + 1 / MK) + 3 * (1 - rr)^2 / (1 + 2 / MK) -
                                        (1 - rr)^3 / (1 + 3 / MK))
  rr <- numeric(n - 1); G <- numeric(n - 1)
  for (o in seq_len(n - 1)) {
    rr[o] <- (1 - r[o + 1])^(FK * SL[o]) / (1 - r[o])^(FK * SL[o])
    G[o] <- prod(rr[1:o])
  }
  YR <- numeric(n); B <- numeric(n)
  for (o in seq_len(n)) {
    YR[o] <- if (o == 1) yfun(r[o], SL[o]) - yfun(r[o + 1], SL[o]) * G[o] else
      if (o == n) yfun(r[o], SL[o]) * G[o - 1] else
        yfun(r[o], SL[o]) * G[o - 1] - yfun(r[o + 1], SL[o]) * G[o]
    B[o] <- if (o < n) (bfun(r[o]) - bfun(r[o + 1])) * SL[o] else bfun(r[o]) * SL[o]
  }
  BB0 <- sum(YR / FM) / sum(B)
  YRs <- sum(YR)
  if (BB0 < 0.25) YRs <- YRs * BB0 / 0.25
  c(BB0 = BB0, YR = YRs)
}

#' LBB reference points and status by year
#'
#' Table of LBB quantities by year from an [lbb_flicc()] result: Linf, Lc
#' (= SL50) and alpha (logistic) or mean and SD (normal), M/K, F/K, F/M, Z/K
#' with 95% intervals (from the fits' standard errors on the log scale), Lopt =
#' 3/(3 + M/K) Linf, Lc_opt = Linf (2 + 3 F/M) / ((1 + F/M)(3 + M/K)), Lmean
#' (above Lc), L95, the share of mature fish, B/B0 and Y/R' from
#' [lbb_bh_flicc()] (intervals by LBB's error propagation), FLicc's EPR (the
#' same quantity from the FLicc engine, see [plot_epr()]) and SPR.
#' Attributes: `BFM1B0` (B/B0 at F = M with Lc = Lc_opt, LBB's proxy for
#' Bmsy) and `YRFM1`.
#' @param x Output of [lbb_flicc()].
#' @return A `data.frame`.
#' @export
lbbref_flicc <- function(x) {
  gaus <- identical(x$sel_fun, "normal")
  L <- as.numeric(dimnames(x$lfd)$len)
  one <- function(f, yr) {
    p <- f$par
    g <- function(nm, i = 1) { r <- p[p$par == nm, , drop = FALSE]; if (nrow(r) < i) c(NA, NA) else c(r$mpd[i], r$se[i]) }
    ci <- function(v) c(exp(v[1]), exp(v[1] - 1.96 * v[2]), exp(v[1] + 1.96 * v[2]))
    linf <- ci(g("log_Linf")); mk <- ci(g("log_Mk"))
    if (anyNA(linf)) linf <- rep(as.numeric(f$report$lhpar["linf"]), 3)
    if (anyNA(mk)) mk <- rep(as.numeric(f$report$lhpar["M"]) / as.numeric(f$report$lhpar["k"]), 3)
    yi <- match(yr, as.character(f$tmb_data$year_names))
    fk <- ci(g("log_Fk", yi))
    s1 <- g("Sm", 1); s2 <- g("Sm", 2)
    sp1 <- exp(s1[1]) * linf[1]; sp2 <- exp(s2[1]) * linf[1]
    obs <- as.numeric(x$lfd[, yr])
    obs[is.na(obs)] <- 0
    keepL <- L < linf[1]
    if (!gaus) {
      Lc <- sp1; alpha <- log(19) * linf[1] / sp2       # relative alpha (per L/Linf)
      Lc_ci <- exp(s1[1] + c(-1.96, 1.96) * s1[2]) * linf[1]
      bh <- lbb_bh_flicc(L[obs > 0], linf[1], mk[1], fk[1], FALSE, Lc / linf[1], alpha)
    } else {
      Lc <- sp1; alpha <- NA; Lc_ci <- exp(s1[1] + c(-1.96, 1.96) * s1[2]) * linf[1]
      bh <- lbb_bh_flicc(L[obs > 0], linf[1], mk[1], fk[1], TRUE, sp1 / linf[1], sp2 / linf[1])
    }
    fm <- fk / mk[1]
    rel <- function(lo, mid, hi) c(lo = (mid - lo) / mid, hi = (hi - mid) / mid)
    rl <- sqrt(rel(fm[2], fm[1], fm[3])^2 + rel(mk[2], mk[1], mk[3])^2 +
                 rel(fk[2], fk[1], fk[3])^2 + rel(linf[2], linf[1], linf[3])^2)
    Lmid <- L + c(diff(L), utils::tail(diff(L), 1)) / 2
    sel <- L >= Lc & keepL
    Lmean <- if (sum(obs[sel]) > 0) sum(L[sel] * obs[sel]) / sum(obs[sel]) else NA
    cs <- cumsum(obs[keepL]) / sum(obs[keepL])
    L95 <- L[keepL][which(cs >= 0.95)[1]]
    data.frame(year = as.numeric(yr), Linf = linf[1], Linf_lcl = linf[2], Linf_ucl = linf[3],
               Lc = Lc, Lc_lcl = Lc_ci[1], Lc_ucl = Lc_ci[2], alpha = alpha,
               sel_sd = if (gaus) sp2 else NA_real_,
               MK = mk[1], MK_lcl = mk[2], MK_ucl = mk[3],
               FK = fk[1], FK_lcl = fk[2], FK_ucl = fk[3],
               FM = fm[1], FM_lcl = fm[2], FM_ucl = fm[3],
               ZK = mk[1] + fk[1],
               Lopt = 3 / (3 + mk[1]) * linf[1],
               Lc_opt = linf[1] * (2 + 3 * fm[1]) / ((1 + fm[1]) * (3 + mk[1])),
               Lmean = Lmean, L95 = L95,
               perc_mat = sum(obs[L > x$Lm50]) / sum(obs),
               BB0 = bh[["BB0"]], BB0_lcl = bh[["BB0"]] * (1 - rl[["lo"]]), BB0_ucl = bh[["BB0"]] * (1 + rl[["hi"]]),
               YR = bh[["YR"]],
               EPR = if (!is.null(f$report$epr)) as.numeric(f$report$epr[, yr]) else NA_real_,
               SPR = as.numeric(f$report$spr[, yr]))
  }
  rows <- if (isTRUE(x$by_year)) {
    do.call(rbind, lapply(names(x$fits), function(y) if (is.null(x$fits[[y]])) NULL else one(x$fits[[y]], y)))
  } else {
    f <- x$fits$joint
    do.call(rbind, lapply(as.character(f$tmb_data$year_names), function(y) one(f, y)))
  }
  if (is.null(rows) || nrow(rows) == 0) return(NULL)
  rows$Lmean_Lopt <- rows$Lmean / rows$Lopt
  rows$Lc_Lc_opt <- rows$Lc / rows$Lc_opt
  # LBB proxy for Bmsy: B/B0 at F = M with Lc = Lc_opt (medians over years)
  Linf.med <- stats::median(rows$Linf); MK.med <- stats::median(rows$MK)
  bfm <- if (!gaus) {
    lbb_bh_flicc(L, Linf.med, MK.med, MK.med, FALSE, 5 / (2 * (3 + MK.med)), stats::median(rows$alpha))
  } else {
    lbb_bh_flicc(L, Linf.med, MK.med, MK.med, TRUE, 3 / (3 + MK.med), stats::median(rows$sel_sd) / Linf.med)
  }
  attr(rows, "BFM1B0") <- bfm[["BB0"]]
  attr(rows, "YRFM1") <- bfm[["YR"]]
  rows
}

#' Plot an LBB analysis
#'
#' Four panels in the style of LBB's summary figure: (a) observed and predicted
#' length frequencies of the last year with Lopt and Linf; (b) Lc, Lc_opt,
#' Lmean and Lopt by year; (c) F/M by year with 95% intervals and F = M;
#' (d) B/B0 by year with intervals, B/B0 at F = M and Lc = Lc_opt (proxy for
#' Bmsy, dashed) and half of it (dotted). B/B0 in LBB is the exploited
#' (vulnerable) biomass relative to the unfished state, not spawning biomass;
#' FLicc's SPR is in `x$table$SPR`.
#' @param x Output of [lbb_flicc()].
#' @param year Year for panel (a) (default the last).
#' @param smooth Apply LBB's 3-year moving average to the time series.
#' @return A patchwork of ggplots.
#' @export
plot_lbb <- function(x, year = NULL, smooth = FALSE) {
  tab <- x$table
  if (is.null(tab)) stop("No LBB results to plot.")
  ma <- function(v) { out <- v; if (length(v) >= 3) for (i in seq_along(v)) out[i] <- mean(v[max(1, i - 2):i]); out }
  if (isTRUE(smooth)) for (nm in c("Lc", "Lmean", "Lopt", "Lc_opt", "FM", "FM_lcl", "FM_ucl", "BB0", "BB0_lcl", "BB0_ucl"))
    tab[[nm]] <- ma(tab[[nm]])
  yr <- if (is.null(year)) as.character(max(tab$year)) else as.character(year)
  f <- if (isTRUE(x$by_year)) x$fits[[yr]] else x$fits$joint
  ob <- as.data.frame(f$report$obslen[[1]][, yr])
  pr <- as.data.frame(f$report$predlen_gear[[1]][, yr])
  d1 <- data.frame(L = as.numeric(as.character(ob$len)), obs = ob$data / sum(ob$data, na.rm = TRUE),
                   pred = pr$data[match(ob$len, pr$len)] / sum(pr$data, na.rm = TRUE))
  tr <- tab[tab$year == as.numeric(yr), ]
  p1 <- ggplot2::ggplot(d1, ggplot2::aes(L)) +
    ggplot2::geom_col(ggplot2::aes(y = obs), fill = "grey75", width = 0.9 * stats::median(diff(d1$L))) +
    ggplot2::geom_line(ggplot2::aes(y = pred), colour = "#D55E00", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = c(tr$Lopt, tr$Linf), linetype = "dashed", colour = "darkgreen") +
    ggplot2::annotate("text", x = c(tr$Lopt, tr$Linf), y = max(d1$obs, na.rm = TRUE) * 1.05,
                      label = c("Lopt", "Linf"), colour = "darkgreen", size = 3) +
    ggplot2::labs(x = "Length", y = "Relative frequency", title = paste0("(a) ", yr, ": Z/K = ", round(tr$ZK, 2))) +
    ggplot2::theme_bw()
  d2 <- rbind(data.frame(year = tab$year, value = tab$Lc, what = "Lc"),
              data.frame(year = tab$year, value = tab$Lc_opt, what = "Lc_opt"),
              data.frame(year = tab$year, value = tab$Lmean, what = "Lmean"),
              data.frame(year = tab$year, value = tab$Lopt, what = "Lopt"))
  p2 <- ggplot2::ggplot(d2, ggplot2::aes(year, value, colour = what, linetype = what)) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 1.3) +
    ggplot2::scale_colour_manual(values = c(Lc = "#0072B2", Lc_opt = "#0072B2", Lmean = "#D55E00", Lopt = "#D55E00")) +
    ggplot2::scale_linetype_manual(values = c(Lc = "solid", Lc_opt = "dashed", Lmean = "solid", Lopt = "dashed")) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(x = "Year", y = "Length", colour = NULL, linetype = NULL, title = "(b) Lengths") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  p3 <- ggplot2::ggplot(tab, ggplot2::aes(year, FM)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = FM_lcl, ymax = FM_ucl), fill = "grey80") +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 1.3) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "darkgreen") +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(x = "Year", y = "F/M", title = "(c) F/M") + ggplot2::theme_bw()
  bfm <- attr(x$table, "BFM1B0")
  p4 <- ggplot2::ggplot(tab, ggplot2::aes(year, BB0)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(BB0_lcl, 0), ymax = BB0_ucl), fill = "grey80") +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 1.3) +
    ggplot2::geom_hline(yintercept = bfm, linetype = "dashed", colour = "darkgreen") +
    ggplot2::geom_hline(yintercept = bfm / 2, linetype = "dotted", colour = "red") +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(x = "Year", y = "B/B0", title = "(d) B/B0 (exploited biomass)") + ggplot2::theme_bw()
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
}

#' @export
print.flicc_lbb <- function(x, ...) {
  pr <- x$priors
  cat("LBB in FLicc (", if (x$by_year) "years fitted separately" else "joint fit", ", ", x$sel_fun, " selectivity)\n", sep = "")
  cat("Linf prior ", round(pr$Linf, 2), " (SD ", signif(pr$Linf.sd, 2), if (pr$user$Linf) ", user" else "",
      "), M/K prior ", pr$MK, " (SD ", pr$MK.sd, if (pr$user$MK) ", user" else "", "), Z/K guess ",
      round(pr$ZK, 2), ", Lmax ", pr$Lmax, ", median Lmax ", pr$Lmax.med, "\n", sep = "")
  if (!is.null(x$table)) {
    tb <- x$table
    print(data.frame(year = tb$year, Linf = round(tb$Linf, 1), Lc = round(tb$Lc, 1), MK = round(tb$MK, 2),
                     FM = round(tb$FM, 2), ZK = round(tb$ZK, 2), BB0 = round(tb$BB0, 2), SPR = round(tb$SPR, 2)),
          row.names = FALSE)
    cat("B/B0 at F = M, Lc = Lc_opt (Bmsy proxy):", round(attr(tb, "BFM1B0"), 2), "\n")
  }
  invisible(x)
}
