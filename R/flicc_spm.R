# FLicc-SPM: surplus production model in the FLicc objective (prototype)
#
# One TMB objective: FLicc length likelihood (gtg, one or two sexes; F as
# fixed effects, a tuning device for SPR) + Pella-Tomlinson state-space surplus production model in
# spawning biomass fitted to catch and CPUE (translated from JABBA3::jb3_nll,
# Selex state) + coupling. r and m come from the per-recruit equilibrium of the
# gtg model with the estimated selectivity (at Fmsy or Fspr, by taped Newton
# iterations), the EB/SB conversion from the same equilibrium on an F grid,
# and SPR by year is linked to SB/SB0 under steepness. See the design log in
# the FLicc - Development project (flicc_spm_log.md).

#' Settings of FLicc-SPM
#'
#' Default settings for the surplus production part of FLicc-SPM, passed as
#' \code{settings$spm} to \code{fiticc()}. With \code{settings$spm} the length
#' model is fitted jointly with a Pella-Tomlinson state-space surplus production
#' model in spawning biomass (as JABBA-Select / JABBA3 Selex) to catch and CPUE.
#'
#' @param catch total catch in weight, a numeric vector named by year (or an
#'   FLQuant); its years are the years of the production model and must include
#'   the length-composition years
#' @param index CPUE indices: a matrix (years x indices, row names years), a
#'   numeric vector named by year (one index), or an FLQuant / FLQuants
#' @param index.se optional log-scale standard errors, same shape as \code{index}
#' @param index.gear per index, the gear whose estimated selectivity defines the
#'   exploitable stock the index measures (a gear name), \code{"fishery"} (joint
#'   selectivity of the fishery, default) or \code{"SB"} (spawning biomass)
#' @param index.units per index, \code{"biomass"} (default) or \code{"numbers"}
#' @param h Beverton-Holt steepness (default \code{lhpar["s"]}, else 0.7)
#' @param shape source of r and m: \code{"fmsy"} (maximum equilibrium yield, default)
#'   or \code{"spr"} (Fspr at \code{spr.tgt})
#' @param spr.tgt SPR target for \code{shape = "spr"}
#' @param r.sd log-scale SD of the penalty log r ~ N(log r(theta), r.sd); 0 fixes r at r(theta)
#' @param m.sd as \code{r.sd} for the shape m (default 0: m = m(theta))
#' @param link \code{"spr"} (default): the depletion implied by the length
#'   model's SPR under steepness, log D_spr ~ N(log SB/SB0, sigma.link), links
#'   the two models; \code{"none"}: no link
#' @param sigma.link SD of log D_spr around log SB/SB0
#' @param link.rho fixed AR1 correlation of the link errors over the
#'   length-data years (default 0, independent errors). It cannot be estimated:
#'   with the link SD fixed, rho runs to 1
#' @param lag production lag in years
#' @param bkinit prior on SB_1/SB0: c(median, log sd, on)
#' @param sdb prior on the process error SD: c(median, log sd, on). With
#'   \code{re = FALSE} the SD is fixed at the median
#' @param re logical. \code{TRUE} (default): the SPM states (log SB_t/SB0) are
#'   random effects, integrated out by the Laplace approximation, and the process
#'   error SD is estimated. \code{FALSE}: penalised likelihood, as SS3 treats
#'   recruitment deviations; the states are fixed effects with the process error
#'   density as a penalty and the SD fixed at \code{sdb[1]} (an estimated SD
#'   collapses to zero without the Laplace approximation). About 10-20 times
#'   faster; for ensembles, sensitivity grids and MSE, with \code{re = TRUE} for
#'   final runs. Without \code{index.se} the additive index SD is then fixed
#'   too (\code{sdi}), since the free states would otherwise fit the index
#'   exactly and the SD would collapse
#' @param sdi additive index SD (log scale), one value or one per index.
#'   \code{NULL} (default): estimated (total SD sqrt(se^2 + sdi^2)), except with
#'   \code{re = FALSE} and no \code{index.se}, where it is fixed at 0.2. A value
#'   fixes it
#' @param K prior on K: \code{"half8"} (default, as JABBA3: flat below 5.7 x max
#'   catch, i.e. 8 x max catch as expected value with CV 1, half-normal on log K
#'   above, log sd 0.833), or c(median, log sd, on) in catch units (lognormal), or
#'   \code{NULL} (no prior). With two sexes K is the female SB0, and the half8 scale is
#'   multiplied by the female share of recruits (\code{sex_ratio})
#' @param xgrid F grid for the EB/SB curves, multiples of M at Linf (first 0)
#' @param newton.it taped Newton iterations for the reference F and m
#' @param P.bound soft bounds on SB/SB0
#' @param sel.weights gear weights of the per-recruit selectivity: \code{"catch"}
#'   (catch shares, default), \code{"F"} (F-weighted joint selectivity) or a vector
#' @param timing biomass per recruit for the SPM (EB/SB curves, Hmsy):
#'   \code{"start"} of the year (default, as the discrete SPM's harvest rate
#'   C_t / EB_t) or \code{"mean"} of the year (FLicc numbers)
#' @param frw random walk on the log harvest rate H_t = C_t / EB_t (SPiCT-like F
#'   process, as JABBA3 \code{frw}): \code{"none"} (default; catch-conditioned,
#'   as JABBA), \code{"index"} (increments from the first index year on) or
#'   \code{"all"} (all years). Its SD is estimated with the prior \code{sdf}
#'   (fixed at \code{sdf[1]} with \code{re = FALSE}). Use with
#'   \code{catch.error = TRUE}
#' @param catch.error logical: catch observed with lognormal error (catch
#'   deviations; random effects with \code{re = TRUE}), log SD \code{sdc}
#' @param sdc log SD of the catch error
#' @param sdf prior on the SD of the F random walk: c(median, log sd, on)
#' @return list of settings
#' @export
spm_settings_flicc <- function(catch = NULL, index = NULL, index.se = NULL, index.gear = "fishery",
                               index.units = "biomass", h = NULL, shape = c("fmsy", "spr"), spr.tgt = 0.4,
                               r.sd = 0.2, m.sd = 0, link = "spr", sigma.link = 0.2, link.rho = 0,
                               lag = 0, bkinit = c(0.9, 0.25, 1), sdb = c(0.15, 0.5, 1), K = "half8",
                               xgrid = c(0, exp(seq(log(0.02), log(40), length.out = 45))),
                               newton.it = 12, P.bound = c(0.02, 1.3), sel.weights = "catch",
                               timing = c("start", "mean"), re = TRUE, sdi = NULL,
                               frw = c("none", "index", "all"), catch.error = FALSE, sdc = 0.1, sdf = c(0.3, 1, 1)) {
  list(catch = catch, index = index, index.se = index.se, index.gear = index.gear, index.units = index.units,
       h = h, shape = match.arg(shape), spr.tgt = spr.tgt, r.sd = r.sd, m.sd = m.sd, link = link,
       sigma.link = sigma.link, link.rho = link.rho, lag = lag, bkinit = bkinit, sdb = sdb, K = K,
       xgrid = xgrid, newton.it = newton.it, P.bound = P.bound, sel.weights = sel.weights,
       timing = match.arg(timing), re = isTRUE(re), sdi = sdi,
       frw = match.arg(frw), catch.error = isTRUE(catch.error), sdc = as.numeric(sdc), sdf = as.numeric(sdf))
}

# gear weights of the per-recruit selectivity: "F" (0: F-weighted joint
# selectivity of the terminal years, default; slow only if F is made a random
# effect with settings$F_re = TRUE), "catch" (catch
# shares of the fishery gears, data) or a numeric vector by gear
.spm_selw <- function(w, td) {
  ng <- td$ngear
  if (is.null(w) || identical(w, "F")) return(0)
  if (identical(w, "catch")) {
    cw <- colMeans(as.matrix(td$catch_wt)); cw[td$gear_type == 1L] <- 0
    if (!all(is.finite(cw)) || sum(cw) <= 0) cw <- as.numeric(td$gear_type == 0L)
    return(as.numeric(cw / sum(cw)))
  }
  w <- as.numeric(w); if (length(w) != ng) stop("settings$spm$sel.weights: one weight per gear", call. = FALSE)
  w / sum(w)
}

# year-named numeric vector from a vector or FLQuant
.spm_series <- function(x, what) {
  if (inherits(x, "FLQuant")) return(stats::setNames(as.numeric(x), dimnames(x)$year))
  if (is.null(names(x))) stop("settings$spm$", what, " must be named by year", call. = FALSE)
  stats::setNames(as.numeric(x), names(x))
}
# years x indices matrix from a vector, matrix, FLQuant or FLQuants
.spm_index <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "FLQuants") || (is.list(x) && !is.data.frame(x))) {
    l <- lapply(x, .spm_series, what = "index")
    yrs <- sort(unique(unlist(lapply(l, names))))
    m <- sapply(l, function(v) v[yrs]); m <- matrix(m, length(yrs), dimnames = list(yrs, names(l)))
    return(m)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    m <- as.matrix(x)
    if (is.null(rownames(m))) stop("settings$spm$index: row names must be years", call. = FALSE)
    if (is.null(colnames(m))) colnames(m) <- paste0("index", seq_len(ncol(m)))
    return(m)
  }
  v <- .spm_series(x, "index")
  matrix(v, ncol = 1, dimnames = list(names(v), "index1"))
}

# dummy SPM data when FLicc-SPM is off (and for fits made before it)
spm_defaults_tmb_flicc <- function(td) {
  d <- list(spm_use = 0L, spm_nyear = 1L, spm_catch = 1, spm_lfd_t = rep(0L, td$nyear), spm_nidx = 0L,
            spm_logI = numeric(0), spm_obs_t = integer(0), spm_obs_i = integer(0), spm_se2 = numeric(0),
            spm_idx_gear = integer(0), spm_idx_units = integer(0), spm_h = 0.7, spm_shape = 0L,
            spm_spr_tgt = 0.4, spm_r_sd = 0, spm_m_sd = 0, spm_link_spr = 0L, spm_sigma_link = 1,
            spm_lagidx = 0L, spm_bk_pr = c(1, 1, 0),
            spm_sdb_pr = c(1, 1, 0), spm_K_pr = c(1, 1, 0), spm_xgrid = c(0, 1), spm_newton_it = 1L,
            spm_Pbound = c(0.02, 1.3), spm_selw = 0, spm_timing = 1L,
            spm_catch_err = 0L, spm_sdc = 0.1, spm_frw_t = integer(0), spm_sdf_pr = c(0.3, 1, 0))
  for (n in names(d)) if (is.null(td[[n]])) td[[n]] <- d[[n]]
  td
}
spm_par_defaults_flicc <- function(pl, td) {
  nidx <- if (is.null(td$spm_nidx)) 0L else td$spm_nidx
  d <- list(spm_logK = 0, spm_logr = 0, spm_logm = log(2), spm_lbk = 0, spm_logsdb = log(0.1),
            spm_logq = rep(0, max(nidx, 1L)), spm_logsdi = rep(log(0.2), max(nidx, 1L)),
            spm_logP = rep(0, if (is.null(td$spm_nyear)) 1L else td$spm_nyear), spm_trho = 0,
            spm_cdev = rep(0, if (is.null(td$spm_nyear)) 1L else td$spm_nyear), spm_logsdf = log(0.3))
  for (n in names(d)) if (is.null(pl[[n]])) pl[[n]] <- d[[n]]
  pl
}

# SPM data, parameters, map and random effects for fiticc_core()
spm_setup_flicc <- function(spm, tmb_data, parameters, map, random, stklen) {
  spm_names <- c("spm_logK", "spm_logr", "spm_logm", "spm_lbk", "spm_logsdb", "spm_logq", "spm_logsdi", "spm_logP", "spm_trho",
                 "spm_cdev", "spm_logsdf")
  if (is.null(spm)) {
    tmb_data <- spm_defaults_tmb_flicc(tmb_data)
    parameters <- spm_par_defaults_flicc(parameters, tmb_data)
    for (n in spm_names) map[[n]] <- factor(rep(NA, length(parameters[[n]])))
    return(list(tmb_data = tmb_data, parameters = parameters, map = map, random = random))
  }
  if (!isTRUE(tmb_data$pop_model == 2L)) stop("FLicc-SPM needs pop_model = \"gtg\"", call. = FALSE)
  sp <- do.call(spm_settings_flicc, spm[intersect(names(spm), names(formals(spm_settings_flicc)))])
  unk <- setdiff(names(spm), names(formals(spm_settings_flicc)))
  if (length(unk)) stop("Unknown settings$spm: ", paste(unk, collapse = ", "), call. = FALSE)
  if (is.null(sp$catch)) stop("settings$spm$catch is required", call. = FALSE)
  catch <- .spm_series(sp$catch, "catch")
  if (any(!is.finite(catch) | catch < 0)) stop("settings$spm$catch: missing or negative values", call. = FALSE)
  yrs <- as.numeric(names(catch))
  if (any(diff(yrs) != 1)) stop("settings$spm$catch: consecutive years needed", call. = FALSE)
  lfd_years <- as.numeric(tmb_data$year_names)
  lt <- match(lfd_years, yrs)
  if (anyNA(lt)) stop("length-composition years outside the catch years: ",
                      paste(lfd_years[is.na(lt)], collapse = ", "), call. = FALSE)
  I <- .spm_index(sp$index)
  nidx <- if (is.null(I)) 0L else ncol(I)
  Im <- matrix(NA_real_, length(yrs), max(nidx, 1)); SEm <- Im
  if (nidx > 0) {
    ri <- match(as.numeric(rownames(I)), yrs)
    if (anyNA(ri)) stop("index years outside the catch years", call. = FALSE)
    Im[ri, seq_len(nidx)] <- I
    if (!is.null(sp$index.se)) { S <- .spm_index(sp$index.se); SEm[match(as.numeric(rownames(S)), yrs), seq_len(nidx)] <- S }
  }
  Im[!is.finite(Im) | Im <= 0] <- NA
  ok <- if (nidx > 0) which(!is.na(Im[, seq_len(nidx), drop = FALSE]), arr.ind = TRUE) else matrix(0L, 0, 2)
  ok <- ok[order(ok[, 1], ok[, 2]), , drop = FALSE]
  se2 <- if (nrow(ok)) { s <- SEm[ok]; ifelse(is.finite(s), s^2 + 0.001^2, 0.1^2) } else numeric(0)
  gn <- tmb_data$gear_names
  ig <- rep_len(as.character(sp$index.gear), max(nidx, 1))[seq_len(nidx)]
  igc <- vapply(ig, function(g) if (g == "SB") -1L else if (g == "fishery") -2L else {
    k <- match(g, gn); if (is.na(k)) stop("settings$spm$index.gear: unknown gear ", g, call. = FALSE); k - 1L }, 0L)
  iu <- rep_len(as.character(sp$index.units), max(nidx, 1))[seq_len(nidx)]
  lh <- tryCatch(attr(stklen, "lhpar"), error = function(e) NULL)
  if (is.null(lh)) lh <- tryCatch(stklen@lhpar, error = function(e) NULL)
  h <- sp$h
  if (is.null(h)) h <- if (!is.null(lh) && "s" %in% dimnames(lh)$params) as.numeric(lh["s"])[1] else 0.7
  mucatch <- mean(catch)
  fsr <- if (isTRUE(tmb_data$nsex == 2L)) as.numeric(tmb_data$sex_ratio) else 1   # K is female SB0 with two sexes
  NT <- length(yrs)
  lk <- c(sp$link)
  if (!all(lk %in% c("spr", "none")))
    stop("settings$spm$link: \"spr\" or \"none\" (the F link was removed; it is kept on branch memory/spm-flink)",
         call. = FALSE)
  tmb_data <- c(tmb_data[!grepl("^spm_", names(tmb_data))], list(
    spm_use = 1L, spm_nyear = NT, spm_catch = pmax(catch / mucatch, 1e-8), spm_lfd_t = as.integer(lt - 1L),
    spm_nidx = as.integer(nidx), spm_logI = log(Im[ok]), spm_obs_t = as.integer(ok[, 1] - 1L),
    spm_obs_i = as.integer(ok[, 2] - 1L), spm_se2 = as.numeric(se2), spm_idx_gear = as.integer(igc),
    spm_idx_units = as.integer(iu == "biomass"), spm_h = as.numeric(h),
    spm_shape = as.integer(sp$shape == "spr"), spm_spr_tgt = as.numeric(sp$spr.tgt),
    spm_r_sd = as.numeric(sp$r.sd), spm_m_sd = as.numeric(sp$m.sd),
    spm_link_spr = as.integer("spr" %in% lk), spm_sigma_link = as.numeric(sp$sigma.link),
    spm_lagidx = as.integer(pmax(seq_len(NT) - 1L - round(sp$lag), 0L)),
    spm_bk_pr = as.numeric(sp$bkinit), spm_sdb_pr = as.numeric(sp$sdb),
    spm_K_pr = if (is.null(sp$K)) c(1, 1, 0) else if (identical(sp$K, "half8"))
                 c(5.7 * max(catch) / mucatch * fsr, 0.833, 2)
               else c(sp$K[1] / mucatch, sp$K[2], sp$K[3]),
    spm_xgrid = as.numeric(sp$xgrid), spm_newton_it = as.integer(sp$newton.it),
    spm_Pbound = as.numeric(sp$P.bound),
    spm_selw = .spm_selw(sp$sel.weights, tmb_data),
    spm_timing = as.integer(sp$timing == "start"),
    spm_catch_err = as.integer(isTRUE(sp$catch.error)), spm_sdc = as.numeric(sp$sdc %||% 0.1),
    spm_frw_t = .spm_frw_t(sp$frw %||% "none", NT, if (nidx) min(ok[, 1]) else NA),
    spm_sdf_pr = as.numeric(sp$sdf %||% c(0.3, 1, 1)),
    spm_info = list(years = yrs, mucatch = mucatch, index.names = if (nidx) colnames(I) else character(0),
                    obs = ok, settings = sp)))
  bk <- min(max(sp$bkinit[1], 0.05), 1.2)
  parameters[spm_names] <- list(log(10), log(0.3), log(2), log(bk), log(sp$sdb[1]),
                                rep(0, max(nidx, 1)), rep(log(0.2), max(nidx, 1)), rep(log(bk), NT), 0,
                                rep(0, NT), log((sp$sdf %||% 0.3)[1]))
  rl <- sp$link.rho
  if (length(rl) > 1) rl <- if ("spr" %in% names(rl)) rl[["spr"]] else rl[[1]]   # old c(spr = , F = )
  rl <- as.numeric(rl)
  if (length(rl) != 1 || !is.finite(rl) || abs(rl) >= 1)
    stop("settings$spm$link.rho: one fixed value in (-1, 1); it cannot be estimated", call. = FALSE)
  parameters$spm_trho <- stats::qlogis((rl + 1) / 2)
  map$spm_trho <- factor(NA)
  if (sp$r.sd <= 0) map$spm_logr <- factor(NA)
  if (sp$m.sd <= 0) map$spm_logm <- factor(NA)
  if (nidx == 0) { map$spm_logq <- factor(NA); map$spm_logsdi <- factor(NA) }
  sdi_fix <- if (!is.null(sp$sdi)) sp$sdi else if (isFALSE(sp$re) && is.null(sp$index.se)) 0.2 else NULL
  if (nidx > 0 && !is.null(sdi_fix)) {
    parameters$spm_logsdi <- log(rep_len(as.numeric(sdi_fix), nidx)); map$spm_logsdi <- factor(rep(NA, nidx))
  }
  if (!isTRUE(sp$catch.error)) map$spm_cdev <- factor(rep(NA, NT))
  if (!length(tmb_data$spm_frw_t)) map$spm_logsdf <- factor(NA)
  if (isFALSE(sp$re)) {          # penalised likelihood: states (and catch deviations) fixed effects, SDs fixed
    map$spm_logsdb <- factor(NA); map$spm_logsdf <- factor(NA)
    return(list(tmb_data = tmb_data, parameters = parameters, map = map, random = random))
  }
  list(tmb_data = tmb_data, parameters = parameters, map = map,
       random = c(random, "spm_logP", if (isTRUE(sp$catch.error)) "spm_cdev"))
}

# starting values: deterministic trajectories over a K grid with r and m of the
# per-recruit model at the initial values; q from the mean log ratio
spm_inits_flicc <- function(obj, tmb_data, parameters) {
  rp <- obj$report()
  r <- rp$r_pr; m <- rp$m_pr
  if (!is.finite(r) || r <= 0) r <- 0.3
  if (!is.finite(m) || m <= 0 || abs(m - 1) < 1e-3) m <- 2
  Pg <- rp$Pg; rg <- rp$rhog
  o <- order(Pg)
  rem <- function(p, col = 1) stats::approx(Pg[o], rg[o, col], xout = p, rule = 2)$y
  C <- tmb_data$spm_catch; NT <- length(C); lag <- tmb_data$spm_lagidx
  bk <- exp(parameters$spm_lbk)
  traj <- function(K) {
    P <- numeric(NT); P[1] <- bk
    for (t in 2:NT) { Pl <- P[lag[t - 1] + 1]
      P[t] <- max(P[t - 1] + r / (m - 1) * Pl * (1 - Pl^(m - 1)) - C[t - 1] / (K * rem(P[t - 1])), 0.02) }
    P
  }
  nidx <- tmb_data$spm_nidx
  score <- function(K) {
    P <- traj(K)
    pen <- 10 * mean(P <= 0.021)
    if (nidx == 0) return(list(s = pen + mean((log(P) - log(0.5))^2) * 0.01, P = P, lq = 0))
    ss <- 0; lq <- numeric(nidx)
    for (i in seq_len(nidx)) {
      w <- tmb_data$spm_obs_i == (i - 1); if (!any(w)) next
      t <- tmb_data$spm_obs_t[w] + 1
      ri <- if (tmb_data$spm_idx_gear[i] == -1) 1 else rem(P[t], 1 + i)
      x <- tmb_data$spm_logI[w] - log(K * P[t] * ri)
      lq[i] <- mean(x); ss <- ss + sum((x - lq[i])^2)
    }
    list(s = ss + pen * 100, P = P, lq = lq)
  }
  Ks <- exp(seq(log(1), log(1000), length.out = 120))
  sc <- lapply(Ks, score)
  b <- which.min(vapply(sc, `[[`, 0, "s"))
  parameters$spm_logK <- log(Ks[b])
  parameters$spm_logP <- log(pmax(sc[[b]]$P, 0.03))
  if (nidx > 0) parameters$spm_logq <- sc[[b]]$lq
  parameters$spm_logr <- log(r); parameters$spm_logm <- log(m)
  parameters
}

#' FLicc-SPM results
#'
#' Time series and reference points of the surplus production part of a
#' FLicc-SPM fit, in the layout of JABBA3 fits: \code{ts} (long data frame:
#' year, quantity, est, lci, uci), \code{refpts} and \code{pars}.
#' @param fit a \code{fiticc()} fit with \code{settings$spm}
#' @return list with \code{ts}, \code{refpts}, \code{pars}, \code{index} (fits)
#' @export
spm_flicc <- function(fit) {
  td <- fit$tmb_data
  if (!isTRUE(td$spm_use == 1L)) stop("not a FLicc-SPM fit", call. = FALSE)
  info <- td$spm_info; mu <- info$mucatch; yrs <- info$years
  rp <- fit$report                                    # fiticc_core(): TMB report; fiticc(): summary report
  if (is.null(rp$P) && !is.null(fit$obj)) rp <- fit$obj$report(fit$obj$env$last.par.best)
  sm <- if (!is.null(fit$rep)) summary(fit$rep, "report") else NULL
  get <- function(nm, est, scale = 1) {
    e <- est; lci <- uci <- rep(NA_real_, length(e))
    if (!is.null(sm)) { k <- rownames(sm) == nm
      if (sum(k) == length(e)) { s <- sm[k, ]; lci <- exp(s[, 1] - 1.96 * s[, 2]) * scale
        uci <- exp(s[, 1] + 1.96 * s[, 2]) * scale } }
    list(e = e, lci = lci, uci = uci)
  }
  q <- list(SB = get("spm_logB", exp(log(rp$P) + log(rp$spm_K)) * mu, mu),
            `SB/SB0` = get("spm_logP", rp$P), `B/Bmsy` = get("spm_logBBmsy", rp$BBmsy),
            `F/Fmsy` = get("spm_logHHmsy", rp$HHmsy))
  ts <- do.call(rbind, lapply(names(q), function(n) data.frame(year = yrs, quantity = n, est = q[[n]]$e,
                                                               lci = q[[n]]$lci, uci = q[[n]]$uci)))
  ts <- rbind(ts, data.frame(year = yrs, quantity = "H", est = rp$H, lci = NA, uci = NA),
              data.frame(year = yrs, quantity = "catch", est = td$spm_catch * mu, lci = NA, uci = NA),
              data.frame(year = as.numeric(td$year_names), quantity = "D_spr", est = rp$D_spr, lci = NA, uci = NA),
              data.frame(year = as.numeric(td$year_names), quantity = "SPR", est = rp$spr_y, lci = NA, uci = NA))
  refpts <- c(K = rp$spm_K * mu, r = rp$spm_r, m = rp$spm_m, BmsyK = rp$BmsyK_spm, Hmsy = rp$Hmsy_spm,
              MSY = rp$MSY_spm * mu, r_pr = rp$r_pr, m_pr = rp$m_pr, BmsyK_pr = rp$BmsyK_pr, Hmsy_pr = rp$Hmsy_pr,
              Fref = rp$Fref, spr_ref = rp$spr_ref, sdb = rp$spm_sdb, sdf = rp$sdf,
              rho_link = rp$rho_link)
  idx <- if (length(td$spm_logI)) data.frame(year = yrs[td$spm_obs_t + 1],
                                             index = info$index.names[td$spm_obs_i + 1],
                                             obs = exp(td$spm_logI), fit = exp(rp$spm_Imean)) else NULL
  list(ts = ts, refpts = refpts, index = idx, sdi = stats::setNames(rp$spm_sdi, info$index.names))
}

# years (0-based t >= 1) with a random-walk increment of log H: none, from the first index year (t0, 1-based row), all
.spm_frw_t <- function(frw, NT, t0) {
  if (identical(frw, "all")) return(as.integer(seq_len(NT - 1)))
  if (identical(frw, "index") && is.finite(t0) && t0 < NT) return(as.integer(seq(t0, NT - 1)))
  integer(0)
}
