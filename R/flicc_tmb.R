#' Build multi-year TMB data list for FLicc
#'
#' @param lfd named list of FLQuant objects, one per gear
#' @param stklen FLStockLen-like object with matching years
#' @param sel_fun selectivity function per gear
#' @param catch_by_gear numeric vector of relative catch by gear (in weight unless settings$catch_units = "numbers")
#' @param settings list with CVL, GL, catch.sd, optional linf.sd, Mk.sd, CVL.sd
#' @param years optional vector of years to include
#'
#' @return list for TMB::MakeADFun
data_tmb_flicc <- function(lfd, stklen, sel_fun, catch_by_gear,
                             settings = list(CVL = 0.1, GL = 50, catch.sd = 0.05),
                             years = NULL) {

  stopifnot(is.list(lfd))
  if (is.null(names(lfd))) {
    names(lfd) <- paste0("gear", seq_along(lfd))
  }
  gear_names <- names(lfd)
  ngear <- length(lfd)

  lens <- as.numeric(dimnames(lfd[[1]])$len)
  yrs_all <- as.numeric(dimnames(lfd[[1]])$year)

  if (is.null(years)) {
    years <- yrs_all
  }
  years <- as.numeric(years)
  nyear <- length(years)
  nlen  <- length(lens)

  # observed length frequencies: obs[len, year, gear]
  obs <- array(NA_real_,
               dim = c(nlen, nyear, ngear),
               dimnames = list(len = lens, year = years, gear = gear_names))

  for (g in seq_len(ngear)) {
    x <- lfd[[g]]
    xyrs <- as.numeric(dimnames(x)$year)
    if (!all(years %in% xyrs)) {
      stop("Not all requested years found in lfd gear: ", gear_names[g])
    }

    for (iy in seq_along(years)) {
      y <- years[iy]
      obs[, iy, g] <- as.numeric(as.data.frame(yearMeans(x[, ac(y)]))$data)
    }
  }

  # stock years
  stk_years <- as.numeric(dimnames(stklen)$year)
  if (!all(years %in% stk_years)) {
    stop("Not all requested years found in stklen")
  }

  stkw <- window(stklen, start = min(years), end = max(years))

  # extract year-varying biology as len x year matrices
  wt_df  <- as.data.frame(catch.wt(stkw))
  mat_df <- as.data.frame(mat(stkw))
  m_df   <- as.data.frame(m(stkw))

  # assume ordering len-fast, then year
  wt_mat  <- matrix(wt_df$data,  nrow = nlen, ncol = nyear)
  mat_mat <- matrix(mat_df$data, nrow = nlen, ncol = nyear)
  m_mat   <- matrix(m_df$data,   nrow = nlen, ncol = nyear)

  # scale M around annual mean so Mk remains the average level
  Mscaler <- sweep(m_mat, 2, colMeans(m_mat, na.rm = TRUE), "/")

  # lengths
  LLB <- lens
  step <- c(diff(LLB), tail(diff(LLB), 1))
  Lmid <- LLB + step / 2

  # selectivity
  if (length(sel_fun) != ngear) {
    stop("length(sel_fun) must equal number of gears")
  }

  sel_type <- unname(c(logistic = 1L, dsnormal = 2L, normal = 3L, binormal = 4L,
                       bidsnormal = 5L)[sel_fun])
  if (any(is.na(sel_type))) {
    stop("Supported selectivity options are logistic, dsnormal, normal, binormal, bidsnormal")
  }

  sm_n <- unname(c(2L, 3L, 2L, 5L, 6L)[sel_type])
  sm_start <- cumsum(c(0L, head(sm_n, -1L)))  # 0-based for C++





  # catch weights by gear
  catch_all <- do.call(
    cbind,
    lapply(catch_by_gear, function(x) as.numeric(x))
  )
  rownames(catch_all) <- dimnames(catch_by_gear[[1]])$year
  # catch may cover more years than the fitted ones (e.g. fiticc(years = ))
  catch_wt <- if (all(as.character(years) %in% rownames(catch_all)) && nrow(catch_all) != nyear)
    catch_all[as.character(years), , drop = FALSE] else catch_all

  catch_wt <- catch_wt / rowSums(catch_wt)


  if (ncol(catch_wt) != ngear) {
    stop("length(catch_by_gear) must equal number of gears")
  }
  if (nrow(catch_wt) != nyear) {
    stop("catch_by_gear must have one value per year")
  }

  catch_sd <- if (!is.null(settings$catch.sd)) as.numeric(settings$catch.sd) else 0.05

  # observation model
  obs_model <- if (!is.null(settings$obs_model)) settings$obs_model else "dm"

  obs_model_code <- switch(obs_model,
                           nb = 1L,
                           mn = 2L,
                           dm = 3L,
                           stop("Supported obs_model options are 'nb', 'mn', 'dm'")
  )

  catch_units <- if (!is.null(settings$catch_units)) settings$catch_units else "weight"
  if (!catch_units %in% c("weight", "numbers"))
    stop("settings$catch_units must be \"weight\" or \"numbers\".")
  dm_theta <- if (!is.null(settings$dm_theta)) settings$dm_theta else "saturating"
  dm_type <- switch(dm_theta, fixed = 0L, linear = 1L, saturating = 2L,
                    stop("settings$dm_theta must be \"fixed\", \"linear\" or \"saturating\"."))

  pop_model <- if (!is.null(settings$pop_model)) settings$pop_model else "gamma"

  pop_model_code <- switch(pop_model,
                           gamma = 1L,
                           gtg   = 2L,
                           stop("Supported pop_model options are 'gamma' and 'gtg'")
  )



  # quadrature
  gl <- statmod::gauss.quad(settings$GL, kind = "laguerre")
  node <- as.numeric(gl$nodes)
  quad_wt <- as.numeric(gl$weights)

  # life history initial values
  lhpar <- stklen@lhpar
  Linf_init <- as.numeric(lhpar["linf"])
  Mk_init   <- as.numeric(lhpar["Mk"])

  CVL_init <- if (!is.null(settings$CVL)) as.numeric(settings$CVL) else 0.1
  Galpha_init <- 1 / CVL_init^2

  # GTG controls
  ngtg  <- if (!is.null(settings$ngtg))  as.integer(settings$ngtg)  else 13L
  maxsd <- if (!is.null(settings$maxsd)) as.numeric(settings$maxsd) else 3
  Mpow  <- if (!is.null(settings$Mpow))  as.numeric(settings$Mpow)  else 0

  # GTG distribution over Linf on a standardised grid. The Linf groups are
  # built inside the TMB template as Linf * (1 + CVL * gtg_z), so they follow
  # the estimated Linf and CVL. recP depends only on the grid.
  # one growth type (e.g. LBB): the single group sits at the mean Linf
  gtg_z <- if (ngtg == 1L) 0 else seq(-maxsd, maxsd, length.out = ngtg)
  recP  <- dnorm(gtg_z)
  recP  <- recP / sum(recP)

  # Bin boundaries for GTG mortality scaling; M/K at each boundary is
  # Mk * (Linf / LBins_mk)^Mpow, computed in the template from the estimates.
  # A zero lower boundary would give infinite M/K when Mpow > 0, so a
  # non-positive boundary is replaced by half the first bin width.
  LBins    <- c(LLB, tail(LLB, 1) + tail(step, 1))
  LBins_mk <- ifelse(LBins > 0, LBins, step[1] / 2)

  # SPR targets whose Fspr is computed at report time (not during fitting)
  spr_ref <- if (!is.null(settings$spr_ref)) as.numeric(settings$spr_ref) else 40
  if (any(!is.finite(spr_ref)) || any(spr_ref <= 0) || any(spr_ref >= 100)) {
    stop("settings$spr_ref must be SPR percentages in (0, 100).")
  }

  # Upper-tail robustness for the composition likelihood
  rob_eps <- if (!is.null(settings$rob_eps)) as.numeric(settings$rob_eps) else 0
  if (length(rob_eps) != 1 || !is.finite(rob_eps) || rob_eps < 0 || rob_eps >= 1) {
    stop("settings$rob_eps must be a single value in [0, 1).")
  }
  plus_bin <- -1L
  if (!is.null(settings$Lplus)) {
    Lplus <- as.numeric(settings$Lplus)
    idx <- which(LLB >= Lplus)
    if (length(idx) == 0) {
      warning("settings$Lplus = ", Lplus, " is above the largest length bin; no plus group used.")
    } else {
      plus_bin <- as.integer(min(idx) - 1L)   # 0-based for C++
    }
  }

  # Tail cut: exclude length bins with lower bound >= tail_cut (an absolute
  # length, as Lplus, so the same bins are used whatever Linf is assumed) from
  # the composition likelihood. Observed and predicted compositions are
  # truncated and the predicted proportions renormalised over the retained bins.
  cut_bin <- -1L
  if (!is.null(settings$tail_cut)) {
    Lcut <- tail_cut_length_flicc(settings$tail_cut, Linf_init)
    idx  <- which(LLB >= Lcut)
    if (length(idx) == 0) {
      warning("settings$tail_cut = ", round(Lcut, 2),
              " is above the largest length bin; no bins cut.")
    } else {
      cut_bin <- as.integer(min(idx) - 1L)    # 0-based first excluded bin
      if (cut_bin < 2L) {
        stop("settings$tail_cut = ", round(Lcut, 2), " would leave fewer than two ",
             "length bins in the likelihood.")
      }
      if (plus_bin >= cut_bin) {
        warning("settings$Lplus lies at or above the tail cut; plus group ignored.")
        plus_bin <- -1L
      }
    }
  }

  # Warn when observations fall in bins the GTG model cannot reach
  # (lower bound at or beyond the largest GTG Linf at the input values).
  if (pop_model == "gtg" && isTRUE(settings$tail_warning)) {
    maxLinf <- Linf_init * (1 + max(gtg_z) * CVL_init)
    beyond  <- LLB >= maxLinf
    if (plus_bin >= 0 && LLB[plus_bin + 1] < maxLinf) beyond <- beyond & seq_along(LLB) <= plus_bin
    if (cut_bin >= 0) beyond <- beyond & seq_along(LLB) <= cut_bin
    n_beyond <- sum(obs[beyond, , , drop = FALSE], na.rm = TRUE)
    if (n_beyond > 0) {
      warning(sprintf(paste0(
        "%s observed fish (%.1f%% of the total) lie in length bins at or above the ",
        "largest GTG Linf (%.1f = Linf * (1 + maxsd * CVL)); the model gives them ",
        "near-zero probability, which can drive CVL and status estimates. ",
        "Consider a larger maxsd, settings$Lplus or settings$rob_eps."),
        format(round(n_beyond, 1)), 100 * n_beyond / sum(obs, na.rm = TRUE), maxLinf),
        call. = FALSE)
    }
  }


  prior_mu <- numeric()
  prior_sd <- numeric()
  prior_code <- integer()

  if (!is.null(settings$linf.sd)) {
    prior_mu   <- c(prior_mu, log(Linf_init))
    prior_sd   <- c(prior_sd, as.numeric(settings$linf.sd))
    prior_code <- c(prior_code, 1L)
  }

  if (!is.null(settings$Mk.sd)) {
    prior_mu   <- c(prior_mu, log(Mk_init))
    prior_sd   <- c(prior_sd, as.numeric(settings$Mk.sd))
    prior_code <- c(prior_code, 3L)
  }

  if (!is.null(settings$CVL.sd)) {
    prior_mu   <- c(prior_mu, -2 * log(CVL_init))
    prior_sd   <- c(prior_sd, 2 * as.numeric(settings$CVL.sd))
    prior_code <- c(prior_code, 2L)
  }

  list(
    nlen       = nlen,
    nyear      = nyear,
    ngear      = ngear,
    Lmid       = as.numeric(Lmid),
    LLB        = as.numeric(LLB),
    obs        = unname(obs),
    Mscaler    = unname(Mscaler),
    wt         = unname(wt_mat),
    mat        = unname(mat_mat),
    catch_wt   = unname(catch_wt),
    catch_sd   = catch_sd,
    pop_model  = as.integer(pop_model_code),
    obs_model  = as.integer(obs_model_code),
    catch_units = as.integer(catch_units == "weight"),
    dm_type    = as.integer(dm_type),
    dm_prior_mu = if (dm_type == 2L) log(10 * pmax(apply(apply(obs, c(2, 3), sum, na.rm = TRUE), 2, max), 1))
                  else rep(log(10), ngear),
    dm_prior_sd = if (!is.null(settings$dm_prior_sd)) as.numeric(settings$dm_prior_sd) else 1,
    pile = if (!is.null(settings$pile)) as.numeric(settings$pile) else 1,
    ngtg       = as.integer(ngtg),
    gtg_z      = as.numeric(gtg_z),
    recP       = as.numeric(recP),
    LBins_mk   = as.numeric(LBins_mk),
    Mpow       = as.numeric(Mpow),
    rob_eps    = as.numeric(rob_eps),
    plus_bin   = as.integer(plus_bin),
    cut_bin    = as.integer(cut_bin),
    # per-recruit block, off during fitting (see brp_tmb_flicc())
    brp_F       = numeric(0),
    brp_spr     = as.numeric(spr_ref) / 100,   # Fspr reported by the fit
    brp_nyears  = 1L,
    brp_k       = as.numeric(lhpar["k"]),
    brp_FMscale = 0L,
    brp_spawn   = 0,
    brp_Fmax    = 5,
    brp_only    = 0L,
    sel_type   = as.integer(sel_type),
    sm_start   = as.integer(sm_start),
    sm_n       = as.integer(sm_n),
    node       = node,
    quad_wt    = quad_wt,
    prior_mu   = as.numeric(prior_mu),
    prior_sd   = as.numeric(prior_sd),
    prior_code = as.integer(prior_code),
    Linf_init   = Linf_init,
    Mk_init     = Mk_init,
    CVL_init    = CVL_init,
    Galpha_init = Galpha_init,
    gear_names  = gear_names,
    year_names  = years
  )
}

#' Tail cut as a length (internal)
#'
#' `tail_cut` is an absolute length (same units as the length bins, as
#' `Lplus`). For backward compatibility a value <= 1.5 is read as a fraction of
#' Linf (the earlier definition), with a warning.
#' @param tail_cut single positive number
#' @param linf input Linf
#' @return the cut length
#' @noRd
tail_cut_length_flicc <- function(tail_cut, linf) {
  tc <- as.numeric(tail_cut)
  if (length(tc) != 1 || !is.finite(tc) || tc <= 0) {
    stop("tail_cut must be a single positive length (e.g. 50), as Lplus.")
  }
  if (tc <= 1.5) {
    warning(sprintf(paste0("tail_cut = %g is read as a fraction of Linf (%.1f); tail_cut is now an ",
                           "absolute length, as Lplus: use tail_cut = %.1f."), tc, linf, tc * linf),
            call. = FALSE)
    tc <- tc * linf
  }
  tc
}

#' Generate initial parameter values for FLicc TMB fits
#'
#' Builds a named parameter list suitable for \code{TMB::MakeADFun()} from a
#' TMB data list produced by \code{data_tmb_flicc()}. Initial values are derived
#' from the data where possible, using observed modal length by gear as the
#' starting point for selectivity location parameters.
#'
#' @param tmb_data A list produced by \code{data_tmb_flicc()}.
#'
#' @details
#' Initial values are currently set as follows:
#' \itemize{
#'   \item \code{log_Linf}: initialized from \code{tmb_data$Linf_init}
#'   \item \code{log_Mk}: initialized from \code{tmb_data$Mk_init}
#'   \item \code{log_Galpha}: initialized to \code{log(20)}
#'   \item \code{log_Fk}: initialized to \code{log(0.5)} for all gears
#'   \item \code{log_phi}: initialized to \code{log(20)}
#'   \item \code{Sm}: selectivity parameters initialized by gear using the
#'     observed modal length and default spread values
#' }
#'
#' Selectivity parameters are stored in a packed vector \code{Sm}, with the
#' meaning depending on selectivity type.
#'
#' For logistic selectivity, the packed parameter block is
#' \code{c(log_SL50_rel, log_dSL_rel)}, where \code{SL50 = exp(log_SL50_rel) * Linf}
#' and \code{SL95 = SL50 + exp(log_dSL_rel) * Linf}.
#'
#' For normal selectivity, the packed parameter block is
#' \code{c(log_mode_rel, log_sd_rel)}, where \code{mode = exp(log_mode_rel) * Linf}
#' and \code{sd = exp(log_sd_rel) * Linf}.
#'
#' For double-sided normal selectivity, the packed parameter block is
#' \code{c(log_mode_rel, log_lsd_rel, log_rsd_rel)}, where
#' \code{mode = exp(log_mode_rel) * Linf}, \code{lsd = exp(log_lsd_rel) * Linf},
#' and \code{rsd = exp(log_rsd_rel) * Linf}.
#'
#' @return A named list of parameter vectors and scalars for use in
#'   \code{TMB::MakeADFun()}.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' tmb_data <- data_tmb_flicc(lfd,stklen,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))
#' par0 <- init_tmb_flicc(tmb_data)
#' str(par0)
#' }
#'
#' @export
#' Initial values for multi-year FLicc TMB
init_tmb_flicc <- function(tmb_data) {

  nlen     <- tmb_data$nlen
  nyear    <- tmb_data$nyear
  ngear    <- tmb_data$ngear
  sel_type <- tmb_data$sel_type
  Lmid     <- tmb_data$Lmid
  obs      <- tmb_data$obs
  Linf     <- tmb_data$Linf_init

  Sm <- numeric(sum(tmb_data$sm_n))
  pos <- 1L

  for (g in seq_len(ngear)) {
    obs_g <- apply(obs[, , g, drop = FALSE], 1, sum, na.rm = TRUE)
    mode_g <- Lmid[which.max(obs_g)]
    mode_rel <- max(mode_g / Linf, 1e-4)

    if (sel_type[g] == 1L) {
      SL50_rel <- mode_rel
      dSL_rel  <- 0.10
      Sm[pos:(pos + 1L)] <- c(log(SL50_rel), log(dSL_rel))
      pos <- pos + 2L

    } else if (sel_type[g] == 2L) {
      lsd_rel <- 0.10
      rsd_rel <- 0.10
      Sm[pos:(pos + 2L)] <- c(log(mode_rel), log(lsd_rel), log(rsd_rel))
      pos <- pos + 3L

    } else if (sel_type[g] == 3L) {
      sd_rel <- 0.10
      Sm[pos:(pos + 1L)] <- c(log(mode_rel), log(sd_rel))
      pos <- pos + 2L

    } else if (sel_type[g] == 4L) {
      # binormal: modes at the 20% and 80% quantiles of the observed lengths
      cs <- cumsum(obs_g) / max(sum(obs_g), 1e-12)
      q1 <- Lmid[which(cs >= 0.2)[1]]; q2 <- Lmid[which(cs >= 0.8)[1]]
      if (!is.finite(q1) || !is.finite(q2) || q2 <= q1) { q1 <- 0.8 * mode_g; q2 <- 1.2 * mode_g }
      Sm[pos:(pos + 4L)] <- c(log(max(q1 / Linf, 1e-4)), log(0.05), log(max((q2 - q1) / Linf, 1e-3)),
                              log(0.05), 0)
      pos <- pos + 5L

    } else if (sel_type[g] == 5L) {
      # bidsnormal: as binormal, second curve double-sided (lsd2, rsd2)
      cs <- cumsum(obs_g) / max(sum(obs_g), 1e-12)
      q1 <- Lmid[which(cs >= 0.2)[1]]; q2 <- Lmid[which(cs >= 0.8)[1]]
      if (!is.finite(q1) || !is.finite(q2) || q2 <= q1) { q1 <- 0.8 * mode_g; q2 <- 1.2 * mode_g }
      Sm[pos:(pos + 5L)] <- c(log(max(q1 / Linf, 1e-4)), log(0.05), log(max((q2 - q1) / Linf, 1e-3)),
                              log(0.05), log(0.10), 0)
      pos <- pos + 6L
    }
  }

  gear_wt <- colMeans(as.matrix(tmb_data$catch_wt))

  if (length(gear_wt) != ngear ||
      any(!is.finite(gear_wt)) ||
      sum(gear_wt) <= 0) {
    gear_wt <- rep(1 / ngear, ngear)
  } else {
    gear_wt <- gear_wt / sum(gear_wt)
  }

  Fk_init <- matrix(
    rep(tmb_data$Mk_init * gear_wt, each = nyear),
    nrow = nyear,
    ncol = ngear,
    dimnames = list(tmb_data$year_names, tmb_data$gear_names)
  )

  log_Fk <- log(pmax(Fk_init, 1e-8))

  list(
    log_Linf   = log(tmb_data$Linf_init),
    log_Galpha = log(tmb_data$Galpha_init),
    log_Mk     = log(tmb_data$Mk_init),
    log_Fk     = log_Fk,
    Sm         = Sm,
    log_phi    = log(10),
    log_dm     = if (!is.null(tmb_data$dm_prior_mu)) as.numeric(tmb_data$dm_prior_mu) else rep(0, length(tmb_data$gear_names)),
    log_sigmaF = tmb_data$prior_sigmaF_mean
  )
}

#' Fit the FLicc model in a single optimizer pass
#'
#' Internal fitting engine for FLicc. This function performs one model fit only
#' and does not apply any automatic retry or stability refit logic.
#'
#' @param lfd A named list of \code{FLQuant} objects containing observed
#'   length-frequency distributions by gear.
#' @param stklen An \code{FLStockLen} object, or compatible FLR object,
#'   containing stock and life-history information.
#' @param sel_fun Character vector giving the selectivity function for each
#'   gear. Supported options include \code{"logistic"}, \code{"normal"},
#'   \code{"dsnormal"}, \code{"binormal"} and \code{"bidsnormal"}.
#' @param catch_by_gear Relative catch by gear, in weight (e.g. landings in
#'   tonnes) unless \code{settings$catch_units = "numbers"}. Can be supplied as a numeric
#'   vector or as an \code{FLQuants} object.
#' @param settings A named list of model settings. Recognized elements are
#'   \code{CVL}, \code{GL}, \code{catch.sd}, and \code{prior_sigmaF}. Missing
#'   elements are filled from internal defaults.
#' @param years Optional numeric vector of years to include in the fit. If
#'   \code{NULL}, all years in \code{lfd} are used.
#' @param iter Integer giving the iteration to fit. Iterations are subset in R
#'   before being passed to TMB.
#' @param compile Logical; if \code{TRUE}, compile the TMB model before fitting.
#' @param silent Logical; passed to \code{TMB::MakeADFun()}.
#' @param dll Character string giving the compiled TMB DLL name.
#' @param sel_fixed Optional object used to fix selectivity parameters.
#' @param start Optional numeric vector of starting values for the estimated
#'   TMB parameters. This is mainly used to restart a fit from
#'   \code{fit$opt$par}. If supplied as a vector, its length must match the
#'   free parameter vector after applying any TMB parameter mapping.
#' @param FLRreport Logical; if \code{TRUE}, convert reported outputs using
#'   \code{as_FLQuants()}.
#' @param newton Maximum number of Newton polishing steps (0 = none).
#'
#' @details
#' \code{fiticc_core()} prepares the TMB data using \code{data_tmb_flicc()},
#' initializes model parameters, applies optional selectivity constraints, builds
#' the TMB objective function with \code{TMB::MakeADFun()}, optimizes it using
#' \code{nlminb()}, and returns both the optimizer output and reported
#' quantities.
#'
#' This function is primarily intended for internal use. Most users should call
#' \code{fiticc()}, which wraps \code{fiticc_core()} and can automatically
#' refit marginal models with a higher \code{GL}.
#'
#' @return
#' An object of class \code{"flicc_tmb_fit"} with optimizer results, sdreport
#' output, parameter summaries, reported quantities, and the processed TMB data.
#'
#' @seealso \code{\link{fiticc}}, \code{\link{need_refit_flicc}}
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc_core(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#' @keywords internal
fiticc_core <- function(lfd, stklen,
                     sel_fun = c("logistic", "normal", "dsnormal"),
                     catch_by_gear,
                     settings = list(
                       CVL = 0.1,
                       GL = 50,
                       catch.sd = 0.05
                     ),
                     years = NULL,
                     iter = 1,
                     compile = FALSE,
                     silent = TRUE,
                     dll = "FLicc",
                     sel_fixed = NULL,
                     start = NULL,
                     n_restart = 3,
                     grad_tol = 1e-3,
                     control = list(),
                     FLRreport=FALSE,
                     se = TRUE,
                     newton = 3) {

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required")
  }


  default_settings <- list(
    CVL = 0.1,
    GL = 50,
    catch.sd = 0.05,
    pop_model = "gtg",
    obs_model = "dm",
    dm_theta = "saturating",
    dm_shared = FALSE,
    catch_units = "weight",
    dm_prior_sd = 1,
    pile = 1,
    ngtg = 13,
    maxsd = 3,
    Mpow = 0,
    prior_sigmaF = c(log(0.5), 0.3, 1),
    linf.sd = NULL,
    Mk.sd = NULL,
    CVL.sd = NULL,
    rob_eps = 0,
    Lplus = NULL,
    F_re = FALSE,
    sigmaF = NULL,
    tail_cut = NULL,
    spr_ref = 40,
    tail_warning = TRUE,
    FM_min = 0.05,
    FM_max = 4,
    FMpen_sd = 0.2
  )
  # ---- merge user settings with defaults ----
  if (is.null(settings)) settings <- list()

  if(isTRUE(any(settings$prior_sigmaF == FALSE, na.rm = TRUE))) settings$prior_sigmaF <- c(NA_real_, NA_real_, 0)

  # Default (FLicc >= 1.0.10): random walk on log F with sigmaF fixed at 0.3
  # (fast, no collapse of sigmaF). sigmaF is estimated only if a prior_sigmaF
  # is given (penalised) or with F_re = TRUE (random effects).
  if (is.null(settings$sigmaF) && is.null(settings$prior_sigmaF) && !isTRUE(settings$F_re)) {
    settings$sigmaF <- 0.3
  }

  # F as random effects: weak default prior on sigmaF (median 0.3, sd 0.5 on
  # the log scale), so sigmaF cannot collapse to 0 when the length data carry
  # little year-to-year information. A user-supplied prior_sigmaF is kept.
  if (isTRUE(settings$F_re) && is.null(settings$prior_sigmaF)) {
    settings$prior_sigmaF <- c(log(0.3), 0.5, 1)
  }

  # fill missing elements
  for (nm in names(default_settings)) {
    if (is.null(settings[[nm]])) {
      settings[[nm]] <- default_settings[[nm]]
    }
  }

  if (!settings$pop_model %in% c("gamma", "gtg")) {
    stop("settings$pop_model must be 'gamma' or 'gtg'")
  }

  if (!settings$obs_model %in% c("nb", "mn", "dm")) {
    stop("settings$obs_model must be 'nb', 'mn', or 'dm'")
  }

  if (!is.null(settings$sigmaF) &&
      (length(settings$sigmaF) != 1 || !is.finite(as.numeric(settings$sigmaF)) || settings$sigmaF <= 0)) {
    stop("settings$sigmaF must be a single positive value (the fixed random-walk SD of log F). ",
         "To switch the random walk off use prior_sigmaF = FALSE.")
  }
  if (isTRUE(settings$prior_sigmaF[3] == 1) && is.null(settings$sigmaF) && !isTRUE(settings$F_re) &&
      !(is.finite(settings$prior_sigmaF[2]) && settings$prior_sigmaF[2] > 0)) {
    warning("Random walk on log F with sigmaF estimated without a prior: sigmaF tends to collapse to 0 ",
            "(flat F). Fix it (sigmaF = 0.3), use F_re = TRUE, or switch the walk off ",
            "(prior_sigmaF = FALSE).", call. = FALSE)
  }

  unknown <- setdiff(names(settings), names(default_settings))
  if (length(unknown) > 0) {
    stop("Unknown settings: ", paste(unknown, collapse = ", "))
  }

  # subset iter in R; do not loop over iter in TMB
  lfd_i <- lapply(lfd, function(x) x[, , iter])
  stk_i <- stklen[, , , iter]

  if (is.null(years)) {
    years <- as.numeric(dimnames(lfd_i[[1]])$year)
  }
  years <- as.numeric(years)

  if(!inherits(catch_by_gear, "FLQuants")){

 catch_by_gear  <-  FLCore::FLQuants(Map(function(x,y){
    flq = FLCore::quantSums(x)
    flq[] <- y
    flq
  },x=lfd_i,y=catch_by_gear))

 }

  tmb_data <- data_tmb_flicc(
    lfd = lfd_i,
    stklen = stk_i,
    sel_fun = sel_fun,
    catch_by_gear = catch_by_gear,
    settings = settings,
    years = years
  )

  # Set priors for random walk on F (if nyear > 0)
  tmb_data$prior_sigmaF_mean <- as.numeric(settings$prior_sigmaF[1])
  tmb_data$prior_sigmaF_sd   <- as.numeric(settings$prior_sigmaF[2])
  tmb_data$prior_sigmaF_use  <- as.integer(settings$prior_sigmaF[3])
  # no prior on sigmaF when its sd is NA / Inf / <= 0 (useful with F_re = TRUE)
  if (!is.finite(tmb_data$prior_sigmaF_sd) || tmb_data$prior_sigmaF_sd <= 0) {
    tmb_data$prior_sigmaF_sd <- 0
  }
  if (!is.finite(tmb_data$prior_sigmaF_mean)) {
    tmb_data$prior_sigmaF_mean <- log(0.5)       # start value only
  }

  # Soft penalty FM
  tmb_data$FM_min    <- as.numeric(settings$FM_min)
  tmb_data$FM_max    <- as.numeric(settings$FM_max)
  tmb_data$FMpen_sd  <- as.numeric(settings$FMpen_sd)

  parameters <- init_tmb_flicc(tmb_data)

  fix_res <- apply_sel_fixed_flicc(
    parameters = parameters,
    tmb_data   = tmb_data,
    sel_fixed  = sel_fixed
  )
  parameters <- fix_res$parameters

  map <- list()

  if (is.null(settings$linf.sd)) {
    map$log_Linf <- factor(NA)
  }
  if (is.null(settings$Mk.sd)) {
    map$log_Mk <- factor(NA)
  }
  if (is.null(settings$CVL.sd)) {
    map$log_Galpha <- factor(NA)
  }
  # NEW: phi only relevant for NB
  if (settings$obs_model != "nb") {
    map$log_phi <- factor(NA)
  }
  # Dirichlet-multinomial precision: estimated by gear (or shared) only for
  # obs_model = "dm" with dm_theta "linear" or "saturating"
  ng <- length(tmb_data$gear_names)
  if (settings$obs_model != "dm" || isTRUE(tmb_data$dm_type == 0L)) {
    map$log_dm <- factor(rep(NA, ng))
  } else if (isTRUE(settings$dm_shared)) {
    map$log_dm <- factor(rep(1, ng))
  }

  rw_on <- length(unique(years)) > 1 && isTRUE(tmb_data$prior_sigmaF_use == 1L)
  if (!rw_on) {
    map$log_sigmaF <- factor(NA)
  }

  # F as a state-space random walk: log_Fk integrated out by the Laplace
  # approximation, so sigmaF is estimated from the marginal likelihood
  # sigmaF fixed at a given value (e.g. estimated once with F_re = TRUE)
  if (!is.null(settings$sigmaF) && rw_on) {
    parameters$log_sigmaF <- log(as.numeric(settings$sigmaF))
    map$log_sigmaF <- factor(NA)
    tmb_data$prior_sigmaF_sd <- 0
  }

  random <- NULL
  if (isTRUE(settings$F_re)) {
    if (rw_on) {
      random <- "log_Fk"
    } else {
      warning("settings$F_re = TRUE ignored: needs more than one year and ",
              "prior_sigmaF[3] = 1 (random walk on).", call. = FALSE)
    }
  }

  # CVL (log_Galpha) is fixed unless CVL.sd is given, for both pop models.
  # For gtg it sets the spread of the Linf groups.
  if (!is.null(fix_res$map_Sm)) {
    map$Sm <- fix_res$map_Sm
  }
  if (length(map) == 0) {
    map <- NULL
  }

  if (compile) {
    TMB::compile("src/FLicc.cpp")
    dyn.load(dynlib("src/FLicc"))
  }

  obj <- TMB::MakeADFun(
    data = tmb_data[setdiff(names(tmb_data),
                            c("gear_names", "year_names",
                              "Linf_init", "Mk_init", "CVL_init", "Galpha_init"))],
    parameters = parameters,
    map = map,
    random = random,
    DLL = dll,
    silent = silent
  )

  # Optional restart from user-supplied starting values.
  # This is useful for reducing gradients by restarting from a previous optimum,
  # e.g. start = fit$opt$par. The supplied vector must match obj$par after
  # mapping, because mapped/fixed parameters are not optimized by nlminb().
  opt_start <- obj$par
  if (!is.null(start)) {
    if (is.list(start)) start <- unlist(start)
    start <- as.numeric(start)

    if (length(start) != length(opt_start)) {
      stop(
        "Length of 'start' (", length(start), ") does not match the number ",
        "of estimated TMB parameters (", length(opt_start), ")."
      )
    }

    names(start) <- names(opt_start)
    opt_start <- start
  }

  ctrl <- modifyList(
    list(
      iter.max = 1000,
      eval.max = 2000,
      rel.tol =  1e-9,
      x.tol =  1e-7
    ),
    control
  )

  opt <- nlminb(
    start = opt_start,
    objective = obj$fn,
    gradient = obj$gr,
    control = ctrl
  )

  # Newton polishing: a few Newton steps from the nlminb solution using the
  # Hessian of the objective (AD for fixed-effects models, finite differences
  # of the gradient with random effects). Steps that do not lower the
  # objective are rejected. Also gives the Hessian check and the Newton
  # decrement 0.5 g'H^-1 g (expected remaining improvement in log-likelihood).
  hess_fun <- function(p) {
    H <- if (is.null(random)) tryCatch(obj$he(p), error = function(e) NULL) else NULL
    if (is.null(H)) H <- tryCatch(stats::optimHess(p, obj$fn, obj$gr), error = function(e) NULL)
    if (is.null(H) || any(!is.finite(H))) return(NULL)
    (H + t(H)) / 2
  }
  newton_state <- function(p) {
    g <- as.numeric(obj$gr(p))
    H <- hess_fun(p)
    R <- if (is.null(H)) NULL else tryCatch(chol(H), error = function(e) NULL)
    if (is.null(R)) return(list(pd = if (is.null(H)) NA else FALSE, dec = NA_real_, step = NULL,
                                max_grad = max(abs(g), na.rm = TRUE)))
    step <- as.numeric(backsolve(R, forwardsolve(t(R), g)))
    list(pd = TRUE, dec = 0.5 * sum(g * step), step = step, max_grad = max(abs(g), na.rm = TRUE))
  }
  polish <- function(opt, nsteps) {
    st <- newton_state(opt$par); k <- 0L
    while (k < nsteps && isTRUE(st$pd) && is.finite(st$dec) && st$dec > 1e-12) {
      p_new <- opt$par - st$step
      f_new <- tryCatch(obj$fn(p_new), error = function(e) NA_real_)
      if (!is.finite(f_new) || f_new > opt$objective + 1e-10) break
      opt$par[] <- p_new; opt$objective <- f_new; k <- k + 1L
      st <- newton_state(opt$par)
    }
    opt$newton <- list(steps = k, hess_pd = st$pd, decrement = st$dec, max_gradient = st$max_grad)
    opt
  }
  newton <- as.integer(newton)
  if (newton > 0L) opt <- polish(opt, newton)

  restart_log <- data.frame(
    restart = 0L,
    objective = opt$objective,
    max_gradient = max(abs(obj$gr(opt$par)), na.rm = TRUE),
    improved = TRUE
  )

  n_restart_used <- 0L

  for (i in seq_len(n_restart)) {

    old_grad <- max(abs(obj$gr(opt$par)), na.rm = TRUE)

    if (old_grad < grad_tol) break

    opt_new <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = ctrl
    )

    new_grad <- max(abs(obj$gr(opt_new$par)), na.rm = TRUE)

    accept <- is.finite(opt_new$objective) &&
      is.finite(new_grad) &&
      new_grad < old_grad &&
      opt_new$objective <= opt$objective + 1e-6

    restart_log <- rbind(
      restart_log,
      data.frame(
        restart = i,
        objective = opt_new$objective,
        max_gradient = new_grad,
        improved  = accept
      )
    )

    if (accept) {
      opt_new$newton <- opt$newton
      opt <- opt_new
      n_restart_used <- n_restart_used + 1L
    } else {
      break
    }
  }

  if (newton > 0L && n_restart_used > 0L) opt <- polish(opt, newton)
  opt$max_gradient <- max(abs(obj$gr(opt$par)), na.rm = TRUE)
  opt$n_restarts <- n_restart_used
  opt$restart_log <- restart_log

  if (is.null(random)) {
    obj$env$last.par.best <- opt$par
  } else {
    # evaluate at the optimum so the random effects sit at their mode;
    # last.par.best then holds fixed and random parameters
    invisible(obj$fn(opt$par))
    obj$env$last.par.best <- obj$env$last.par
  }

  if (isTRUE(se)) {
    rep <- TMB::sdreport(obj)
    adrep <- summary(rep, "fixed")
  } else {
    # no standard errors requested (e.g. inside an MSE): skip sdreport
    rep <- NULL
    adrep <- cbind(Estimate = opt$par, `Std. Error` = NA_real_)
  }
  report <- obj$report(obj$env$last.par.best)
  par_tab <- data.frame(
    par = rownames(adrep),
    mpd = adrep[, 1],
    se  = adrep[, 2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )


  fit <- structure(
    list(
      opt = opt,
      obj = obj,
      rep = rep,
      par = par_tab,
      report = report,
      stklen= stklen,
      tmb_data = tmb_data
    ),
    class = "flicc_tmb_fit"
  )

  fit$settings <- settings
  fit$pop_model <- settings$pop_model
  fit$F_re <- !is.null(random)
  fit$obs_model <- settings$obs_model
  fit$start_used <- !is.null(start)
  fit$max_gradient <- opt$max_gradient
  fit$n_restarts <- opt$n_restarts
  fit$restart_log <- opt$restart_log
  fit$sel_fixed <- sel_fixed

  fit$logLik <- LLflicc(fit)
  fit$diagnostics <- diagnostics_flicc(fit)
  # a polished fit with a positive definite Hessian and a negligible Newton
  # decrement is converged, whatever nlminb's stopping code
  if (identical(fit$diagnostics$conv_tier, "ok") && !identical(as.integer(fit$opt$convergence), 0L)) {
    fit$opt$convergence_nlminb <- fit$opt$convergence
    fit$opt$convergence <- 0L
  }


  if(FLRreport){fit$report <- as_FLQuants(fit,stklen)} else {fit$report$FLReport=FALSE}
  return(fit)

}

#' Fit FLicc length-based model using TMB
#'
#' Fits the FLicc length-based equilibrium model to length-frequency data
#' using Template Model Builder (TMB). The model separates the population
#' process from the observation process, allowing flexible combinations of
#' population structure and likelihood assumptions.
#'
#' @param lfd A named list of \code{FLQuant} objects containing observed
#'   length-frequency data by gear. Values can be counts or proportions.
#'
#' @param stklen An \code{FLStockLen}-like object containing life-history
#'   parameters (e.g. Linf, k, M) and matching dimensions to \code{lfd}.
#'
#' @param sel_fun Character vector (named by gear, or in the order of
#'   \code{lfd}) giving the selectivity function of each gear; parameters are
#'   on the length scale:
#'   \describe{
#'     \item{\code{"logistic"}}{Asymptotic: \code{SL50}, \code{SL95}.}
#'     \item{\code{"dsnormal"}}{Double-sided normal (dome): \code{mode},
#'       \code{lsd} (left), \code{rsd} (right); \code{rsd = 1000} in
#'       \code{sel_fixed} makes it asymptotic.}
#'     \item{\code{"normal"}}{Symmetric dome: \code{mode}, \code{sd}.}
#'     \item{\code{"binormal"}}{Two size modes, the sum of two normal curves
#'       scaled to a maximum of 1: \code{mode1}, \code{sd1} (first curve),
#'       \code{mode2 > mode1}, \code{sd2} (second curve) and \code{h}, the
#'       height of the second curve relative to the first. For gears whose
#'       catch has two size groups, e.g. juveniles and adults, or a gillnet
#'       fleet with two mesh groups.}
#'     \item{\code{"bidsnormal"}}{Two size modes with a double-sided second
#'       curve: \code{mode1}, \code{sd1}, \code{mode2}, \code{lsd2},
#'       \code{rsd2}, \code{h}. The right limb of the large-fish mode is
#'       separate; \code{rsd2 = 1000} in \code{sel_fixed} gives an asymptotic
#'       second mode.}
#'   }
#'   Use \code{plot_selfun()} to look at the shapes over the data before
#'   fitting. The two-mode functions have five or six parameters and need
#'   enough data in both modes.
#'
#' @param catch_by_gear Numeric vector of relative catch proportions by gear,
#'   in weight (e.g. landings in tonnes) unless \code{settings$catch_units = "numbers"}.
#'   Must sum to 1.
#'
#' @param settings Named list of model settings. Key options include:
#'   \describe{
#'     \item{\code{pop_model}}{Population model. Options are:
#'       \code{"gtg"} (default), a growth-type-group approximation inspired
#'       by LBSPR, and \code{"gamma"}, the original FLicc/BLICC continuous
#'       growth formulation.}
#'     \item{\code{obs_model}}{Observation likelihood. Options are:
#'       \code{"dm"} (default, Dirichlet-multinomial, with the precision
#'       estimated by default, see \code{dm_theta}; give raw counts or an
#'       effective sample size of at least 500), \code{"mn"} (multinomial, as
#'       LBSPR) and \code{"nb"} (negative binomial, as fishblicc).}
#'     \item{\code{CVL}}{Coefficient of variation in length-at-age.}
#'     \item{\code{GL}}{Number of Gauss–Laguerre quadrature nodes (used only
#'       when \code{pop_model = "gamma"}).}
#'     \item{\code{catch.sd}}{Standard deviation for catch penalty.}
#'     \item{\code{catch_units}}{Units of \code{catch_by_gear}: \code{"weight"}
#'       (default; predicted catch-at-length times weight at length is compared
#'       with the input shares, e.g. landings in tonnes) or \code{"numbers"}
#'       (predicted catch in numbers, as in fishblicc).}
#'     \item{\code{dm_theta}}{Precision of the Dirichlet-multinomial
#'       (\code{obs_model = "dm"}): \code{"saturating"} (default;
#'       \eqn{\alpha_0 = \beta} estimated; the effective
#'       sample size \eqn{n(1+\beta)/(n+\beta)} levels off at \eqn{\beta + 1}, so
#'       the fit becomes insensitive to the input sample size),
#'       \code{"linear"} (\eqn{\alpha_0 = \theta n}, \eqn{\theta} estimated) or
#'       \code{"fixed"} (\eqn{\alpha_0 = n}, effective sample size about
#'       \eqn{n/2}). Thorson et al. (2017, Fisheries Research 192: 84-93).}
#'     \item{\code{dm_shared}}{Logical; one precision parameter for all gears
#'       instead of one per gear (default \code{FALSE}).}
#'     \item{\code{dm_prior_sd}}{SD of a weak normal prior on the log precision
#'       (default 1; 0 = none), centred on \eqn{\beta = 10\max n} (saturating)
#'       or \eqn{\theta = 10} (linear), i.e. close to a multinomial. It only
#'       matters when the data show no overdispersion (input sample size below
#'       what the data support), where it keeps the precision finite.}
#'     \item{\code{pile}}{GTG only: predicted length composition from numbers
#'       integrated over each length class (\code{1}, default) or from the
#'       survivorship at the class bound (\code{0}, LBB without the pile-up
#'       correction, Pile = 0); values in between mix the two as LBB.
#'       Per-recruit quantities (SPR, EPR) always use the integrated numbers.}
#'     \item{\code{ngtg}}{Number of growth-type groups for
#'       \code{pop_model = "gtg"}. Default is 13, as in LBSPR; with \code{maxsd = 3}
#'       more groups (e.g. 19) give practically the same fit at higher cost.}
#'     \item{\code{maxsd}}{Range (in standard deviations) of the GTG Linf
#'       distribution. Default is 3 (LBSPR: 2). The groups cover
#'       \code{Linf * (1 +/- maxsd * CVL)}: no fish in the model grows beyond
#'       \code{Linf * (1 + maxsd * CVL)}, and observed fish above it get a
#'       near-zero probability. Truncation also narrows the actual Linf spread
#'       (0.99 * CVL at 3, about 0.88-0.91 * CVL at 2). With
#'       \code{maxsd = 2}, a few large fish can cause steps in the likelihood
#'       when Linf or CVL change (sensitivity runs, Monte Carlo). Use
#'       \code{maxsd = 2} to reproduce LBSPR. Default changed from 2 in FLicc
#'       1.0.11.}
#'     \item{\code{tail_warning}}{Logical. If \code{TRUE}, warn when observed
#'       fish lie in length bins at or above the largest GTG Linf, where the
#'       model gives them near-zero probability. Default \code{TRUE} (from FLicc 1.0.11).}
#'     \item{\code{Mpow}}{Optional length-dependent mortality scaling
#'       exponent (GTG only). Default is 0.}
#'     \item{\code{sigmaF}, default}{By default log F follows a random walk by gear
#'       with sigmaF fixed at 0.3 (fast and stable). sigmaF is estimated only
#'       when \code{prior_sigmaF} is given (penalised fixed effect with that
#'       prior) or with \code{F_re = TRUE} (random effects). The random walk is
#'       switched off with \code{prior_sigmaF = FALSE}.}
#'     \item{\code{prior_sigmaF}}{Random walk on log F by gear and its
#'       standard deviation sigmaF: \code{c(mean, sd, use)}, a normal prior on
#'       log(sigmaF) with \code{use = 1} switching the random walk on. With
#'       \code{sd = NA} there is no prior on sigmaF.}
#'     \item{\code{F_re}}{Logical. If \code{TRUE}, log F by year and gear are
#'       random effects (state-space random walk) integrated out by the Laplace
#'       approximation, so sigmaF is estimated from the marginal likelihood and
#'       standard errors include the F process. If \code{FALSE} (default), log F
#'       are fixed effects with a random-walk penalty, where sigmaF is mainly
#'       set by its prior. With \code{F_re = TRUE} and no \code{prior_sigmaF}
#'       supplied, a weak prior is used: \code{c(log(0.3), 0.5, 1)}.}
#'     \item{\code{sigmaF}}{Optional fixed value of sigmaF (random-walk SD of
#'       log F). Fits are then as fast as the penalised default. A practical
#'       use is to estimate sigmaF once with \code{F_re = TRUE} on the full
#'       data, then fix it for repeated fits (e.g. in an MSE); the yearly
#'       estimates closely match the full random-effects fit.}
#'     \item{\code{rob_eps}}{Robustness of the composition likelihood to
#'       sparse or outlying bins. Predicted proportions become
#'       \code{(1 - rob_eps) * p + rob_eps / nbins}, which bounds the cost of a
#'       few fish in bins the model barely reaches (e.g. 0.001). Default 0 (off).}
#'     \item{\code{Lplus}}{Optional plus-group length. Observed and predicted
#'       counts in bins with lower bound \code{>= Lplus} are pooled into one
#'       bin before the likelihood (tail compression). Reduces the influence
#'       of the binning and sampling of the largest fish. Default \code{NULL}
#'       (off). Predicted compositions in the report stay unpooled.}
#'     \item{\code{tail_cut}}{Optional tail cut, an absolute length in the
#'       units of the length bins (as \code{Lplus}), e.g. \code{50}. Length
#'       bins with lower bound \code{>= tail_cut} are excluded from the
#'       composition likelihood: observed and predicted compositions are
#'       truncated and the predicted proportions renormalised over the retained
#'       bins, so the fit conditions on fish below the cut. Reduces the
#'       influence of the largest fish, which are most sensitive to
#'       \eqn{L_\infty}, growth variability (CVL), dome-shaped selectivity and
#'       the memory of past F. As a fixed length, the same bins are used
#'       whatever \eqn{L_\infty} is assumed (sensitivity runs, MSE
#'       iterations). \code{LBIspr()} and \code{LBImean()} apply the same cut by
#'       default. Values <= 1.5 are read as a fraction of \eqn{L_\infty} (earlier
#'       definition) with a warning. Default \code{NULL} (off).}
#'     \item{\code{spr_ref}}{SPR target(s) in percent (default \code{40}).
#'       The fit reports the corresponding Fspr (\code{fit$report$Fspr}) and
#'       per-recruit numbers-at-length, computed in TMB after fitting at no
#'       cost to the optimisation; \code{fspr_flicc()} and \code{LBIspr()}
#'       reuse them for these targets.}
#'     \item{\code{FM_min}, \code{FM_max}, \code{FMpen_sd}}{Bounds on
#'       apical F/M (F at the most selected length over M at that length) by
#'       year: outside [\code{FM_min}, \code{FM_max}] (defaults 0.05 and 4) a
#'       steep one-sided penalty on log F/M applies; inside there is no
#'       penalty. \code{FMpen_sd > 0} (default 0.2) switches the bounds on,
#'       \code{0} off. The lower bound stops F sliding to 0 (SPR = 1) when no
#'       gear is asymptotic; raise \code{FM_max} (e.g. 20) for heavily fished
#'       stocks if F/M sits at 4.}
#'     \item{\code{linf.sd, Mk.sd, CVL.sd}}{Optional penalties on life-history
#'       parameters. Supplying one makes that parameter estimated (otherwise it
#'       is fixed at its input value). For \code{pop_model = "gtg"}, the
#'       estimated Linf, M/K and CVL all enter GTG survival: the Linf groups are
#'       \code{Linf * (1 + CVL * z)} on a fixed standardised grid \code{z}, and
#'       M/K is \code{Mk * (Linf / L)^Mpow}.}
#'   }
#'
#' @param years Optional vector of years to fit (default all years in \code{lfd}).
#' @param iter Iteration of \code{lfd} and \code{stklen} to fit.
#' @param compile,dll Compile the TMB template (\code{FALSE}, default, uses
#'   the compiled package library \code{dll}).
#' @param silent Logical; suppress TMB tracing output.
#' @param sel_fixed Optional named list (by gear) of selectivity parameters
#'   fixed at given values on the length scale, with the names used by
#'   \code{sel_fun}: \code{SL50}, \code{SL95} (logistic); \code{mode},
#'   \code{lsd}, \code{rsd} (dsnormal); \code{mode}, \code{sd} (normal);
#'   \code{mode1}, \code{sd1}, \code{mode2}, \code{sd2}, \code{h}
#'   (binormal); \code{mode1}, \code{sd1}, \code{mode2}, \code{lsd2},
#'   \code{rsd2}, \code{h} (bidsnormal). E.g.
#'   \code{list(TZGN = list(rsd = 12), TZLN = list(rsd2 = 1000))}. Unnamed
#'   parameters are estimated. \code{selpars_flicc(fit)} gives values in the
#'   same form.
#' @param FLRreport Logical; return reported quantities as \code{FLQuants}.
#' @param safe_fit,refit_GL For \code{pop_model = "gamma"}: refit with
#'   \code{refit_GL} quadrature nodes if the gradient check fails.
#' @param start Optional numeric vector of starting values for the estimated
#'   TMB parameters. A common use is \code{start = fit$opt$par} to restart
#'   from a previous fit and reduce residual gradients. Currently this is used
#'   for ordinary multi-year fits; it is ignored for \code{by_year = TRUE}.
#' @param n_restart Integer. Number of automatic optimizer restarts from the
#'   current optimum. Restarts can reduce residual gradients when \code{nlminb()}
#'   stops by relative convergence before reaching a stationary point. Default
#'   is 3. Use 0 to disable automatic restarts.
#' @param grad_tol Numeric. Maximum absolute gradient tolerance used to stop
#'   automatic restarts and by the optional safe-fit check. Default is
#'   \code{1e-3}.
#' @param control Optional list of control arguments passed to \code{nlminb()}.
#'   Values supplied here override the FLicc defaults.
#' @param se Logical. If \code{FALSE}, \code{TMB::sdreport()} is skipped: no
#'   standard errors, but faster (useful in an MSE). The Hessian is still
#'   checked by the Newton polishing.
#'   Default \code{TRUE}.
#' @param newton Integer. Maximum number of Newton steps after \code{nlminb()}
#'   (default 3; 0 = none). Each step uses the Hessian of the objective and is
#'   kept only if it lowers the objective. This drives the gradient of a
#'   well-defined optimum to near zero, checks that the Hessian is positive
#'   definite and gives the Newton decrement (see \code{\link{flicc_convergence}}).
#' @param verbose Logical. If \code{TRUE} (default), print a message when the
#'   convergence check gives a warning or failure, or when an estimated
#'   Dirichlet-multinomial precision is at its upper bound. The same
#'   information is always stored in \code{fit$diagnostics}.
#' @param by_year Logical. If \code{TRUE}, the model is fitted independently
#'   for each year (LBSPR-style). Only the final fit object is retained, while
#'   annual report quantities are updated and combined across years.
#'   Default is \code{FALSE}.
#'
#' @details
#' The FLicc model separates:
#'
#' \itemize{
#'   \item \strong{Population model} (\code{pop_model}):
#'     \itemize{
#'       \item \code{"gtg"}: discrete growth-type-group approximation where
#'         variability in asymptotic length is represented by a finite set
#'         of groups with associated recruitment proportions.
#'       \item \code{"gamma"}: continuous formulation based on a gamma
#'         distribution of growth, solved using Gauss–Laguerre quadrature.
#'     }
#'
#'   \item \strong{Observation model} (\code{obs_model}):
#'     \itemize{
#'       \item \code{"mn"}: multinomial likelihood for length compositions.
#'         Effective sample size is taken directly from the observed LFD
#'         (sum of counts), consistent with LBSPR and fishblicc.
#'       \item \code{"nb"}: negative binomial likelihood for counts with
#'         overdispersion.
#'       \item \code{"dm"} (default): Dirichlet-multinomial likelihood for
#'         overdispersed compositions, with the precision estimated
#'         (\code{dm_theta = "saturating"}) so that the data set their own
#'         weight; use raw counts or an effective sample size of at least 500.
#'     }
#' }
#'
#' The default model (\code{pop_model = "gtg"}, \code{obs_model = "dm"} with
#' estimated precision, all years jointly with a random walk on F) is FLicc's
#' recommended configuration. LBSPR corresponds to \code{obs_model = "mn"} with
#' Linf, M/k and CVL fixed and each year fitted separately
#' (\code{by_year = TRUE}); LBB to \code{lbb_flicc()}; fishblicc to
#' \code{pop_model = "gamma"}, \code{obs_model = "nb"}.
#'
#' @return An object of class \code{"flicc_tmb_fit"} containing:
#' \describe{
#'   \item{\code{opt}}{Optimizer output.}
#'   \item{\code{report}}{Reported model quantities from TMB.}
#'   \item{\code{par}}{Estimated parameters.}
#'   \item{\code{tmb_data}}{Input data passed to TMB.}
#'   \item{\code{pop_model}}{Population model used.}
#'   \item{\code{obs_model}}{Observation model used.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("logistic", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#' }
#'
#' @export
fiticc <- function(lfd, stklen,
                   sel_fun = c("logistic", "normal", "dsnormal"),
                   catch_by_gear,
                   settings = list(
                     CVL = 0.1,
                     GL = 30,
                     catch.sd = 0.05,
                     pop_model = "gtg",
                     obs_model = "dm",
                     ngtg = 13,
                     maxsd = 3,
                     Mpow = 0,
                     linf.sd = NULL,
                     Mk.sd = NULL,
                     CVL.sd = NULL
                   ),
                   years = NULL,
                   iter = 1,
                   compile = FALSE,
                   silent = TRUE,
                   dll = "FLicc",
                   start = NULL,
                   by_year = FALSE,
                   sel_fixed = NULL,
                   FLRreport = TRUE,
                   safe_fit = TRUE,
                   refit_GL = 100,
                   n_restart = 3,
                   grad_tol = 1e-3,
                   control = list(),
                   se = TRUE,
                   newton = 3,
                   verbose = TRUE) {


     if(by_year){
      if (!is.null(start)) {
        warning("Argument 'start' is currently ignored when by_year = TRUE.")
      }
      yrs <- dimnames(stklen)$year
      if(!is.null(years)){
       yrs <-  yrs[an(yrs)%in%years]
      }

      nyrs <- length(yrs)
      # Loops through years
      fity <- lapply(yrs,function(y){

           if (inherits(catch_by_gear, "FLQuants")){
             cbg <- window(catch_by_gear,start= yrs[nyrs],end=yrs[nyrs])
           } else {
             cbg <- catch_by_gear
           }

         fiticc_core(
         lfd = window(lfd,start= y,end=y),
         stklen = window(stklen,start= y,y),
         sel_fun = sel_fun,
         catch_by_gear =cbg,
         settings = settings,
         years = NULL,
         iter = iter,
         compile = compile,
         silent = silent,
         dll = dll,
         sel_fixed = sel_fixed,
         n_restart = n_restart,
         se = se,
         grad_tol = grad_tol,
         control = control,
         FLRreport = TRUE,
         newton = newton
         )}
      )
      repy <-  fity[[1]]$report
      fit_final = fity[[nyrs]]
      dn <-  names(fity[[1]]$report)[-c(3)]
      qs <- dn[c(1:14)]

      for(y in 2:nyrs){
       for(i in seq(qs)){
        flq <- repy[[qs[i]]]
        if(!inherits(repy[[qs[i]]], "FLQuants")){
        repy[[qs[i]]] <- append(repy[[qs[i]]],fity[[y]]$report[[qs[i]]])
        } else {
        for(f in 1:length(repy[[qs[i]]])){
          repy[[qs[i]]][[f]] <- append(repy[[qs[i]]][[f]],fity[[y]]$report[[qs[i]]][[f]])
        }}
       }
      }
      fit_final <- fity[[nyrs]]
      fit_final$stklen <- stklen
      fit_final$report[qs] <- repy[qs]
      if (!is.null(fity[[1]]$report$epr)) {
        ep <- fity[[1]]$report$epr
        for (y in 2:nyrs) ep <- append(ep, fity[[y]]$report$epr)
        fit_final$report$epr <- ep
      }
  }


  if(!by_year){
  fit1 <- fiticc_core(
    lfd = lfd,
    stklen = stklen,
    sel_fun = sel_fun,
    catch_by_gear = catch_by_gear,
    settings = settings,
    years = years,
    iter = iter,
    compile = compile,
    silent = silent,
    dll = dll,
    sel_fixed = sel_fixed,
    start = start,
    n_restart = n_restart,
    se = se,
    grad_tol = grad_tol,
    control = control,
    FLRreport = FLRreport,
    newton = newton
  )

  if(pop_model_flicc(fit1)%in%"gtg")  safe_fit = FALSE

  refit_used <- FALSE
  fit_final <- fit1

  if (isTRUE(safe_fit) && need_refit_flicc(fit1, grad_tol = grad_tol)) {

    settings2 <- settings
    settings2$GL <- refit_GL

    fit2 <- tryCatch(
      fiticc_core(
        lfd = lfd,
        stklen = stklen,
        sel_fun = sel_fun,
        catch_by_gear = catch_by_gear,
        settings = settings2,
        years = years,
        iter = iter,
        compile = compile,
        silent = silent,
        dll = dll,
        sel_fixed = sel_fixed,
        start = fit1$opt$par,
        n_restart = n_restart,
        se = se,
        grad_tol = grad_tol,
        control = control,
        FLRreport = FLRreport,
        newton = newton
      ),
      error = function(e) NULL
    )

    if (!is.null(fit2)) {
      fit_final <- fit2
      refit_used <- TRUE
    }
  }

  fit_final$safe_fit <- list(
    enabled = safe_fit,
    refit_used = refit_used,
    initial_GL = settings$GL,
    refit_GL = if (refit_used) refit_GL else NA_real_,
    final_needs_refit = need_refit_flicc(fit_final, grad_tol = grad_tol)
  )
  }

  if (isTRUE(verbose)) messages_flicc(fit_final)
  fit_final
}


#' Check whether a FLicc fit should be refitted
#'
#' Applies a small set of basic convergence checks to a fitted FLicc model and
#' returns whether a conditional refit is recommended.
#'
#' @param fit A fitted object of class \code{"flicc_tmb_fit"}.
#'
#' @details
#' The current checks are intentionally simple and are designed for stable
#' routine use in \code{fiticc()}. A refit is recommended when:
#' \itemize{
#'   \item the fit object is \code{NULL},
#'   \item the optimizer convergence code is not zero,
#'   \item the objective function is not finite, or
#'   \item the Hessian reported by \code{TMB::sdreport()} is not positive
#'   definite.
#' }
#'
#' These checks are used by \code{fiticc()} to decide whether to rerun the model
#' once with a higher \code{GL} value.
#'
#' @return Logical; \code{TRUE} if a refit is recommended, otherwise
#'   \code{FALSE}.
#'
#' @seealso \code{\link{fiticc}}, \code{\link{fiticc_core}}
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#'
#' need_refit_flicc(fit)

need_refit_flicc <- function(fit, grad_tol = 0.02) {

  if (is.null(fit))
    return(TRUE)

  if (is.null(fit$opt$convergence) || fit$opt$convergence != 0)
    return(TRUE)

  if (is.null(fit$opt$objective) || !is.finite(fit$opt$objective))
    return(TRUE)

  if (!is.null(fit$rep$pdHess) && !isTRUE(fit$rep$pdHess))
    return(TRUE)

  # Hessian check from the Newton polishing (available without sdreport)
  if (!is.null(fit$diagnostics$hess_pd) && isFALSE(fit$diagnostics$hess_pd))
    return(TRUE)



  FALSE
}

#' Apply fixed selectivity settings to FLicc parameters
#'
#' Updates the packed selectivity parameter vector \code{Sm} and builds a TMB
#' mapping object for fixing selected selectivity parameters by gear. This
#' provides a user-facing interface for constraining selectivity parameters on
#' the natural length scale while retaining the internal relative-to-\code{Linf}
#' parameterization used by the TMB model.
#'
#' @param parameters A named parameter list, typically produced by
#'   \code{init_tmb_flicc()}, containing at least the packed selectivity vector
#'   \code{Sm}.
#' @param tmb_data A TMB data list produced by \code{data_tmb_flicc()}.
#' @param sel_fixed An optional named list specifying fixed selectivity values by
#'   gear. Names may be either gear names or character representations of gear
#'   indices. Each element should itself be a list of fixed parameter values on
#'   the natural length scale.
#'
#' @details
#' The function allows users to fix selectivity parameters by gear without
#' needing to work directly with the internal packed \code{Sm} vector.
#'
#' Supported user-facing selectivity inputs are:
#'
#' \strong{Logistic selectivity} (\code{sel_type == 1})
#' \itemize{
#'   \item \code{SL50}: length at 50 percent selectivity
#'   \item \code{SL95}: length at 95 percent selectivity
#'   \item \code{dSL}: difference \code{SL95 - SL50}
#' }
#'
#' \strong{Normal selectivity} (\code{sel_type == 3})
#' \itemize{
#'   \item \code{mode}: modal length of selectivity
#'   \item \code{sd}: standard deviation of the symmetric dome
#' }
#'
#' \strong{Double-sided normal selectivity} (\code{sel_type == 2})
#' \itemize{
#'   \item \code{mode}: modal length of selectivity
#'   \item \code{lsd}: left-side standard deviation
#'   \item \code{rsd}: right-side standard deviation
#' }
#'
#' Internally, all selectivity parameters are stored relative to \code{Linf} on
#' the log scale. The function converts natural-scale user inputs to this
#' internal parameterization and updates the TMB map so that selected values are
#' held fixed during optimization.
#'
#' This makes it possible, for example, to:
#' \itemize{
#'   \item fix a logistic selectivity curve using \code{SL50} and \code{SL95},
#'   \item fix a symmetric dome using \code{mode} and \code{sd},
#'   \item or approximate a one-sided shoulder by fixing a very large
#'     \code{rsd} for a \code{dsnormal} gear.
#' }
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{parameters}}{The updated parameter list with modified
#'       selectivity values in \code{Sm}.}
#'     \item{\code{map_Sm}}{A factor vector suitable for inclusion in the TMB
#'       \code{map} argument. Elements set to \code{NA} are fixed during
#'       optimization. Returns \code{NULL} if no selectivity parameters are
#'       fixed.}
#'   }
#'
#' @examples
#' \dontrun{
#' tmb_data <- data_tmb_flicc(dl)
#' parameters <- init_tmb_flicc(tmb_data)
#'
#' # Fix logistic selectivity using SL50 and SL95
#' sel_fixed <- list(
#'   "Gill net" = list(SL50 = 22, SL95 = 27)
#' )
#'
#' fix_res <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed)
#'
#' # Fix right tail of a double-sided normal to be very broad
#' sel_fixed2 <- list(
#'   "Gill net" = list(rsd = 1000)
#' )
#'
#' fix_res2 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed2)
#'
#' # Fix a symmetric normal dome
#' sel_fixed3 <- list(
#'   "Marine set bagnet" = list(mode = 24, sd = 4)
#' )
#'
#' fix_res3 <- apply_sel_fixed_flicc(parameters, tmb_data, sel_fixed3)
#' }
#'
#' @export

apply_sel_fixed_flicc <- function(parameters, tmb_data, sel_fixed = NULL) {
  if (is.null(sel_fixed) || length(sel_fixed) == 0) {
    return(list(parameters = parameters, map_Sm = NULL))
  }

  Sm_map <- factor(seq_along(parameters$Sm))

  gear_names <- tmb_data$gear_names
  sel_type   <- tmb_data$sel_type
  sm_start   <- tmb_data$sm_start
  sm_n       <- tmb_data$sm_n
  Linf       <- tmb_data$Linf_init

  for (nm in names(sel_fixed)) {
    spec <- sel_fixed[[nm]]
    if (is.null(spec)) next

    if (nm %in% gear_names) {
      g <- match(nm, gear_names)
    } else if (grepl("^[0-9]+$", nm)) {
      g <- as.integer(nm)
      if (g < 1 || g > length(gear_names)) {
        stop("Gear index out of range in sel_fixed: ", nm)
      }
    } else {
      stop("Unknown gear in sel_fixed: ", nm)
    }

    start <- sm_start[g] + 1L

    if (sel_type[g] == 1L) {
      # logistic: c(log_SL50_rel, log_dSL_rel)

      current_SL50 <- exp(parameters$Sm[start]) * Linf
      current_dSL  <- exp(parameters$Sm[start + 1L]) * Linf

      if (!is.null(spec$SL50)) {
        current_SL50 <- as.numeric(spec$SL50)
        if (current_SL50 <= 0) stop("SL50 must be > 0")
        parameters$Sm[start] <- log(current_SL50 / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$SL95)) {
        dSL <- as.numeric(spec$SL95) - current_SL50
        if (dSL <= 0) stop("For logistic selectivity, SL95 must be > SL50")
        parameters$Sm[start + 1L] <- log(dSL / Linf)
        Sm_map[start + 1L] <- NA
      }

      if (!is.null(spec$dSL)) {
        dSL <- as.numeric(spec$dSL)
        if (dSL <= 0) stop("For logistic selectivity, dSL must be > 0")
        parameters$Sm[start + 1L] <- log(dSL / Linf)
        Sm_map[start + 1L] <- NA
      }

    } else if (sel_type[g] == 2L) {
      # dsnormal: c(log_mode_rel, log_lsd_rel, log_rsd_rel)

      if (!is.null(spec$mode)) {
        mode <- as.numeric(spec$mode)
        if (mode <= 0) stop("mode must be > 0")
        parameters$Sm[start] <- log(mode / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$lsd)) {
        lsd <- as.numeric(spec$lsd)
        if (lsd <= 0) stop("lsd must be > 0")
        parameters$Sm[start + 1L] <- log(lsd / Linf)
        Sm_map[start + 1L] <- NA
      }

      if (!is.null(spec$rsd)) {
        rsd <- as.numeric(spec$rsd)
        if (rsd <= 0) stop("rsd must be > 0")
        parameters$Sm[start + 2L] <- log(rsd / Linf)
        Sm_map[start + 2L] <- NA
      }

    } else if (sel_type[g] == 3L) {
      # normal: c(log_mode_rel, log_sd_rel)

      if (!is.null(spec$mode)) {
        mode <- as.numeric(spec$mode)
        if (mode <= 0) stop("mode must be > 0")
        parameters$Sm[start] <- log(mode / Linf)
        Sm_map[start] <- NA
      }

      if (!is.null(spec$sd)) {
        sd <- as.numeric(spec$sd)
        if (sd <= 0) stop("sd must be > 0")
        parameters$Sm[start + 1L] <- log(sd / Linf)
        Sm_map[start + 1L] <- NA
      }

    } else if (sel_type[g] == 4L) {
      # binormal: c(log_mode1_rel, log_sd1_rel, log_dmode_rel, log_sd2_rel, log_h)
      cur_m1 <- exp(parameters$Sm[start]) * Linf
      cur_m2 <- cur_m1 + exp(parameters$Sm[start + 2L]) * Linf
      if (!is.null(spec$mode1)) {
        cur_m1 <- as.numeric(spec$mode1)
        if (cur_m1 <= 0) stop("mode1 must be > 0")
        parameters$Sm[start] <- log(cur_m1 / Linf)
        Sm_map[start] <- NA
      }
      if (!is.null(spec$sd1)) {
        if (spec$sd1 <= 0) stop("sd1 must be > 0")
        parameters$Sm[start + 1L] <- log(as.numeric(spec$sd1) / Linf)
        Sm_map[start + 1L] <- NA
      }
      if (!is.null(spec$mode2)) {
        dm <- as.numeric(spec$mode2) - cur_m1
        if (dm <= 0) stop("For binormal selectivity, mode2 must be > mode1")
        parameters$Sm[start + 2L] <- log(dm / Linf)
        Sm_map[start + 2L] <- NA
      }
      if (!is.null(spec$sd2)) {
        if (spec$sd2 <= 0) stop("sd2 must be > 0")
        parameters$Sm[start + 3L] <- log(as.numeric(spec$sd2) / Linf)
        Sm_map[start + 3L] <- NA
      }
      if (!is.null(spec$h)) {
        if (spec$h <= 0) stop("h must be > 0")
        parameters$Sm[start + 4L] <- log(as.numeric(spec$h))
        Sm_map[start + 4L] <- NA
      }

    } else if (sel_type[g] == 5L) {
      # bidsnormal: c(log_mode1_rel, log_sd1_rel, log_dmode_rel, log_lsd2_rel, log_rsd2_rel, log_h)
      cur_m1 <- exp(parameters$Sm[start]) * Linf
      if (!is.null(spec$mode1)) {
        cur_m1 <- as.numeric(spec$mode1)
        if (cur_m1 <= 0) stop("mode1 must be > 0")
        parameters$Sm[start] <- log(cur_m1 / Linf)
        Sm_map[start] <- NA
      }
      fixpos <- function(val, off, nm, f = function(v) log(v / Linf)) {
        if (is.null(val)) return(invisible())
        if (val <= 0) stop(nm, " must be > 0")
        parameters$Sm[start + off] <<- f(as.numeric(val))
        Sm_map[start + off] <<- NA
      }
      fixpos(spec$sd1, 1L, "sd1")
      if (!is.null(spec$mode2)) {
        dm <- as.numeric(spec$mode2) - cur_m1
        if (dm <= 0) stop("For bidsnormal selectivity, mode2 must be > mode1")
        parameters$Sm[start + 2L] <- log(dm / Linf)
        Sm_map[start + 2L] <- NA
      }
      fixpos(spec$lsd2, 3L, "lsd2")
      fixpos(spec$rsd2, 4L, "rsd2")
      fixpos(spec$h, 5L, "h", function(v) log(v))

    } else {
      stop("Unsupported selectivity type for gear ", g)
    }

    if ((start + sm_n[g] - 1L) > length(parameters$Sm)) {
      stop("Selectivity indexing exceeded Sm length for gear ", g)
    }
  }

  # drop the levels of fixed entries: TMB sizes the parameter vector by the
  # number of levels, so unused levels would add phantom NA parameters
  Sm_map <- factor(as.character(Sm_map), levels = unique(as.character(Sm_map[!is.na(Sm_map)])))
  list(parameters = parameters, map_Sm = Sm_map)
}


#' Summarise convergence diagnostics for a FLicc fit
#'
#' Provides a compact, lay-readable convergence summary for a fitted FLicc/TMB
#' model. The function checks whether the optimiser terminated normally, whether
#' the maximum absolute gradient is small, whether the Hessian is positive
#' definite, whether standard errors are available, and whether the numerical
#' Hessian has positive eigenvalues.
#'
#' @param fit A fitted FLicc model object, typically returned by `fiticc()`.
#'   The object is expected to contain `fit$opt`, `fit$obj`, `fit$par`, and
#'   either `fit$rep$pdHess` or `fit$sdrep$pdHess`.
#' @param grad_tol Numeric. Tolerance used to judge whether the maximum absolute
#'   gradient is sufficiently small. Default is `1e-3`.
#' @param eig_tol Numeric. Tolerance used to judge whether Hessian eigenvalues
#'   are positive. Default is `1e-8`.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{check}{Name of the convergence diagnostic.}
#'   \item{value}{Observed diagnostic value.}
#'   \item{passed}{Logical flag indicating whether the check passed.}
#'   \item{interpretation}{Plain-language explanation of the check.}
#' }
#'
#' @details
#' Optimiser convergence alone does not guarantee a reliable fit. A model may
#' stop with optimiser code `0` but still have a large gradient, indicating that
#' the solution is not close to a stationary point. Conversely, standard errors
#' and a positive-definite Hessian indicate that local curvature is estimable,
#' but should still be interpreted together with the gradient check.
#'
#' The Hessian positive-definiteness check first looks for `fit$rep$pdHess`,
#' then for `fit$sdrep$pdHess`. The eigenvalue check recomputes a numerical
#' Hessian using `stats::optimHess()`.
#'
#' The Newton decrement \eqn{\frac12 g'H^{-1}g} (from the Newton polishing in
#' [fiticc()]) is the expected further improvement in log-likelihood and does
#' not depend on the scale of the objective or of the parameters, unlike the
#' gradient. With a positive definite Hessian and a decrement available, the
#' overall tier is "ok" when the decrement is below 0.001 and "warning" when it
#' is up to 0.1, whatever the gradient; without a decrement (`newton = 0`) the
#' gradient is compared with `grad_tol` instead. The tier is "fail" when the
#' Hessian is not positive definite or the decrement is above 0.1.
#'
#' @examples
#' \dontrun{
#' fit <- fiticc(lfds, stock, lhpar)
#' flicc_convergence(fit)
#' }
#'
#' @export
flicc_convergence <- function(fit, grad_tol = 1e-2, eig_tol = 1e-8,optimizer.code=FALSE) {

  opt_code <- fit$opt$convergence
  opt_msg  <- fit$opt$message

  max_grad <- if (!is.null(fit$opt$gradient)) {
    max(abs(fit$opt$gradient), na.rm = TRUE)
  } else if (!is.null(fit$obj$gr)) {
    max(abs(fit$obj$gr(fit$opt$par)), na.rm = TRUE)
  } else {
    NA_real_
  }

  pdHess <- if (!is.null(fit$rep$pdHess)) {
    fit$rep$pdHess
  } else if (!is.null(fit$sdrep$pdHess)) {
    fit$sdrep$pdHess
  } else {
    NA
  }

  n_par <- if (!is.null(fit$par)) nrow(fit$par) else length(fit$opt$par)

  n_nan_se <- if (!is.null(fit$par) && "se" %in% names(fit$par)) {
    sum(is.nan(fit$par$se))
  } else {
    NA_integer_
  }

  nan_pars <- if (!is.null(fit$par) && "se" %in% names(fit$par)) {
    unique(fit$par$par[is.nan(fit$par$se)])
  } else {
    character(0)
  }

  eig <- try({
    H <- stats::optimHess(fit$opt$par, fit$obj$fn, fit$obj$gr)
    eigen(H, symmetric = TRUE, only.values = TRUE)$values
  }, silent = TRUE)

  min_eig <- if (!inherits(eig, "try-error")) min(eig, na.rm = TRUE) else NA_real_

  dg <- fit$diagnostics
  dec <- if (!is.null(dg$decrement)) dg$decrement else NA_real_
  tier <- conv_tier_flicc(fit$opt$objective,
                          if (is.na(pdHess)) (if (is.finite(min_eig)) min_eig > eig_tol else NA) else isTRUE(pdHess),
                          dec, max_grad, grad_tol = grad_tol)

  out <- data.frame(
    check = c(
      "Optimizer finished normally",
      "Gradient is acceptable",
      "Hessian is positive definite",
      "Standard errors available",
      "Hessian eigenvalues are positive",
      "Newton decrement is negligible",
      "Overall"
    ),
    value = c(
      paste0(opt_code, " - ", opt_msg),
      signif(max_grad, 4),
      as.character(pdHess),
      paste0(n_par - n_nan_se, " / ", n_par, " parameters"),
      signif(min_eig, 4),
      signif(dec, 3),
      tier
    ),
    passed = c(
      identical(opt_code, 0L),
      isTRUE(max_grad < grad_tol),
      isTRUE(pdHess),
      isTRUE(n_nan_se == 0),
      isTRUE(min_eig > eig_tol),
      if (is.na(dec)) NA else dec <= 1e-3,
      tier != "fail"
    ),
    interpretation = c(
      if (identical(opt_code, 0L)) "The optimizer found a numerical optimum." else
        "The optimizer stopped without reporting convergence.",
      if (isTRUE(max_grad < grad_tol)) "The fit is close to a stationary point." else
        "The gradient is not small: the fit may not be at an optimum.",
      if (isTRUE(pdHess)) "The model curvature is well behaved." else
        "The Hessian is not positive definite: a saddle point or a flat direction.",
      if (isTRUE(n_nan_se == 0)) "Uncertainty could be estimated for all parameters." else
        paste(n_nan_se, "parameter(s) without standard error: see fit$par[!is.finite(fit$par$se), ]."),
      if (isTRUE(min_eig > eig_tol)) "No flat or negative curvature directions were detected." else
        "Flat or negative curvature: some parameter(s) are not identified by the data.",
      "Expected further gain in log-likelihood (0.5 g'H^-1 g) is below 0.001; NA without Newton polishing.",
      "ok: converged; warning: usable, but check the gradient or decrement; fail: do not use."
    ),
    stringsAsFactors = FALSE
  )

  if(!optimizer.code) out <- out[-1,]
  out
}

#' Random-walk standard deviation of log F from an FLicc fit
#'
#' Returns sigmaF, the standard deviation of the random walk on log F, as
#' estimated (or fixed) in a fit. Typical use: estimate it once with
#' \code{settings$F_re = TRUE}, then fix it in repeated fits (e.g. in an MSE)
#' with \code{settings$sigmaF = sigmaF_flicc(fit)}.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @return Numeric scalar sigmaF, or \code{NA} if the fit has no random walk on F.
#' @examples
#' \dontrun{
#' fit_re <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
#'                  settings = list(F_re = TRUE))
#' sig <- sigmaF_flicc(fit_re)
#' fit <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
#'               settings = list(sigmaF = sig))
#' }
#' @export
sigmaF_flicc <- function(fit) {
  if (!is.null(fit$par) && "log_sigmaF" %in% fit$par$par) {
    return(exp(fit$par$mpd[fit$par$par == "log_sigmaF"][1]))
  }
  if (!is.null(fit$settings$sigmaF)) return(as.numeric(fit$settings$sigmaF))
  pl <- tryCatch(fit$obj$env$parList(par = fit$obj$env$last.par.best), error = function(e) NULL)
  if (!is.null(pl) && isTRUE(fit$tmb_data$prior_sigmaF_use == 1L) &&
      length(unique(fit$tmb_data$year_names)) > 1) return(exp(pl$log_sigmaF))
  NA_real_
}
