#' Extract FLicc population model name
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#'
#' @return Character string, one of \code{"gamma"} or \code{"gtg"}.
#' @export
gtg_inputs_flicc <- function(fit) {
  # GTG Linf groups, recruitment weights and M/K at the estimates.
  # Taken from the fit report (FLicc >= 1.0.6), else from the TMB object,
  # else from tmb_data for fits made with earlier versions.
  pick <- function(src) {
    if (!is.null(src$gtgLinfs) && !is.null(src$recP) && !is.null(src$MKMat)) {
      list(gtgLinfs = as.numeric(src$gtgLinfs),
           recP     = as.numeric(src$recP),
           MKMat    = as.matrix(src$MKMat))
    } else NULL
  }
  out <- pick(fit$report)
  if (is.null(out) && !is.null(fit$obj))
    out <- pick(fit$obj$report(fit$obj$env$last.par.best))
  if (is.null(out)) out <- pick(fit$tmb_data)
  if (is.null(out)) stop("GTG inputs (gtgLinfs, recP, MKMat) not found in fit.")
  out
}

pop_model_flicc <- function(fit) {

  if (!is.null(fit$pop_model)) {
    return(fit$pop_model)
  }

  if (!is.null(fit$settings$pop_model)) {
    return(fit$settings$pop_model)
  }

  code <- fit$tmb_data$pop_model

  if (is.null(code)) {
    return("gamma")
  }

  switch(as.character(code),
         "1" = "gamma",
         "2" = "gtg",
         stop("Unknown pop_model code in fit$tmb_data$pop_model")
  )
}


#' FLicc helper functions for equilibrium numbers-, catch-, and per-recruit analyses
#'
#' Helper functions for extracting fitted reference quantities from a
#' \code{"flicc_tmb_fit"} object and computing equilibrium numbers-at-length,
#' catch-at-length, and per-recruit quantities under arbitrary fishing
#' intensity.
#'
#' The functions in this file are intended as building blocks for reference
#' point calculations such as YPR, SBPR, SPR, \code{Fspr40}, and
#' length-based indicators such as \code{LBIspr}. They mirror the core
#' FLicc/TMB population recursion while remaining usable directly from R.
#'
#' @name fliccit
NULL

#' Equilibrium numbers-at-length per recruit using FLicc recursion
#'
#' Computes equilibrium numbers-at-length per recruit from the FLicc recursion,
#' conditional on length-specific total mortality and growth parameters.
#'
#' @param node Numeric vector of Gauss-Laguerre quadrature nodes.
#' @param quad_wt Numeric vector of Gauss-Laguerre quadrature weights.
#' @param Len Numeric vector of lower-bound or representative length classes.
#' @param Zki Numeric vector of total mortality at length.
#' @param Galpha Positive scalar growth parameter.
#' @param Gbeta Positive scalar growth parameter.
#' @param return_surv Logical. If \code{TRUE}, also returns the survival
#'   function used in the recursion.
#'
#' @return If \code{return_surv = FALSE}, a numeric vector of equilibrium
#'   numbers-at-length per recruit. If \code{TRUE}, a list with elements:
#'   \describe{
#'     \item{\code{NI}}{Numbers-at-length per recruit.}
#'     \item{\code{surv}}{Survival function by length class.}
#'   }
#'
#' @export

n_f_flicc <- function(node, quad_wt, Len, Zki, Galpha, Gbeta,
                      return_surv = FALSE) {

  nv <- length(node)
  LN <- length(Len)

  if (length(quad_wt) != nv) {
    stop("'quad_wt' must have same length as 'node'.")
  }
  if (length(Zki) != LN) {
    stop("'Zki' must have same length as 'Len'.")
  }
  if (any(Zki <= 0)) {
    stop("'Zki' must be strictly positive.")
  }
  if (!is.numeric(Galpha) || length(Galpha) != 1 || Galpha <= 0) {
    stop("'Galpha' must be a positive scalar.")
  }
  if (!is.numeric(Gbeta) || length(Gbeta) != 1 || Gbeta <= 0) {
    stop("'Gbeta' must be a positive scalar.")
  }

  lgamma_Galpha <- lgamma(Galpha)
  Galpha_1 <- Galpha - 1

  surv <- numeric(LN)
  NI <- numeric(LN)
  Zii <- numeric(LN)

  x_beta <- node / Gbeta
  log_x_beta <- log(x_beta)

  Zii[1] <- -Zki[1]
  if (LN >= 3) {
    Zii[2:(LN - 1)] <- Zki[1:(LN - 2)] - Zki[2:(LN - 1)]
  }

  Len1 <- Len[1]
  tmp1 <- node + Gbeta * Len1
  surv[1] <- sum(exp(log(tmp1) * Galpha_1 - Gbeta * Len1 - lgamma_Galpha) * quad_wt)

  if (LN >= 2) {
    for (Li in 2:LN) {
      Ln <- Len[Li]
      acc <- 0

      for (i in seq_len(nv)) {
        lim <- 0

        if (Li >= 2) {
          Lrange <- Ln - Len[1:(Li - 1)]
          lim <- sum(log(x_beta[i] + Lrange) * Zii[1:(Li - 1)])
        }

        acc <- acc + exp(
          lim +
            log_x_beta[i] * Zki[Li - 1] +
            log(node[i] + Gbeta * Ln) * Galpha_1 -
            (Gbeta * Ln + lgamma_Galpha)
        ) * quad_wt[i]
      }

      surv[Li] <- acc
    }
  }

  if (LN >= 2) {
    NI[1:(LN - 1)] <- (surv[1:(LN - 1)] - surv[2:LN]) / Zki[1:(LN - 1)]
  }
  NI[LN] <- surv[LN] / Zki[LN]

  if (return_surv) {
    list(NI = NI, surv = surv)
  } else {
    NI
  }
}


#' Equilibrium numbers-at-length per recruit using GTG (LBSPR) recursion
#'
#' @param Len Numeric vector of lower length bin boundaries.
#' @param gtgLinfs Numeric vector of Linf values for growth-type groups.
#' @param ZKLMat Numeric matrix of total mortality over K, dimension
#'   \code{(length(Len)+1) x length(gtgLinfs)}.
#' @param recP Numeric vector of recruitment proportions for each GTG.
#' @param return_surv Logical. If \code{TRUE}, also returns boundary survival.
#'
#' @return Numeric vector of equilibrium numbers-at-length per recruit, or a
#'   list when \code{return_surv = TRUE}.
#'
#' @export
n_f_gtg_flicc <- function(Len, gtgLinfs, ZKLMat, recP, return_surv = FALSE) {

  nlen <- length(Len)
  ngtg <- length(gtgLinfs)

  if (!is.matrix(ZKLMat)) {
    stop("'ZKLMat' must be a matrix.")
  }
  if (nrow(ZKLMat) != (nlen + 1)) {
    stop("'ZKLMat' must have nlen + 1 rows.")
  }
  if (ncol(ZKLMat) != ngtg) {
    stop("'ZKLMat' must have one column per GTG.")
  }
  if (length(recP) != ngtg) {
    stop("'recP' must have same length as 'gtgLinfs'.")
  }
  if (any(ZKLMat <= 0)) {
    stop("'ZKLMat' must be strictly positive.")
  }

  # build upper boundary
  if (nlen > 1) {
    step_last <- Len[nlen] - Len[nlen - 1]
  } else {
    step_last <- 1
  }
  LBins <- c(Len, Len[nlen] + step_last)

  NPRFished <- matrix(0, nrow = nlen + 1, ncol = ngtg)
  NPRFished[1, ] <- recP

  for (g in seq_len(ngtg)) {
    for (l in 2:(nlen + 1)) {
      num <- gtgLinfs[g] - LBins[l]
      den <- gtgLinfs[g] - LBins[l - 1]

      num <- max(num, 1e-12)
      den <- max(den, 1e-12)

      NPRFished[l, g] <-
        NPRFished[l - 1, g] * (num / den)^ZKLMat[l - 1, g]
    }
  }

  NI <- numeric(nlen)
  for (l in seq_len(nlen)) {
    NI[l] <- sum((NPRFished[l, ] - NPRFished[l + 1, ]) / ZKLMat[l, ])
  }

  if (!return_surv) {
    return(NI)
  }

  list(
    NI = NI,
    surv = NPRFished
  )
}

#' Equilibrium numbers-at-length per recruit from an FLicc fit
#'
#' Computes equilibrium numbers-at-length per recruit from a fitted
#' \code{"flicc_tmb_fit"} object using the FLicc recursion, conditional on
#' fishing mortality-at-length and total mortality-at-length derived from
#' \code{calc_Z_l()}.
#'
#' This is a low-level helper that returns unscaled numbers-per-recruit.
#' Higher-level wrappers such as \code{nf_flicc()} can then apply recruitment
#' scaling (for example via \code{R0}).
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param F Optional numeric scalar apical fishing mortality.
#' @param FM Optional numeric scalar apical \code{F/M} multiplier.
#' @param nyears Integer. Number of terminal years to average over when
#'   constructing biological reference quantities.
#' @param Sel Optional selectivity override passed to \code{calc_Z_l()}.
#'   If \code{NULL}, the terminal-year average selectivity from the fit is used.
#' @param scale_sel Logical. If \code{TRUE}, selectivity is scaled to an
#'   apical maximum of 1 when extracting reference parameters.
#' @param return_surv Logical. If \code{TRUE}, also return survival-at-length
#'   and associated mortality objects.
#'
#' @return If \code{return_surv = FALSE}, an \code{FLQuant} of equilibrium
#'   numbers-at-length per recruit.
#'
#' If \code{return_surv = TRUE}, a list with elements:
#' \describe{
#'   \item{\code{NI}}{An \code{FLQuant} of numbers-at-length per recruit.}
#'   \item{\code{surv}}{An \code{FLQuant} of survival-at-length.}
#'   \item{\code{Z_l}}{An \code{FLQuant} of total mortality-at-length.}
#'   \item{\code{F_l}}{An \code{FLQuant} of fishing mortality-at-length.}
#'   \item{\code{sel}}{An \code{FLQuant} of apical-scaled selectivity-at-length.}
#' }
#'
#' @details
#' The mortality-at-length objects are first constructed using
#' \code{calc_Z_l()}. Total mortality is then converted to the internal FLicc
#' recursion scale by dividing annual mortality by the fitted growth parameter
#' \code{k}. The recursion itself is evaluated by \code{n_f_flicc()}.
#'
#' This function is intended as a low-level engine for:
#' \itemize{
#'   \item \code{nf_flicc()} for scaled equilibrium numbers,
#'   \item \code{cf_flicc()} for catch-at-length,
#'   \item \code{pr_flicc()} for YPR and SBPR calculations.
#' }
#'
#' @seealso
#' \code{\link{calc_Z_l}},
#' \code{\link{flicc_refpars}},
#' \code{\link{n_f_flicc}},
#' \code{\link{nf_flicc}}
#'
#' @examples
#' \dontrun{
#' # raw numbers-per-recruit at FM = 1
#' nf0 <- nf_from_flicc(fit, FM = 1)
#'
#' # also return survival and mortality objects
#' nf1 <- nf_from_flicc(fit, FM = 1, return_surv = TRUE)
#' nf1$NI
#' nf1$surv
#' nf1$Z_l
#' }
#'
#' @export
nf_from_flicc <- function(fit,
                          F = NULL, FM = NULL,
                          nyears = 1,
                          Sel = NULL,
                          scale_sel = TRUE,
                          return_surv = FALSE) {

  zobj <- calc_Z_l(
    fit = fit,
    Sel = Sel,
    F = F,
    FM = FM,
    nyears = nyears
  )

  ref <- flicc_refpars(
    fit = fit,
    nyears = nyears,
    scale_sel = scale_sel
  )

  pop_model <- pop_model_flicc(fit)

  if (pop_model == "gamma") {

    res <- n_f_flicc(
      node = ref$node,
      quad_wt = ref$quad_wt,
      Len = ref$Len,
      Zki = as.numeric(zobj$Z_l) / ref$k,
      Galpha = ref$Galpha,
      Gbeta = ref$Gbeta,
      return_surv = return_surv
    )

  } else if (pop_model == "gtg") {

    td <- fit$tmb_data
    gtg <- gtg_inputs_flicc(fit)

    nlen <- length(ref$Len)
    ngtg <- length(gtg$gtgLinfs)

    # total mortality on the gtg scale
    # natural mortality comes from MKMat * yearly Mscaler
    # fishing mortality enters through zobj$F_l / k
    Mscale <- as.numeric(apply(td$Mscaler[, tail(seq_len(ncol(td$Mscaler)), nyears), drop = FALSE], 1, mean))

    Fki <- as.numeric(zobj$F_l) / ref$k

    ZKLMat <- matrix(NA_real_, nrow = nlen + 1, ncol = ngtg)

    for (g in seq_len(ngtg)) {
      for (l in seq_len(nlen)) {
        ZKLMat[l, g] <- gtg$MKMat[l, g] * Mscale[l] + Fki[l]
      }
      ZKLMat[nlen + 1, g] <- gtg$MKMat[nlen + 1, g] * Mscale[nlen] + Fki[nlen]
    }

    res <- n_f_gtg_flicc(
      Len = ref$Len,
      gtgLinfs = gtg$gtgLinfs,
      ZKLMat = ZKLMat,
      recP = gtg$recP,
      return_surv = return_surv
    )

  } else {
    stop("Unsupported pop_model: ", pop_model)
  }
  NI <- zobj$Z_l
  FLCore::units(NI) <- "1000"

  if (!return_surv) {
    NI[] <- res / 1000
    return(NI)
  }

  NI[] <- res$NI / 1000

  surv <- zobj$Z_l
  FLCore::units(surv) <- ""
  surv[] <- res$surv

  list(
    NI = NI,
    surv = surv,
    Z_l = zobj$Z_l,
    F_l = zobj$F_l,
    sel = zobj$sel
  )
}



#' Validate fishing input for FLicc helpers
#'
#' @param F Optional numeric scalar fully-selected fishing mortality.
#' @param FM Optional numeric scalar F/M multiplier.
#' @param default_FM Optional default value for \code{FM} if both \code{F}
#'   and \code{FM} are \code{NULL}.
#'
#' @return A list with elements \code{F} and \code{FM}.
#' @export
resolve_f_input_flicc <- function(F = NULL, FM = NULL, default_FM = NULL) {

  if (is.null(F) && is.null(FM)) {
    if (is.null(default_FM)) {
      stop("Provide one of 'F' or 'FM'.")
    }
    FM <- default_FM
  }

  if (!is.null(F) && !is.null(FM)) {
    stop("Provide only one of 'F' or 'FM', not both.")
  }

  list(F = F, FM = FM)
}

#' Calculate fishing and total mortality at length from an FLicc fit
#'
#' Computes annual natural mortality, fishing mortality, and total mortality
#' at length from a fitted \code{"flicc_tmb_fit"} object, averaged over the
#' selected terminal years and returned as \code{FLQuants}.
#'
#' Natural mortality is converted from the internal \code{Mk} scale to annual
#' mortality using the fitted growth parameter \code{k}. Fishing mortality is
#' then constructed either from an apical fishing mortality \code{F} or an
#' \code{F/M} multiplier \code{FM}, applied to the apical-scaled selectivity
#' curve.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param Sel Optional \code{FLQuant} or numeric vector of selectivity at
#'   length. If \code{NULL}, the terminal-year average of
#'   \code{fit$report$sel} is used.
#' @param F Optional numeric scalar apical fishing mortality.
#' @param FM Optional numeric scalar apical \code{F/M} multiplier.
#' @param M Optional \code{FLQuant} of annual natural mortality at length.
#'   If \code{NULL}, it is derived from \code{fit$report$Mk} and the fitted
#'   growth parameter \code{k}.
#' @param nyears Integer. Number of terminal years to average over.
#'
#' @return An \code{FLQuants} object with components:
#' \describe{
#'   \item{\code{F_l}}{Fishing mortality at length.}
#'   \item{\code{Z_l}}{Total mortality at length.}
#'   \item{\code{sel}}{Apical-scaled selectivity at length.}
#' }
#'
#' @export

calc_Z_l <- function(fit, Sel = NULL, F = NULL, FM = NULL, M = NULL, nyears=1) {
  if(is.null(F) & is.null(FM)) FM= 1
  years = an(dimnames(fit$report$N)$year)
  lens = an(an(dimnames(fit$report$Mk)$len))
  nyears <- min(length(years),nyears)
  years <- tail(years,nyears)
  if (is.null(M)) {
    M <- yearMeans(fit$report$M[, ac(years)])
  }
  units(M) <- "m"
  if(is.null(Sel)) Sel <- yearMeans(fit$report$sel[,ac(years)])

  nlen <- length(lens)
  fin <- resolve_f_input_flicc(F = F, FM = FM)

  sel_apical <- Sel/max(Sel,na.rm=T)

  smax <- max(Sel, na.rm = TRUE)
  if (!is.finite(smax) || smax <= 0) {
    stop("Combined selectivity must have a positive maximum.")
  }

  if (is.null(fin$F) && is.null(M)) {
    stop("Provide scalar 'M' when using 'FM'.")
  }

  F_l <- if (!is.null(fin$F)) {
    fin$F * sel_apical
  } else {
    fin$FM * M * sel_apical
  }

  Z_l <- M + F_l

  FLQuants(
    F_l = F_l,
    Z_l = Z_l,
    sel = sel_apical
  )
}



#' Extract reference biological and numerical quantities from an FLicc fit
#'
#' Extracts terminal-year or recent-year average biological parameters,
#' selectivity, weight-at-length, and numerical integration settings from a
#' fitted \code{"flicc_tmb_fit"} object.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param nyears Integer. Number of terminal years to average over.
#' @param scale_sel Logical. If \code{TRUE}, selectivity is scaled to an
#'   apical maximum of 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{year}}{Selected year labels.}
#'   \item{\code{Len}}{Length midpoints.}
#'   \item{\code{Mk}}{Average natural mortality on the internal Mk scale.}
#'   \item{\code{M}}{Scalar annual natural mortality.}
#'   \item{\code{k}}{Growth parameter used for annual scaling.}
#'   \item{\code{Sel}}{Average selectivity-at-length.}
#'   \item{\code{wl}}{Average weight-at-length.}
#'   \item{\code{Galpha, Gbeta}}{Growth parameters for the FLicc recursion.}
#'   \item{\code{node, quad_wt}}{Gauss-Laguerre quadrature settings.}
#' }
#'
#' @export

flicc_refpars <- function(fit, nyears = 1, scale_sel = TRUE) {

  yrs <- dimnames(fit$report$Mk)$year
  nyears <- min(length(yrs),nyears)
  year = ac(tail(yrs ,nyears))

  Len <- as.numeric(dimnames(fit$report$N)$len)
  Mk  <- as.numeric(yearMeans(fit$report$Mk[, year]))
  Sel <- as.numeric(yearMeans(fit$report$sel[, year]))
  wl  <- as.numeric(yearMeans(fit$report$w_at_len[, year]))

  if (scale_sel) {
    smax <- max(Sel, na.rm = TRUE)
    if (is.finite(smax) && smax > 0) {
      Sel <- Sel / smax
    }
  }

  out <- list(
    year = year,
    Len = Len,
    Mk = Mk,
    M = as.numeric(fit$report$lhpar["M"]),
    k = as.numeric(fit$report$lhpar["k"]),
    Sel = Sel,
    wl = wl,
    nyears= 1
  )
  if(fit$pop_model=="gamma"){
    out$Galpha = as.numeric(fit$report$pars["Galpha"])
    out$Gbeta  = as.numeric(fit$report$pars["Gbeta"])
    out$node = as.numeric(fit$report$node)
    out$quad_wt = as.numeric(fit$report$quad_wt)
  }
  return(out)
}

#' Fast per-recruit quantities from an FLicc fit, computed in TMB
#'
#' Evaluates numbers, catch and fishing mortality at length, yield and
#' spawning biomass per recruit for a vector of fishing levels, and the F
#' giving target SPR values, in one compiled TMB evaluation at the fitted
#' parameters (double precision, no AD taping). It reproduces
#' \code{pr_flicc()} / \code{spr_flicc()} for \code{Sel = NULL} and
#' \code{scale_sel = TRUE}, and is used internally by \code{fspr_flicc()},
#' \code{spr_flicc()}, \code{pr_flicc()}, \code{nf_flicc()} and
#' \code{prbrp_flicc()} (so also \code{eqstklen()}, \code{LBIspr()} and
#' \code{flicc2FLStockR()}). Set \code{options(FLicc.brp_tmb = FALSE)} to use
#' the previous R implementation instead.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object (FLicc >= 1.0.9).
#' @param Fseq Numeric vector of apical F (or F/M if \code{FM = TRUE}).
#' @param spr Numeric vector of target SPR in percent for Fspr.
#' @param nyears Number of terminal years averaged.
#' @param FM Logical; \code{Fseq} and the returned Fspr are F/M multipliers
#'   of natural mortality at length.
#' @param spawn_time Spawning time fraction.
#' @param Fmax Upper bound of the Fspr search.
#'
#' @return A list with \code{F}, \code{SBPR0}, \code{SBPR}, \code{YPR},
#'   matrices \code{N}, \code{C}, \code{Fl} (length x F) and \code{Fspr}
#'   (\code{NA} if not reached), or \code{NULL} if the fit cannot be
#'   evaluated this way (e.g. fits from earlier FLicc versions).
#' @export
brp_tmb_flicc <- function(fit, Fseq = numeric(0), spr = numeric(0), nyears = 1,
                          FM = FALSE, spawn_time = 0, Fmax = 5) {
  td <- fit$tmb_data
  if (is.null(fit$obj) || is.null(fit$opt$par) || is.null(td$gtg_z)) return(NULL)
  k <- as.numeric(fit$stklen@lhpar["k"])
  td$brp_F       <- as.numeric(Fseq)
  td$brp_spr     <- as.numeric(spr) / 100
  td$brp_nyears  <- as.integer(nyears)
  td$brp_k       <- k
  td$brp_FMscale <- as.integer(isTRUE(FM))
  td$brp_spawn   <- as.numeric(spawn_time)
  td$brp_Fmax    <- as.numeric(Fmax)
  td <- td[setdiff(names(td), c("gear_names", "year_names", "Linf_init",
                                "Mk_init", "CVL_init", "Galpha_init"))]
  pl <- fit$obj$env$parList(fit$opt$par)
  # all parameters fixed at the estimates: nothing is recorded on the AD
  # tape, and report() runs one double evaluation including the brp block
  map <- lapply(pl, function(x) factor(rep(NA, length(x))))
  obj <- TMB::MakeADFun(data = td, parameters = pl, map = map,
                        DLL = fit$obj$env$DLL, silent = TRUE)
  rep <- obj$report()
  Fspr <- as.numeric(rep$brp_Fspr)
  Fspr[Fspr < 0] <- NA_real_
  list(F = as.numeric(Fseq), SBPR0 = rep$brp_SBPR0,
       SBPR = as.numeric(rep$brp_SBPR), YPR = as.numeric(rep$brp_YPR),
       N = rep$brp_N, C = rep$brp_C, Fl = rep$brp_Fl, Fspr = Fspr)
}

brp_fast_ok <- function(fit, Sel = NULL, scale_sel = TRUE) {
  is.null(Sel) && isTRUE(scale_sel) && !isFALSE(getOption("FLicc.brp_tmb", TRUE)) &&
    !is.null(fit$obj) && !is.null(fit$tmb_data$gtg_z) && !is.null(fit$report$M)
}

# FLQuant template with the same dimensions as calc_Z_l() output
brp_tmpl_flicc <- function(fit, nyears = 1) {
  years <- dimnames(fit$report$N)$year
  years <- tail(years, min(length(years), nyears))
  FLCore::yearMeans(fit$report$M[, years])
}

#' Equilibrium numbers-at-length from an FLicc fit
#'
#' Computes equilibrium numbers-at-length from a fitted
#' \code{"flicc_tmb_fit"} object under a specified fishing mortality.
#'
#' Fishing can be supplied either as apical annual fishing mortality
#' (\code{F}) or as an \code{F/M} multiplier (\code{FM}). Numbers are scaled
#' so that recruitment at the smallest length class equals \code{R0} under
#' unfished conditions.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param nyears Integer. Number of terminal years to average over.
#' @param F Optional numeric scalar apical fishing mortality.
#' @param FM Optional numeric scalar apical \code{F/M} multiplier.
#' @param Sel Optional selectivity override.
#' @param scale_sel Logical. If \code{TRUE}, selectivity is scaled to an
#'   apical maximum of 1.
#' @param return_surv Logical. If \code{TRUE}, also return survival and
#'   intermediate mortality objects.
#' @param R0 Numeric recruitment scaling constant.
#'
#' @return An \code{FLQuant} of numbers-at-length, or a list if
#'   \code{return_surv = TRUE}.
#'
#' @export

nf_flicc <- function(fit, nyears = 1, F = NULL, FM = NULL,
                     Sel = NULL, scale_sel = TRUE,
                     return_surv = FALSE,
                     R0 = 1000) {

  fin <- resolve_f_input_flicc(F = F, FM = FM, default_FM = 1)

  if (!return_surv && brp_fast_ok(fit, Sel, scale_sel)) {
    b <- brp_tmb_flicc(fit, Fseq = if (is.null(fin$F)) fin$FM else fin$F,
                       nyears = nyears, FM = is.null(fin$F))
    if (!is.null(b)) {
      out <- brp_tmpl_flicc(fit, nyears)
      out[] <- b$N[, 1] * R0 / 1000
      FLCore::units(out) <- "1000"
      return(out)
    }
  }

  ref <- flicc_refpars(fit, nyears = nyears, scale_sel = scale_sel)

  if (!is.null(Sel)) {
    ref$Sel <- as.numeric(Sel)
  }

  res0 <- nf_from_flicc(
    fit,
    F = 0,
    FM = NULL,
    nyears = nyears,
    return_surv = FALSE
  )

  scale0 <- an(res0[1, 1])

  res <- nf_from_flicc(
    fit,
    F = fin$F,
    FM = fin$FM,
    nyears = nyears,
    return_surv = return_surv
  )

  if (!return_surv) {
    return(res / an(scale0) * R0 / 1000)
  }

  res$NI <- res$NI / an(scale0) * R0
  res$surv <- res$surv / scale0 * R0
  res
}

#' Equilibrium catch-at-length from an FLicc fit
#'
#' Computes equilibrium catch-at-length from a fitted
#' \code{"flicc_tmb_fit"} object under a specified fishing mortality,
#' using the Baranov catch equation.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param nyears Integer. Number of terminal years to average over.
#' @param F Optional numeric scalar apical fishing mortality.
#' @param FM Optional numeric scalar apical \code{F/M} multiplier.
#' @param Sel Optional selectivity override.
#' @param scale_sel Logical. If \code{TRUE}, selectivity is scaled to an
#'   apical maximum of 1.
#' @param biomass Logical. If \code{TRUE}, return catch biomass-at-length
#'   instead of catch numbers-at-length.
#' @param return_all Logical. If \code{TRUE}, return a list of intermediate
#'   quantities.
#' @param R0 Numeric recruitment scaling constant.
#'
#' @return An \code{FLQuant} of catch-at-length, or a list when
#'   \code{return_all = TRUE}.
#'
#' @export

cf_flicc <- function(fit, nyears = 1, F = NULL, FM = NULL,
                     Sel = NULL, scale_sel = TRUE,
                     biomass = FALSE, return_all = FALSE,
                     R0 = 1000) {

  year <- dimnames(fit$report$N)$year[]
  yrs <- tail(year,min(nyears,length(year)))


  ref <- flicc_refpars(fit, scale_sel = scale_sel,nyears=nyears)

  if (!is.null(F) && is.null(FM)) {
    FM <- F / ref$M
    F <- NULL
  }

  fin <- resolve_f_input_flicc(F = F, FM = FM, default_FM = 1)

  if (!is.null(Sel)) {
    ref$Sel <- as.numeric(Sel)
  }

  zobj <- calc_Z_l(
    fit,
    F = fin$F,
    FM = fin$FM,
  )

  Nl <- nf_flicc(
    fit = fit,
    nyears = nyears,
    F = fin$F,
    FM = fin$FM,
    scale_sel = FALSE,
    return_surv = FALSE,
    R0 = R0
  )

  Cl <- Nl * zobj$F_l / zobj$Z_l * (1 - exp(-zobj$Z_l))
  Cl[!is.finite(Cl)] <- 0
  units(Cl) <- "1000"

  wl <- fit$report$w_at_len[,yrs]
  out <- if(biomass) Cl * wl else Cl

  if (!return_all) {
    return(out)
  }

  list(
    C =  Cl * wl,
    Cn = Cl,
    N = Nl,
    Z = zobj$Z_l,
    F = zobj$F_l,
    Sel = zobj$sel,
    year = an(ref$year),
    Len = an(ref$Len),
    Finput = fin$F,
    FMinput = fin$FM,
    R0 = R0
  )
}

#' Per-recruit quantities from an FLicc fit
#'
#' Computes yield-per-recruit and spawning biomass-per-recruit from a fitted
#' \code{"flicc_tmb_fit"} object under a specified fishing mortality.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param nyears Integer. Number of terminal years to average over.
#' @param F Optional numeric scalar apical fishing mortality.
#' @param FM Optional numeric scalar apical \code{F/M} multiplier.
#' @param Sel Optional selectivity override.
#' @param scale_sel Logical. If \code{TRUE}, selectivity is scaled to an
#'   apical maximum of 1.
#' @param spawn_time Numeric fraction of the year at which spawning occurs.
#' @param R0 Numeric recruitment scaling constant.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{YPR}}{Yield per recruit.}
#'   \item{\code{SBPR}}{Spawning biomass per recruit.}
#'   \item{\code{Cn}}{Catch-at-length per recruit.}
#'   \item{\code{N}}{Numbers-at-length per recruit.}
#'   \item{\code{Z}}{Total mortality-at-length.}
#'   \item{\code{F}}{Fishing mortality-at-length.}
#'   \item{\code{Len}}{Length midpoints.}
#'   \item{\code{year}}{Reference year label.}
#'   \item{\code{R0}}{Recruitment scaling used.}
#' }
#'
#' @export

pr_flicc <- function(fit, nyears = 1, F = NULL, FM = NULL,
                     Sel = NULL, scale_sel = TRUE,
                     spawn_time = 0,
                     R0 = 1000) {

  fin <- resolve_f_input_flicc(F = F, FM = FM, default_FM = 1)

  if (brp_fast_ok(fit, Sel, scale_sel)) {
    b <- brp_tmb_flicc(fit, Fseq = if (is.null(fin$F)) fin$FM else fin$F,
                       nyears = nyears, FM = is.null(fin$F),
                       spawn_time = spawn_time)
    if (!is.null(b)) {
      sc <- R0 / 1000
      Mq <- brp_tmpl_flicc(fit, nyears)
      Fq <- Mq; Fq[] <- b$Fl[, 1]
      Zq <- Mq + Fq
      Nl <- Mq; Nl[] <- b$N[, 1] * sc; FLCore::units(Nl) <- "1000"
      Cn <- Mq; Cn[] <- b$C[, 1] * sc; FLCore::units(Cn) <- "1000"
      ref <- flicc_refpars(fit, nyears = nyears, scale_sel = scale_sel)
      return(list(YPR = b$YPR * sc, SBPR = b$SBPR * sc, Cn = Cn, N = Nl,
                  Z = Zq, F = Fq, Len = ref$Len, year = ref$year, R0 = R0))
    }
  }

  ref <- flicc_refpars(fit, nyears = nyears, scale_sel = scale_sel)

  if (!is.null(Sel)) {
    ref$Sel <- as.numeric(Sel)
  }

  zobj <- calc_Z_l(
    fit,
    F = fin$F,
    FM = fin$FM,
    nyears = nyears
  )

  Nl <- nf_flicc(
    fit = fit,
    nyears = nyears,
    F = fin$F,
    FM = fin$FM,
    Sel = ref$Sel,
    scale_sel = FALSE,
    return_surv = FALSE,
    R0 = R0
  )

  year_chr <- ref$year
  matl <- (FLCore::mat(fit$stklen)[, year_chr])

  Cn <- Nl * zobj$F_l / zobj$Z_l * (1 - exp(-zobj$Z_l))
  Cn[!is.finite(Cn)] <- 0
  units(Cn) <- "1000"
  wl <- yearMeans(fit$report$w_at_len)
  wl[] <- ref$w
  units(wl) <- "kg"

  YPR <- sum(Cn * wl, na.rm = TRUE)
  SSN <- Nl * matl * exp(-zobj$Z_l * spawn_time)
  SBPR <- sum(SSN * wl, na.rm = TRUE)

  list(
    YPR = YPR,
    SBPR = SBPR,
    Cn = Cn,
    N = Nl,
    Z = zobj$Z_l,
    F = zobj$F_l,
    Len = ref$Len,
    year = year_chr,
    R0 = R0
  )
}

#' Spawning potential ratio from an FLicc fit
#'
#' Computes spawning potential ratio (SPR) as the ratio of spawning biomass
#' per recruit under a specified fishing mortality to unfished spawning
#' biomass per recruit.
#'
#' @inheritParams pr_flicc
#'
#' @return Numeric or \code{FLQuant}-like SPR value, depending on upstream
#'   object classes.
#'
#' @export

spr_flicc <- function(fit, nyears = 1, F = NULL, FM = NULL,
                      Sel = NULL, scale_sel = TRUE,
                      spawn_time = 0) {

  fin <- resolve_f_input_flicc(F = F, FM = FM)

  if (brp_fast_ok(fit, Sel, scale_sel)) {
    b <- brp_tmb_flicc(fit, Fseq = if (is.null(fin$F)) fin$FM else fin$F,
                       nyears = nyears, FM = is.null(fin$F),
                       spawn_time = spawn_time)
    if (!is.null(b)) {
      if (!is.finite(b$SBPR0) || b$SBPR0 <= 0) return(NA_real_)
      return(b$SBPR / b$SBPR0)
    }
  }

  sb0 <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = 0,
    FM = NULL,
    Sel = Sel,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )$SBPR

  sbf <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = fin$F,
    FM = fin$FM,
    Sel = Sel,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )$SBPR

  if (!is.finite(sb0) || sb0 <= 0) {
    return(NA_real_)
  }

  sbf / sb0
}

#' Solve for fishing mortality giving a target SPR
#'
#' Solves for the fishing mortality that gives a requested target spawning
#' potential ratio.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param spr Numeric SPR target in percent, for example \code{40} for
#'   \code{SPR = 0.4}.
#' @param nyears Integer. Number of terminal years to average over.
#' @param Sel Optional selectivity override.
#' @param scale_sel Logical. If \code{TRUE}, selectivity is scaled to an
#'   apical maximum of 1.
#' @param spawn_time Numeric spawning timing fraction.
#' @param interval Numeric interval passed to \code{stats::uniroot()}.
#' @param input Character string indicating whether the solution should be
#'   returned on the \code{"FM"} or \code{"F"} scale.
#'
#' @return Numeric fishing mortality corresponding to the requested SPR.
#'
#' @export

fspr_flicc <- function(fit, spr = 40, nyears = 1,
                       Sel = NULL, scale_sel = TRUE,
                       spawn_time = 0,
                       interval = c(1e-8, 5),
                       input = c( "F","FM")) {

  input <- match.arg(input)

  if (brp_fast_ok(fit, Sel, scale_sel)) {
    b <- brp_tmb_flicc(fit, spr = spr, nyears = nyears, FM = input == "FM",
                       spawn_time = spawn_time, Fmax = interval[2])
    if (!is.null(b)) {
      if (anyNA(b$Fspr)) {
        stop("Target SPR not reached for F in (1e-8, ", interval[2], ").")
      }
      return(b$Fspr)
    }
  }

  target <- spr / 100

  fobj <- function(x) {
    if (input == "FM") {
      spr_flicc(
        fit = fit,
        nyears = nyears,
        FM = x,
        Sel = Sel,
        scale_sel = scale_sel,
        spawn_time = spawn_time
      ) - target
    } else {
      spr_flicc(
        fit = fit,
        nyears = nyears,
        F = x,
        Sel = Sel,
        scale_sel = scale_sel,
        spawn_time = spawn_time
      ) - target
    }
  }

  stats::uniroot(fobj, interval = interval)$root
}

#' Biomass-equivalent reference points from SPR under steepness
#'
#' Under Beverton-Holt recruitment with steepness \code{s}, equilibrium
#' spawning biomass is linear in SPR, \code{B/B0 = (4 s SPR - (1 - s)) / (5 s - 1)},
#' so biomass falls faster than SPR as the stock is depleted. These helpers
#' convert between SPR and equilibrium biomass relative to the biomass at a
#' reference SPR (e.g. SPR40 as an Bmsy proxy, or SPRmsy), without refitting.
#'
#' \code{fbrel_flicc()} returns the fishing mortality at which equilibrium
#' biomass is \code{brel} times the reference biomass, e.g. \code{brel = 0.25}
#' for a Blim at a quarter of the B40 proxy. It solves for the equivalent SPR
#' in closed form and then calls \code{fspr_flicc()}, so it is as fast as
#' \code{fspr_flicc()}.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param brel Equilibrium biomass relative to the reference biomass
#'   (e.g. \code{0.25}, \code{0.5}, \code{1}). Vectorised in the helpers.
#' @param s Beverton-Holt steepness.
#' @param sprref Reference SPR in percent, defining the reference biomass
#'   (default \code{40}). Use \code{100 * SPRmsy} to express results relative
#'   to Bmsy.
#' @param spr SPR as a proportion (e.g. \code{fit$report$spr} or
#'   \code{sprcur_flicc(fit)}); numeric or \code{FLQuant}.
#' @param ... Further arguments passed to \code{fspr_flicc()}, e.g.
#'   \code{nyears}, \code{input}.
#'
#' @return \code{fbrel_flicc()}: fishing mortality (scale set by \code{input}).
#'   \code{sprbrel_flicc()}: SPR in percent giving \code{brel}.
#'   \code{spr2brel()}: equilibrium biomass relative to the reference
#'   biomass, truncated at zero (stock collapse).
#'
#' @examples
#' \dontrun{
#' sprbrel_flicc(c(1, 0.5, 0.25), s = 0.75)   # 40.0 24.2 16.2 (percent)
#' fbrel_flicc(fit, brel = 0.25, s = 0.75)    # F at 0.25 * B40
#' spr2brel(fit$report$spr, s = 0.75)         # status as B / B40
#' }
#' @export
fbrel_flicc <- function(fit, brel = 0.25, s = 0.75, sprref = 40, ...) {
  if (length(brel) != 1) stop("fbrel_flicc() takes a single 'brel'.")
  fspr_flicc(fit, spr = sprbrel_flicc(brel, s = s, sprref = sprref), ...)
}

#' @rdname fbrel_flicc
#' @export
sprbrel_flicc <- function(brel, s = 0.75, sprref = 40) {
  check_steep_flicc(s, sprref)
  ref <- sprref / 100
  bref <- 4 * s * ref - (1 - s)
  100 * ((1 - s) + brel * bref) / (4 * s)
}

#' @rdname fbrel_flicc
#' @export
spr2brel <- function(spr, s = 0.75, sprref = 40) {
  check_steep_flicc(s, sprref)
  ref <- sprref / 100
  bref <- 4 * s * ref - (1 - s)
  b <- 4 * s * spr - (1 - s)
  b[] <- pmax(as.numeric(b), 0)
  b / bref
}

check_steep_flicc <- function(s, sprref) {
  if (length(s) != 1 || !is.finite(s) || s <= 0.2 || s > 1) {
    stop("'s' must be a single steepness value in (0.2, 1].")
  }
  if (4 * s * sprref / 100 <= 1 - s) {
    stop("The reference SPR (", sprref, "%) is at or below the collapse SPR ",
         "for s = ", s, " (", round(100 * (1 - s) / (4 * s), 1), "%).")
  }
  invisible(TRUE)
}
#' Yield-per-recruit from an FLicc fit
#'
#' @inheritParams pr_flicc
#' @return Numeric YPR.
#' @export

ypr_flicc <- function(fit, nyears = 1, F = NULL, FM = NULL,
                      Sel = NULL, scale_sel = TRUE,
                      spawn_time = 0) {

  pr_flicc(
    fit = fit,
    nyears = nyears,
    F = F,
    FM = FM,
    Sel = Sel,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )$YPR
}


