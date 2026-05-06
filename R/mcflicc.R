#' Propagate FLQuant objects in a FLicc report
#'f
#' @param report A `fit$report` list from `fiticc()`.
#' @param nsim Number of Monte Carlo iterations.
#'
#' @return A report list where all `FLQuant` objects have an expanded `iter` dimension.
prop_report <- function(report, nsim) {

  out <- report

  is_flq <- vapply(out, inherits, logical(1), "FLQuant")

  for (nm in names(out)[is_flq]) {
    out[[nm]] <- FLCore::propagate(out[[nm]], nsim)
  }

  out
}


#' Fill one iteration of a propagated FLicc report
#'
#' @param mcrep A propagated report from `prop_report()`.
#' @param report_i A report list from one Monte Carlo fit.
#' @param i Iteration number.
#'
#' @return Updated Monte Carlo report.
fill_report_iter <- function(mcrep, report_i, i) {

  for (nm in names(mcrep)) {

    if (
      inherits(mcrep[[nm]], "FLQuant") &&
      !is.null(report_i[[nm]]) &&
      inherits(report_i[[nm]], "FLQuant")
    ) {
      mcrep[[nm]][, , , , , i] <- report_i[[nm]]
    }
  }

  mcrep
}

#' Default FishLife-style life-history correlation matrix
#'
#' Returns a default correlation matrix for drawing correlated life-history
#' parameters in \code{mc_flicc()}.
#'
#' @return A numeric correlation matrix.
#'
#' @export
fishlife_corr <- function() {

  corr <- matrix(
    c(
      1.0000000, -0.7648994, -0.6878538,  0.9643831,
      -0.7648994,  1.0000000,  0.8588491, -0.6939758,
      -0.6878538,  0.8588491,  1.0000000, -0.6859586,
      0.9643831, -0.6939758, -0.6859586,  1.0000000
    ),
    nrow = 4,
    byrow = TRUE
  )

  dimnames(corr) <- list(
    c("linf", "k", "M", "L50"),
    c("linf", "k", "M", "L50")
  )

  corr
}

#' Extract life-history parameters for Monte Carlo sampling
#'
#' Extracts the life-history parameters required for FLicc Monte Carlo
#' uncertainty propagation from a fitted FLicc object.
#'
#' @param fit A fitted FLicc object returned by \code{fiticc()}.
#'
#' @return A named numeric vector or \code{FLPar} containing life-history
#'   parameters used by \code{draw_lh_mc()}.
#'
#' @export
get_lhpar_mc <- function(fit) {

  lh <- fit$stklen@lhpar

  list(
    linf = as.numeric(lh["linf"]),
    k    = as.numeric(lh["k"]),
    M    = as.numeric(lh["M"]),
    Mk   = as.numeric(lh["Mk"]),
    L50  = as.numeric(lh["L50"]),
    L95  = as.numeric(lh["L95"]),
    CVL  = fit$settings$CVL
  )
}




#' Draw Monte Carlo life-history parameters
#'
#' Generates Monte Carlo draws of uncertain life-history parameters for use in
#' \code{mc_flicc()}. Parameters are drawn around the fitted values using
#' user-specified coefficients of variation and an optional correlation matrix.
#'
#' @param fit A fitted FLicc object returned by \code{fiticc()}.
#' @param nsim Number of Monte Carlo draws.
#' @param cv Named list of coefficients of variation for parameters to be
#'   varied. Common entries are \code{linf}, \code{k}, \code{M}, \code{L50},
#'   and \code{CVL}.
#' @param corr Correlation matrix for correlated draws. Defaults are usually
#'   supplied by \code{fishlife_corr()}.
#'
#' @return A data frame of Monte Carlo life-history parameter draws, with one
#'   row per simulation.
#'
#' @export
draw_lh_mc <- function(
    fit,
    nsim,
    cv = list(
      linf = 0.10,
      k    = 0.15,
      M    = 0.25,
      L50  = 0.10,
      CVL  = 0.20
    ),
    corr = fishlife_corr(),
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  lh0 <- get_lhpar_mc(fit)

  cv <- unlist(cv)

  req <- c("linf", "k", "M", "L50", "CVL")
  if (!all(req %in% names(cv))) {
    stop("`cv` must contain: ", paste(req, collapse = ", "))
  }

  mu <- log(c(
    linf = lh0$linf,
    k    = lh0$k,
    M    = lh0$M,
    L50  = lh0$L50
  ))

  logsd <- sqrt(log(1 + cv[c("linf", "k", "M", "L50")]^2))

  Sigma <- (logsd %*% t(logsd)) * corr

  z <- MASS::mvrnorm(
    n = nsim,
    mu = mu,
    Sigma = Sigma
  )

  draws <- as.data.frame(exp(z))
  names(draws) <- c("linf", "k", "M", "L50")

  draws$Mk <- draws$M / draws$k

  L95_ratio <- lh0$L95 / lh0$L50
  draws$L95 <- draws$L50 * L95_ratio

  cvl_logsd <- sqrt(log(1 + cv["CVL"]^2))
  cvl_logmu <- log(lh0$CVL) - 0.5 * cvl_logsd^2

  draws$CVL <- stats::rlnorm(
    nsim,
    meanlog = cvl_logmu,
    sdlog = cvl_logsd
  )

  draws
}

#' Monte Carlo refits for FLicc uncertainty analysis
#'
#' Runs Monte Carlo refits of a fitted FLicc model to propagate uncertainty in
#' externally specified life-history and model inputs. For each simulation,
#' uncertain parameters such as \eqn{L_\infty}, \eqn{k}, \eqn{M}, maturity
#' length and length-at-age variation are drawn from user-defined distributions.
#' The length-structured stock object is rebuilt and the FLicc model is refitted.
#'
#' This function is intended to support uncertainty analysis and ensemble
#' evaluation in length-based assessments. The resulting object can be used to
#' examine parameter sensitivity, likelihood support, structural uncertainty,
#' and derived quantities such as fishing mortality and spawning potential ratio.
#'
#' @param fit A fitted FLicc object returned by \code{fiticc()}.
#' @param nsim Number of Monte Carlo simulations.
#' @param cv Named list of coefficients of variation for uncertain inputs.
#'   Common entries are \code{linf}, \code{k}, \code{M}, \code{L50}, and
#'   \code{CVL}. Parameters not included are held fixed.
#' @param corr Correlation matrix used when drawing correlated life-history
#'   parameters. See \code{fishlife_corr()}.
#' @param seed Optional random seed for reproducibility.
#' @param verbose Logical. If \code{TRUE}, print progress and elapsed time.
#' @param parallel Logical. If \code{TRUE}, run Monte Carlo refits in parallel
#'   using \code{foreach}, \code{doFuture}, and \code{future}.
#' @param workers Optional number of parallel workers. If \code{NULL}, the
#'   current \code{future} plan is used.
#' @param n_restart Number of optimizer restarts passed to \code{fiticc()}.
#'   The default \code{0} is recommended for Monte Carlo refits because each
#'   run starts from the original optimum.
#' @param ... Additional arguments passed to internal fitting routines.
#'
#' @details
#' The Monte Carlo procedure is designed for uncertainty propagation rather than
#' formal hypothesis testing. In particular, it can be used to explore
#' structural uncertainty arising from alternative plausible life-history
#' inputs, growth assumptions, maturity schedules, or length-composition
#' variation. The original fitted model is stored in \code{mc$fit}, allowing
#' direct comparison between the reference fit and the Monte Carlo ensemble.
#'
#' Likelihood values stored in \code{mc$logLik} can be used to compute
#' \eqn{\Delta AIC} values when all refits use the same number of estimated
#' parameters. A practical ensemble filter is to retain simulations with
#' \eqn{\Delta AIC \le 20} as a broad plausible set for structural uncertainty,
#' while using \eqn{\Delta AIC \le 10} as a more strongly supported subset.
#'
#' @return An object of class \code{mcflicc}. The object contains propagated
#'   FLR quantities, Monte Carlo life-history draws, fitted parameter values,
#'   selectivity parameters, log-likelihoods, convergence flags, and the
#'   original fitted report stored in \code{mc$fit}.
#'
#' @examples
#' \dontrun{
#' # Basic Monte Carlo uncertainty analysis
#' mc <- mc_flicc(
#'   fit,
#'   nsim = 500,
#'   cv = list(linf = 0.10, k = 0.10, M = 0.15, L50 = 0.10, CVL = 0.10),
#'   seed = 123,
#'   verbose = TRUE
#' )
#'
#' # Run in parallel
#' # Check available cores/workers before running in parallel
#' parallel::detectCores()
#' parallelly::availableCores()
#'
#' # Start conservatively, especially for large FLR/TMB objects.
#' # For example, on a machine with many cores, 6--10 workers is often safer
#' # than using all available cores.

#' # Parallel refits
#' mc <- mc_flicc(
#'   fit,
#'   nsim = 500,
#'   parallel = TRUE,
#'   workers = 10,
#'   n_restart = 0
#' )
#'
#' # Likelihood-supported ensemble
#' nll <- -as.numeric(mc$logLik)
#' dAIC <- 2 * (nll - min(nll, na.rm = TRUE))
#' keep <- dAIC <= 20
#'
#' summary(an(mc$spr[, ac(2024)])[keep])
#'
#' # Diagnostic plots
#' plot_mcpars(mc)
#' plot_mccor(mc, dAIC_cut = 20)
#' }
#'
#' @export
mc_flicc <- function(
    fit,
    nsim = 100,
    cv = list(linf = 0.10, k = 0.1, M = 0.15, L50 = 0.10, CVL = 0.1),
    corr = fishlife_corr(),
    seed = 123,
    verbose = TRUE,
    parallel = FALSE,
    workers = NULL,
    n_restart = 0,
    ...
) {

  if (!is.null(seed)) set.seed(seed)

  # 1. Draw life-history parameters
  draws <- draw_lh_mc(fit, nsim = nsim, cv = cv, corr = corr)

  # 2. Initialise MC object from report
  mc <- fit$report
  choose <- c("Fap","spr","logLik", "lhpar",   "selpars", "pars")
  mc <- mc[names(mc)%in%choose]

  is_flq <- vapply(mc, inherits, logical(1), "FLQuant")
  is_flp <- vapply(mc, inherits, logical(1), "FLPar")
  #is_flqs <- vapply(mc, inherits, logical(1), "FLQuants")
  for (nm in names(mc)[is_flq]) {
    mc[[nm]] <- FLCore::propagate(mc[[nm]], nsim)
  }

  for (nm in names(mc)[is_flp]) {
    mc[[nm]] <- FLCore::propagate(mc[[nm]], nsim)
  }
  #for (nm in names(mc)[is_flqs]) {
  #  mc[[nm]] <- FLCore::propagate(mc[[nm]], nsim)
  #}
  mc$selpars <- FLCore::propagate(  mc$selpars, nsim)

  # 3. Propagate lhpar
  # add CVL
  lhpar <- fit$stklen@lhpar
  lhpar <- rbind(lhpar,FLPar(CVL=fit$settings$CVL))
  mc$lhpar <- FLCore::propagate( lhpar, nsim)

  # 4. Convergence flag
  mc$ok <- rep(FALSE, nsim)

  # 5. Paralell computing
  if (parallel) {
    requireNamespace("foreach")
    requireNamespace("doFuture")
    requireNamespace("future")
    requireNamespace("progressr")

    doFuture::registerDoFuture()

    if (!is.null(workers)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = workers)
    }
  }

  # 6. Monte Carlo loop function

  if (!parallel && verbose) {
    cat("Monte Carlo FLicc: ")
    pb <- utils::txtProgressBar(
      min = 0,
      max = nsim,
      style = 3,
      char = "><>"
    )
    on.exit(close(pb), add = TRUE)
  }

  # Run i loop
  run_one <- function(i) {

    stk_i <- fit$stklen
    lh <- lhpar

    lh["linf"] <- draws$linf[i]
    lh["k"]    <- draws$k[i]
    lh["M"]    <- draws$M[i]
    lh["Mk"]   <- draws$Mk[i]
    lh["L50"]  <- draws$L50[i]
    lh["L95"]  <- draws$L95[i]
    lh["CVL"]  <- draws$CVL[i]

    stk_i <- stocklen(
      fit$report$obslen,
      lh,
      m_model = fit$stklen@m_model
    )

    attr(stk_i, "lhpar") <- lh

    settings_i <- modifyList(
      fit$settings,
      list(CVL = draws$CVL[i])
    )

    sel_fun <- c("logistic", "dsnormal", "normal")[fit$tmb_data$sel_type]

    fit_i <- try(
      fiticc(
        fit$report$obslen,
        stk_i,
        sel_fun = sel_fun,
        catch_by_gear = fit$report$catch_by_gear,
        settings = settings_i,
        start = fit$opt$par,
        n_restart = n_restart
      ),
      silent = TRUE
    )

    if (inherits(fit_i, "try-error")) {
      return(list(i = i, ok = FALSE))
    }

    list(
      i = i,
      ok = TRUE,
      lhpar = lh,
      pars = fit_i$report$pars,
      selpars = fit_i$report$selpars,
      report = fit_i$report[names(fit_i$report)%in%choose]
    )
  }


  # Run iterations
  if (verbose) t0 <- Sys.time()

  if (parallel) {

    if (verbose) {
      message("Running Monte Carlo FLicc on ",
              future::nbrOfWorkers(), " workers")
    }

    res <- progressr::with_progress({
      p <- progressr::progressor(steps = nsim)

      foreach::foreach(
        i = seq_len(nsim),
        .options.future = list(seed = TRUE)
      ) %dofuture% {

        suppressWarnings(
          suppressPackageStartupMessages({

            out <- run_one(i)

            if (verbose) p()

            out
          })
        )
      }
    })

  } else {

    if (verbose) {
      cat("Monte Carlo FLicc: ")
      pb <- utils::txtProgressBar(
        min = 0,
        max = nsim,
        style = 3,
        char = "><>"
      )
      on.exit(close(pb), add = TRUE)
    }

    res <- vector("list", nsim)

    for (i in seq_len(nsim)) {
      res[[i]] <- run_one(i)
      if (verbose) utils::setTxtProgressBar(pb, i)
    }
  }


  # Fill results
  for (ri in res) {

    i <- ri$i

    if (!isTRUE(ri$ok)) next

    mc$ok[i] <- TRUE
    mc$logLik[i] <-  ri$report$logLik
    iter(mc$lhpar, i) <- ri$lhpar
    iter(mc$pars, i)  <- ri$pars


    for (g in seq_along(mc$selpars)) {
      iter(mc$selpars[[g]], i) <- ri$selpars[[g]]
    }

    for (nm in names(mc)[is_flq]) {
      if (!is.null(ri$report[[nm]]) &&
          inherits(ri$report[[nm]], "FLQuant")) {
        mc[[nm]][, , , , , i] <- ri$report[[nm]]
      }
    }

    #for (nm in names(mc)[is_flqs]) {
    #  if (!is.null(ri$report[[nm]]) &&
    ###      inherits(ri$report[[nm]], "FLQuants")) {
    #
    #    for (g in seq_along(mc[[nm]])) {
    #      mc[[nm]][[g]][, , , , , i] <- ri$report[[nm]][[g]]
    #    }
    #  }
    #}
  }
  # Timer
  if (verbose) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    msg <- if (elapsed < 60) {
      paste0(round(elapsed, 1), " sec")
    } else {
      paste0(floor(elapsed / 60), " min ", round(elapsed %% 60, 1), " sec")
    }

    cat("\nMonte Carlo FLicc completed in ", msg, "\n", sep = "")
  }
  mc$fit = fit$report[names(fit$report)%in%choose]
  mc$fit$lhpar <- rbind(mc$fit$lhpar,FLPar(CVL = fit$settings$CVL))

  class(mc) <- c("mcflicc", class(mc))
  mc
}

#' Compute delta AIC for FLicc Monte Carlo fits
#'
#' Computes relative AIC differences for Monte Carlo refits stored in an
#' object of class \code{mcflicc}. The calculation assumes that all Monte Carlo
#' refits have the same number of estimated parameters, so that
#' \code{delta AIC = 2 * delta NLL}.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#'
#' @return A numeric vector of delta AIC values, one per Monte Carlo iteration.
#'
#' @keywords internal
mc_dAIC <- function(mc) {
  nll <- -as.numeric(mc$logLik)
  2 * (nll - min(nll, na.rm = TRUE))
}



#' Compute Akaike weights for FLicc Monte Carlo fits
#'
#' Computes delta AIC values and Akaike weights for Monte Carlo refits stored in
#' an object returned by \code{mc_flicc()}. The calculation assumes that all
#' refits have the same number of estimated parameters, so that
#' \code{delta AIC = 2 * delta NLL}.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param dAIC_cut Optional delta AIC cutoff. If supplied, iterations with
#'   \code{dAIC > dAIC_cut} receive zero weight before weights are normalized.
#'
#' @return A data frame with iteration number, log-likelihood, delta AIC,
#'   normalized Akaike weight and retained/filtered flag.
#'
#' @export
mc_aic_weights <- function(mc, dAIC_cut = NULL) {

  dAIC <- mc_dAIC(mc)

  w <- exp(-0.5 * dAIC)

  if (!is.null(dAIC_cut)) {
    w[dAIC > dAIC_cut] <- 0
  }

  w <- w / sum(w, na.rm = TRUE)

  data.frame(
    iter = seq_along(dAIC),
    logLik = as.numeric(mc$logLik),
    dAIC = dAIC,
    weight = w,
    keep = if (is.null(dAIC_cut)) TRUE else dAIC <= dAIC_cut
  )

   data.frame(
    iter = seq_along(w),
    logLik = as.numeric(mc$logLik),
    dAIC = dAIC,
    weight = w,
    keep = if (is.null(dAIC_cut)) TRUE else dAIC <= dAIC_cut
  )
}

#' Weighted quantiles
#'
#' Computes weighted empirical quantiles for a numeric vector. This helper is
#' used internally by FLicc Monte Carlo summary functions to calculate
#' likelihood-weighted uncertainty intervals.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of non-negative weights with the same length as
#'   \code{x}.
#' @param probs Numeric vector of probabilities in \code{[0, 1]}. Default is
#'   \code{c(0.025, 0.5, 0.975)}.
#' @param na.rm Logical. If \code{TRUE}, remove non-finite values in \code{x}
#'   and \code{w}. Default is \code{TRUE}.
#'
#' @return A numeric vector of weighted quantiles with length equal to
#'   \code{length(probs)}.
#'
#' @keywords internal
wquantile <- function(x, w, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) {

  if (na.rm) {
    ok <- is.finite(x) & is.finite(w)
    x <- x[ok]
    w <- w[ok]
  }

  if (length(x) == 0 || sum(w) <= 0)
    return(rep(NA_real_, length(probs)))

  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w[o])

  cw <- cumsum(w)

  sapply(probs, function(p) {
    x[which(cw >= p)[1]]
  })
}

#' Summarize likelihood-weighted SPR trajectories
#'
#' Computes likelihood-weighted summaries of spawning potential ratio (SPR)
#' across Monte Carlo refits. Akaike weights are calculated from relative
#' likelihood support among Monte Carlo iterations.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param years Optional vector of years to summarize. If \code{NULL}, all years
#'   in \code{mc$spr} are used.
#' @param dAIC_cut Optional delta AIC cutoff used before calculating weights.
#'   If \code{NULL}, all iterations are weighted by their relative AIC support.
#' @param probs Numeric vector of weighted quantiles to compute. Default is
#'   \code{c(0.025, 0.5, 0.975)}.
#'
#' @return A data frame with year, weighted mean, lower quantile, median, upper
#'   quantile, number of simulations and effective sample size.
#'
#' @export
mc_spr_summary <- function(mc, years = NULL, dAIC_cut = NULL,
                           probs = c(0.025, 0.5, 0.975)) {

  wtab <- mc_aic_weights(mc, dAIC_cut = dAIC_cut)
  w <- wtab$weight

  yrs <- dimnames(mc$spr)$year
  if (is.null(years)) years <- yrs
  years <- as.character(years)

  out <- data.frame()

  for (yr in years) {

    spr <- as.numeric(an(mc$spr[, ac(yr)]))

    q <- wquantile(spr, w, probs = probs)

    tmp <- data.frame(
      year = yr,
      mean = weighted.mean(spr, w, na.rm = TRUE),
      q025 = q[1],
      median = q[2],
      q975 = q[3],
      n = length(spr),
      n_eff = 1 / sum(w^2, na.rm = TRUE)
    )

    out <- rbind(out, tmp)
  }

  rownames(out) <- NULL
  out
}



#' Plot likelihood-weighted SPR distributions
#'
#' Produces right-sided density plots of likelihood-weighted spawning potential
#' ratio (SPR) distributions across years. The red point and vertical line show
#' the weighted median and uncertainty interval, while the blue star shows the
#' original fitted model stored in \code{mc$fit}.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param years Optional vector of years to plot. If \code{NULL}, all years in
#'   \code{mc$spr} are used.
#' @param dAIC_cut Optional delta AIC cutoff used before calculating Akaike
#'   weights. If \code{NULL}, all iterations are weighted by their relative
#'   AIC support.
#' @param scale Width scaling factor for density polygons. Default is
#'   \code{0.45}.
#' @param n Number of points used to evaluate each density. Default is
#'   \code{300}.
#' @param bw Bandwidth passed to \code{stats::density()}. Setting this explicitly
#'   avoids warnings because automatic bandwidth selection does not use weights.
#'   Default is \code{0.04}.
#' @param trim Relative density cutoff used to remove very small density tails.
#'   For example, \code{trim = 0.02} only draws parts of the density greater
#'   than two percent of the yearly maximum density. Default is \code{0.01}.
#' @param show_fit Logical. If \code{TRUE}, show the original fitted SPR values
#'   from \code{mc$fit$spr}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_mcsprdist <- function(mc, years = NULL, dAIC_cut = 20,
                           scale = 0.45, n = 300,
                           bw = 0.04, trim = 0.01,
                           show_fit = TRUE)  {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  wtab <- mc_aic_weights(mc, dAIC_cut = dAIC_cut)
  w <- wtab$weight

  yrs <- dimnames(mc$spr)$year
  if (is.null(years)) years <- yrs
  years <- as.character(years)

  dens_df <- data.frame()
  summ_df <- mc_spr_summary(mc, years = years, dAIC_cut = dAIC_cut)

  for (j in seq_along(years)) {

    yr <- years[j]
    spr <- as.numeric(an(mc$spr[, ac(yr)]))

    ok <- is.finite(spr) & is.finite(w) & w > 0

    if (sum(ok) < 3) next

    dd <- stats::density(
      spr[ok],
      weights = w[ok] / sum(w[ok]),
      from = 0,
      to = 1,
      n = n,
      bw = bw,
      na.rm = TRUE
    )

    # remove very small density tails
    keep_dens <- dd$y > max(dd$y, na.rm = TRUE) * trim

    x <- dd$x[keep_dens]
    y <- dd$y[keep_dens]

    width <- y / max(y, na.rm = TRUE) * scale

    tmp <- data.frame(
      year = yr,
      year_id = j,
      spr = c(x, rev(x)),
      x = c(rep(j, length(x)), j + rev(width))
    )

    dens_df <- rbind(dens_df, tmp)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = dens_df,
      ggplot2::aes(x = x, y = spr, group = year_id),
      fill = "grey70",
      colour = "grey35",
      alpha = 0.75,
      linewidth = 0.3
    ) +
    ggplot2::geom_linerange(
      data = summ_df,
      ggplot2::aes(x = seq_along(years), ymin = q025, ymax = q975),
      colour = "red",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = summ_df,
      ggplot2::aes(x = seq_along(years), y = median),
      colour = "red",
      size = 2
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(years),
      labels = years
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Year",
      y = "Likelihood-weighted SPR"
    )

  if (show_fit && !is.null(mc$fit) && !is.null(mc$fit$spr)) {

    fit_df <- data.frame(
      year = years,
      year_id = seq_along(years),
      spr = as.numeric(an(mc$fit$spr[, ac(years)]))
    )

    p <- p +
      ggplot2::geom_point(
        data = fit_df,
        ggplot2::aes(x = year_id, y = spr),
        shape = 8,
        colour = "blue",
        size = 2.5,
        stroke = 1
      )
  }

  p
}

#' Plot Monte Carlo life-history parameter draws
#'
#' Produces faceted histograms of life-history parameters sampled during a
#' FLicc Monte Carlo uncertainty run. Red dashed vertical lines show the
#' corresponding parameter values from the original fitted model stored in
#' \code{mc$fit$lhpar}.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param pars Character vector of life-history parameters to plot. Defaults to
#'   \code{c("linf", "Mk", "k", "M", "L50", "CVL")}.
#' @param bins Number of histogram bins. Default is \code{30}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mc <- mc_flicc(fit, nsim = 500)
#' plot_mcpars(mc)
#' }
plot_mcpars <- function(mc,
                          pars = c("linf", "Mk", "k", "M", "L50", "CVL"),
                          bins = 30) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  df <- as.data.frame(mc$lhpar)
  df <- df[df$params %in% pars, ]
  df$value <- as.numeric(df$data)
  df$params <- factor(df$params, levels = pars)

  fit_lh <- mc$fit$lhpar

  ref <- data.frame(
    params = pars,
    value = NA_real_
  )

  for (p in pars) {
    ref$value[ref$params == p] <- if (p == "CVL") {
      fit$settings$CVL
    } else {
      as.numeric(fit_lh[p])
    }
  }

  ref$params <- factor(ref$params, levels = pars)

  ggplot2::ggplot(df, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = "grey75",
      colour = "white"
    ) +
    ggplot2::geom_vline(
      data = ref,
      ggplot2::aes(xintercept = value),
      colour = "blue",
      linewidth = 0.8,
      linetype = 2
    ) +
    ggplot2::facet_wrap(~params, scales = "free", ncol = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Monte Carlo draw",
      y = "Frequency"
    )
}

#' Plot correlations among Monte Carlo life-history parameter draws
#'
#' Produces faceted scatter plots for selected pairs of life-history parameters
#' sampled during a FLicc Monte Carlo uncertainty run. Points are classified by
#' delta AIC support. Fits with \code{dAIC <= dAIC_cut} are shown using a colour
#' gradient, while less-supported fits are shown as open grey circles. The
#' original fitted parameter combination stored in \code{mc$fit$lhpar} is shown
#' as a blue star.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param pairs A list of character vectors of length two defining parameter
#'   pairs to plot. Defaults to \code{list(c("linf", "k"), c("M", "k"),
#'   c("linf", "L50"), c("linf", "CVL"))}.
#' @param dAIC_cut Delta AIC threshold used to define the retained plausible
#'   ensemble. Default is \code{20}.
#' @param ncol Number of columns in the faceted plot. Default is \code{2}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mc <- mc_flicc(fit, nsim = 500)
#' plot_mccor(mc, dAIC_cut = 20)
#' }
plot_mccor <- function(mc,
                          pairs = list(
                            c("linf", "k"),
                            c("M", "k"),
                            c("linf", "L50"),
                            c("linf", "CVL")
                          ),
                          dAIC_cut = 20,
                          ncol = 2) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  mc_dAIC <- function(mc) {
    nll <- -as.numeric(mc$logLik)
    2 * (nll - min(nll, na.rm = TRUE))
  }

  df <- as.data.frame(mc$lhpar)
  df$value <- as.numeric(df$data)

  wide <- reshape(
    df[, c("iter", "params", "value")],
    idvar = "iter",
    timevar = "params",
    direction = "wide"
  )

  names(wide) <- sub("^value\\.", "", names(wide))

  wide$dAIC <- mc_dAIC(mc)
  wide$keep <- wide$dAIC <= dAIC_cut

  pd <- data.frame()
  panel_levels <- character(length(pairs))

  for (j in seq_along(pairs)) {
    z <- pairs[[j]]
    panel_levels[j] <- paste(z[1], "vs", z[2])

    if (!all(z %in% names(wide))) {
      warning("Skipping pair ", paste(z, collapse = " vs "),
              ": parameter not found in mc$lhpar.")
      next
    }

    pd <- rbind(pd, data.frame(
      iter = wide$iter,
      panel = panel_levels[j],
      x = wide[[z[1]]],
      y = wide[[z[2]]],
      dAIC = wide$dAIC,
      keep = wide$keep
    ))
  }

  pd$panel <- factor(pd$panel, levels = panel_levels)

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(
      data = pd[!pd$keep, ],
      shape = 1,
      colour = "grey55",
      alpha = 0.75,
      size = 1.6,
      stroke = 0.45
    ) +
    ggplot2::geom_point(
      data = pd[pd$keep, ],
      ggplot2::aes(colour = dAIC),
      size = 2.1,
      alpha = 0.95
    ) +
    ggplot2::scale_colour_gradientn(
      colours = c("darkred", "orangered", "orange", "gold"),
      limits = c(0, dAIC_cut),
      name = expression(Delta*AIC)
    ) +
    ggplot2::facet_wrap(~panel, scales = "free", ncol = ncol) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL)


   fit_lh <- mc$fit$lhpar


    get_fit_par <- function(z) {
      if (z == "CVL") {
        fit$settings$CVL
      } else {
        as.numeric(fit_lh[z])
      }
    }

    fit_pd <- data.frame()

    for (j in seq_along(pairs)) {
      z <- pairs[[j]]
      fit_pd <- rbind(fit_pd, data.frame(
        panel = panel_levels[j],
        x = get_fit_par(z[1]),
        y = get_fit_par(z[2])
      ))
    }

    fit_pd$panel <- factor(fit_pd$panel, levels = panel_levels)

    p <- p +
      ggplot2::geom_point(
        data = fit_pd,
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        shape = 8,
        colour = "blue",
        size = 3.2,
        stroke = 1.1
      )


  p
}




#' Plot Monte Carlo likelihood profile over SPR
#'
#' Produces a likelihood-profile diagnostic plot for a FLicc Monte Carlo
#' uncertainty run. The x-axis shows spawning potential ratio (SPR) for a
#' selected year and the y-axis shows delta AIC relative to the best-supported
#' Monte Carlo refit. Points within the selected delta AIC threshold are shown
#' using a likelihood-support colour gradient, while less-supported runs are
#' shown as open grey circles. The original fitted model stored in
#' \code{mc$fit} is shown as a blue star.
#'
#' This plot is useful for identifying likelihood-supported regions, weakly
#' supported boundary solutions, and potential confounding between fishing
#' mortality, selectivity and life-history inputs.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param year Year for which SPR should be plotted. If \code{NULL}, the final
#'   year in \code{mc$spr} is used.
#' @param dAIC_cut Delta AIC threshold used to define the retained plausible
#'   ensemble. Default is \code{20}.
#' @param ylim Numeric vector of length two giving the y-axis limits for
#'   delta AIC. Default is \code{c(0, 200)}.
#' @param ref Logical. If \code{TRUE}, show the original fitted model stored in
#'   \code{mc$fit}. Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mc <- mc_flicc(fit, nsim = 500, parallel = TRUE, workers = 8)
#'
#' # Likelihood profile for the final year
#' plot_mcprofile(mc)
#'
#' # Profile for a specific year
#' plot_mcprofile(mc, year = 2024, dAIC_cut = 20)
#' }
plot_mcprofile <- function(mc,
                       year = NULL,
                       dAIC_cut = 20,
                       ylim = c(0, 100),
                       ref = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  yrs <- dimnames(mc$spr)$year

  if (is.null(year)) year <- tail(yrs, 1)
  year <- as.character(year)

  if (!year %in% yrs)
    stop("year not found in mc$spr.")

  nll <- -as.numeric(mc$logLik)
  dAIC <- 2 * (nll - min(nll, na.rm = TRUE))

  spr <- as.numeric(an(mc$spr[, ac(year)]))

  df <- data.frame(
    spr = spr,
    dAIC = dAIC,
    keep = dAIC <= dAIC_cut
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = spr, y = dAIC)) +
    ggplot2::geom_point(
      data = df[!df$keep, ],
      shape = 1,
      colour = "grey55",
      alpha = 0.75,
      size = 1.6,
      stroke = 0.45
    ) +
    ggplot2::geom_point(
      data = df[df$keep, ],
      ggplot2::aes(colour = dAIC),
      size = 2.1,
      alpha = 0.95
    ) +
    ggplot2::scale_colour_gradientn(
      colours = c("darkred", "orangered", "orange", "gold"),
      limits = c(0, dAIC_cut),
      name = expression(Delta*AIC)
    ) +
    ggplot2::geom_hline(
      yintercept = dAIC_cut,
      colour = "red",
      linetype = 2,
      linewidth = 0.7
    ) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste("SPR", year),
      y = expression(Delta*AIC)
    )

  if (ref && !is.null(mc$fit) && !is.null(mc$fit$spr)) {

    fit_spr <- as.numeric(an(mc$fit$spr[, ac(year)]))
    fit_nll <- -as.numeric(mc$fit$logLik)
    fit_dAIC <- 2 * (fit_nll - min(nll, na.rm = TRUE))

    p <- p +
      ggplot2::geom_point(
        data = data.frame(spr = fit_spr, dAIC = fit_dAIC),
        ggplot2::aes(x = spr, y = dAIC),
        inherit.aes = FALSE,
        shape = 8,
        colour = "blue",
        size = 3.4,
        stroke = 1.1
      )
  }

  p
}


#' Summarize Monte Carlo life-history priors
#'
#' Creates a compact summary table of the life-history parameter draws used in a
#' FLicc Monte Carlo uncertainty run. The table includes arithmetic mean,
#' standard deviation, coefficient of variation, log-scale mean and log-scale
#' standard deviation.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param pars Character vector of life-history parameters to summarize.
#'   Defaults to \code{c("linf", "Mk", "k", "M", "L50", "CVL")}.
#'
#' @return A data frame with one row per parameter.
#'
#' @export
mc_lhtab <- function(mc,
                           pars = c("linf", "Mk", "k", "M", "L50", "CVL")) {

  df <- as.data.frame(mc$lhpar)
  df$value <- as.numeric(df$data)
  df <- df[df$params %in% pars, ]

  out <- data.frame()

  for (p in pars) {

    x <- df$value[df$params == p]
    lx <- log(x[x > 0])

    tmp <- data.frame(
      params = p,
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      cv = stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE),
      log.mean = mean(lx, na.rm = TRUE),
      log.sd = stats::sd(lx, na.rm = TRUE)
    )

    out <- rbind(out, tmp)
  }

  rownames(out) <- NULL
  out
}



#' Summarize Monte Carlo SPR uncertainty
#'
#' Produces a year-specific summary table of spawning potential ratio (SPR)
#' uncertainty from a FLicc Monte Carlo ensemble. The table includes the
#' reference fit estimate, unweighted summaries of retained ensemble members,
#' Akaike-weighted summaries, the number of retained realizations, and the
#' effective ensemble size.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param dAIC_cut Delta AIC threshold used to define the retained ensemble.
#'   Default is \code{20}.
#' @param years Optional vector of years to summarize. If \code{NULL}, all years
#'   in \code{mc$spr} are used.
#'
#' @return A data frame with columns \code{year}, \code{mle}, \code{mean_uw},
#'   \code{cv_uw}, \code{mean_w}, \code{cv_w}, \code{n_retained}, and
#'   \code{n_eff}. Attributes \code{dAIC_cut} and \code{n_total} store the AIC
#'   cutoff and total number of Monte Carlo realizations.
#'
#' @details
#' The unweighted mean and coefficient of variation summarize the spread of
#' retained ensemble members after delta AIC filtering. The weighted mean and
#' coefficient of variation use normalized Akaike weights and therefore reflect
#' relative likelihood support. The effective ensemble size is calculated as
#' \deqn{n_{\mathrm{eff}} = 1 / \sum_i w_i^2,}
#' where \eqn{w_i} are the normalized Akaike weights.
#'
#' @export
mc_spr_tab <- function(mc,
                                 dAIC_cut = 20,
                                 years = NULL) {

  if (is.null(years))
    years <- dimnames(mc$spr)$year

  years <- as.character(years)

  wtab <- mc_aic_weights(mc, dAIC_cut = dAIC_cut)

  w <- wtab$weight
  keep <- wtab$keep

  n_eff <- 1 / sum(w^2, na.rm = TRUE)

  out <- vector("list", length(years))

  for (i in seq_along(years)) {

    yr <- years[i]

    spr <- as.numeric(an(mc$spr[, ac(yr)]))

    spr_keep <- spr[keep]
    w_keep   <- w[keep]

    # normalize retained weights
    w_keep <- w_keep / sum(w_keep)

    mle <- if (!is.null(mc$fit))
      as.numeric(an(mc$fit$spr[, ac(yr)]))
    else
      NA_real_

    mean_uw <- mean(spr_keep, na.rm = TRUE)

    cv_uw <- stats::sd(spr_keep, na.rm = TRUE) /
      mean_uw

    mean_w <- weighted.mean(
      spr_keep,
      w_keep,
      na.rm = TRUE
    )

    sd_w <- sqrt(
      weighted.mean(
        (spr_keep - mean_w)^2,
        w_keep,
        na.rm = TRUE
      )
    )

    cv_w <- sd_w / mean_w

    out[[i]] <- data.frame(
      year = yr,
      mle = mle,
      mean_uw = mean_uw,
      cv_uw = cv_uw,
      mean_w = mean_w,
      cv_w = cv_w,
      n_retained = sum(keep),
      n_eff = n_eff
    )
  }

  out <- do.call(rbind, out)
  attr(out, "dAIC_cut") <- dAIC_cut
  attr(out, "n_total") <- length(w)
  out
}

#' Plot weighted and unweighted SPR densities
#'
#' Produces faceted density plots comparing unweighted and Akaike-weighted
#' spawning potential ratio (SPR) distributions across years. The unweighted
#' density is calculated from the retained Monte Carlo ensemble, while the
#' weighted density uses normalized Akaike weights. The reference FLicc fit
#' stored in \code{mc$fit} is shown as a blue dashed vertical line.
#'
#' @param mc An object returned by \code{mc_flicc()}.
#' @param years Optional vector of years to plot. If \code{NULL}, all years in
#'   \code{mc$spr} are used.
#' @param dAIC_cut Delta AIC threshold used to define the retained ensemble.
#'   Default is \code{20}.
#' @param n Number of points used to evaluate each density. Default is
#'   \code{300}.
#' @param bw Bandwidth passed to \code{stats::density()}. Default is
#'   \code{0.04}.
#' @param show_fit Logical. If \code{TRUE}, show the reference fit stored in
#'   \code{mc$fit}. Default is \code{TRUE}.
#' @param ncol Number of columns in the faceted plot. Default is \code{3}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' This plot is intended as a diagnostic for assessing how likelihood weighting
#' changes the retained structural ensemble. The unweighted density represents
#' the spread of all retained realizations after the \code{dAIC_cut} filter,
#' while the weighted density emphasizes the subset of realizations with higher
#' relative support. Large differences between the two densities indicate that
#' the likelihood is strongly concentrating support within the plausible
#' structural ensemble.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mc <- mc_flicc(fit, nsim = 500, parallel = TRUE, workers = 8)
#'
#' plot_mcspr_weights(mc)
#'
#' plot_mcsprwts(mc, years = 2020:2024, dAIC_cut = 20)
#' }
plot_mcsprwts<- function(mc, years = NULL, dAIC_cut = 20,
                                 n = 300, bw = 0.04,
                                 show_fit = TRUE, ncol = 3) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  yrs <- dimnames(mc$spr)$year
  if (is.null(years)) years <- yrs
  years <- as.character(years)

  wtab <- mc_aic_weights(mc, dAIC_cut = dAIC_cut)
  w <- wtab$weight
  keep <- wtab$keep

  dens_df <- data.frame()

  for (yr in years) {

    spr <- as.numeric(an(mc$spr[, ac(yr)]))
    ok <- is.finite(spr) & keep

    if (sum(ok) < 3) next

    dd_uw <- stats::density(
      spr[ok],
      from = 0, to = 1,
      n = n, bw = bw,
      na.rm = TRUE
    )

    dens_df <- rbind(dens_df, data.frame(
      year = yr,
      spr = dd_uw$x,
      density = dd_uw$y,
      type = "Unweighted"
    ))

    ok_w <- ok & is.finite(w) & w > 0

    if (sum(ok_w) >= 3) {
      dd_w <- stats::density(
        spr[ok_w],
        weights = w[ok_w] / sum(w[ok_w]),
        from = 0, to = 1,
        n = n, bw = bw,
        na.rm = TRUE
      )

      dens_df <- rbind(dens_df, data.frame(
        year = yr,
        spr = dd_w$x,
        density = dd_w$y,
        type = "Weighted"
      ))
    }
  }

  dens_df$year <- factor(dens_df$year, levels = years)
  dens_df$type <- factor(dens_df$type, levels = c("Unweighted", "Weighted"))

  p <- ggplot2::ggplot(
    dens_df,
    ggplot2::aes(x = spr, y = density, colour = type, fill = type)
  ) +
    ggplot2::geom_area(alpha = 0.35, position = "identity") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(
      values = c(Unweighted = "grey35", Weighted = "grey8"),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c(Unweighted = "grey75", Weighted = "grey50"),
      name = NULL
    ) +
    ggplot2::facet_wrap(~year, ncol = ncol) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "SPR",
      y = "Density"
    )

  if (show_fit && !is.null(mc$fit) && !is.null(mc$fit$spr)) {

    fit_df <- data.frame(
      year = factor(years, levels = years),
      spr = as.numeric(an(mc$fit$spr[, ac(years)]))
    )

    p <- p +
      ggplot2::geom_vline(
        data = fit_df,
        ggplot2::aes(xintercept = spr),
        inherit.aes = FALSE,
        colour = "blue",
        linetype = 2,
        linewidth = 0.9
      )
  }

  p
}
