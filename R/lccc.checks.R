# Length-converted catch-curve checks for FLicc
#
# This file contains lightweight diagnostic helpers for estimating empirical
# length-converted catch-curve (LCCC) total mortality, Z, from FLQuant length
# compositions. The functions are intended as checks/diagnostics and should not
# replace the full FLicc likelihood-based fits.

#' Extract a life-history parameter as numeric
#'
#' Small internal helper used by the LCCC functions. It works with named numeric
#' vectors, FLPar-like objects, and simple lists.
#'
#' @param lhpar Named life-history parameter object.
#' @param name Character. Parameter name to extract.
#' @param default Optional numeric default used when the parameter is absent.
#'
#' @return A numeric scalar.
#' @keywords internal
.get_lhpar <- function(lhpar, name, default = NULL) {

  out <- NULL

  if (is.list(lhpar) && !is.null(lhpar[[name]])) {
    out <- lhpar[[name]]
  } else {
    nms <- try(names(lhpar), silent = TRUE)
    if (!inherits(nms, "try-error") && name %in% nms) {
      out <- lhpar[name]
    }
  }

  if (is.null(out)) {
    if (!is.null(default)) return(as.numeric(default)[1])
    stop("Parameter '", name, "' not found in 'lhpar'.", call. = FALSE)
  }

  out <- as.numeric(out)
  out <- out[is.finite(out)]

  if (length(out) == 0) {
    if (!is.null(default)) return(as.numeric(default)[1])
    stop("Parameter '", name, "' is missing or not finite.", call. = FALSE)
  }

  out[1]
}

#' Extract length classes from an FLQuant length object
#'
#' @param lfd FLQuant with length classes in the first dimension. The function
#'   first looks for a dimension named `len`; if not found, it uses the first
#'   dimension.
#'
#' @return Numeric vector of length-class labels.
#' @keywords internal
.get_lens <- function(lfd) {

  dmn <- dimnames(lfd)

  if (!is.null(dmn$len)) {
    len <- as.numeric(dmn$len)
  } else {
    len <- as.numeric(dmn[[1]])
  }

  if (any(!is.finite(len))) {
    stop("Could not extract finite numeric length classes from 'lfd'.", call. = FALSE)
  }

  len
}

#' Estimate an empirical lower cutoff for a length-converted catch curve
#'
#' Estimates the lower length cutoff used to fit the descending limb of a
#' length-converted catch curve. The cutoff is based on the empirical modal
#' length class. If the distribution is bimodal, the second local maximum in
#' length order can be selected, which is often preferable when the first peak
#' represents partial recruitment, juveniles, or a recruitment pulse.
#'
#' The function collapses the length composition over all non-length dimensions,
#' converts the result to relative frequencies, optionally smooths the relative
#' frequency curve, finds local maxima, and returns `lc = lc_mult * Lpeak`.
#'
#' @param lfd FLQuant length-frequency object.
#' @param smooth Logical. Should relative frequencies be smoothed before finding
#'   peaks? Default is `TRUE`.
#' @param span Numeric smoothing span passed to `stats::loess` when
#'   `smooth = TRUE`.
#' @param second_if_bimodal Logical. If two or more valid local maxima are
#'   detected, should the second peak in length order be used? Default is `TRUE`.
#' @param peak_rel_min Numeric between 0 and 1. Minimum relative height for a
#'   local maximum to be considered a peak.
#' @param lc_mult Numeric multiplier applied to the selected peak to define the
#'   lower cutoff. The default is `0.9`, i.e. the catch curve uses length classes
#'   larger than 90 percent of the selected peak length.
#'
#' @return A list of class `flicc_empirical_lc` with elements `lc`, `Lpeak`,
#'   `peak_index`, `len`, `observed`, `smoothed`, and `peaks`.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' elc <- empirical_lc(lfd)
#' elc$lc
#' plot_empirical_lc(elc)
#' }
#'
#' @export
empirical_lc <- function(lfd,
                         smooth = TRUE,
                         span = 0.35,
                         second_if_bimodal = TRUE,
                         peak_rel_min = 0.20,
                         lc_mult = 0.90) {

  len <- .get_lens(lfd)

  arr <- as.array(lfd)
  y <- apply(arr, 1, sum, na.rm = TRUE)

  ok <- is.finite(len) & is.finite(y) & y > 0
  len <- len[ok]
  y <- y[ok]

  if (length(y) < 4) {
    stop("Too few positive length bins to estimate an empirical cutoff.", call. = FALSE)
  }

  y <- y / max(y, na.rm = TRUE)

  if (smooth && length(y) >= 6) {
    ys <- try(
      as.numeric(stats::predict(stats::loess(y ~ len, span = span, degree = 2),
                                newdata = data.frame(len = len))),
      silent = TRUE
    )

    if (inherits(ys, "try-error") || length(ys) != length(y)) {
      ys <- y
    }

    ys[!is.finite(ys)] <- y[!is.finite(ys)]
    ys <- pmax(ys, 0)
  } else {
    ys <- y
  }

  ## Local maxima. This is intentionally simple and transparent.
  is_peak <- c(FALSE, diff(sign(diff(ys))) < 0, FALSE)
  peaks <- which(is_peak & ys >= peak_rel_min)

  if (length(peaks) == 0) {
    peak <- which.max(ys)
  } else {
    peaks <- peaks[order(len[peaks])]

    if (second_if_bimodal && length(peaks) >= 2) {
      peak <- peaks[2]
    } else {
      peak <- peaks[which.max(ys[peaks])]
    }
  }

  Lpeak <- len[peak]
  lc <- lc_mult * Lpeak

  out <- list(
    lc = as.numeric(lc),
    Lpeak = as.numeric(Lpeak),
    peak_index = as.integer(peak),
    len = as.numeric(len),
    observed = as.numeric(y),
    smoothed = as.numeric(ys),
    peaks = as.numeric(len[peaks]),
    lc_mult = as.numeric(lc_mult),
    peak_rel_min = as.numeric(peak_rel_min),
    second_if_bimodal = isTRUE(second_if_bimodal)
  )

  class(out) <- "flicc_empirical_lc"
  out
}

#' Plot empirical selection of the LCCC lower cutoff
#'
#' Visualises the relative length-frequency distribution used by
#' [empirical_lc()]. The selected peak (`Lpeak`) and the actual catch-curve
#' cutoff (`lc`) are shown as vertical lines. The catch curve uses length bins
#' larger than `lc`, not larger than `Lpeak`.
#'
#' @param x Object returned by [empirical_lc()].
#' @param show_peaks Logical. Should all detected local peaks be shown as points?
#' @param alpha_col Numeric transparency for the observed relative frequency bars.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' elc <- empirical_lc(lfd)
#' plot_empirical_lc(elc)
#' }
#'
#' @export
plot_empirical_lc <- function(x, show_peaks = TRUE, alpha_col = 0.35) {

  if (!inherits(x, "flicc_empirical_lc")) {
    stop("'x' must be the output from empirical_lc().", call. = FALSE)
  }

  dat <- data.frame(
    len = as.numeric(x$len),
    observed = as.numeric(x$observed),
    smoothed = as.numeric(x$smoothed)
  )

  dat <- dat[is.finite(dat$len) & is.finite(dat$observed), ]

  vdat <- data.frame(
    xintercept = c(as.numeric(x$Lpeak), as.numeric(x$lc)),
    label = c("selected peak", "catch-curve cutoff")
  )

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = len)) +
    ggplot2::geom_col(ggplot2::aes(y = observed), alpha = alpha_col) +
    ggplot2::geom_line(ggplot2::aes(y = smoothed), linewidth = 1) +
    ggplot2::geom_vline(
      data = vdat,
      ggplot2::aes(xintercept = xintercept, linetype = label),
      linewidth = 0.8
    ) +
    ggplot2::labs(
      x = "Length",
      y = "Relative frequency",
      linetype = NULL,
      title = "Empirical catch-curve cutoff",
      subtitle = "Lpeak is the selected mode; lc is the lower cutoff used for the LCCC"
    ) +
    ggplot2::theme_bw()

  if (show_peaks && !is.null(x$peaks) && length(x$peaks) > 0) {
    peak_dat <- data.frame(len = as.numeric(x$peaks))
    peak_dat$smoothed <- stats::approx(
      x = dat$len,
      y = dat$smoothed,
      xout = peak_dat$len,
      rule = 2
    )$y

    p <- p +
      ggplot2::geom_point(
        data = peak_dat,
        ggplot2::aes(x = len, y = smoothed),
        size = 2
      )
  }

  p
}


#' Build regression data for a length-converted catch curve
#'
#' Constructs the regression data used by a Pauly-style length-converted
#' catch curve. The response is log(N / Delta t), where N is the observed
#' length-frequency count and Delta t is the time required to grow through
#' the length interval. The predictor is relative age from the inverse
#' von Bertalanffy curve. The intercept is arbitrary; Z is estimated from
#' the negative slope.
#'
#' @param lfd An FLQuant with length in the first dimension.
#' @param lhpar Life-history parameters containing at least \code{linf}
#'   and \code{k}.
#' @param lc Optional lower length cutoff. If \code{NULL}, estimated with
#'   \code{empirical_lc()}.
#' @param empirical_lc_args Optional list of arguments passed to
#'   \code{empirical_lc()}.
#' @param max_rel_linf Maximum allowed length relative to Linf.
#'
#' @return A data.frame with columns \code{len}, \code{year}, \code{unit},
#'   \code{season}, \code{area}, \code{iter}, \code{n}, \code{dt},
#'   \code{y}, and \code{x}.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' dat <- zlcc_data(lfd, lhpar)
#' head(dat)
#' plot_zlcc(lfd, lhpar)
#' }
#'
#' @export
zlcc_data <- function(lfd,
                      lhpar,
                      lc = NULL,
                      empirical_lc_args = list(),
                      max_rel_linf = 0.999) {

  linf <- an(lhpar["linf"])
  k <- an(lhpar["k"])

  if (is.null(lc)) {
    elc <- do.call(empirical_lc, c(list(lfd = lfd), empirical_lc_args))
    lc <- elc$lc
  } else {
    elc <- NULL
  }

  len0 <- as.numeric(dimnames(lfd)$len)
  keep <- len0 > lc & len0 < linf * max_rel_linf

  lfd <- lfd[keep,]

  if (length(as.numeric(dimnames(lfd)$len)) < 4) {
    stop("Too few length bins after applying lc and linf filters.", call. = FALSE)
  }

  ## Midpoint lengths for retained bins
  l <- lfd
  l[] <- as.numeric(dimnames(l)[[1]])

  ## Adjacent lower and upper midpoint lengths
  l1 <- lfd[-dim(lfd)[1]]
  l1[] <- as.numeric(dimnames(l1)[[1]])

  l2 <- lfd[-1]
  l2[] <- as.numeric(dimnames(l2)[[1]])

  ## Time needed to grow from one midpoint to the next.
  ## t0 cancels because this is a relative-age difference.
  dt <- -log(1 - l2 / linf) / k +
    log(1 - l1 / linf) / k

  ## Relative age at length midpoint.
  t <- -log(1 - l / linf) / k

  ## Align counts with dt: omit final length bin.
  n <- lfd[-dim(lfd)[1]]

  ## Build model.frame before taking logs.
  dat <- model.frame(FLQuants(
    n  = n,
    dt = dt,
    x  = t[-dim(t)[1]]
  ))

  ## Force plain numeric vectors.
  dat$n <- as.numeric(dat$n)
  dat$dt <- as.numeric(dat$dt)
  dat$x <- as.numeric(dat$x)

  ## Keep only usable positive observations.
  dat <- dat[
    is.finite(dat$n)  & dat$n  > 0 &
      is.finite(dat$dt) & dat$dt > 0 &
      is.finite(dat$x),
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) {
    stop(
      "No positive length-frequency observations remain after lc/Linf filtering. ",
      "Try lowering lc, setting lc manually, or checking empirical_lc(lfd).",
      call. = FALSE
    )
  }

  dat$y <- log(dat$n / dat$dt)

  dat <- dat[is.finite(dat$y) & is.finite(dat$x), , drop = FALSE]

  if (nrow(dat) == 0) {
    stop("No finite LCCC regression observations after log(N / dt).", call. = FALSE)
  }

  if ("year" %in% names(dat)) dat$year <- factor(dat$year)
  if ("iter" %in% names(dat)) dat$iter <- factor(dat$iter)

  if (!"year" %in% names(dat)) dat$year <- factor("all")
  if (!"iter" %in% names(dat)) dat$iter <- factor("1")

  attr(dat, "lc") <- as.numeric(lc)
  attr(dat, "empirical_lc") <- elc
  attr(dat, "lhpar") <- lhpar
  attr(dat, "linf") <- linf
  attr(dat, "k") <- k

  dat
}

#' Fit length-converted catch-curve regressions
#'
#' Fits one linear regression by year and iteration using data returned by
#' [zlcc_data()]. Total mortality is estimated as `Z = -slope`.
#'
#' @param dat Data.frame returned by [zlcc_data()].
#' @param min_n Integer. Minimum number of points required to fit each regression.
#'
#' @return Data.frame with columns `year`, `iter`, `intercept`, `slope`, `Z`,
#'   `r2`, and `n`.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' dat <- zlcc_data(lfd, lhpar)
#' zlcc_fits(dat)
#' }
#'
#' @export
zlcc_fits <- function(dat, min_n = 3) {

  if (!all(c("x", "y", "year", "iter") %in% names(dat))) {
    stop("'dat' must contain x, y, year and iter columns.", call. = FALSE)
  }

  spl <- split(dat, interaction(dat$year, dat$iter, drop = TRUE))

  fits <- do.call(rbind, lapply(spl, function(d) {
    if (nrow(d) < min_n) return(NULL)

    fit <- stats::lm(y ~ x, data = d)
    cf <- stats::coef(fit)

    data.frame(
      year = as.character(d$year[1]),
      iter = as.character(d$iter[1]),
      intercept = unname(cf[1]),
      slope = unname(cf[2]),
      Z = -unname(cf[2]),
      r2 = summary(fit)$r.squared,
      n = nrow(d)
    )
  }))

  if (is.null(fits) || nrow(fits) == 0) {
    stop("No valid LCCC regressions could be fitted.", call. = FALSE)
  }

  fits$year <- factor(fits$year, levels = levels(dat$year))
  fits$iter <- factor(fits$iter, levels = levels(dat$iter))

  fits
}

#' Estimate total mortality from length-converted catch curves
#'
#' Estimates total mortality, Z, by fitting Pauly-style length-converted catch
#' curves independently by year and iteration. This is intended as a diagnostic
#' check against FLicc fits, not as a replacement for the full FLicc model.
#'
#' @param lfd FLQuant length-frequency object.
#' @param lhpar Named life-history parameter object containing at least `linf`
#'   and `k`.
#' @param lc Optional numeric lower cutoff. If `NULL`, [empirical_lc()] is used.
#' @param empirical_lc_args Optional list of arguments passed to [empirical_lc()]
#'   when `lc = NULL`.
#' @param min_n Integer. Minimum number of points required for each regression.
#' @param max_rel_linf Numeric. Length classes greater than this fraction of
#'   `linf` are excluded.
#'
#' @return FLQuant containing Z by year and iteration.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' z <- zlcc(lfd, lhpar)
#' plot(z)
#' plot_zlcc_year(z)
#' }
#'
#' @export
zlcc <- function(lfd, lhpar, lc = NULL, empirical_lc_args = list()) {

  dat <- zlcc_data(
    lfd = lfd,
    lhpar = lhpar,
    lc = lc,
    empirical_lc_args = empirical_lc_args
  )

  if (nrow(dat) == 0) {
    stop("No valid data available for length-converted catch curve.", call. = FALSE)
  }

  spl <- split(dat, interaction(dat$year, dat$iter, drop = TRUE))

  z <- do.call(rbind, lapply(spl, function(d) {

    if (nrow(d) < 3) return(NULL)

    fit <- lm(y ~ x, data = d)

    data.frame(
      year = d$year[1],
      iter = d$iter[1],
      data = -coef(fit)[["x"]]
    )
  }))

  if (is.null(z) || nrow(z) == 0) {
    stop("No valid regressions could be fitted.", call. = FALSE)
  }

  z$year <- factor(z$year)
  z$iter <- factor(z$iter)

  as.FLQuant(z)
}

#' Plot length-converted catch-curve regressions
#'
#' Produces diagnostic plots for the linear regressions used to estimate
#' total mortality, Z, from a length-converted catch curve. The response is
#' log(N / Delta t), and the predictor is relative age from the inverse
#' von Bertalanffy growth curve. Z is estimated as the negative slope of
#' the regression.
#'
#' @param lfd An FLQuant with length in the first dimension.
#' @param lhpar Life-history parameters containing at least \code{linf}
#'   and \code{k}.
#' @param lc Optional lower length cutoff for the descending limb. If
#'   \code{NULL}, this is estimated with \code{empirical_lc()}.
#' @param empirical_lc_args Optional list of arguments passed to
#'   \code{empirical_lc()}.
#' @param min_n Minimum number of observations required to fit a regression.
#' @param max_rel_linf Maximum allowed length relative to Linf.
#' @param facet Logical. Should regressions be shown in facets?
#' @param show_points Logical. Show observed regression points.
#' @param show_lm Logical. Show fitted linear regression line.
#' @param show_eq Logical. Add labels with Z, R-squared and n.
#' @param ncol Optional number of facet columns.
#' @param return_data Logical. If TRUE, return a list with the plot,
#'   regression data and fitted Z table.
#'
#' @return A ggplot object, or a list if \code{return_data = TRUE}.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' plot_zlcc(lfd, lhpar)
#'
#' elc <- empirical_lc(lfd)
#' plot_zlcc(lfd, lhpar, lc = elc$lc)
#'
#' out <- plot_zlcc(lfd, lhpar, return_data = TRUE)
#' out$fits
#' }
#'
#' @export
plot_zlcc <- function(lfd,
                      lhpar,
                      lc = NULL,
                      empirical_lc_args = list(),
                      min_n = 3,
                      max_rel_linf = 0.999,
                      facet = TRUE,
                      show_points = TRUE,
                      show_lm = TRUE,
                      show_eq = TRUE,
                      ncol = NULL,
                      return_data = FALSE) {

  dat <- zlcc_data(
    lfd = lfd,
    lhpar = lhpar,
    lc = lc,
    empirical_lc_args = empirical_lc_args,
    max_rel_linf = max_rel_linf
  )

  fits <- zlcc_fits(dat, min_n = min_n)

  if (is.null(fits) || nrow(fits) == 0) {
    stop("No valid LCCC regressions could be fitted.", call. = FALSE)
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y))

  if (show_points) {
    p <- p +
      ggplot2::geom_point(alpha = 0.65, size = 1.8)
  }

  if (show_lm) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      )
  }

  p <- p +
    ggplot2::labs(
      x = "Relative age from inverse von Bertalanffy",
      y = expression(log(N / Delta*t)),
      title = "Length-converted catch curve",
      subtitle = "Z is estimated as the negative regression slope"
    ) +
    ggplot2::theme_bw()

  if (facet) {
    if ("iter" %in% names(dat) && length(unique(dat$iter)) > 1) {
      p <- p +
        ggplot2::facet_wrap(
          ~ year + iter,
          scales = "free_y",
          ncol = ncol
        )
    } else {
      p <- p +
        ggplot2::facet_wrap(
          ~ year,
          scales = "free_y",
          ncol = ncol
        )
    }
  }

  if (show_eq) {
    fits$label <- paste0(
      "Z = ", round(fits$Z, 3),
      "\nR² = ", round(fits$r2, 2),
      "\nn = ", fits$n
    )

    p <- p +
      ggplot2::geom_label(
        data = fits,
        ggplot2::aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.1,
        size = 3,
        linewidth = 0.15,
        alpha = 0.85
      )
  }

  if (return_data) {
    return(list(
      plot = p,
      data = dat,
      fits = fits,
      lc = attr(dat, "lc"),
      empirical_lc = attr(dat, "empirical_lc")
    ))
  }

  p
}

#' Plot length-converted catch-curve Z by year
#'
#' Plots annual LCCC estimates of total mortality, Z. The input can either be an
#' FLQuant returned by [zlcc()] or a data.frame returned by [zlcc_fits()].
#'
#' @param z FLQuant returned by [zlcc()] or data.frame returned by [zlcc_fits()].
#' @param points Logical. Should points be drawn?
#' @param line Logical. Should lines be drawn?
#' @param facet_iter Logical. Should iterations be faceted when multiple iters
#'   are present?
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lfd <- lfd_alfonsino$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' z <- zlcc(lfd, lhpar)
#' plot_zlcc_year(z)
#' }
#'
#' @export
plot_zlcc_year <- function(z,
                           points = TRUE,
                           line = TRUE,
                           facet_iter = TRUE) {

  if (inherits(z, "FLQuant")) {
    dat <- as.data.frame(z)
    names(dat)[names(dat) == "data"] <- "Z"
  } else if (is.data.frame(z) && "Z" %in% names(z)) {
    dat <- z
  } else {
    stop("'z' must be an FLQuant from zlcc() or a data.frame from zlcc_fits().",
         call. = FALSE)
  }

  if (!"year" %in% names(dat)) stop("No 'year' column found.", call. = FALSE)
  if (!"iter" %in% names(dat)) dat$iter <- factor("1")

  dat$year_num <- suppressWarnings(as.numeric(as.character(dat$year)))
  if (any(!is.finite(dat$year_num))) {
    dat$year_num <- seq_len(nrow(dat))
  }

  dat$Z <- as.numeric(dat$Z)
  dat$iter <- factor(dat$iter)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = year_num, y = Z, group = iter))

  if (line) {
    p <- p + ggplot2::geom_line(linewidth = 0.8)
  }

  if (points) {
    p <- p + ggplot2::geom_point(size = 2)
  }

  p <- p +
    ggplot2::labs(
      x = "Year",
      y = "Total mortality, Z",
      title = "Length-converted catch-curve Z by year"
    ) +
    ggplot2::theme_bw()

  if (facet_iter && length(unique(dat$iter)) > 1) {
    p <- p + ggplot2::facet_wrap(~ iter)
  }

  p
}


#' Compare empirical LCCC Z with FLicc model-implied Z
#'
#' Compares empirical total mortality estimated from length-converted
#' catch-curve regressions with the FLicc model-implied total mortality
#' proxy M + Fap. The empirical Z is gear-specific if a gear-specific
#' length-frequency object is supplied.
#'
#' @param fit A fitted FLicc object containing \code{fit$stklen} and
#'   \code{fit$report$Fap}.
#' @param lfd Gear-specific FLQuant length-frequency data, e.g.
#'   \code{lfd_alfonsino$Gillnet}.
#' @param lc Optional lower length cutoff for LCCC. If \code{NULL}, estimated
#'   using \code{empirical_lc()}.
#' @param empirical_lc_args Optional list of arguments passed to
#'   \code{empirical_lc()}.
#' @param max_rel_linf Maximum allowed length relative to Linf.
#' @param min_n Minimum number of observations for each LCCC regression.
#'
#' @return A list containing empirical Z, FLicc Z, merged comparison table,
#'   regression data, lc and empirical_lc.
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#'
#' lfds <- lfd_alfonsino
#' lfd <- lfds$Gillnet
#' lhpar <- lhpar_alfonsino
#'
#' stklen <- stocklen(lfds, lhpar)
#' fit <- fiticc(
#'   lfds, stklen,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#'
#' stkfit <- flicc_stklen(fit)
#'
#' cmp <- compare_zlcc_flicc(fit, lfd =lfd)
#' cmp$table
#' plot_compare_zlcc_flicc(cmp, quantity = "Z")+ylim(0,NA)
#' }
#'
#' @export
compare_zlcc_flicc <- function(fit,
                               lfd,
                               lc = NULL,
                               empirical_lc_args = list(),
                               max_rel_linf = 0.999,
                               min_n = 3,
                               weight = c("lfd", "equal")) {

  weight <- match.arg(weight)

  if (is.null(fit$stklen)) {
    stop("fit$stklen not found.", call. = FALSE)
  }

  lhpar <- fit$stklen@lhpar
  linf <- an(lhpar["linf"])

  dat <- zlcc_data(
    lfd = lfd,
    lhpar = lhpar,
    lc = lc,
    empirical_lc_args = empirical_lc_args,
    max_rel_linf = max_rel_linf
  )

  z_emp <- zlcc_fits(dat, min_n = min_n)

  z_emp$year <- as.character(z_emp$year)
  z_emp$iter <- as.character(z_emp$iter)

  lc_used <- attr(dat, "lc")

  ## FLicc model-implied Z-at-length
  fitstk <- flicc_stklen(fit)
  zmod <- z(fitstk)

  len <- as.numeric(dimnames(zmod)$len)
  keep <- len > lc_used & len < linf * max_rel_linf

  zmod <- zmod[keep]

  ## Match observed LFD to same range for optional weighting
  lfd_use <- lfd[as.numeric(dimnames(lfd)$len) %in% as.numeric(dimnames(zmod)$len)]

  zdat <- model.frame(FLQuants(
    Z_flicc_len = zmod,
    w = lfd_use
  ))

  zdat$year <- as.character(zdat$year)
  zdat$iter <- as.character(zdat$iter)
  zdat$Z_flicc_len <- as.numeric(zdat$Z_flicc_len)
  zdat$w <- as.numeric(zdat$w)

  zdat <- zdat[is.finite(zdat$Z_flicc_len), , drop = FALSE]

  spl <- split(zdat, interaction(zdat$year, zdat$iter, drop = TRUE))

  z_flicc <- do.call(rbind, lapply(spl, function(d) {

    if (weight == "lfd" && sum(d$w, na.rm = TRUE) > 0) {
      zz <- weighted.mean(d$Z_flicc_len, w = d$w, na.rm = TRUE)
    } else {
      zz <- mean(d$Z_flicc_len, na.rm = TRUE)
    }

    data.frame(
      year = d$year[1],
      iter = d$iter[1],
      Z_flicc = zz
    )
  }))

  tab <- merge(
    z_emp,
    z_flicc,
    by = c("year", "iter"),
    all.x = TRUE
  )

  M <- an(lhpar["M"])

  tab$F_lccc <- pmax(tab$Z - M, 0)
  tab$F_flicc <- pmax(tab$Z_flicc - M, 0)
  tab$Z_ratio <- tab$Z / tab$Z_flicc
  tab$F_ratio <- tab$F_lccc / tab$F_flicc

  list(
    empirical = z_emp,
    model = z_flicc,
    table = tab,
    data = dat,
    z_len = zdat,
    lc = lc_used,
    empirical_lc = attr(dat, "empirical_lc"),
    weight = weight
  )
}


#' Plot empirical LCCC Z/F against FLicc model-implied Z/F
#'
#' @param x Output from \code{compare_zlcc_flicc()}.
#' @param quantity Either \code{"Z"} or \code{"F"}.
#'
#' @return A ggplot object.
#'
#' @export
plot_compare_zlcc_flicc <- function(x, quantity = c("Z", "F")) {

  quantity <- match.arg(quantity)

  tab <- x$table
  tab$year <- as.numeric(as.character(tab$year))

  if (quantity == "Z") {

    pdat <- rbind(
      data.frame(
        year = tab$year,
        value = tab$Z,
        source = "LCCC empirical Z"
      ),
      data.frame(
        year = tab$year,
        value = tab$Z_flicc,
        source = "FLicc length-range Z"
      )
    )

    ylab <- expression(Z~(yr^{-1}))
    title <- "Empirical LCCC Z versus FLicc model-implied Z"

  } else {

    pdat <- rbind(
      data.frame(
        year = tab$year,
        value = tab$F_lccc,
        source = "LCCC implied F = max(Z - M, 0)"
      ),
      data.frame(
        year = tab$year,
        value = tab$F_flicc,
        source = "FLicc length-range implied F"
      )
    )

    ylab <- expression(F~(yr^{-1}))
    title <- "Empirical LCCC-implied F versus FLicc model-implied F"
  }

  ggplot2::ggplot(pdat, ggplot2::aes(x = year, y = value, linetype = source)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Year",
      y = ylab,
      linetype = NULL,
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

