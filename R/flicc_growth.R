# Growth diagnostics: growth-type groups (plot_growth), a MULTIFAN-type phase
# fit of the cohort modes (growth_phase_flicc) and cohort tracks over the
# length compositions (plot_growth_lfd)

# growth parameters from a fit, an FLStockLen or an FLPar
growth_pars_flicc <- function(x, CVL = NULL, maxsd = NULL, ngtg = NULL) {
  getp <- function(lh, nm) {
    pn <- dimnames(lh)[[1]]
    if (nm %in% pn) as.numeric(lh[nm]) else NA_real_
  }
  lfd <- NULL
  if (is.list(x) && !is.null(x$tmb_data) && !is.null(x$report)) {
    lh <- x$report$lhpar
    linf <- getp(lh, "linf"); k <- getp(lh, "k"); L50 <- getp(lh, "L50")
    Mk <- getp(lh, "Mk")
    if (is.null(CVL)) CVL <- tryCatch(as.numeric(x$report$pars["CVL"]), error = function(e) NA_real_)
    if (!is.finite(CVL)) CVL <- x$settings$CVL
    z <- as.numeric(x$tmb_data$gtg_z)
    if (!is.null(maxsd) || !is.null(ngtg)) {
      ms <- if (is.null(maxsd)) max(abs(z)) else maxsd
      ng <- if (is.null(ngtg)) length(z) else ngtg
      z <- if (ng == 1) 0 else seq(-ms, ms, length.out = ng)
    }
    lfd <- x$report$obslen
  } else {
    lh <- if (methods::is(x, "FLStockLen")) x@lhpar else x
    if (!methods::is(lh, "FLPar")) stop("x must be a fitted FLicc object, an FLStockLen or an FLPar.")
    linf <- getp(lh, "linf"); k <- getp(lh, "k"); L50 <- getp(lh, "L50")
    Mk <- getp(lh, "Mk")
    if (!is.finite(Mk)) Mk <- getp(lh, "M") / k
    if (is.null(CVL)) CVL <- 0.1
    ms <- if (is.null(maxsd)) 3 else maxsd
    ng <- if (is.null(ngtg)) 13 else ngtg
    z <- if (ng == 1) 0 else seq(-ms, ms, length.out = ng)
  }
  if (!is.finite(linf) || !is.finite(k)) stop("linf and k are needed.")
  w <- stats::dnorm(z); w <- w / sum(w)
  list(linf = linf, k = k, L50 = L50, Mk = Mk, CVL = CVL, z = z, w = w,
       linf_g = linf * (1 + CVL * z), cv_real = CVL * sqrt(sum(w * z^2)), lfd = lfd)
}

# counts by length (rows) and year (columns) of one FLQuant
lfd_matrix_flicc <- function(q) {
  a <- q@.Data
  m <- apply(a, 1:2, sum, na.rm = TRUE)
  L <- as.numeric(dimnames(q)[[1]])
  dimnames(m) <- list(len = L, year = dimnames(q)$year)
  m
}

# largest length class (lower bound) with fish, by gear
max_len_flicc <- function(lfd) {
  vapply(names(lfd), function(g) {
    m <- lfd_matrix_flicc(lfd[[g]])
    L <- as.numeric(rownames(m)); n <- rowSums(m)
    if (any(n > 0)) max(L[n > 0]) else NA_real_
  }, numeric(1))
}

#' Growth of the growth-type groups
#'
#' Draws the von Bertalanffy growth curves of the growth-type groups (GTG) of
#' FLicc's population model, as in Fig. 1 of Hordyk et al. (2015): one curve
#' per group, \eqn{L_g(a) = L_{\infty,g}(1 - e^{-k(a - t_0)})} with
#' \eqn{L_{\infty,g} = L_\infty(1 + \mathrm{CVL}\, z_g)} on the grid
#' `seq(-maxsd, maxsd, length = ngtg)`, the line shade showing each group's share
#' of the recruits. Reference lines: Linf, the Linf of the largest group (no
#' fish in the model grows beyond it), L50 and, when length data are available,
#' the largest observed fish of each gear. Fish above the largest group cannot
#' be produced by the model (see `maxsd` in [fiticc()]).
#'
#' With `composition = TRUE` a side panel shows the unfished length composition
#' of each group and their sum, on the same length axis (Hordyk et al. 2015,
#' Fig. 1, right column), for constant M/k: the density of group \eqn{g} is
#' proportional to \eqn{(1 - L/L_{\infty,g})^{M/k - 1}}.
#'
#' Use it before the fit to check the growth inputs (Linf, CVL, maxsd) against
#' the largest fish, and after the fit to see the groups of the fitted model
#' (estimated Linf and CVL when `linf.sd` or `CVL.sd` are used).
#'
#' @param x A fitted FLicc object, an `FLStockLen` or an `FLPar` with `linf`,
#'   `k` and optionally `L50`, `M` or `Mk`.
#' @param CVL,maxsd,ngtg Growth variability and GTG grid. For a fit, the
#'   fitted values; otherwise defaults 0.1, 3 and 13 (FLicc's defaults).
#' @param maxage Largest age shown; default the age at which the central group
#'   reaches 0.95 Linf, \eqn{t_0 + \ln(20)/k}.
#' @param t0 Age at length zero (only shifts the age axis; FLicc does not need it).
#' @param lfd Optional `FLQuants` of length frequencies for the largest observed
#'   fish; for a fit, the fitted data.
#' @param composition Logical; add the unfished length composition (needs
#'   patchwork).
#' @return A `ggplot` (or a `patchwork` with `composition = TRUE`).
#' @references Hordyk, A., Ono, K., Valencia, S., Loneragan, N. and Prince, J.
#'   (2015). A novel length-based empirical estimation method of spawning
#'   potential ratio (SPR), and tests of its performance, for small-scale,
#'   data-poor fisheries. ICES Journal of Marine Science 72: 217-231.
#' @examples
#' \dontrun{
#' plot_growth(lhpar, CVL = 0.1, lfd = lfd)
#' plot_growth(fit, composition = TRUE)
#' plot_growth(lhpar, maxsd = 2, ngtg = 13, lfd = lfd)   # LBSPR's grid
#' }
#' @export
plot_growth <- function(x, CVL = NULL, maxsd = NULL, ngtg = NULL, maxage = NULL, t0 = 0,
                        lfd = NULL, composition = FALSE) {
  gp <- growth_pars_flicc(x, CVL, maxsd, ngtg)
  if (is.null(lfd)) lfd <- gp$lfd
  if (is.null(maxage)) maxage <- t0 + log(20) / gp$k
  age <- seq(0, maxage, length.out = 200)
  ng <- length(gp$z)
  d <- do.call(rbind, lapply(seq_len(ng), function(g)
    data.frame(group = g, age = age, L = pmax(gp$linf_g[g] * (1 - exp(-gp$k * (age - t0))), 0),
               w = gp$w[g])))
  d$shade <- d$w / max(d$w)
  central <- which.min(abs(gp$z))
  ref <- data.frame(what = c("Linf", "largest group", "L50"),
                    L = c(gp$linf, max(gp$linf_g), gp$L50),
                    col = c("black", "#D6604D", "darkgreen"), lty = c("dashed", "dashed", "dotted"))
  if (!is.null(lfd)) {
    ml <- max_len_flicc(lfd)
    ref <- rbind(ref, data.frame(what = paste("largest fish,", names(ml)), L = ml,
                                 col = "#2166AC", lty = "solid"))
  }
  ref <- ref[is.finite(ref$L), ]
  ymax <- max(c(ref$L, d$L), na.rm = TRUE) * 1.05
  p <- ggplot2::ggplot(d, ggplot2::aes(age, L, group = group)) +
    ggplot2::geom_line(ggplot2::aes(alpha = shade), colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_line(data = d[d$group == central, ], colour = "black", linewidth = 1.2) +
    ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = L), colour = ref$col,
                        linetype = ref$lty, inherit.aes = FALSE) +
    ggplot2::annotate("text", x = 0, y = ref$L, label = ref$what, hjust = 0, vjust = -0.3,
                      size = 3, colour = ref$col) +
    ggplot2::scale_alpha(range = c(0.2, 1), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::labs(x = "Age (years)", y = "Length",
                  title = sprintf("Growth-type groups: Linf %.1f, k %.3f, CVL %.2f, %d groups over +/- %.1f SD",
                                  gp$linf, gp$k, gp$CVL, ng, max(abs(gp$z)))) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                                         plot.title = ggplot2::element_text(size = 10))
  if (!composition) return(p)
  if (!is.finite(gp$Mk)) stop("composition = TRUE needs M/k (M and k, or Mk, in the parameters).")
  Lg <- seq(0, max(gp$linf_g), length.out = 300)
  dc <- do.call(rbind, lapply(seq_len(ng), function(g) {
    f <- ifelse(Lg < gp$linf_g[g], (1 - Lg / gp$linf_g[g])^(gp$Mk - 1) * gp$Mk / gp$linf_g[g], 0)
    data.frame(group = g, L = Lg, dens = gp$w[g] * f)
  }))
  tot <- stats::aggregate(dens ~ L, data = dc, FUN = sum)
  pc <- ggplot2::ggplot(dc, ggplot2::aes(L, dens, group = group)) +
    ggplot2::geom_line(colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_line(data = tot, ggplot2::aes(L, dens), colour = "black", linewidth = 1, inherit.aes = FALSE) +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = L), colour = ref$col,
                        linetype = ref$lty, inherit.aes = FALSE) +
    ggplot2::coord_flip(xlim = c(0, ymax)) +
    ggplot2::labs(x = NULL, y = "Unfished density", title = sprintf("M/k = %.2f", gp$Mk)) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                                         axis.text.x = ggplot2::element_blank(),
                                         plot.title = ggplot2::element_text(size = 10))
  if (!requireNamespace("patchwork", quietly = TRUE)) return(list(growth = p, composition = pc))
  patchwork::wrap_plots(p, pc, widths = c(3, 1))
}

#' Phase of the cohort modes and an estimate of t0
#'
#' Fits the length compositions as mixtures of age components whose mean
#' lengths follow von Bertalanffy growth, as in MULTIFAN (Fournier et al.
#' 1990): component \eqn{j = 0, 1, \dots} has mean
#' \eqn{\mu_j = L_\infty(1 - e^{-k(j + \phi)})} and standard deviation
#' \eqn{\sqrt{(\mathrm{CV}\,\mu_j)^2 + sd_0^2}}, with CV the realised spread of
#' Linf among the growth-type groups and \eqn{sd_0} an extra spread (spawning
#' over part of the year, measurement). The weights of the components are free
#' in each composition, so year-class strength and selectivity do not enter;
#' only the positions and widths of the modes have to fit. The phase \eqn{\phi}
#' (in years) is common to all compositions and is profiled over `phase`.
#'
#' A clear maximum of the profile means that the modes follow the growth curve;
#' a flat profile that they do not, or that the compositions show no modes
#' (continuous recruitment, data pooled over seasons). With a vector of `k`
#' values the profile is two-dimensional and checks the growth rate as well
#' (Linf fixed).
#'
#' **From phase to t0.** With the compositions sampled \eqn{\tau} years after
#' spawning (`sample - spawn`, as fractions of the year),
#' \eqn{t_0 = \tau - \phi - m}. The integer \eqn{m} is not identified by the
#' lengths alone (is the first mode age 0 or age 1?); the candidates
#' \eqn{m = -1, 0, 1} are returned with the implied lengths at age 1 and 2, to be
#' chosen with knowledge of the biology. In SS3 terms, the length at age 1 (or
#' at `Amin`) and CVL as the CV of the young fish.
#'
#' Partly selected small fish shift the first mode upwards; exclude them with
#' `Lmin`. Compositions pooled over months smear the modes; use compositions
#' from a short sampling period where possible.
#'
#' @param x A fitted FLicc object (fitted Linf, k, CVL and data), an
#'   `FLStockLen` or an `FLPar`.
#' @param lfd `FLQuants` of length frequencies (default: the data of the fit).
#' @param CVL,maxsd Growth variability (default: the fit's, or 0.1 and 3).
#' @param k Optional vector of k values for a joint profile (default: the
#'   input k).
#' @param phase Grid of phases in [0, 1).
#' @param sd0 Extra standard deviation of length at age; default the width of
#'   one length class.
#' @param years,gear Optional selection of years and gears.
#' @param Lmin Optional smallest length used.
#' @param ess Each composition is rescaled to this sample size for the
#'   likelihood (the profile interval is approximate).
#' @param spawn,sample Spawning and sampling time as fractions of the year.
#' @param maxage Largest age component (default: age at 0.99 Linf).
#' @return A list of class `"flicc_growth_phase"`: `profile` (k, phase,
#'   logLik, dLL), `best` (k, phase, logLik), `t0` (candidates with lengths at
#'   age 1 and 2), the inputs, and `contrast` (largest minus smallest logLik of
#'   the profile; near zero means no information on the phase).
#' @references Fournier, D.A., Sibert, J.R., Majkowski, J. and Hampton, J.
#'   (1990). MULTIFAN a likelihood-based method for estimating growth parameters
#'   and age composition from multiple length frequency data sets. Canadian
#'   Journal of Fisheries and Aquatic Sciences 47: 301-317.
#' @examples
#' \dontrun{
#' gp <- growth_phase_flicc(fit)
#' gp$best; gp$t0
#' plot_growth_phase(gp)
#' gp2 <- growth_phase_flicc(fit, k = seq(0.06, 0.12, 0.01))   # check k as well
#' }
#' @export
growth_phase_flicc <- function(x, lfd = NULL, CVL = NULL, maxsd = NULL, k = NULL,
                               phase = seq(0, 0.98, by = 0.02), sd0 = NULL,
                               years = NULL, gear = NULL, Lmin = NULL, ess = 200,
                               spawn = 0, sample = 0, maxage = NULL) {
  gp <- growth_pars_flicc(x, CVL, maxsd)
  if (is.null(lfd)) lfd <- gp$lfd
  if (is.null(lfd)) stop("lfd is needed (FLQuants of length frequencies).")
  gears <- if (is.null(gear)) names(lfd) else if (is.numeric(gear)) names(lfd)[gear] else gear
  ks <- if (is.null(k)) gp$k else k
  # compositions
  comps <- list()
  for (g in gears) {
    m <- lfd_matrix_flicc(lfd[[g]])
    L <- as.numeric(rownames(m)); step <- if (length(L) > 1) L[2] - L[1] else 1
    keep <- if (is.null(Lmin)) rep(TRUE, length(L)) else L >= Lmin
    yrs <- colnames(m); if (!is.null(years)) yrs <- intersect(yrs, as.character(years))
    for (y in yrs) {
      o <- m[keep, y]
      if (sum(o) <= 0) next
      comps[[length(comps) + 1]] <- list(gear = g, year = y, lo = L[keep], hi = L[keep] + step,
                                         o = o / sum(o) * (if (is.null(ess)) sum(o) else ess))
    }
  }
  if (!length(comps)) stop("No compositions with data.")
  step <- comps[[1]]$hi[1] - comps[[1]]$lo[1]
  if (is.null(sd0)) sd0 <- step
  ll_comp <- function(o, P) {
    keepc <- colSums(P) > 1e-10
    P <- P[, keepc, drop = FALSE]
    if (!ncol(P)) return(-Inf)
    P <- sweep(P, 2, colSums(P), "/")
    w <- rep(1 / ncol(P), ncol(P)); N <- sum(o); ll0 <- -Inf
    for (it in 1:300) {
      q <- as.numeric(P %*% w); q[q < 1e-300] <- 1e-300
      ll <- sum(o * log(q))
      if (ll - ll0 < 1e-6) break
      ll0 <- ll
      w <- w * as.numeric(crossprod(P, o / q)) / N
    }
    ll
  }
  prof <- expand.grid(k = ks, phase = phase)
  prof$logLik <- NA_real_
  for (i in seq_len(nrow(prof))) {
    kk <- prof$k[i]; ph <- prof$phase[i]
    amax <- if (is.null(maxage)) log(100) / kk else maxage
    mu <- gp$linf * (1 - exp(-kk * (seq(0, ceiling(amax)) + ph)))
    sdv <- sqrt((gp$cv_real * mu)^2 + sd0^2)
    prof$logLik[i] <- sum(vapply(comps, function(cc) {
      P <- sapply(seq_along(mu), function(j)
        stats::pnorm(cc$hi, mu[j], sdv[j]) - stats::pnorm(cc$lo, mu[j], sdv[j]))
      ll_comp(cc$o, matrix(P, nrow = length(cc$lo)))
    }, numeric(1)))
  }
  prof$dLL <- prof$logLik - max(prof$logLik)
  b <- prof[which.max(prof$logLik), ]
  tau <- sample - spawn
  m <- -1:1
  t0 <- tau - b$phase - m
  t0tab <- data.frame(m = m, t0 = round(t0, 3),
                      L_age1 = round(gp$linf * (1 - exp(-b$k * (1 - t0))), 2),
                      L_age2 = round(gp$linf * (1 - exp(-b$k * (2 - t0))), 2))
  out <- list(profile = prof, best = b, t0 = t0tab, contrast = diff(range(prof$logLik)),
              linf = gp$linf, k = b$k, cv_real = gp$cv_real, sd0 = sd0, spawn = spawn,
              sample = sample, ess = ess, gears = gears, lfd = lfd, CVL = gp$CVL)
  class(out) <- c("flicc_growth_phase", "list")
  out
}

#' Profile of the cohort-mode phase
#'
#' @param x Output of [growth_phase_flicc()].
#' @return A `ggplot`: the profile log-likelihood (relative to its maximum) over
#'   the phase, with the approximate 95% limit (dLL = -1.92); a heat map when
#'   several k were profiled.
#' @export
plot_growth_phase <- function(x) {
  pr <- x$profile
  if (length(unique(pr$k)) > 1) {
    return(ggplot2::ggplot(pr, ggplot2::aes(phase, k, fill = pmax(dLL, -50))) +
             ggplot2::geom_tile() +
             ggplot2::geom_point(data = x$best, ggplot2::aes(phase, k), inherit.aes = FALSE,
                                 shape = 4, size = 3) +
             ggplot2::scale_fill_viridis_c(name = "dLL") +
             ggplot2::labs(x = "Phase (years)", y = "k", title = "Profile over phase and k") +
             ggplot2::theme_bw())
  }
  ggplot2::ggplot(pr, ggplot2::aes(phase, dLL)) +
    ggplot2::geom_hline(yintercept = -1.92, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::geom_vline(xintercept = x$best$phase, colour = "#D6604D") +
    ggplot2::coord_cartesian(ylim = c(max(min(pr$dLL), -50), 0)) +
    ggplot2::labs(x = "Phase (years)", y = "Log-likelihood relative to the maximum",
                  title = sprintf("Phase of the cohort modes: %.2f (k = %.3f)", x$best$phase, x$best$k)) +
    ggplot2::theme_bw()
}

#' Cohort tracks over the length compositions
#'
#' Ridgeline plot of the length compositions by year (one panel per gear), with
#' von Bertalanffy cohort tracks at the phase of [growth_phase_flicc()] and,
#' optionally, bands of \eqn{\pm 2} standard deviations of length at age
#' (realised CVL of the growth-type groups and `sd0`). Tracks that follow the
#' modes from year to year support the growth parameters; modes that move
#' faster or slower than the tracks point to k or Linf. Each track is a year
#' class; their spacing is the annual growth increment
#' \eqn{(L_\infty - L)(1 - e^{-k})}, which does not depend on t0.
#'
#' @param x Output of [growth_phase_flicc()], or a fit / `FLStockLen` / `FLPar`
#'   (then [growth_phase_flicc()] is run with `...`).
#' @param phase Optional phase to draw instead of the estimate.
#' @param band Logical; draw the \eqn{\pm 2} SD bands.
#' @param scale Width of the largest composition, in years.
#' @param maxtracks Largest age of the tracks drawn (default: age at 0.95 Linf).
#' @param ... Passed to [growth_phase_flicc()] when `x` is not its output.
#' @return A `ggplot`.
#' @examples
#' \dontrun{
#' gp <- growth_phase_flicc(fit)
#' plot_growth_lfd(gp)
#' }
#' @export
plot_growth_lfd <- function(x, phase = NULL, band = TRUE, scale = 0.9, maxtracks = NULL, ...) {
  if (!inherits(x, "flicc_growth_phase")) x <- growth_phase_flicc(x, ...)
  ph <- if (is.null(phase)) x$best$phase else phase
  kk <- x$k; linf <- x$linf
  dd <- do.call(rbind, lapply(x$gears, function(g) {
    m <- lfd_matrix_flicc(x$lfd[[g]])
    L <- as.numeric(rownames(m)); step <- if (length(L) > 1) L[2] - L[1] else 1
    p <- sweep(m, 2, pmax(colSums(m), 1e-12), "/")
    do.call(rbind, lapply(colnames(m), function(y) {
      if (sum(m[, y]) <= 0) return(NULL)
      data.frame(gear = g, year = as.numeric(y), L = L + step / 2, p = p[, y])
    }))
  }))
  dd$gear <- factor(dd$gear, levels = x$gears)
  dd$x <- dd$year + scale * dd$p / max(dd$p)
  yrs <- sort(unique(dd$year))
  amax <- if (is.null(maxtracks)) log(20) / kk else maxtracks
  tr <- do.call(rbind, lapply(seq(min(yrs) - ceiling(amax), max(yrs)), function(b) {
    yy <- seq(min(yrs), max(yrs) + 1, by = 0.1)
    age <- yy - b + ph                 # year class b is in component 0 in year b
    ok <- age >= 0 & age <= amax
    if (!any(ok)) return(NULL)
    mu <- linf * (1 - exp(-kk * age[ok]))
    s <- sqrt((x$cv_real * mu)^2 + x$sd0^2)
    data.frame(cohort = b, year = yy[ok], mu = mu, lo = mu - 2 * s, hi = mu + 2 * s)
  }))
  p <- ggplot2::ggplot()
  if (band) p <- p + ggplot2::geom_ribbon(data = tr, ggplot2::aes(year, ymin = lo, ymax = hi, group = cohort),
                                          fill = "#D6604D", alpha = 0.06)
  p <- p +
    ggplot2::geom_polygon(data = do.call(rbind, lapply(split(dd, list(dd$gear, dd$year), drop = TRUE), function(z) {
      z <- z[order(z$L), ]
      rbind(z, transform(z[rev(seq_len(nrow(z))), ], x = year))
    })), ggplot2::aes(x, L, group = interaction(gear, year)), fill = "#6BAED6", colour = NA, alpha = 0.8) +
    ggplot2::geom_line(data = tr, ggplot2::aes(year, mu, group = cohort), colour = "#B2182B",
                       linetype = "dashed", linewidth = 0.4) +
    ggplot2::facet_wrap(~gear, ncol = 1) +
    ggplot2::coord_cartesian(ylim = c(0, max(dd$L[dd$p > 0], na.rm = TRUE) * 1.05)) +
    ggplot2::scale_x_continuous(breaks = yrs) +
    ggplot2::labs(x = "Year", y = "Length",
                  title = sprintf("Cohort tracks: Linf %.1f, k %.3f, phase %.2f", linf, kk, ph)) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  p
}
