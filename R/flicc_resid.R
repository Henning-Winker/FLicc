# Residual diagnostics and approximate posterior predictive checks for the
# length compositions

# Observed counts and predicted proportions as used in the likelihood: bins at
# or beyond the tail cut removed, bins from the plus group pooled, predicted
# proportions renormalised over the retained bins. One element per gear-year
# with data.
comp_data_flicc <- function(fit) {
  td <- fit$tmb_data
  rp <- fit$obj$report(fit$obj$env$last.par.best)
  obs <- td$obs
  plen <- rp$plen
  nlen <- dim(obs)[1]; nyear <- dim(obs)[2]; ngear <- dim(obs)[3]
  LLB <- as.numeric(td$LLB); Lmid <- as.numeric(td$Lmid)
  cut <- if (!is.null(td$cut_bin)) td$cut_bin else -1L
  plus <- if (!is.null(td$plus_bin)) td$plus_bin else -1L
  nkeep <- if (cut > 0 && cut < nlen) cut else nlen
  nb <- if (plus >= 0 && plus < nkeep) plus + 1 else nkeep
  grp <- pmin(seq_len(nkeep), nb)                      # bin -> retained/pooled bin
  Lb <- LLB[seq_len(nb)]
  Lm <- vapply(seq_len(nb), function(k) mean(Lmid[seq_len(nkeep)][grp == k]), numeric(1))
  out <- list()
  for (g in seq_len(ngear)) for (y in seq_len(nyear)) {
    o <- as.numeric(tapply(obs[seq_len(nkeep), y, g], grp, sum))
    p <- as.numeric(tapply(plen[seq_len(nkeep), y, g], grp, sum))
    n <- sum(o)
    if (!(n > 0) || !(sum(p) > 0)) next
    p <- pmax(p / sum(p), 1e-12); p <- p / sum(p)
    a0 <- if (!is.null(rp$dm_alpha0) && rp$dm_alpha0[y, g] > 0) rp$dm_alpha0[y, g] else n
    out[[length(out) + 1]] <- list(gear = td$gear_names[g], year = td$year_names[y],
                                   len = Lb, lmid = Lm, obs = o, p = p, n = n, a0 = a0)
  }
  phi <- tryCatch(exp(fit$obj$env$parList()$log_phi), error = function(e) NA_real_)
  list(comp = out, obs_model = fit$obs_model, phi = phi)
}

# integer year axis
year_axis_flicc <- function() ggplot2::scale_x_continuous(breaks = function(x) unique(round(pretty(x))))

# variance inflation of the Dirichlet-multinomial with precision alpha0
dm_factor_flicc <- function(n, a0 = n) (n + a0) / (1 + a0)

# beta-binomial cdf at x for N trials, shape a, b (vectorised over x)
pbetabinom_flicc <- function(x, N, a, b) {
  if (N <= 0) return(rep(1, length(x)))
  k <- 0:N
  lp <- lchoose(N, k) + lbeta(k + a, N - k + b) - lbeta(a, b)
  cdf <- cumsum(exp(lp - max(lp))); cdf <- cdf / cdf[length(cdf)]
  ifelse(x < 0, 0, ifelse(x >= N, 1, cdf[pmin(floor(x), N) + 1]))
}

#' Residuals of the length compositions
#'
#' Computes residuals of the observed length compositions by gear, year and
#' length bin, consistent with the likelihood of the fit (effective sample
#' size, observation model, tail cut and plus group).
#'
#' * `type = "pearson"`: \eqn{(o - n p) / \sqrt{V}}, with the variance of the
#'   observation model: \eqn{V = n p (1 - p)} (multinomial), the same multiplied
#'   by \eqn{(n + \alpha_0)/(1 + \alpha_0)} with the precision \eqn{\alpha_0}
#'   of the fit (Dirichlet-multinomial; \eqn{\alpha_0 = n} unless `dm_theta` is
#'   estimated), or \eqn{\mu + \mu^2/\phi} with \eqn{\mu = n p}
#'   (negative binomial). Pearson residuals of a composition are correlated
#'   (they sum to about zero) and are not normal for small expected counts;
#'   they are best read in bubble plots.
#' * `type = "osa"`: one-step-ahead (randomized quantile, PIT) residuals. The
#'   composition is decomposed along increasing length into conditional
#'   binomials (multinomial) or beta-binomials (Dirichlet-multinomial); each
#'   observed count gives \eqn{u = F(o-1) + v [F(o) - F(o-1)]}, \eqn{v \sim U(0,1)},
#'   and \eqn{r = \Phi^{-1}(u)}. Under a correct model the residuals are
#'   independent standard normal (Trijoulet et al. 2023, Fisheries Research
#'   257: 106487). The last bin of each composition is determined by the others
#'   and has no residual. Counts are rounded to integers for this step. For the
#'   negative binomial (independent bins) the PIT residual of each bin is used.
#'
#' @param fit A fitted FLicc object.
#' @param type `"pearson"` or `"osa"`.
#' @param seed Random seed for the randomization of the OSA residuals.
#' @return A data frame with gear, year, len (lower bound of the bin), lmid,
#'   obs, exp (expected count n p), n (sample size of the composition) and
#'   resid.
#' @examples
#' \dontrun{
#' r <- resid_flicc(fit, type = "osa")
#' qqnorm(r$resid); qqline(r$resid)
#' }
#' @export
resid_flicc <- function(fit, type = c("pearson", "osa"), seed = 1) {
  type <- match.arg(type)
  cd <- comp_data_flicc(fit)
  om <- cd$obs_model
  if (type == "osa") {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  rows <- lapply(cd$comp, function(z) {
    o <- z$obs; p <- z$p; n <- z$n; K <- length(o)
    if (type == "pearson") {
      V <- switch(om,
        nb = n * p + (n * p)^2 / cd$phi,
        dm = n * p * (1 - p) * dm_factor_flicc(n, z$a0),
        n * p * (1 - p))
      r <- (o - n * p) / sqrt(V)
    } else {
      oi <- round(o); N <- sum(oi); r <- rep(NA_real_, K)
      if (om == "nb") {
        mu <- n * p
        lo <- stats::pnbinom(oi - 1, mu = mu, size = cd$phi); hi <- stats::pnbinom(oi, mu = mu, size = cd$phi)
        r <- stats::qnorm(lo + stats::runif(K) * (hi - lo))
      } else {
        rem_n <- N; rem_p <- 1
        for (k in seq_len(K - 1)) {
          if (rem_n <= 0 || rem_p <= 1e-12) break
          q <- min(max(p[k] / rem_p, 1e-12), 1 - 1e-12)
          if (om == "dm") {
            a0 <- z$a0
            a <- a0 * p[k]; b <- a0 * max(rem_p - p[k], 1e-12)
            lo <- pbetabinom_flicc(oi[k] - 1, rem_n, a, b); hi <- pbetabinom_flicc(oi[k], rem_n, a, b)
          } else {
            lo <- stats::pbinom(oi[k] - 1, rem_n, q); hi <- stats::pbinom(oi[k], rem_n, q)
          }
          u <- lo + stats::runif(1) * (hi - lo)
          r[k] <- stats::qnorm(min(max(u, 1e-12), 1 - 1e-12))
          rem_n <- rem_n - oi[k]; rem_p <- rem_p - p[k]
        }
      }
    }
    data.frame(gear = z$gear, year = z$year, len = z$len, lmid = z$lmid, obs = o,
               exp = n * p, n = n, resid = r)
  })
  out <- do.call(rbind, rows)
  out$gear <- factor(out$gear, levels = fit$tmb_data$gear_names)
  attr(out, "type") <- type
  out
}

#' Bubble plot of length-composition residuals
#'
#' Residuals by length (y) and year (x), one panel per gear; bubble area is
#' proportional to the absolute residual, colour shows the sign. Optional cohort
#' tracks (von Bertalanffy growth curves of a few year classes, with
#' \eqn{t_0 = 0}) help to recognise cohorts growing through the compositions:
#' runs of same-signed residuals along the tracks indicate year-class strength,
#' which an equilibrium model cannot follow. Vertical bands at fixed lengths in
#' all years point instead to selectivity, misfit in the largest fish in all
#' years to Linf, CVL or doming.
#'
#' @param x Output of [resid_flicc()] or a fitted FLicc object.
#' @param type Residual type when `x` is a fit.
#' @param cohorts Logical; draw cohort tracks.
#' @param cohort_by Years between the cohort tracks.
#' @param fit The fit, needed for the cohort tracks when `x` is a residual
#'   data frame.
#' @param max_size Maximum bubble size.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' plot_resid_bubble(fit)
#' plot_resid_bubble(fit, type = "osa", cohorts = TRUE)
#' }
#' @export
plot_resid_bubble <- function(x, type = c("pearson", "osa"), cohorts = FALSE, cohort_by = 2,
                              fit = NULL, max_size = 6) {
  type <- match.arg(type)
  if (!is.data.frame(x)) { fit <- x; x <- resid_flicc(fit, type = type) }
  type <- if (!is.null(attr(x, "type"))) attr(x, "type") else type
  d <- x[is.finite(x$resid), ]
  d$year <- as.numeric(as.character(d$year))
  d$sign <- factor(ifelse(d$resid >= 0, "positive", "negative"), levels = c("positive", "negative"))
  p <- ggplot2::ggplot(d, ggplot2::aes(year, lmid))
  if (isTRUE(cohorts)) {
    if (is.null(fit)) stop("cohort tracks need the fit: plot_resid_bubble(r, fit = fit, cohorts = TRUE)")
    linf <- as.numeric(fit$report$lhpar["linf"]); k <- as.numeric(fit$report$lhpar["k"])
    yrs <- range(d$year)
    a_min <- -log(1 - min(min(d$lmid) / linf, 0.99)) / k      # age at the smallest length
    a_max <- -log(1 - min(0.9, max(d$lmid) / linf)) / k         # age at 0.9 Linf (or the largest length)
    births <- seq(floor(yrs[1] - a_max), ceiling(yrs[2] - a_min), by = cohort_by)
    tr <- do.call(rbind, lapply(births, function(b) {
      yy <- seq(yrs[1], yrs[2], by = 0.25)
      data.frame(year = yy, lmid = linf * (1 - exp(-k * (yy - b))), cohort = b)
    }))
    tr <- tr[tr$lmid > min(d$lmid) & tr$lmid < max(d$lmid) & tr$year >= tr$cohort, ]
    p <- p + ggplot2::geom_line(data = tr, ggplot2::aes(year, lmid, group = cohort),
                                colour = "grey70", linewidth = 0.3)
  }
  p + ggplot2::geom_point(ggplot2::aes(size = abs(resid), fill = sign), shape = 21,
                          colour = "grey30", stroke = 0.2, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(positive = "#2166AC", negative = "#D6604D"),
                               labels = c(positive = "observed > expected", negative = "observed < expected"),
                               name = NULL) +
    ggplot2::scale_size_area(max_size = max_size, name = "|residual|") +
    ggplot2::facet_wrap(~gear, scales = "free_y") + year_axis_flicc() +
    ggplot2::labs(x = "Year", y = "Length",
                  title = paste0(if (type == "osa") "OSA" else "Pearson", " residuals of the length compositions")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
}

#' QQ plot of OSA residuals by gear
#'
#' @param x Output of [resid_flicc()] with `type = "osa"`, or a fitted FLicc object.
#' @return A `ggplot` object.
#' @export
plot_resid_qq <- function(x) {
  if (!is.data.frame(x)) x <- resid_flicc(x, type = "osa")
  d <- x[is.finite(x$resid), ]
  ggplot2::ggplot(d, ggplot2::aes(sample = resid)) +
    ggplot2::stat_qq(size = 0.8, alpha = 0.5) + ggplot2::stat_qq_line(colour = "#D6604D") +
    ggplot2::facet_wrap(~gear) +
    ggplot2::labs(x = "Theoretical N(0,1) quantiles", y = "Residual quantiles",
                  title = "QQ plot of the residuals") +
    ggplot2::theme_bw()
}

# Wald-Wolfowitz runs test on the signs of x (one-sided: too few runs, i.e.
# positive autocorrelation); returns the p-value
runs_test_flicc <- function(x) {
  s <- sign(x[is.finite(x)]); s <- s[s != 0]
  n1 <- sum(s > 0); n2 <- sum(s < 0); n <- n1 + n2
  if (n1 == 0 || n2 == 0) return(NA_real_)
  runs <- 1 + sum(diff(s) != 0)
  mu <- 2 * n1 * n2 / n + 1
  v <- 2 * n1 * n2 * (2 * n1 * n2 - n) / (n^2 * (n - 1))
  if (v <= 0) return(NA_real_)
  stats::pnorm((runs - mu) / sqrt(v))
}

#' Mean-length residuals by year with runs test
#'
#' For each gear and year, compares the observed mean length of the
#' composition with its expected value under the fit (Francis 2011, TA1.8).
#' The standard error of the expected mean length follows from the
#' observation model: \eqn{\sqrt{\mathrm{Var}(L)/n}} for the multinomial,
#' multiplied by \eqn{\sqrt{(n+\alpha_0)/(1+\alpha_0)}} for the
#' Dirichlet-multinomial. The standardized residuals are tested for
#' randomness with a one-sided runs test (too few runs = systematic
#' patterns), and limits of three times their standard deviation are drawn,
#' as in the runs tests of `ss3diags` and JABBA: a green background means the
#' test is passed (p >= 0.05), red that it fails. Runs of years with longer or
#' shorter fish than expected are the clearest sign of strong year classes or
#' changes that the equilibrium model does not follow.
#'
#' @param fit A fitted FLicc object.
#' @return A list with `table` (gear, year, obs, exp, se, z) and `runs`
#'   (gear, p-value, pass), and a `ggplot` in `plot`.
#' @examples
#' \dontrun{
#' lm <- lmean_resid_flicc(fit)
#' lm$runs
#' lm$plot
#' }
#' @export
lmean_resid_flicc <- function(fit) {
  cd <- comp_data_flicc(fit)
  tab <- do.call(rbind, lapply(cd$comp, function(z) {
    L <- z$lmid; n <- z$n
    ex <- sum(L * z$p); vl <- sum(L^2 * z$p) - ex^2
    f <- switch(cd$obs_model, dm = dm_factor_flicc(n, z$a0),
                nb = 1 + n * sum(z$p^2 * (L - ex)^2) / max(cd$phi * vl, 1e-12), 1)
    se <- sqrt(vl / n * f)
    ob <- sum(L * z$obs) / n
    data.frame(gear = z$gear, year = as.numeric(as.character(z$year)), obs = ob, exp = ex, se = se,
               z = (ob - ex) / se)
  }))
  tab$gear <- factor(tab$gear, levels = fit$tmb_data$gear_names)
  runs <- do.call(rbind, lapply(split(tab, tab$gear), function(d) {
    if (!nrow(d)) return(NULL)
    pv <- runs_test_flicc(d$z)
    data.frame(gear = d$gear[1], p_value = round(pv, 3), pass = is.na(pv) || pv >= 0.05,
               sd_z = round(stats::sd(d$z), 2))
  }))
  lim <- do.call(rbind, lapply(split(tab, tab$gear), function(d) {
    if (!nrow(d)) return(NULL)
    s3 <- 3 * stats::sd(d$z)
    data.frame(gear = d$gear[1], lo = -s3, hi = s3, xmin = min(d$year) - 0.5, xmax = max(d$year) + 0.5)
  }))
  lim <- merge(lim, runs[, c("gear", "pass")], by = "gear")
  tab$out <- abs(tab$z) > 3 * stats::ave(tab$z, tab$gear, FUN = stats::sd)
  p <- ggplot2::ggplot(tab, ggplot2::aes(year, z)) +
    ggplot2::geom_rect(data = lim, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = lo, ymax = hi, fill = pass),
                       inherit.aes = FALSE, alpha = 0.25) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#4DAF4A", `FALSE` = "#E41A1C"),
                               labels = c(`TRUE` = "runs test passed", `FALSE` = "runs test failed"),
                               name = NULL, drop = FALSE) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_segment(ggplot2::aes(xend = year, yend = 0), colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = out), size = 2) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#E41A1C"), guide = "none") +
    ggplot2::facet_wrap(~gear) + year_axis_flicc() +
    ggplot2::labs(x = "Year", y = "Standardized residual of mean length",
                  title = "Mean-length residuals and runs test") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  list(table = tab, runs = runs, plot = p)
}

#' Approximate posterior predictive checks of the length compositions
#'
#' Simulates length compositions from the fitted observation model
#' (multinomial, Dirichlet-multinomial with the fitted precision, or negative
#' binomial) at the fitted expected compositions and sample sizes, and compares
#' summary statistics of the observed compositions with their predictive
#' distribution by gear and year. With `par_draws = TRUE` the parameters are
#' first drawn from their approximate normal distribution (estimates and
#' covariance from `sdreport`), so the predictive distribution also includes
#' estimation uncertainty; otherwise the estimates are used (plug-in).
#'
#' Statistics: `Lmean` (mean length), `L95` (95th percentile of length),
#' `Pmat` (proportion of fish at or above the length at 50% maturity) and
#' `Lmax` (lower bound of the largest bin with fish). If far more than 5% of
#' the observed values fall outside the central 95% predictive interval, the
#' observation model or the effective sample size understate the variability
#' of the data; runs of years outside the interval on the same side point to
#' systematic misfit (e.g. cohorts).
#'
#' @param fit A fitted FLicc object.
#' @param nsim Number of simulated data sets.
#' @param par_draws Logical; draw parameters from their normal approximation.
#' @param seed Random seed.
#' @return A list with `table` (gear, year, stat, obs, lo, med, hi, ppp =
#'   predictive p-value \eqn{P(sim \ge obs)}, outside), `coverage` (share of
#'   observed values outside the 95% interval by gear and statistic) and
#'   `plot`.
#' @examples
#' \dontrun{
#' pc <- ppc_flicc(fit, nsim = 500)
#' pc$coverage
#' pc$plot
#' }
#' @export
ppc_flicc <- function(fit, nsim = 500, par_draws = FALSE, seed = 1) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  cd <- comp_data_flicc(fit)
  om <- cd$obs_model
  lmat <- if ("L50" %in% dimnames(fit$report$lhpar)[[1]]) as.numeric(fit$report$lhpar["L50"]) else NA_real_
  stats_of <- function(cnt, L, lb) {
    n <- sum(cnt); if (n <= 0) return(c(Lmean = NA, L95 = NA, Pmat = NA, Lmax = NA))
    cp <- cumsum(cnt) / n
    c(Lmean = sum(L * cnt) / n, L95 = L[which(cp >= 0.95)[1]],
      Pmat = if (is.finite(lmat)) sum(cnt[lb >= lmat]) / n else NA,
      Lmax = lb[max(which(cnt > 0))])
  }
  # expected proportions per draw
  p_draws <- if (isTRUE(par_draws)) {
    if (is.null(fit$rep) || is.null(fit$rep$cov.fixed)) stop("par_draws = TRUE needs a fit with se = TRUE.")
    est <- fit$opt$par; V <- fit$rep$cov.fixed
    Lc <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(Lc)) stop("Covariance of the parameters is not positive definite.")
    lapply(seq_len(nsim), function(i) {
      th <- est + as.numeric(stats::rnorm(length(est)) %*% Lc)
      fit$obj$report(th)$plen
    })
  } else NULL
  cd_list <- cd$comp
  td <- fit$tmb_data
  nlen <- dim(td$obs)[1]
  cut <- if (!is.null(td$cut_bin)) td$cut_bin else -1L
  plus <- if (!is.null(td$plus_bin)) td$plus_bin else -1L
  nkeep <- if (cut > 0 && cut < nlen) cut else nlen
  nb <- if (plus >= 0 && plus < nkeep) plus + 1 else nkeep
  grp <- pmin(seq_len(nkeep), nb)
  rows <- lapply(cd_list, function(z) {
    gi <- match(z$gear, td$gear_names); yi <- match(z$year, td$year_names)
    n <- z$n; K <- length(z$p)
    sims <- vapply(seq_len(nsim), function(i) {
      p <- if (is.null(p_draws)) z$p else {
        pp <- as.numeric(tapply(p_draws[[i]][seq_len(nkeep), yi, gi], grp, sum)); pmax(pp / sum(pp), 1e-12)
      }
      cnt <- switch(om,
        nb = stats::rnbinom(K, mu = n * p, size = cd$phi),
        dm = { a <- stats::rgamma(K, shape = z$a0 * p); stats::rmultinom(1, round(n), a / sum(a))[, 1] },
        stats::rmultinom(1, round(n), p)[, 1])
      stats_of(cnt, z$lmid, z$len)
    }, numeric(4))
    ob <- stats_of(z$obs, z$lmid, z$len)
    data.frame(gear = z$gear, year = as.numeric(as.character(z$year)), stat = names(ob), obs = ob,
               lo = apply(sims, 1, stats::quantile, 0.025, na.rm = TRUE),
               med = apply(sims, 1, stats::median, na.rm = TRUE),
               hi = apply(sims, 1, stats::quantile, 0.975, na.rm = TRUE),
               ppp = rowMeans(sims >= ob, na.rm = TRUE), row.names = NULL)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[is.finite(tab$obs), ]
  tab$outside <- tab$obs < tab$lo | tab$obs > tab$hi
  tab$gear <- factor(tab$gear, levels = td$gear_names)
  tab$stat <- factor(tab$stat, levels = c("Lmean", "L95", "Pmat", "Lmax"))
  cov <- stats::aggregate(outside ~ gear + stat, data = tab, FUN = function(v) round(mean(v), 2))
  names(cov)[3] <- "share_outside_95"
  p <- ggplot2::ggplot(tab, ggplot2::aes(year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(y = med), colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(y = obs, colour = outside), size = 1.8) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#E41A1C"),
                                 labels = c(`FALSE` = "inside 95%", `TRUE` = "outside 95%"), name = NULL) +
    ggplot2::facet_grid(stat ~ gear, scales = "free_y") + year_axis_flicc() +
    ggplot2::labs(x = "Year", y = NULL,
                  title = paste0("Posterior predictive check (", if (isTRUE(par_draws)) "with parameter uncertainty" else "plug-in", ")")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  list(table = tab, coverage = cov, plot = p)
}

#' Effective sample sizes of the length compositions
#'
#' For each gear and year: the input sample size \eqn{n} (after `lfdess()`,
#' within the retained bins), the Dirichlet-multinomial precision
#' \eqn{\alpha_0} and the effective sample size it implies,
#' \eqn{n_{eff} = n(1 + \alpha_0)/(n + \alpha_0)}. With `dm_theta = "fixed"`
#' (\eqn{\alpha_0 = n}) \eqn{n_{eff} \approx n/2}; with `"saturating"` it
#' levels off at \eqn{\beta + 1}. For the multinomial \eqn{n_{eff} = n}.
#'
#' @param fit A fitted FLicc object.
#' @return A data frame with gear, year, n, alpha0 and neff.
#' @examples
#' \dontrun{
#' fit <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
#'               settings = list(obs_model = "dm", dm_theta = "saturating"))
#' neff_flicc(fit)
#' }
#' @export
neff_flicc <- function(fit) {
  cd <- comp_data_flicc(fit)
  out <- do.call(rbind, lapply(cd$comp, function(z) {
    a0 <- if (identical(cd$obs_model, "dm")) z$a0 else NA_real_
    data.frame(gear = z$gear, year = as.numeric(as.character(z$year)), n = z$n, alpha0 = a0,
               neff = if (identical(cd$obs_model, "dm")) z$n * (1 + a0) / (z$n + a0) else z$n)
  }))
  out$gear <- factor(out$gear, levels = fit$tmb_data$gear_names)
  out
}
