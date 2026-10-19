#' Combine several FLicc fits into an ensemble
#'
#' Integrates a list of fitted FLicc models (e.g. the dome levels of
#' [dome_sens_flicc()], alternative Linf or tail settings) into one estimate of
#' SPR by year with an interval that includes both the uncertainty within each
#' fit and the differences between fits. For each year, the ensemble is the
#' mixture of the members' normal distributions of log SPR (estimate and
#' standard error from `sdreport`, as in [spr_ci_flicc()]):
#' \deqn{p(\log SPR) = \sum_i w_i\, N(\mu_i, se_i^2),}
#' with total variance \eqn{\sum_i w_i se_i^2 + \sum_i w_i (\mu_i - \bar\mu)^2}
#' (within + between fits). The median and interval are quantiles of the
#' mixture, so a skewed or bimodal ensemble is represented as it is.
#'
#' Weights:
#' * `"equal"` (default): all fits count the same. Appropriate for structural
#'   assumptions the data cannot distinguish (e.g. slight doming).
#' * `"aic"`: Akaike weights \eqn{\exp(-\Delta AIC/2)}. Only valid when all
#'   fits use the same data and likelihood (same observations, observation
#'   model, tail cut, plus group and robust mixing); otherwise the function
#'   stops. Akaike weights tend to collapse onto one fit when the effective
#'   sample size is large.
#' * a numeric vector: user weights (normalised to sum 1).
#'
#' With steepness `s`, the quantiles are also given for spawning biomass
#' relative to the SPR target under Beverton-Holt recruitment,
#' \eqn{B/B_{tgt} = (4s\,SPR - (1-s)) / (4s\,SPR_{tgt} - (1-s))}, a monotone
#' transform of SPR. Apical F has no standard error in the fits; the ensemble
#' reports its weighted mean and the range over the members.
#'
#' @param fits A named list of fitted FLicc objects (from [fiticc()]), fitted
#'   to the same years. Fits without standard errors (`se = FALSE`) enter with
#'   their point estimate only.
#' @param weights `"equal"`, `"aic"` or a numeric vector, one per fit.
#' @param level Interval level.
#' @param s Optional Beverton-Holt steepness for B/Btgt.
#' @param spr_tgt SPR target (proportion) for B/Btgt; default from the first
#'   fit's `spr_ref` or 0.4.
#' @return A list of class `"flicc_ensemble"` with `table` (year, spr = mixture
#'   median, lwr, upr, mean_log, sd_within, sd_between, sd_total, Fap_mean,
#'   Fap_min, Fap_max and, with `s`, brel, brel_lwr, brel_upr), `members`
#'   (name, weight, AIC, dAIC, last-year SPR and its se; dAIC is `NA` when the
#'   fits' data or likelihood differ, as AIC is then not comparable) and `fits`.
#' @examples
#' \dontrun{
#' ds <- dome_sens_flicc(fit, gear = "Gillnet")
#' ens <- ensemble_flicc(ds$fits[-1])          # the dome levels, equal weights
#' ens$members
#' plot_ensemble_flicc(ens)
#' ensemble_flicc(list(Linf_lo = f1, Linf = f2, Linf_hi = f3), s = 0.75)
#' }
#' @export
ensemble_flicc <- function(fits, weights = "equal", level = 0.95, s = NULL, spr_tgt = NULL) {
  if (!is.list(fits) || length(fits) < 2) stop("fits must be a list of at least two fitted FLicc objects.")
  if (is.null(names(fits)) || any(names(fits) == "")) names(fits) <- paste0("fit", seq_along(fits))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  K <- length(fits)
  yrs <- lapply(fits, function(f) as.character(f$tmb_data$year_names))
  if (!all(vapply(yrs, identical, logical(1), yrs[[1]])))
    stop("All fits must cover the same years.")
  yrs <- yrs[[1]]

  # AIC and comparability
  aic <- vapply(fits, function(f) tryCatch(stats::AIC(LLflicc(f)), error = function(e) NA_real_), numeric(1))
  same_lik <- function() {
    key <- function(f) {
      td <- f$tmb_data
      list(obs = td$obs, om = f$obs_model, cut = td$cut_bin, plus = td$plus_bin,
           rob = if (is.null(td$rob_eps)) 0 else td$rob_eps, dm = td$dm_type)
    }
    k1 <- key(fits[[1]])
    all(vapply(fits[-1], function(f) isTRUE(all.equal(key(f), k1)), logical(1)))
  }
  w <- if (is.numeric(weights)) {
    if (length(weights) != K || any(weights < 0) || sum(weights) <= 0) stop("numeric weights: one non-negative value per fit.")
    weights / sum(weights)
  } else {
    weights <- match.arg(weights, c("equal", "aic"))
    if (weights == "equal") rep(1 / K, K) else {
      if (!same_lik())
        stop("AIC weights need the same data and likelihood in all fits (observations, obs_model, ",
             "tail cut, plus group, rob_eps, DM precision type). Use equal or user weights.")
      if (any(!is.finite(aic))) stop("AIC not available for all fits.")
      d <- aic - min(aic); ww <- exp(-d / 2); ww / sum(ww)
    }
  }
  names(w) <- names(fits)

  # log SPR estimates and standard errors by year and fit
  est <- lapply(fits, function(f) {
    ci <- spr_ci_flicc(f)
    if (is.null(ci)) {
      data.frame(year = yrs, mu = log(as.numeric(f$report$spr)), se = 0)
    } else data.frame(year = ci$year, mu = log(ci$spr), se = ci$se_log)
  })
  if (any(vapply(est, function(e) all(e$se == 0), logical(1))))
    warning("Some fits have no standard errors (se = FALSE); they enter with their point estimates only.", call. = FALSE)
  MU <- sapply(est, `[[`, "mu"); SE <- sapply(est, `[[`, "se")
  if (!is.matrix(MU)) { MU <- matrix(MU, nrow = 1); SE <- matrix(SE, nrow = 1) }

  # mixture quantile on the log scale
  qmix <- function(p, mu, se, w) {
    se <- pmax(se, 1e-8)
    f <- function(x) sum(w * stats::pnorm((x - mu) / se)) - p
    lo <- min(mu - 10 * se); hi <- max(mu + 10 * se)
    stats::uniroot(f, c(lo, hi), tol = 1e-10)$root
  }
  a <- (1 - level) / 2
  tab <- do.call(rbind, lapply(seq_along(yrs), function(j) {
    mu <- MU[j, ]; se <- SE[j, ]
    m <- sum(w * mu)
    data.frame(year = as.numeric(yrs[j]),
               spr = exp(qmix(0.5, mu, se, w)),
               lwr = exp(qmix(a, mu, se, w)), upr = exp(qmix(1 - a, mu, se, w)),
               mean_log = m,
               sd_within = sqrt(sum(w * se^2)), sd_between = sqrt(sum(w * (mu - m)^2)),
               sd_total = sqrt(sum(w * se^2) + sum(w * (mu - m)^2)))
  }))
  # apical F: weighted mean and range of the members
  Fap <- sapply(fits, function(f) as.numeric(f$report$Fap))
  if (!is.matrix(Fap)) Fap <- matrix(Fap, nrow = 1)
  tab$Fap_mean <- as.numeric(Fap %*% w)
  tab$Fap_min <- apply(Fap, 1, min); tab$Fap_max <- apply(Fap, 1, max)

  if (!is.null(s)) {
    check_steep_flicc(s, 40)
    if (is.null(spr_tgt)) {
      sr <- fits[[1]]$report$spr_ref
      spr_tgt <- if (!is.null(sr)) sr[1] / 100 else 0.4
    }
    br <- function(x) pmax(4 * s * x - (1 - s), 0) / (4 * s * spr_tgt - (1 - s))
    tab$brel <- br(tab$spr); tab$brel_lwr <- br(tab$lwr); tab$brel_upr <- br(tab$upr)
  }
  nl <- length(yrs)
  comparable <- same_lik()
  members <- data.frame(name = names(fits), weight = round(w, 3),
                        AIC = round(aic, 2),
                        dAIC = if (comparable) round(aic - min(aic, na.rm = TRUE), 2) else NA_real_,
                        SPR_last = round(exp(MU[nl, ]), 3), se_log_last = round(SE[nl, ], 3),
                        row.names = NULL)
  structure(list(table = tab, members = members, fits = fits, weights = w, level = level, s = s),
            class = "flicc_ensemble")
}

#' Plot an FLicc ensemble
#'
#' SPR by year: the ensemble median and interval (ribbon) with the member fits
#' as thin lines (line width proportional to the weight). With `what = "brel"`
#' (needs `s` in [ensemble_flicc()]) spawning biomass relative to the target.
#' With fewer than four years the plot shows error bars instead: each member's
#' estimate and interval (labelled with its weight) and the ensemble median and
#' interval, one panel per year.
#'
#' @param x Output of [ensemble_flicc()].
#' @param what `"spr"` or `"brel"`.
#' @param members Logical; draw the member fits.
#' @param ref Reference lines (proportions of SPR, or of the target for brel).
#' @return A `ggplot` object.
#' @export
plot_ensemble_flicc <- function(x, what = c("spr", "brel"), members = TRUE, ref = NULL) {
  what <- match.arg(what)
  tab <- x$table
  if (what == "brel") {
    if (is.null(tab$brel)) stop("what = \"brel\" needs s in ensemble_flicc().")
    tab$y <- tab$brel; tab$lo <- tab$brel_lwr; tab$hi <- tab$brel_upr
    if (is.null(ref)) ref <- c(1, 0.5)
    ylab <- "B / B target"
  } else {
    tab$y <- tab$spr; tab$lo <- tab$lwr; tab$hi <- tab$upr
    if (is.null(ref)) ref <- c(0.4, 0.2)
    ylab <- "SPR"
  }
  tg_of <- function(f) if (!is.null(f$report$spr_ref)) f$report$spr_ref[1] / 100 else 0.4
  to_brel <- function(v, f) pmax(4 * x$s * v - (1 - x$s), 0) / (4 * x$s * tg_of(f) - (1 - x$s))
  ttl <- paste0("Ensemble of ", length(x$fits), " fits (", round(100 * x$level), "% interval)")
  # few years: estimate and interval of each member and of the ensemble as error bars
  if (length(unique(tab$year)) < 4) {
    md <- do.call(rbind, lapply(names(x$fits), function(nm) {
      f <- x$fits[[nm]]
      ci <- tryCatch(spr_ci_flicc(f, level = x$level), error = function(e) NULL)
      v <- as.numeric(f$report$spr)
      lo <- if (is.null(ci)) v else ci$lwr; hi <- if (is.null(ci)) v else ci$upr
      if (what == "brel") { v <- to_brel(v, f); lo <- to_brel(lo, f); hi <- to_brel(hi, f) }
      data.frame(year = as.numeric(f$tmb_data$year_names), y = v, lo = lo, hi = hi,
                 member = nm, w = x$weights[[nm]], type = "member")
    }))
    if (!isTRUE(members)) md <- md[0, ]
    ed <- data.frame(year = tab$year, y = tab$y, lo = tab$lo, hi = tab$hi,
                     member = "Ensemble", w = 1, type = "ensemble")
    d <- rbind(md, ed)
    d$member <- factor(d$member, levels = c(names(x$fits), "Ensemble"))
    d$lab <- if (isTRUE(members)) sprintf("w = %.2f", d$w) else ""
    d$lab[d$type == "ensemble"] <- ""
    p <- ggplot2::ggplot(d, ggplot2::aes(member, y, colour = type)) +
      ggplot2::geom_hline(yintercept = ref, linetype = c("solid", rep("dashed", length(ref) - 1)),
                          colour = c("darkgreen", rep("orange", length(ref) - 1))) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.25) +
      ggplot2::geom_point(ggplot2::aes(size = type)) +
      ggplot2::geom_text(ggplot2::aes(y = hi, label = lab), vjust = -0.6, size = 3, colour = "grey30") +
      ggplot2::scale_colour_manual(values = c(member = "grey35", ensemble = "#08519C"), guide = "none") +
      ggplot2::scale_size_manual(values = c(member = 1.8, ensemble = 3), guide = "none") +
      ggplot2::labs(x = NULL, y = ylab, title = ttl) + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (length(unique(d$year)) > 1) p <- p + ggplot2::facet_wrap(~year)
    return(p)
  }
  p <- ggplot2::ggplot(tab, ggplot2::aes(year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#9ECAE1", alpha = 0.6)
  if (isTRUE(members)) {
    md <- do.call(rbind, lapply(names(x$fits), function(nm) {
      f <- x$fits[[nm]]; v <- as.numeric(f$report$spr)
      if (what == "brel") {
        v <- to_brel(v, f)
      }
      data.frame(year = as.numeric(f$tmb_data$year_names), y = v, member = nm, w = x$weights[[nm]])
    }))
    p <- p + ggplot2::geom_line(data = md, ggplot2::aes(year, y, group = member, linewidth = w),
                                colour = "grey35", alpha = 0.7) +
      ggplot2::scale_linewidth(range = c(0.2, 1.2), guide = "none")
  }
  p + ggplot2::geom_hline(yintercept = ref, linetype = c("solid", rep("dashed", length(ref) - 1)),
                          colour = c("darkgreen", rep("orange", length(ref) - 1))) +
    ggplot2::geom_line(ggplot2::aes(y = y), colour = "#08519C", linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(y = y), colour = "#08519C", size = 1.6) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(x = "Year", y = ylab,
                  title = ttl) +
    ggplot2::theme_bw()
}
