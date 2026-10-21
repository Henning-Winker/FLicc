#' Selectivity at length for given parameters (natural length scale)
#'
#' Evaluates the FLicc selectivity functions with parameters on the length
#' scale (as in `sel_fixed`):
#' `"logistic"` (`SL50`, `SL95`), `"dsnormal"` (`mode`, `lsd`, `rsd`),
#' `"normal"` (`mode`, `sd`), `"binormal"` (`mode1`, `sd1`, `mode2`, `sd2`,
#' `h`) and `"bidsnormal"` (`mode1`, `sd1`, `mode2`, `lsd2`, `rsd2`, `h`); the
#' two-mode functions are scaled to a maximum of 1 over `L`.
#'
#' @param L Lengths (e.g. bin mid-points).
#' @param sel_fun `"logistic"`, `"dsnormal"`, `"normal"`, `"binormal"` or `"bidsnormal"`.
#' @param pars Named list or vector of parameters (single values).
#' @return Numeric vector of selectivity (maximum 1 for domes, asymptote 1 for logistic).
#' @examples
#' sel_curve_flicc(seq(10, 60, 5), "dsnormal", list(mode = 30, lsd = 4, rsd = 15))
#' @export
sel_curve_flicc <- function(L, sel_fun = c("logistic", "dsnormal", "normal", "binormal", "bidsnormal"), pars) {
  sel_fun <- match.arg(sel_fun)
  p <- as.list(pars)
  need <- switch(sel_fun, logistic = c("SL50", "SL95"), dsnormal = c("mode", "lsd", "rsd"),
                 normal = c("mode", "sd"), binormal = c("mode1", "sd1", "mode2", "sd2", "h"),
                 bidsnormal = c("mode1", "sd1", "mode2", "lsd2", "rsd2", "h"))
  miss <- setdiff(need, names(p))
  if (length(miss)) stop("Missing selectivity parameter(s) for ", sel_fun, ": ", paste(miss, collapse = ", "))
  switch(sel_fun,
    logistic = 1 / (1 + exp(-log(19) * (L - p$SL50) / (p$SL95 - p$SL50))),
    dsnormal = { d <- L - p$mode; exp(-d^2 / (2 * ifelse(L < p$mode, p$lsd, p$rsd)^2)) },
    normal   = exp(-(L - p$mode)^2 / (2 * p$sd^2)),
    binormal = { s <- exp(-(L - p$mode1)^2 / (2 * p$sd1^2)) + p$h * exp(-(L - p$mode2)^2 / (2 * p$sd2^2))
                 s / max(s) },
    bidsnormal = { d2 <- L - p$mode2
                   s <- exp(-(L - p$mode1)^2 / (2 * p$sd1^2)) +
                     p$h * exp(-d2^2 / (2 * ifelse(L < p$mode2, p$lsd2, p$rsd2)^2))
                   s / max(s) })
}

#' Explore selectivity curves over a range of parameter values
#'
#' Plots the selectivity function for one or more values of its parameters,
#' e.g. a sequence of right-hand standard deviations of a dome (from 1000 =
#' flat, no dome, to moderate doming), optionally over the observed length
#' composition and with reference lengths (length at maturity, mean length).
#' Useful for choosing fixed values for `sel_fixed` and the scenarios of
#' [dome_sens_flicc()].
#'
#' @param sel_fun `"logistic"` (`SL50`, `SL95`), `"dsnormal"` (`mode`, `lsd`,
#'   `rsd`), `"normal"` (`mode`, `sd`), `"binormal"` (`mode1`, `sd1`, `mode2`,
#'   `sd2`, `h`) or `"bidsnormal"` (`mode1`, `sd1`, `mode2`, `lsd2`, `rsd2`,
#'   `h`); see [fiticc()].
#' @param pars Named list of parameters on the length scale; elements with
#'   several values are crossed (`expand.grid`) and each combination is drawn.
#' @param len Lengths at which to evaluate; default from `lfd` or `0:80`.
#' @param lfd Optional `FLQuant` (one gear) or `FLQuants` (pooled over gears)
#'   of length frequencies, drawn as grey bars (scaled to a maximum of 1,
#'   pooled over years).
#' @param ref Optional named numeric of reference lengths drawn as vertical
#'   lines, e.g. `c(L50 = 31, Lmean = 36, Linf = 55.7)`.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' plot_selfun("dsnormal", list(mode = 32, lsd = 4, rsd = c(1000, 35, 20, 12)),
#'             lfd = lfd$Gillnet, ref = c(L50 = 31.2, Linf = 55.7))
#' plot_selfun("logistic", list(SL50 = c(26, 30, 34), SL95 = 38))
#' }
#' # two size modes (juveniles at 10 cm, adults at 25 cm), width of the adult curve
#' plot_selfun("binormal", list(mode1 = 10, sd1 = 1, mode2 = 25, sd2 = c(4, 8), h = 0.3))
#' # second mode from flat (rsd2 = 1000) to domed, and its relative height
#' plot_selfun("bidsnormal", list(mode1 = 10, sd1 = 1, mode2 = 25, lsd2 = 5,
#'                                rsd2 = c(1000, 15, 8), h = c(0.3, 1)))
#' @export
plot_selfun <- function(sel_fun = c("dsnormal", "logistic", "normal", "binormal", "bidsnormal"), pars, len = NULL,
                        lfd = NULL, ref = NULL) {
  sel_fun <- match.arg(sel_fun)
  obs <- NULL
  if (!is.null(lfd)) {
    q <- if (inherits(lfd, "FLQuants")) Reduce(`+`, lapply(lfd, function(x) x)) else lfd
    L0 <- as.numeric(dimnames(q)$len)
    n <- apply(q@.Data, 1, sum, na.rm = TRUE)
    step <- if (length(L0) > 1) diff(L0)[1] else 1
    obs <- data.frame(L = L0 + step / 2, p = n / max(n))
    if (is.null(len)) len <- seq(min(L0), max(L0) + step, length.out = 300)
  }
  if (is.null(len)) len <- seq(0, 80, length.out = 300)
  grid <- expand.grid(lapply(pars, unlist), KEEP.OUT.ATTRS = FALSE)
  vary <- names(pars)[vapply(pars, function(x) length(unlist(x)) > 1, logical(1))]
  lab <- if (length(vary)) apply(grid[, vary, drop = FALSE], 1, function(r)
    paste(paste0(vary, " = ", ifelse(r >= 999, paste0(r, " (no dome)"), signif(r, 3))), collapse = ", ")) else "selectivity"
  d <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    data.frame(L = len, sel = sel_curve_flicc(len, sel_fun, grid[i, , drop = FALSE]),
               scenario = lab[i])))
  d$scenario <- factor(d$scenario, levels = unique(lab))
  p <- ggplot2::ggplot()
  if (!is.null(obs)) p <- p + ggplot2::geom_col(data = obs, ggplot2::aes(L, p), fill = "grey85", width = 0.9)
  p <- p + ggplot2::geom_line(data = d, ggplot2::aes(L, sel, colour = scenario), linewidth = 0.9)
  if (!is.null(ref)) {
    rd <- data.frame(L = as.numeric(ref), name = names(ref))
    p <- p + ggplot2::geom_vline(data = rd, ggplot2::aes(xintercept = L), linetype = 2, colour = "grey40") +
      ggplot2::geom_text(data = rd, ggplot2::aes(x = L, y = 1.04, label = name), size = 3, colour = "grey30")
  }
  p + ggplot2::labs(x = "Length", y = "Selectivity", colour = NULL,
                    title = paste0(sel_fun, " selectivity")) +
    ggplot2::coord_cartesian(ylim = c(0, 1.06)) + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
}

#' Dome-shaped selectivity sensitivity
#'
#' Refits a fitted model with one gear's selectivity changed to a
#' double-sided normal whose right-hand standard deviation is fixed at a
#' sequence of values, from very large (`rsd = 1000`: flat right limb, no
#' doming) to mild and moderate doming; the mode and left-hand standard
#' deviation are estimated. All other settings, data and gears are those of
#' `fit`. Recommended as a standard sensitivity when the gear is assumed
#' logistic: an unmodelled slight dome reduces the number of large fish in the
#' catch and is otherwise read as higher F.
#'
#' @param fit A fitted FLicc object (from [fiticc()], with `FLRreport = TRUE`).
#' @param gear Gear name (or index) whose selectivity is varied.
#' @param rsd Right-hand standard deviations (length units). Default
#'   `c(1000, 0.85, 0.7, 0.55, 0.4) * Linf` (flat to
#'   moderate doming in even steps).
#' @param fixed Optional further fixed parameters for this gear, e.g.
#'   `list(lsd = 4)`.
#' @param se Passed to [fiticc()]; `TRUE` gives intervals for [plot_spr()].
#' @param ... Further arguments passed to [fiticc()].
#' @return A list with `fits` (named list: `base` = `fit`, then one fit per
#'   `rsd`), `table` (data frame: scenario, rsd, convergence, log-likelihood,
#'   AIC, last-year SPR, Fspr, apical F in the last year), `gear`, `rsd` and
#'   `errors` (messages of refits that failed).
#' @examples
#' \dontrun{
#' ds <- dome_sens_flicc(fit, gear = "Gillnet")
#' ds$table
#' plot_dome_sens(ds)
#' }
#' @export
dome_sens_flicc <- function(fit, gear, rsd = NULL, fixed = NULL, se = TRUE, ...) {
  gears <- fit$tmb_data$gear_names
  g <- if (is.numeric(gear)) gear else match(gear, gears)
  if (is.na(g)) stop("Unknown gear: ", gear, ". Gears: ", paste(gears, collapse = ", "))
  if (is.null(fit$report$obslen)) stop("fit needs FLRreport = TRUE (the default of fiticc()).")
  linf <- as.numeric(fit$stklen@lhpar["linf"])
  if (is.null(rsd)) rsd <- c(1000, round(c(0.85, 0.7, 0.55, 0.4) * linf, 1))
  sel_fun <- c("logistic", "dsnormal", "normal", "binormal", "bidsnormal")[fit$tmb_data$sel_type]
  sel_fun[g] <- "dsnormal"
  cbg <- fit$report$catch_by_gear
  fits <- list(base = fit)
  errs <- character()
  for (r in rsd) {
    sf <- list(c(list(rsd = r), fixed)); names(sf) <- gears[g]
    nm <- paste0("rsd", r)
    fits[[nm]] <- tryCatch(
      fiticc(fit$report$obslen, fit$stklen, sel_fun = sel_fun, catch_by_gear = cbg,
             settings = fit$settings, sel_fixed = sf, se = se, ...),
      error = function(e) {
        errs[[nm]] <<- conditionMessage(e)
        warning("rsd = ", r, ": ", conditionMessage(e), call. = FALSE); NULL })
  }
  fits <- fits[!vapply(fits, is.null, logical(1))]
  tab <- do.call(rbind, lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    ll <- LLflicc(f)
    data.frame(scenario = nm, rsd = if (nm == "base") NA else as.numeric(sub("rsd", "", nm)),
               conv = f$opt$convergence, logLik = round(as.numeric(ll), 2), AIC = round(stats::AIC(ll), 2),
               SPR_last = round(tail(c(f$report$spr), 1), 3), Fspr = round(as.numeric(fspr_flicc(f)), 4),
               Fap_last = round(tail(c(f$report$Fap), 1), 4))
  }))
  tab$dAIC <- round(tab$AIC - min(tab$AIC), 2)
  list(fits = fits, table = tab, gear = gears[g], rsd = rsd, errors = errs)
}

#' Plot a dome-shaped selectivity sensitivity
#'
#' @param x Output of [dome_sens_flicc()].
#' @param what `"sel"`: the varied gear's fitted selectivity per scenario;
#'   `"spr"`: SPR trajectories (via [plot_spr()]).
#' @param ... Passed to [plot_spr()] (e.g. `ci = TRUE`).
#' @return A `ggplot` object.
#' @export
plot_dome_sens <- function(x, what = c("sel", "spr"), ...) {
  what <- match.arg(what)
  if (what == "spr") return(plot_spr(x$fits, ...))
  d <- do.call(rbind, lapply(names(x$fits), function(nm) {
    s <- x$fits[[nm]]$report$sel_gear[[x$gear]]
    data.frame(L = as.numeric(dimnames(s)$len), sel = as.numeric(s), scenario = nm)
  }))
  d$scenario <- factor(d$scenario, levels = names(x$fits))
  ggplot2::ggplot(d, ggplot2::aes(L, sel, colour = scenario)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(x = "Length", y = "Selectivity", colour = NULL,
                  title = paste0(x$gear, ": fitted selectivity by scenario")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

#' Compare estimated selectivity across FLicc fits
#'
#' Plots the fitted selectivity-at-length curves of several FLicc fits, one
#' panel per gear and one colour per fit. Useful for sensitivity runs
#' (e.g. dome levels, tail cuts, Linf or CVL settings, data subsets).
#'
#' @param fits A named list of fitted FLicc objects (from [fiticc()]). Unnamed
#'   fits are labelled `fit1`, `fit2`, ...
#' @param gear Optional character vector of gears to show (default all).
#' @param colours Optional colours, one per fit (named by fit or in fit order).
#' @param linewidth Line width.
#' @param ncol Number of panel columns (default: one row per gear pair).
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' plot_compare_sel(list(base = fit, cut47 = fit_cut, dome = fit_dome))
#' plot_compare_sel(ds$fits, gear = "Gillnet")   # dome sensitivity runs
#' }
#' @export
plot_compare_sel <- function(fits, gear = NULL, colours = NULL, linewidth = 0.8, ncol = NULL) {
  if (inherits(fits, "flicc_tmb_fit")) fits <- list(fit = fits)
  if (!is.list(fits) || length(fits) < 1) stop("fits must be a list of fitted FLicc objects.")
  if (is.null(names(fits)) || any(names(fits) == "")) names(fits) <- paste0("fit", seq_along(fits))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  d <- do.call(rbind, lapply(names(fits), function(nm) {
    s <- fits[[nm]]$report$sel_gear
    if (is.null(s)) stop("Fit '", nm, "' has no report$sel_gear (needs FLRreport = TRUE).")
    df <- as.data.frame(s)
    if ("year" %in% names(df) && length(unique(df$year)) > 1)
      df <- df[df$year == max(df$year), , drop = FALSE]          # last year
    data.frame(L = as.numeric(as.character(df$len)), sel = as.numeric(df$data),
               gear = as.character(df$qname), fit = nm)
  }))
  if (!is.null(gear)) {
    miss <- setdiff(gear, unique(d$gear))
    if (length(miss)) stop("Unknown gear(s): ", paste(miss, collapse = ", "))
    d <- d[d$gear %in% gear, , drop = FALSE]
  }
  d$fit <- factor(d$fit, levels = names(fits))
  d$gear <- factor(d$gear, levels = unique(d$gear))
  p <- ggplot2::ggplot(d, ggplot2::aes(L, sel, colour = fit)) +
    ggplot2::geom_line(linewidth = linewidth) +
    ggplot2::facet_wrap(~gear, ncol = ncol) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02)) +
    ggplot2::labs(x = "Length", y = "Selectivity", colour = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_rect(fill = "grey92"))
  if (!is.null(colours)) {
    if (is.null(names(colours))) colours <- stats::setNames(rep(colours, length.out = length(fits)), names(fits))
    p <- p + ggplot2::scale_colour_manual(values = colours)
  } else if (length(fits) <= 8) {
    p <- p + ggplot2::scale_colour_manual(values = stats::setNames(
      c("#000000", "#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999")[seq_along(fits)],
      names(fits)))
  }
  p
}
