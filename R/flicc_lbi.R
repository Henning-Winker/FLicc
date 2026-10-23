# Empirical length-based indicators (no model fit needed)

# weighted quantile of lengths (class mid-points, counts), linear within classes
lbi_quantile <- function(lo, hi, n, p) {
  cs <- cumsum(n) / sum(n)
  vapply(p, function(pp) {
    i <- which(cs >= pp)[1]
    c0 <- if (i > 1) cs[i - 1] else 0
    lo[i] + (hi[i] - lo[i]) * (pp - c0) / max(cs[i] - c0, 1e-12)
  }, numeric(1))
}

# indicators of one composition
lbi_one <- function(lo, hi, n, linf, L50, Lopt, k, wt) {
  mid <- (lo + hi) / 2
  N <- sum(n)
  if (!(N > 0)) return(NULL)
  # Lc: first class with at least half the modal abundance (WKLIFE)
  ic <- which(n >= 0.5 * max(n))[1]
  Lc <- mid[ic]
  sel <- seq_along(n) >= ic
  Lmean <- sum(mid[sel] * n[sel]) / sum(n[sel])
  Lbar <- sum(mid * n) / N
  q <- lbi_quantile(lo, hi, n, c(0.25, 0.95))
  # mean length of the largest 5% of fish
  ord <- rev(seq_along(n)); cum <- 0; s <- 0; target <- 0.05 * N
  for (i in ord) {
    if (!(n[i] > 0)) next
    take <- min(n[i], target - cum)
    if (take <= 0) break
    s <- s + take * mid[i]; cum <- cum + take
  }
  Lmax5 <- s / cum
  LmaxY <- mid[which.max(n * wt)]
  LFM <- 0.75 * Lc + 0.25 * linf
  Pmat <- sum(n[mid >= L50]) / N
  Popt <- sum(n[mid >= 0.9 * Lopt & mid <= 1.1 * Lopt]) / N
  Pmega <- sum(n[mid > 1.1 * Lopt]) / N
  # Beverton-Holt total mortality from the mean length above the lower bound of the Lc class
  Zk <- (linf - Lmean) / (Lmean - lo[ic])
  data.frame(n = N, Lc = Lc, L25 = q[1], Lbar = Lbar, Lmean = Lmean, LmaxY = LmaxY,
             L95 = q[2], Lmax5 = Lmax5, LFM = LFM, Lopt = Lopt, L50 = L50, Linf = linf,
             Pmat = Pmat, Popt = Popt, Pmega = Pmega, Pobj = Pmat + Popt + Pmega,
             Z_k = Zk, Z = if (is.finite(k)) k * Zk else NA_real_)
}

# Cope and Punt (2009) decision tree: selectivity pattern from Pobj, then
# whether spawning biomass is at or above the target (0.4 B0) from Pmat or Popt
cope_punt_flicc <- function(Pmat, Popt, Pmega, Lmat_Lopt, tol = 0.02) {
  Pobj <- Pmat + Popt + Pmega
  # Pmat threshold interpolated between Lmat = 0.75 Lopt and Lmat = 0.9 Lopt
  w <- rep_len(pmin(pmax((Lmat_Lopt - 0.75) / 0.15, 0), 1), length(Pobj))
  out <- data.frame(Pobj = Pobj, cp_selectivity = NA_character_, cp_test = NA_character_,
                    cp_status = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_along(Pobj)) {
    if (Pobj[i] < 1) {
      if (Popt[i] + Pmega[i] <= tol) {
        sel <- "small, immature fish"; thr <- 0.25; v <- Pmat[i]; nm <- "Pmat"
      } else {
        sel <- "small and optimal-size fish"; thr <- 0.4 + (0.25 - 0.4) * w[i]; v <- Pmat[i]; nm <- "Pmat"
      }
      ok <- v > thr
    } else if (Pobj[i] < 2 - tol) {
      sel <- "maturity ogive"; thr <- 0.95 + (0.90 - 0.95) * w[i]; v <- Pmat[i]; nm <- "Pmat"
      ok <- v > thr
    } else if (Popt[i] < 1 - tol) {
      sel <- "optimal-size and larger fish"; thr <- 0.65; v <- Popt[i]; nm <- "Popt"
      ok <- v < thr
    } else {
      sel <- "optimal-size fish only"; nm <- NA; ok <- NA
    }
    out$cp_selectivity[i] <- sel
    out$cp_test[i] <- if (is.na(nm)) "no status information" else
      sprintf("%s %s %.2f", nm, if (nm == "Popt") "<" else ">", thr)
    out$cp_status[i] <- if (is.na(ok)) "no information" else if (ok) "SB >= target" else "SB < target"
  }
  out
}

# indicator ratios, reference levels and properties (ICES WKLIFE; Kell et al. 2022)
lbi_table_def <- function() {
  data.frame(
    indicator = c("Lc/L50", "L25/L50", "Lmax5/Linf", "L95/Linf", "Pmega", "LmaxY/Lopt", "Lmean/Lopt", "Lmean/LFM"),
    num = c("Lc", "L25", "Lmax5", "L95", "Pmega", "LmaxY", "Lmean", "Lmean"),
    den = c("L50", "L50", "Linf", "Linf", NA, "Lopt", "Lopt", "LFM"),
    ref = c(1, 1, 0.8, 0.8, 0.3, 0.9, 0.9, 1),
    lim = c(0.75, 0.75, 0.7, 0.7, 0.15, 0.75, 0.75, 0.9),
    property = c("Conservation (immatures)", "Conservation (immatures)",
                 "Conservation (large fish)", "Conservation (large fish)", "Conservation (large fish)",
                 "Optimal yield", "Optimal yield", "MSY proxy"),
    stringsAsFactors = FALSE)
}

#' Empirical length-based indicators by gear and year
#'
#' Computes the length-based indicators (LBIs) used in data-limited assessments
#' directly from the length frequencies of each gear and year, without fitting a
#' model: the ICES WKLIFE indicators (ICES 2015, 2018) and those of Froese
#' (2004), with the reference levels used by ICES and evaluated by Kell et al.
#' (2022), and a Beverton-Holt estimate of total mortality from the mean length.
#'
#' Indicators (lengths are class mid-points):
#' * `Lc`: first length class with at least half the modal abundance;
#'   `L25`, `L95`: 25th and 95th percentiles; `Lbar`: mean length;
#'   `Lmean`: mean length of fish from `Lc`; `LmaxY`: class with the largest catch
#'   biomass; `Lmax5`: mean length of the largest 5% of the fish.
#' * `Pmat`: proportion mature (\eqn{\ge} L50); `Popt`: proportion within
#'   \eqn{\pm}10% of Lopt; `Pmega`: proportion above 1.1 Lopt; `Pobj = Pmat +
#'   Popt + Pmega` (Cope and Punt 2009).
#' * Cope and Punt (2009) decision tree (`cp_selectivity`, `cp_test`,
#'   `cp_status`): `Pobj` identifies the selectivity pattern (< 1: small and
#'   immature fish; 1-2: selection close to the maturity ogive; 2: only mature fish
#'   of optimal size and larger), and `Pmat` or `Popt` then indicates whether
#'   spawning biomass is at or above the target of 0.4 B0. The `Pmat` thresholds
#'   depend on Lmat/Lopt (0.75 and 0.9 in Cope and Punt; interpolated between,
#'   held constant outside). A fishery taking only fish of optimal size gives no
#'   status information. In the plots, the row "Pobj (Cope-Punt)" shows Pobj,
#'   green if SB >= target, red if not, grey if no information.
#' * Reference lengths: `LFM = 0.75 Lc + 0.25 Linf` (length at F = M), `Lopt`
#'   (length of maximum cohort biomass: `2/3 Linf`, or `Linf 3 / (3 + M/k)` with
#'   `Lopt = "mk"`).
#' * `Z_k`, `Z`: Beverton-Holt mortality, \eqn{Z/k = (L_\infty - \bar L)/(\bar L - L')},
#'   with \eqn{\bar L} the mean length above \eqn{L'}, the lower bound of the Lc class
#'   (assumes equilibrium and knife-edge selection at \eqn{L'}).
#'
#' Ratios, reference levels and limits (`$ratios`): `status` is green at or
#' above the reference level (ICES WKLIFE), orange between the limit and the
#' reference level, red below the limit; `pass` = green.
#'
#' | Ratio | Reference (green) | Limit (red below) | Property |
#' |---|---|---|---|
#' | Lc/L50, L25/L50 | 1 | 0.75 | conservation of immatures |
#' | Lmax5/Linf, L95/Linf | 0.8 | 0.7 | conservation of large fish |
#' | Pmega | 0.3 | 0.15 | conservation of large fish |
#' | LmaxY/Lopt, Lmean/Lopt | 0.9 | 0.75 | optimal yield |
#' | Lmean/LFM | 1 | 0.9 | MSY proxy |
#'
#' The reference levels are those of ICES; the limits are pragmatic choices
#' (about 10-25% below the reference, half the reference for Pmega) that mark
#' a clear departure rather than a marginal one, and can be changed with
#' `limits`. For Lmean/LFM, 0.9 corresponds to a mean length about 10% below
#' the length expected at F = M. The Cope-Punt row has only two colours.
#'
#' The indicators of the upper tail (Lmax5, L95, Pmega) identify overfishing
#' best; those of the left limb (Lc, L25) are sensitive to recruitment
#' variability (Kell et al. 2022). Indicators follow trends better than they
#' classify status, and the reference levels depend on the assumed life history.
#' With `ref_years`, each indicator is also given relative to its mean over a
#' reference period when the stock is believed to have been in a good state
#' (`$relative`), which avoids proxies based on uncertain parameters (Kell and
#' Sharma 2025).
#'
#' @param x `FLQuants` of length frequencies (one element per gear; lengths as
#'   lower class bounds), or a fitted FLicc object (its data and life history,
#'   with the fitted Linf).
#' @param lhpar `FLPar` with `linf` and `L50` (and `k`, `M` for `Z` and `Lopt =
#'   "mk"`; `a`, `b` for `LmaxY`). Taken from the fit if `x` is a fit.
#' @param gear,years Optional selection.
#' @param pool Logical; add the gears pooled into one composition as gear
#'   `"pooled"`. The length frequency of each gear and year is converted to
#'   proportions, raised to catch numbers with the gear's relative catch
#'   (`catch_by_gear`; in weight divided by the mean weight of a sampled fish)
#'   and summed over gears ([combine_lfd()]), so each gear contributes in
#'   proportion to its catch rather than to its sampling effort. Default `TRUE`
#'   when `catch_by_gear` is given.
#' @param catch_by_gear Relative catch by gear for pooling: a numeric vector
#'   (named by gear, or in the order of the data) or `FLQuants` by gear with the
#'   catch by year; only relative values matter. For a fit, the fit's catch by
#'   gear is used when not given. Without it, gears are pooled with equal catch
#'   shares (with a message).
#' @param catch_units Units of `catch_by_gear`: `"weight"` (default; needs `a`,
#'   `b` in `lhpar`) or `"numbers"`. For a fit, the fit's setting.
#' @param ref_years Optional reference years for indicators relative to a
#'   historical reference period.
#' @param Lopt `"wklife"` (2/3 Linf, default) or `"mk"` (Linf 3/(3 + M/k)).
#' @param limits Optional named numeric vector replacing limits (orange/red
#'   boundary) by ratio name, e.g. `c("Lmean/LFM" = 0.85, Pmega = 0.1)`.
#' @return A list of class `"flicc_lbi"`: `values` (one row per gear and year),
#'   `ratios` (long format: gear, year, indicator, value, ref, lim, status, pass,
#'   property),
#'   `relative` (with `ref_years`), and the inputs.
#' @references
#' Cope, J.M. and Punt, A.E. (2009). Length-based reference points for
#' data-limited situations. Marine and Coastal Fisheries 1: 169-186.
#'
#' Froese, R. (2004). Keep it simple: three indicators to deal with
#' overfishing. Fish and Fisheries 5: 86-91.
#'
#' ICES (2015). Report of the Fifth Workshop on the Development of Quantitative
#' Assessment Methodologies based on Life-history Traits (WKLIFE V).
#'
#' Kell, L.T., Minto, C. and Gerritsen, H.D. (2022). Evaluation of the skill of
#' length-based indicators to identify stock status and trends. ICES Journal of
#' Marine Science 79: 1202-1216.
#'
#' Kell, L.T. and Sharma, R. (2025). An evaluation of the robustness of
#' length-based stock assessment approaches for sustainable fisheries
#' management in data and capacity limited situations. Sustainability 17: 4791.
#' @examples
#' \dontrun{
#' lbi <- lbi_flicc(lfd, lhpar)
#' lbi$values
#' plot_lbi(lbi)                    # traffic light
#' plot_lbi(lbi, type = "trend")
#' lbi_flicc(lfd, lhpar, catch_by_gear = c(Trawl = 0.7, Gillnet = 0.3), ref_years = 2013:2015)
#' }
#' @export
lbi_flicc <- function(x, lhpar = NULL, gear = NULL, years = NULL, pool = !is.null(catch_by_gear),
                      catch_by_gear = NULL, catch_units = NULL,
                      ref_years = NULL, Lopt = c("wklife", "mk"), limits = NULL) {
  Lopt <- match.arg(Lopt)
  is_fit <- is.list(x) && !is.null(x$tmb_data) && !is.null(x$report)
  if (is_fit) {
    lfd <- x$report$obslen
    if (is.null(lhpar)) lhpar <- x$report$lhpar
    if (pool && is.null(catch_by_gear)) catch_by_gear <- x$report$catch_by_gear
    if (is.null(catch_units)) catch_units <- if (is.null(x$settings$catch_units)) "weight" else x$settings$catch_units
  } else lfd <- x
  if (is.null(catch_units)) catch_units <- "weight"
  catch_units <- match.arg(catch_units, c("weight", "numbers"))
  if (!inherits(lfd, "FLQuants")) stop("x must be FLQuants of length frequencies or a fitted FLicc object.")
  if (is.null(lhpar)) stop("lhpar (FLPar with linf and L50) is needed.")
  getp <- function(nm) if (nm %in% dimnames(lhpar)[[1]]) as.numeric(lhpar[nm]) else NA_real_
  linf <- getp("linf"); L50 <- getp("L50"); k <- getp("k")
  Mk <- getp("Mk"); if (!is.finite(Mk)) Mk <- getp("M") / k
  if (!is.finite(linf)) stop("lhpar needs linf.")
  if (!is.finite(L50)) { L50 <- 2 / 3 * linf; message("lbi_flicc: no L50 in lhpar; using 2/3 Linf.") }
  Lo <- if (Lopt == "mk") {
    if (!is.finite(Mk)) stop('Lopt = "mk" needs M and k (or Mk) in lhpar.')
    linf * 3 / (3 + Mk)
  } else 2 / 3 * linf
  if (pool) {
    if (is.null(catch_by_gear)) {
      catch_by_gear <- rep(1, length(lfd))
      message("lbi_flicc: no catch_by_gear; gears pooled with equal catch shares.")
    }
    if (catch_units == "weight" && !all(c("a", "b") %in% dimnames(lhpar)[[1]]))
      stop('catch_units = "weight" needs a and b in lhpar (or use catch_units = "numbers").')
    pooled <- combine_lfd(lfd, catch_by_gear, lhpar = lhpar, ess.n = NULL, name = "pooled",
                          catch_units = catch_units)
    lfd <- FLCore::FLQuants(c(as.list(lfd), as.list(pooled)))
  }
  gears <- if (is.null(gear)) names(lfd) else if (is.numeric(gear)) names(lfd)[gear] else gear
  a <- getp("a"); b <- getp("b")
  rows <- list()
  for (g in gears) {
    m <- lfd_matrix_flicc(lfd[[g]])
    L <- as.numeric(rownames(m)); step <- if (length(L) > 1) L[2] - L[1] else 1
    wt <- if (is.finite(a) && is.finite(b)) a * (L + step / 2)^b else (L + step / 2)^3
    yrs <- colnames(m); if (!is.null(years)) yrs <- intersect(yrs, as.character(years))
    for (y in yrs) {
      r <- lbi_one(L, L + step, m[, y], linf, L50, Lo, k, wt)
      if (!is.null(r)) rows[[length(rows) + 1]] <- cbind(gear = g, year = as.numeric(y), r)
    }
  }
  if (!length(rows)) stop("No compositions with data.")
  vals <- do.call(rbind, rows)
  vals$gear <- factor(vals$gear, levels = gears)
  cp <- cope_punt_flicc(vals$Pmat, vals$Popt, vals$Pmega, L50 / Lo)
  vals$cp_selectivity <- cp$cp_selectivity
  vals$cp_test <- cp$cp_test
  vals$cp_status <- cp$cp_status
  def <- lbi_table_def()
  if (!is.null(limits)) {
    bad <- setdiff(names(limits), def$indicator)
    if (length(bad)) stop("Unknown indicator(s) in limits: ", paste(bad, collapse = ", "))
    def$lim[match(names(limits), def$indicator)] <- as.numeric(limits)
  }
  rat <- do.call(rbind, lapply(seq_len(nrow(def)), function(i) {
    v <- vals[[def$num[i]]] / (if (is.na(def$den[i])) 1 else vals[[def$den[i]]])
    st <- ifelse(v >= def$ref[i], "green", ifelse(v >= def$lim[i], "orange", "red"))
    data.frame(gear = vals$gear, year = vals$year, indicator = def$indicator[i], value = v,
               ref = def$ref[i], lim = def$lim[i], status = st, pass = st == "green",
               property = def$property[i])
  }))
  cp_ok <- ifelse(vals$cp_status == "no information", NA, vals$cp_status == "SB >= target")
  rat <- rbind(rat, data.frame(gear = vals$gear, year = vals$year, indicator = "Pobj (Cope-Punt)",
                               value = vals$Pobj, ref = NA_real_, lim = NA_real_,
                               status = ifelse(is.na(cp_ok), NA, ifelse(cp_ok, "green", "red")),
                               pass = cp_ok, property = "Spawning biomass vs 0.4 B0"))
  rat$indicator <- factor(rat$indicator, levels = c(def$indicator, "Pobj (Cope-Punt)"))
  rel <- NULL
  if (!is.null(ref_years)) {
    cols <- c("Lc", "L25", "Lbar", "Lmean", "LmaxY", "L95", "Lmax5", "Pmat", "Popt", "Pmega", "Z")
    rel <- do.call(rbind, lapply(split(vals, vals$gear), function(d) {
      ref <- d[d$year %in% ref_years, cols, drop = FALSE]
      if (!nrow(ref)) return(NULL)
      out <- d[, c("gear", "year")]
      for (cc in cols) out[[cc]] <- d[[cc]] / mean(ref[[cc]], na.rm = TRUE)
      out
    }))
  }
  out <- list(values = vals, ratios = rat, relative = rel, ref_years = ref_years,
              Lopt = Lo, linf = linf, L50 = L50)
  class(out) <- c("flicc_lbi", "list")
  out
}

#' Plot length-based indicators
#'
#' @param x Output of [lbi_flicc()].
#' @param type `"traffic"` (default): one tile per indicator ratio and year, green
#'   where the reference level is met, red where not, with the value; one panel
#'   per gear. `"trend"`: the ratios over years with their reference levels
#'   (dashed), one panel per indicator, one colour per gear. `"relative"`:
#'   indicators relative to the reference period (needs `ref_years` in
#'   [lbi_flicc()]).
#' @param indicators Optional subset of ratio names (e.g. `c("Lmax5/Linf",
#'   "Pmega", "Lmean/LFM")`).
#' @return A `ggplot`.
#' @export
plot_lbi <- function(x, type = c("traffic", "trend", "relative"), indicators = NULL) {
  type <- match.arg(type)
  r <- x$ratios
  if (!is.null(indicators)) r <- r[r$indicator %in% indicators, ]
  if (type == "traffic") {
    r$lab <- formatC(r$value, digits = 2, format = "f")
    r$indicator <- factor(r$indicator, levels = rev(levels(droplevels(r$indicator))))
    r$status <- factor(r$status, levels = c("green", "orange", "red"))
    return(ggplot2::ggplot(r, ggplot2::aes(factor(year), indicator, fill = status)) +
             ggplot2::geom_tile(colour = "white") +
             ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.6) +
             ggplot2::scale_fill_manual(values = c(green = "#7FBF7B", orange = "#FDB863", red = "#EF8A62"),
                                        labels = c(green = "meets reference level", orange = "between limit and reference",
                                                   red = "below limit"),
                                        name = NULL, na.value = "grey85", drop = FALSE) +
             ggplot2::facet_wrap(~gear, ncol = 1) +
             ggplot2::labs(x = "Year", y = NULL, title = "Length-based indicators") +
             ggplot2::theme_bw() +
             ggplot2::theme(legend.position = "bottom", panel.grid = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)))
  }
  if (type == "trend") {
    refs <- unique(r[, c("indicator", "ref")])
    refs <- refs[is.finite(refs$ref), ]
    if ("Pobj (Cope-Punt)" %in% r$indicator)
      refs <- rbind(refs, data.frame(indicator = "Pobj (Cope-Punt)", ref = c(1, 2)))
    lims <- unique(r[is.finite(r$lim), c("indicator", "lim")])
    return(ggplot2::ggplot(r, ggplot2::aes(year, value, colour = gear)) +
             ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = ref), linetype = "dashed",
                                 colour = "darkgreen") +
             ggplot2::geom_hline(data = lims, ggplot2::aes(yintercept = lim), linetype = "dotted",
                                 colour = "red") +
             ggplot2::geom_line() + ggplot2::geom_point(size = 1.2) +
             ggplot2::facet_wrap(~indicator, scales = "free_y") + year_axis_flicc() +
             ggplot2::labs(x = "Year", y = "Indicator ratio", colour = "Gear",
                           title = "Length-based indicators against reference levels") +
             ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom"))
  }
  if (is.null(x$relative)) stop('type = "relative" needs ref_years in lbi_flicc().')
  rl <- stats::reshape(x$relative, direction = "long", varying = setdiff(names(x$relative), c("gear", "year")),
                       v.names = "value", timevar = "indicator",
                       times = setdiff(names(x$relative), c("gear", "year")), idvar = c("gear", "year"))
  rl <- rl[rl$indicator %in% c("Lmean", "L95", "Lmax5", "Pmega", "Pmat", "Z"), ]
  yr <- range(x$ref_years)
  ggplot2::ggplot(rl, ggplot2::aes(year, value, colour = gear)) +
    ggplot2::annotate("rect", xmin = yr[1] - 0.5, xmax = yr[2] + 0.5, ymin = -Inf, ymax = Inf,
                      fill = "grey85", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~indicator, scales = "free_y") + year_axis_flicc() +
    ggplot2::labs(x = "Year", y = "Relative to the reference period", colour = "Gear",
                  title = "Indicators relative to the reference period (shaded)") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
