# Per-recruit surfaces over gear selectivity and F (YPR, SPR isopleths)

# Fast per-recruit engine at the fitted life history: returns a function of
# fishing mortality at length (annual, length nlen) giving SPR and YPR. Mirrors
# the reference-point block of the TMB template (GTG or gamma recursion, catch
# N F / Z (1 - exp(-Z)), SBPR = sum(N mat w)).
pr_engine_flicc <- function(fit, nyears = 1) {
  rp <- fit$report
  td <- fit$tmb_data
  yrs <- dimnames(rp$N)$year
  yrs <- utils::tail(yrs, min(nyears, length(yrs)))
  Len <- as.numeric(dimnames(rp$N)$len)
  nlen <- length(Len)
  k   <- as.numeric(rp$lhpar["k"])
  M   <- as.numeric(FLCore::yearMeans(rp$M[, yrs]))
  w   <- as.numeric(FLCore::yearMeans(rp$w_at_len[, yrs]))
  mat <- as.numeric(FLCore::yearMeans(FLCore::mat(fit$stklen)[, yrs]))
  pm  <- pop_model_flicc(fit)

  if (pm == "gtg") {
    gtg <- gtg_inputs_flicc(fit)
    ms <- td$Mscaler
    Mscale <- as.numeric(rowMeans(ms[, utils::tail(seq_len(ncol(ms)), length(yrs)), drop = FALSE]))
    step <- if (nlen > 1) Len[nlen] - Len[nlen - 1] else 1
    LB <- c(Len, Len[nlen] + step)
    Linfs <- gtg$gtgLinfs
    ngtg <- length(Linfs)
    lnum <- outer(LB[-1], Linfs, function(L, Li) pmax(Li - L, 1e-12))
    lden <- outer(LB[-(nlen + 1)], Linfs, function(L, Li) pmax(Li - L, 1e-12))
    lr <- log(lnum / lden)                                  # nlen x ngtg, independent of F
    MK <- gtg$MKMat * c(Mscale, Mscale[nlen])               # (nlen + 1) x ngtg
    recP <- gtg$recP
    popN <- function(Fl) {
      ZK <- MK + c(Fl, Fl[nlen]) / k
      lS <- apply(ZK[seq_len(nlen), , drop = FALSE] * lr, 2, cumsum)
      if (!is.matrix(lS)) lS <- matrix(lS, nrow = nlen)
      P <- rbind(recP, sweep(exp(lS), 2, recP, "*"))
      rowSums((P[seq_len(nlen), , drop = FALSE] - P[-1, , drop = FALSE]) / ZK[seq_len(nlen), , drop = FALSE])
    }
  } else {
    ref <- flicc_refpars(fit, nyears = nyears)
    popN <- function(Fl) n_f_flicc(ref$node, ref$quad_wt, ref$Len, (M + Fl) / k, ref$Galpha, ref$Gbeta)
  }
  N0 <- popN(rep(0, nlen))
  s0 <- N0[1]
  SB0 <- sum(N0 / s0 * mat * w)
  function(Fl, detail = FALSE) {
    N <- popN(Fl) / s0
    Z <- M + Fl
    C <- ifelse(Z > 0, N * Fl / Z * (1 - exp(-Z)), 0)
    out <- c(SPR = sum(N * mat * w) / SB0, YPR = sum(C * w))
    if (detail) return(list(pr = out, Cw = C * w, Fl = Fl))
    out
  }
}

# selectivity parameters (length scale) and current F by gear
gear_state_flicc <- function(fit, nyears = 1) {
  rp <- fit$report
  td <- fit$tmb_data
  gears <- td$gear_names
  type <- c("logistic", "dsnormal", "normal", "binormal", "bidsnormal")[td$sel_type]
  k <- as.numeric(rp$lhpar["k"])
  yrs <- utils::tail(dimnames(rp$N)$year, nyears)
  pars <- lapply(gears, function(g) { p <- rp$selpars[[g]]; stats::setNames(as.list(as.numeric(p)), dimnames(p)[[1]]) })
  Fg <- vapply(gears, function(g) k * mean(as.numeric(rp$Fk[[g]][, yrs])), numeric(1))
  names(pars) <- names(type) <- gears
  list(gears = gears, type = type, pars = pars, F = Fg, Lmid = as.numeric(td$Lmid))
}

shift_sel_flicc <- function(type, pars, loc) {
  if (type == "logistic") {
    w <- pars$SL95 - pars$SL50
    pars$SL50 <- loc; pars$SL95 <- loc + w
  } else if (type %in% c("binormal", "bidsnormal")) {
    d <- loc - pars$mode1
    pars$mode1 <- loc; pars$mode2 <- pars$mode2 + d
  } else {
    pars$mode <- loc
  }
  pars
}

sel_loc_flicc <- function(type, pars) if (type == "logistic") pars$SL50 else if (type %in% c("binormal", "bidsnormal")) pars$mode1 else pars$mode

#' Yield and spawning-potential surfaces over selectivity, F and F shares
#'
#' Computes equilibrium yield per recruit (YPR) and SPR over a grid for each
#' gear, in one of three views:
#'
#' * `view = "gear"` (default): the gear's apical F (x) against its selectivity
#'   location (y; `SL50` for logistic, the mode for dome-shaped selectivity,
#'   the curve shifted along the length axis with its shape kept). The other
#'   gears stay at their current F and selectivity, so the catch split changes
#'   along the F axis.
#' * `view = "fleet"`: total F (x; sum of the apical F of all gears) with the
#'   gears' F shares fixed at their current values, against the gear's
#'   selectivity location (y). The other gears keep their selectivity.
#' * `view = "share"`: total F (x) against the gear's share of the total F (y,
#'   0-100%); the rest of the total F is split among the other gears in their
#'   current proportions. All selectivities are as fitted.
#'
#' Life history, growth-type groups and natural mortality are those of the fit;
#' current values are means of the last `nyears`.
#'
#' @param fit A fitted FLicc object (from [fiticc()], with `FLRreport = TRUE`).
#' @param gear Gear names (or indices) to compute; default all gears.
#' @param view `"gear"`, `"fleet"` or `"share"` (see Details).
#' @param n Grid size in each dimension.
#' @param Fmax Upper limit of the F axis; default `2 * max(current F, F at
#'   SPR 20%)` along the view's F direction. A single value or one per gear.
#' @param Lrange Range of the selectivity location; default from the first
#'   length bin to Linf (not used for `view = "share"`).
#' @param nyears Number of terminal years for current F, M, weight and maturity.
#' @return A list of class `"flicc_contour"` with `surface` (data frame: gear,
#'   type, F = x-axis value, loc = y-axis value, SPR, YPR, YPRrel), `current`
#'   (one row per gear: the current state on the same axes, with the gear's
#'   current F share `Fshare` and catch share `Cshare`), `view`, `Lmat` and
#'   `linf`.
#' @examples
#' \dontrun{
#' plot_contour_flicc(contour_flicc(fit))                     # one gear changes
#' plot_contour_flicc(contour_flicc(fit, view = "fleet"))     # total F, shares fixed
#' plot_contour_flicc(contour_flicc(fit, view = "share"))     # total F x F share
#' }
#' @export
contour_flicc <- function(fit, gear = NULL, view = c("fleet", "gear", "share"), n = 40,
                          Fmax = NULL, Lrange = NULL, nyears = 1) {
  view <- match.arg(view)
  if (is.null(fit$report$selpars)) stop("fit needs FLRreport = TRUE (the default of fiticc()).")
  eng <- pr_engine_flicc(fit, nyears)
  st <- gear_state_flicc(fit, nyears)
  gears <- st$gears
  gsel <- if (is.null(gear)) gears else if (is.numeric(gear)) gears[gear] else gear
  bad <- setdiff(gsel, gears)
  if (length(bad)) stop("Unknown gear(s): ", paste(bad, collapse = ", "), ". Gears: ", paste(gears, collapse = ", "))
  if (view == "share" && length(gears) < 2) stop("view = \"share\" needs at least two gears.")
  if (!is.null(Fmax)) Fmax <- rep_len(Fmax, length(gsel))
  L <- st$Lmid
  S <- lapply(gears, function(g) sel_curve_flicc(L, st$type[[g]], st$pars[[g]]))
  names(S) <- gears
  linf <- as.numeric(fit$report$lhpar["linf"])
  if (is.null(Lrange)) Lrange <- c(as.numeric(fit$tmb_data$LLB[1]), linf)
  Fcur <- unlist(st$F)
  Ftot <- sum(Fcur)
  share <- Fcur / Ftot

  # current catch share of each gear (weight)
  dcur <- eng(Reduce(`+`, Map(function(f, sg) f * sg, Fcur, S)), detail = TRUE)
  Cshare <- vapply(gears, function(g) {
    part <- ifelse(dcur$Fl > 0, Fcur[[g]] * S[[g]] / dcur$Fl, 0)
    sum(dcur$Cw * part) / sum(dcur$Cw)
  }, numeric(1))

  surf <- list(); cur <- list()
  for (i in seq_along(gsel)) {
    g <- gsel[i]
    oth <- setdiff(gears, g)
    Sg_at <- function(loc) sel_curve_flicc(L, st$type[[g]], shift_sel_flicc(st$type[[g]], st$pars[[g]], loc))
    loc0 <- sel_loc_flicc(st$type[[g]], st$pars[[g]])
    if (view == "gear") {
      Foth <- Reduce(`+`, lapply(oth, function(h) Fcur[[h]] * S[[h]]), rep(0, length(L)))
      pr_at <- function(y, x) eng(Foth + x * Sg_at(y))
      x0 <- Fcur[[g]]; y0 <- loc0
    } else if (view == "fleet") {
      Soth <- Reduce(`+`, lapply(oth, function(h) share[[h]] * S[[h]]), rep(0, length(L)))
      pr_at <- function(y, x) eng(x * (share[[g]] * Sg_at(y) + Soth))
      x0 <- Ftot; y0 <- loc0
    } else {
      wo <- if (sum(Fcur[oth]) > 0) Fcur[oth] / sum(Fcur[oth]) else rep(1 / length(oth), length(oth))
      Smix <- Reduce(`+`, Map(function(wh, h) wh * S[[h]], wo, oth))
      pr_at <- function(y, x) eng(x * (y / 100 * S[[g]] + (1 - y / 100) * Smix))
      x0 <- Ftot; y0 <- 100 * share[[g]]
    }
    Fmx <- if (!is.null(Fmax)) Fmax[i] else {
      f20 <- tryCatch(stats::uniroot(function(f) pr_at(y0, f)[["SPR"]] - 0.2,
                                     c(0, max(50 * x0, 5)))$root, error = function(e) NA_real_)
      2 * max(x0, f20, na.rm = TRUE)
    }
    ys <- if (view == "share") seq(0, 100, length.out = n) else seq(Lrange[1], Lrange[2], length.out = n)
    grid <- expand.grid(loc = ys, F = seq(0, Fmx, length.out = n))
    v <- t(mapply(pr_at, grid$loc, grid$F))
    d <- data.frame(gear = g, type = st$type[[g]], grid, SPR = v[, "SPR"], YPR = v[, "YPR"])
    ymax <- max(d$YPR, na.rm = TRUE)
    d$YPRrel <- d$YPR / ymax
    v0 <- pr_at(y0, x0)
    surf[[g]] <- d
    cur[[g]] <- data.frame(gear = g, type = st$type[[g]], loc = y0, F = x0,
                           SPR = v0[["SPR"]], YPR = v0[["YPR"]], YPRrel = v0[["YPR"]] / ymax,
                           Fshare = share[[g]], Cshare = Cshare[[g]])
  }
  Lmat <- if ("L50" %in% dimnames(fit$report$lhpar)[[1]]) as.numeric(fit$report$lhpar["L50"]) else NA_real_
  structure(list(surface = do.call(rbind, surf), current = do.call(rbind, cur), view = view,
                 Lmat = Lmat, linf = linf, nyears = nyears),
            class = "flicc_contour")
}


#' FLSelex colour palette
#'
#' Subset of `rich.colors` (Arni Magnusson), as used by FLSelex for isopleth
#' plots. Reversed (`rev()`), it runs from dark blue (high) to red (low).
#'
#' @param n Number of colours.
#' @param alpha Transparency.
#' @return Character vector of colour codes.
#' @keywords internal
r4col_flicc <- function(n, alpha = 1) {
  x <- seq(0, 1, length.out = n)
  r <- 1 / (1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- stats::dnorm(x, 0.25, 0.15) / max(stats::dnorm(x, 0.25, 0.15))
  grDevices::rgb(r, g, b, alpha = alpha)
}

# Contour line labels: one label per contour line, placed at the point that is
# farthest from the panel edges and from labels already placed, so labels do
# not stack up on steep line ends or get clipped at the borders.
contour_labels_flicc <- function(d, z, levels, fmt = function(x) x, placed = NULL) {
  x <- sort(unique(d$loc)); y <- sort(unique(d$F))
  zm <- matrix(d[[z]][order(d$F, d$loc)], nrow = length(x))
  cl <- grDevices::contourLines(x, y, zm, levels = levels)
  if (!length(cl)) return(NULL)
  rx <- range(x); ry <- range(y)
  nx <- function(v) (v - rx[1]) / diff(rx)
  ny <- function(v) (v - ry[1]) / diff(ry)
  out <- NULL
  for (l in cl) {
    xs <- nx(l$x); ys <- ny(l$y)
    score <- pmin(xs, 1 - xs, ys, 1 - ys)
    prev <- rbind(placed, out)
    if (!is.null(prev) && nrow(prev)) {
      dmin <- apply(cbind(xs, ys), 1, function(p)
        min(sqrt((nx(prev$loc) - p[1])^2 + (ny(prev$F) - p[2])^2)))
      score <- pmin(score, dmin)
    }
    j <- which.max(score)
    out <- rbind(out, data.frame(loc = l$x[j], F = l$y[j], label = fmt(l$level)))
  }
  out
}

#' Contour plots of yield and spawning potential over selectivity, F and F shares
#'
#' Draws, for each gear (rows), yield (left) and spawning potential (right) on
#' the axes of the [contour_flicc()] view: the gear's F against its
#' selectivity location (`"gear"`), total F with fixed F shares against the
#' gear's selectivity location (`"fleet"`), or total F against the gear's share
#' of total F (`"share"`). Values are in %:
#'
#' * without `s`: YPR% = 100 x YPR / max(YPR over the grid) and SPR% = 100 x SPR;
#' * with steepness `s` (Beverton-Holt): Yield% = 100 x Y / max(Y over the
#'   grid), with Y = YPR x R/R0, and SSB% = 100 x B/B0, with
#'   R/R0 = (4 s SPR - (1 - s)) / ((5 s - 1) SPR) and
#'   B/B0 = (4 s SPR - (1 - s)) / (5 s - 1) (both 0 when negative).
#'
#' Each column has its own colour scale. Spawning potential uses a diverging
#' scale anchored on the reference levels: red below the limit (`ref[2]`),
#' orange to pale yellow between limit and target (`ref[1]`), blue above the
#' target. Yield uses the sequential viridis scale stretched over `ylim` (yield
#' below `ylim[1]` takes the lowest colour), with more resolution near the
#' maximum. Grey labelled isopleths, the reference levels as black lines on
#' both panels, the current state as a point with dotted guides to the axes,
#' and the length at 50% maturity (dotted, when the y-axis is a length).
#'
#' @param x Output of [contour_flicc()] or a fitted FLicc object (then passed to
#'   [contour_flicc()] with `...`, e.g. `view = "share"`).
#' @param gear Gears to show (rows); default all in `x`.
#' @param s Optional steepness of a Beverton-Holt stock-recruitment
#'   relationship; gives Yield% and SSB% instead of YPR% and SPR%.
#' @param ref Target and limit (%) of SPR (or SSB with `s`), drawn as black
#'   lines (target solid, limit dashed) and anchoring the colour scale.
#' @param levels Isopleth levels (%) of the spawning-potential panels.
#' @param ylevels Isopleth levels (%) of the yield panels.
#' @param ylim Range (%) of the yield colour scale; default from the 5%
#'   quantile of the yield shown (rounded down to 10) to 100.
#' @param ... Passed to [contour_flicc()] when `x` is a fit.
#' @return A `patchwork` object (rows = gears, columns = yield, spawning
#'   potential), or a list of `ggplot`s if patchwork is not installed.
#' @examples
#' \dontrun{
#' plot_contour_flicc(fit)
#' plot_contour_flicc(fit, view = "fleet")
#' plot_contour_flicc(fit, view = "share", gear = "Gillnet")
#' plot_contour_flicc(fit, s = 0.75)
#' }
#' @export
plot_contour_flicc <- function(x, gear = NULL, s = NULL, ref = c(40, 20),
                               levels = c(10, 30, 50, 60, 80),
                               ylevels = c(50, 70, 80, 90, 95), ylim = NULL, ...) {
  if (!inherits(x, "flicc_contour")) x <- contour_flicc(x, gear = gear, ...)
  view <- if (is.null(x$view)) "gear" else x$view
  gears <- if (is.null(gear)) unique(x$surface$gear) else gear
  bad <- setdiff(gears, unique(x$surface$gear))
  if (length(bad)) stop("Gear(s) not in the surface: ", paste(bad, collapse = ", "))
  if (!is.null(s)) check_steep_flicc(s, 40)
  ref <- sort(ref, decreasing = TRUE)
  tgt <- ref[1]; lim <- if (length(ref) > 1) ref[2] else tgt / 2

  # yield and spawning potential in % (0-100)
  to_pct <- function(d, ymax) {
    if (is.null(s)) {
      d$left  <- 100 * d$YPR / ymax
      d$right <- 100 * d$SPR
    } else {
      rr0 <- (4 * s * d$SPR - (1 - s)) / ((5 * s - 1) * d$SPR)
      rr0[!is.finite(rr0) | rr0 < 0] <- 0
      d$left  <- 100 * d$YPR * rr0 / ymax
      d$right <- 100 * pmax(4 * s * d$SPR - (1 - s), 0) / (5 * s - 1)
    }
    d
  }
  pct <- lapply(gears, function(g) {
    d <- x$surface[x$surface$gear == g, ]
    cp <- x$current[x$current$gear == g, ]
    ymax <- max(to_pct(d, 1)$left, na.rm = TRUE) / 100
    list(d = to_pct(d, ymax), cp = to_pct(cp, ymax))
  })
  names(pct) <- gears
  nm_left  <- if (is.null(s)) "YPR%" else "Yield%"
  nm_right <- if (is.null(s)) "SPR%" else "SSB%"
  ref_nm   <- if (is.null(s)) "SPR" else "SSB"
  ref_lt   <- c("solid", rep("22", max(0, length(ref) - 1)))

  # spawning potential: diverging, anchored on limit and target
  stops <- c(0, lim / 2, lim, (lim + tgt) / 2, tgt, (tgt + 100) / 2, 100)
  fill_right <- ggplot2::scale_fill_gradientn(
    colours = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FFFFBF", "#74ADD1", "#313695"),
    values = stops / 100, limits = c(0, 100), breaks = sort(unique(c(0, lim, tgt, 60, 80, 100))),
    oob = scales::squish, name = nm_right)
  # yield: viridis stretched over the range shown, more resolution near the top
  if (is.null(ylim)) {
    lo <- stats::quantile(unlist(lapply(pct, function(z) z$d$left)), 0.05, na.rm = TRUE)
    ylim <- c(max(0, min(floor(lo / 10) * 10, 80)), 100)
  }
  ystops <- ylim[1] + (ylim[2] - ylim[1]) * c(0, 0.4, 0.6, 0.75, 0.85, 0.92, 1)
  fill_left <- ggplot2::scale_fill_gradientn(
    colours = scales::viridis_pal()(7), values = (ystops - ylim[1]) / diff(ylim),
    limits = ylim, oob = scales::squish, name = nm_left,
    breaks = unique(round(c(ylim[1], ystops[3:7]))))

  lab_layer <- function(lab, size = 2.6, face = "plain") {
    if (is.null(lab)) return(NULL)
    ggplot2::geom_label(data = lab, ggplot2::aes(F, loc, label = label),
                        size = size, fontface = face, label.size = 0,
                        fill = "white", alpha = 0.75,
                        label.padding = ggplot2::unit(0.06, "lines"))
  }

  panel <- function(g, what) {
    d <- pct[[g]]$d; cp <- pct[[g]]$cp
    d$z <- if (what == "left") d$left else d$right
    lev <- setdiff(if (what == "left") ylevels else levels, ref)
    ylab <- switch(view,
      share = paste0("Share of total F: ", g, " (%)"),
      if (d$type[1] == "logistic") paste0("SL50, ", g) else paste0("Selectivity mode, ", g))
    xlab <- if (view == "gear") paste0("F (", g, ")") else
      if (view == "fleet") "Total F (F shares as current)" else "Total F"

    p <- ggplot2::ggplot(d, ggplot2::aes(F, loc)) +
      ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = TRUE) +
      (if (what == "left") fill_left else fill_right) +
      ggplot2::geom_contour(ggplot2::aes(z = z), breaks = lev,
                            colour = grDevices::grey(0.35), linewidth = 0.35)
    lab <- contour_labels_flicc(d, "z", lev)

    ref_lab <- NULL
    for (j in seq_along(ref)) {
      p <- p + ggplot2::geom_contour(ggplot2::aes(z = right), breaks = ref[j], colour = "black",
                                     linewidth = if (j == 1) 0.8 else 0.55, linetype = ref_lt[j])
      rl <- contour_labels_flicc(d, "right", ref[j], function(v) paste0(ref_nm, v),
                                 placed = rbind(lab, ref_lab))
      if (!is.null(rl)) ref_lab <- rbind(ref_lab, rl[1, ])
    }
    p <- p + lab_layer(lab) + lab_layer(ref_lab, size = 2.7, face = "bold")

    if (view != "share" && is.finite(x$Lmat)) {
      p <- p + ggplot2::geom_hline(yintercept = x$Lmat, linetype = "dotted", colour = "grey20") +
        ggplot2::annotate("text", x = 0, y = x$Lmat, label = "Lmat",
                          hjust = -0.15, vjust = -0.4, size = 2.7)
    }

    sub <- if (what == "left") sprintf("current %s = %.0f%%", sub("%", "", nm_left), cp$left)
           else sprintf("current %s = %.0f%%", sub("%", "", nm_right), cp$right)
    if (view == "share" && !is.null(cp$Cshare))
      sub <- sprintf("%s; %s: %.0f%% of F, %.0f%% of catch", sub, g, 100 * cp$Fshare, 100 * cp$Cshare)
    p + ggplot2::annotate("segment", x = cp$F, xend = cp$F, y = min(d$loc), yend = cp$loc,
                          linetype = "dotted", linewidth = 0.3) +
      ggplot2::annotate("segment", x = 0, xend = cp$F, y = cp$loc, yend = cp$loc,
                        linetype = "dotted", linewidth = 0.3) +
      ggplot2::annotate("point", x = cp$F, y = cp$loc, size = 2.5) +
      ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
      ggplot2::labs(x = xlab, y = ylab,
                    title = paste0(if (what == "left") nm_left else nm_right, ": ", g),
                    subtitle = sub) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11),
                     plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
                     axis.title = ggplot2::element_text(size = 10),
                     panel.grid = ggplot2::element_blank(),
                     legend.key.height = ggplot2::unit(1, "cm"),
                     legend.key.width = ggplot2::unit(0.6, "cm"),
                     legend.text = ggplot2::element_text(size = 7),
                     legend.title = ggplot2::element_text(size = 9),
                     plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5))
  }

  plots <- unlist(lapply(gears, function(g) list(panel(g, "left"), panel(g, "right"))), recursive = FALSE)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(plots)
  # one colour bar per column, collected into two legends
  patchwork::wrap_plots(plots, ncol = 2, byrow = TRUE) +
    patchwork::plot_layout(guides = "collect")
}
