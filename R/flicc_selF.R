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
  function(Fl) {
    N <- popN(Fl) / s0
    Z <- M + Fl
    C <- ifelse(Z > 0, N * Fl / Z * (1 - exp(-Z)), 0)
    c(SPR = sum(N * mat * w) / SB0, YPR = sum(C * w))
  }
}

# selectivity parameters (length scale) and current F by gear
gear_state_flicc <- function(fit, nyears = 1) {
  rp <- fit$report
  td <- fit$tmb_data
  gears <- td$gear_names
  type <- c("logistic", "dsnormal", "normal")[td$sel_type]
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
  } else {
    pars$mode <- loc
  }
  pars
}

sel_loc_flicc <- function(type, pars) if (type == "logistic") pars$SL50 else pars$mode

#' Yield- and spawning-potential-per-recruit surfaces over selectivity and F
#'
#' For each gear, computes YPR and SPR at equilibrium over a grid of the gear's
#' selectivity location (`SL50` for logistic, the mode for dome-shaped
#' selectivity; the shape is kept, i.e. the curve is shifted along the length
#' axis) and the gear's apical fishing mortality. The other gears are held at
#' their fitted selectivity and current F (mean of the last `nyears`). Life
#' history, growth-type groups and natural mortality are those of the fit.
#'
#' @param fit A fitted FLicc object (from [fiticc()], with `FLRreport = TRUE`).
#' @param gear Gear names (or indices) to compute; default all gears.
#' @param n Grid size in each dimension.
#' @param Fmax Upper limit of F; default per gear `2 * max(Fcur, F at SPR 20%)`.
#'   A single value or one per gear.
#' @param Lrange Range of the selectivity location; default from the first
#'   length bin to Linf.
#' @param nyears Number of terminal years for current F, M, weight and maturity.
#' @return A list of class `"flicc_contour"` with `surface` (data frame: gear,
#'   type, loc, F, SPR, YPR, YPRrel = YPR relative to the maximum over the
#'   gear's grid), `current` (gear, type, loc, F, SPR, YPR, YPRrel) and `Lmat`
#'   (length at 50% maturity, if available).
#' @examples
#' \dontrun{
#' cf <- contour_flicc(fit)
#' plot_contour_flicc(cf)
#' plot_contour_flicc(cf, s = 0.75)
#' }
#' @export
contour_flicc <- function(fit, gear = NULL, n = 40, Fmax = NULL, Lrange = NULL, nyears = 1) {
  if (is.null(fit$report$selpars)) stop("fit needs FLRreport = TRUE (the default of fiticc()).")
  eng <- pr_engine_flicc(fit, nyears)
  st <- gear_state_flicc(fit, nyears)
  gears <- st$gears
  gsel <- if (is.null(gear)) gears else if (is.numeric(gear)) gears[gear] else gear
  bad <- setdiff(gsel, gears)
  if (length(bad)) stop("Unknown gear(s): ", paste(bad, collapse = ", "), ". Gears: ", paste(gears, collapse = ", "))
  if (!is.null(Fmax)) Fmax <- rep_len(Fmax, length(gsel))
  L <- st$Lmid
  S <- lapply(gears, function(g) sel_curve_flicc(L, st$type[[g]], st$pars[[g]]))
  names(S) <- gears
  linf <- as.numeric(fit$report$lhpar["linf"])
  if (is.null(Lrange)) Lrange <- c(as.numeric(fit$tmb_data$LLB[1]), linf)

  surf <- list(); cur <- list()
  for (i in seq_along(gsel)) {
    g <- gsel[i]
    Foth <- Reduce(`+`, lapply(setdiff(gears, g), function(h) st$F[[h]] * S[[h]]), rep(0, length(L)))
    loc0 <- sel_loc_flicc(st$type[[g]], st$pars[[g]])
    pr_at <- function(loc, Fg) eng(Foth + Fg * sel_curve_flicc(L, st$type[[g]], shift_sel_flicc(st$type[[g]], st$pars[[g]], loc)))
    Fmx <- if (!is.null(Fmax)) Fmax[i] else {
      f20 <- tryCatch(stats::uniroot(function(f) pr_at(loc0, f)[["SPR"]] - 0.2,
                                     c(0, max(50 * st$F[[g]], 5)))$root, error = function(e) NA_real_)
      2 * max(st$F[[g]], f20, na.rm = TRUE)
    }
    locs <- seq(Lrange[1], Lrange[2], length.out = n)
    Fs <- seq(0, Fmx, length.out = n)
    grid <- expand.grid(loc = locs, F = Fs)
    v <- t(mapply(pr_at, grid$loc, grid$F))
    d <- data.frame(gear = g, type = st$type[[g]], grid, SPR = v[, "SPR"], YPR = v[, "YPR"])
    ymax <- max(d$YPR, na.rm = TRUE)
    d$YPRrel <- d$YPR / ymax
    v0 <- pr_at(loc0, st$F[[g]])
    surf[[g]] <- d
    cur[[g]] <- data.frame(gear = g, type = st$type[[g]], loc = loc0, F = st$F[[g]],
                           SPR = v0[["SPR"]], YPR = v0[["YPR"]], YPRrel = v0[["YPR"]] / ymax)
  }
  Lmat <- if ("L50" %in% dimnames(fit$report$lhpar)[[1]]) as.numeric(fit$report$lhpar["L50"]) else NA_real_
  structure(list(surface = do.call(rbind, surf), current = do.call(rbind, cur), Lmat = Lmat,
                 linf = linf, nyears = nyears),
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

#' Contour plots of yield and spawning potential over selectivity and F
#'
#' Draws, for each gear (rows), yield (left) and spawning potential (right) as
#' functions of the gear's apical F and its selectivity location (`SL50` or
#' mode), with the other gears at their current state. Both panels are on the
#' same 0-100% scale:
#'
#' * without `s`: YPR% = 100 x YPR / max(YPR over the grid) and SPR% = 100 x SPR;
#' * with steepness `s` (Beverton-Holt): Yield% = 100 x Y / max(Y over the
#'   grid), with Y = YPR x R/R0, and SSB% = 100 x B/B0, with
#'   R/R0 = (4 s SPR - (1 - s)) / ((5 s - 1) SPR) and
#'   B/B0 = (4 s SPR - (1 - s)) / (5 s - 1) (both 0 when negative).
#'
#' Styled after the FLSelex isopleth plots: a continuous `rich.colors` fill
#' (dark blue = high, red = low), labelled grey isopleths, reference levels of
#' the right-hand quantity (`ref`, e.g. SPR 40% or SSB 40%) as black lines on
#' both panels, the current state as a black point with dotted guides to the
#' axes, and the length at 50% maturity as a dotted horizontal line.
#'
#' @param x Output of [contour_flicc()] or a fitted FLicc object (then passed to
#'   [contour_flicc()] with `...`).
#' @param gear Gears to show (rows); default all in `x`.
#' @param s Optional steepness of a Beverton-Holt stock-recruitment
#'   relationship; gives Yield% and SSB% instead of YPR% and SPR%.
#' @param ref Reference levels (%) of SPR (or SSB with `s`) drawn as black lines;
#'   the first solid, others dashed.
#' @param levels Isopleth levels (%), the same on both panels.
#' @param ... Passed to [contour_flicc()] when `x` is a fit.
#' @return A `patchwork` object (rows = gears, columns = yield, spawning
#'   potential), or a list of `ggplot`s if patchwork is not installed.
#' @examples
#' \dontrun{
#' plot_contour_flicc(fit)
#' plot_contour_flicc(fit, s = 0.75)
#' plot_contour_flicc(fit, gear = "Gillnet", n = 60)
#' }
#' @export
plot_contour_flicc <- function(x, gear = NULL, s = NULL, ref = c(40, 20),
                               levels = c(10, 20, 30, 50, 60, 80, 90, 95), ...) {
  if (!inherits(x, "flicc_contour")) x <- contour_flicc(x, gear = gear, ...)
  gears <- if (is.null(gear)) unique(x$surface$gear) else gear
  bad <- setdiff(gears, unique(x$surface$gear))
  if (length(bad)) stop("Gear(s) not in the surface: ", paste(bad, collapse = ", "))
  if (!is.null(s)) check_steep_flicc(s, 40)

  # yield and spawning potential in % (0-100) on a common scale
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
  nm_left  <- if (is.null(s)) "YPR%" else "Yield%"
  nm_right <- if (is.null(s)) "SPR%" else "SSB%"
  ref_nm   <- if (is.null(s)) "SPR" else "SSB"
  ref_lt <- c("solid", rep("22", max(0, length(ref) - 1)))

  # one FLSelex fill scale for both panels: uneven breaks give more resolution
  colbr <- c(0, 5, seq(10, 100, 10))
  fill_pct <- ggplot2::scale_fill_gradientn(
    colours = rev(r4col_flicc(length(colbr) + 4))[-c(4:7)],
    limits = c(-3, 100), breaks = colbr, name = "%")

  lab_layer <- function(lab, size = 2.6, face = "plain") {
    if (is.null(lab)) return(NULL)
    ggplot2::geom_label(data = lab, ggplot2::aes(F, loc, label = label),
                        size = size, fontface = face, label.size = 0,
                        fill = "white", alpha = 0.75,
                        label.padding = ggplot2::unit(0.06, "lines"))
  }

  panel <- function(g, what) {
    d <- x$surface[x$surface$gear == g, ]
    cp <- x$current[x$current$gear == g, ]
    # maximum yield over the grid (YPR without s, YPR x R/R0 with s)
    ymax <- max(to_pct(d, 1)$left, na.rm = TRUE) / 100
    d <- to_pct(d, ymax); cp <- to_pct(cp, ymax)
    d$z <- if (what == "left") d$left else d$right
    ylab <- if (d$type[1] == "logistic") "SL50 (length at 50% selectivity)" else "Selectivity mode"

    p <- ggplot2::ggplot(d, ggplot2::aes(F, loc)) +
      ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = TRUE) + fill_pct +
      ggplot2::geom_contour(ggplot2::aes(z = z), breaks = setdiff(levels, ref),
                            colour = grDevices::grey(0.4), linewidth = 0.35)
    lab <- contour_labels_flicc(d, "z", setdiff(levels, ref))

    # reference levels of the right-hand quantity as black lines on both panels
    ref_lab <- NULL
    for (j in seq_along(ref)) {
      p <- p + ggplot2::geom_contour(ggplot2::aes(z = right), breaks = ref[j], colour = "black",
                                     linewidth = if (j == 1) 0.8 else 0.55, linetype = ref_lt[j])
      rl <- contour_labels_flicc(d, "right", ref[j], function(v) paste0(ref_nm, v),
                                 placed = rbind(lab, ref_lab))
      if (!is.null(rl)) ref_lab <- rbind(ref_lab, rl[1, ])
    }
    p <- p + lab_layer(lab) + lab_layer(ref_lab, size = 2.7, face = "bold")

    if (is.finite(x$Lmat)) {
      p <- p + ggplot2::geom_hline(yintercept = x$Lmat, linetype = "dotted", colour = "grey20") +
        ggplot2::annotate("text", x = 0, y = x$Lmat, label = "Lmat",
                          hjust = -0.15, vjust = -0.4, size = 2.7)
    }

    sub <- if (what == "left") sprintf("current %s = %.0f%%", sub("%", "", nm_left), cp$left)
           else sprintf("current %s = %.0f%%", sub("%", "", nm_right), cp$right)
    # current state: black point with dotted guides to the axes (as FLSelex)
    p + ggplot2::annotate("segment", x = cp$F, xend = cp$F, y = min(d$loc), yend = cp$loc,
                          linetype = "dotted", linewidth = 0.3) +
      ggplot2::annotate("segment", x = 0, xend = cp$F, y = cp$loc, yend = cp$loc,
                        linetype = "dotted", linewidth = 0.3) +
      ggplot2::annotate("point", x = cp$F, y = cp$loc, size = 2.5) +
      ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
      ggplot2::labs(x = paste0("Fishing Mortality (", g, ")"), y = ylab,
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
  # all panels share one colour bar, collected into a single legend
  patchwork::wrap_plots(plots, ncol = 2, byrow = TRUE) +
    patchwork::plot_layout(guides = "collect")
}
