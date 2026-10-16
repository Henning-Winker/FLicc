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
#' @return A list of class `"flicc_selF"` with `surface` (data frame: gear,
#'   type, loc, F, SPR, YPR, YPRrel = YPR relative to the maximum over the
#'   gear's grid), `current` (gear, type, loc, F, SPR, YPR, YPRrel) and `Lmat`
#'   (length at 50% maturity, if available).
#' @examples
#' \dontrun{
#' sf <- selF_surface_flicc(fit)
#' plot_selF_flicc(sf)
#' }
#' @export
selF_surface_flicc <- function(fit, gear = NULL, n = 40, Fmax = NULL, Lrange = NULL, nyears = 1) {
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
            class = "flicc_selF")
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

#' Plot YPR and SPR isopleths over selectivity and F by gear
#'
#' Draws, for each gear (rows), yield per recruit (left) and spawning potential
#' ratio (right) as functions of the gear's apical F and its selectivity
#' location (`SL50` or mode), with the other gears at their current state.
#' Styled after the FLSelex isopleth plots: a continuous `rich.colors` fill
#' (dark blue = high, red = low), labelled grey isopleths, the SPR reference
#' levels (`spr_ref`) as black lines, the current state as a black point with
#' dotted guides to the axes, and the length at 50% maturity as a dotted
#' horizontal line.
#'
#' @param x Output of [selF_surface_flicc()] or a fitted FLicc object (then
#'   passed to [selF_surface_flicc()] with `...`).
#' @param gear Gears to show (rows); default all in `x`.
#' @param spr_ref SPR reference levels (percent) drawn as black lines; the
#'   first is drawn solid (target), others dashed.
#' @param breaks SPR isopleth levels (proportions); 0 and 1 are ignored.
#' @param ypr_levels Relative YPR isopleth levels.
#' @param ypr_fill `"ypr"` (default): left panel shaded by relative YPR;
#'   `"spr"`: shaded by SPR.
#' @param ... Passed to [selF_surface_flicc()] when `x` is a fit.
#' @return A `patchwork` object (rows = gears, columns = YPR, SPR), or a list
#'   of `ggplot`s if patchwork is not installed.
#' @examples
#' \dontrun{
#' plot_selF_flicc(fit)
#' plot_selF_flicc(fit, gear = "Gillnet", n = 60)
#' }
#' @export
plot_selF_flicc <- function(x, gear = NULL, spr_ref = c(40, 20),
                            breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8),
                            ypr_levels = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95),
                            ypr_fill = c("ypr", "spr"), ...) {
  ypr_fill <- match.arg(ypr_fill)
  if (!inherits(x, "flicc_selF")) x <- selF_surface_flicc(x, gear = gear, ...)
  gears <- if (is.null(gear)) unique(x$surface$gear) else gear
  bad <- setdiff(gears, unique(x$surface$gear))
  if (length(bad)) stop("Gear(s) not in the surface: ", paste(bad, collapse = ", "))

  spr_lev <- breaks[breaks > 0 & breaks < 1]
  ref_lev <- spr_ref / 100
  ref_lt <- c("solid", rep("22", max(0, length(ref_lev) - 1)))

  # FLSelex fill scales: uneven breaks give more colour resolution near the top
  colbr_y <- c(0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
  fill_ypr <- ggplot2::scale_fill_gradientn(
    colours = rev(r4col_flicc(length(colbr_y) + 3))[-c(10:11, 13)],
    limits = c(-0.03, 1), breaks = colbr_y, name = "YPR")
  colbr_s <- c(0, 0.05, seq(0.1, 1, 0.1))
  fill_spr <- ggplot2::scale_fill_gradientn(
    colours = rev(r4col_flicc(length(colbr_s) + 4))[-c(4:7)],
    limits = c(-0.03, 1), breaks = colbr_s, name = "SPR")

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
    ylab <- if (d$type[1] == "logistic") "SL50 (length at 50% selectivity)" else "Selectivity mode"
    fill_ypr_here <- what == "YPR" && ypr_fill == "ypr"

    d$zfill <- if (fill_ypr_here) d$YPRrel else d$SPR

    p <- ggplot2::ggplot(d, ggplot2::aes(F, loc)) +
      ggplot2::geom_raster(ggplot2::aes(fill = zfill), interpolate = TRUE) +
      (if (fill_ypr_here) fill_ypr else fill_spr)

    # grey isopleths with labels
    if (what == "YPR") {
      p <- p + ggplot2::geom_contour(ggplot2::aes(z = YPRrel), breaks = ypr_levels,
                                     colour = grDevices::grey(0.4), linewidth = 0.35)
      lab <- contour_labels_flicc(d, "YPRrel", ypr_levels)
    } else {
      p <- p + ggplot2::geom_contour(ggplot2::aes(z = SPR), breaks = spr_lev,
                                     colour = grDevices::grey(0.4), linewidth = 0.35)
      lab <- contour_labels_flicc(d, "SPR", setdiff(spr_lev, ref_lev))
    }

    # SPR reference levels as black lines, labelled away from the grey labels
    ref_lab <- NULL
    for (j in seq_along(ref_lev)) {
      p <- p + ggplot2::geom_contour(ggplot2::aes(z = SPR), breaks = ref_lev[j], colour = "black",
                                     linewidth = if (j == 1) 0.8 else 0.55, linetype = ref_lt[j])
      rl <- contour_labels_flicc(d, "SPR", ref_lev[j], function(v) paste0("SPR", 100 * v),
                                 placed = rbind(lab, ref_lab))
      if (!is.null(rl)) ref_lab <- rbind(ref_lab, rl[1, ])
    }
    p <- p + lab_layer(lab) + lab_layer(ref_lab, size = 2.7, face = "bold")

    if (is.finite(x$Lmat)) {
      p <- p + ggplot2::geom_hline(yintercept = x$Lmat, linetype = "dotted", colour = "grey20") +
        ggplot2::annotate("text", x = max(d$F), y = x$Lmat, label = "Lmat",
                          hjust = 1.1, vjust = -0.4, size = 2.7)
    }

    # current state: black point with dotted guides to the axes (as FLSelex)
    p + ggplot2::annotate("segment", x = cp$F, xend = cp$F, y = min(d$loc), yend = cp$loc,
                          linetype = "dotted", linewidth = 0.3) +
      ggplot2::annotate("segment", x = 0, xend = cp$F, y = cp$loc, yend = cp$loc,
                        linetype = "dotted", linewidth = 0.3) +
      ggplot2::annotate("point", x = cp$F, y = cp$loc, size = 2.5) +
      ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
      ggplot2::labs(x = paste0("Fishing Mortality (", g, ")"), y = ylab,
                    title = paste0(what, ": ", g),
                    subtitle = if (what == "YPR") sprintf("current YPR = %.0f%% of max", 100 * cp$YPRrel)
                               else sprintf("current SPR = %.0f%%", 100 * cp$SPR)) +
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

  plots <- unlist(lapply(gears, function(g) list(panel(g, "YPR"), panel(g, "SPR"))), recursive = FALSE)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(plots)
  # identical colour bars across gears are merged into one YPR and one SPR legend
  patchwork::wrap_plots(plots, ncol = 2, byrow = TRUE) +
    patchwork::plot_layout(guides = "collect")
}
