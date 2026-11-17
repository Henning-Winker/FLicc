# Diagnostics of the age-length keys (conditional age at length) of an FLicc fit with growth from ALKs

# observed counts and expected proportions of the ALKs of a fit, one row per key, length bin and age
alk_cells_flicc <- function(fit) {
  if (is.null(fit$alk) || !isTRUE(fit$tmb_data$alk_on == 1L))
    stop("the fit has no age-length keys (fiticc(..., alk = )).")
  td <- fit$tmb_data
  obs <- td$alk_obs; pex <- fit$alk_pexp
  if (is.null(pex)) stop("the fit has no expected ALK proportions: refit with this FLicc version.")
  LLB <- as.numeric(td$LLB)
  bw <- if (length(LLB) > 1) stats::median(diff(LLB)) else 1
  ages <- td$alk_age
  do.call(rbind, lapply(seq_along(fit$alk), function(j) {
    k <- fit$alk[[j]]
    nl <- colSums(obs[, , j, drop = FALSE][, , 1])
    do.call(rbind, lapply(which(nl > 0), function(l)
      data.frame(key = j, gears = paste(k$gears, collapse = "+"), year = k$year, timing = k$timing,
                 len = LLB[l], lmid = LLB[l] + bw / 2, age = ages, obs = obs[, l, j], pexp = pex[, l, j],
                 n = nl[l])))
  }))
}

#' Residuals of the age-length keys
#'
#' Residuals of the conditional age-at-length data of an FLicc fit with growth from age-length keys
#' (`fiticc(..., alk = )`). Each length bin of a key is a multinomial sample over ages, as in the likelihood.
#'
#' * `type = "osa"` (default): one-step-ahead (randomized quantile) residuals, decomposing each length bin
#'   along increasing age into conditional binomials (Trijoulet et al. 2023, Fisheries Research 257: 106487).
#'   Under a correct model they are independent standard normal; the oldest age of each bin has no residual.
#' * `type = "pearson"`: \eqn{(o - n p) / \sqrt{n p (1 - p)}}.
#'
#' @param fit A fitted FLicc object with age-length keys.
#' @param type `"osa"` or `"pearson"`.
#' @param seed Random seed for the OSA randomization.
#' @param min_exp OSA residuals of cells with no fish observed and an expected count (n p) below `min_exp` are
#'   set to `NA`: the randomized residual of a cell that can hardly be observed (old ages in small length bins,
#'   young ages in large ones) is pure noise. 0 keeps all (as in ss3diags2).
#' @return A data frame with key, gears, year, len (lower bound), lmid, age, obs, exp (n p), pexp, n (aged fish in
#'   the length bin) and resid.
#' @export
resid_alk_flicc <- function(fit, type = c("osa", "pearson"), seed = 1, min_exp = 0.01) {
  type <- match.arg(type)
  d <- alk_cells_flicc(fit)
  if (type == "osa") {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  d$resid <- NA_real_
  for (g in split(seq_len(nrow(d)), paste(d$key, d$len))) {
    o <- round(d$obs[g]); p <- pmax(d$pexp[g], 1e-12); p <- p / sum(p); n <- sum(o)
    if (type == "pearson") {
      d$resid[g] <- (o - n * p) / sqrt(n * p * (1 - p))
    } else {
      rem_n <- n; rem_p <- 1; r <- rep(NA_real_, length(o))
      for (a in seq_len(length(o) - 1)) {
        if (rem_n <= 0 || rem_p <= 1e-12) break
        q <- min(max(p[a] / rem_p, 1e-12), 1 - 1e-12)
        lo <- stats::pbinom(o[a] - 1, rem_n, q); hi <- stats::pbinom(o[a], rem_n, q)
        r[a] <- stats::qnorm(min(max(lo + stats::runif(1) * (hi - lo), 1e-12), 1 - 1e-12))
        rem_n <- rem_n - o[a]; rem_p <- rem_p - p[a]
      }
      d$resid[g] <- r
    }
  }
  d$exp <- d$n * d$pexp
  if (type == "osa" && min_exp > 0) d$resid[d$obs == 0 & d$exp < min_exp] <- NA
  attr(d, "type") <- type
  d
}

# mean length at age of the aged fish: observed, and expected under the same numbers aged per length bin
alk_means_flicc <- function(d) {
  do.call(rbind, lapply(split(d, list(d$key, d$age), drop = TRUE), function(z) {
    no <- sum(z$obs); ne <- sum(z$n * z$pexp)
    if (no < 1 || ne < 0.5) return(NULL)                 # ages hardly expected in the sample: no mean
    data.frame(key = z$key[1], gears = z$gears[1], year = z$year[1], age = z$age[1], nobs = no,
               obs = sum(z$lmid * z$obs) / no, exp = sum(z$lmid * z$n * z$pexp) / ne)
  }))
}

# naive von Bertalanffy fit to the age-length pairs of one key (selection ignored), t0 fixed
vb_naive_flicc <- function(z, t0, timing) {
  a <- rep(z$age + timing, z$obs); L <- rep(z$lmid, z$obs)
  if (length(unique(a)) < 3) return(c(linf = NA, k = NA))
  f <- try(stats::nls(L ~ linf * (1 - exp(-k * (a - t0))), start = list(linf = max(L), k = 0.3),
                      algorithm = "port", lower = c(0.5 * max(L), 0.01), upper = c(3 * max(L), 5)), silent = TRUE)
  if (inherits(f, "try-error")) c(linf = NA, k = NA) else stats::coef(f)
}

# growth curve of the fit at the sampling time, with +/- 2 SD of length at age (CVL)
alk_curve_flicc <- function(fit, amax, timing) {
  g <- fit$growth; cvl <- fit$alk_cvl %||% 0.1
  a <- seq(0, amax + 0.5, by = 0.05)
  L <- g[["linf"]] * (1 - exp(-g[["k"]] * (a + timing - g[["t0"]])))
  data.frame(age = a, L = pmax(L, 0), sd = cvl * pmax(L, 0))
}

#' Plot the fit to the age-length keys
#'
#' Observed age-length data of each key (bubbles: fish aged by age and length bin) with
#'
#' * the observed mean length at age of the aged fish (diamonds) and the mean length at age FLicc expects for
#'   the same numbers aged per length bin (bars): these compare with each other, as both carry the sampling;
#' * the population growth curve estimated by FLicc (solid, with \eqn{\pm} 2 SD of length at age);
#' * a naive von Bertalanffy curve fitted directly to the age-length pairs (dashed, t0 fixed as in the fit),
#'   as an external growth fit would give. Where the aged fish come from a selective gear or are sampled by
#'   length, the naive curve follows the sample rather than the population.
#'
#' @param fit A fitted FLicc object with age-length keys.
#' @param naive Add the naive von Bertalanffy curve.
#' @param band Shade \eqn{\pm} 2 SD of length at age around the FLicc curve.
#' @return A `ggplot` object.
#' @export
plot_alk_fit <- function(fit, naive = TRUE, band = TRUE) {
  d <- alk_cells_flicc(fit)
  mn <- alk_means_flicc(d)
  lab <- function(k) paste0(k$year, " (", k$gears, ")")
  keys <- unique(d[, c("key", "gears", "year", "timing")])
  keys$panel <- lab(keys)
  d$panel <- keys$panel[match(d$key, keys$key)]; mn$panel <- keys$panel[match(mn$key, keys$key)]
  t0 <- fit$growth[["t0"]]
  cur <- do.call(rbind, lapply(seq_len(nrow(keys)), function(j) {
    cc <- alk_curve_flicc(fit, max(d$age[d$obs > 0]), keys$timing[j]); cc$panel <- keys$panel[j]; cc$what <- "FLicc"
    if (naive) {
      vb <- vb_naive_flicc(d[d$key == keys$key[j], ], t0, keys$timing[j])
      if (all(is.finite(vb))) {
        nv <- cc; nv$L <- vb[["linf"]] * (1 - exp(-vb[["k"]] * (nv$age + keys$timing[j] - t0)))
        nv$sd <- NA; nv$what <- sprintf("naive VB (Linf %.1f, k %.2f)", vb[["linf"]], vb[["k"]])
        cc <- rbind(cc, nv)
      }
    }
    cc
  }))
  lmax <- max(d$lmid[d$obs > 0]) * 1.15
  cur <- cur[cur$L <= lmax, ]
  g <- fit$growth
  flab <- sprintf("FLicc (Linf %.1f, k %.2f)", g[["linf"]], g[["k"]])
  cur$what[cur$what == "FLicc"] <- flab
  cols <- stats::setNames(c("#08519c", rep("#D6604D", length(unique(cur$what)) - 1)), c(flab, setdiff(unique(cur$what), flab)))
  lts <- stats::setNames(c(1, rep(2, length(cols) - 1)), names(cols))
  p <- ggplot2::ggplot(d[d$obs > 0, ], ggplot2::aes(age, lmid))
  if (band) p <- p + ggplot2::geom_ribbon(data = cur[cur$what == flab, ],
                                          ggplot2::aes(age, ymin = pmax(L - 2 * sd, 0), ymax = L + 2 * sd),
                                          inherit.aes = FALSE, fill = "#08519c", alpha = 0.1)
  p + ggplot2::geom_point(ggplot2::aes(size = obs), shape = 21, fill = "grey70", colour = "grey35", alpha = 0.7) +
    ggplot2::geom_line(data = cur, ggplot2::aes(age, L, colour = what, linetype = what), inherit.aes = FALSE,
                       linewidth = 0.8) +
    ggplot2::geom_segment(data = mn, ggplot2::aes(x = age - 0.2, xend = age + 0.2, y = exp, yend = exp),
                          inherit.aes = FALSE, linewidth = 0.7) +
    ggplot2::geom_point(data = mn, ggplot2::aes(age, obs), inherit.aes = FALSE, shape = 23, fill = "gold",
                        size = 2, stroke = 0.4) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_linetype_manual(values = lts, name = NULL) +
    ggplot2::scale_size_area(max_size = 5, name = "Fish aged") +
    ggplot2::facet_wrap(~panel) +
    ggplot2::labs(x = "Age (years)", y = "Length (cm)", title = "Fit to the age-length keys",
                  caption = paste("diamonds: observed, bars: expected mean length at age of the aged fish;",
                                  "band: +/- 2 SD of length at age")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}

#' Plot residuals of the age-length keys
#'
#' OSA (default) or Pearson residuals of the conditional age-at-length data by age (x) and length bin (y), one
#' panel per key: bubbles sized by the absolute residual, blue for observed > expected, red for observed <
#' expected. The growth curve of the fit (\eqn{\pm} 2 SD of length at age) and the observed (diamonds) and
#' expected (bars) mean length at age of the aged fish are drawn on top. Runs of same-signed residuals along
#' the curve point to misfit in growth (Linf, k or CVL); with `qq = TRUE` a QQ plot of the residuals instead.
#'
#' @param x A fitted FLicc object with age-length keys, or the output of [resid_alk_flicc()].
#' @param type `"osa"` or `"pearson"` (when `x` is a fit).
#' @param fit The fit, for the growth curve, when `x` is a residual data frame.
#' @param means Add observed and expected mean length at age.
#' @param cap Residuals capped at this absolute value for the bubble size.
#' @param qq Draw a QQ plot of the residuals instead of the bubbles.
#' @return A `ggplot` object.
#' @export
plot_alk_resid <- function(x, type = c("osa", "pearson"), fit = NULL, means = TRUE, cap = 4, qq = FALSE) {
  type <- match.arg(type)
  if (!is.data.frame(x)) { fit <- x; x <- resid_alk_flicc(x, type = type) }
  type <- attr(x, "type") %||% type
  keys <- unique(x[, c("key", "gears", "year", "timing")])
  keys$panel <- paste0(keys$year, " (", keys$gears, ")")
  x$panel <- keys$panel[match(x$key, keys$key)]
  d <- x[is.finite(x$resid), ]
  if (qq) return(ggplot2::ggplot(d, ggplot2::aes(sample = resid)) +
                   ggplot2::stat_qq(size = 0.8, alpha = 0.6) + ggplot2::stat_qq_line(colour = "#D6604D") +
                   ggplot2::facet_wrap(~panel) +
                   ggplot2::labs(x = "Theoretical N(0,1) quantiles", y = "Residual quantiles",
                                 title = paste("QQ plot,", type, "residuals of the age-length keys")) +
                   ggplot2::theme_bw())
  d$r <- pmin(abs(d$resid), cap)
  d$sign <- ifelse(d$resid > 0, "obs > exp", "obs < exp")
  p <- ggplot2::ggplot(d, ggplot2::aes(age, lmid))
  if (!is.null(fit)) {
    cur <- do.call(rbind, lapply(seq_len(nrow(keys)), function(j) {
      cc <- alk_curve_flicc(fit, max(x$age[x$obs > 0]), keys$timing[j]); cc$panel <- keys$panel[j]; cc }))
    cur <- cur[cur$L <= max(x$lmid) * 1.1, ]
    p <- p + ggplot2::geom_ribbon(data = cur, ggplot2::aes(age, ymin = pmax(L - 2 * sd, 0), ymax = L + 2 * sd),
                                  inherit.aes = FALSE, fill = "grey50", alpha = 0.12) +
      ggplot2::geom_line(data = cur, ggplot2::aes(age, L), inherit.aes = FALSE, colour = "grey40")
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(size = r, fill = sign), shape = 21, colour = "grey30", stroke = 0.2,
                               alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(`obs > exp` = "#3C78B4", `obs < exp` = "#D2604B"), name = NULL) +
    ggplot2::scale_size_area(max_size = 4, limits = c(0, cap), name = "|residual|")
  if (means) {
    mn <- alk_means_flicc(x); mn$panel <- keys$panel[match(mn$key, keys$key)]
    p <- p + ggplot2::geom_segment(data = mn, ggplot2::aes(x = age - 0.18, xend = age + 0.18, y = exp, yend = exp),
                                   inherit.aes = FALSE, linewidth = 0.6) +
      ggplot2::geom_point(data = mn, ggplot2::aes(age, obs), inherit.aes = FALSE, shape = 23, fill = "gold",
                          size = 1.6, stroke = 0.4)
  }
  p + ggplot2::facet_wrap(~panel) +
    ggplot2::labs(x = "Age (years)", y = "Length (cm)", title = "Residuals of the age-length keys",
                  subtitle = paste(type, "residuals, conditional age at length"),
                  caption = paste(c(if (!is.null(fit)) "band: +/- 2 SD of length at age",
                                    if (means) "diamonds: observed, bars: expected mean length at age of the aged fish"),
                                  collapse = "; ")) +
    ggplot2::theme_bw()
}
