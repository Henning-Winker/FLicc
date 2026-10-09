
#' plot_hcrflicc
#'
#' Plot an SPR-based harvest control rule (HCR) for FLicc.
#'
#' The x-axis shows current SPR and the y-axis shows the relative change
#' in TAC or effort. The horizontal zero line indicates no change.
#'
#' Stock status zones are defined independently from the HCR triggers:
#'   - below blim: critical
#'   - blim to bthr: cautious / rebuilding
#'   - above bthr: healthyW
#'
#' The HCR is defined by four trigger points:
#'   - below b1: fixed low change (dlow)
#'   - b1 to b2: linear increase from dlow to dopt
#'   - b2 to b3: constant at dopt
#'   - b3 to b4: linear increase from dopt to dup
#'   - above b4: constant at dup
#'
#' @param b1 Lower HCR trigger for reduced effort/TAC.
#' @param b2 HCR trigger where no-change plateau starts.
#' @param b3 HCR trigger where upper slope starts.
#' @param b4 HCR trigger where maximum increase is reached.
#' @param blim Lower biological limit for background zone.
#' @param bthr Biological threshold separating cautious and healthy zones.
#' @param dlow Relative change below b1, e.g. -0.20 for -20%.
#' @param dopt Relative change in the target range, usually 0.
#' @param dup Relative change above b4, e.g. 0.10 for +10%.
#' @param metric Character string for y-axis label: "Effort" or "TAC".
#' @param spr_target Optional target SPR shown as a vertical dashed line.
#' @param show_spr Logical, if TRUE annotate the SPR target.
#' @param obs Optional numeric vector of observed SPR values.
#' @param years Optional year labels for obs.
#' @param xmax Maximum x-axis value.
#' @param ymin Minimum y-axis value.
#' @param ymax Maximum y-axis value.
#' @param alpha Transparency for background shading.
#' @param triggers Logical, if TRUE annotate trigger points and levels.
#' @param refpts Logical, if TRUE annotate trigger points and levels.
#' @param status.text Logical, if TRUE add stock-zone labels.
#' @param line_colour Colour of HCR curve.
#' @param line_size Line width of HCR curve.
#'
#' @return A ggplot object.
#' @examples
#'
#' plot_hcrflicc(
#'   spr_target = 0.4,
#'   show_target = TRUE,
#'   spr_target_label = "SPR target"
#' )
#' # example code
#'
#' @export
plot_hcrspr <- function(
    b1 = 0.1,
    b2 = 0.35,
    b3 = 0.45,
    b4 = 0.6,
    blim = 0.10,
    bthr = 0.2,
    dlow = -0.20,
    dopt =  0.00,
    dup  =  0.15,
    metric = c("TAC","Effort"),
    spr_target = 0.4,
    show_target = TRUE,
    spr_target_label = "SPR target",
    obs = NULL,
    years = NULL,
    xmax = 0.9,
    ymin = -0.30,
    ymax =  0.30,
    alpha = 0.18,
    triggers = TRUE,
    refpts = TRUE,
    status.text = FALSE,
    line_colour = "blue",
    line_size = 1.2
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  metric <- match.arg(metric)

  stopifnot(
    is.numeric(c(b1, b2, b3, b4, blim, bthr, dlow, dopt, dup)),
    b1 < b2,
    b2 <= b3,
    b3 < b4,
    blim < bthr
  )

  # background zones based on stock status, not HCR triggers
  zone_dat <- data.frame(
    xmin = c(0, blim, bthr),
    xmax = c(blim, bthr, xmax),
    ymin = ymin,
    ymax = ymax,
    zone = factor(
      c("Critical", "Cautious", "Healthy"),
      levels = c("Critical", "Cautious", "Healthy")
    )
  )

  # piecewise HCR curve
  hcr_fun <- function(x) {
    ifelse(
      x <= b1, dlow,
      ifelse(
        x < b2, dlow + (dopt - dlow) * (x - b1) / (b2 - b1),
        ifelse(
          x <= b3, dopt,
          ifelse(
            x < b4, dopt + (dup - dopt) * (x - b3) / (b4 - b3),
            dup
          )
        )
      )
    )
  }

  xseq <- seq(0, xmax, length.out = 500)
  hcr_dat <- data.frame(
    spr = xseq,
    change = hcr_fun(xseq)
  )

  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_rect(
      data = zone_dat,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone),
      alpha = alpha,
      colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Critical" = "indianred2",
        "Cautious" = "wheat3",
        "Healthy"  = "lightblue3"
      ),
      guide = "none"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    ggplot2::geom_line(
      data = hcr_dat,
      ggplot2::aes(x = spr, y = change),
      colour = line_colour,
      linewidth = line_size
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, xmax),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(ymin, ymax),
      expand = c(0, 0),
      labels = function(z) paste0(round(z * 100), "%")
    ) +
    ggplot2::xlab(expression(SPR[curr])) +
    ggplot2::ylab(paste("Relative change in", metric)) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )

  p <- p + ggplot2::annotate("text", x = xmax * 0.98, y = dlow, label = paste0(round(100 * dlow), "%"), hjust = 1, vjust = -0.2) +
    ggplot2::annotate("text", x = xmax * 0.98, y = dopt, label = paste0(round(100 * dopt), "%"), hjust = 1, vjust = -0.2) +
    ggplot2::annotate("text", x = xmax * 0.98, y = dup,  label = paste0("+", round(100 * dup), "%"), hjust = 1, vjust = -0.2)

  if (!is.null(spr_target)) {
    p <- p +
      ggplot2::geom_segment(
        aes(
          x = spr_target, xend = spr_target,
          y = dopt, yend = dup*0.75
        ),
        linetype = 5,
        linewidth = 0.7,
        colour = "blue"
      )


      p <- p +
        ggplot2::annotate(
          "text",
          x = spr_target,
          y = dup*0.75 + 0.01 * (ymax - ymin),
          label = spr_target_label,
          colour = "blue",
          hjust = 0.5,
          vjust = 0,
          size = 3.5
        )

  }

  if (triggers) {
    p <- p +
      ggplot2::annotate("text", x = b1, y = ymin + 0.03 * (ymax - ymin), label = "b1", hjust = 1.1) +
      ggplot2::annotate("text", x = b2, y = ymin + 0.03 * (ymax - ymin), label = "b2", hjust = 1.1) +
      ggplot2::annotate("text", x = b3, y = ymin + 0.03 * (ymax - ymin), label = "b3", hjust = 1.1) +
      ggplot2::annotate("text", x = b4, y = ymin + 0.03 * (ymax - ymin), label = "b4", hjust = 1.1) +
      ggplot2::geom_vline(xintercept = c(b1, b2, b3, b4), linewidth = 0.4, colour = "blue",linetype=2)
  }
  if(refpts) {
    p <- p + ggplot2::geom_vline(xintercept = c(blim, bthr), linewidth = 0.6, colour = "black")
  }

  if (status.text) {
    p <- p +
      ggplot2::annotate("text", x = blim / 2, y = ymax * 0.80, label = "Critical", angle = 90) +
      ggplot2::annotate("text", x = (blim + bthr) / 2, y = ymax * 0.80, label = "Cautious", angle = 90) +
      ggplot2::annotate("text", x = (bthr + xmax) / 2, y = ymax * 0.80, label = "Healthy", angle = 90)
  }

  if (!is.null(obs)) {
    obs <- as.numeric(obs)
    obs_dat <- data.frame(
      spr = obs,
      change = hcr_fun(obs),
      idx = seq_along(obs)
    )

    p <- p +
      ggplot2::geom_path(
        data = obs_dat,
        ggplot2::aes(x = spr, y = change, group = 1),
        colour = "grey40",
        alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = obs_dat,
        ggplot2::aes(x = spr, y = change),
        shape = 21,
        fill = c(rep("grey80", max(0, nrow(obs_dat) - 1)), "blue"),
        colour = "black",
        size = c(rep(2, max(0, nrow(obs_dat) - 1)), 3)
      )

    if (!is.null(years) && length(years) == length(obs)) {
      obs_dat$year <- years
      p <- p +
        ggplot2::geom_text(
          data = obs_dat,
          ggplot2::aes(x = spr, y = change, label = year),
          nudge_y = 0.02,
          size = 3
        )
    }
  }

  return(p)
}


#' hcr_spr (flicc)
#'
#' SPR-based harvest control rule for FLicc.
#'
#' Returns a multiplicative TAC/effort adjustment factor based on current SPR.
#' The rule follows the same piecewise structure as plot_hcrflicc():
#'   - spr <= b1: dlow
#'   - b1 < spr < b2: linear slope from dlow to dopt
#'   - b2 <= spr <= b3: dopt
#'   - b3 < spr < b4: linear slope from dopt to dhi
#'   - spr >= b4: dhi
#'
#' Relative changes are converted to multipliers as:
#'   multiplier = 1 + change
#'
#' @param spr Numeric SPR value(s).
#' @param b1 Lower trigger point.
#' @param b2 Trigger where no-change plateau starts.
#' @param b3 Trigger where upper slope starts.
#' @param b4 Trigger where maximum increase is reached.
#' @param dlow Relative change below b1, e.g. -0.20.
#' @param dopt Relative change between b2 and b3, usually 0.
#' @param dhi Relative change above b4, e.g. 0.10.
#' @param multiplier Logical, if TRUE return multiplier (default),
#'   otherwise return relative change.
#'
#' @return Numeric vector of multipliers or relative changes.
#'
#' @examples
#' # example code
#'
#' spr <- seq(0,1,0.02)
#' change <- spr_rule(spr=seq(0,1,0.02),
#' b1 = 0.1,
#' b2 = 0.35,
#' b3 = 0.45,
#' b4 = 0.6,
#' dlow = -0.20,
#' dopt =  0.00,
#' dhi  =  0.15)
#' plot(spr,change-1,type="l",col=4,lwd=2,ylab=("Change TAC"),xlab="SPR")
#' abline(0,0,lty=2)
#' plot(spr,change,type="l",col=4,lwd=2,ylab=("Multiplier"),xlab="SPR")
#' abline(0,0,lty=2)
#'
#' @export
spr_rule <- function(
    spr,
    b1 = 0.10,
    b2 = 0.35,
    b3 = 0.45,
    b4 = 0.6,
    dlow = -0.20,
    dopt =  0.00,
    dhi  =  0.15,
    multiplier = TRUE
) {
  stopifnot(
    is.numeric(spr),
    is.numeric(c(b1, b2, b3, b4, dlow, dopt, dhi)),
    b1 < b2,
    b2 <= b3,
    b3 < b4
  )

  change <- ifelse(
    spr <= b1, dlow,
    ifelse(
      spr < b2,
      dlow + (dopt - dlow) * (spr - b1) / (b2 - b1),
      ifelse(
        spr <= b3, dopt,
        ifelse(
          spr < b4,
          dopt + (dhi - dopt) * (spr - b3) / (b4 - b3),
          dhi
        )
      )
    )
  )

  if (multiplier) {
    return(1 + change)
  } else {
    return(change)
  }
}
#' Biomass ratio trend rule (2 over 3)
#'
#' Calculate a simple trend ratio from a biomass-related index time series.
#'
#' This rule compares the mean of the most recent \code{n1} years of an index
#' with the mean of the preceding \code{n2} years, returning the ratio:
#'
#' \deqn{
#' \bar{I}_{recent} / \bar{I}_{previous}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\bar{I}_{recent}} is the mean over the last \code{n1} years
#'   \item \eqn{\bar{I}_{previous}} is the mean over the \code{n2} years
#'         immediately before that
#' }
#'
#' Values greater than 1 indicate that the recent biomass index is higher than
#' the earlier period, while values below 1 indicate a declining trend.
#'
#' @param idx An FLQuant-like index object, for example one element of
#'   \code{LBIspr(fit)}.
#' @param n1 Number of most recent years used to calculate the current mean.
#'   Default is 2.
#' @param n2 Number of preceding years used as the comparison period.
#'   Default is 3.
#'
#' @return A numeric scalar returned via \code{an()}, giving the ratio of the
#'   recent mean index to the previous mean index.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item extracting the index values with \code{index(idx)}
#'   \item taking the last \code{n1} years and calculating their mean
#'   \item taking the \code{n2} years immediately before those and calculating
#'         their mean
#'   \item returning the ratio of recent mean to previous mean
#' }
#'
#' This can be used as a simple decision rule component for biomass or
#' recruitment indicators.
#'
#' @examples
#' # Biomass ratio trend rule
#' #
#' # idx <- LBIspr(fit)[["Gillnet"]]
#' # r_rule(idx)
#' # r_rule(idx, n1 = 3, n2 = 5)
#'
#' @export
r_rule <- function(idx, n1 = 2, n2 = 3) {
  an(
    yearMeans(tail(index(idx), n1)) /
      yearMeans(head(tail(index(idx), n1 + n2), n2))
  )
}

#' Combined SPR and index trend harvest control rule
#'
#' Apply a combined harvest control rule using current SPR and a recent
#' biomass-ratio trend indicator.
#'
#' The rule first calculates a gear-specific biomass trend ratio using
#' \code{r_rule()} applied to \code{LBIspr(fit, thresh = thresh)[[gear]]}.
#' It then calculates an SPR-based change multiplier using \code{spr_rule()}.
#'
#' The SPR multiplier is applied only when the biomass trend and SPR rule point
#' in the same direction:
#' \itemize{
#'   \item if \code{r < 1} and \code{s < 1}, apply \code{s}
#'   \item if \code{r > 1} and \code{s > 1}, apply \code{s}
#'   \item otherwise, return \code{1}, meaning no change
#' }
#'
#' @param fit A fitted FLicc object.
#' @param gear Numeric or character index selecting the gear-specific indicator
#'   from \code{LBIspr(fit, thresh = thresh)}.
#' @param thresh Threshold passed to \code{LBIspr()}.
#' @param b1 Lower SPR trigger point.
#' @param b2 SPR trigger where the no-change plateau starts.
#' @param b3 SPR trigger where the upper slope starts.
#' @param b4 SPR trigger where the maximum increase is reached.
#' @param dlow Relative change below \code{b1}, e.g. \code{-0.20} for a
#'   20% reduction.
#' @param dopt Relative change between \code{b2} and \code{b3}, usually
#'   \code{0}.
#' @param dhi Relative change above \code{b4}, e.g. \code{0.15} for a
#'   15% increase.
#' @param nyrs Number of most recent years of \code{fit$report$spr} used to
#'   calculate mean current SPR. Default is 1.
#' @param n1 Number of recent years used by \code{r_rule()} for the numerator
#'   mean. Default is 2.
#' @param n2 Number of preceding years used by \code{r_rule()} for the
#'   denominator mean. Default is 3.
#'
#' @return A numeric multiplier. Values below 1 imply a reduction, 1 implies
#'   no change, and values above 1 imply an increase.
#'
#' @examples
#' # Default combined SPR-ratio rule
#' # hcr_sprr(fit, gear = 1)
#'
#' # Custom SPR control points
#' # hcr_sprlbi(
#' #   fit, gear = 1,
#' #   b1 = 0.10, b2 = 0.35, b3 = 0.50, b4 = 0.75,
#' #   dlow = -0.20, dopt = 0, dhi = 0.15
#' # )
#'
#' # Use a three-year mean SPR and alternative trend windows
#' # hcr_sprr(fit, gear = 1, nyrs = 3, n1 = 3, n2 = 5)
#'
#' @export
hcr_sprlbi <- function(fit, gear,
                     thresh = 0.75,
                     b1 = 0.10,
                     b2 = 0.35,
                     b3 = 0.50,
                     b4 = 0.75,
                     dlow = -0.20,
                     dopt =  0.00,
                     dhi  =  0.15,
                     nyrs = 1,
                     n1 = 2,
                     n2 = 3) {

  idx <- LBIspr(fit, thresh = thresh)[[gear]]

  r <- r_rule(idx, n1 = n1, n2 = n2)

  spr <- mean(tail(fit$report$spr, nyrs))

  s <- spr_rule(
    spr,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b4 = b4,
    dlow = dlow,
    dopt = dopt,
    dhi = dhi
  )

  mult <- ifelse((r < 1 && s < 1) | (r > 1 && s > 1), s, 1)

  return(mult)
}




#' Call fiticc() inside the mp() function
#'
#' This function provides an interface to FLicc::fiticc() to be used inside
#' the mp() function of the mse package, directly analogous to sca.sa()/
#' sep.sa() (FLa4a's ADMB-based sca(), above) but for FLicc's TMB-based,
#' length-only assessment.
#'
#' Unlike sca(), fiticc() fits one FLR iteration at a time (its own `iter`
#' argument defaults to 1) rather than being vectorised across the iter
#' dimension -- so flicc.sa() runs its own loop over args$it and collects
#' the per-iteration fits, optionally in parallel via the same foreach/
#' doFuture/future idiom FLicc::mc_flicc() already uses for Monte Carlo
#' refits.
#'
#' @param stk The FLStock input object. flicc.sa() reads the sampled
#'   length-frequency data from attr(stk, "lfd"), and the per-gear catch
#'   split (with its own reporting error) from attr(stk, "catch_gear") --
#'   both attached by sampling.lfd.oem() -- not from stk itself, which is
#'   returned unchanged.
#' @param idx Unused; present for mpCtrl/mseCtrl dispatch compatibility.
#' @param args The mse per-cycle args used by mp() (ay, dy, it, ...).
#' @param tracking mp()'s tracking object.
#' @param lhpar An FLPar of life-history parameters (linf, k, M or Mk, L50,
#'   a, b, ...) passed to FLicc::stocklen() to build the length-bin/weight/
#'   maturity/M template. Static across cycles -- pass the OM's own lhpar.
#' @param sel_fun Character vector of per-gear selectivity function types
#'   ("logistic"/"normal"/"dsnormal"), matching gear_cfg$type per gear.
#' @param catch_by_gear Fallback only, used when attr(stk, "catch_gear") is
#'   absent (e.g. an OEM that doesn't attach it): a static numeric vector of
#'   per-gear catch/effort shares, or an FLQuants of shares over time.
#'   Normally left NULL -- the OEM's own sampled catch split is preferred,
#'   since it (unlike a fixed share) carries genuine observation error and
#'   varies by iteration and cycle.
#' @param m_model Natural-mortality-at-length model passed to stocklen():
#'   one of "constant", "inverse", "Lorenzen", "Gislason".
#' @param settings List of fiticc() model settings (pop_model, obs_model,
#'   ngtg, tail_cut, ...); see ?fiticc. A `tail_cut` (fraction of Linf) is
#'   applied in the likelihood and, by default, in LBIspr() and LBImean().
#' @param spr,thresh Target SPR (%) for LBIspr() and Ztgt, and the cumulative
#'   threshold defining LBIspr()'s Lref.
#' @param lbi Which length indicators to compute: "spr" (LBIspr(), share of
#'   catch above Lref relative to the SPR-target reference) and/or "mean"
#'   (LBImean(), mean length relative to the reference at F = M).
#' @param pool How gear indices are pooled into `LBIspr` / `LBImean`:
#'   `"none"` (default, only per-gear indices are returned), `"equal"`,
#'   `"catch"` (reported catch shares by year), or a named numeric vector of
#'   weights by gear (normalised). Gears differ in how informative their
#'   length compositions are (e.g. a gear selecting large fish responds more
#'   to changes in survival), so catch weighting is not necessarily the best
#'   choice for an indicator.
#' @param lbimean_args List of arguments passed to LBImean(): `ref`
#'   ("FM" = reference at F = FM x M, the ICES-type L_F=M; "spr" = reference
#'   at F_SPRx), `FM` (F/M ratio, default 1), `lc` ("sel50" = length at 50%
#'   gear selectivity, "none", or numeric by gear), and optionally `spr`
#'   (defaults to flicc.sa()'s `spr`), `nyears`, `scale_sel`.
#' @param ess.g Numeric scalar or vector (gear order) of effective sample
#'   sizes passed to lfdess(); each gear-year LFD is rescaled to this total.
#'   `NULL` keeps the sampled counts (e.g. the OEM's `ess_len`). A scalar
#'   gives every gear the same weight regardless of how many fish were
#'   measured.
#' @param drop_nonconv Logical; if TRUE, non-converged fits are not used
#'   (their indicators are NA, i.e. mult = 1 in lbi.hcr()).
#' @param spr_max Fits with any SPR >= `spr_max` are flagged as boundary
#'   solutions and not used; `NULL` disables the check.
#' @param refit_sel_fun Optional alternative `sel_fun` tried for boundary
#'   fits; the refit is kept only if it is no longer at the boundary.
#' @param debug_dir Optional folder where the inputs of boundary fits are
#'   saved (one .rds per assessment year) for offline diagnosis.
#' @param n_restart,grad_tol,compile,silent,dll Passed through to fiticc().
#' @param parallel,workers If TRUE, fit iterations in parallel via foreach/
#'   doFuture/future::multisession. Automatically disabled (with a warning)
#'   whenever outer_parallel = TRUE -- see outer_parallel below.
#' @param outer_parallel Set TRUE whenever mp() itself is called with
#'   parallel = TRUE. The two parallel loops (mp()'s own, and this
#'   function's) must never both be live at once -- a fresh
#'   future::multisession plan started inside every already-parallel mp()
#'   worker is exactly the nested-parallel-backend collision that caused an
#'   earlier browser()/worker crash in this project. Setting outer_parallel
#'   = TRUE forces parallel = FALSE here regardless of what was passed, so
#'   the two settings can never accidentally both be on. Default FALSE.
#' @param verbose If TRUE, show a progress bar (sequential) or `progressr`
#'   progress (parallel) and don't suppress fiticc()'s own console output.
#'   Default FALSE -- inside mp(), which already prints its own per-year
#'   progress, per-iteration chatter from n_restart optimizer runs is noise
#'   rather than useful signal.
#' @param nyrs If not NULL, only the trailing `nyrs` years ending at the
#'   current data year (args$dy) are used to fit -- both attr(stk, "lfd")
#'   and (when it is an FLQuants time series, i.e. the OEM's own sampled
#'   catch split) catch_by_gear are windowed to the same year range before
#'   fiticc() is called. Default NULL fits the full accumulated history in
#'   attr(stk, "lfd") (y0:dy), matching the previous behaviour.
#' @param ... Additional arguments passed to fiticc().
#'
#' @return A list with the (unchanged) stk, the indicators `ind` and the
#'   tracking object. `ind` is an FLQuants (year x iter; NA for unusable
#'   iterations) with
#'   \describe{
#'     \item{LBIspr.<gear>, LBIspr}{LBIspr() index by gear (= 1 at the SPR
#'       target), and the pooled index if `pool` is not "none".}
#'     \item{LBImean.<gear>, LBImean}{LBImean() index by gear (= 1 at the
#'       reference, default F = M), and the pooled index if `pool` is not
#'       "none".}
#'     \item{SPR}{Fitted SPR (absolute, 0-1).}
#'     \item{Z, Ztgt, Zrel}{Z = Fap + M, Ztgt = F_SPRx + M (apical F, fitted
#'       M), and Zrel = Z / Ztgt (> 1: fishing above target).}
#'   }
#'   Tracking adds conv.est, sprbound.est, refit.est, used.est, the LBIspr
#'   cut-off Lref.<gear> and the LBImean reference Lmeanref.<gear>, so that
#'   drift of the references between refits can be monitored.
#'
#' @name flicc.sa
#' @rdname flicc.sa
#' @keywords classes
#' @export
flicc.sa <- function(stk, idx = NULL, args, tracking,
                     lhpar, sel_fun, catch_by_gear = NULL,
                     m_model = "constant",
                     settings = list(
                       CVL = 0.1, GL = 30, catch.sd = 0.05,
                       pop_model = "gtg", obs_model = "dm",
                       ngtg = 13, maxsd = 2, Mpow = 0
                     ),
                     spr = 40, thresh = 0.75, ess.g = 150,
                     lbi = c("spr", "mean"),  # LBIspr() and/or LBImean()
                     lbimean_args = list(ref = "FM", FM = 1, lc = "sel50"),
                     pool = "none",           # "none" | "equal" | "catch" | named weights
                     nyrs = NULL,
                     drop_nonconv = TRUE,     # non-converged fits -> mult = 1
                     spr_max = 0.99,          # SPR >= spr_max flagged as boundary; NULL = off
                     refit_sel_fun = NULL,    # e.g. c("logistic","logistic"); NULL = no refit
                     debug_dir = NULL,        # folder to save flagged-fit inputs; NULL = off
                     n_restart = 3, grad_tol = 1e-3,
                     compile = FALSE, silent = TRUE, dll = "FLicc",
                     parallel = FALSE, workers = NULL,
                     outer_parallel = FALSE, verbose = FALSE, ...) {

  # --- LFD: sum areas, downweight -----------------------------------------

  #browser()
  lfd <- attr(stk, "lfd")
  if (is.null(lfd))
    stop("attr(stk, 'lfd') is NULL -- flicc.sa() must be paired with a ",
         "length-sampling OEM (e.g. sampling.lfd.oem()) that attaches lfd to stk.")
  lfd <- FLQuants(lapply(lfd, unitSums))
  if (!is.null(ess.g)) lfd <- FLicc::lfdess(lfd, ess.g = ess.g)

  # --- Catch split: OEM-sampled shares preferred over static fallback ------
  cg_attr <- attr(stk, "catch_gear")
  if (!is.null(cg_attr) && length(cg_attr) > 0) {
    cg_attr  <- FLQuants(lapply(cg_attr, unitSums))
    cg_total <- Reduce(`+`, cg_attr)
    catch_by_gear <- FLQuants(lapply(cg_attr, function(x) {
      out <- x / cg_total
      out[!is.finite(out)] <- 0
      out
    }))
  } else if (is.null(catch_by_gear)) {
    stop("Neither attr(stk, 'catch_gear') nor a fallback catch_by_gear was ",
         "supplied -- flicc.sa() needs one or the other.")
  }

  # outer_parallel always wins: never let both loops be live at once
  if (isTRUE(outer_parallel) && isTRUE(parallel)) {
    warning("flicc.sa(): parallel = TRUE ignored because outer_parallel = TRUE ",
            "(mp() itself is running parallel) -- nested parallel backends are unsafe.")
    parallel <- FALSE
  }

  gears <- names(lfd)
  it <- args$it
  ay <- args$ay
  dy <- args$dy

  # --- Optional rolling window --------------------------------------------
  if (!is.null(nyrs)) {
    fit_years <- ac(seq(as.numeric(dy) - nyrs + 1, as.numeric(dy)))
    lfd <- FLQuants(setNames(lapply(lfd, function(x) {
      yrs <- intersect(fit_years, dimnames(x)$year)
      if (length(yrs) == 0)
        stop("flicc.sa(): nyrs = ", nyrs, " leaves no overlapping years in ",
             "attr(stk, 'lfd') at ay = ", ay, " (dy = ", dy, ").")
      x[, yrs]
    }), gears))
    if (inherits(catch_by_gear, "FLQuants")) {
      catch_by_gear <- FLQuants(setNames(lapply(catch_by_gear, function(x) {
        x[, intersect(fit_years, dimnames(x)$year)]
      }), names(catch_by_gear)))
    }
  }

  # --- Fit one iteration, with diagnostics and optional refit -------------
  fit_one <- function(i) {
    lfd_i    <- iter(FLQuants(lfd), i)
    stklen_i <- stocklen(lfd_i, lhpar = lhpar, m_model = m_model)
    cbg_i    <- iter(catch_by_gear, i)

    do_fit <- function(sf) {
      try(fiticc(
        lfd = lfd_i, stklen = stklen_i, sel_fun = sf,
        catch_by_gear = cbg_i, settings = settings,
        iter = 1, compile = compile, silent = !verbose, dll = dll,
        n_restart = n_restart, grad_tol = grad_tol, ...
      ), silent = TRUE)
    }

    check <- function(fit) {
      if (inherits(fit, "try-error"))
        return(list(ok = FALSE, conv = NA_real_, bound = NA))
      conv  <- as.numeric(!need_refit_flicc(fit))
      bound <- if (is.null(spr_max)) FALSE else
        any(c(fit$report$spr) >= spr_max, na.rm = TRUE)
      list(ok = TRUE, conv = conv, bound = bound)
    }

    fit   <- do_fit(sel_fun)
    chk   <- check(fit)
    refit <- FALSE

    # boundary solution: try alternative selectivity, keep only if it helps
    if (!is.null(refit_sel_fun) && chk$ok && isTRUE(chk$bound)) {
      fit2 <- do_fit(refit_sel_fun)
      chk2 <- check(fit2)
      if (chk2$ok && !isTRUE(chk2$bound)) {
        fit <- fit2; chk <- chk2; refit <- TRUE
      }
    }

    list(i = i, ok = chk$ok, fit = if (chk$ok) fit else NULL,
         conv = chk$conv, bound = chk$bound, refit = refit)
  }

  # --- Fit across iterations ----------------------------------------------
  if (parallel) {
    requireNamespace("foreach"); requireNamespace("doFuture"); requireNamespace("future")
    doFuture::registerDoFuture()
    if (!is.null(workers)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = workers)
    }
    run <- function(i) {
      if (verbose) fit_one(i) else
        suppressWarnings(suppressPackageStartupMessages(suppressMessages(fit_one(i))))
    }
    res <- foreach::foreach(
      i = seq_len(it), .options.future = list(seed = TRUE)
    ) %dofuture% { run(i) }
  } else {
    res <- lapply(seq_len(it), fit_one)
  }

  # --- Classify iterations ------------------------------------------------
  ok    <- vapply(res, `[[`, logical(1), "ok")
  conv  <- vapply(res, function(r) as.numeric(r$conv),  numeric(1))
  bound <- vapply(res, function(r) as.numeric(r$bound), numeric(1))
  refit <- vapply(res, function(r) as.numeric(r$refit), numeric(1))

  use <- ok & bound %in% 0
  if (isTRUE(drop_nonconv)) use <- use & conv %in% 1

  # save inputs of boundary fits for offline diagnosis
  if (!is.null(debug_dir) && any(bound %in% 1)) {
    dir.create(debug_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(ay = ay, dy = dy, iters = which(bound %in% 1),
                 lfd = lfd, catch_by_gear = catch_by_gear,
                 lhpar = lhpar, sel_fun = sel_fun, m_model = m_model,
                 settings = settings),
            file.path(debug_dir, paste0("flicc_bound_ay", ay, ".rds")))
  }

  if (!any(use))
    warning("flicc.sa(): no usable fits at ay = ", ay,
            " -- all iterations get NA indicators (mult = 1 in lbi.hcr)")

  # --- Indicators on usable iterations -----------------------------------
  lbi <- match.arg(lbi, several.ok = TRUE)

  ## catch weights by gear for pooled indices: from the OEM's reported catch
  ## (yearly shares), or the catch_by_gear fallback -- FLQuants (shares or
  ## catches, by year) or a static numeric vector; always normalised
  gear_w <- function(i) {
    if (identical(pool, "equal"))
      return(as.list(setNames(rep(1 / length(gears), length(gears)), gears)))
    if (is.numeric(pool)) {
      if (is.null(names(pool)) || !all(gears %in% names(pool)))
        stop("flicc.sa(): numeric 'pool' must be named by gear: ", paste(gears, collapse = ", "))
      return(as.list(pool[gears] / sum(pool[gears])))
    }
    if (inherits(catch_by_gear, "FLQuants")) {
      ## shares or catches: normalise to sum to 1 per year
      w   <- lapply(setNames(nm = gears), function(g) iter(catch_by_gear[[g]], i))
      tot <- Reduce(`+`, w)
      lapply(w, function(x) { s <- x / tot; s[!is.finite(s)] <- 0; s })
    } else {
      w <- catch_by_gear
      if (is.null(names(w))) names(w) <- gears
      as.list(w[gears] / sum(w[gears]))
    }
  }

  make_ind <- function(fit, i) {
    out <- list()
    w <- if (identical(pool, "none")) NULL else gear_w(i)
    if ("spr" %in% lbi) {
      li <- LBIspr(fit, gear = gears, spr = spr, thresh = thresh)
      for (g in gears) out[[paste0("LBIspr.", g)]] <- index(li[[g]])
      if (!identical(pool, "none"))
        out$LBIspr <- Reduce(`+`, lapply(gears, function(g)
          out[[paste0("LBIspr.", g)]] * w[[g]]))
      attr(out, "Lref") <- attr(li, "Lref")
    }
    if ("mean" %in% lbi) {
      la <- utils::modifyList(list(spr = spr), lbimean_args)   # flicc.sa's spr unless set
      lm <- do.call(LBImean, c(list(fit = fit, gear = gears), la))
      for (g in gears) out[[paste0("LBImean.", g)]] <- index(lm[[g]])
      if (!identical(pool, "none"))
        out$LBImean <- Reduce(`+`, lapply(gears, function(g)
          out[[paste0("LBImean.", g)]] * w[[g]]))
      attr(out, "Lmean_ref") <- attr(lm, "Lmean_ref")
    }
    out$SPR  <- fit$report$spr
    M        <- c(fit$report$lhpar["M"])                 # fitted M (Mk * k)
    out$Z    <- fit$report$Fap + M                       # apical F + M
    out$Ztgt <- out$Z
    out$Ztgt[] <- fspr_flicc(fit, spr = spr) + M         # F_SPRx + M
    out$Zrel <- out$Z / out$Ztgt
    out
  }

  ind_i <- vector("list", it)
  ind_i[use] <- lapply(which(use), function(i) make_ind(res[[i]]$fit, i))

  ## NA templates so unusable iterations never inherit another's values
  if (any(use)) {
    first <- ind_i[[which(use)[1]]]
    ind <- lapply(first, function(x) { q <- propagate(x, it); q[] <- NA; q })
    for (i in which(use))
      for (nm in names(ind)) iter(ind[[nm]], i) <- ind_i[[i]][[nm]]
  } else {
    pooled <- !identical(pool, "none")
    nms <- c(if ("spr" %in% lbi) c(paste0("LBIspr.", gears), if (pooled) "LBIspr"),
             if ("mean" %in% lbi) c(paste0("LBImean.", gears), if (pooled) "LBImean"),
             "SPR", "Z", "Ztgt", "Zrel")
    na_q <- FLQuant(NA, dimnames = list(year = dimnames(lfd[[1]])$year,
                                        iter = seq_len(it)))
    ind <- setNames(rep(list(na_q), length(nms)), nms)
  }
  ind <- FLQuants(ind)

  # --- Tracking -----------------------------------------------------------
  track(tracking, "conv.est",     ac(ay)) <- conv           # 1 / 0 / NA (error)
  track(tracking, "sprbound.est", ac(ay)) <- bound          # 1 = SPR at boundary
  track(tracking, "refit.est",    ac(ay)) <- refit          # 1 = rescued by refit
  track(tracking, "used.est",     ac(ay)) <- as.numeric(use) # 0 -> mult = 1
  if ("spr" %in% lbi) {
    for (g in gears) {
      lr <- rep(NA_real_, it)
      lr[use] <- vapply(ind_i[use], function(x) unname(attr(x, "Lref")[g]), numeric(1))
      track(tracking, paste0("Lref.", g), ac(ay)) <- lr   # LBIspr cut-off drift
    }
  }
  if ("mean" %in% lbi) {
    for (g in gears) {
      lm_ref <- rep(NA_real_, it)
      lm_ref[use] <- vapply(ind_i[use], function(x)
        unname(attr(x, "Lmean_ref")[g]), numeric(1))
      track(tracking, paste0("Lmeanref.", g), ac(ay)) <- lm_ref  # reference drift
    }
  }

  list(stk = stk, ind = ind, tracking = tracking)
}

#' mp()-compatible trend x status harvest control rule
#'
#' Ports hcr_sprlbi()'s own combined logic -- a recent-trend ratio on a
#' length indicator (r_rule()'s "recent n1 years vs prior n2 years" ratio),
#' combined with an SPR-based status check (spr_rule()), applied only when
#' both point the same direction -- onto flicc.sa()'s own `ind` output
#' (FLQuants: one element per gear plus "SPR"), returning a *relative*
#' multiplicative adjustment to the previous catch/TAC rather than an
#' absolute output level. Broadly the shape used by trend-driven empirical
#' MPs (e.g. CCSBT's own Bali Procedure): NOT a reproduction of any
#' particular stock's tuned gain/trigger constants -- b1-b4/dlow/dhi below
#' are spr_rule()'s own defaults and need calibrating against this OM.
#'
#' @param gear Name of the `ind` element (from [flicc.sa()]) supplying the
#'   trend ratio r: `"SPR"` (default), a gear's length indicator
#'   `"LBIspr.<gear>"` or `"LBImean.<gear>"` (e.g. `"LBIspr.Trawl"`,
#'   `"LBIspr.Gillnet"`), or a pooled `"LBIspr"` / `"LBImean"` (only when
#'   flicc.sa() is run with `pool` other than "none"). Plain gear names
#'   (e.g. "Trawl") are no longer valid. A numeric index is also accepted.
#' @param output "catch" (TAC, the default) or "effort" -- whichever
#'   fwd.om's projection method expects a relative-to-previous adjustment
#'   applied to. NOT "fbar" -- this rule has no F to set.
#' @param spr_trigger SPR below which `decision.hcr` is set to 3, so that
#'   flicc.is() applies its `dtaclow`/`dtacupp`/`Cmax` limits. Default 0.1;
#'   0 disables the limits.
#' @param ... Passed through, unused (dispatch compatibility).
#' @examples
#' \dontrun{
#' hcr <- mseCtrl(method = lbi.hcr, args = list(gear = "LBIspr.Gillnet"))
#' hcr <- mseCtrl(method = lbi.hcr, args = list(gear = "LBImean.Trawl"))
#' }
#' @export
lbi.hcr <- function(stk, ind, args, tracking, gear = "SPR",
                    b1 = 0.10, b2 = 0.3, b3 = 0.50, b4 = 0.75,
                    dlow = -0.20, dopt = 0.00, dhi = 0.15,
                    nyrs = 1, n1 = 2, n2 = 3, output = "catch",
                    spr_trigger = 0.1, ...) {

  FLCore::spread(args)

  # trend ratio -- r_rule()'s own logic, applied to the FLQuant ind[[gear]]
  # (flicc.sa() already unwraps LBIspr()/LBImean() from FLIndexBiomass)
  if (is.character(gear) && !gear %in% names(ind))
    stop("lbi.hcr(): '", gear, "' not in ind (",
         paste(names(ind), collapse = ", "), "). Gear indicators are named ",
         "e.g. 'LBIspr.Trawl' or 'LBImean.Trawl'.")
  idx <- ind[[gear]]
  r <- c(yearMeans(tail(idx, n1)) / yearMeans(head(tail(idx, n1 + n2), n2)))
  #browser()
  # SPR-based status rule, unchanged from spr_rule()
  spr_cur <- c(yearMeans(tail(ind[["SPR"]], nyrs)))
  s <- spr_rule(spr_cur, b1 = b1, b2 = b2, b3 = b3, b4 = b4,
                dlow = dlow, dopt = dopt, dhi = dhi)


  # apply s only when r and s agree on direction, exactly as hcr_sprlbi()
  mult <- ifelse((r < 1 & s < 1) | (r > 1 & s > 1), s, 1)
  #browser()
  # iterations without a usable fit (NA indicators from flicc.sa()) keep
  # the previous TAC: mult = 1
  mult[!is.finite(mult)] <- 1

  # decision.hcr trigger flag, read by flicc.is() the same way tacspm.is():
  # 3 = TAC-change limits apply (SPR below spr_trigger), 1 = otherwise
  decision <- ifelse(is.finite(spr_cur) & spr_cur < spr_trigger, 3, 1)

  track(tracking, "trend.hcr", ac(ay)) <- FLQuant(c(r), dimnames = list(iter = seq_along(decision)))
  track(tracking, "metric.hcr", ac(ay)) <- ind[["SPR"]][, ac(dy)]
  track(tracking, "decision.hcr", ac(ay)) <-
    FLQuant(decision, dimnames = list(iter = seq_along(decision)))
  track(tracking, "mult.hcr", ac(ay)) <-
    FLQuant(c(mult), dimnames = list(iter = seq_along(mult)))

  #browser()
  ctrl <- fwdControl(list(year = ay, quant = output, value = c(mult)))

  list(ctrl = ctrl, tracking = tracking)
}


#' mp()-compatible trend x status harvest control rule
#'
#' Ports hcr_sprlbi()'s own combined logic -- a recent-trend ratio on a
#' length indicator (r_rule()'s "recent n1 years vs prior n2 years" ratio),
#' combined with an SPR-based status check (spr_rule()), applied only when
#' both point the same direction -- onto flicc.sa()'s own `ind` output
#' (FLQuants: one element per gear plus "SPR"), returning a *relative*
#' multiplicative adjustment to the previous catch/TAC rather than an
#' absolute output level. Broadly the shape used by trend-driven empirical
#' MPs (e.g. CCSBT's own Bali Procedure): NOT a reproduction of any
#' particular stock's tuned gain/trigger constants -- b1-b4/dlow/dhi below
#' are spr_rule()'s own defaults and need calibrating against this OM.
#'
#' @param gear Name of the `ind` element (from [flicc.sa()]) supplying the
#'   trend ratio r: `"SPR"` (default), a gear's length indicator
#'   `"LBIspr.<gear>"` or `"LBImean.<gear>"` (e.g. `"LBIspr.Trawl"`,
#'   `"LBIspr.Gillnet"`), or a pooled `"LBIspr"` / `"LBImean"` (only when
#'   flicc.sa() is run with `pool` other than "none"). Plain gear names
#'   (e.g. "Trawl") are no longer valid. A numeric index is also accepted.
#' @param output "catch" (TAC, the default) or "effort" -- whichever
#'   fwd.om's projection method expects a relative-to-previous adjustment
#'   applied to. NOT "fbar" -- this rule has no F to set.
#' @param spr_trigger SPR below which `decision.hcr` is set to 3, so that
#'   flicc.is() applies its `dtaclow`/`dtacupp`/`Cmax` limits. Default 0.1;
#'   0 disables the limits.
#' @param ... Passed through, unused (dispatch compatibility).
#' @examples
#' \dontrun{
#' hcr <- mseCtrl(method = lbi.hcr, args = list(gear = "LBIspr.Gillnet"))
#' hcr <- mseCtrl(method = lbi.hcr, args = list(gear = "LBImean.Trawl"))
#' }
#' @export
lbi.hcr <- function(stk, ind, args, tracking, gear = "SPR",
                    b1 = 0.10, b2 = 0.3, b3 = 0.50, b4 = 0.75,
                    dlow = -0.20, dopt = 0.00, dhi = 0.15,
                    nyrs = 1, n1 = 2, n2 = 3, output = "catch",
                    spr_trigger = 0.1, ...) {

  FLCore::spread(args)

  # trend ratio -- r_rule()'s own logic, applied to the FLQuant ind[[gear]]
  # (flicc.sa() already unwraps LBIspr()/LBImean() from FLIndexBiomass)
  if (is.character(gear) && !gear %in% names(ind))
    stop("lbi.hcr(): '", gear, "' not in ind (",
         paste(names(ind), collapse = ", "), "). Gear indicators are named ",
         "e.g. 'LBIspr.Trawl' or 'LBImean.Trawl'.")
  idx <- ind[[gear]]
  r <- c(yearMeans(tail(idx, n1)) / yearMeans(head(tail(idx, n1 + n2), n2)))
  #browser()
  # SPR-based status rule, unchanged from spr_rule()
  spr_cur <- c(yearMeans(tail(ind[["SPR"]], nyrs)))
  s <- spr_rule(spr_cur, b1 = b1, b2 = b2, b3 = b3, b4 = b4,
                dlow = dlow, dopt = dopt, dhi = dhi)


  # apply s only when r and s agree on direction, exactly as hcr_sprlbi()
  mult <- ifelse((r < 1 & s < 1) | (r > 1 & s > 1), s, 1)
  #browser()
  # iterations without a usable fit (NA indicators from flicc.sa()) keep
  # the previous TAC: mult = 1
  mult[!is.finite(mult)] <- 1

  # decision.hcr trigger flag, read by flicc.is() the same way tacspm.is():
  # 3 = TAC-change limits apply (SPR below spr_trigger), 1 = otherwise
  decision <- ifelse(is.finite(spr_cur) & spr_cur < spr_trigger, 3, 1)

  track(tracking, "trend.hcr", ac(ay)) <- FLQuant(c(r), dimnames = list(iter = seq_along(decision)))
  track(tracking, "metric.hcr", ac(ay)) <- ind[["SPR"]][, ac(dy)]
  track(tracking, "decision.hcr", ac(ay)) <-
    FLQuant(decision, dimnames = list(iter = seq_along(decision)))
  track(tracking, "mult.hcr", ac(ay)) <-
    FLQuant(c(mult), dimnames = list(iter = seq_along(mult)))

  #browser()
  ctrl <- fwdControl(list(year = ay, quant = output, value = c(mult)))

  list(ctrl = ctrl, tracking = tracking)
}

#' rfb-type harvest control rule for flicc.sa() output
#'
#' ICES category-3 style rule adapted to the indicators returned by
#' [flicc.sa()] (`LBIspr`, `LBImean` pooled and by gear, `SPR`, `Zrel`):
#' \deqn{A_{y+1} = A_y \, r^{\gamma} \, f \, b \, m}
#' returned as a TAC multiplier for [flicc.is()].
#'
#' \describe{
#'   \item{r}{Trend ratio of the `index` series: mean of the last `n1`
#'     years over the mean of the preceding `n2` years (ICES 2-over-3).
#'     `index` is any element of `ind`, e.g. `"LBIspr.Gillnet"` (default),
#'     `"LBImean.Trawl"`, a pooled `"LBIspr"`/`"LBImean"` (if flicc.sa() was
#'     run with `pool`), or `"SPR"`. A trend cancels a
#'     constant scale bias, so it is the component least sensitive to
#'     misspecified M, Linf or selectivity.}
#'   \item{f}{Fishing-pressure component. `NULL` (default) sets f = 1 (an
#'     rb rule). Otherwise an element of `ind` whose *level* is used, e.g.
#'     `"LBImean"` (= 1 at F = M, the ICES-type f), `"LBIspr"` (= 1 at the
#'     SPR target) or `"Zrel"` (used as 1/Zrel). Levels inherit the model's
#'     scale bias.}
#'   \item{b}{Safeguard \eqn{\min(1, SPR_y / SPR_{trigger})}. A low-biased
#'     SPR estimate makes it more precautionary.}
#'   \item{m}{Precautionary multiplier (ICES rfb: 0.95 for k < 0.2).}
#' }
#'
#' Iterations without usable indicators (NA from flicc.sa()) keep the
#' previous TAC (mult = 1). `decision.hcr` is 3 when b = 1, so that
#' flicc.is() applies its stability limits (`dtacupp`, `dtaclow`) only when
#' the safeguard is not active, as in the ICES rfb rule.
#'
#' @param stk,ind,args,tracking Standard mp() arguments; `ind` from flicc.sa().
#' @param index Element of `ind` for the trend r.
#' @param f Level component, taken in the latest year: `NULL` (f = 1, rb
#'   rule); `"LBImean"` or `"LBImean.<gear>"` (Lmean / L_F=M, the ICES-type
#'   f); `"Zrel"` (used as 1 / Zrel = Ztgt / Z); or `"LBIspr"` /
#'   `"LBIspr.<gear>"`. In all cases f < 1 indicates fishing above the
#'   reference.
#' @param n1,n2 Recent and reference windows of the trend ratio.
#' @param gamma Exponent on r (1 = ICES; < 1 damps the response).
#' @param spr_trigger SPR below which b < 1.
#' @param nspr Years of SPR averaged for b.
#' @param m Precautionary multiplier.
#' @param output "catch" (default) or "effort".
#' @param ... Unused.
#' @return List with `ctrl` (fwdControl holding the TAC multiplier) and
#'   `tracking` (r.hcr, f.hcr, b.hcr, mult.hcr, decision.hcr).
#' @examples
#' \dontrun{
#' ctrl <- mpCtrl(list(
#'   est  = mseCtrl(method = flicc.sa, args = est_args),
#'   hcr  = mseCtrl(method = rfb.flicc.hcr,
#'                  args = list(index = "LBIspr.Gillnet", spr_trigger = 0.25)),
#'   isys = mseCtrl(method = flicc.is,
#'                  args = list(initac = initac, dtacupp = 1.2, dtaclow = 0.7))
#' ))
#' }
#' @export
rfb.flicc.hcr <- function(stk, ind, args, tracking,
                          index = "LBIspr.Gillnet", f = NULL, n1 = 2, n2 = 3, gamma = 1,
                          spr_trigger = 0.25, nspr = 1, m = 0.95,
                          output = "catch", ...) {

  FLCore::spread(args)

  get_series <- function(nm) {
    if (!nm %in% names(ind))
      stop("rfb.flicc.hcr(): '", nm, "' not in ind (",
           paste(names(ind), collapse = ", "), ")")
    ind[[nm]]
  }

  ## r: trend ratio (n1 over n2)
  I <- get_series(index)
  r <- c(yearMeans(tail(I, n1)) / yearMeans(head(tail(I, n1 + n2), n2)))^gamma

  ## f: optional level component
  fv <- if (is.null(f)) rep(1, length(r)) else c(tail(get_series(f), 1))
  if (identical(f, "Zrel")) fv <- 1 / fv              # Z above target -> f < 1

  ## b: SPR safeguard
  spr_cur <- c(yearMeans(tail(ind[["SPR"]], nspr)))
  b <- pmin(1, spr_cur / spr_trigger)

  mult <- r * fv * b * m
  mult[!is.finite(mult)] <- 1                         # no usable fit: keep TAC
  b[!is.finite(b)] <- 1

  ## stability limits in flicc.is() only when the safeguard is inactive
  decision <- ifelse(b >= 1, 3, 1)

  its <- seq_along(mult)
  track(tracking, "r.hcr", ac(ay))        <- FLQuant(r,        dimnames = list(iter = its))
  track(tracking, "f.hcr", ac(ay))        <- FLQuant(fv,       dimnames = list(iter = its))
  track(tracking, "b.hcr", ac(ay))        <- FLQuant(b,        dimnames = list(iter = its))
  track(tracking, "mult.hcr", ac(ay))     <- FLQuant(mult,     dimnames = list(iter = its))
  track(tracking, "decision.hcr", ac(ay)) <- FLQuant(decision, dimnames = list(iter = its))

  ctrl <- fwdControl(list(year = ay, quant = output, value = c(mult)))
  list(ctrl = ctrl, tracking = tracking)
}


#' mp()-compatible implementation system for lbi.hcr() and rfb.flicc.hcr
#'
#' Converts the TAC multiplier returned by lbi.hcr() into a TAC,
#' \eqn{TAC_y = TAC_{y-1} \times mult}, applies optional limits where
#' `decision.hcr > 2`, and optionally splits the TAC across units (areas).
#'
#' @param stk,ctrl,args,tracking Standard mp() arguments.
#' @param output Quantity controlled, default "catch".
#' @param dtaclow,dtacupp Lower/upper multipliers on the previous TAC.
#' @param Cmax Maximum TAC.
#' @param initac Previous TAC for the first cycle (ay == iy).
#' @param catch_area Optional TAC shares by unit, named by unit or in the
#'   order of `dimnames(stk)$unit` (excluding "combined").
#' @return List with `ctrl` and `tracking` (TAC tracked as "tac.is").
#' @export
flicc.is <- function(stk, ctrl, args, tracking, output = "catch",
                     dtaclow = NA, dtacupp = NA, Cmax = NA,
                     initac = NULL,catch_area=NULL) {

  FLCore::spread(args)

  # EXTRACT the multiplier lbi.hcr() put in ctrl (NOT an absolute TAC --
  # see the note in lbi.hcr() above)
  mult <- c(ctrl@iters[1, 2, ])

  # GET previous TAC -- same tracking convention as tacspm.is()
  if (ay == iy) {
    if (is.null(initac))
      stop("flicc.is(): initac must be supplied for the first mp() cycle (ay == iy).")
    prev_tac <- rep(initac, args$it)
  } else {
    prev_tac <- c(tracking[metric == "tac.is" & year == ay - frq, data])
  }

  #browser()
  TAC <- prev_tac * mult

  # ID iters where the hcr's decision crossed its own trigger -- unchanged
  # convention from tacspm.is()

  id <- tracking[metric  == "decision.hcr" & year == ay, data > 2]

  # APPLY upper/lower TAC limits relative to previous TAC, only for id
  # iters -- unchanged from tacspm.is()
  if (!is.na(dtacupp))
    TAC[id] <- pmin(TAC[id], prev_tac[id] * dtacupp)
  if (!is.na(dtaclow))
    TAC[id] <- pmax(TAC[id], prev_tac[id] * dtaclow)
  if (!is.na(Cmax))
    TAC[id] <- pmin(TAC[id], Cmax)

  ctrl <- fwdControl(list(year = ay + management_lag, quant = output, value = TAC))

  if(!is.null(catch_area)){
    ### get some dimensions
    units <- dimnames(stk)$unit
    units <- setdiff(units, "combined")
    names(units) <- units ### named vector needed for subsetting later

    ### create fwdControl object with units
    ctrl <- fwdControl(list(year = args$ay + args$management_lag,
                            quant = output,
                            unit = units,
                            value = NA))
    n_units <- length(units)
    n_its <- dims(stk)$iter
    iters_array <- array(NA, dim = c(n_units, 3, n_its),
                         dimnames = list(row = units,
                                         value = c("min", "value", "max"),
                                         iter = seq(n_its)))
    ### split catch advice into units based on smpf
    shares <- if (!is.null(names(catch_area))) catch_area[units] else catch_area
    for(i in seq_along(units)){
      iters_array[units[i], "value", ] <- TAC * shares[i]
    }

    ### intert into ctrl
    ctrl@iters <- iters_array
  }

  track(tracking, "tac.is", ac(ay)) <-
    FLQuant(c(TAC), dimnames = list(iter = seq_along(TAC)))

  list(ctrl = ctrl, tracking = tracking)
}

## flicc.is() reads mse's data.table tracking with data.table syntax
## (metric, year, data); make the FLicc namespace data.table-aware
.datatable.aware <- TRUE
utils::globalVariables(c("metric", "year", "data"))
