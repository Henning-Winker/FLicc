
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
#' @section Multiple areas (units):
#' When the OEM delivers length samples and catch by gear with areas in the
#' `unit` dimension (e.g. a spatial OM with areas as units), `unit_fit`
#' decides how they are assessed:
#' \describe{
#'   \item{"pooled"}{(default) compositions and catches are summed over units
#'     and one stock-wide model is fitted per iteration. Indicator names as
#'     below.}
#'   \item{"each"}{one model per unit and iteration, using only the gears
#'     with length samples in that unit (gears absent from a unit get NA
#'     indicators there). Every indicator and tracking metric gets the unit
#'     name as suffix, e.g. `SPR.a`, `LBIspr.trawl.a`, `used.est.a`.}
#'   \item{"both"}{the pooled fit and the per-unit fits.}
#' }
#' With a single unit, `unit_fit` is ignored. Per-unit fits multiply the
#' number of fits by the number of units with data.
#'
#' @param stk The FLStock input object. flicc.sa() reads the sampled
#'   length-frequency data from attr(stk, "lfd"), and the per-gear catch
#'   (with its own reporting error) from attr(stk, "catch_gear") --
#'   both attached by sampling.lfd.oem() -- not from stk itself, which is
#'   returned unchanged.
#' @param idx Unused; present for mpCtrl/mseCtrl dispatch compatibility.
#' @param args The mse per-cycle args used by mp() (ay, dy, it, ...).
#' @param tracking mp()'s tracking object.
#' @param lhpar An FLPar of life-history parameters (linf, k, M or Mk, L50,
#'   a, b, ...) passed to FLicc::stocklen() to build the length-bin/weight/
#'   maturity/M template. Static across cycles -- pass the OM's own lhpar.
#' @param sel_fun Character vector of per-gear selectivity function types
#'   ("logistic"/"normal"/"dsnormal"), named by gear or in the order of
#'   `names(attr(stk, "lfd"))`. Subset by gear for per-unit fits.
#' @param catch_by_gear Fallback only, used when attr(stk, "catch_gear") is
#'   absent (e.g. an OEM that doesn't attach it): a static numeric vector of
#'   per-gear catch/effort shares, or an FLQuants of shares over time.
#'   Normally left NULL -- the OEM's own sampled catch split is preferred,
#'   since it (unlike a fixed share) carries genuine observation error and
#'   varies by iteration and cycle.
#' @param m_model Natural-mortality-at-length model passed to stocklen():
#'   one of "constant", "inverse", "Lorenzen", "Gislason".
#' @param settings List of fiticc() model settings (pop_model, obs_model,
#'   ngtg, ...); see ?fiticc. `spr` is added to `settings$spr_ref`, so that
#'   Fspr at the target is reported by every fit.
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
#' @param ess.g Numeric scalar or vector (named by gear or in gear order) of
#'   effective sample sizes passed to lfdess(); each gear-year LFD is rescaled
#'   to this total. `NULL` keeps the sampled counts (e.g. the OEM's
#'   `ess_len`).
#' @param unit_fit `"pooled"` (default), `"each"` or `"both"`; see the
#'   section on multiple areas.
#' @param drop_nonconv Logical; if TRUE, non-converged fits are not used
#'   (their indicators are NA, i.e. mult = 1 in the HCR).
#' @param spr_max Fits with any SPR >= `spr_max` are flagged as boundary
#'   solutions and not used; `NULL` disables the check.
#' @param refit_sel_fun Optional alternative `sel_fun` (named by gear or in
#'   gear order) tried for boundary fits; the refit is kept only if it is no
#'   longer at the boundary.
#' @param debug_dir Optional folder where the inputs of boundary fits are
#'   saved (one .rds per assessment year and unit) for offline diagnosis.
#' @param n_restart,grad_tol,compile,silent,dll Passed through to fiticc().
#' @param parallel,workers If TRUE, fit iterations in parallel via foreach/
#'   doFuture/future::multisession. Automatically disabled (with a warning)
#'   whenever outer_parallel = TRUE -- see outer_parallel below.
#' @param outer_parallel Set TRUE whenever mp() itself is called with
#'   parallel = TRUE. The two parallel loops (mp()'s own, and this
#'   function's) must never both be live at once. Setting outer_parallel
#'   = TRUE forces parallel = FALSE here regardless of what was passed.
#'   Default FALSE.
#' @param verbose If TRUE, don't suppress fiticc()'s own console output.
#'   Default FALSE.
#' @param nyrs If not NULL, only the trailing `nyrs` years ending at the
#'   current data year (args$dy) are used to fit -- both attr(stk, "lfd")
#'   and (when it is an FLQuants time series) the catch split are windowed
#'   to the same year range before fiticc() is called. Default NULL fits the
#'   full accumulated history in attr(stk, "lfd") (y0:dy).
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
#'   For per-unit fits, the same names with `.<unit>` appended.
#'   Tracking adds conv.est, sprbound.est, refit.est, used.est, the LBIspr
#'   cut-off Lref.<gear> and the LBImean reference Lmeanref.<gear> (with
#'   `.<unit>` appended for per-unit fits), so that drift of the references
#'   between refits can be monitored.
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
                     unit_fit = c("pooled", "each", "both"),
                     nyrs = NULL,
                     drop_nonconv = TRUE,     # non-converged fits -> mult = 1
                     spr_max = 0.99,          # SPR >= spr_max flagged as boundary; NULL = off
                     refit_sel_fun = NULL,    # e.g. c("logistic","logistic"); NULL = no refit
                     debug_dir = NULL,        # folder to save flagged-fit inputs; NULL = off
                     n_restart = 3, grad_tol = 1e-3,
                     compile = FALSE, silent = TRUE, dll = "FLicc",
                     parallel = FALSE, workers = NULL,
                     outer_parallel = FALSE, verbose = FALSE, ...) {

  unit_fit <- match.arg(unit_fit)
  lbi      <- match.arg(lbi, several.ok = TRUE)
  pooled   <- !identical(pool, "none")

  # --- Inputs from the OEM -------------------------------------------------
  lfd_all <- attr(stk, "lfd")
  if (is.null(lfd_all))
    stop("attr(stk, 'lfd') is NULL -- flicc.sa() must be paired with a ",
         "length-sampling OEM (e.g. sampling.lfd.oem()) that attaches lfd to stk.")
  ## annual compositions: pool seasons if the OEM kept them by season
  annual <- function(x) if (dim(x)[4] > 1) seasonSums(x) else x
  lfd_all   <- FLQuants(lapply(lfd_all, annual))
  gears_all <- names(lfd_all)
  units_all <- setdiff(dimnames(lfd_all[[1]])$unit, "combined")
  if (length(units_all) <= 1) unit_fit <- "pooled"

  cg_all <- attr(stk, "catch_gear")
  if (is.null(cg_all) || length(cg_all) == 0) {
    cg_all <- NULL
    if (is.null(catch_by_gear))
      stop("Neither attr(stk, 'catch_gear') nor a fallback catch_by_gear was ",
           "supplied -- flicc.sa() needs one or the other.")
  } else {
    cg_all <- FLQuants(lapply(cg_all, annual))
  }

  ## Fspr at the target is reported by every fit (no per-recruit search in
  ## LBIspr() / fspr_flicc())
  settings$spr_ref <- unique(c(spr, settings$spr_ref))

  # outer_parallel always wins: never let both loops be live at once
  if (isTRUE(outer_parallel) && isTRUE(parallel)) {
    warning("flicc.sa(): parallel = TRUE ignored because outer_parallel = TRUE ",
            "(mp() itself is running parallel) -- nested parallel backends are unsafe.")
    parallel <- FALSE
  }

  it <- args$it
  ay <- args$ay
  dy <- args$dy

  ## per-gear settings, named by gear
  sel_fun <- .by_gear(sel_fun, gears_all, "sel_fun")
  if (!is.null(refit_sel_fun))
    refit_sel_fun <- .by_gear(refit_sel_fun, gears_all, "refit_sel_fun")
  if (!is.null(ess.g) && length(ess.g) > 1)
    ess.g <- .by_gear(ess.g, gears_all, "ess.g")

  if (parallel) {
    requireNamespace("foreach"); requireNamespace("doFuture"); requireNamespace("future")
    doFuture::registerDoFuture()
    if (!is.null(workers)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = workers)
    }
  }

  ## catch shares by gear and year (sum to one over gears)
  shares <- function(cg) {
    tot <- Reduce(`+`, cg)
    FLQuants(lapply(cg, function(x) {
      s <- x / tot
      s[!is.finite(s)] <- 0
      s
    }))
  }

  ## indicator names returned for every assessment, for all gears, so that
  ## names are the same in every unit (absent gears are NA)
  ind_names <- c(if ("spr" %in% lbi) c(paste0("LBIspr.", gears_all), if (pooled) "LBIspr"),
                 if ("mean" %in% lbi) c(paste0("LBImean.", gears_all), if (pooled) "LBImean"),
                 "SPR", "Z", "Ztgt", "Zrel")

  # --- One assessment (pooled or one unit) over all iterations -------------
  assess <- function(lfd, cbg, sfx, tracking) {

    nm  <- function(x) if (sfx == "") x else paste0(x, ".", sfx)
    lab <- if (sfx == "") "" else paste0(" (unit ", sfx, ")")

    ## optional rolling window
    if (!is.null(nyrs)) {
      fit_years <- ac(seq(as.numeric(dy) - nyrs + 1, as.numeric(dy)))
      lfd <- FLQuants(setNames(lapply(lfd, function(x) {
        yrs <- intersect(fit_years, dimnames(x)$year)
        if (length(yrs) == 0)
          stop("flicc.sa(): nyrs = ", nyrs, " leaves no overlapping years in ",
               "attr(stk, 'lfd') at ay = ", ay, " (dy = ", dy, ").")
        x[, yrs]
      }), names(lfd)))
      if (inherits(cbg, "FLQuants")) {
        cbg <- FLQuants(setNames(lapply(cbg, function(x) {
          x[, intersect(fit_years, dimnames(x)$year)]
        }), names(cbg)))
      }
    }

    ## gears with length samples here; in a unit, the others are not fished
    has_data <- vapply(lfd, function(x) sum(x, na.rm = TRUE) > 0, logical(1))
    gears    <- names(lfd)[has_data]
    yrs_out  <- dimnames(lfd[[1]])$year

    conv  <- rep(NA_real_, it)
    bound <- rep(NA_real_, it)
    refit <- rep(0, it)
    use   <- rep(FALSE, it)
    ind_i <- vector("list", it)

    if (length(gears) > 0) {

      lfd <- lfd[gears]
      if (!is.null(ess.g)) {
        eg  <- if (length(ess.g) > 1) unname(ess.g[gears]) else ess.g
        lfd <- FLicc::lfdess(lfd, ess.g = eg)
        lfd <- FLQuants(lapply(lfd, function(x) { x[!is.finite(x)] <- 0; x }))
      }

      ## catch split of the gears assessed here, renormalised
      if (inherits(cbg, "FLQuants")) {
        cbg <- shares(cbg[gears])
      } else {
        w <- cbg
        if (is.null(names(w))) names(w) <- gears_all
        cbg <- w[gears] / sum(w[gears])
      }

      sf  <- unname(sel_fun[gears])
      rsf <- if (is.null(refit_sel_fun)) NULL else unname(refit_sel_fun[gears])

      ## fit one iteration, with diagnostics and optional refit
      fit_one <- function(i) {
        lfd_i    <- iter(FLQuants(lfd), i)
        stklen_i <- stocklen(lfd_i, lhpar = lhpar, m_model = m_model)
        cbg_i    <- if (inherits(cbg, "FLQuants")) iter(cbg, i) else cbg

        do_fit <- function(s) {
          try(fiticc(
            lfd = lfd_i, stklen = stklen_i, sel_fun = s,
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

        fit   <- do_fit(sf)
        chk   <- check(fit)
        rf    <- FALSE

        ## boundary solution: try alternative selectivity, keep only if it helps
        if (!is.null(rsf) && chk$ok && isTRUE(chk$bound)) {
          fit2 <- do_fit(rsf)
          chk2 <- check(fit2)
          if (chk2$ok && !isTRUE(chk2$bound)) {
            fit <- fit2; chk <- chk2; rf <- TRUE
          }
        }

        list(i = i, ok = chk$ok, fit = if (chk$ok) fit else NULL,
             conv = chk$conv, bound = chk$bound, refit = rf)
      }

      if (parallel) {
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

      ## classify iterations
      ok    <- vapply(res, `[[`, logical(1), "ok")
      conv  <- vapply(res, function(r) as.numeric(r$conv),  numeric(1))
      bound <- vapply(res, function(r) as.numeric(r$bound), numeric(1))
      refit <- vapply(res, function(r) as.numeric(r$refit), numeric(1))

      use <- ok & bound %in% 0
      if (isTRUE(drop_nonconv)) use <- use & conv %in% 1

      ## save inputs of boundary fits for offline diagnosis
      if (!is.null(debug_dir) && any(bound %in% 1)) {
        dir.create(debug_dir, showWarnings = FALSE, recursive = TRUE)
        saveRDS(list(ay = ay, dy = dy, unit = sfx, iters = which(bound %in% 1),
                     lfd = lfd, catch_by_gear = cbg,
                     lhpar = lhpar, sel_fun = sf, m_model = m_model,
                     settings = settings),
                file.path(debug_dir, paste0("flicc_bound_ay", ay,
                                            if (sfx != "") paste0("_", sfx), ".rds")))
      }

      ## catch weights by gear for pooled indices: from the OEM's reported
      ## catch (yearly shares) or the catch_by_gear fallback; always normalised
      gear_w <- function(i) {
        if (identical(pool, "equal"))
          return(as.list(setNames(rep(1 / length(gears), length(gears)), gears)))
        if (is.numeric(pool)) {
          if (is.null(names(pool)) || !all(gears %in% names(pool)))
            stop("flicc.sa(): numeric 'pool' must be named by gear: ",
                 paste(gears_all, collapse = ", "))
          return(as.list(pool[gears] / sum(pool[gears])))
        }
        if (inherits(cbg, "FLQuants")) {
          w   <- lapply(setNames(nm = gears), function(g) iter(cbg[[g]], i))
          tot <- Reduce(`+`, w)
          lapply(w, function(x) { s <- x / tot; s[!is.finite(s)] <- 0; s })
        } else {
          as.list(cbg[gears] / sum(cbg[gears]))
        }
      }

      make_ind <- function(fit, i) {
        out <- list()
        w <- if (pooled) gear_w(i) else NULL
        if ("spr" %in% lbi) {
          li <- LBIspr(fit, gear = gears, spr = spr, thresh = thresh)
          for (g in gears) out[[paste0("LBIspr.", g)]] <- index(li[[g]])
          if (pooled)
            out$LBIspr <- Reduce(`+`, lapply(gears, function(g)
              out[[paste0("LBIspr.", g)]] * w[[g]]))
          attr(out, "Lref") <- attr(li, "Lref")
        }
        if ("mean" %in% lbi) {
          la <- utils::modifyList(list(spr = spr), lbimean_args)  # flicc.sa's spr unless set
          lm <- do.call(LBImean, c(list(fit = fit, gear = gears), la))
          for (g in gears) out[[paste0("LBImean.", g)]] <- index(lm[[g]])
          if (pooled)
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

      ind_i[use] <- lapply(which(use), function(i) make_ind(res[[i]]$fit, i))
    }

    ## (a unit without length samples is simply not assessed: no warning)
    if (!any(use) && length(gears) > 0)
      warning("flicc.sa(): no usable fits at ay = ", ay, lab,
              " -- all iterations get NA indicators (mult = 1 in the HCR)")

    ## NA templates so unusable iterations never inherit another's values;
    ## names absent from the fits (gears not fished here) stay NA
    if (any(use)) {
      first <- ind_i[[which(use)[1]]]
      tmpl  <- lapply(first, function(x) { q <- propagate(x, it); q[] <- NA; q })
      ind   <- setNames(lapply(ind_names, function(n)
        if (!is.null(tmpl[[n]])) tmpl[[n]] else tmpl$SPR), ind_names)
      for (i in which(use))
        for (n in names(ind_i[[i]])) iter(ind[[n]], i) <- ind_i[[i]][[n]]
    } else {
      na_q <- FLQuant(NA, dimnames = list(year = yrs_out, iter = seq_len(it)))
      ind  <- setNames(rep(list(na_q), length(ind_names)), ind_names)
    }
    names(ind) <- nm(names(ind))

    ## tracking
    track(tracking, nm("conv.est"),     ac(ay)) <- conv              # 1 / 0 / NA (error)
    track(tracking, nm("sprbound.est"), ac(ay)) <- bound             # 1 = SPR at boundary
    track(tracking, nm("refit.est"),    ac(ay)) <- refit             # 1 = rescued by refit
    track(tracking, nm("used.est"),     ac(ay)) <- as.numeric(use)   # 0 -> mult = 1
    if ("spr" %in% lbi) {
      for (g in gears) {
        lr <- rep(NA_real_, it)
        lr[use] <- vapply(ind_i[use], function(x) unname(attr(x, "Lref")[g]), numeric(1))
        track(tracking, nm(paste0("Lref.", g)), ac(ay)) <- lr        # LBIspr cut-off drift
      }
    }
    if ("mean" %in% lbi) {
      for (g in gears) {
        lm_ref <- rep(NA_real_, it)
        lm_ref[use] <- vapply(ind_i[use], function(x)
          unname(attr(x, "Lmean_ref")[g]), numeric(1))
        track(tracking, nm(paste0("Lmeanref.", g)), ac(ay)) <- lm_ref # reference drift
      }
    }

    list(ind = ind, tracking = tracking)
  }

  # --- Pooled and/or per-unit assessments ----------------------------------
  ind <- list()

  if (unit_fit %in% c("pooled", "both")) {
    lfd_p <- FLQuants(lapply(lfd_all, unitSums))
    cbg_p <- if (!is.null(cg_all)) shares(FLQuants(lapply(cg_all, unitSums))) else catch_by_gear
    out <- assess(lfd_p, cbg_p, "", tracking)
    ind <- c(ind, out$ind)
    tracking <- out$tracking
  }

  if (unit_fit %in% c("each", "both")) {
    for (u in units_all) {
      lfd_u <- FLQuants(lapply(lfd_all, function(x) x[, , u]))
      cbg_u <- if (!is.null(cg_all))
        shares(FLQuants(lapply(cg_all, function(x) x[, , u]))) else catch_by_gear
      out <- assess(lfd_u, cbg_u, u, tracking)
      ind <- c(ind, out$ind)
      tracking <- out$tracking
    }
  }

  list(stk = stk, ind = FLQuants(ind), tracking = tracking)
}

## per-gear argument as a vector named by gear: by name if all gears are
## named, otherwise in gear order (or recycled if a single value)
.by_gear <- function(x, gears, what) {
  if (!is.null(names(x)) && all(gears %in% names(x)))
    return(x[gears])
  if (length(x) == length(gears))
    return(setNames(unname(x), gears))
  if (length(x) == 1)
    return(setNames(rep(x, length(gears)), gears))
  stop("flicc.sa(): '", what, "' needs one value per gear (",
       paste(gears, collapse = ", "), "), named by gear or in gear order.")
}

## indicator / tracking name for a unit ("" = stock-wide)
.unit_name <- function(x, u) if (identical(u, "")) x else paste0(x, ".", u)

## ".<unit>" suffix for tracking metrics
.unit_sfx <- function(u) if (identical(u, "")) "" else paste0(".", u)

#' mp()-compatible trend x status harvest control rule
#'
#' Ports hcr_sprlbi()'s own combined logic -- a recent-trend ratio on a
#' length indicator (r_rule()'s "recent n1 years vs prior n2 years" ratio),
#' combined with an SPR-based status check (spr_rule()), applied only when
#' both point the same direction -- onto flicc.sa()'s own `ind` output,
#' returning a *relative* multiplicative adjustment to the previous
#' catch/TAC rather than an absolute output level. b1-b4/dlow/dhi need
#' calibrating against the OM.
#'
#' @section Multiple areas (units):
#' With `units` (e.g. `c("a", "b")`), the rule is applied separately to the
#' per-unit indicators of [flicc.sa()] run with `unit_fit = "each"` or
#' `"both"`: for unit `u` the trend uses `ind[[paste0(gear, ".", u)]]` and
#' the status `ind[[paste0("SPR.", u)]]`. The multipliers are tracked as
#' `mult.hcr.<unit>` (with `trend.hcr.<unit>`, `metric.hcr.<unit>`,
#' `decision.hcr.<unit>`) and read by [flicc.is()] with
#' `area_tac = "each"`. The returned `ctrl` then holds the mean of the unit
#' multipliers, for information only.
#'
#' @param stk,ind,args,tracking Standard mp() arguments; `ind` from flicc.sa().
#' @param gear Name of the `ind` element supplying the trend ratio r:
#'   `"SPR"` (default), a gear's length indicator `"LBIspr.<gear>"` or
#'   `"LBImean.<gear>"`, or a pooled `"LBIspr"` / `"LBImean"` (only when
#'   flicc.sa() is run with `pool` other than "none"). With `units`, the
#'   name without the unit suffix. A numeric index is also accepted when
#'   `units` is NULL.
#' @param b1,b2,b3,b4,dlow,dopt,dhi SPR rule, see [spr_rule()].
#' @param nyrs Years of SPR averaged for the status.
#' @param n1,n2 Recent and reference windows of the trend ratio.
#' @param output "catch" (TAC, the default) or "effort".
#' @param spr_trigger SPR below which `decision.hcr` is set to 3, so that
#'   flicc.is() applies its `dtaclow`/`dtacupp`/`Cmax` limits. Default 0.1;
#'   0 disables the limits.
#' @param units NULL (default, stock-wide rule) or the units (areas) to which
#'   the rule is applied separately.
#' @param ... Passed through, unused (dispatch compatibility).
#' @examples
#' \dontrun{
#' hcr <- mseCtrl(method = lbi.hcr, args = list(gear = "LBIspr.Gillnet"))
#' # per area, with flicc.sa(unit_fit = "each") and flicc.is(area_tac = "each")
#' hcr <- mseCtrl(method = lbi.hcr, args = list(gear = "SPR", units = c("a", "b")))
#' }
#' @export
lbi.hcr <- function(stk, ind, args, tracking, gear = "SPR",
                    b1 = 0.10, b2 = 0.3, b3 = 0.50, b4 = 0.75,
                    dlow = -0.20, dopt = 0.00, dhi = 0.15,
                    nyrs = 1, n1 = 2, n2 = 3, output = "catch",
                    spr_trigger = 0.1, units = NULL, ...) {

  FLCore::spread(args)

  us <- if (is.null(units)) "" else as.character(units)
  if (!is.null(units) && is.numeric(gear))
    stop("lbi.hcr(): use a character 'gear' with 'units'.")

  one <- function(u) {
    gnm <- if (is.numeric(gear)) gear else .unit_name(gear, u)
    snm <- .unit_name("SPR", u)
    for (n in list(gnm, snm)) {
      if (is.character(n) && !n %in% names(ind))
        stop("lbi.hcr(): '", n, "' not in ind (",
             paste(names(ind), collapse = ", "), "). Gear indicators are named ",
             "e.g. 'LBIspr.Trawl'; per-unit indicators need flicc.sa(unit_fit = 'each').")
    }
    # trend ratio -- r_rule()'s own logic, applied to the FLQuant ind[[gear]]
    idx <- ind[[gnm]]
    r <- c(yearMeans(tail(idx, n1)) / yearMeans(head(tail(idx, n1 + n2), n2)))
    # SPR-based status rule, unchanged from spr_rule()
    spr_cur <- c(yearMeans(tail(ind[[snm]], nyrs)))
    s <- spr_rule(spr_cur, b1 = b1, b2 = b2, b3 = b3, b4 = b4,
                  dlow = dlow, dopt = dopt, dhi = dhi)
    # apply s only when r and s agree on direction, exactly as hcr_sprlbi()
    mult <- ifelse((r < 1 & s < 1) | (r > 1 & s > 1), s, 1)
    # iterations without a usable fit (NA indicators) keep the previous TAC
    mult[!is.finite(mult)] <- 1
    # decision.hcr: 3 = TAC-change limits apply in flicc.is() (SPR below
    # spr_trigger), 1 = otherwise
    decision <- ifelse(is.finite(spr_cur) & spr_cur < spr_trigger, 3, 1)
    list(r = r, mult = mult, decision = decision, spr_dy = ind[[snm]][, ac(dy)])
  }

  res <- lapply(us, one)

  for (k in seq_along(us)) {
    sx  <- .unit_sfx(us[k])
    its <- seq_along(res[[k]]$mult)
    track(tracking, paste0("trend.hcr", sx), ac(ay)) <-
      FLQuant(c(res[[k]]$r), dimnames = list(iter = its))
    track(tracking, paste0("metric.hcr", sx), ac(ay)) <- res[[k]]$spr_dy
    track(tracking, paste0("decision.hcr", sx), ac(ay)) <-
      FLQuant(res[[k]]$decision, dimnames = list(iter = its))
    track(tracking, paste0("mult.hcr", sx), ac(ay)) <-
      FLQuant(c(res[[k]]$mult), dimnames = list(iter = its))
  }

  ## stock-wide: the multiplier; per unit: mean of the unit multipliers
  ## (information only -- flicc.is(area_tac = "each") reads mult.hcr.<unit>)
  mult <- Reduce(`+`, lapply(res, `[[`, "mult")) / length(res)
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
#'     run with `pool`), or `"SPR"`. A trend cancels a constant scale bias,
#'     so it is the component least sensitive to misspecified M, Linf or
#'     selectivity.}
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
#' @section Multiple areas (units):
#' With `units`, the rule is applied per unit to `paste0(index, ".", u)`,
#' `paste0(f, ".", u)` and `paste0("SPR.", u)` (flicc.sa() with
#' `unit_fit = "each"` or `"both"`), and the components are tracked with
#' the suffix `.<unit>` (`mult.hcr.<unit>` is read by
#' `flicc.is(area_tac = "each")`).
#'
#' @param stk,ind,args,tracking Standard mp() arguments; `ind` from flicc.sa().
#' @param index Element of `ind` for the trend r (without unit suffix when
#'   `units` is given).
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
#' @param units NULL (default, stock-wide rule) or the units (areas) to which
#'   the rule is applied separately.
#' @param ... Unused.
#' @return List with `ctrl` (fwdControl holding the TAC multiplier; the mean
#'   over units when `units` is given) and `tracking` (r.hcr, f.hcr, b.hcr,
#'   mult.hcr, decision.hcr, with `.<unit>` for per-unit rules).
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
                          output = "catch", units = NULL, ...) {

  FLCore::spread(args)

  us <- if (is.null(units)) "" else as.character(units)

  get_series <- function(nm) {
    if (!nm %in% names(ind))
      stop("rfb.flicc.hcr(): '", nm, "' not in ind (",
           paste(names(ind), collapse = ", "), ")")
    ind[[nm]]
  }

  one <- function(u) {
    ## r: trend ratio (n1 over n2)
    I <- get_series(.unit_name(index, u))
    r <- c(yearMeans(tail(I, n1)) / yearMeans(head(tail(I, n1 + n2), n2)))^gamma

    ## f: optional level component
    fv <- if (is.null(f)) rep(1, length(r)) else c(tail(get_series(.unit_name(f, u)), 1))
    if (identical(f, "Zrel")) fv <- 1 / fv              # Z above target -> f < 1

    ## b: SPR safeguard
    spr_cur <- c(yearMeans(tail(get_series(.unit_name("SPR", u)), nspr)))
    b <- pmin(1, spr_cur / spr_trigger)

    mult <- r * fv * b * m
    mult[!is.finite(mult)] <- 1                         # no usable fit: keep TAC
    b[!is.finite(b)] <- 1

    ## stability limits in flicc.is() only when the safeguard is inactive
    decision <- ifelse(b >= 1, 3, 1)
    list(r = r, f = fv, b = b, mult = mult, decision = decision)
  }

  res <- lapply(us, one)

  for (k in seq_along(us)) {
    sx  <- .unit_sfx(us[k])
    its <- seq_along(res[[k]]$mult)
    for (v in c("r", "f", "b", "mult", "decision"))
      track(tracking, paste0(v, ".hcr", sx), ac(ay)) <-
        FLQuant(res[[k]][[v]], dimnames = list(iter = its))
  }

  mult <- Reduce(`+`, lapply(res, `[[`, "mult")) / length(res)
  ctrl <- fwdControl(list(year = ay, quant = output, value = c(mult)))
  list(ctrl = ctrl, tracking = tracking)
}


#' mp()-compatible TAC implementation system for lbi.hcr() and rfb.flicc.hcr()
#'
#' Turns the relative change returned by the HCR into a TAC:
#' \deqn{TAC_y = TAC_{y-1} \times mult,}
#' applies optional limits to the change and a cap, sets the TAC for the
#' management years of the cycle, and optionally divides it among units
#' (areas).
#'
#' @details
#' Steps in each management cycle:
#' \enumerate{
#'   \item \strong{Multiplier.} Read from the HCR's `ctrl` (stock-wide) or,
#'     with `area_tac = "each"`, from the tracked `mult.hcr.<unit>`.
#'   \item \strong{Previous TAC.} `initac` in the first cycle (`ay == iy`),
#'     afterwards the TAC tracked as `tac.is` (or `tac.is.<unit>`) at
#'     `ay - frq`.
#'   \item \strong{New TAC.} Previous TAC times multiplier.
#'   \item \strong{Limits.} For iterations with `decision.hcr > 2` (or
#'     `decision.hcr.<unit> > 2`) the TAC is bounded to
#'     `[dtaclow, dtacupp]` times the previous TAC and capped at `Cmax`.
#'     `lbi.hcr()` flags iterations with SPR below its `spr_trigger`;
#'     `rfb.flicc.hcr()` flags iterations with an inactive safeguard (ICES
#'     rfb convention).
#'   \item \strong{Years.} The TAC is set for every management year of the
#'     cycle (`args$mys`, i.e. `ay + management_lag` to
#'     `ay + management_lag + frq - 1`), so it is held constant between
#'     assessments when `frq > 1`.
#'   \item \strong{Areas.} Optionally split among units (see below).
#' }
#'
#' @section Multiple areas (units):
#' \describe{
#'   \item{`area_tac = "split"`}{(default) one stock TAC from the stock-wide
#'     multiplier, limits on the stock TAC, then divided by `catch_area`:
#'     fixed shares (named by unit or in unit order), or `"recent"` for the
#'     observed catch shares by unit over the last `area_years` data years
#'     (from `attr(stk, "catch_gear")`, else `catch(stk)`), by iteration.
#'     Without `catch_area` the TAC is set for the stock as a whole.}
#'   \item{`area_tac = "each"`}{one TAC per unit, each with its own
#'     multiplier from an HCR run with `units` (falls back to the stock-wide
#'     multiplier for units without `mult.hcr.<unit>`) and its own previous
#'     TAC and limits. The first-cycle TACs are `initac` by unit (named
#'     vector) or `initac` times `catch_area`. `Cmax` caps the sum of the
#'     area TACs, scaling them down proportionally.}
#' }
#' Areas with a zero share or zero initial TAC stay closed (zero TAC) in
#' both modes. The control then has one row per management year and unit;
#' the OM projection method must take catch targets by unit.
#'
#' @param stk,ctrl,args,tracking Standard mp() arguments.
#' @param output Quantity controlled, default "catch".
#' @param dtaclow,dtacupp Lower/upper multipliers on the previous TAC.
#' @param Cmax Maximum TAC (stock total).
#' @param initac Previous TAC for the first cycle (ay == iy): a single
#'   stock total, or (with `area_tac = "each"`) a vector named by unit.
#' @param catch_area Optional TAC shares by unit, named by unit or in the
#'   order of `dimnames(stk)$unit` (excluding "combined"), or `"recent"`.
#' @param area_tac `"split"` (default) or `"each"`; see the section on
#'   multiple areas.
#' @param area_years Data years used for `catch_area = "recent"`.
#' @param ... Unused (dispatch compatibility).
#' @return List with `ctrl` and `tracking` (stock TAC tracked as "tac.is",
#'   area TACs as "tac.is.<unit>").
#' @examples
#' \dontrun{
#' # stock TAC, +/-20% limits when flagged by the HCR
#' isys <- mseCtrl(method = flicc.is,
#'                 args = list(initac = 1000, dtaclow = 0.8, dtacupp = 1.2))
#' # stock TAC split 30:70:0 among areas a, b, c (c closed)
#' isys <- mseCtrl(method = flicc.is,
#'                 args = list(initac = 1000, catch_area = c(a = 0.3, b = 0.7, c = 0)))
#' # area TACs from area-specific HCR multipliers
#' isys <- mseCtrl(method = flicc.is,
#'                 args = list(initac = c(a = 300, b = 700, c = 0), area_tac = "each"))
#' }
#' @export
flicc.is <- function(stk, ctrl, args, tracking, output = "catch",
                     dtaclow = NA, dtacupp = NA, Cmax = NA,
                     initac = NULL, catch_area = NULL,
                     area_tac = c("split", "each"), area_years = 3, ...) {

  FLCore::spread(args)
  area_tac <- match.arg(area_tac)

  it  <- args$it
  # management years of this cycle (mse >= 2.x passes args$mys)
  mys <- if (!is.null(args$mys)) args$mys else ay + management_lag

  # multiplier from the HCR (stock-wide)
  mult <- c(ctrl@iters[1, "value", ])

  # units (areas)
  units <- setdiff(dimnames(stk)$unit, "combined")
  names(units) <- units

  if (area_tac == "each" && length(units) < 2)
    stop("flicc.is(): area_tac = 'each' needs a stock with more than one unit.")

  # tracked value at year yr, one per iteration (NULL if not tracked)
  tracked <- function(metric_name, yr) {
    v <- tracking[metric == metric_name & year == yr, data]
    if (length(v) == 0) NULL else c(v)
  }

  # TAC limits for the iterations flagged by the HCR
  limit <- function(TAC, prev, id) {
    if (length(id) == 0) return(TAC)
    if (!is.na(dtacupp)) TAC[id] <- pmin(TAC[id], prev[id] * dtacupp)
    if (!is.na(dtaclow)) TAC[id] <- pmax(TAC[id], prev[id] * dtaclow)
    TAC
  }

  # area shares: units x iterations
  area_shares <- function() {
    if (identical(catch_area, "recent")) {
      cg <- attr(stk, "catch_gear")
      cu <- if (!is.null(cg) && length(cg) > 0) Reduce(`+`, cg) else catch(stk)
      yrs <- intersect(ac(seq(dy - area_years + 1, dy)), dimnames(cu)$year)
      cu  <- yearSums(cu[, yrs])
      sh  <- t(vapply(units, function(u) c(cu[, , u]), numeric(it)))
      tot <- colSums(sh)
      sh  <- sweep(sh, 2, tot, "/")
      sh[!is.finite(sh)] <- 0
      return(matrix(sh, nrow = length(units), dimnames = list(units, NULL)))
    }
    sh <- if (!is.null(names(catch_area))) catch_area[units] else catch_area
    if (length(sh) != length(units) || anyNA(sh))
      stop("flicc.is(): 'catch_area' needs one share per unit (",
           paste(units, collapse = ", "), ").")
    sh <- sh / sum(sh)
    matrix(rep(sh, it), nrow = length(units), dimnames = list(units, NULL))
  }

  if (area_tac == "split") {

    # stock TAC
    if (ay == iy) {
      if (is.null(initac))
        stop("flicc.is(): initac must be supplied for the first mp() cycle (ay == iy).")
      prev_tac <- rep(sum(initac), it)
    } else {
      prev_tac <- tracked("tac.is", ay - frq)
    }

    TAC <- prev_tac * mult
    id  <- tracking[metric == "decision.hcr" & year == ay, data > 2]
    TAC <- limit(TAC, prev_tac, id)
    if (!is.na(Cmax)) TAC[id] <- pmin(TAC[id], Cmax)

    tac_u <- if (!is.null(catch_area)) area_shares() * rep(TAC, each = length(units)) else NULL

  } else {

    # first-cycle area TACs
    if (ay == iy) {
      if (is.null(initac))
        stop("flicc.is(): initac must be supplied for the first mp() cycle (ay == iy).")
      init_u <- if (!is.null(names(initac)) && all(units %in% names(initac))) {
        initac[units]
      } else if (length(initac) == 1 && !is.null(catch_area)) {
        initac * area_shares()[, 1]
      } else if (length(initac) == length(units)) {
        setNames(initac, units)
      } else {
        stop("flicc.is(): with area_tac = 'each', initac must be named by unit ",
             "or a stock total with catch_area.")
      }
    }

    tac_u <- t(vapply(units, function(u) {
      prev <- if (ay == iy) rep(init_u[[u]], it) else tracked(paste0("tac.is.", u), ay - frq)
      mu   <- tracked(paste0("mult.hcr.", u), ay)
      if (is.null(mu)) mu <- mult                    # no area multiplier: stock-wide
      mu[!is.finite(mu)] <- 1
      tac  <- prev * mu
      id   <- tracking[metric == paste0("decision.hcr.", u) & year == ay, data > 2]
      limit(tac, prev, id)
    }, numeric(it)))
    tac_u <- matrix(tac_u, nrow = length(units), dimnames = list(units, NULL))

    TAC <- colSums(tac_u)
    if (!is.na(Cmax)) {
      sc    <- pmin(1, Cmax / TAC)
      sc[!is.finite(sc)] <- 1
      tac_u <- sweep(tac_u, 2, sc, "*")
      TAC   <- colSums(tac_u)
    }
  }

  # control: the same TAC in every management year of the cycle
  if (is.null(tac_u)) {
    ctrl <- fwdControl(lapply(mys, function(y)
      list(year = y, quant = output, value = TAC)))
  } else {
    ctrl <- fwdControl(unlist(lapply(mys, function(y)
      lapply(units, function(u)
        list(year = y, quant = output, unit = u, value = tac_u[u, ]))),
      recursive = FALSE))
  }

  # tracking
  track(tracking, "tac.is", ac(ay)) <-
    FLQuant(c(TAC), dimnames = list(iter = seq_along(TAC)))
  if (!is.null(tac_u)) {
    for (u in units)
      track(tracking, paste0("tac.is.", u), ac(ay)) <-
        FLQuant(tac_u[u, ], dimnames = list(iter = seq_len(ncol(tac_u))))
  }

  list(ctrl = ctrl, tracking = tracking)
}

## flicc.is() reads mse's data.table tracking with data.table syntax
## (metric, year, data); make the FLicc namespace data.table-aware
.datatable.aware <- TRUE
utils::globalVariables(c("metric", "year", "data"))
