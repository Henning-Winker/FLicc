
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
    b2 = 0.3,
    b3 = 0.5,
    b4 = 0.8,
    blim = 0.10,
    bthr = 0.3,
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
#' b1 = 0.20,
#' b2 = 0.35,
#' b3 = 0.6,
#' b4 = 0.8,
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
    b2 = 0.3,
    b3 = 0.50,
    b4 = 0.75,
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
#'   ngtg, ...); see ?fiticc. If \code{settings$spr_ref} is not given it is
#'   set to \code{spr}, so each fit reports Fspr at the indicator target and
#'   \code{LBIspr()} reuses it.
#' @param spr,thresh Passed to LBIspr() to compute the length-based
#'   indicator (target SPR%, cumulative threshold defining Lref).
#' @param ess.g Numeric scalar or vector of effective sample sizes by gear.
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
#' @return A list containing the (unchanged) stk, the length-based
#'   indicator ind -- an FLQuants with one element per gear (LBIspr()'s own
#'   indicator, unwrapped from FLIndexBiomass) plus a final "SPR" element
#'   (the fitted SPR time series) -- and the tracking object with a
#'   per-iteration "conv.est" convergence flag.
#'
#' @name flicc.sa
#' @rdname flicc.sa
#' @keywords classes
flicc.sa <- function(stk, idx = NULL, args, tracking,
                     lhpar, sel_fun, catch_by_gear = NULL,
                     m_model = "constant",
                     settings = list(
                       CVL = 0.1, GL = 30, catch.sd = 0.05,
                       pop_model = "gtg", obs_model = "dm",
                       ngtg = 13, maxsd = 2, Mpow = 0
                     ),
                     spr = 40, thresh = 0.75, ess.g = 150,
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

  # Fspr at the indicator target is reported by each fit, so LBIspr()
  # reuses it instead of solving for it (see settings$spr_ref in ?fiticc)
  if (is.null(settings$spr_ref)) settings$spr_ref <- spr

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

  # --- LBIspr() on usable iterations --------------------------------------
  ind_i <- vector("list", it)
  ind_i[use] <- lapply(res[use], function(r)
    LBIspr(r$fit, gear = gears, spr = spr, thresh = thresh))

  # NA templates so unusable iterations never inherit another's values
  if (any(use)) {
    template <- ind_i[[which(use)[1]]]
    gear_q <- setNames(lapply(gears, function(g) {
      q <- propagate(index(template[[g]]), it); q[] <- NA; q
    }), gears)
    spr_q <- propagate(res[[which(use)[1]]]$fit$report$spr, it)
    spr_q[] <- NA
  } else {
    na_q <- FLQuant(NA, dimnames = list(year = dimnames(lfd[[1]])$year,
                                        iter = seq_len(it)))
    gear_q <- setNames(rep(list(na_q), length(gears)), gears)
    spr_q  <- na_q
  }

  for (i in which(use)) {
    for (g in gears)
      iter(gear_q[[g]], i) <- index(ind_i[[i]][[g]])
    iter(spr_q, i) <- res[[i]]$fit$report$spr
  }
  ind <- FLQuants(c(gear_q, list(SPR = spr_q)))

  # --- Tracking -----------------------------------------------------------
  track(tracking, "conv.est",     ac(ay)) <- conv           # 1 / 0 / NA (error)
  track(tracking, "sprbound.est", ac(ay)) <- bound          # 1 = SPR at boundary
  track(tracking, "refit.est",    ac(ay)) <- refit          # 1 = rescued by refit
  track(tracking, "used.est",     ac(ay)) <- as.numeric(use) # 0 -> mult = 1

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
#' @param gear Character name or numeric index of the ind element (a gear,
#'   e.g. "Trawl") supplying the trend ratio r.
#' @param output "catch" (TAC, the default) or "effort" -- whichever
#'   fwd.om's projection method expects a relative-to-previous adjustment
#'   applied to. NOT "fbar" -- this rule has no F to set.
#' @param ... Passed through, unused (dispatch compatibility).
lbi.hcr <- function(stk, ind, args, tracking, gear = 1,
                    b1 = 0.10, b2 = 0.3, b3 = 0.50, b4 = 0.75,
                    dlow = -0.20, dopt = 0.00, dhi = 0.15,
                    nyrs = 1, n1 = 2, n2 = 3, output = "catch", ...) {

  FLCore::spread(args)

  # trend ratio -- r_rule()'s own logic, applied directly to the FLQuant
  # (ind[[gear]] already IS the indicator, unlike LBIspr(fit)[[gear]] in
  # hcr_sprlbi(), which needed index() to unwrap it from FLIndexBiomass)
  idx <- ind[[gear]]
  r <- c(yearMeans(tail(idx, n1)) / yearMeans(head(tail(idx, n1 + n2), n2)))
  #browser()
  # SPR-based status rule, unchanged from spr_rule()
  spr_cur <- c(yearMeans(tail(ind[["SPR"]], nyrs)))
  s <- spr_rule(spr_cur, b1 = b1, b2 = b2, b3 = b3, b4 = b4,
                dlow = dlow, dopt = dopt, dhi = dhi)

  # apply s only when r and s agree on direction, exactly as hcr_sprlbi()
  mult <- ifelse((r < 1 & s < 1) | (r > 1 & s > 1), s, 1)

  # decision.hcr trigger flag, read by flicc.is() the same way tacspm.is()

  decision <- ifelse(spr_cur < 0., 3, 1)

  track(tracking, "metric.hcr", ac(ay)) <- ind[["SPR"]][, ac(dy)]
  track(tracking, "decision.hcr", ac(ay)) <-
    FLQuant(decision, dimnames = list(iter = seq_along(decision)))
  track(tracking, "mult.hcr", ac(ay)) <-
    FLQuant(c(mult), dimnames = list(iter = seq_along(decision)))

  #browser()
  ctrl <- fwdControl(list(year = ay, quant = output, value = c(mult)))

  list(ctrl = ctrl, tracking = tracking)
}

#' mp()-compatible implementation system for lbi.hcr()
#'
#' Turns the relative multiplier returned by \code{lbi.hcr()} into a TAC
#' (or effort) advice: the previous TAC (from tracking, or \code{initac} in
#' the first cycle) times the multiplier, with optional limits on the change
#' and an optional split of the advice across units (areas).
#'
#' @param stk The FLStock passed by \code{mp()}.
#' @param ctrl The \code{fwdControl} returned by \code{lbi.hcr()}, holding
#'   the multiplier.
#' @param args,tracking \code{mp()}'s per-cycle args and tracking object.
#' @param output "catch" (default) or "effort".
#' @param dtaclow,dtacupp Optional lower and upper limits on the TAC change
#'   relative to the previous TAC (e.g. 0.85, 1.15), applied to iterations
#'   flagged by the HCR decision.
#' @param Cmax Optional absolute maximum TAC.
#' @param initac TAC used as the previous TAC in the first cycle (required).
#' @param catch_area Optional vector of shares to split the TAC across units.
#'
#' @return A list with \code{ctrl} (a \code{fwdControl} for year
#'   \code{ay + management_lag}) and \code{tracking} (with "tac.is").
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
    for(i in 1:length(units)){
      iters_array[units[i], "value", ] <- TAC * catch_area[i]
    }

    ### intert into ctrl
    ctrl@iters <- iters_array
  }

  track(tracking, "tac.is", ac(ay)) <-
    FLQuant(c(TAC), dimnames = list(iter = seq_along(TAC)))

  list(ctrl = ctrl, tracking = tracking)
}


