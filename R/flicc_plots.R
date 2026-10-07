
#' Plot spawning potential ratio by year
#'
#' Plots spawning potential ratio (SPR) by year from a fitted
#' `flicc_tmb_fit` object using `ggplot2`.
#'
#' SPR values are taken from `fit$report$spr`, which is expected to be an
#' `FLQuant` indexed by year.
#'
#' @param fit A fitted `flicc_tmb_fit` object, list of fits or *FLQuants*.
#' @param title Optional plot title.
#' @param line_colour Colour of the SPR time-series line.
#' @param point_colour Colour of the SPR points.
#' @param linewidth Line width for the SPR time-series line.
#' @param point_size Point size for the SPR observations.
#' @param ymax Optional upper y-axis limit. If `NULL`, this is set to the
#'   maximum of the data, the reference levels, and `0.6`.
#' @param reflines Logical. If `TRUE`, horizontal reference lines are added.
#' @param ref_levels Numeric vector of SPR reference levels.
#' @param ref_colours Character vector of colours for the SPR reference lines.
#' @param ref_linetype Line type for the SPR reference lines.
#'
#' @details
#' Reference lines can be toggled on or off using `reflines`. Their levels and
#' colours can be controlled with `ref_levels` and `ref_colours`.
#'
#' The function expects the fitted object to contain:
#' \describe{
#'   \item{`fit$report$spr`}{An `FLQuant` of spawning potential ratio by year.}
#' }
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_len()]
#'
#' @examples
#' \dontrun{
#'
#' fit <- fiticc(lfd_alfonsino,stklen_alfonsino,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))
#' plot_spr(fit)
#' plot_spr(fit, reflines = FALSE)
#' plot_spr(fit, ref_levels = c(0.7, 0.4, 0.2))
#' }
#'
#' @export
plot_spr <- function(fit,
                     title = NULL,
                     by_year = 2,
                     line_colour = "black",
                     point_colour = "black",
                     linewidth = 0.6,
                     point_size = 2,
                     ymax = NULL,
                     reflines = TRUE,
                     ref_levels = c(0.4, 0.2, 0.1),
                     ref_colours = c("darkgreen", "orange", "red"),
                     ref_linetype = 2) {



  if (inherits(fit, "FLQuants")){
    if(is.null(names(fit))) names <- paste0("run",1:length(fit))
    res <- FLQuants(lapply(fit,function(x){
        FLQuant(x,quant="all")
    }))
  }

  if (inherits(fit, "list")&&inherits(fit[[1]], "flicc_tmb_fit") ){
    if(is.null(names(fit))) names <- paste0("run",1:length(fit))
    res <- FLQuants(lapply(fit,function(x){
      x <- x$report$spr
      FLQuant(x,quant="all")
    }))
  }

  if (inherits(fit, "flicc_tmb_fit")){
      res <- fit$report$spr
    }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  df <- as.data.frame(res)

  yrs <- sort(unique(df$year))
  yr_breaks <- seq(min(yrs), max(yrs), by = by_year)

  if (!all(c("year", "data") %in% names(df))) {
    stop("Could not coerce fit$report$spr to a data.frame with columns 'year' and 'data'.")
  }

  df <- df[!is.na(df$data), ]
  df$year_num <- as.numeric(as.character(df$year))

  if (is.null(ymax)) {
    ymax <- max(df$data, ref_levels, 0.6, na.rm = TRUE)
  }

  if (is.null(title)) {
    title <- "Spawning potential ratio"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = year_num, y = data))
  if(inherits(res, "FLQuant")){
   p <- p+ ggplot2::geom_line(
      colour = line_colour,
      linewidth = linewidth
    ) +
    ggplot2::geom_point(
      colour = point_colour,
      size = point_size
    )
  } else {
    p <- p+ ggplot2::geom_line(
      aes(colour = qname),
      linewidth = linewidth
    ) +
      ggplot2::geom_point(
        aes(colour = qname),
        size = point_size
      )
   }
  p <- p+ ggplot2::scale_x_continuous(breaks = yr_breaks) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::labs(
      x = "Year",
      y = "SPR",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (isTRUE(reflines)) {
    for (i in seq_along(ref_levels)) {
      p <- p + ggplot2::geom_hline(
        yintercept = ref_levels[i],
        colour = ref_colours[i],
        linetype = ref_linetype,linewidth=0.7
      )
    }
  }

  return(p)
}


#' Plot observed and predicted length compositions
#'
#' Plots observed and predicted length compositions from a fitted
#' `flicc_tmb_fit` object using `ggplot2`. Observed length frequencies are
#' taken from `fit$report$obslen` and predicted length compositions from
#' `fit$report$predlen_gear`.
#'
#' The function can subset by year and gear, and can either display separate
#' year-by-gear panels or aggregate across selected gears using the gear
#' weights in `fit$tmb_data$catch_wt`.
#'
#' @param fit A fitted `flicc_tmb_fit` object.
#' @param year Optional year selection. Can be `NULL` for all years, a vector
#'   of year values matching the year dimnames, or numeric indices.
#' @param gear Optional gear selection. Can be `NULL` for all gears, a vector
#'   of gear names, or numeric indices.
#' @param by_gear Logical. If `TRUE`, panels are shown by year and gear. If
#'   `FALSE`, selected gears are combined using `fit$tmb_data$catch_wt`.
#' @param title Optional plot title.
#' @param observed_fill Fill colour for observed length composition bars.
#' @param observed_colour Border colour for observed length composition bars.
#' @param fitted_colour Colour for fitted length composition lines.
#' @param fitted_linewidth Line width for fitted length composition lines.
#' @param len_by Integer giving the spacing of x-axis tick marks for length.
#'
#' @details
#' Observed length frequencies are standardized within each year and gear to
#' proportions before plotting. When `by_gear = FALSE`, the selected gears are
#' combined using the relative gear weights stored in `fit$tmb_data$catch_wt`.
#'
#' The function expects the fitted object to contain:
#' \describe{
#'   \item{`fit$report$obslen`}{An `FLQuants` object of observed length
#'   frequencies by gear.}
#'   \item{`fit$report$predlen_gear`}{An `FLQuants` object of predicted length
#'   compositions by gear.}
#' }
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_spr()]
#'
#' @examples
#' \dontrun{
#' fit <- fiticc(lfd_alfonsino,stklen_alfonsino,sel_fun=c("dsnormal","logistic"),catch_by_gear = c(0.7,0.3))
#' plot_len(fit)
#' plot_len(fit, year = 2022:2024)
#' plot_len(fit, gear = "Gillnet", by_gear = TRUE)
#' plot_len(fit, gear = c("Trawl", "Gillnet"), year = 2022:2024, by_gear = TRUE)
#' }
#'
#' @export
plot_len <- function(fit,
                      year = NULL,
                      gear = NULL,
                      by_gear = FALSE,
                      title = NULL,
                      observed_fill = "grey75",
                      observed_colour = "grey60",
                      fitted_colour = "blue",
                      fitted_linewidth = 1,
                      len_by = 10) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (is.null(fit$report$obslen)) {
    stop("fit$report$obslen not found.")
  }
  if (is.null(fit$report$predlen_gear)) {
    stop("fit$report$predlen_gear not found.")
  }

  obs_flqs <- fit$report$obslen
  fit_flqs <- fit$report$predlen_gear

  gear_labels <- names(obs_flqs)
  if (is.null(gear_labels) || any(gear_labels == "")) {
    stop("fit$report$obslen must be a named FLQuants object.")
  }

  #-----------------------------
  # resolve gear selection
  #-----------------------------
  if (is.null(gear)) {
    gear_sel <- gear_labels
  } else if (is.numeric(gear) || is.integer(gear)) {
    if (!all(gear %in% seq_along(gear_labels))) {
      stop("Numeric 'gear' must be valid gear indices.")
    }
    gear_sel <- gear_labels[gear]
  } else {
    gear_sel <- as.character(gear)
    if (!all(gear_sel %in% gear_labels)) {
      stop("Character 'gear' must match names(fit$report$obslen).")
    }
  }

  #-----------------------------
  # helper: FLQuant -> data.frame
  #-----------------------------
  flq_to_df <- function(flq, gear_name, value_name) {
    df <- as.data.frame(flq)

    len_col <- intersect(c("len", "quant"), names(df))
    if (length(len_col) == 0) {
      stop("Could not find length column ('len' or 'quant') in FLQuant data.")
    }
    len_col <- len_col[1]

    if (!all(c("year", "data") %in% names(df))) {
      stop("FLQuant coercion did not return expected columns 'year' and 'data'.")
    }

    df <- df[, c(len_col, "year", "data")]
    names(df) <- c("len", "year", value_name)

    df$len  <- as.numeric(as.character(df$len))
    df$year <- as.character(df$year)
    df$gear <- gear_name
    df
  }

  #-----------------------------
  # build long data frame from FLQuants
  #-----------------------------
  obs_list <- lapply(gear_sel, function(g) {
    flq_to_df(obs_flqs[[g]], gear_name = g, value_name = "observed")
  })

  fit_list <- lapply(gear_sel, function(g) {
    flq_to_df(fit_flqs[[g]], gear_name = g, value_name = "fitted")
  })

  obs_df <- do.call(rbind, obs_list)
  fit_df <- do.call(rbind, fit_list)

  df <- merge(obs_df, fit_df, by = c("len", "year", "gear"), all = TRUE)

  #-----------------------------
  # resolve year selection
  #-----------------------------
  year_labels <- sort(unique(df$year))

  if (is.null(year)) {
    year_sel <- year_labels
  } else {
    year_chr <- as.character(year)

    if (is.numeric(year) || is.integer(year)) {
      if (all(year_chr %in% year_labels)) {
        year_sel <- year_chr
      } else if (all(year %in% seq_along(year_labels))) {
        year_sel <- year_labels[year]
      } else {
        stop("Requested 'year' not found.")
      }
    } else {
      if (!all(year_chr %in% year_labels)) {
        stop("Requested 'year' not found.")
      }
      year_sel <- year_chr
    }
  }

  df <- df[df$year %in% year_sel, , drop = FALSE]

  # preserve selected order
  df$year <- factor(df$year, levels = year_sel)
  df$gear <- factor(df$gear, levels = gear_sel)

  #-----------------------------
  # standardize observed within year x gear
  #-----------------------------
  split_id <- interaction(df$year, df$gear, drop = TRUE)
  obs_sum <- ave(df$observed, split_id, FUN = function(x) sum(x, na.rm = TRUE))
  df$observed <- ifelse(obs_sum > 0, df$observed / obs_sum, NA_real_)

  #-----------------------------
  # optional aggregate over gears
  #-----------------------------
  if (!by_gear) {
    gear_wt_all <- fit$tmb_data$catch_wt
    names(gear_wt_all) <- fit$tmb_data$gear_names

    gear_wt <- gear_wt_all[gear_sel]
    gear_wt <- gear_wt / sum(gear_wt)

    df$w <- gear_wt[as.character(df$gear)]
    df$observed <- df$observed * df$w
    df$fitted   <- df$fitted * df$w

    df <- stats::aggregate(
      cbind(observed, fitted) ~ len + year,
      data = df,
      FUN = sum,
      na.rm = TRUE
    )
  }

  #-----------------------------
  # x-axis breaks
  #-----------------------------
  lens_all <- sort(unique(df$len))
  x_breaks <- lens_all
  if (!is.null(len_by) && length(lens_all) > 1) {
    x_breaks <- lens_all[seq(1, length(lens_all), by = len_by)]
  }

  #-----------------------------
  # default title
  #-----------------------------
  if (is.null(title)) {
    title <- "Observed and fitted length compositions"
  }

  #-----------------------------
  # plot
  #-----------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = len)) +
    ggplot2::geom_col(
      ggplot2::aes(y = observed),
      fill = observed_fill,
      colour = observed_colour
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = fitted, group = 1),
      colour = fitted_colour,
      linewidth = fitted_linewidth
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::labs(
      x = "Length",
      y = "Proportion",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (by_gear) {
    p <- p + ggplot2::facet_grid(year ~ gear, scales = "free_y")
  } else {
    p <- p + ggplot2::facet_wrap(~ year, scales = "free_y")
  }

  return(p)
}



#' Plot annual LBIspr trajectories by gear
#'
#' Plots annual \code{LBIspr} ratios by gear, with the option to also show
#' annual \code{SPR/SPRtarget} in an additional panel. A horizontal reference
#' line at 1 indicates the target level.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param gear Optional character vector of gear names to include. Defaults to
#'   \code{NULL}, meaning all available gears are plotted.
#' @param spr Numeric SPR target percentage used in \code{LBIspr}, default
#'   \code{40}.
#' @param thresh Numeric cumulative threshold used inside \code{LBIspr},
#'   default \code{0.8}.
#' @param nyears Integer number of years used in the \code{LBIspr}
#'   calculation, default \code{1}.
#' @param plot.spr Logical. If \code{TRUE}, an additional panel is added for
#'   annual \code{SPR / (spr/100)}.
#' @param smoother Logical. If \code{TRUE}, adds a loess smoother to each
#'   panel.
#' @param scale_sel Logical passed to \code{LBIspr()}.
#' @param facets Logical. If \code{TRUE}, facet by series name. If
#'   \code{FALSE}, plot all series in a single panel.
#' @param year_angle Numeric angle for year labels on the x-axis, default
#'   \code{90}.
#' @param line_width Numeric line width for the annual series.
#' @param smooth_width Numeric line width for the smoother.
#' @param colours Optional named character vector of colours. If \code{NULL},
#'   colours are assigned automatically using the same generic palette logic as
#'   \code{plot_lfd()}.
#' @param by_year Integer spacing between x-axis year labels.
#' @param point_pch Plotting symbol for annual observations.
#' @param point_size Point size for annual observations.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' plot_LBIspr(fit, smoother = FALSE)
#' plot_LBIspr(fit, smoother = TRUE)
#' plot_LBIspr(fit, plot.spr = TRUE, facets = TRUE)
#' }
#'
#' @export
plot_LBIspr <- function(fit, gear = NULL, spr = 40, thresh = 0.8,
                        nyears = 1, plot.spr = TRUE,
                        smoother = TRUE, scale_sel = TRUE,
                        facets = TRUE, year_angle = 90,
                        line_width = 0.8, smooth_width = 1.1,
                        colours = NULL, by_year = 2,
                        point_pch = 15, point_size = 1.8) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  idx <- LBIspr(
    fit,
    gear = gear,
    spr = spr,
    thresh = thresh,
    nyears = nyears,
    scale_sel = scale_sel
  )

  dfi <- as.data.frame(idx)
  dfi <- dfi[dfi$slot %in% "index", -c(1:2), drop = FALSE]

  dfi$year <- as.numeric(as.character(dfi$year))
  dfi$cname <- as.character(dfi$cname)

  if (isTRUE(plot.spr)) {
    spry <- as.data.frame(fit$report$spr / (spr / 100))
    spry$cname <- paste0("SPR/SPR", spr)
    spry$year <- as.numeric(as.character(spry$year))
    spry <- spry[, c("year", "data", "cname")]

    dfi <- rbind(dfi[, c("year", "data", "cname")], spry)
  } else {
    dfi <- dfi[, c("year", "data", "cname")]
  }

  dfi$cname <- factor(dfi$cname, levels = unique(dfi$cname))
  qlevels <- levels(dfi$cname)

  if (is.null(colours)) {
    base_cols <- c(
      "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#999999"
    )
    colours <- stats::setNames(
      rep(base_cols, length.out = length(qlevels)),
      qlevels
    )
  } else {
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching cname levels")
    }
    colours <- colours[qlevels]
  }

  yrs <- sort(unique(dfi$year))
  yr_breaks <- seq(min(yrs), max(yrs), by = by_year)

  p <- ggplot2::ggplot(
    dfi,
    ggplot2::aes(x = year, y = data, colour = cname, group = cname)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5) +
    ggplot2::scale_x_continuous(breaks = yr_breaks) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::labs(
      x = "Year",
      y = "LBI ratio",
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = year_angle, vjust = 0.5, hjust = 1
      ),
      axis.title = ggplot2::element_text(face = "bold")
    )

  ymax <- max(dfi$data, 1.4, na.rm = TRUE)

  p <- p + ggplot2::coord_cartesian(ylim = c(0, ymax))

  if (isTRUE(smoother)) {
    p <- p +
      ggplot2::geom_smooth(
        se = FALSE,
        method = "loess",
        formula = y ~ x,
        linewidth = smooth_width,
        alpha = 0.35
      ) +
      ggplot2::geom_point(size = point_size, pch = point_pch)
  } else {
    p <- p +
      ggplot2::geom_line(linewidth = line_width, alpha = 0.95) +
      ggplot2::geom_point(size = point_size, pch = point_pch)
  }

  if (isTRUE(facets)) {
    p <- p + ggplot2::facet_wrap(~ cname, scales = "fixed")
  }

  p
}


#' Plot observed length-frequency data
#'
#' Plots raw or scaled length-frequency distributions (LFDs) from an
#' \code{FLQuants} object, faceted by year and coloured by fishery or gear.
#'
#' The function supports plotting raw frequencies or values scaled relative
#' to the maximum within each \code{FLQuant} (e.g. per gear), allowing
#' comparison of distribution shapes independent of absolute magnitude.
#'
#' @param lfd An \code{FLQuants} object containing observed length-frequency data.
#' @param iter Numeric. Iteration to plot (default = 1).
#' @param len_by Numeric. Step size used to generate x-axis breaks (default = 5).
#' @param yticks Optional numeric vector of y-axis breaks.
#' @param scales Character. Passed to \code{facet_wrap()} (e.g. "fixed", "free_y").
#' @param type Character. Either:
#' \describe{
#'   \item{"raw"}{Plot observed frequencies (default).}
#'   \item{"relmax"}{Scale each \code{FLQuant} by its maximum value.}
#' }
#' @param colours Optional named vector of colours corresponding to \code{qname}.
#' If \code{NULL}, colours are assigned automatically.
#' @param legend_title Character. Title for the legend.
#' @param ncol Optional number of columns in facet layout.
#' @param lwd Numeric. Line width (default = 0.5).
#'
#' @details
#' Length bins are extracted from the \code{FLQuant} dimension names, ensuring
#' consistency with the underlying FLR structure. The function converts the data
#' to a \code{data.frame} for plotting with \pkg{ggplot2}.
#'
#' When \code{type = "relmax"}, each \code{FLQuant} is normalised by its maximum
#' value across the length dimension, preserving relative distribution shapes.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{FLQuants}}, \code{\link{FLQuant}}
#'
#' @author Henning Winker
#'
#' @examples
#' \dontrun{
#' library(FLicc)
#'
#' data(alfonsino)
#'
#' # Raw observed LFDs
#' plot_lfd(lfd_alfonsino)
#'
#' # Relative distributions (scaled to max per gear)
#' plot_lfd(lfd_alfonsino, type = "relmax")
#'
#' # Custom axis ticks
#' plot_lfd(
#'   lfd_alfonsino,
#'   len_by = 20,
#'   yticks = seq(0, 80, by = 20)
#' )
#'
#' # Free y-scale for better visual comparison
#' plot_lfd(lfd_alfonsino, scales = "free_y")
#'
#' # Custom colours
#' plot_lfd(
#'   lfd_alfonsino,
#'   colours = c(Trawl = "steelblue", Gillnet = "tomato")
#' )
#' }
#'
#' @export
plot_lfd <- function(lfd,
                     type=c("raw","relmax"),
                     iter = 1,
                     len_by = 10,
                     yticks = NULL,
                     scales = "fixed",
                     colours = NULL,
                     legend_title = "Fishery",
                     ncol = NULL,
                     lwd = 0.5) {

  type <- match.arg(type)

  L <- as.numeric(dimnames(lfd[[1]])$len)

  lfd <- FLCore::iter(lfd, iter)
  if(type=="relmax")
    lfd <- lapply(lfd, function(x){
      x%/%apply(x,2:5,max)
    })
  lfd <- as.data.frame(lfd)

  x_breaks <- seq(
    from = floor(min(L, na.rm = TRUE) / len_by) * len_by,
    to   = ceiling(max(L, na.rm = TRUE) / len_by) * len_by,
    by   = len_by
  )

  qlevels <- unique(as.character(lfd$qname))

  if (is.null(colours)) {
    base_cols <- c(
      "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#999999"
    )
    colours <- stats::setNames(
      rep(base_cols, length.out = length(qlevels)),
      qlevels
    )
  } else {
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching qname levels")
    }
    colours <- colours[qlevels]
  }

  p <- ggplot2::ggplot(lfd) +
    ggplot2::geom_line(
      ggplot2::aes(x = len, y = data, colour = qname),
      linewidth = lwd
    ) +
    ggplot2::facet_wrap(~year, scales = scales, ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Length", y = "Frequency", colour = legend_title) +
    ggplot2::theme(
      legend.position = "right",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(yticks)) {
    p <- p + ggplot2::scale_y_continuous(breaks = yticks)
  }

  p
}


##' Plot LBSPR-style equilibrium curves
#'
#' Plots equilibrium spawning potential ratio (SPR), relative spawning biomass
#' (\code{SSB/SSB0}), and relative equilibrium yield against apical fishing
#' mortality from an equilibrium \code{FLStockLen} object returned by
#' \code{eqstklen()}.
#'
#' The plot is shown in a single LBSPR-style panel with coloured equilibrium
#' curves, optional horizontal SPR reference lines, an optional vertical line
#' at \code{Fspr}, and optional points showing the current stock state.
#'
#' Current stock state points are taken from \code{eqstk@refpts}, using:
#' \itemize{
#'   \item \code{Fcur} for the x-coordinate
#'   \item \code{SPRcur} for the SPR curve
#'   \item \code{Ycur} for the relative yield curve
#'   \item \code{SSBcur} for the \code{SSB/SSB0} curve
#' }
#'
#' Equilibrium curves are extracted from the \code{eqstk} object as:
#' \itemize{
#'   \item apical fishing mortality from \code{fapl(eqstk)}
#'   \item SPR from \code{eqstk@stock}
#'   \item equilibrium spawning biomass from \code{ssbl(eqstk)}
#'   \item equilibrium yield from \code{eqstk@catch}
#' }
#'
#' Relative spawning biomass is calculated as \code{SSB/SSB0}, where
#' \code{SSB0} is taken as the maximum equilibrium spawning biomass across the
#' fishing mortality sequence. Relative yield is calculated as equilibrium yield
#' divided by the maximum equilibrium yield.
#'
#' Optionally, the x-axis can be trimmed to the last fishing mortality value
#' where relative yield exceeds a small tolerance, plus a user-defined buffer.
#' This helps focus the plot on the informative part of the equilibrium curve.
#'
#' @param eqstk An equilibrium \code{FLStockLen} object returned by
#'   \code{eqstklen()}.
#' @param reflines Numeric vector of SPR-style reference levels to draw as
#'   horizontal dashed lines. Defaults to \code{c(0.4, 0.2, 0.1)}.
#' @param refcols Character vector of colours for \code{reflines}. Defaults to
#'   \code{c("darkgreen", "orange", "red")}.
#' @param show_current Logical. If \code{TRUE}, add current-state reference
#'   points from \code{eqstk@refpts}.
#' @param show_legend Logical. If \code{TRUE}, show the legend for the
#'   equilibrium curves.
#' @param show_s Logical. If \code{TRUE}, annotate the plot with the steepness
#'   value.
#' @param s Optional Beverton-Holt steepness value to label. If \code{NULL},
#'   the function first tries \code{attr(eqstk, "lhpar")["s"]}, then
#'   \code{attr(eqstk, "steepness")}.
#' @param xlab Character x-axis label. Defaults to \code{"Apical F"}.
#' @param ylab Optional y-axis label. If \code{NULL}, defaults to
#'   \code{"Proportion"}.
#' @param title Optional plot title. If \code{NULL}, defaults to
#'   \code{"Equilibrium curves"}.
#' @param line_width Numeric line width for the equilibrium curves.
#' @param point_size Numeric point size for current-state reference points.
#' @param colours Optional named character vector giving colours for the three
#'   curves. Names should match \code{c("SPR", "SSB/SSB0", "Relative Yield")}.
#' @param trim_x Logical. If \code{TRUE}, trim the x-axis to the last fishing
#'   mortality value where relative yield exceeds \code{yield_tol}, plus a
#'   buffer.
#' @param trim_buffer Numeric proportional buffer added to the trimmed upper
#'   x-limit. The default \code{0.2} corresponds to 20%.
#' @param yield_tol Numeric tolerance used to determine where relative yield is
#'   effectively zero when trimming the x-axis.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{eqstklen}}, \code{\link{prbrp_flicc}}
#'
#' @examples
#' \dontrun{
#' eqstk <- eqstklen(fit)
#'
#' plot_eqcurves(eqstk)
#'
#' plot_eqcurves(
#'   eqstk,
#'   reflines = c(0.5, 0.4, 0.2),
#'   refcols = c("grey40", "darkgreen", "orange")
#' )
#'
#' plot_eqcurves(
#'   eqstk,
#'   trim_x = TRUE,
#'   trim_buffer = 0.15,
#'   title = "Equilibrium SPR, biomass and yield"
#' )
#' }
#'
#' @export

plot_eqcurves <- function(eqstk,
                          reflines = c(0.4, 0.2, 0.1),
                          refcols = c("darkgreen", "orange", "red"),
                          show_current = TRUE,
                          show_legend = TRUE,
                          show_s = TRUE,
                          s = NULL,
                          xlab = "Apical F",
                          ylab = NULL,
                          title = NULL,
                          line_width = 1.1,
                          point_size = 3,
                          colours = NULL,
                          trim_x = TRUE,
                          trim_buffer = 0.2,
                          yield_tol = 1e-6) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (!methods::is(eqstk, "FLStockLen")) {
    stop("'eqstk' must be an FLStockLen object returned by eqstklen().")
  }

  if (is.null(eqstk@refpts) || length(eqstk@refpts) == 0) {
    stop("eqstk@refpts is empty. Current reference points are required.")
  }

  rp <- eqstk@refpts
  rp_names <- dimnames(rp)[[1]]

  get_rp <- function(x) {
    if (!x %in% rp_names) return(NA_real_)
    as.numeric(rp[x])
  }

  Fcur   <- get_rp("Fcur")
  SPRcur <- get_rp("SPRcur")
  Ycur   <- get_rp("Ycur")
  SSBcur <- get_rp("SSBcur")
  Fspr   <- get_rp("Fspr")

  # equilibrium curves
  Fap <- as.numeric(fapl(eqstk))
  SPR <- as.numeric(eqstk@stock)
  SSB <- as.numeric(ssbl(eqstk))
  Y   <- as.numeric(eqstk@catch)

  if (length(Fap) == 0 || length(SPR) == 0 || length(SSB) == 0 || length(Y) == 0) {
    stop("eqstk does not contain the expected equilibrium quantities.")
  }

  SSB0 <- max(SSB, na.rm = TRUE)
  Ymax <- max(Y, na.rm = TRUE)

  plotdf <- rbind(
    data.frame(
      Fap = Fap,
      value = SPR,
      metric = "SPR"
    ),
    data.frame(
      Fap = Fap,
      value = SSB / SSB0,
      metric = "SSB/SSB0"
    ),
    data.frame(
      Fap = Fap,
      value = Y / Ymax,
      metric = "Relative Yield"
    )
  )

  plotdf$metric <- factor(
    plotdf$metric,
    levels = c("SPR", "SSB/SSB0", "Relative Yield")
  )

  if (is.null(colours)) {
    colours <- c(
      "SPR" = "#F8766D",
      "SSB/SSB0" = "#00BA38",
      "Relative Yield" = "#619CFF"
    )
  }

  if (is.null(title)) {
    title <- "Equilibrium curves"
  }

  p <- ggplot2::ggplot(
    plotdf,
    ggplot2::aes(x = Fap, y = value, colour = metric)
  ) +
    ggplot2::geom_line(linewidth = line_width, show.legend = show_legend) +
    ggplot2::scale_colour_manual(
      values = colours,
      labels = c("SPR", "SSB/SSB0", "Relative yield")
    )+
    ggplot2::labs(
      x = xlab,
      y = if (is.null(ylab)) "Proportion" else ylab,
      title = title,
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "top" else "none",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  # SPR reference lines only in SPR panel
  if (length(reflines) > 0) {
    refdf <- data.frame(
      metric = "SPR",
      yint = reflines,
      col = rep_len(refcols, length(reflines))
    )

    p <- p +
      ggplot2::geom_hline(
        data = refdf,
        ggplot2::aes(yintercept = yint),
        colour = refdf$col,
        linetype = 2,
        linewidth = 0.7,
        inherit.aes = FALSE
      )
  }

  # optional Fspr vertical line in all panels
  if (is.finite(Fspr)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = Fspr,
        linetype = 2,
        colour = "grey40",
        linewidth = 0.6
      )
  }

  if (isTRUE(show_current)) {
    curdf <- rbind(
      data.frame(Fap = Fcur, value = SPRcur, metric = "SPR"),
      data.frame(Fap = Fcur, value = SSBcur, metric = "SSB/SSB0"),
      data.frame(Fap = Fcur, value = Ycur, metric = "Relative Yield")
    )

    curdf <- curdf[is.finite(curdf$Fap) & is.finite(curdf$value), , drop = FALSE]

    if (nrow(curdf) > 0) {
      curdf$metric <- factor(curdf$metric, levels = levels(plotdf$metric))

      p <- p +
        ggplot2::geom_point(
          data = curdf,
          ggplot2::aes(x = Fap, y = value),
          inherit.aes = FALSE,
          size = point_size,
          colour = "black"
        )
    }
  }

  if (is.null(s) && !is.null(attr(eqstk, "lhpar"))) {
    lh <- attr(eqstk, "lhpar")
    if ("s" %in% rownames(lh)) {
      s <- as.numeric(lh["s"])
    }
  }

  if (is.null(s)) {
    s <- attr(eqstk, "steepness")
  }

  if (isTRUE(show_s) && !is.null(s) && is.finite(s)) {
    sprdf <- plotdf[plotdf$metric == "SPR", , drop = FALSE]
    xmax <- max(sprdf$Fap, na.rm = TRUE)
    ymax <- max(sprdf$value, na.rm = TRUE)



    if (isTRUE(trim_x)) {
      Y <- as.numeric(eqstk@catch)
      Fap <- as.numeric(fapl(eqstk))
      Yrel <- Y / max(Y, na.rm = TRUE)

      keep <- which(Yrel > yield_tol)

      if (length(keep) > 0) {
        xmax <- max(Fap[keep], na.rm = TRUE) * (1 + trim_buffer)

        if (isTRUE(show_current) && is.finite(Fcur)) {
          xmax <- max(xmax, Fcur * 1.05, na.rm = TRUE)
        }

        xmax <- min(xmax, max(Fap, na.rm = TRUE))

        p <- p + ggplot2::coord_cartesian(xlim = c(0, xmax))
      }
    }
    p <- p +
      ggplot2::annotate(
        "text",
        x = xmax,
        y = ymax,
        label = paste0("s = ", format(round(s, 2), nsmall = 2)),
        hjust = 1.02,
        vjust = -0.3,
        size = 4
      )
  }
  p <- p +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 11),
      legend.key.width = grid::unit(1.5, "cm")
    )

  p
}


#' Plot FAO-style stock status from SPR relative to target
#'
#' Plots annual spawning potential ratio relative to a target SPR level,
#' \code{SPR / SPRtarget}, from a fitted \code{flicc_tmb_fit} object using
#' \code{ggplot2}. The background is shaded using FAO-style stock status bands.
#'
#' By default, the status zones are:
#' \itemize{
#'   \item \code{0.0 - 0.8}: overfished
#'   \item \code{0.8 - 1.2}: sustainably fished
#'   \item \code{> 1.2}: underfished
#' }
#'
#' The two upper zones are shown with the same blue fill and separated by a
#' white dashed line, following the intended FAO/SOSI visual style.
#'
#' @param fit A fitted \code{flicc_tmb_fit} object.
#' @param spr.tgt Numeric SPR target in percent. Default is \code{40}, so the
#'   plotted quantity is \code{SPR / 0.4}.
#' @param title Optional plot title. Defaults to \code{"FAO Status"}.
#' @param line_colour Colour of the time-series line.
#' @param point_colour Colour of the annual points.
#' @param linewidth Line width for the time-series line.
#' @param point_size Point size for the annual observations.
#' @param side_alpha Transparency of the side status bar.
#' @param ymax Optional upper y-axis limit. If \code{NULL}, this is set to the
#'   maximum of the scaled SPR series and \code{1.5}.
#' @param overfished_col Fill colour for the overfished zone.
#' @param sustainable_col Fill colour for the sustainable and underfished
#'   zones.
#' @param band_alpha Transparency of the background status bands.
#' @param separator_colour Colour of the dashed separator lines between FAO
#'   status zones.
#' @param separator_linetype Line type for the dashed separator lines.
#' @param separator_linewidth Line width for the dashed separator lines.
#' @param y_breaks Numeric vector of y-axis breaks.
#'
#' @return A \code{ggplot} object.
#' @export
plot_lbfao <- function(fit,
                       spr.tgt = 40,
                       title = "FAO Status",
                       line_colour = "black",
                       point_colour = "black",
                       linewidth = 0.6,
                       point_size = 3,
                       side_alpha = 0.18,
                       ymax = NULL,
                       overfished_col = "orange",
                       sustainable_col = "dodgerblue3",
                       band_alpha = 0.8,
                       separator_colour = "white",
                       separator_linetype = 3,
                       separator_linewidth = 0.7,
                       y_breaks = seq(0, 10, 0.2)) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (is.null(fit$report$spr)) {
    stop("fit$report$spr not found.")
  }

  df <- as.data.frame(fit$report$spr)

  if (!all(c("year", "data") %in% names(df))) {
    stop("Could not coerce fit$report$spr to a data.frame with columns 'year' and 'data'.")
  }

  df <- df[!is.na(df$data), , drop = FALSE]
  df$year_num <- as.numeric(as.character(df$year))
  df$status_ratio <- df$data / (spr.tgt / 100)

  if (nrow(df) == 0) {
    stop("No non-missing SPR values found in fit$report$spr.")
  }

  if (is.null(ymax)) {
    ymax <- max(df$status_ratio, 1.5, na.rm = TRUE)
  }

  xmin <- min(df$year_num, na.rm = TRUE) - 0.5
  xmax <- max(df$year_num, na.rm = TRUE) + 0.5

  band_df <- data.frame(
    xmin = c(xmin, xmin, xmin),
    xmax = c(xmax, xmax, xmax),
    ymin = c(0, 0.8, 1.2),
    ymax = c(0.8, 1.2, ymax),
    fill = factor(
      c("overfished", "sustainable", "underfished"),
      levels = c("overfished", "sustainable", "underfished", "unsustainable")
    )
  )

  status_cut <- 0.8
  bar_xmin <- xmax + 0.15
  bar_xmax <- xmax + 0.50

  status_bar <- data.frame(
    xmin = c(bar_xmin, bar_xmin),
    xmax = c(bar_xmax, bar_xmax),
    ymin = c(0, status_cut),
    ymax = c(status_cut, ymax),
    fill = factor(
      c("unsustainable", "sustainable"),
      levels = c("overfished", "sustainable", "underfished", "unsustainable")
    )
  )

  fill_vals <- c(
    overfished = overfished_col,
    sustainable = sustainable_col,
    underfished = sustainable_col,
    unsustainable = overfished_col
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = year_num, y = status_ratio)) +
    ggplot2::geom_rect(
      data = band_df,
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = fill
      ),
      inherit.aes = FALSE,
      alpha = band_alpha,
      colour = NA
    ) +
    ggplot2::geom_rect(
      data = status_bar,
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = fill
      ),
      inherit.aes = FALSE,
      alpha = side_alpha,
      colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = fill_vals,
      guide = "none"
    ) +
    ggplot2::geom_hline(
      yintercept = 1.2,
      colour = separator_colour,
      linetype = separator_linetype,
      linewidth = separator_linewidth
    ) +
    ggplot2::geom_line(
      colour = line_colour,
      linewidth = linewidth
    ) +
    ggplot2::geom_point(
      shape = 21,
      fill = "white",
      colour = point_colour,
      size = point_size,
      stroke = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x = (bar_xmin + bar_xmax) / 2,
      y = status_cut / 2,
      label = "UNSUSTAINABLE",
      angle = 90,
      size = 3,
      fontface = "bold",
      colour = "grey30"
    ) +
    ggplot2::annotate(
      "text",
      x = (bar_xmin + bar_xmax) / 2,
      y = status_cut + (ymax - status_cut) / 2,
      label = "SUSTAINABLE",
      angle = 90,
      size = 3,
      fontface = "bold",
      colour = "grey30"
    ) +
    ggplot2::scale_x_continuous(
      breaks = df$year_num,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks
    ) +
    ggplot2::coord_cartesian(
      xlim = c(xmin, bar_xmax + 0.25),
      ylim = c(0, ymax),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::labs(
      x = "Year",
      y = paste0("SPR / SPR", spr.tgt),
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(5.5, 35, 5.5, 5.5)
    )

  p
}

#' Plot advice diagnostics for length-based per-recruit analysis
#'
#' Plots annual trajectories of spawning potential ratio (SPR) and fishing
#' mortality from simplified \code{FLStock} or \code{FLStockR} objects,
#' including uncertainty intervals and optional reference points.
#'
#' The function is designed for outputs derived from length-based per-recruit
#' workflows such as \code{flicc2FLStockR()}, where:
#'
#' \itemize{
#'   \item \code{ssb()} is interpreted as SPR or relative SPR,
#'   \item \code{fbar()} is interpreted as applied fishing mortality or
#'         relative fishing mortality.
#' }
#'
#' If several stocks are supplied in an \code{FLStocks} object, the function
#' overlays them and plots median and interval summaries for each metric.
#'
#' @param object An object of class \code{FLStock}, \code{FLStockR}, or
#'   \code{FLStocks}.
#' @param panels Integer vector indicating which panels to plot:
#'   \code{1} for SPR and \code{2} for fishing mortality. The default is
#'   \code{c(1, 2)}.
#' @param plotrefs Logical. If \code{TRUE}, reference points stored in
#'   \code{object@refpts} are added to the plot.
#' @param probs Numeric vector of probabilities used for plotting quantile
#'   intervals. The default is \code{c(0.05, 0.2, 0.50, 0.8, 0.95)}.
#' @param by_year Integer spacing between x-axis year labels.
#' @param colour Fill colour used for uncertainty intervals when plotting a
#'   single stock. The default is \code{"dodgerblue"}.
#' @param ncol Number of columns in the facet layout. Currently set internally
#'   to 1 in the function.
#' @param label.size Numeric size of the reference point labels.
#' @param ssbQ Integer spawning season or quarter used when extracting
#'   \code{ssb()} for seasonal models. Default is \code{1}.
#' @param recQ Included for consistency with related plotting functions. Not
#'   currently used in the function.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The plotted metrics are:
#'
#' \itemize{
#'   \item \strong{SPR}: computed from \code{ssb(x)} using
#'         \code{apply(ssb(x)[,,,ssbQ], c(2,6), sum)}
#'   \item \strong{F}: computed from \code{fbar(x)} using
#'         \code{apply(fbar(x), c(2,6), mean)}
#' }
#'
#' For a single stock, shaded quantile bands are drawn using the selected
#' \code{colour}. For multiple stocks, interval bands and median lines are
#' overlaid by stock.
#'
#' Reference points are grouped by names starting with \code{"B"} or
#' \code{"F"} and coloured heuristically according to common advice categories,
#' for example limit, precautionary, target, or boundary reference points.
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#'
#' stk <- flicc2FLStockR(fit)
#' stkr <- flicc2FLStockR(fit,rel=TRUE)
#'
#' plot_LBAdvice(stk)
#' plot_LBAdvice(stkr)
#' plot_LBAdvice(stkr, panels = 1)
#' plot_LBAdvice(stkr, panels = 2, plotrefs = FALSE)
#'
#' @export
plot_LBAdvice <- function(object,panels=c(1,2),plotrefs=TRUE,probs=c(0.05,0.2,0.50,0.8,0.95),by_year = 2,colour="dodgerblue",ncol=NULL,label.size=2.5,ssbQ = 1,recQ=1){

  old_lifecycle <- getOption("lifecycle_verbosity")
  options(lifecycle_verbosity = "quiet")
  on.exit(options(lifecycle_verbosity = old_lifecycle), add = TRUE)

  if(class(object)%in%c("FLStock","FLStockR"))
    object = FLStocks(object)


  if(is.null(ncol)&&length(panels)>1){
    ncol=2
  } else {
    ncol=1
  }

  stks <- FLStocks(lapply(object,function(stock){
    landings = stock@landings
    stock@landings =computeLandings(stock)
    stock@landings[is.na(stock@landings)] = landings[is.na(stock@landings)]
    stock@discards = computeDiscards(stock)
    stock@catch = computeCatch(stock)

    stock
  }))



  rp <- stks[[1]]@refpts
  year <- an(dimnames(stks[[1]])$year)

  styr = endyr = NULL
  for(i in 1:length(stks)){
    styr = min(styr, range(stks[[i]])["minyear"])
    endyr = max(endyr, range(stks[[i]])["maxyear"])
  }
  iv = ceiling(length(styr:endyr)/8)



  if(length(stks)==1) stks = stks[[1]]


  leg = theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
                legend.key.height = unit(0.5, 'cm'), #change legend key height
                legend.key.width = unit(0.3, 'cm'), #change legend key width
                legend.title = element_blank(), #change legend title font size
                legend.text = element_text(size=7),
                axis.text.x= element_text(size=7)) #change legend text font size


      ### PLOT
      p = suppressMessages(ggplotFL::plot(stks,
                         metrics=list(SPR=function(x)apply(ssb(x)[,,,ssbQ],c(2,6),sum),F=function(x)apply(fbar(x),c(2,6),mean))[panels])+
        ylim(c(0, NA))+ theme_bw()+leg+
        xlab("Year")+ facet_wrap(~qname, scales="free",ncol=ncol))

      if(length(stks)==1){
        p = p +ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(1,3,5)], alpha=0.2) +
          ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(2,3,4)], alpha=0.4)
      }  else {
        p = p +ggplotFL::geom_flquantiles(colour=NA, probs=probs[c(1,3,5)], alpha=0.2) +
          ggplotFL::geom_flquantiles(colour=NA, probs=probs[c(2,3,4)], alpha=0.4)+
          ggplotFL::geom_flquantiles( probs=probs[c(3)])

      }

  if(class(stks)%in%c("FLStockR","FLStock")){
    stk=stks
  } else {
    stk = stks[[1]]
  }

  xy =quantile(dims(stk)$minyear:dims(stk)$maxyear,c(0.2,0.45,0.75,0.6,0.3,0.5,0.1))


    Fs = FLPar(an(rp[substr(rownames(rp),1,1)%in%"F"]),
               params=rownames(rp)[substr(rownames(rp),1,1)%in%"F"])
    nf = length(Fs)
    Bs = FLPar(an(rp[substr(rownames(rp),1,1)%in%"B"]),
               params=paste0(rownames(rp)[substr(rownames(rp),1,1)%in%"B"],""))

    nb = length(Bs)

  if(plotrefs){

    Bsc = data.frame(red =  rownames(Bs)%in%c("Bcrit","Blim"),
                     orange = rownames(Bs)%in%c("Bpa","Bthr"),
                     green = !rownames(Bs)%in%c("B0","SPR0","Bpa","Bthr","Btri","Btrigger","Bcrit","Blim"),
                     blue = rownames(Bs)%in%c("B0","SPR0","Btrigger","Btri")
    )
    Fsc = data.frame(red =  rownames(Fs)%in%c("Flim","Fext","Fcrash","Fcrit"),
                     orange = rownames(Fs)%in%c("Fpa","Fthr","Fthresh","Ftri","Fp0.5"),
                     green = !rownames(Fs)%in%c("Fpa","Fthr","Fcrit","Flim","Fext","Fcrash","Fp0.5","Flower","Fupper","Fmys.low","Fmsy.up"),
                     blue = rownames(Fs)%in%c("Flower","Fupper","Fmsy.low","Fmsy.up")
    )


    bcol = NULL
    for(i in 1:length(Bs))  bcol = c(bcol,c("red","orange","darkgreen","blue")[which(Bsc[i,]==TRUE)])
    fcol = NULL
    for(i in 1:length(Fs))  fcol = c(fcol,c("red","darkorange","darkgreen","blue")[which(Fsc[i,]==TRUE)])

   colo <- unname(unlist(list(b=bcol,f=fcol)[panels]))
   posx <- unname(unlist(list(b=xy[1:length(Bs)],f=xy[1:length(Fs)])[panels]))


    qn = c("SPR","F")[panels]

      #posx = posx[-length(posx)] #if(any(!is.na(Ys))){
      #colo = colo[-length(colo)]
      #posx = posx[-c(8)]
      #colo = colo[-c(8)]
      #} else {
      #  posx = posx[-c(7:8)]
      #  colo = colo[-c(7:8)]

      #}
      p = p+facet_wrap(~qname,scales="free_y",ncol= ncol)

    fps =FLPars(SSB  =Bs,
                F    =Fs,
                )[panels]
    fps@names = qn
    ggp = ggplotFL::geom_flpar(data=fps,x=posx,colour=colo)

    ggp[[2]]$aes_params$size=label.size
    p = p +ggp
  }

  yr_breaks <- rev(seq(max(year),min(year),-by_year))
  p <-  p+ggplot2::scale_x_continuous(breaks = yr_breaks)
  return(p)
}
# }}}


#' Plot selectivity curves by gear
#'
#' Plots selectivity-at-length either from a fitted \code{"flicc_tmb_fit"}
#' object or from an \code{FLStockLen} object.
#'
#' For fitted objects, selectivity is taken from \code{fit$report$sel_gear}.
#' For \code{FLStockLen} objects, selectivity is approximated by scaling
#' \code{harvest()} within each year to a maximum of 1.
#'
#' @param object A \code{"flicc_tmb_fit"} or \code{FLStockLen} object.
#' @param year Optional year or vector of years to include. Defaults to
#'   \code{NULL}, meaning all available years are shown.
#' @param gear Optional character vector of gears to include. Defaults to
#'   \code{NULL}, meaning all available gears are shown.
#' @param facets Logical. If \code{TRUE}, facet by year for \code{FLStockLen}
#'   objects. Ignored for fitted objects unless annual selectivity is available.
#' @param colours Optional named character vector of colours. If \code{NULL},
#'   colours are assigned automatically using the same generic palette logic as
#'   \code{plot_lfd()}.
#' @param line_width Numeric line width for selectivity curves.
#' @param len_by Numeric. Step size used to generate x-axis breaks (default = 5).
#'
#' @return A \code{ggplot2} object.
#'
#' @export
plot_sel <- function(object, year = NULL, gear = NULL,
                     facets = FALSE, colours = NULL,
                     line_width = 0.7, len_by = 5) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  base_cols <- c(
    "#D55E00", "#0072B2", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#999999"
  )

  base_cols <- rep(base_cols,10)

  #--------------------------------------------------
  # flicc_tmb_fit
  #--------------------------------------------------
  if (inherits(object, "flicc_tmb_fit")) {

    if (is.null(object$report$sel_gear)) {
      stop("object$report$sel_gear not found.")
    }

    sel <- object$report$sel_gear
    df <- as.data.frame(sel)
    lens <- sort(unique(df$len))
    len_breaks <- seq(min(lens, na.rm = TRUE), max(lens, na.rm = TRUE), by = len_by)

    # expected FLQuants names are gear names
    df$gear <- as.character(df$qname)

    if (!is.null(gear)) {
      df <- df[df$gear %in% gear, , drop = FALSE]
    }

    glevels <- unique(df$gear)

    if (is.null(colours)) {
      colours <- stats::setNames(
        rep(base_cols, length.out = length(glevels)),
        glevels
      )
    } else {
      if (is.null(names(colours))) {
        stop("'colours' must be a named vector with names matching gear names")
      }
      colours <- colours[glevels]
    }

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = len, y = data, colour = gear, group = gear)
    ) +
      ggplot2::geom_line(linewidth = line_width) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(breaks = len_breaks) +
      ggplot2::labs(x = "Length", y = "Selectivity", colour = NULL) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(face = "bold")
      )

    return(p)
  }

  #--------------------------------------------------
  # FLStockLen
  #--------------------------------------------------
  if (inherits(object, "FLStockLen")) {


      h <- harvest(object)
      out <- h %/% apply(h, 2:6, max, na.rm = TRUE)


    if (!is.null(year)) {
      out <- out[, ac(year)]
    }

    if (!is.null(gear)) {
      keep <- names(out) %in% gear
      out <- out[keep]
    }

    df <- as.data.frame(out)
    df$year <- ac(df$year)
    lens <- sort(unique(df$len))
    len_breaks <- seq(min(lens, na.rm = TRUE), max(lens, na.rm = TRUE), by = len_by)



    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = len, y = data, colour = year)
    ) +
      ggplot2::geom_line(linewidth = line_width) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(breaks = len_breaks) +
      ggplot2::labs(x = "Length", y = "Selectivity", colour = NULL) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(face = "bold")
      )

    if (isTRUE(facets)) {
      p <- p + ggplot2::facet_wrap(~ year)
    }

    return(p)
  }

  stop("object must be of class 'flicc_tmb_fit' or 'FLStockLen'.")
}

#' Plot natural mortality-at-length curves
#'
#' Plots mean natural mortality-at-length, \eqn{M(L)}, from an
#' \code{FLStockLen} or \code{FLStocks} object. If multiple years are present,
#' the plot uses the mean over the last \code{nyears}.
#'
#' If a single \code{FLStockLen} object is supplied, it is internally coerced
#' to an \code{FLStocks} object named \code{M_l}.
#'
#' @param object An \code{FLStockLen} or \code{FLStocks} object.
#' @param nyears Integer number of most recent years over which to average
#'   \code{m()}. Default is \code{1}.
#' @param len_by Numeric spacing for length-axis tick marks. Default is
#'   \code{5}.
#' @param colours Optional named character vector of colours. If \code{NULL},
#'   colours are assigned automatically using the same generic palette logic as
#'   other FLicc plotting functions.
#' @param line_width Numeric line width for mortality curves. Default is
#'   \code{0.8}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' plot_m(stklen)
#' plot_m(stklen, nyears = 3)
#' plot_m(stklen, len_by = 5)
#' plot_m(stklen, colours = c(M_l = "black"))
#' }
#'
#' @export
plot_m <- function(object, nyears = 1, len_by = 5,
                   colours = NULL, line_width = 0.8) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (inherits(object, "FLStockLen")) {
    object <- FLCore::FLStocks(M_l = object)
  }

  year <- dimnames(object[[1]])$year
  yrs <- tail(year, min(nyears, length(year)))

  out <- FLCore::FLQuants(lapply(object, function(x) {
    FLCore::yearMeans(FLCore::m(x)[, yrs])
  }))

  df <- as.data.frame(out)

  qlevels <- unique(as.character(df$qname))

  if (is.null(colours)) {
    base_cols <- c(
      "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#999999"
    )
    colours <- stats::setNames(
      rep(base_cols, length.out = length(qlevels)),
      qlevels
    )
  } else {
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching series names")
    }
    colours <- colours[qlevels]
  }

  lens <- sort(unique(df$len))
  len_breaks <- seq(min(lens, na.rm = TRUE), max(lens, na.rm = TRUE), by = len_by)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = len, y = data, colour = qname)) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::scale_x_continuous(breaks = len_breaks) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::ylab("M(L)") +
    ggplot2::xlab("Length") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (length(out) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p + ggplot2::scale_y_continuous(n.breaks = 8)

  return(p)
}


#' Plot maturity-at-length curves
#'
#' Plots mean maturity-at-length from an \code{FLStockLen} or
#' \code{FLStocks} object. If multiple years are present, the plot uses the
#' mean over the last \code{nyears}.
#'
#' If a single \code{FLStockLen} object is supplied, it is internally coerced
#' to an \code{FLStocks} object named \code{mat_l}.
#'
#' @param object An \code{FLStockLen} or \code{FLStocks} object.
#' @param nyears Integer number of most recent years over which to average
#'   \code{mat()}. Default is \code{1}.
#' @param len_by Numeric spacing for length-axis tick marks. Default is
#'   \code{5}.
#' @param colours Optional named character vector of colours. If \code{NULL},
#'   colours are assigned automatically using the same generic palette logic as
#'   other FLicc plotting functions.
#' @param line_width Numeric line width for maturity curves. Default is
#'   \code{0.8}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' plot_mat(stklen)
#' plot_mat(stklen, nyears = 3)
#' plot_mat(stklen, len_by = 5)
#' plot_mat(stklen, colours = c(mat_l = "black"))
#' }
#'
#' @export
plot_mat <- function(object, nyears = 1, len_by = 5,
                     colours = NULL, line_width = 0.8) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (inherits(object, "FLStockLen")) {
    object <- FLCore::FLStocks(mat_l = object)
  }

  year <- dimnames(object[[1]])$year
  yrs <- tail(year, min(nyears, length(year)))

  out <- FLCore::FLQuants(lapply(object, function(x) {
    FLCore::yearMeans(mat(x)[, yrs])
  }))

  df <- as.data.frame(out)

  qlevels <- unique(as.character(df$qname))

  if (is.null(colours)) {
    base_cols <- c(
      "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#999999"
    )
    colours <- stats::setNames(
      rep(base_cols, length.out = length(qlevels)),
      qlevels
    )
  } else {
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching series names")
    }
    colours <- colours[qlevels]
  }

  lens <- sort(unique(df$len))
  len_breaks <- seq(min(lens, na.rm = TRUE), max(lens, na.rm = TRUE), by = len_by)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = len, y = data, colour = qname)) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::scale_x_continuous(breaks = len_breaks) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Maturity") +
    ggplot2::xlab("Length (cm)") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )+geom_hline(yintercept = 0.5,linetype=2)

  if (length(out) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}


#' Plot catch weight-at-length curves
#'
#' Plots mean catch weight-at-length from an \code{FLStockLen} or
#' \code{FLStocks} object. If multiple years are present, the plot uses the
#' mean over the last \code{nyears}.
#'
#' If a single \code{FLStockLen} object is supplied, it is internally coerced
#' to an \code{FLStocks} object named \code{catch.wt_l}.
#'
#' @param object An \code{FLStockLen} or \code{FLStocks} object.
#' @param nyears Integer number of most recent years over which to average
#'   \code{catch.wt()}. Default is \code{1}.
#' @param len_by Numeric spacing for length-axis tick marks. Default is
#'   \code{5}.
#' @param colours Optional named character vector of colours. If \code{NULL},
#'   colours are assigned automatically using the same generic palette logic as
#'   other FLicc plotting functions.
#' @param line_width Numeric line width for catch weight curves. Default is
#'   \code{0.8}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' plot_lw(stklen)
#' plot_lw(stklen, nyears = 3)
#' plot_lw(stklen, len_by = 5)
#' plot_lw(stklen, colours = c(catch.wt_l = "black"))
#' }
#'
#' @export
plot_lw <- function(object, nyears = 1, len_by = 5,
                    colours = NULL, line_width = 0.8) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (inherits(object, "FLStockLen")) {
    object <- FLCore::FLStocks(catch.wt_l = object)
  }

  year <- dimnames(object[[1]])$year
  yrs <- tail(year, min(nyears, length(year)))

  out <- FLCore::FLQuants(lapply(object, function(x) {
    FLCore::yearMeans(FLCore::catch.wt(x)[, yrs])
  }))

  df <- as.data.frame(out)

  qlevels <- unique(as.character(df$qname))

  if (is.null(colours)) {
    base_cols <- c(
      "#D55E00", "#0072B2", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#999999"
    )
    colours <- stats::setNames(
      rep(base_cols, length.out = length(qlevels)),
      qlevels
    )
  } else {
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching series names")
    }
    colours <- colours[qlevels]
  }

  lens <- sort(unique(df$len))
  len_breaks <- seq(min(lens, na.rm = TRUE), max(lens, na.rm = TRUE), by = len_by)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = len, y = data, colour = qname)) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::scale_x_continuous(breaks = len_breaks) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Weight  (kg)") +
    ggplot2::xlab("Length (vm)") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (length(out) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p <- p + ggplot2::scale_y_continuous(n.breaks = 6)

  return(p)
}


#' Plot observed vs expected length composition at Fspr
#'
#' A plotting helper for the length-based indicator workflow. For each selected
#' gear and year, the function overlays the observed length-frequency
#' distribution with the expected vulnerable length distribution at Fspr, and
#' shows the Lspr threshold used in the indicator calculation.
#'
#' Panels are arranged with years in rows and gears in columns.
#'
#' @param fit fitted FLicc object
#' @param year optional years to plot; defaults to all available years
#' @param gear optional gears to plot; defaults to all available gears
#' @param spr target SPR percentage, default 40
#' @param thresh cumulative threshold used to define Lspr, default 0.75
#' @param nyears number of terminal years used when computing Fspr
#' @param scale_sel logical; passed to nf_flicc()
#' @param observed_fill fill colour for observed bars
#' @param observed_colour outline colour for observed bars
#' @param expected_colour line colour for expected Fspr distribution
#' @param threshold_colour colour of vertical threshold line
#' @param threshold_linetype linetype of threshold line
#' @param threshold_linewidth linewidth of threshold line
#' @param free_y logical; if TRUE, panels use free y scales
#' @param title optional plot title
#'
#' @return A ggplot object
#' @export
plot_LBIp <- function(fit,
                      year = NULL,
                      gear = NULL,
                      spr = 40,
                      thresh = 0.75,
                      nyears = 1,
                      scale_sel = TRUE,
                      observed_fill = "grey85",
                      observed_colour = "black",
                      expected_colour = "blue",
                      threshold_colour = "black",
                      threshold_linetype = 2,
                      threshold_linewidth = 0.8,
                      free_y = FALSE,
                      title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  `%||%` <- function(x, y) if (is.null(x)) y else x

  if (is.null(fit$report$obslen) || is.null(fit$report$sel_gear)) {
    stop("fit$report$obslen and fit$report$sel_gear are required.")
  }

  gear_names <- names(fit$report$obslen)
  if (is.null(gear)) {
    gear <- gear_names
  }
  gear <- as.character(gear)

  if (!all(gear %in% gear_names)) {
    stop("Unknown gear name(s): ", paste(setdiff(gear, gear_names), collapse = ", "))
  }

  # preserve requested gear order in facet columns
  gear <- unique(gear)

  # helper: pull a length-vector from array / matrix / FLQuant-like object
  get_slice <- function(x, yi = 1L) {
    dx <- dim(x)
    if (is.null(dx)) {
      return(as.numeric(x))
    }
    if (length(dx) == 1L) {
      return(as.numeric(x))
    }
    if (length(dx) == 2L) {
      yi <- max(1L, min(yi, dx[2]))
      return(as.numeric(x[, yi, drop = TRUE]))
    }
    if (length(dx) >= 6L) {
      yi <- max(1L, min(yi, dx[2]))
      return(as.numeric(x[, yi, 1, 1, 1, 1, drop = TRUE]))
    }
    idx <- rep(list(1), length(dx))
    idx[[1]] <- seq_len(dx[1])
    if (length(dx) >= 2L) {
      yi <- max(1L, min(yi, dx[2]))
      idx[[2]] <- yi
    }
    return(as.numeric(do.call(`[`, c(list(x), idx, list(drop = TRUE)))))
  }

  yrs_all <- dimnames(fit$report$obslen[[gear[1]]])$year
  if (is.null(yrs_all)) {
    dx_obs <- dim(fit$report$obslen[[gear[1]]])
    yrs_all <- if (!is.null(dx_obs) && length(dx_obs) >= 2L) seq_len(dx_obs[2]) else 1L
  }

  if (is.null(year)) {
    year <- yrs_all
  }
  year <- as.character(year)

  if (!all(year %in% as.character(yrs_all))) {
    stop("Unknown year(s): ", paste(setdiff(year, as.character(yrs_all)), collapse = ", "))
  }
  year_idx <- match(year, as.character(yrs_all))

  Ftgt <- fspr_flicc(fit, spr = spr, nyears = nyears, input = "F")
  Nref <- nf_flicc(fit, nyears = nyears, F = Ftgt, scale_sel = scale_sel)

  len_chr <- dimnames(Nref)$len %||% dimnames(fit$report$obslen[[gear[1]]])$len
  if (is.null(len_chr)) {
    dx <- dim(fit$report$obslen[[gear[1]]])
    len_chr <- seq_len(dx[1])
  }
  len_num <- as.numeric(len_chr)

  yrs_ref <- dimnames(Nref)$year
  if (is.null(yrs_ref)) {
    dx_ref <- dim(Nref)
    yrs_ref <- if (!is.null(dx_ref) && length(dx_ref) >= 2L) seq_len(dx_ref[2]) else 1L
  }

  out <- vector("list", length(gear) * length(year))
  k <- 1L

  for (g in gear) {
    obs_g_all <- fit$report$obslen[[g]]
    sel_g_all <- fit$report$sel_gear[[g]]

    yrs_sel <- dimnames(sel_g_all)$year
    if (is.null(yrs_sel)) {
      dx_sel <- dim(sel_g_all)
      yrs_sel <- if (!is.null(dx_sel) && length(dx_sel) >= 2L) seq_len(dx_sel[2]) else 1L
    }

    for (yy in seq_along(year)) {
      ylab <- year[yy]
      yi_obs <- year_idx[yy]

      yi_sel <- match(ylab, as.character(yrs_sel))
      if (is.na(yi_sel)) yi_sel <- min(yy, length(yrs_sel))

      yi_ref <- match(ylab, as.character(yrs_ref))
      if (is.na(yi_ref)) yi_ref <- min(yy, length(yrs_ref))

      obs_g  <- get_slice(obs_g_all, yi_obs)
      sel_g  <- get_slice(sel_g_all, yi_sel)
      nref_y <- get_slice(Nref, yi_ref)

      if (length(obs_g) != length(len_num) ||
          length(sel_g) != length(len_num) ||
          length(nref_y) != length(len_num)) {
        stop("Length mismatch in year ", ylab, ", gear ", g,
             ". Check dimensions of obslen, sel_gear, and nf_flicc output.")
      }

      if (sum(obs_g, na.rm = TRUE) <= 0 || sum(nref_y * sel_g, na.rm = TRUE) <= 0) {
        next
      }

      exp_vuln <- nref_y * sel_g
      obs_prop <- obs_g / sum(obs_g, na.rm = TRUE)
      exp_prop <- exp_vuln / sum(exp_vuln, na.rm = TRUE)

      keep <- is.finite(len_num) & is.finite(obs_prop) & is.finite(exp_prop)
      len2 <- len_num[keep]
      obs2 <- obs_prop[keep]
      exp2 <- exp_prop[keep]

      if (length(len2) < 2L) next

      len_thr <- len2[-1]
      exp_thr <- exp2[-1]
      cums <- cumsum(exp_thr)
      target <- sum(exp_thr, na.rm = TRUE) * thresh
      Lref <- len_thr[which.min((target - cums)^2)]

      p_obs_above <- sum(obs2[len2 >= Lref], na.rm = TRUE)
      p_exp_above <- sum(exp2[len2 >= Lref], na.rm = TRUE)

      out[[k]] <- data.frame(
        gear = g,
        year = ylab,
        len = len2,
        observed = obs2,
        expected = exp2,
        Lref = Lref,
        pobs = p_obs_above,
        pref = p_exp_above,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  out <- Filter(Negate(is.null), out)
  if (!length(out)) {
    stop("No valid data available for the requested gears/years.")
  }

  dat <- do.call(rbind, out)
  dat$gear <- factor(dat$gear, levels = gear)
  dat$year <- factor(dat$year, levels = year)

  vdat <- unique(dat[c("gear", "year", "Lref", "pobs", "pref")])
  vdat$label <- sprintf("Lspr = %.1f\nObs = %.2f\nExp = %.2f",
                        vdat$Lref, vdat$pobs, vdat$pref)

  ytop <- aggregate(pmax(observed, expected) ~ gear + year, data = dat, FUN = max)
  names(ytop)[3] <- "ymax"
  vdat <- merge(vdat, ytop, by = c("gear", "year"), all.x = TRUE, sort = FALSE)
  vdat$gear <- factor(vdat$gear, levels = gear)
  vdat$year <- factor(vdat$year, levels = year)
  vdat$ypos <- vdat$ymax * 0.95

  w <- if (length(unique(dat$len)) > 1L) {
    stats::median(diff(sort(unique(dat$len)))) * 0.9
  } else {
    0.9
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = len)) +
    ggplot2::geom_col(
      ggplot2::aes(y = observed),
      fill = observed_fill,
      colour = observed_colour,
      width = w
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = expected),
      colour = expected_colour,
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      data = vdat,
      ggplot2::aes(xintercept = Lref),
      colour = threshold_colour,
      linetype = threshold_linetype,
      linewidth = threshold_linewidth,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = vdat,
      ggplot2::aes(x = Lref, y = ypos, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1,
      size = 3
    ) +
    ggplot2::facet_grid(
      year ~ gear,
      scales = if (free_y) "free_y" else "fixed"
    ) +
    ggplot2::labs(
      x = "Length",
      y = "Proportion",
      title = if (is.null(title)) {
        paste0("Observed vs expected length composition at Fspr", spr)
      } else {
        title
      },
      subtitle = paste0(
        "Expected curve from nf_flicc() at SPR", spr,
        "; dashed line shows Lspr threshold (thresh = ", thresh, ")"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95"),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}
