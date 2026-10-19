% Please edit documentation in R/flicc_plots.R
\name{plot_epr}
\alias{plot_epr}
\title{Plot the exploitable potential ratio (EPR)}
\usage{
plot_epr(
  fit,
  title = NULL,
  by_year = 2,
  line_colour = "black",
  point_colour = "black",
  linewidth = 0.6,
  point_size = 2,
  ymax = NULL,
  reflines = TRUE,
  ref_levels = c(0.4, 0.2),
  ref_colours = c("darkgreen", "red"),
  ref_linetype = 2,
  ci = FALSE,
  ci_level = 0.95,
  ci_type = c("bars", "ribbon")
)
}
\arguments{
\item{fit}{A fitted FLicc object, a named list of fits, or \code{FLQuants}.}

\item{title}{Optional plot title.}

\item{by_year}{Spacing of the year axis breaks.}

\item{line_colour, point_colour}{Colours of the line and points (single fit).}

\item{linewidth, point_size}{Line width and point size.}

\item{ymax}{Optional upper y-axis limit.}

\item{reflines}{Logical; draw reference lines.}

\item{ref_levels}{EPR reference levels (default 0.4 and 0.2). For LBB's reference use \code{attr(lbb$table, "BFM1B0")} (B/B0 at F = M and Lc = Lc_opt) and half of it.}

\item{ref_colours}{Colours of the reference lines.}

\item{ref_linetype}{Line type of the reference lines.}

\item{ci}{Logical; add confidence intervals from the standard errors of log EPR.}

\item{ci_level}{Confidence level.}

\item{ci_type}{\code{"bars"} or \code{"ribbon"}.}
}
\value{
A \code{ggplot} object.
}
\description{
EPR is the biomass vulnerable to the fishery (length weights times the combined
selectivity of all gears, i.e. F at length relative to apical F) per recruit at
the year's F, relative to the unfished state. For one gear it is the quantity
LBB reports as B/B0 (Froese et al. 2018). Same features as \code{plot_spr()}.
}
\examples{
\dontrun{
plot_epr(fit, ci = TRUE)
plot_epr(list(base = fit, cut = fit_cut), ci = TRUE, ci_type = "ribbon")
}
}
