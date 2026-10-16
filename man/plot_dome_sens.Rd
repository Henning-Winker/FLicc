\name{plot_dome_sens}
\alias{plot_dome_sens}
\title{Plot a dome-shaped selectivity sensitivity}
\usage{
plot_dome_sens(x, what = c("sel", "spr"), ...)
}
\arguments{
\item{x}{Output of \code{dome_sens_flicc()}.}
\item{what}{\code{"sel"}: the varied gear's fitted selectivity per scenario; \code{"spr"}: SPR trajectories via \code{plot_spr()}.}
\item{...}{Passed to \code{plot_spr()} (e.g. \code{ci = TRUE}).}
}
\value{A \code{ggplot} object.}
\description{Plots the fitted selectivity or SPR of each scenario of \code{dome_sens_flicc()}.}
