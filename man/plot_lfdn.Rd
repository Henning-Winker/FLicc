% Generated manually (roxygen2 style)
\name{plot_lfdn}
\alias{plot_lfdn}
\title{Sample size of length-frequency data by year and gear}
\usage{
plot_lfdn(lfd, ess = NULL, bins = TRUE)
}
\arguments{
\item{lfd}{\code{FLQuants} of length frequencies by gear (raw counts).}
\item{ess}{Optional \code{FLQuants} after \code{lfdess()}, or a scalar / named vector of effective sample sizes per gear (passed to \code{lfdess()}).}
\item{bins}{Logical; label bars with the number of non-empty length bins.}
}
\value{A \code{ggplot} object.}
\description{
Number of fish measured per year and gear, stacked so the bar height is the combined sample size, optionally next to the effective sample sizes after \code{lfdess()}. Gear-years without samples show as gaps; labels give the number of non-empty bins.
}
\examples{
\dontrun{
plot_lfdn(lfd)
plot_lfdn(lfd, ess = c(Trawl = 200, Gillnet = 150))
}
}
