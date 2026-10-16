\name{dome_sens_flicc}
\alias{dome_sens_flicc}
\title{Dome-shaped selectivity sensitivity}
\usage{
dome_sens_flicc(fit, gear, rsd = NULL, fixed = NULL, se = TRUE, ...)
}
\arguments{
\item{fit}{A fitted FLicc object (from \code{fiticc()}, with \code{FLRreport = TRUE}).}
\item{gear}{Gear name (or index) whose selectivity is varied.}
\item{rsd}{Right-hand standard deviations (length units). Default \code{c(1000, 0.6, 0.4, 0.25) * Linf}; 1000 gives a flat right limb (no dome).}
\item{fixed}{Optional further fixed parameters for this gear, e.g. \code{list(lsd = 4)}.}
\item{se}{Passed to \code{fiticc()}.}
\item{...}{Further arguments passed to \code{fiticc()}.}
}
\value{A list with \code{fits} (\code{base} = \code{fit}, then one fit per \code{rsd}), \code{table} (scenario, rsd, convergence, log-likelihood, AIC, dAIC, last-year SPR, Fspr, last-year apical F), \code{gear}, \code{rsd} and \code{errors} (messages of refits that failed).}
\description{Refits \code{fit} with one gear's selectivity as a double-sided normal whose right-hand standard deviation is fixed at a sequence of values, from no doming to mild and moderate doming; mode and left-hand SD are estimated, all other settings and data as in \code{fit}. Recommended as a standard sensitivity for gears assumed logistic.}
\examples{
\dontrun{
ds <- dome_sens_flicc(fit, gear = "Gillnet")
ds$table
plot_dome_sens(ds)
plot_dome_sens(ds, "spr", ci = TRUE)
}
}
