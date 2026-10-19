\name{sel_curve_flicc}
\alias{sel_curve_flicc}
\title{Selectivity at length for given parameters}
\usage{
sel_curve_flicc(L, sel_fun = c("logistic", "dsnormal", "normal", "binormal", "bidsnormal"), pars)
}
\arguments{
\item{L}{Lengths (e.g. bin mid-points).}
\item{sel_fun}{\code{"logistic"} (\code{SL50}, \code{SL95}), \code{"dsnormal"} (\code{mode}, \code{lsd}, \code{rsd}), \code{"normal"} (\code{mode}, \code{sd}) \code{"binormal"} (\code{mode1}, \code{sd1}, \code{mode2}, \code{sd2}, \code{h}) or \code{"bidsnormal"} (\code{mode1}, \code{sd1}, \code{mode2}, \code{lsd2}, \code{rsd2}, \code{h}); the two-mode functions are scaled to a maximum of 1 over \code{L}.}
\item{pars}{Named list or vector of parameters on the length scale (as in \code{sel_fixed}).}
}
\value{Numeric vector of selectivity.}
\description{Evaluates the FLicc selectivity functions with parameters on the natural length scale.}
\examples{
sel_curve_flicc(seq(10, 60, 5), "dsnormal", list(mode = 30, lsd = 4, rsd = 15))
}
