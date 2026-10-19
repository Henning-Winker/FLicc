\name{sel_curve_flicc}
\alias{sel_curve_flicc}
\title{Selectivity at length for given parameters}
\usage{
sel_curve_flicc(L, sel_fun = c("logistic", "dsnormal", "normal"), pars)
}
\arguments{
\item{L}{Lengths (e.g. bin mid-points).}
\item{sel_fun}{\code{"logistic"} (\code{SL50}, \code{SL95}), \code{"dsnormal"} (\code{mode}, \code{lsd}, \code{rsd}) or \code{"normal"} (\code{mode}, \code{sd}).}
\item{pars}{Named list or vector of parameters on the length scale (as in \code{sel_fixed}).}
}
\value{Numeric vector of selectivity.}
\description{Evaluates the FLicc selectivity functions with parameters on the natural length scale.}
\examples{
sel_curve_flicc(seq(10, 60, 5), "dsnormal", list(mode = 30, lsd = 4, rsd = 15))
}
