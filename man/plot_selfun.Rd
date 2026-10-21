\name{plot_selfun}
\alias{plot_selfun}
\title{Explore selectivity curves over a range of parameter values}
\usage{
plot_selfun(sel_fun = c("dsnormal", "logistic", "normal", "binormal", "bidsnormal"), pars, len = NULL,
  lfd = NULL, ref = NULL)
}
\arguments{
\item{sel_fun}{\code{"logistic"} (\code{SL50}, \code{SL95}), \code{"dsnormal"} (\code{mode}, \code{lsd},
\code{rsd}), \code{"normal"} (\code{mode}, \code{sd}), \code{"binormal"} (\code{mode1}, \code{sd1}, \code{mode2},
\code{sd2}, \code{h}) or \code{"bidsnormal"} (\code{mode1}, \code{sd1}, \code{mode2}, \code{lsd2}, \code{rsd2},
\code{h}); see \code{\link[=fiticc]{fiticc()}}.}
\item{pars}{Named list of parameters on the length scale; elements with several values are crossed and each combination is drawn.}
\item{len}{Lengths at which to evaluate; default from \code{lfd} or \code{0:80}.}
\item{lfd}{Optional \code{FLQuant} (one gear) or \code{FLQuants} (pooled) of length frequencies, drawn as grey bars scaled to a maximum of 1.}
\item{ref}{Optional named numeric of reference lengths drawn as vertical lines, e.g. \code{c(L50 = 31, Lmean = 36)}.}
}
\value{A \code{ggplot} object.}
\description{Plots selectivity for one or more parameter values (e.g. a sequence of right-hand standard deviations of a dome, from 1000 = no dome to moderate doming), optionally over the observed length composition. Useful for choosing \code{sel_fixed} values and the scenarios of \code{dome_sens_flicc()}.}
\examples{
\dontrun{
plot_selfun("dsnormal", list(mode = 32, lsd = 4, rsd = c(1000, 35, 20, 12)),
            lfd = lfd$Gillnet, ref = c(L50 = 31.2, Linf = 55.7))
plot_selfun("logistic", list(SL50 = c(26, 30, 34), SL95 = 38))
}
# two size modes (juveniles at 10 cm, adults at 25 cm), width of the adult curve
plot_selfun("binormal", list(mode1 = 10, sd1 = 1, mode2 = 25, sd2 = c(4, 8), h = 0.3))
# second mode from flat (rsd2 = 1000) to domed, and its relative height
plot_selfun("bidsnormal", list(mode1 = 10, sd1 = 1, mode2 = 25, lsd2 = 5,
                               rsd2 = c(1000, 15, 8), h = c(0.3, 1)))
}
