% Please edit documentation in R/flicc_plots.R
\name{epr_ci_flicc}
\alias{epr_ci_flicc}
\title{Confidence intervals of the exploitable potential ratio}
\usage{
epr_ci_flicc(fit, level = 0.95)
}
\arguments{
\item{fit}{A fitted FLicc object (with standard errors).}

\item{level}{Confidence level.}
}
\value{
A \code{data.frame} with year, epr, se_log, lwr and upr, or \code{NULL} without standard errors.
}
\description{
As \code{spr_ci_flicc()} for EPR: intervals from the standard errors of log EPR
(\code{TMB::sdreport()}), back-transformed.
}
