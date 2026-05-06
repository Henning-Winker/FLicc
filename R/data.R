
#' Example alfonsino data
#'
#' Example data objects for alfonsino, including life-history parameters,
#' length-frequency data, and an example \code{FLStockLen} object.
#'
#' @name alfonsino
#' @docType data
#' @usage data(alfonsino)
#' @format The dataset loads three objects:
#' \describe{
#'   \item{\code{lhpar_alfonsino}}{An \code{FLPar} object with life-history parameters.}
#'   \item{\code{lfd_alfonsino}}{An \code{FLQuants} object with simulated length-frequency data.}
#'   \item{\code{stklen_alfonsino}}{An \code{FLStockLen} example object.}
#' }
#' @source FLife simulation
#' @author Henning Winker
NULL


#' Example multi-gear length-frequency data from fishblicc (Bombay duck)
#'
#' Example dataset adapted from the \pkg{fishblicc} package, representing
#' a multi-gear artisanal fishery targeting Bombay duck
#' (\emph{Harpadon nehereus}). The data consist of length-frequency
#' distributions and associated life-history parameters, formatted for use
#' in \pkg{FLicc}.
#'
#' The dataset includes length-frequency data for three artisanal gears.
#' Note that this is a simplified example and does not include trawl catches,
#' so it should be treated as illustrative rather than a complete assessment
#' dataset.
#'
#' @name fishblicc_example
#'
#' @docType data
#'
#' @format A list containing:
#' \describe{
#'   \item{lfd_fishblicc}{An \code{FLQuants} object of length-frequency data
#'   by gear (counts-at-length).}
#'   \item{lhpar_fishblicc}{An \code{FLPar} object containing life-history
#'   parameters (e.g. \eqn{L_\infty}, \eqn{k}, maturity, weight-length, etc.).}
#' }
#'
#' @details
#' The data originate from the example provided in the
#' \pkg{fishblicc} repository (Paul Medley), illustrating a
#' multi-gear length-based catch-curve analysis.
#'
#' The dataset has been reformatted to:
#' \itemize{
#'   \item conform to \pkg{FLR} object classes (\code{FLQuants}, \code{FLPar})
#'   \item support direct use in \code{fit_flicc()}
#'   \item provide a reproducible example for testing and demonstration
#' }
#'
#' Users should note that:
#' \itemize{
#'   \item the dataset is incomplete (missing trawl component)
#'   \item results should not be interpreted as a full stock assessment
#'   \item the example is intended for method illustration only
#' }
#'
#' @source
#' Adapted from:
#' \url{https://github.com/PaulAHMedley/fishblicc}
#'
#' Original description:
#' Bombay duck (\emph{Harpadon nehereus}) multi-gear artisanal fishery.
#'
#' @references
#' Medley, P. A. (2025). A Bayesian catch-curve stock assessment model for
#' the analysis of length data from multi-gear fisheries.
#'
#' @examples
#' data(fishblicc_example)
#'
#' stklen <- stocklen(
#'   lfd_fishblicc,
#'   lhpar_fishblicc,
#'   m_model = "inverse"
#' )
#'
#' \dontrun{
#' fit <- fiticc(
#'   lfd_fishblicc,
#'   stklen,
#'   sel_fun = c("dsnormal", "dsnormal", "dsnormal"),
#'   catch_by_gear = c(0.1802070, 0.2101353, 0.6096577),
#'   settings = list(
#'     pop_model = "gamma",
#'     obs_model = "nb",
#'     CVL = 0.15,
#'     CVL.sd = 0.1,
#'     linf.sd = 2 / 40,
#'     Mk.sd = 0.1
#'   )
#' )
#' }
#'
NULL


#' Example FLicc Monte Carlo ensemble
#'
#' Example Monte Carlo ensemble object generated from a fitted FLicc model.
#' The object is used to illustrate structural uncertainty diagnostics,
#' likelihood/AIC weighting, and SPR ensemble summaries.
#'
#' The ensemble contains Monte Carlo refits based on alternative biologically
#' plausible life-history parameter draws. These draws are typically generated
#' from life-history priors and correlations, then refitted with \code{fiticc()}.
#' The object also includes the reference fit quantities in \code{mc$fit}, which
#' allows comparison between the original fitted model and the retained or
#' likelihood-weighted ensemble.
#'
#' @format A list containing, among others:
#' \describe{
#'   \item{\code{lhpar}}{Monte Carlo life-history parameter draws, typically
#'   including \code{linf}, \code{k}, \code{M}, \code{Mk}, \code{L50}, and
#'   \code{CVL}.}
#'   \item{\code{spr}}{SPR estimates from each Monte Carlo refit by year.}
#'   \item{\code{ll}}{Log-likelihood or negative log-likelihood information
#'   from the Monte Carlo refits, used to compute delta AIC and Akaike weights.}
#'   \item{\code{fit}}{Reference fit quantities extracted from the original
#'   FLicc fit, including life-history parameters and SPR.}
#' }
#'
#' @details
#' This data set is intended for examples and package testing. It supports
#' functions such as \code{plot_mcpars()}, \code{plot_mccor()},
#' \code{plot_mcprofile()}, \code{plot_mcsprdist()},
#' \code{plot_mcsprwts()}, and \code{mc_spr_tab()}.
#'
#' @examples
#' data(flicc_mc_example)
#'
#' plot_mcpars(flicc_mc_example)
#' plot_mccor(flicc_mc_example)
#' plot_mcprofile(flicc_mc_example)
#' plot_mcsprdist(flicc_mc_example)
#' plot_mcsprwts(flicc_mc_example)
#'
#' mc_spr_tab(flicc_mc_example)
#'
#' @name flicc_mc_example
#' @docType data
#' @keywords datasets
NULL

