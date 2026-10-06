#' Build an FLStockLen BRP object from an FLicc fit
#'
#' Constructs an \code{FLStockLen}-like object over a grid of fishing levels.
#' The BRP grid is stored on the year dimension, with pseudo-years
#' \code{1:length(Fseq)}.
#'
#' Slots are filled with equilibrium quantities at length:
#' \itemize{
#'   \item \code{stock.n}: equilibrium numbers-at-length
#'   \item \code{harvest}: fishing mortality-at-length
#'   \item \code{catch.n}: equilibrium catch-at-length in numbers
#'   \item \code{stock.wt}, \code{catch.wt}: weight-at-length
#'   \item \code{mat}: maturity-at-length
#'   \item \code{m}: natural mortality-at-length
#' }
#'
#' Optional Beverton-Holt steepness \code{s} is stored as an attribute and used
#' later by helper functions for equilibrium recruitment and biomass ratios.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param year Assessment year used as the biological reference. Defaults to
#'   the terminal year.
#' @param Fseq Numeric vector of fishing levels.
#' @param input Character string, one of \code{"FM"} or \code{"F"}, indicating
#'   the scale of \code{Fseq}.
#' @param Sel Optional numeric selectivity vector to override the default joint
#'   selectivity from the fit.
#' @param scale_sel Logical. If \code{TRUE}, the default joint selectivity is
#'   scaled to a maximum of 1.
#' @param spawn_time Fraction of the interval elapsed before spawning.
#' @param s Optional Beverton-Holt steepness.
#'
#' @return An object of class \code{FLStockLen}, with attributes
#'   \code{Fseq}, \code{input}, \code{steepness}, \code{spawn_time},
#'   and \code{brp_year}.
#'
#' @export
eqstklen <- function(fit, nyears = 1,
                             Fseq =NULL,
                             Sel = NULL, scale_sel = TRUE,
                             spawn_time = 0,
                             s = 0.7,spr.tgt=40) {

  fcur <- fcur_flicc(fit)
  if(is.null(Fseq)) Fseq = seq(0,max(an(fit$report$lhpar["M"])*5,fcur*1.05) , length.out = 101)

  prbrp <- prbrp_flicc(fit, nyears = nyears, Fseq = Fseq,
        scale_sel = TRUE, spawn_time = spawn_time)


  year <- dimnames(fit$report$N)$year
  lens <- dimnames(fit$report$N)$len
  nlen <- length(lens)
  nyear <- length(year)
  yrs <- tail(year,min(nyear,nyears))
  SPR <- FLQuant(prbrp$SBPR/prbrp$spr0,quant="len",unit="")
  if(!is.null(s)){
  recs <- rr0_bh(SPR,s=s)
  } else {
  recs <- rr0_bh(SPR,s=0.7)
  recs[] <- 1
  }
  # template from original stklen to preserve class/range structure
  stk <- FLCore::FLStockLen(stock.n=prbrp$NPR%*%recs,
                    catch.n=prbrp$CPR%*%recs,
                    landings.n =prbrp$CPR,
                    stock.wt = prbrp$wl,
                    catch.wt = prbrp$wl,
                    landings.wt = prbrp$wl,
                    discards.wt = prbrp$wl,
                    harvest = prbrp$FPR,
                    m = prbrp$m,
                    mat=prbrp$mat,
                    stock =   SPR
                    )

  FLCore::discards.n(stk)[] <- 0
  FLCore::m.spwn(stk)[] <- spawn_time
  FLCore::harvest.spwn(stk)[] <- spawn_time
  FLCore::catch(stk) <- FLCore::computeCatch(stk)
  FLCore::landings(stk) <- FLCore::computeLandings(stk)
  FLCore::discards(stk) <- FLCore::computeDiscards(stk)
  #------ Calculate refpts-------#


  Fcur <- fcur_flicc(fit,nyears = nyears)
  Fspr = fspr_flicc(fit,input="F",spr=spr.tgt)
  Ycur = yeql(fit,F=Fcur,s=s,nyears=nyears)/max(catch(stk))
  SSBcur = ssbeql(fit,F=Fcur,s=s,nyears=nyears)/prbrp$spr0
  SPRcur = spreql(fit,F=Fcur,s=s,nyears=nyears)
  # Life history
  lhpar <- fit$report$lhpar
  # refpts
  refpts <- FLPar(
  Fcur =  Fcur ,
  Fspr =  Fspr,
  Ycur =  Ycur,
  SSBcur = SSBcur ,
  SPRcur =  SPRcur ,
  SPR0 = prbrp$spr0
  )

  if("s"%in%rownames(lhpar)){
    lhpar["s"] <- s
  }
  attr(stk,"lhpar") <- lhpar
  attr(stk,"selpars") <-fit$report$selpars
  attr(stk,"pars") <- fit$report$pars
  attr(stk,"refpts") <-   refpts
  return(stk)
}

#' Extract BRP plotting quantities from an FLStockLen BRP object
#'
#' @param stk An object returned by \code{brp_stklen_flicc()}.
#'
#' @return An \code{FLQuants} object with \code{F}, \code{Rec}, \code{Catch},
#'   and \code{SSB}.
#' @export
brp_flqs_flicc <- function(stk) {
  Rec <- recl(stk)
  SSB <- ssbl(stk)
  Fq  <- fapl(stk)
  Cq  <- FLCore::catch(stk)

  FLCore::FLQuants(F = Fq, Rec = Rec, Catch = Cq, SSB = SSB)
}

#' SPR-style reference points from an FLStockLen BRP object
#'
#' @param stk An object returned by \code{brp_stklen_flicc()}.
#' @param sprx Numeric vector of SPR reference levels as proportions.
#' @param labels Character labels for the reference points.
#'
#' @return An \code{FLPar} with SPR-based fishing reference points.
#' @export
refpts_flicc <- function(stk,fit=0,
                         spr = c(40, 20, 10),
                         labels = paste0("Fspr",spr )) {

  if (length(spr) != length(labels))
    stop("'spr' and 'labels' must have the same length.")

  Fseq <- an(fapl(stk))
  SSB  <- an(ssbl(stk))
  SPR  <- an(SSB / an(SSB[1]))

  if(!is.null(fit)){
  rps <- do.call(c,lapply(spr, function(x){
    fspr_flicc(fit,spr=x,input="F")}))
  }

  out <- FLCore::FLPar(setNames(rps, labels))

  out
}



#' Build per-recruit BRP object from FLicc fit
#'
#' Computes per-recruit quantities over a sequence of fishing mortalities
#' using \code{pr_flicc()}, and stores the results in an FLR-friendly list
#' with \code{FLQuant} components. This object serves as the core input
#' for equilibrium (stock-recruit) calculations and the construction of
#' \code{FLStockLen} objects.
#'
#' The function evaluates per-recruit quantities once over \code{Fseq} and
#' returns:
#' \itemize{
#'   \item Yield-per-recruit (\code{YPR})
#'   \item Spawning biomass-per-recruit (\code{SBPR})
#'   \item Numbers-per-recruit at length (\code{NPR})
#'   \item Catch-per-recruit at length (\code{CPR})
#'   \item Fishing mortality-at-length (\code{FPR})
#' }
#'
#' These outputs can then be combined with a stock-recruit relationship
#' (e.g. Beverton–Holt) to compute equilibrium quantities such as
#' yield, SSB, and recruitment.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param nyears Integer. Number of terminal years used to construct
#'   reference parameters (default = 1).
#' @param Fseq Numeric vector of fishing mortality values on the
#'   apical (annual F) scale. If \code{NULL}, a default sequence from
#'   0 to \code{10 * M} is used.
#' @param scale_sel Logical. Passed to \code{pr_flicc()} to control
#'   selectivity scaling (default = TRUE).
#' @param spawn_time Numeric. Spawning time as fraction of the year,
#'   passed to \code{pr_flicc()} (default = 0).
#'
#' @return A list containing per-recruit quantities:
#' \describe{
#'   \item{Fseq}{\code{FLQuant}. Fishing mortality sequence.}
#'   \item{Len}{Numeric vector of length classes.}
#'   \item{years}{Reference year used for parameters.}
#'   \item{YPR}{\code{FLQuant}. Yield per recruit.}
#'   \item{SBPR}{\code{FLQuant}. Spawning biomass per recruit.}
#'   \item{NPR}{\code{FLQuant}. Numbers per recruit at length.}
#'   \item{CPR}{\code{FLQuant}. Catch per recruit at length.}
#'   \item{FPR}{\code{FLQuant}. Fishing mortality at length.}
#'   \item{wl}{\code{FLQuant}. Weight-at-length (kg).}
#'   \item{mat}{\code{FLQuant}. Maturity-at-length.}
#'   \item{m}{\code{FLQuant}. Natural mortality-at-length.}
#'   \item{M}{Numeric. Scalar natural mortality.}
#'   \item{k}{Numeric. Growth parameter used for scaling.}
#'   \item{spr0}{Unfished spawning biomass per recruit.}
#'   \item{scale_sel}{Logical. Selectivity scaling flag used.}
#'   \item{spawn_time}{Numeric. Spawning timing used.}
#' }
#'
#' @details
#' This function separates the per-recruit calculations from stock-recruit
#' dynamics. It is intended to be used as input to higher-level functions
#' such as:
#' \itemize{
#'   \item \code{yieldrp()} – equilibrium yield
#'   \item \code{ssbrp()} – equilibrium spawning biomass
#'   \item \code{fmsy_flicc()} – MSY reference point estimation
#'   \item \code{brp_stklen_flicc()} – construction of equilibrium stock objects
#' }
#'
#' @examples
#' \dontrun{
#' prb <- prbrp_flicc(fit)
#'
#' # plot SPR curve
#' SPR <- prb$SBPR / prb$spr0
#' par(mfrow=c(1,2))
#' plot(as.numeric(prb$Fseq), as.numeric(SPR), type = "l",ylab="SPR",xlab="F")
#'
#' # check YPR
#' plot(as.numeric(prb$Fseq), as.numeric(prb$YPR), type = "l",ylab="YPR",xlab="F")
#' }
#'
#' @export
prbrp_flicc <- function(fit, nyears = 1, Fseq = NULL,
                        scale_sel = TRUE, spawn_time = 0) {

  if (is.null(Fseq)) {
    Fseq <- seq(0, as.numeric(fit$report$lhpar["M"]) * 5, length.out = 101)
  }


  ref <- flicc_refpars(fit, nyears = nyears, scale_sel = scale_sel)
  lens = dimnames(fit$report$sel)$len
  nlen <- length(ref$Len)

  # one compiled TMB evaluation for all F values (see brp_tmb_flicc())
  b <- if (brp_fast_ok(fit, NULL, scale_sel)) {
    brp_tmb_flicc(fit, Fseq = Fseq, nyears = nyears, spawn_time = spawn_time)
  } else NULL

  if (!is.null(b)) {
    spr0 <- b$SBPR0
    YPR  <- FLQuant(b$YPR,  quant = "all", unit = "kg")
    SBPR <- FLQuant(b$SBPR, quant = "all", unit = "kg")
    nf <- unname(b$N); cf <- unname(b$C); ff <- unname(b$Fl)
  } else {

  pr_list <- lapply(Fseq, function(ff) {
    pr_flicc(
      fit = fit,
      nyears = nyears,
      F = ff,
      scale_sel = scale_sel,
      spawn_time = spawn_time
    )
  })

  spr0 <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = 0,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )$SBPR

  YPR  <- FLQuant(vapply(pr_list, function(x) as.numeric(x$YPR), numeric(1)), quant = "all",unit="kg")
  SBPR <- FLQuant(vapply(pr_list, function(x) as.numeric(x$SBPR), numeric(1)), quant = "all",unit="kg")

  nf <- array(NA_real_, dim = c(nlen, length(Fseq)))
  cf <- array(NA_real_, dim = c(nlen, length(Fseq)))
  ff <- array(NA_real_, dim = c(nlen, length(Fseq)))

  for (i in seq_along(pr_list)) {
    nf[, i] <- as.numeric(pr_list[[i]]$N)
    cf[, i] <- as.numeric(pr_list[[i]]$Cn)
    ff[, i] <- as.numeric(pr_list[[i]]$F)
  }
  }

  nf <- FLQuant(nf,quant="len",unit="1000",dimnames=list(len=lens))
  cf <- FLQuant(cf,quant="len",unit="1000",dimnames=list(len=lens))
  ff <- FLQuant(ff,quant="len",unit="f",dimnames=list(len=lens))
  wl <- FLQuant(nf, unit = "kg", dimnames = list(len = lens))
  wl[] <- ref$wl
  mat <- FLQuant(wl,unit="")
  mat[] <- mat(fit$stklen)[, ref$year]
  m <- FLQuant(wl,unit="m")
  m[] <- m(fit$stklen)[, ref$year]


  fap <- FLQuant(Fseq,quant="len",unit="f")
  list(
    Fseq = fap,
    Len = ref$Len,
    years = ref$year,
    YPR = YPR,
    SBPR = SBPR,
    NPR = nf,
    CPR = cf,
    FPR = ff,
    wl = wl,
    mat = mat,
    m = m,
    M = ref$M,
    k = ref$k,
    spr0 = spr0,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )
}

#' Beverton-Holt equilibrium recruitment ratio
#'
#' Computes equilibrium recruitment relative to unfished recruitment
#' (\eqn{R/R_0}) under a Beverton-Holt stock-recruitment relationship,
#' expressed as a function of spawning potential ratio (\code{SPR}) and
#' steepness (\code{s}).
#'
#' The function is intended as a lightweight helper for converting
#' per-recruit spawning biomass ratios into equilibrium recruitment ratios
#' when constructing stock-recruit based BRP quantities such as equilibrium
#' yield, spawning biomass, or stock objects.
#'
#' The Beverton-Holt equilibrium recruitment ratio is:
#' \deqn{
#' R/R_0 = \frac{4 s \cdot SPR - (1 - s)}{(5 s - 1)\cdot SPR}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{SPR = SBPR(F) / SBPR(0)}
#'   \item \eqn{s} is Beverton-Holt steepness
#' }
#'
#' Non-finite values are set to \code{NA}, and negative recruitment ratios
#' are truncated to zero.
#'
#' @param SPR Numeric vector, \code{FLQuant}, or object coercible to numeric,
#'   giving spawning potential ratio values.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#'
#' @return An \code{FLQuant} containing equilibrium recruitment ratio
#'   (\code{R/R0}). Values are dimensioned consistently with the input after
#'   coercion through \code{FLQuant()}.
#'
#' @details
#' This helper does not scale by absolute recruitment. It returns recruitment
#' relative to unfished recruitment, so it is typically combined with
#' per-recruit quantities such as:
#' \itemize{
#'   \item \code{yieldrp()} for equilibrium yield,
#'   \item \code{ssbrp()} for equilibrium spawning biomass,
#'   \item \code{brp_stklen_flicc()} for updating equilibrium stock numbers.
#' }
#'
#' For example, if \code{NPR} is numbers-per-recruit, equilibrium numbers are
#' obtained as:
#' \code{N_eq = NPR * rr0_bh(SPR, s)}.
#'
#' @examples
#' \dontrun{
#' # SPR curve from a PR object
#' prb <- prbrp_flicc(fit)
#' SPR <- prb$SBPR / prb$spr0
#'
#' # Beverton-Holt equilibrium recruitment ratio
#' RR0 <- rr0_bh(SPR, s = 0.7)
#'
#' plot(as.numeric(prb$Fseq), as.numeric(RR0), type = "l",
#'      xlab = "F", ylab = "R/R0")
#' }
#'
#' @export
rr0_bh <- function(SPR, s = 0.7) {
  rr0 <- (4 * s * SPR - (1 - s)) / ((5 * s - 1) * SPR)
  rr0[!is.finite(rr0)] <- NA_real_
  rr0[rr0 < 0] <- 0
  FLQuant(rr0, unit = "")
}
#' Equilibrium yield from a per-recruit BRP object
#'
#' Computes equilibrium yield as the product of yield-per-recruit
#' (\code{YPR}) and equilibrium recruitment relative to unfished recruitment
#' (\code{R/R0}) under a Beverton-Holt stock-recruit relationship.
#'
#' @param prbrp A per-recruit BRP object returned by \code{prbrp_flicc()}.
#'   It must contain components \code{YPR}, \code{SBPR}, and \code{spr0}.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#'
#' @return An \code{FLQuant} of equilibrium yield over the fishing mortality
#'   sequence stored in \code{prbrp}.
#'
#' @examples
#' \dontrun{
#' prb <- prbrp_flicc(fit)
#' Yeq <- yieldf(prb, s = 0.7)
#' plot(as.numeric(prb$Fseq), as.numeric(Yeq), type = "l")
#' }
#'
#' @export
yieldf <- function(prbrp, s = 0.7) {
  SPR <- prbrp$SBPR / prbrp$spr0
  RR0 <- rr0_bh(SPR, s = s)
  prbrp$YPR * RR0
}

#' Equilibrium spawning biomass from a per-recruit BRP object
#'
#' Computes equilibrium spawning biomass as the product of
#' spawning-biomass-per-recruit (\code{SBPR}) and equilibrium recruitment
#' relative to unfished recruitment (\code{R/R0}) under a Beverton-Holt
#' stock-recruit relationship.
#'
#' @param prbrp A per-recruit BRP object returned by \code{prbrp_flicc()}.
#'   It must contain components \code{SBPR} and \code{spr0}.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#'
#' @return An \code{FLQuant} of equilibrium spawning biomass over the fishing
#'   mortality sequence stored in \code{prbrp}.
#'
#' @examples
#' \dontrun{
#' prb <- prbrp_flicc(fit)
#' SSBeq <- ssbf(prb, s = 0.7)
#' plot(as.numeric(prb$Fseq), as.numeric(SSBeq), type = "l")
#' }
#'
#' @export
ssbf <- function(prbrp, s = 0.7) {
  SPR <- prbrp$SBPR / prbrp$spr0
  RR0 <- rr0_bh(SPR, s = s)
  prbrp$SBPR * RR0
}

#' Relative equilibrium spawning biomass from a per-recruit BRP object
#'
#' Computes equilibrium spawning biomass relative to unfished spawning biomass
#' (\code{SSB/SSB0}) under a Beverton-Holt stock-recruit relationship.
#'
#' This is calculated as the product of spawning potential ratio
#' (\code{SPR = SBPR / spr0}) and equilibrium recruitment ratio
#' (\code{R/R0}).
#'
#' @param prbrp A per-recruit BRP object returned by \code{prbrp_flicc()}.
#'   It must contain components \code{SBPR} and \code{spr0}.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#'
#' @return An \code{FLQuant} of equilibrium spawning biomass relative to
#'   unfished spawning biomass.
#'
#' @examples
#' \dontrun{
#' prb <- prbrp_flicc(fit)
#' SSBrel <- ssbrel(prb, s = 0.7)
#' plot(as.numeric(prb$Fseq), as.numeric(SSBrel), type = "l")
#' abline(h = c(0.4, 0.2, 0.1), lty = 2)
#' }
#'
#' @export
ssbrel <- function(prbrp, s = 0.7) {
  SPR <- prbrp$SBPR / prbrp$spr0
  RR0 <- rr0_bh(SPR, s = s)
  prbrp$SBPR / prbrp$spr0 * RR0
}


#' Estimate Fmsy from a per-recruit BRP object
#'
#' Estimates the fishing mortality corresponding to maximum equilibrium yield
#' (\eqn{F_{MSY}}) under a Beverton-Holt stock-recruit relationship.
#'
#' The function uses a per-recruit BRP object produced by
#' \code{prbrp_flicc()} and computes equilibrium yield from
#' \code{yieldf(prbrp, s = s)}. By default, \eqn{F_{MSY}} is approximated as
#' the fishing mortality in \code{prbrp$Fseq} that maximizes equilibrium yield.
#'
#' If the fitted \code{"flicc_tmb_fit"} object is also supplied through
#' \code{fit}, the function can optionally refine the estimate using a local
#' one-dimensional optimization around the discrete maximum.
#'
#' @param prbrp A per-recruit BRP object returned by \code{prbrp_flicc()}.
#'   It must contain at least \code{Fseq}, \code{YPR}, \code{SBPR}, and
#'   \code{spr0}.
#' @param fit Optional fitted \code{"flicc_tmb_fit"} object. If supplied,
#'   \eqn{F_{MSY}} can be refined using a local root/optimization search around
#'   the discrete maximum.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#' @param nyears Integer. Number of terminal years passed to the refining
#'   \code{pr_flicc()} call when \code{fit} is supplied.
#' @param scale_sel Logical. Passed to \code{pr_flicc()} during refinement.
#' @param spawn_time Numeric. Spawning timing passed to \code{pr_flicc()}
#'   during refinement.
#' @param refine Logical. If \code{TRUE} and \code{fit} is supplied,
#'   refine the discrete estimate using local numerical optimization.
#'   Default is \code{!is.null(fit)}.
#' @param return_all Logical. If \code{FALSE}, return only \eqn{F_{MSY}}.
#'   If \code{TRUE}, return a list with the full equilibrium BRP quantities.
#'
#' @return
#' If \code{return_all = FALSE}, a numeric estimate of \eqn{F_{MSY}}.
#'
#' If \code{return_all = TRUE}, a list with elements:
#' \describe{
#'   \item{\code{Fmsy}}{Estimated fishing mortality at MSY.}
#'   \item{\code{MSY}}{Maximum equilibrium yield.}
#'   \item{\code{Fmsy_approx}}{Discrete-grid approximation to \eqn{F_{MSY}}.}
#'   \item{\code{MSY_approx}}{Maximum equilibrium yield on the discrete grid.}
#'   \item{\code{Yeq}}{Equilibrium yield over \code{Fseq}.}
#'   \item{\code{SSB}}{Equilibrium spawning biomass over \code{Fseq}.}
#'   \item{\code{SSBrel}}{Relative equilibrium spawning biomass over \code{Fseq}.}
#'   \item{\code{Fseq}}{Fishing mortality sequence.}
#'   \item{\code{s}}{Steepness used in the Beverton-Holt relationship.}
#'   \item{\code{optimize}}{Result of the local \code{optimize()} call, or
#'   \code{NULL} if refinement was not used.}
#' }
#'
#' @details
#' The discrete approximation is obtained from:
#' \preformatted{
#' Yeq <- yieldf(prbrp, s = s)
#' Fmsy_approx <- prbrp$Fseq[which.max(Yeq)]
#' }
#'
#' If \code{fit} is provided and \code{refine = TRUE}, a local optimization is
#' carried out between the neighbouring \code{Fseq} values around the discrete
#' maximum. This provides a smoother estimate of \eqn{F_{MSY}} without having
#' to recompute the full BRP grid.
#'
#' This function keeps the standard separation between:
#' \itemize{
#'   \item per-recruit reference points (e.g. \code{Fspr40})
#'   \item stock-recruit based equilibrium reference points (e.g. \eqn{F_{MSY}})
#' }
#'
#' @seealso
#' \code{\link{prbrp_flicc}},
#' \code{\link{yieldf}},
#' \code{\link{ssbf}},
#' \code{\link{ssbrel}},
#' \code{\link{rr0_bh}}
#'
#' @examples
#' \dontrun{
#' # Build per-recruit BRP object
#' prb <- prbrp_flicc(fit)
#'
#' # Approximate Fmsy from the BRP grid
#' fmsy_flicc(prb, s = 0.7)
#'
#' # Refined estimate using the original fit
#' fmsy_flicc(prb, fit = fit, s = 0.7)
#'
#' # Return full equilibrium BRP output
#' msy <- fmsy_flicc(prb, fit = fit, s = 0.7, return_all = TRUE)
#' msy$Fmsy
#' msy$MSY
#'
#' plot(as.numeric(msy$Fseq), as.numeric(msy$Yeq), type = "l",
#'      xlab = "F", ylab = "Equilibrium yield")
#' abline(v = msy$Fmsy, col = 2, lty = 2)
#' }
#'
#' @export
fmsy_flicc <- function(prbrp, fit = NULL, s = 0.7,
                       nyears = 1,
                       scale_sel = TRUE,
                       spawn_time = 0,
                       refine = !is.null(fit),
                       return_all = FALSE) {

  Yeq <- yieldf(prbrp, s = s)
  Fseq <- as.numeric(prbrp$Fseq)

  i_max <- which.max(as.numeric(Yeq))
  Fmsy_approx <- Fseq[i_max]
  MSY_approx <- as.numeric(Yeq[,i_max])

  Fmsy <- Fmsy_approx
  MSY <- MSY_approx
  opt <- NULL

  if (isTRUE(refine) && !is.null(fit) && length(Fseq) >= 3) {

    i_lo <- max(1, i_max - 1)
    i_hi <- min(length(Fseq), i_max + 1)

    lower <- Fseq[i_lo]
    upper <- Fseq[i_hi]

    if (upper > lower) {
      obj <- function(f) {
        -yield_eq_flicc(
          fit = fit,
          F = f,
          nyears = nyears,
          s = s,
          scale_sel = scale_sel,
          spawn_time = spawn_time
        )
      }

      opt <- stats::optimize(obj, interval = c(lower, upper))
      Fmsy <- opt$minimum
      MSY <- -opt$objective
    }
  }

  if (!return_all) {
    return(Fmsy)
  }

  list(
    Fmsy = Fmsy,
    MSY = MSY,
    Fmsy_approx = Fmsy_approx,
    MSY_approx = MSY_approx,
    Yeq = Yeq,
    SSB = ssbf(prbrp, s = s),
    SSBrel = ssbrel(prbrp, s = s),
    Fseq = prbrp$Fseq,
    s = s,
    optimize = opt
  )
}


#' Equilibrium yield at a single fishing mortality
#'
#' Computes equilibrium yield for a single apical fishing mortality value by
#' combining yield-per-recruit from \code{pr_flicc()} with equilibrium
#' recruitment relative to unfished recruitment under a Beverton-Holt
#' stock-recruit relationship.
#'
#' This helper is intended primarily for single-point evaluation, for example
#' when refining \eqn{F_{MSY}} estimates by numerical optimization.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param F Numeric scalar apical annual fishing mortality.
#' @param nyears Integer. Number of terminal years used to construct biological
#'   reference inputs in \code{pr_flicc()}.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#' @param scale_sel Logical. Passed to \code{pr_flicc()} to control selectivity
#'   scaling.
#' @param spawn_time Numeric. Spawning time as fraction of the year.
#'
#' @return Numeric scalar equilibrium yield at fishing mortality \code{F}.
#'
#' @seealso \code{\link{pr_flicc}}, \code{\link{rr0_bh}},
#'   \code{\link{fmsy_flicc}}
#'
#' @examples
#' \dontrun{
#' yeql(fit, F = 0.2)
#' }
#'
#' @export
yeql <- function(fit, F, nyears = 1, s = 0.7,
                 scale_sel = TRUE, spawn_time = 0) {

  pr0 <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = 0,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )

  pr <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = F,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )

  SPR <- pr$SBPR / pr0$SBPR
  RR0 <- rr0_bh(SPR, s = s)

  as.numeric(pr$YPR * RR0)
}

#' Equilibrium spawning biomass at a single fishing mortality
#'
#' Computes equilibrium spawning biomass for a single apical fishing mortality
#' value by combining spawning-biomass-per-recruit from \code{pr_flicc()} with
#' equilibrium recruitment relative to unfished recruitment under a
#' Beverton-Holt stock-recruit relationship.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param F Numeric scalar apical annual fishing mortality.
#' @param nyears Integer. Number of terminal years used to construct biological
#'   reference inputs in \code{pr_flicc()}.
#' @param s Numeric Beverton-Holt steepness, default \code{0.7}.
#' @param scale_sel Logical. Passed to \code{pr_flicc()} to control selectivity
#'   scaling.
#' @param spawn_time Numeric. Spawning time as fraction of the year.
#'
#' @return Numeric scalar equilibrium spawning biomass at fishing mortality
#'   \code{F}.
#'
#' @seealso \code{\link{pr_flicc}}, \code{\link{rr0_bh}},
#'   \code{\link{ssbf}}
#'
#' @examples
#' \dontrun{
#' ssbeql(fit, F = 0.2)
#' }
#'
#' @export
ssbeql <- function(fit, F, nyears = 1, s = 0.7,
                   scale_sel = TRUE, spawn_time = 0) {

  pr0 <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = 0,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )

  pr <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = F,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )

  SPR <- pr$SBPR / pr0$SBPR
  RR0 <- rr0_bh(SPR, s = s)

  as.numeric(pr$SBPR * RR0)
}

#' Equilibrium spawning potential ratio at a single fishing mortality
#'
#' Computes spawning potential ratio (SPR) for a single apical fishing
#' mortality value as the ratio of spawning-biomass-per-recruit under
#' fishing to unfished spawning-biomass-per-recruit.
#'
#' This helper returns the classic per-recruit SPR and does not apply
#' stock-recruit scaling.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object.
#' @param F Numeric scalar apical annual fishing mortality.
#' @param nyears Integer. Number of terminal years used to construct biological
#'   reference inputs in \code{pr_flicc()}.
#' @param s Numeric Beverton-Holt steepness. Included for interface
#'   consistency; not used in the SPR calculation.
#' @param scale_sel Logical. Passed to \code{pr_flicc()} to control selectivity
#'   scaling.
#' @param spawn_time Numeric. Spawning time as fraction of the year.
#'
#' @return Numeric scalar spawning potential ratio at fishing mortality
#'   \code{F}.
#'
#' @seealso \code{\link{pr_flicc}}, \code{\link{spr_flicc}}
#'
#' @examples
#' \dontrun{
#' spreql(fit, F = 0.2)
#' }
#'
#' @export
spreql <- function(fit, F, nyears = 1, s = 0.7,
                   scale_sel = TRUE, spawn_time = 0) {

  pr0 <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = 0,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )

  pr <- pr_flicc(
    fit = fit,
    nyears = nyears,
    F = F,
    scale_sel = scale_sel,
    spawn_time = spawn_time
  )

  SPR <- pr$SBPR / pr0$SBPR

  as.numeric(SPR)
}

#' Current fishing mortality from an FLicc fit
#'
#' Extracts the terminal-year or recent-year average current fishing mortality
#' from a fitted \code{"flicc_tmb_fit"} object.
#'
#' Depending on \code{type}, the function returns either apical fishing
#' mortality (\code{"F"}) from \code{fit$report$Fap}, or fishing mortality
#' relative to natural mortality (\code{"FM"}) from \code{fit$report$FM}.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object with
#'   \code{fit$report$Fap} and \code{fit$report$FM}.
#' @param nyears Integer. Number of terminal years to average over.
#' @param type Character string, either \code{"F"} or \code{"FM"}.
#'
#' @return Numeric scalar current fishing mortality on the requested scale.
#'
#' @details
#' This helper assumes that \code{as_FLQuants()} has added \code{Fy} and
#' \code{FM} to the fitted report object.
#'
#' @seealso \code{\link{as_FLQuants}}
#'
#' @examples
#' \dontrun{
#' fcur_flicc(fit)
#' fcur_flicc(fit, type = "FM")
#' fcur_flicc(fit, nyears = 3, type = "F")
#' }
#'
#' @export
fcur_flicc <- function(fit, nyears = 1, type = c("F", "FM")) {

  type <- match.arg(type)

  year <- dimnames(fit$report$Fap)$year
  yrs <- tail(year, min(nyears, length(year)))

  if (type == "F") {
    an(yearMeans(fit$report$Fap[, yrs]))
  } else {
    an(yearMeans((fit$report$F/fit$report$M)[, yrs]))
  }
}

sprcur_flicc <- function(fit, nyears = 1) {


  year <- dimnames(fit$report$Fap)$year
  yrs <- tail(year, min(nyears, length(year)))


  an(yearMeans(fit$report$spr[, yrs]))

}
