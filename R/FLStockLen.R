

#' Construct an FLStockLen object from length bins and life-history parameters
#'
#' Builds an \code{FLStockLen} object from the length bins in a length-frequency
#' \code{FLQuants} together with scalar life-history inputs in an \code{FLPar}
#' object. The constructor derives midpoint lengths, weight-at-length,
#' maturity-at-length, and natural mortality-at-length under a chosen mortality
#' model.
#'
#' @param lfds An \code{FLQuants} containing length-frequency data by gear. The
#'   \code{len} dimension is used to define the length bins.
#' @param lhpar An \code{FLPar} object containing life-history parameters.
#'   Expected entries include \code{linf}, \code{L50}, \code{a}, and \code{b},
#'   and optionally \code{k}, \code{M}, \code{MK}, \code{L95}, and \code{s}.
#'   If \code{s} is not supplied, a default value of \code{0.7} is used.
#' @param m_model Character string specifying the natural mortality-at-length
#'   model. Supported options are \code{"constant"}, \code{"inverse"},
#'   \code{"Lorenzen"}, and \code{"Gislason"}.
#' @param reflen Optional reference length used by mortality models that require
#'   it, such as \code{"inverse"} and \code{"Lorenzen"}.
#'
#' @details
#' The \code{len} dimension of \code{lfd} is interpreted as lower length-bin
#' boundaries. Biological schedules are evaluated at bin midpoints derived from
#' these lower bounds.
#'
#' Weight at length is calculated as:
#' \deqn{W_L = aL^b}
#'
#' Maturity at length is represented by a logistic ogive parameterized by
#' \code{L50} and \code{L95}. If \code{L95} is not supplied, a default value is
#' constructed internally. If \code{s} is not supplied in \code{lhpar}, it
#' defaults to \code{0.7}.
#'
#' Natural mortality at length depends on \code{m_model}:
#' \itemize{
#'   \item \code{"constant"}: constant natural mortality at all lengths
#'   \item \code{"inverse"}: inverse-length mortality referenced to
#'     \code{reflen}
#'   \item \code{"Lorenzen"}: Lorenzen-style allometric mortality
#'   \item \code{"Gislason"}: Gislason length-based empirical mortality
#' }
#'
#' If \code{M} is not supplied directly but both \code{MK} and \code{k} are
#' available, then \code{M} is derived as \code{M = MK * k}. Likewise,
#' \code{MK} may be derived from \code{M / k} when needed.
#'
#' @return An \code{FLStockLen} object containing at least:
#'   \describe{
#'     \item{\code{stock.wt}}{Weight-at-length schedule}
#'     \item{\code{mat}}{Maturity-at-length schedule}
#'     \item{\code{m}}{Natural mortality-at-length schedule}
#'   }
#'
#' @examples
#' \dontrun{
#' data(alfonsino)
#' lhpar <- FLPar(
#'   linf = 55.7,
#'   k    = 0.08,
#'   M    = 0.162,
#'   L50  = 31.1859,
#'   a    = 0.004721956,
#'   b    = 3.146168
#' )
#'
#' stklen <- stocklen(
#'   lfd_alfonsino,
#'   lhpar = lhpar,
#'   m_model = "constant"
#' )
#'  stklen@lhpar
#' }
#'
#' @export

stocklen <- function(lfds, lhpar, m_model = c("constant", "inverse", "Lorenzen", "Gislason"),reflen=NULL) {
  m_model <- match.arg(m_model)

  pars <- dimnames(lhpar)$params


  if (!inherits(lfds, "FLQuants")) {
    stop("lfds must be class FLQuants")
  }
  if (!inherits(lhpar, "FLPar")) {
    stop("lhpar must be an FLPar")
  }


  getpar <- function(x, default = NA_real_) {
    if (x %in% pars) as.numeric(lhpar[x]) else default
  }

  linf   <- getpar("linf")
  k      <- getpar("k")
  M      <- getpar("M")
  Mk     <- getpar("Mk")
  L50    <- getpar("L50")
  L95    <- getpar("L95")
  a      <- getpar("a")
  b      <- getpar("b")
  s      <- getpar("s")
  if(is.null(reflen)) reflen <- getpar("reflen", L50)

  if (is.na(M) && !is.na(Mk) && !is.na(k)) M <- Mk * k
  if (is.na(Mk) && !is.na(M) && !is.na(k)) Mk <- M / k

  if(is.na(s)) s <- 0.7
  if (is.na(linf)) stop("lhpar must contain 'linf'")
  if (is.na(L50)) stop("lhpar must contain 'L50'")
  if (is.na(a) || is.na(b)) stop("lhpar must contain 'a' and 'b'")

  lens <- an(dimnames(lfds[[1]])$len)
  year <- an(dimnames(lfds[[1]])$year)
  its <- an(dimnames(lfds[[1]])$iter)
  if (anyNA(lens)) stop("lfds must have numeric 'len' dimnames")

  step <- c(diff(lens), tail(diff(lens), 1))
  midL <- lens + step / 2



  if (is.na(L95)) {
    L95 <- L50 + 0.2 * linf
  }
  if (L95 <= L50) stop("L95 must be > L50")

  # weight at length
  wt <- a * midL^b

  # maturity at length using L50/L95 parameterization
  dL <- L95 - L50
  mat <- 1 / (1 + exp(-log(19) * (midL - L50) / dL))

  # natural mortality at length
  if (is.na(M) && m_model %in% c("constant", "inverse")) {
    stop("Need M (or Mk and k) for selected m_model")
  }



  if (m_model %in% c("inverse", "Lorenzen") && is.null(reflen)) {
    stop("reflen must be supplied for m_model = '", m_model, "'")
  }
  if (m_model == "Gislason" && is.na(k)) {
    stop("k must be supplied in lhpar for m_model = 'Gislason'")
  }

  m <- switch(
    m_model,
    constant = rep(M, length(midL)),
    inverse  = M * reflen / midL,
    Lorenzen = m_lorenzen(midL, Mref = M, reflen = reflen, a = a, b = b),
    Gislason = m_gislason(midL, linf = linf, k = k)
  )

  # make FLQuants
  dn <- list(
    len = ac(lens),
    year = ac(year),
    unit = "unique",
    season = "all",
    area = "unique",
    iter = ac(its)
  )

  stock.n <- FLCore::FLQuant(rep(NA_real_, length(lens)), dimnames = dn)
  stock.wt <- FLCore::FLQuant(wt, dimnames = dn)
  mat.q <- FLCore::FLQuant(mat, dimnames = dn)
  m.q <- FLCore::FLQuant(m, dimnames = dn)
  catch.n <- FLCore::FLQuant(rep(NA_real_, length(lens)), dimnames = dn)
  catch.wt <- stock.wt
  flq0 <-   FLCore::FLQuant(rep(0, length(lens)), dimnames = dn)

  out <- FLCore::FLStockLen(
    stock.n = stock.n,
    catch.n = catch.n,
    catch.wt = catch.wt,
    stock.wt = stock.wt,
    mat = mat.q,
    m = m.q,
    harvest.spwn = flq0,
    m.spwn = flq0,
    discards.n = flq0,
    discards.wt = catch.wt,
    name = "FLicc stocklen"
  )

  lhpar.out <- FLPar(linf=linf,k=k,M=M,Mk=Mk,L50=L50,L95=L95,a=a,b=b,s=0.7)
  units(out) <- standardUnits(out)
  attr(out,"lhpar") <- lhpar.out
  return(out)
}

#' Spawning biomass by BRP year from an FLStockLen BRP object
#'
#' @param stk An object returned by \code{brp_stklen_flicc()}.
#' @param spawn_time Fraction of the interval elapsed before spawning.
#'
#' @return An \code{FLQuant} over the BRP year dimension.
#' @export
ssbl <- function(stk, spawn_time = attr(stk, "spawn_time", exact = TRUE) %||% 0) {
  z <- FLCore::harvest(stk) + FLCore::m(stk)
  FLCore::quantSums(FLCore::stock.n(stk) * FLCore::stock.wt(stk) *
                      FLCore::mat(stk) * exp(-z * spawn_time))
}

#' Fully selected fishing mortality from an FLStockLen BRP object
#'
#' @param stk An object returned by \code{brp_stklen_flicc()}.
#'
#' @return An \code{FLQuant} over the BRP year dimension.
#' @export
fapl <- function(stk) {
  apply(FLCore::harvest(stk) ,c(2:6), max, na.rm = TRUE)

}


#' Recruitment proxy from an FLStockLen object
#'
#' Returns the abundance in the first length class, used here as a recruitment
#' proxy for FLicc BRP objects.
#'
#' @param stk An object of class \code{FLStockLen}.
#'
#' @return An \code{FLQuant} over year.
#' @export
recl <- function(stk) {
  FLCore::stock.n(stk)[1]
}

#' Compute catch numbers from stock numbers and mortalities
#'
#' Computes catch-at-length or catch-at-age in numbers using the Baranov
#' catch equation.
#'
#' @param stock.n An \code{FLQuant} of stock numbers.
#' @param harvest An \code{FLQuant} of fishing mortality.
#' @param m An \code{FLQuant} of natural mortality.
#'
#' @return An \code{FLQuant} of catch numbers.
#' @export
compute_catch.n <- function(stock.n, harvest, m) {

  z <- harvest + m
  cn <- stock.n * harvest / z * (1 - exp(-z))
  cn[z == 0] <- 0

  return(cn)
}


#' Mean fishing mortality over a length range
#'
#' @param stk An \code{FLStockLen} object.
#' @param lrange Optional numeric vector of length 2 giving min and max length.
#'
#' @return An \code{FLQuant}.
#' @export
fbarl <- function(stk, lrange = NULL) {

  L <- as.numeric(dimnames(harvest(stk))$len)

  idx <- if (is.null(lrange)) {
    seq_along(L)
  } else {
    which(L >= lrange[1] & L <= lrange[2])
  }

  FLCore::quantMeans(harvest(stk)[idx,,,,,])
}

#' Update an FLStockLen object with FLicc fitted outputs
#'
#' Fills an existing \code{FLStockLen} object with equilibrium
#' length-based population and fishery quantities derived from a fitted
#' FLicc model.
#'
#' The function uses the fitted report components to reconstruct
#' numbers-at-length, fishing mortality, natural mortality, and catch,
#' and updates the corresponding slots of the supplied \code{FLStockLen}
#' template.
#'
#' @param fit A fitted \code{"flicc_tmb_fit"} object containing at least
#'   \code{report$N}, \code{report$FM}, and \code{report$Mk}.
#' @param year Optional year selection (currently not used explicitly but
#'   reserved for future extensions).
#' @param R0 Numeric recruitment scaling factor. Numbers-at-length are scaled
#'   such that the smallest length class equals \code{R0}.
#'
#' @return An \code{FLStockLen} object with updated length-based population and
#'   fishery quantities derived from the fitted FLicc outputs.
#'
#' @details
#' This function assumes that the fitted report contains at least the following
#' elements:
#' \describe{
#'   \item{\code{N}}{Numbers-at-length on a relative scale.}
#'   \item{\code{FM}}{Fishing mortality ratio on the \eqn{F/M} scale.}
#'   \item{\code{Mk}}{Natural mortality on the \eqn{M/k} scale.}
#' }
#'
#' It also assumes that the life-history parameter slot
#' \code{stklen@lhpar} contains \code{k}, which is used to convert
#' \code{Mk} and \code{FM} back to annual mortality rates:
#' \deqn{
#' M = Mk \cdot k \quad \text{and} \quad F = FM \cdot M
#' }
#'
#' Numbers-at-length are scaled relative to the smallest length class and
#' multiplied by \code{R0}. Catch-at-length is then computed using the
#' Baranov catch equation via \code{compute_catch.n()}.
#'
#' Discards are currently set to zero, and landings are assumed equal to catch.
#'
#' @seealso \code{\link{fiticc}}, \code{\link{compute_catch.n}}
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#'
#' stkl <- flicc_stklen(fit, stklen_alfonsino)
#'
#' stock.n(stkl)
#' harvest(stkl)
#' catch.n(stkl)
#'
#' @export
flicc_stklen <- function(fit, year = NULL, R0 = 1000){

  stklen <- fit$stklen
  # -----------------------------
  # checks
  # -----------------------------
  stopifnot(methods::is(stklen, "FLStockLen"))

  rep <- fit$report
  lhpar <- stklen@lhpar

  if (!all(c("N", "F", "M") %in% names(rep))) {
    stop("fit$report must contain 'N', 'FM', and 'Mk'.")
  }

  if (!"k" %in% rownames(lhpar)) {
    stop("stklen@lhpar must contain 'k'.")
  }

  # -----------------------------
  # clone template
  # -----------------------------
  stky <- stklen

  # -----------------------------
  # numbers-at-length (scaled)
  # -----------------------------
  Nrel <- rep$N %/% rep$N[1, ]   # relative to first length
  stock.n(stky) <- FLQuant(Nrel * R0/1000,unit="1000")

  # -----------------------------
  # mortality
  # -----------------------------
  k <- lhpar["k"]

  harvest(stky) <- rep$F
  m(stky)       <- rep$M           # M

  # -----------------------------
  # catch numbers
  # -----------------------------
  catch.n(stky) <- compute_catch.n(
    stock.n(stky),
    harvest(stky),
    m(stky)
  )

  landings.n(stky) <- catch.n(stky)
  discards.n(stky) <- catch.n(stky) * 0

  # -----------------------------
  # biomass quantities
  # -----------------------------
  stock(stky)    <- rep$spr
  catch(stky)    <- FLCore::computeCatch(stky)
  landings(stky) <- catch(stky)
  discards(stky) <- catch(stky) * 0

  # -----------------------------
  # units
  # -----------------------------
  units(stky) <- standardUnits(stky)

  return(stky)
}



#' Convert a FLicc fit and FLStockLen template to a simplified FLStockR
#'
#' Builds a simplified one-age \code{FLStockR} object from a fitted
#' \code{FLicc} model and an \code{FLStockLen} template.
#'
#' The function first updates the supplied \code{stklen} object with fitted
#' quantities using \code{flicc_stklen()}. It then constructs a simplified
#' annual \code{FLStockR} object in which:
#'
#' \itemize{
#'   \item \code{ssb()} represents spawning potential ratio (SPR),
#'   \item \code{fbar()} represents applied fishing mortality from
#'         \code{fapl(stk)},
#'   \item reference points are derived from \code{eqstklen()} using the
#'         target SPR level.
#' }
#'
#' If \code{rel = TRUE}, the returned object is rescaled to relative reference
#' point units:
#'
#' \itemize{
#'   \item \code{ssb()} is expressed as \eqn{SPR / SPR_{target}},
#'   \item \code{fbar()} is expressed as \eqn{F / F_{SPR}},
#'   \item \code{Fspr = 1} and \code{Bspr = 1}.
#' }
#'
#' This format is intended for plotting, diagnostics, and simple harvest
#' control rule applications using relative SPR-based metrics.
#'
#' @param fit A fitted \code{FLicc} object.
#' @param spr.btgt Numeric SPR target on the proportional scale. The default
#'   is \code{0.4}, corresponding to \code{Fspr40}.
#' @param bpa Numeric precautionary biomass threshold on the same SPR scale.
#'   The default is \code{0.2}.
#' @param blim Numeric biomass limit on the same SPR scale. The default is
#'   \code{0.1}.
#' @param rel Logical. If \code{FALSE} (default), the output is on the absolute
#'   SPR scale. If \code{TRUE}, SPR and fishing mortality are rescaled relative
#'   to the SPR target and \code{Fspr}, respectively.
#'
#' @return An object of class \code{FLStockR}.
#'
#' @details
#' The returned object is a simplified single-age representation. Quantities
#' such as catch, landings, discards, and harvest are carried forward as annual
#' fitted values, while the maturity slot is used to store SPR so that
#' \code{ssb()} returns the SPR metric directly.
#'
#' Reference points are computed from:
#' \preformatted{
#' eqstklen(fit, F = c(0, mean(stk@m)), spr.tgt = spr.btgt * 100)
#' }
#'
#' where \code{SPR0} and \code{Fspr} are extracted from the resulting reference
#' point table.
#'
#' @seealso \code{\link{flicc_stklen}}, \code{\link{eqstklen}},
#'   \code{\link[FLCore]{FLStock}}, \code{\link[FLasher]{FLStockR}}
#'
#' @examples
#' data(alfonsino)
#'
#' fit <- fiticc(
#'   lfd_alfonsino,
#'   stklen_alfonsino,
#'   sel_fun = c("dsnormal", "logistic"),
#'   catch_by_gear = c(0.7, 0.3)
#' )
#'
#' stkr <- flicc2FLStockR(fit)
#' ssb(stkr)
#' fbar(stkr)
#' refpts(stkr)
#'
#' stkr_rel <- flicc2FLStockR(fit, rel = TRUE)
#' ssb(stkr_rel)
#' fbar(stkr_rel)
#' refpts(stkr_rel)
#'
#' @export
flicc2FLStockR <- function(fit, spr.btgt = 0.4, bpa = 0.2, blim = 0.1, rel = FALSE) {

  stklen <- fit$stklen

  stk <- flicc_stklen(fit)

  year <- dimnames(stklen)$year

  N = as.FLQuant(data.frame(age=1,year=year,unit="unique",
                            season="all",area="unique",iter=1,data=1))
  B = N
  H = N
  C = N
  D = N
  L = N
  B[] <- FLQuant(stock(stk),quant="age")
  H[] <- FLQuant(fapl(stk),quant="age")
  C[] <- FLQuant(catch(stk),quant="age")
  L[] <-  FLQuant(landings(stk),quant="age")
  D[] <- FLQuant(discards(stk),quant="age")

  eqs<- eqstklen(fit,F=c(0,mean(stk@m)),spr.tgt = spr.btgt*100)
  spr0 <- an(eqs@refpts["SPR0"])
  fspr <-  an(eqs@refpts["Fspr"])

  if(rel){
  H <- H/ fspr
  B <- B/ spr.btgt
  }


  stkr = FLStockR(
    stock.n=N,
    catch.n = C,
    landings.n = C,
    discards.n = FLQuant(0, dimnames=list(age="1", year = (year))),
    stock.wt=FLQuant(1, dimnames=list(age="1", year = (year))),
    landings.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    discards.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    catch.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    mat=B,
    m=FLQuant(0.0001, dimnames=list(age="1", year = year)),
    harvest = H,
    m.spwn = FLQuant(0, dimnames=list(age="1", year = year)),
    harvest.spwn = FLQuant(0.0, dimnames=list(age="1", year = year))
  )


  units(stkr) = standardUnits(stkr)
  stkr@catch = computeCatch(stkr)
  stkr@landings = computeLandings(stkr)
  stkr@discards = computeDiscards(stkr)
  stkr@stock = B

  stkr@refpts = FLPar(
    Fspr =fspr ,
    Bspr = spr.btgt,
    Bpa= bpa,
    Blim= blim,
  )
  if(rel){
    stkr@refpts = FLPar(
      Fspr =1,
      Bspr = 1,
      Bpa= bpa/spr.btgt,
      Blim= blim/spr.btgt,
    )

  }

  return(stkr)
}
