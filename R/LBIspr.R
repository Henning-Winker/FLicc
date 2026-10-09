#' Gear-specific length-based indicator relative to Fspr reference
#'
#' Computes a gear-specific length-based indicator by comparing the observed
#' proportion above a reference threshold length to the expected proportion
#' above that threshold under equilibrium numbers-at-length at Fsprx,
#' filtered through gear-specific selectivity.
#'
#' The threshold is defined separately for each gear from:
#'   Nf(Fsprx) * sel_gear
#'
#' @param fit fitted flicc_tmb_fit object
#' @param gear optional character vector of gears to include; defaults to all
#' @param spr target SPR percentage, default 40
#' @param thresh cumulative threshold used to define Lref, default 0.75
#' @param nyears number of terminal years to average
#' @param scale_sel logical; passed to nf_flicc
#'
#' @return *FLIndices* with *FLIndexBiomass*
#' @export
LBIspr<- function(fit, gear = NULL, spr = 40, thresh = 0.75,
                         nyears = 1, scale_sel = TRUE) {


  gear_names <- names(fit$report$sel_gear)
  if (is.null(gear)) {
    gear <- gear_names
  }
  gear <- as.character(gear)
  yrs <- dimnames(fit$report$obslen[[1]])$year
  lens <- dimnames(fit$report$obslen[[1]])$len

  sel.pattern <- FLQuant(
    dimnames = list(
      len = lens,
      year = yrs,
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )

  flqs = FLQuants(lapply(fit$report$Fk,function(x){
      x <-  FLQuant(x,quant='age')
      x[] <- NA
      units(x) = "lbi"
      x
  }))

  if (!all(gear %in% gear_names)) {
    stop("Unknown gear name(s): ", paste(setdiff(gear, gear_names), collapse = ", "))
  }

  # reference F on F scale
  Ftgt <- fspr_flicc(fit, spr = spr, nyears = nyears, input = "F")

  # equilibrium numbers-at-length at Fsprx
  Nref <- nf_flicc(fit, nyears = nyears, F = Ftgt, scale_sel = scale_sel)
  Len  <- as.numeric(dimnames(Nref)$len)



  vals <- numeric(length(gear))
  Lref_out <- numeric(length(gear))
  pref_out <- numeric(length(gear))
  pobs_out <- numeric(length(gear))

  LFDobs <- fit$report$obslen
  LFDref <- fit$report$obslen
  Lref <- NULL

  for (i in seq_along(gear)) {
    g <- gear[i]
    # gear-specific selectivity
    sg <- fit$report$sel_gear[[g]]


    # reference vulnerable numbers for this gear
    vref <- Nref * sg

    if (sum(vref, na.rm = TRUE) <= 0) {
      vals[i] <- NA_real_
      Lref_out[i] <- NA_real_
      pref_out[i] <- NA_real_
      pobs_out[i] <- NA_real_
      next
    }

    # define threshold length from cumulative vulnerable numbers
    vref2 <- vref[-1]
    Len2  <- Len[-1]

    cums = apply(vref2, 2:6, cumsum)
    n_thresh <- sum(vref2, na.rm = TRUE) * thresh
    Lthresh<- Len2[which.min((n_thresh - cums)^2)]
    Li <- ac(Len2[Len2>= Lthresh])

    LFDref[[g]][] <- vref
    Lref <- c(Lref,Lthresh)


    # expected proportion above threshold
    pref <- sum(vref[Li,], na.rm = TRUE) / sum(vref, na.rm = TRUE)

    for(y in seq(yrs)){
    # observed gear-specific length composition
    obs_g <- fit$report$obslen[[g]][,y]

    if (sum(obs_g, na.rm = TRUE) <= 0) {
      next
    }

    pobs <- sum(obs_g[Li,], na.rm = TRUE) / sum(obs_g, na.rm = TRUE)

    flqs[[g]][,y] <- pobs / pref
  }
  }

  out <- FLIndices(Map(function(x,y){
    idx <-FLIndexBiomass(index=x)
    idx@range[c("startf","endf")] =c(0.4,0.6)
    idx
    },x=flqs,y= fit$report$sel_gear))

 names(Lref) <- gear
 attr(out,"LFDref") <- LFDref
 attr(out,"Lref") <- Lref

 return(out)

}


#' Gear-specific mean-length indicator relative to a model reference
#'
#' Model-based analogue of the ICES length indicator
#' \eqn{f = \bar L / L_{F=M}}. For each gear, the reference mean length is
#' computed from the fitted model's equilibrium numbers-at-length at a
#' reference fishing mortality, filtered through the gear's fitted
#' selectivity:
#' \deqn{\bar L^{ref}_g = \frac{\sum_l \bar l\, N_l(F_{ref})\, s_g(l)}
#'                            {\sum_l N_l(F_{ref})\, s_g(l)},}
#' and the indicator is \eqn{\bar L^{obs}_{g,y} / \bar L^{ref}_g}, so that
#' 1 means the observed catch is as large as expected at the reference.
#'
#' Unlike the ICES approximation
#' \eqn{L_{F=M} = (L_\infty + 2\cdot1.5\,L_c)/(1 + 2\cdot1.5)}, which assumes
#' knife-edge (logistic) selectivity and a single growth curve, the
#' reference here uses the fitted growth variability (GTG), natural
#' mortality and gear selectivity, including dome-shaped curves.
#'
#' @param fit Fitted FLicc object.
#' @param gear Optional character vector of gears; defaults to all.
#' @param ref `"FM"` (default): reference at \eqn{F = FM \times M}
#'   (ICES-type \eqn{L_{F=M}} with `FM = 1`); `"spr"`: reference at
#'   \eqn{F_{\mathrm{SPR}x}}.
#' @param FM F/M ratio of the reference when `ref = "FM"`. Default 1.
#' @param spr Target SPR in percent when `ref = "spr"`. Default 40.
#' @param lc Length cut-off for both observed and reference means:
#'   `"sel50"` (default) uses the first length where the gear's selectivity
#'   reaches 0.5, which excludes poorly selected small fish as in the ICES
#'   \eqn{L_c} convention; `"none"` uses all lengths (except the first bin);
#'   or a numeric value (single, or named by gear).
#' @param nyears,scale_sel Passed to [nf_flicc()] / [fspr_flicc()].
#'
#' @return *FLIndices* of *FLIndexBiomass*, one per gear, with
#'   \eqn{\bar L^{obs}/\bar L^{ref}} by year. Attributes `"Lmean_ref"`
#'   (reference mean length) and `"Lc"` (cut-off) are named by gear.
#'
#' @details
#' Like [LBIspr()], the level of the indicator depends on the model's scale
#' (M, Linf, growth variability, selectivity) through the reference; its
#' trend over time does not, because a constant reference cancels in a
#' ratio of years.
#'
#' @seealso [LBIspr()], [nf_flicc()], [fspr_flicc()]
#'
#' @examples
#' \dontrun{
#' lbm <- LBImean(fit)                 # Lmean / Lmean(F = M), per gear
#' lbm_spr <- LBImean(fit, ref = "spr", spr = 40)
#' attr(lbm, "Lmean_ref")
#' }
#' @export
LBImean <- function(fit, gear = NULL, ref = c("FM", "spr"), FM = 1, spr = 40,
                    lc = "sel50", nyears = 1, scale_sel = TRUE) {

  ref <- match.arg(ref)

  gear_names <- names(fit$report$sel_gear)
  if (is.null(gear)) gear <- gear_names
  gear <- as.character(gear)
  if (!all(gear %in% gear_names))
    stop("Unknown gear name(s): ", paste(setdiff(gear, gear_names), collapse = ", "))

  yrs <- dimnames(fit$report$obslen[[1]])$year

  ## equilibrium numbers-at-length at the reference F
  Nref <- if (ref == "FM") {
    nf_flicc(fit, nyears = nyears, FM = FM, scale_sel = scale_sel)
  } else {
    Ftgt <- fspr_flicc(fit, spr = spr, nyears = nyears, input = "F")
    nf_flicc(fit, nyears = nyears, F = Ftgt, scale_sel = scale_sel)
  }

  Len <- as.numeric(dimnames(Nref)$len)
  bin <- if (length(Len) > 1) stats::median(diff(Len)) else 1
  mid <- Len + bin / 2

  flqs <- FLQuants(lapply(setNames(nm = gear), function(g) {
    FLQuant(NA_real_, dimnames = list(age = "all", year = yrs), units = "lbi")
  }))
  Lmean_ref <- setNames(rep(NA_real_, length(gear)), gear)
  Lc_out    <- setNames(rep(NA_real_, length(gear)), gear)

  for (g in gear) {
    sg <- fit$report$sel_gear[[g]]
    v  <- Nref * sg
    if (dim(v)[2] > 1)  v  <- yearMeans(v)       # one reference length profile
    if (dim(sg)[2] > 1) sg <- yearMeans(sg)
    vref <- c(v)
    sgv  <- c(sg)

    ## cut-off
    Lc <- if (is.numeric(lc)) {
      if (!is.null(names(lc))) lc[[g]] else lc[1]
    } else if (lc == "sel50") {
      Len[which(sgv / max(sgv, na.rm = TRUE) >= 0.5)[1]]
    } else {
      Len[2]                              # drop the first (pooled) bin only
    }
    keep <- Len >= Lc
    Lc_out[g] <- Lc

    if (sum(vref[keep], na.rm = TRUE) <= 0) next
    Lmean_ref[g] <- sum(mid[keep] * vref[keep], na.rm = TRUE) /
      sum(vref[keep], na.rm = TRUE)

    for (y in seq_along(yrs)) {
      obs <- c(fit$report$obslen[[g]][, y])
      if (sum(obs[keep], na.rm = TRUE) <= 0) next
      Lmean_obs <- sum(mid[keep] * obs[keep], na.rm = TRUE) /
        sum(obs[keep], na.rm = TRUE)
      flqs[[g]][, y] <- Lmean_obs / Lmean_ref[g]
    }
  }

  out <- FLIndices(lapply(flqs, function(x) {
    idx <- FLIndexBiomass(index = x)
    idx@range[c("startf", "endf")] <- c(0.4, 0.6)
    idx
  }))
  attr(out, "Lmean_ref") <- Lmean_ref
  attr(out, "Lc") <- Lc_out
  out
}

