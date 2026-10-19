#' SPR with maturity at the same relative size in each growth-type group
#'
#' Recomputes the spawning potential ratio of a GTG fit with the maturity
#' ogive scaled to each group's asymptotic length, as in LBSPR's GTG model
#' (Hordyk et al. 2016; `LBSPRsim_()`: "maturity at same relative size"):
#' \eqn{L_{50,g} = L_{50} L_{\infty,g} / L_\infty} and likewise for
#' \eqn{L_{95}}, so fast-growing groups mature at larger lengths and
#' slow-growing groups at smaller ones. FLicc's own SPR uses one
#' maturity-at-length ogive for all groups (`stocklen()`).
#'
#' The calculation is done after the fit from the fitted quantities (group
#' Linf and recruitment shares, M/K by group, F at length by year), with the
#' same survival recursion as the FLicc template; maturity does not enter the
#' likelihood, so selectivity, F and the fit are unchanged. Because maturity
#' is a fixed input, the result is applied as a correction factor
#' `adj = spr_gtg / spr` to FLicc's SPR and its confidence interval, with the
#' standard error of log SPR unchanged. With `mat = "length"` the function
#' reproduces FLicc's SPR (`adj = 1`), which serves as a check.
#'
#' @param fit A fitted FLicc object with `pop_model = "gtg"`.
#' @param mat `"relative"` (default; LBSPR) or `"length"` (FLicc's ogive).
#' @param sel_lag Integer; number of length classes by which F at length in the
#'   survival recursion lags the fitted selectivity. `0` (default) is FLicc's
#'   model. `1` mimics LBSPR, whose GTG model (`LBSPRsim_()`) evaluates
#'   selectivity in the survival recursion one class below the class midpoint
#'   used for the predicted catch; with `mat = "relative", sel_lag = 1` the
#'   result approximates LBSPR's reported SPR for the same SL50, SL95 and F/M.
#' @param level Confidence level of the interval.
#' @return A `data.frame` by year with `spr` (FLicc), `spr_gtg` (recomputed
#'   with `mat`), the factor `adj = spr_gtg / spr`, and `lwr`, `upr`: FLicc's
#'   SPR interval multiplied by `adj` (`NA` without standard errors).
#' @references Hordyk, A.R., Ono, K., Prince, J.D. and Walters, C.J. (2016). A
#'   simple length-structured model based on life history ratios and
#'   incorporating size-dependent selectivity: application to spawning
#'   potential ratios for data-poor stocks. Canadian Journal of Fisheries and
#'   Aquatic Sciences 73(12): 1787-1799.
#' @examples
#' \dontrun{
#' spr_gtg_flicc(fit)                  # LBSPR's maturity at relative size
#' spr_gtg_flicc(fit, mat = "length")  # reproduces fit$report$spr
#' spr_gtg_flicc(fit, sel_lag = 1)     # as LBSPR's reported SPR
#' }
#' @export
spr_gtg_flicc <- function(fit, mat = c("relative", "length"), sel_lag = 0, level = 0.95) {
  mat <- match.arg(mat)
  sel_lag <- as.integer(sel_lag)
  td <- fit$tmb_data
  if (!identical(as.integer(td$pop_model), 2L)) stop("spr_gtg_flicc() needs a fit with pop_model = \"gtg\".")
  rp <- fit$obj$report(fit$obj$env$last.par.best)   # raw template report
  LLB <- as.numeric(td$LLB); nl <- length(LLB)
  LB <- c(LLB, LLB[nl] + if (nl > 1) LLB[nl] - LLB[nl - 1] else 1)
  Lmid <- (LB[-1] + LB[-(nl + 1)]) / 2
  linf_g <- as.numeric(rp$gtgLinfs); recP <- as.numeric(rp$recP); ng <- length(linf_g)
  MK <- as.matrix(rp$MKMat); Ms <- as.matrix(td$Mscaler); Fl <- as.matrix(rp$Fk_l)
  wt <- as.matrix(td$wt); matL <- as.matrix(td$mat)
  ny <- ncol(Fl)
  # LBSPR's survival recursion uses selectivity sel_lag classes below the midpoint
  if (sel_lag > 0) Fl <- rbind(matrix(0, sel_lag, ny), Fl[seq_len(nl - sel_lag), , drop = FALSE])
  # maturity of each group
  lh <- fit$stklen@lhpar
  pn <- dimnames(lh)$params
  L50 <- as.numeric(lh["L50"])
  L95 <- if ("L95" %in% pn) as.numeric(lh["L95"]) else NA_real_
  if (!is.finite(L95)) L95 <- L50 + 0.2 * as.numeric(lh["linf"])   # as stocklen()
  Linf <- as.numeric(rp$Linf)
  mat_g <- function(g, y) {
    if (mat == "length") return(matL[, y])
    s <- linf_g[g] / Linf
    1 / (1 + exp(-log(19) * (Lmid - L50 * s) / ((L95 - L50) * s)))
  }
  # numbers in each length class of group g for total Z/K by class (template recursion)
  num_g <- function(g, zk) {
    s <- numeric(nl + 1); s[1] <- recP[g]
    for (l in seq_len(nl)) {
      a <- max(linf_g[g] - LB[l + 1], 1e-12); b <- max(linf_g[g] - LB[l], 1e-12)
      s[l + 1] <- s[l] * (a / b)^zk[l]
    }
    (s[seq_len(nl)] - s[-1]) / zk[seq_len(nl)]
  }
  yrs <- as.character(td$year_names)
  out <- vapply(seq_len(ny), function(y) {
    sf <- s0 <- 0
    for (g in seq_len(ng)) {
      z0 <- MK[seq_len(nl), g] * Ms[, 1]
      zf <- MK[seq_len(nl), g] * Ms[, y] + Fl[, y]
      m <- mat_g(g, y)
      s0 <- s0 + sum(num_g(g, z0) * wt[, y] * m)
      sf <- sf + sum(num_g(g, zf) * wt[, y] * m)
    }
    sf / s0
  }, numeric(1))
  spr <- as.numeric(rp$spr_y)
  adj <- out / spr
  ci <- tryCatch(ratio_ci_flicc(fit, "spr", level), error = function(e) NULL)
  lwr <- upr <- rep(NA_real_, ny)
  if (!is.null(ci) && nrow(ci) == ny) { lwr <- ci$lwr * adj; upr <- ci$upr * adj }
  data.frame(year = as.numeric(yrs), spr = spr, spr_gtg = out, adj = adj,
             lwr = lwr, upr = upr)
}
