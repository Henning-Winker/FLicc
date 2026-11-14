# Sex-structured FLicc: input helpers

#' Stock objects by sex for a sex-structured FLicc fit
#'
#' Builds one `FLStockLen` per sex with [stocklen()], for a sex-structured
#' fit with `fiticc(lfd, stklen = stocklen_sex(lfd, lhpar), ...)`. The female
#' stock is the reference: spawning potential (SPR) is female spawning output,
#' and F is scaled by the female k. Males have their own Linf, k, M, maturity
#' and weight at length; recruits are split by `settings$sex_ratio`
#' (proportion female, default 0.5), and fishing mortality and selectivity
#' are shared by the sexes.
#'
#' **Length data by sex.** Sex is carried in the FLR `unit` dimension of each
#' gear's length frequencies, with levels `"F"`, `"M"` and `"U"` (unsexed;
#' `"unique"` is read as unsexed too). As in Stock Synthesis, a gear-year can
#' hold
#' * an unsexed (combined) composition (`U`), compared with the predicted
#'   composition of females and males together (SS3 gender code 0);
#' * females only (`F`; code 1) or males only (`M`; code 2): the composition
#'   within that sex;
#' * both sexes (`F` and `M`; code 3): one joint composition over length and
#'   sex, so the sex ratio at length is informative;
#'
#' and combined and sexed samples together (partly sexed samples). Gears,
#' and years within a gear, can differ. With unsexed data only, the model
#' still splits the catch by sex from the life history and sex ratio; male
#' Linf is then not estimable and stays at its input value.
#'
#' @param lfd Length frequencies (`FLQuants` by gear), sexed or not; used
#'   for the length classes and years.
#' @param lhpar Life history by sex: an `FLPar` with a `unit` dimension
#'   containing `"F"` and `"M"` (e.g. from `ss3diags2::ss3_lhpar()`), or a
#'   list `list(F = FLPar, M = FLPar)`. Each needs the parameters used by
#'   [stocklen()] (`linf`, `k`, `M`, `L50`, `a`, `b`, ...).
#' @param ... Passed to [stocklen()] (e.g. `m_model`).
#' @return A list with elements `F` and `M` (`FLStockLen`).
#' @examples
#' \dontrun{
#' lh <- FLPar(array(c(129, 0.10, 0.32, 65, 1.8e-5, 2.95,
#'                     78,  0.20, 0.32, 65, 1.8e-5, 2.95), dim = c(6, 2, 1),
#'                   dimnames = list(params = c("linf", "k", "M", "L50", "a", "b"),
#'                                   unit = c("F", "M"), iter = 1)))
#' stk <- stocklen_sex(lfd, lh, m_model = "constant")
#' fit <- fiticc(lfd, stk, sel_fun = c("logistic", "dsnormal"), catch_by_gear = cbg)
#' }
#' @export
stocklen_sex <- function(lfd, lhpar, ...) {
  lh <- sex_lhpar_flicc(lhpar)
  # length classes and years only: sex units collapsed
  if (inherits(lfd, "FLQuants") || is.list(lfd)) {
    lfd1 <- FLCore::FLQuants(lapply(lfd, function(x)
      if (dim(x)[3] > 1) FLCore::unitSums(x) else x))
  } else {
    lfd1 <- if (dim(lfd)[3] > 1) FLCore::unitSums(lfd) else lfd
  }
  list(F = stocklen(lfd1, lh$F, ...), M = stocklen(lfd1, lh$M, ...))
}

# split life history by sex into list(F = FLPar, M = FLPar) (internal)
sex_lhpar_flicc <- function(lhpar) {
  if (is.list(lhpar) && !isS4(lhpar)) {
    if (!all(c("F", "M") %in% names(lhpar))) stop("lhpar list needs elements F and M.")
    return(list(F = lhpar$F, M = lhpar$M))
  }
  dn <- dimnames(lhpar)
  if (is.null(dn$unit) || !all(c("F", "M") %in% dn$unit))
    stop("lhpar needs a unit dimension with \"F\" and \"M\" (params x unit x iter), ",
         "or use a list(F = , M = ).")
  one <- function(u) {
    a <- lhpar@.Data[, u, , drop = FALSE]
    FLCore::FLPar(array(a, dim = c(dim(a)[1], dim(a)[3]),
                        dimnames = list(params = dn$params, iter = dn$iter)))
  }
  list(F = one("F"), M = one("M"))
}

# fill the sex-structure data of the template for objects made before it
# (one-sex model) (internal)
sex_defaults_tmb_flicc <- function(td) {
  if (!is.null(td$nsex)) return(td)
  nl <- td$nlen; ny <- td$nyear; ng <- td$ngear
  td$nsex <- 1L
  td$obs_sex <- array(0, dim = c(nl, ny, ng, 2))
  td$sex_type <- matrix(0L, ny, ng)
  td$sex_ratio <- 1; td$k_ratio <- 1; td$mk_ratio <- 1
  td$Mscaler_m <- td$Mscaler; td$wt_m <- td$wt; td$mat_m <- td$mat
  td
}
