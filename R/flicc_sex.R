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
#' **Life history by sex.** Linf, M/k and CVL are separate parameters for
#' each sex. Each is fixed at its input value unless a penalty (prior SD on
#' the log scale) is given: `linf.sd`, `Mk.sd`, `CVL.sd` for females and
#' `linf_m.sd`, `Mk_m.sd`, `CVL_m.sd` for males (male CVL input
#' `settings$CVL_m`, default `settings$CVL`). One sex can therefore be fixed
#' while the other is estimated, or both estimated. Sexed compositions,
#' especially joint F + M samples, carry the information (see
#' [plot_sexratio()]); with unsexed data only the male parameters are weakly
#' determined and should stay fixed or under informative penalties.
#'
#' **Surveys.** Gears named in `settings$surveys` are surveys: their
#' composition is the population seen through the survey's selectivity, with
#' no fishing mortality and no share of the catch (`catch_by_gear` is then
#' given for the fishery gears only). Surveys work in one-sex and
#' sex-structured fits; sexed survey samples are often the main source of
#' information on the sex ratio at length (see [plot_sexratio()]).
#'
#' **Per-recruit quantities.** Spawning biomass per recruit is female
#' spawning output and yield per recruit sums both sexes, both per total
#' recruit; Fspr and SPR therefore refer to females.
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
  if (is.null(td$gear_type)) td$gear_type <- rep(0L, td$ngear)
  if (!is.null(td$nsex)) return(td)
  nl <- td$nlen; ny <- td$nyear; ng <- td$ngear
  td$nsex <- 1L
  td$obs_sex <- array(0, dim = c(nl, ny, ng, 2))
  td$sex_type <- matrix(0L, ny, ng)
  td$sex_ratio <- 1; td$k_ratio <- 1; td$mk_ratio <- 1
  td$Mscaler_m <- td$Mscaler; td$wt_m <- td$wt; td$mat_m <- td$mat
  td
}

# gear types from settings$surveys (names or indices of survey gears):
# 0 fishery, 1 survey (internal)
survey_gears_flicc <- function(surveys, gear_names) {
  gt <- rep(0L, length(gear_names))
  if (is.null(surveys) || !length(surveys)) return(gt)
  idx <- if (is.numeric(surveys)) as.integer(surveys) else match(surveys, gear_names)
  if (any(is.na(idx)) || any(idx < 1 | idx > length(gear_names)))
    stop("settings$surveys must name gears of lfd (", paste(gear_names, collapse = ", "), ").")
  gt[idx] <- 1L
  gt
}

#' Observed and predicted sex ratio at length
#'
#' For a sex-structured fit, plots the proportion of females at length in the
#' sexed compositions (points, sized by the number of sexed fish; lines: 95%
#' binomial intervals) against the model prediction (line), by gear. Only
#' gear-years with both sexes sampled (SS3 gender code 3) carry information
#' on the sex ratio and are used. With `by_year = FALSE` (default) years are
#' pooled: observed counts are summed and the predicted proportions are
#' weighted by each year's sexed sample size.
#'
#' **Reading the plot.** The predicted sex ratio at length follows from the
#' growth (Linf, k) and natural mortality of each sex and the sex ratio at
#' recruitment. A rising share of females at large lengths is the usual
#' signature of females growing larger. Observed proportions systematically
#' above or below the line over a range of lengths point to male Linf or k,
#' M by sex, or `settings$sex_ratio`; a misfit only at small lengths often
#' reflects small fish that could not be sexed (sexing depends on length).
#'
#' @param fit A sex-structured FLicc fit.
#' @param gear Optional gears (names or indices); default: all with sexed data.
#' @param by_year Logical; one panel per year and gear.
#' @param years Optional years.
#' @param min_n Minimum number of sexed fish in a length bin to show the
#'   observed proportion.
#' @return A `ggplot` object; the data are in attribute `"data"`.
#' @export
plot_sexratio <- function(fit, gear = NULL, by_year = FALSE, years = NULL, min_n = 5) {
  td <- fit$tmb_data
  if (is.null(td$nsex) || td$nsex != 2L) stop("plot_sexratio() needs a sex-structured fit (stocklen_sex()).")
  rp <- fit$obj$report()
  ps <- rp$plen_sex
  gn <- td$gear_names; yn <- as.character(td$year_names)
  L <- td$Lmid
  sel_g <- seq_along(gn)
  if (!is.null(gear)) sel_g <- if (is.numeric(gear)) gear else match(gear, gn)
  sel_y <- seq_along(yn)
  if (!is.null(years)) sel_y <- which(yn %in% as.character(years))
  d <- do.call(rbind, lapply(sel_g, function(g) do.call(rbind, lapply(sel_y, function(y) {
    if (td$sex_type[y, g] != 3L) return(NULL)
    oF <- td$obs_sex[, y, g, 1]; oM <- td$obs_sex[, y, g, 2]
    data.frame(gear = gn[g], year = yn[y], len = L, oF = oF, oM = oM, n = sum(oF + oM),
               pF = ps[, y, g, 1], pM = ps[, y, g, 2])
  }))))
  if (is.null(d) || !nrow(d)) stop("No gear-years with both sexes sampled.")
  grp <- if (by_year) c("gear", "year", "len") else c("gear", "len")
  # predicted: sum over years of n_y * proportion within the year's sexed composition
  d$wF <- d$n * d$pF / ave(d$pF + d$pM, d$gear, d$year, FUN = sum)
  d$wM <- d$n * d$pM / ave(d$pF + d$pM, d$gear, d$year, FUN = sum)
  a <- stats::aggregate(d[, c("oF", "oM", "wF", "wM")], d[grp], sum)
  a$nobs <- a$oF + a$oM
  a$obs <- ifelse(a$nobs >= min_n, a$oF / pmax(a$nobs, 1), NA)
  a$pred <- a$wF / (a$wF + a$wM)
  ci <- t(mapply(function(x, n) if (n >= min_n) stats::binom.test(round(x), round(n))$conf.int else c(NA, NA),
                 a$oF, a$nobs))
  a$lo <- ci[, 1]; a$hi <- ci[, 2]
  a$gear <- factor(a$gear, levels = gn)
  p <- ggplot2::ggplot(a, ggplot2::aes(len)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), colour = "grey60", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = obs, size = nobs), colour = "#C2185B", alpha = 0.7, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = pred), colour = "grey15", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_size_area(max_size = 3, name = "Sexed fish") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Length", y = "Proportion female",
                  title = "Sex ratio at length: observed (points, 95% CI) and predicted (line)",
                  subtitle = if (by_year) NULL else "Years pooled") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  p <- p + if (by_year) ggplot2::facet_grid(year ~ gear) else ggplot2::facet_wrap(~gear)
  attr(p, "data") <- a
  p
}

# life-history plots: a sex-structured stklen (list F, M) or fit becomes
# FLStocks(Female, Male); anything else is returned unchanged (internal)
sex_stocks_flicc <- function(object) {
  if (is.list(object) && !isS4(object) && !is.null(object$tmb_data)) {
    if (!is.null(object$stklen_m)) return(FLCore::FLStocks(Female = object$stklen, Male = object$stklen_m))
    return(object$stklen)
  }
  if (is.list(object) && !isS4(object) && all(c("F", "M") %in% names(object)))
    return(FLCore::FLStocks(Female = object$F, Male = object$M))
  object
}

# length frequencies of a fit with the sex units used in the likelihood:
# unsexed counts as "U" (or "unique" for one-sex fits) and sexed counts as "F",
# "M"; units without data are dropped (internal)
lfd_from_fit_flicc <- function(fit) {
  r <- fit$report
  if (is.null(r$obslen_sex)) return(r$obslen)
  out <- lapply(names(r$obslen), function(g) {
    u <- r$obslen[[g]]; sx <- r$obslen_sex[[g]]
    keep <- c(U = sum(u, na.rm = TRUE) > 0, F = sum(sx[, , "F"], na.rm = TRUE) > 0,
              M = sum(sx[, , "M"], na.rm = TRUE) > 0)
    if (!any(keep)) keep["U"] <- TRUE
    d <- dimnames(u); d$unit <- names(keep)[keep]
    q <- FLCore::FLQuant(0, dimnames = d, units = units(u))
    if (keep["U"]) q[, , "U"] <- u
    for (s in c("F", "M")) if (keep[s]) q[, , s] <- sx[, , s]
    q
  })
  FLCore::FLQuants(stats::setNames(out, names(r$obslen)))
}
