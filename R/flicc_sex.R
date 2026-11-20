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
#'   containing `"F"` and `"M"` (e.g. from [lhpar_sex()] or
#'   `ss3diags2::ss3_lhpar()`), or a list `list(F = FLPar, M = FLPar)`. Each
#'   needs the parameters used by [stocklen()] (`linf`, `k`, `M`, `L50`, `a`,
#'   `b`, ...). A one-sex `FLPar` is taken as the female life history, with
#'   the male values given in `male`.
#' @param male Optional male values for a one-sex `lhpar`: a named list or
#'   vector (e.g. `list(linf = 45, k = 0.25)`) or an `FLPar`; parameters not
#'   given are copied from the females (see [lhpar_sex()]).
#' @param ... Passed to [stocklen()] (e.g. `m_model`).
#' @return A list with elements `F` and `M` (`FLStockLen`).
#' @examples
#' \dontrun{
#' # one-sex life history plus the male values that differ
#' stk <- stocklen_sex(lfd, lhpar, male = list(linf = 78, k = 0.20), m_model = "constant")
#'
#' # or a full FLPar by sex
#' lh <- FLPar(array(c(129, 0.10, 0.32, 65, 1.8e-5, 2.95,
#'                     78,  0.20, 0.32, 65, 1.8e-5, 2.95), dim = c(6, 2, 1),
#'                   dimnames = list(params = c("linf", "k", "M", "L50", "a", "b"),
#'                                   unit = c("F", "M"), iter = 1)))
#' stk <- stocklen_sex(lfd, lh, m_model = "constant")
#' fit <- fiticc(lfd, stk, sel_fun = c("logistic", "dsnormal"), catch_by_gear = cbg)
#' }
#' @export
stocklen_sex <- function(lfd, lhpar, male = NULL, ...) {
  if (!is.null(male)) lhpar <- lhpar_sex(lhpar, male)
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

#' Life-history parameters by sex
#'
#' Builds the life history by sex for [stocklen_sex()] from a one-sex
#' `FLPar` (taken as the females) and only the male values that differ:
#' all other male parameters are copied from the females. The result is an
#' `FLPar` with a `unit` dimension (`"F"`, `"M"`), the same format as
#' `ss3diags2::ss3_lhpar()`.
#'
#' @param lhpar One-sex `FLPar` of life-history parameters (`linf`, `k`,
#'   `M`, `L50`, `a`, `b`, ...), used for females; iterations are kept.
#' @param male Male values that differ: a named list or vector, or an
#'   `FLPar` (its parameters replace the female values; one value or one per
#'   iteration). Male-only parameters not in `lhpar` are added (set to `NA`
#'   for females).
#' @param ... Male values as named arguments, alternative to `male`
#'   (e.g. `linf = 45, k = 0.25`).
#' @return An `FLPar` with dimensions params x unit (F, M) x iter.
#' @seealso [stocklen_sex()]
#' @examples
#' data(alfonsino)
#' lh <- lhpar_sex(lhpar_alfonsino, linf = 48, k = 0.10)
#' lh
#' # same as
#' lhpar_sex(lhpar_alfonsino, male = list(linf = 48, k = 0.10))
#' @export
lhpar_sex <- function(lhpar, male = NULL, ...) {
  dots <- list(...)
  if (!is.null(male) && length(dots)) stop("Give the male values in 'male' or as named arguments, not both.")
  if (is.null(male)) male <- dots
  if (!is.null(dimnames(lhpar)$unit))
    stop("lhpar already has a unit dimension; give a one-sex FLPar (used for females).")
  dn <- dimnames(lhpar)
  pf <- dn$params; it <- dn$iter; ni <- length(it)
  if (inherits(male, "FLPar")) {
    mv <- lapply(dimnames(male)$params, function(p) as.numeric(male[p]))
    names(mv) <- dimnames(male)$params
  } else {
    mv <- as.list(male)
  }
  if (length(mv) && (is.null(names(mv)) || any(names(mv) == "")))
    stop("Male values must be named (e.g. linf = 45).")
  pars <- union(pf, names(mv))
  a <- array(NA_real_, dim = c(length(pars), 2, ni),
             dimnames = list(params = pars, unit = c("F", "M"), iter = it))
  a[pf, "F", ] <- lhpar@.Data
  a[pf, "M", ] <- lhpar@.Data
  for (p in names(mv)) {
    v <- as.numeric(mv[[p]])
    if (!(length(v) %in% c(1, ni))) stop("Male value '", p, "' needs 1 or ", ni, " values.")
    a[p, "M", ] <- rep(v, length.out = ni)
  }
  FLCore::FLPar(a, units = "NA")
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

# number of leading length bins with sexes combined in sexed samples, by gear,
# from settings$sex_minlen (named by gear, or one value per gear) (internal)
sex_comb_flicc <- function(sex_minlen, gear_names, lens) {
  cb <- rep(0L, length(gear_names))
  if (is.null(sex_minlen) || !length(sex_minlen)) return(cb)
  v <- as.numeric(sex_minlen)
  if (!is.null(names(sex_minlen))) {
    bad <- setdiff(names(sex_minlen), gear_names)
    if (length(bad)) stop("settings$sex_minlen names must be gears of lfd; unknown: ", paste(bad, collapse = ", "))
    idx <- match(names(sex_minlen), gear_names)
  } else {
    if (length(v) == 1) v <- rep(v, length(gear_names))
    if (length(v) != length(gear_names)) stop("settings$sex_minlen: one value per gear, or named by gear.")
    idx <- seq_along(gear_names)
  }
  for (j in seq_along(idx)) if (is.finite(v[j]) && v[j] > 0) cb[idx[j]] <- as.integer(sum(lens < v[j]))
  cb
}

# fill the sex-structure data of the template for objects made before it
# (one-sex model) (internal)
sex_defaults_tmb_flicc <- function(td) {
  if (is.null(td$gear_type)) td$gear_type <- rep(0L, td$ngear)
  if (is.null(td$sex_comb)) td$sex_comb <- rep(0L, td$ngear)
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
    cb <- if (!is.null(td$sex_comb)) td$sex_comb[g] else 0L
    if (td$sex_type[y, g] != 3L && !(cb > 0 && td$sex_type[y, g] > 0)) return(NULL)
    oF <- td$obs_sex[, y, g, 1]; oM <- td$obs_sex[, y, g, 2]
    if (cb > 0) { oF[seq_len(cb)] <- 0; oM[seq_len(cb)] <- 0 }    # sexes combined below sex_minlen: no ratio
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

# observed length compositions of all fish measured per gear: unsexed plus
# sexed samples (sexes summed); the one-sex observed compositions otherwise
# (internal; used by the model-based indicators)
obslen_all_flicc <- function(fit) {
  r <- fit$report
  if (is.null(r$obslen_sex)) return(r$obslen)
  FLCore::FLQuants(stats::setNames(lapply(names(r$obslen), function(g) {
    z <- FLCore::unitSums(r$obslen_sex[[g]]); dimnames(z)$unit <- "unique"
    r$obslen[[g]] + z
  }), names(r$obslen)))
}


# ---------------------------------------------------------------------------
# Mirrored length compositions by sex (r4ss style): females above zero (red),
# males below (blue), unsexed in their own panel column (grey)
# ---------------------------------------------------------------------------

# colours used for the sexes in all sex views
sex_cols_flicc <- function(cols = NULL) {
  def <- c(F = "#D7191C", M = "#2C7BB6", U = "grey45", C = "#7B3294")
  if (!is.null(cols)) def[names(cols)] <- cols
  def
}

# long table of compositions by gear, year and sex from an FLQuants with units
# (F, M, U or unique); panel = gear for unsexed, "gear sexed" for F/M
# scale: "raw" counts, "prop" proportions within panel-year (females and males
# jointly where both were sampled, so the sex ratio is kept), "relmax" divided
# by the panel-year maximum
lfd_sex_long_flicc <- function(lfd, scale = c("raw", "prop", "relmax"), iter = 1, sex_minlen = NULL) {
  scale <- match.arg(scale)
  out <- list()
  cbl <- if (is.null(sex_minlen)) NULL else
    stats::setNames(if (is.null(names(sex_minlen))) rep(as.numeric(sex_minlen), length.out = length(lfd)) else
      as.numeric(sex_minlen)[match(names(lfd), names(sex_minlen))], names(lfd))
  for (g in names(lfd)) {
    q <- lfd[[g]]
    if (dim(q)[6] > 1) q <- q[, , , , , iter]
    sx <- sex_code_flicc(dimnames(q)$unit)
    df <- as.data.frame(q)
    df$len <- as.numeric(as.character(df$len))
    df$year <- as.character(df$year)
    df$sex <- sx[match(as.character(df$unit), dimnames(q)$unit)]
    if (!is.null(cbl) && is.finite(cbl[g]))      # sexes combined below sex_minlen
      df$sex[df$sex %in% c("F", "M") & df$len < cbl[g]] <- "C"
    df$data[!is.finite(df$data)] <- 0
    df <- stats::aggregate(data ~ len + year + sex, data = df, FUN = sum)
    df$gear <- g
    df$panel <- ifelse(df$sex == "U", g, paste(g, "sexed"))
    out[[g]] <- df
  }
  df <- do.call(rbind, out)
  tot <- stats::ave(df$data, df$panel, df$year, FUN = sum)
  df <- df[tot > 0, , drop = FALSE]
  if (scale == "prop") {
    df$data <- df$data / stats::ave(df$data, df$panel, df$year, FUN = sum)
  } else if (scale == "relmax") {
    df$data <- df$data / stats::ave(df$data, df$panel, df$year, FUN = max)
  }
  df
}

# observed and predicted compositions of a fit by gear, year and sex, as
# compared in the likelihood: unsexed within gear; females and males jointly
# (sum to 1 over both) where both were sampled, within sex where one was
lenfit_sex_long_flicc <- function(fit) {
  rp <- fit$report
  gears <- names(rp$obslen)
  rows <- list()
  todf <- function(q, nm) {
    d <- as.data.frame(q)
    data.frame(len = as.numeric(as.character(d$len)), year = as.character(d$year),
               v = d$data, stringsAsFactors = FALSE)
  }
  for (g in gears) {
    o <- todf(rp$obslen[[g]]); p <- todf(rp$predlen_gear[[g]])
    s <- stats::ave(o$v, o$year, FUN = function(x) sum(x, na.rm = TRUE))
    keep <- s > 0
    if (any(keep))
      rows[[length(rows) + 1]] <- data.frame(gear = g, panel = g, sex = "U",
        len = o$len[keep], year = o$year[keep], observed = o$v[keep] / s[keep],
        fitted = p$v[keep], stringsAsFactors = FALSE)
    if (is.null(rp$obslen_sex)) next
    oF <- todf(rp$obslen_sex[[g]][, , "F"]); oM <- todf(rp$obslen_sex[[g]][, , "M"])
    pF <- todf(rp$predlen_sex[[g]][, , "F"]); pM <- todf(rp$predlen_sex[[g]][, , "M"])
    gi <- match(g, gears)
    cb <- if (!is.null(fit$tmb_data$sex_comb)) fit$tmb_data$sex_comb[gi] else 0L
    for (y in unique(oF$year)) {
      i <- oF$year == y
      sF <- sum(oF$v[i], na.rm = TRUE); sM <- sum(oM$v[i], na.rm = TRUE)
      if (!(sF + sM > 0)) next
      if (cb > 0) {                   # sexes combined in the first cb bins, joint above
        lv <- sort(unique(oF$len[i])); lc <- lv[seq_len(min(cb, length(lv)))]
        iF <- which(i); ic <- iF[oF$len[iF] %in% lc]; iu <- iF[!(oF$len[iF] %in% lc)]
        tot <- sF + sM; ptot <- sum(pF$v[i] + pM$v[i])
        # combined (F + M) expected over all lengths, observed below the cut-off
        rows[[length(rows) + 1]] <- data.frame(gear = g, panel = paste(g, "sexed"), sex = "C",
          len = oF$len[c(ic, iu)], year = y, observed = c((oF$v[ic] + oM$v[ic]) / tot, rep(NA, length(iu))),
          fitted = (pF$v[c(ic, iu)] + pM$v[c(ic, iu)]) / ptot, stringsAsFactors = FALSE)
        # female and male expected lines over all lengths (below the cut-off
        # the model's split of the combined cell), so they run continuously;
        # observed sexes only above the cut-off
        rows[[length(rows) + 1]] <- data.frame(gear = g, panel = paste(g, "sexed"), sex = "F",
          len = oF$len[c(ic, iu)], year = y, observed = c(rep(NA, length(ic)), oF$v[iu] / tot),
          fitted = pF$v[c(ic, iu)] / ptot, stringsAsFactors = FALSE)
        rows[[length(rows) + 1]] <- data.frame(gear = g, panel = paste(g, "sexed"), sex = "M",
          len = oM$len[c(ic, iu)], year = y, observed = c(rep(NA, length(ic)), oM$v[iu] / tot),
          fitted = pM$v[c(ic, iu)] / ptot, stringsAsFactors = FALSE)
        next
      }
      if (sF > 0 && sM > 0) {         # joint length x sex composition
        dF <- oF$v[i] / (sF + sM); dM <- oM$v[i] / (sF + sM)
        fF <- pF$v[i]; fM <- pM$v[i]
      } else {                        # one sex only: within that sex
        dF <- if (sF > 0) oF$v[i] / sF else NA; dM <- if (sM > 0) oM$v[i] / sM else NA
        fF <- if (sF > 0) pF$v[i] / sum(pF$v[i]) else NA
        fM <- if (sM > 0) pM$v[i] / sum(pM$v[i]) else NA
      }
      if (sF > 0) rows[[length(rows) + 1]] <- data.frame(gear = g, panel = paste(g, "sexed"),
        sex = "F", len = oF$len[i], year = y, observed = dF, fitted = fF, stringsAsFactors = FALSE)
      if (sM > 0) rows[[length(rows) + 1]] <- data.frame(gear = g, panel = paste(g, "sexed"),
        sex = "M", len = oM$len[i], year = y, observed = dM, fitted = fM, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  df$fitted[!is.finite(df$fitted)] <- NA
  df
}

# widths of length bins with lower bounds L (last bin: the previous width) (internal)
bin_width_flicc <- function(L) {
  u <- sort(unique(L)); w <- if (length(u) > 1) c(diff(u), utils::tail(diff(u), 1)) else 1
  w[match(L, u)]
}

# the mirrored plot: df with panel, year, sex, len, observed and optionally fitted
plot_sex_mirror_flicc <- function(df, cols, len_by = 10, ylab = "Proportion",
                                  title = NULL, ncol_years = NULL, scales = "free_y",
                                  fitted_U = "black", linewidth = 0.8) {
  # bins of unequal width: proportions per unit length, so that wider bins
  # do not show as steps
  df$w <- bin_width_flicc(df$len)
  if (length(unique(round(df$w, 6))) > 1) {
    df$observed <- df$observed / df$w
    if (!is.null(df$fitted)) df$fitted <- df$fitted / df$w
    if (ylab == "Proportion") ylab <- "Proportion per cm" else ylab <- paste(ylab, "per cm")
  }
  df$mid <- df$len + df$w / 2
  sg <- ifelse(df$sex == "M", -1, 1)
  df$obs_s <- sg * df$observed
  if (!is.null(df$fitted)) df$fit_s <- sg * df$fitted
  df$sex <- factor(df$sex, levels = c("F", "M", "U", "C"))
  lab <- c(F = "Female", M = "Male", U = "Unsexed", C = "Sexes combined")
  lens <- sort(unique(df$len))
  x_breaks <- seq(floor(min(lens) / len_by) * len_by, ceiling(max(lens) / len_by) * len_by, by = len_by)
  df$year <- factor(df$year, levels = sort(unique(df$year)))
  df$panel <- factor(df$panel, levels = unique(df$panel))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = len)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = len, xmax = len + w, ymin = 0, ymax = obs_s, fill = sex),
                       alpha = 0.45, colour = NA, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.3)
  if (!is.null(df$fitted)) {
    lcol <- cols; lcol["U"] <- fitted_U
    p <- p + ggplot2::geom_line(ggplot2::aes(x = mid, y = fit_s, colour = sex,
                                             group = interaction(sex, panel, year)),
                                linewidth = linewidth, na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = lcol, labels = lab, drop = TRUE, name = "Fitted")
  }
  p <- p +
    ggplot2::scale_fill_manual(values = cols, labels = lab, drop = TRUE, name = "Observed") +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_y_continuous(labels = function(x) format(abs(x), trim = TRUE, drop0trailing = TRUE)) +
    ggplot2::labs(x = "Length", y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey90"),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  if (length(levels(df$panel)) > 1) {
    p + ggplot2::facet_grid(year ~ panel, scales = scales)
  } else {
    p + ggplot2::facet_wrap(~ year, scales = scales, ncol = ncol_years)
  }
}
