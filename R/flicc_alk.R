# Growth from age-length keys (conditional age at length)

#' Age-length key for growth estimation in FLicc
#'
#' Aged fish of one gear and year, used in [fiticc()] (`alk = `) to estimate von Bertalanffy growth (Linf, k)
#' jointly with the length compositions. The likelihood is conditional on length, P(age | length), so fish can be
#' selected for ageing at random, stratified by length class or opportunistically by length: the selection drops out.
#' The predicted P(age | length) accounts for length-based selectivity (it cancels at a given length) and for the
#' removal of fast growers by fishing (Rosa Lee effect) through the survival along the length classes. It carries
#' growth information, not information on total mortality.
#'
#' @param data data frame with `age` and `len` (one row per fish), or `age`, `len` and `n` (counts). `len` in the
#'   units of the length compositions; fish are assigned to the length classes of the fit.
#' @param gears gear(s) the fish were sampled from (for diagnostics; P(age | length) does not depend on the gear
#'   when selectivity is length-based).
#' @param year year of the sample (one of the fitted years).
#' @param timing sampling time within the year (fraction of the year; 0.5 = mid-year, the mean of year-round sampling).
#' @return an object of class `flicc_alk`.
#' @export
alk_flicc <- function(data, gears, year, timing = 0.5) {
  data <- as.data.frame(data)
  if (!all(c("age", "len") %in% names(data))) stop("data needs columns 'age' and 'len'.")
  if (is.null(data$n)) data$n <- 1
  data <- data[is.finite(data$age) & is.finite(data$len) & data$n > 0, c("age", "len", "n")]
  if (!nrow(data)) stop("no aged fish.")
  if (length(year) != 1) stop("one year per key: use a list of keys for several years or gears.")
  structure(list(data = data, gears = gears, year = as.numeric(year), timing = as.numeric(timing)),
            class = "flicc_alk")
}

# TMB data of the ALK block (defaults when there is no ALK)
alk_tmb_data_flicc <- function(alk, LLB, years, lhpar, settings) {
  none <- list(alk_on = 0L, alk_mkk = 0L, alk_logM = 0, alk_t0 = 0, alk_maxsd = 3,
               alk_gl_x = 0, alk_gl_w = 2, alk_age = 0, alk_year = 0L, alk_tau = 0.5,
               alk_obs = array(0, c(1, 1, 1)))
  if (is.null(alk)) return(none)
  if (inherits(alk, "flicc_alk")) alk <- list(alk)
  if (!all(vapply(alk, inherits, TRUE, "flicc_alk"))) stop("alk must be an alk_flicc() object or a list of them.")
  nlen <- length(LLB)
  bw <- if (nlen > 1) LLB[nlen] - LLB[nlen - 1] else 1
  ub <- c(LLB[-1], LLB[nlen] + bw)
  amax <- max(vapply(alk, function(k) max(k$data$age), 0)) + 2
  ages <- 0:amax
  obs <- array(0, c(length(ages), nlen, length(alk)))
  yi <- integer(length(alk))
  for (j in seq_along(alk)) {
    d <- alk[[j]]$data
    li <- findInterval(d$len, LLB)
    li[d$len >= ub[nlen]] <- 0L
    ok <- li >= 1
    if (any(!ok)) warning(sum(d$n[!ok]), " aged fish outside the length classes ignored (key ", j, ").", call. = FALSE)
    for (r in which(ok)) obs[d$age[r] + 1, li[r], j] <- obs[d$age[r] + 1, li[r], j] + d$n[r]
    yi[j] <- match(alk[[j]]$year, years) - 1L
    if (is.na(yi[j])) stop("ALK year ", alk[[j]]$year, " is not a fitted year.")
  }
  gq <- statmod::gauss.quad(if (!is.null(settings$alk_nq)) settings$alk_nq else 5, "legendre")
  t0 <- if (!is.null(settings$alk_t0)) settings$alk_t0 else {
    warning("settings$alk_t0 not given: t0 = 0 assumed for the ALK.", call. = FALSE); 0
  }
  fixM <- if (is.null(settings$alk_fixM)) TRUE else isTRUE(settings$alk_fixM)
  list(alk_on = 1L, alk_mkk = as.integer(fixM), alk_logM = log(as.numeric(lhpar["M"])),
       alk_t0 = as.numeric(t0), alk_maxsd = if (!is.null(settings$maxsd)) as.numeric(settings$maxsd) else 3,
       alk_gl_x = as.numeric(gq$nodes), alk_gl_w = as.numeric(gq$weights), alk_age = as.numeric(ages),
       alk_year = as.integer(yi), alk_tau = vapply(alk, function(k) k$timing, 0), alk_obs = obs)
}

# defaults for fits made before the ALK option
alk_defaults_tmb_flicc <- function(td) {
  d <- alk_tmb_data_flicc(NULL)
  for (nm in names(d)) if (is.null(td[[nm]])) td[[nm]] <- d[[nm]]
  td
}

#' Fit of the age-length keys
#'
#' Observed and predicted mean age by length class for each key of a fit with `alk`.
#' @param fit a [fiticc()] fit with `alk`.
#' @return data frame: key, length, observed and predicted mean age, number aged.
#' @export
alk_fit_flicc <- function(fit) {
  td <- fit$tmb_data
  if (is.null(td$alk_on) || td$alk_on != 1L) stop("fit has no ALK.")
  growth <- fit$growth
  Linf <- growth[["linf"]]; k <- growth[["k"]]; CVL <- fit$settings$CVL
  out <- list()
  for (j in seq_len(dim(td$alk_obs)[3])) {
    o <- td$alk_obs[, , j]
    n <- colSums(o)
    obs_mean <- ifelse(n > 0, colSums(o * td$alk_age) / n, NA)
    # predicted mean age at length (no mortality): from the growth curve and the Linf distribution
    g <- 1 - exp(-k * (td$alk_age + td$alk_tau[j] - td$alk_t0))
    Lm <- td$LLB + (td$LLB[2] - td$LLB[1]) / 2
    p <- sapply(Lm, function(L) { d <- stats::dnorm(L / pmax(g, 1e-4), Linf, CVL * Linf) / pmax(g, 1e-4); d / sum(d) })
    out[[j]] <- data.frame(key = j, len = td$LLB, n = n, obs = obs_mean,
                           pred_nomort = colSums(p * td$alk_age))
  }
  do.call(rbind, out)
}
