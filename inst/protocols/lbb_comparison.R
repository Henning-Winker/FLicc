# FLicc vs LBB on the LBB example stocks
# ---------------------------------------------------------------------------
# Runs lbb_flicc() on each stock of lbb_examples.rds with the settings of
# LBB's Example_ID.csv (Years.user, Linf.user, MK.user, Lc.user, Lstart.user,
# Lcut.user, GausSel, Lm50) and two data weightings:
#   "dm"     : FLicc default, raw data with Dirichlet-multinomial precision estimated
#   "lbb"    : LBB-like, counts rescaled to ESS 20000 with the DM precision fixed
#              (n_eff ~ 10^4, as LBB's Dirichlet concentration of 10 x 1000)
# Output: one row per stock, year and weighting, with the LBB quantities
# (Linf, Lc, alpha, M/K, F/K, F/M, Z/K, Lopt, Lc_opt, Lmean/Lopt, B/B0, Y/R')
# plus FLicc's EPR and SPR. Merge with LBB's own results by Stock and Year.
#
# Usage:
#   library(FLicc)
#   source(system.file("protocols", "lbb_comparison.R", package = "FLicc"))
#   lbbx <- readRDS("lbb_examples.rds")
#   res  <- run_lbb_examples(lbbx)                    # all stocks
#   res1 <- run_lbb_examples(lbbx, stocks = "Ille_coi_AD")
#   fit  <- lbb_example_fit(lbbx, "Ille_coi_AD")      # one stock, full object
#   plot_lbb(fit)
#   write.csv(res, "flicc_lbb_examples.csv", row.names = FALSE)
# ---------------------------------------------------------------------------

# LBB settings of one stock as lbb_flicc() arguments
lbb_id_args <- function(lbbx, stock) {
  id <- lbbx$id[lbbx$id$Stock == stock, , drop = FALSE]
  if (nrow(id) != 1) stop("Stock not found in lbbx$id: ", stock)
  num <- function(v) { v <- suppressWarnings(as.numeric(v)); if (length(v) == 0 || is.na(v)) NULL else v }
  yrs <- id$Years.user
  years <- if (!is.na(yrs) && nzchar(yrs) && yrs != "NA") as.numeric(strsplit(yrs, ",")[[1]]) else NULL
  if (is.null(years)) {
    sy <- num(id$StartYear); ey <- num(id$EndYear)
    ally <- sort(unique(lbbx$lfd$year[lbbx$lfd$stock == stock]))
    if (!is.null(sy)) ally <- ally[ally >= sy]
    if (!is.null(ey)) ally <- ally[ally <= ey]
    years <- ally
  }
  list(Linf = num(id$Linf.user), MK = num(id$MK.user), Lc = num(id$Lc.user),
       Lstart = num(id$Lstart.user), Lcut = num(id$Lcut.user), Lm50 = num(id$Lm50),
       sel_fun = if (isTRUE(as.logical(id$GausSel))) "normal" else "logistic",
       years = years)
}

# FLQuants of one stock from the long data
lbb_example_lfd <- function(lbbx, stock, raised = TRUE) {
  d <- lbbx$lfd[lbbx$lfd$stock == stock, ]
  FLQuantLen(lfd_long_to_wide(data.frame(len = d$len, year = d$year,
                                         data = if (raised) d$data1000 else d$data, qname = d$qname)),
             unit = "cm", midL = FALSE)
}

# lbb_flicc() for one stock with its LBB settings
lbb_example_fit <- function(lbbx, stock, weighting = c("dm", "lbb"), ...) {
  weighting <- match.arg(weighting)
  a <- lbb_id_args(lbbx, stock)
  lfd <- lbb_example_lfd(lbbx, stock)
  do.call(lbb_flicc, c(list(lfd = lfd), a,
                       list(ess = if (weighting == "lbb") 20000 else NULL, verbose = FALSE), list(...)))
}

# All (or selected) stocks, both weightings, one table
run_lbb_examples <- function(lbbx, stocks = NULL, weightings = c("dm", "lbb")) {
  if (is.null(stocks)) stocks <- unique(lbbx$id$Stock)
  stocks <- intersect(stocks, unique(lbbx$lfd$stock))
  out <- list()
  for (st in stocks) for (w in weightings) {
    t0 <- Sys.time()
    r <- tryCatch(suppressMessages(lbb_example_fit(lbbx, st, w)), error = function(e) e)
    sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    if (inherits(r, "error") || is.null(r$table)) {
      msg <- if (inherits(r, "error")) conditionMessage(r) else "no fits"
      message(st, " [", w, "]: failed: ", msg)
      out[[length(out) + 1]] <- data.frame(Stock = st, weighting = w, Year = NA, error = msg, sec = sec)
      next
    }
    tb <- r$table
    conv <- vapply(as.character(tb$year), function(y) {
      f <- r$fits[[y]]; if (is.null(f)) NA_character_ else f$diagnostics$conv_tier }, character(1))
    out[[length(out) + 1]] <- data.frame(
      Stock = st, weighting = w, Year = tb$year, conv = conv,
      Linf = tb$Linf, Lc = tb$Lc, alpha = tb$alpha, MK = tb$MK, FK = tb$FK, FM = tb$FM,
      FM_lcl = tb$FM_lcl, FM_ucl = tb$FM_ucl, ZK = tb$ZK, Lopt = tb$Lopt, Lc_opt = tb$Lc_opt,
      Lmean_Lopt = tb$Lmean_Lopt, Lc_Lc_opt = tb$Lc_Lc_opt, L95 = tb$L95, perc_mat = tb$perc_mat,
      BB0 = tb$BB0, BB0_lcl = tb$BB0_lcl, BB0_ucl = tb$BB0_ucl, YR = tb$YR,
      BFM1B0 = attr(tb, "BFM1B0"), EPR = tb$EPR, SPR = tb$SPR,
      Linf_prior = r$priors$Linf, MK_prior = r$priors$MK, error = NA_character_, sec = sec)
  }
  res <- do.call(rbind, lapply(out, function(x) { for (nm in setdiff(names(out[[which.max(sapply(out, ncol))]]), names(x))) x[[nm]] <- NA; x }))
  res[, names(out[[which.max(sapply(out, ncol))]])]
}
