# Illex (LBB example): FLicc variants against LBB's per-year results
library(FLicc)
options(width = 220)
od <- if (dir.exists("out/tests")) "out/tests" else "."
f <- file.path(od, "lbb_examples.rds")
if (!file.exists(f)) { cat("no lbb_examples.rds\n") } else {
  lbbx <- readRDS(f)
  lbb <- data.frame(Year = c(1996, 1999, 2001, 2002, 2009, 2014, 2016),
                    Lc = c(12.71, 5.23, 3.64, 4.33, 7.74, 5.32, 5.56),
                    MK = c(1.82, 1.66, 1.47, 1.49, 1.66, 1.41, 1.50),
                    ZK = c(2.32, 3.12, 2.57, 3.72, 3.15, 5.24, 6.85),
                    FM = c(0.265, 0.878, 0.732, 1.467, 0.910, 2.713, 3.570),
                    BB0 = c(0.6986, 0.2668, 0.2932, 0.1384, 0.2942, 0.0648, 0.0438))
  yrs <- lbb$Year
  d0 <- lbbx$lfd[lbbx$lfd$stock == "Ille_coi_AD", ]
  mk <- function(shift) FLQuantLen(lfd_long_to_wide(data.frame(len = d0$len + shift, year = d0$year,
                                   data = d0$data1000, qname = "Illex")), unit = "cm", midL = FALSE)
  run <- function(lab, shift = 0, ess = 20000, ...) {
    r <- suppressMessages(lbb_flicc(mk(shift), Linf = 23, Lm50 = 16, years = yrs, ess = ess, verbose = FALSE, ...))
    tb <- r$table
    data.frame(variant = lab, Year = tb$year, Lc = tb$Lc, MK = tb$MK, ZK = tb$ZK, FM = tb$FM, BB0 = tb$BB0,
               conv = vapply(as.character(tb$year), function(y) r$fits[[y]]$diagnostics$conv_tier, ""))
  }
  # slope fixed at LBB's estimated alpha per year: dSL = SL95 - SL50 = log(19) Linf / alpha
  lbb_alpha <- c(20.4, 36.7, 35.8, 42.7, 25.7, 32.4, 41.8)
  runa <- function(lab, ...) do.call(rbind, lapply(seq_along(yrs), function(i) {
    r <- suppressMessages(lbb_flicc(mk(0), Linf = 23, Lm50 = 16, years = yrs[i], ess = 20000, verbose = FALSE,
                                    sel_fixed = list(Illex = list(dSL = log(19) * 23 / lbb_alpha[i])), ...))
    tb <- r$table
    data.frame(variant = lab, Year = tb$year, Lc = tb$Lc, MK = tb$MK, ZK = tb$ZK, FM = tb$FM, BB0 = tb$BB0,
               conv = r$fits[[1]]$diagnostics$conv_tier)
  }))
  res <- rbind(
    run("A lbb weighting, no F/M cap"),
    run("B + lengths shifted -0.25 (half bin)", shift = -0.25),
    run("C + lengths shifted -0.5 (one bin)", shift = -0.5),
    run("D default weighting (estimated DM)", ess = NULL),
    run("E lbb weighting, M/K free (Mk.sd 0.5)", settings = list(Mk.sd = 0.5)),
    run("F lbb weighting, Pile = 0", pile = 0),
    runa("G lbb weighting, slope fixed at LBB alpha"),
    runa("H lbb weighting, Pile = 0, slope fixed at LBB alpha", pile = 0))
  cmp <- merge(res, lbb, by = "Year", suffixes = c("", "_LBB"))
  cmp <- cmp[order(cmp$variant, cmp$Year), ]
  print(format(cmp[, c("variant", "Year", "conv", "Lc", "Lc_LBB", "MK", "MK_LBB", "ZK", "ZK_LBB", "FM", "FM_LBB", "BB0", "BB0_LBB")], digits = 3), row.names = FALSE)
  sm <- aggregate(cbind(dZK = ZK - ZK_LBB, dFM = FM - FM_LBB, dBB0 = BB0 - BB0_LBB, dLc = Lc - Lc_LBB) ~ variant, cmp,
                  function(v) round(mean(v), 3))
  print(sm, row.names = FALSE)
  utils::write.csv(cmp, file.path(od, "illex_flicc_vs_lbb.csv"), row.names = FALSE)
}
