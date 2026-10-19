# Build lbb_examples.rds: LBB example data (github.com/SISTA16/LBB) in long format for FLicc
library(FLicc)
od <- if (dir.exists("out/tests")) "out/tests" else "."
base <- "https://raw.githubusercontent.com/SISTA16/LBB/master/"
dat <- tryCatch(utils::read.csv(paste0(base, "Example_Dat.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
id  <- tryCatch(utils::read.csv(paste0(base, "Example_ID.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(dat) || is.null(id)) {
  cat("LBB example data not reachable; skipped.\n")
} else {
  dat <- dat[!is.na(dat$CatchNo), ]
  # LBB: lengths in the catch file are always in mm; with mm.user = FALSE the analysis is in cm
  mm <- stats::setNames(as.logical(id$mm.user), id$Stock)
  dat$len <- ifelse(!is.na(mm[dat$Stock]) & mm[dat$Stock], dat$Length, dat$Length / 10)
  lfd <- aggregate(CatchNo ~ Stock + Year + len, data = dat, FUN = sum)
  names(lfd) <- c("stock", "year", "len", "data")
  lfd <- lfd[order(lfd$stock, lfd$year, lfd$len), ]
  tot <- ave(lfd$data, lfd$stock, lfd$year, FUN = sum)
  lfd$data1000 <- lfd$data / tot * 1000                     # raised to 1000 per stock and year
  lfd$qname <- lfd$stock
  lfd$species <- id$Species[match(lfd$stock, id$Stock)]
  lfd$unit <- ifelse(!is.na(mm[lfd$stock]) & mm[lfd$stock], "mm", "cm")
  rownames(lfd) <- NULL
  idx <- id[, c("Stock", "Species", "StartYear", "EndYear", "Years.user", "Year.select", "Lcut.user", "Lc.user",
                "Lstart.user", "Linf.user", "MK.user", "mm.user", "GausSel", "MergeLF", "Pile", "Lm50", "Comment")]
  lbb_examples <- list(lfd = lfd, id = idx,
                       source = "LBB example data, https://github.com/SISTA16/LBB (Example_Dat.csv, Example_ID.csv); Froese et al. 2018, ICES JMS 75(6): 2004-2015",
                       notes = paste("len: lower bound of the length class in cm (mm for stocks with mm.user = TRUE);",
                                     "data: CatchNo as in LBB (summed over duplicate records); data1000: data raised to 1000 per stock and year;",
                                     "id: LBB settings per stock (Years.user etc.); qname = stock for lfd_long_to_wide()."))
  saveRDS(lbb_examples, file.path(od, "lbb_examples.rds"))
  print(aggregate(data ~ stock, lfd, function(x) length(x)))
  print(table(lfd$unit))
  # conversion check: one stock to FLQuants, raised to 1000
  d <- lfd[lfd$stock == "Ille_coi_AD", ]
  il <- FLQuantLen(lfd_long_to_wide(data.frame(len = d$len, year = d$year, data = d$data1000, qname = d$qname)),
                   unit = "cm", midL = FALSE)
  print(il); print(apply(il[[1]], 2, sum))
  stopifnot(all(abs(apply(il[[1]]@.Data, 2, sum, na.rm = TRUE) - 1000) < 1e-6))
}
