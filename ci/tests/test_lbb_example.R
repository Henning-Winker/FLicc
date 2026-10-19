# FLicc vs LBB protocol on all LBB example stocks (needs lbb_examples.rds from test_lbb_data.R)
library(FLicc)
options(width = 250)
od <- if (dir.exists("out/tests")) "out/tests" else "."
f <- file.path(od, "lbb_examples.rds")
if (!file.exists(f)) {
  cat("lbb_examples.rds not available; skipped.\n")
} else {
  source("inst/protocols/lbb_comparison.R")
  lbbx <- readRDS(f)
  t0 <- system.time(res <- run_lbb_examples(lbbx))
  cat("total seconds", t0[["elapsed"]], "\n")
  utils::write.csv(res, file.path(od, "flicc_lbb_examples.csv"), row.names = FALSE)
  # overview: per stock and weighting, median over years
  ok <- res[is.na(res$error), ]
  agg <- aggregate(cbind(Linf, Lc, MK, FM, BB0, EPR, SPR) ~ Stock + weighting, ok, median)
  agg$n_years <- aggregate(Year ~ Stock + weighting, ok, length)$Year
  agg$n_ok <- aggregate(conv ~ Stock + weighting, ok, function(v) sum(v == "ok", na.rm = TRUE))$conv
  print(format(agg, digits = 3), row.names = FALSE)
  print(res[!is.na(res$error), c("Stock", "weighting", "error")])
  for (st in c("Ille_coi_AD", "Micr_pou_AD", "ANK_GSA22")) {
    fit <- tryCatch(suppressMessages(lbb_example_fit(lbbx, st)), error = function(e) NULL)
    if (!is.null(fit)) ggplot2::ggsave(file.path(od, paste0("lbbx_", st, ".png")), plot_lbb(fit), width = 10, height = 7, dpi = 90)
  }
}
