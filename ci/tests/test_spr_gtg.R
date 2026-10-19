# spr_gtg_flicc(): mat = "length" must reproduce the template SPR; "relative" vs LBSPR
library(FLicc)
if ("lbspr_test" %in% data(package = "FLicc")$results[, "Item"]) data("lbspr_test", package = "FLicc") else
  load(file.path(if (dir.exists("out/tests")) "out/tests" else ".", "lbspr_test.rda"))
set_lbspr <- list(CVL = 0.1, ngtg = 13, maxsd = 2, obs_model = "mn")
res <- lapply(dimnames(lfd_lbspr[[1]])$year, function(y) {
  lfd_y <- FLQuants(lapply(lfd_lbspr, function(x) x[, y]))
  f <- fiticc(lfd_y, stocklen(lfd_y, lhpar_lbspr, m_model = "constant"), sel_fun = "logistic",
              catch_by_gear = 1, settings = set_lbspr, verbose = FALSE)
  a <- spr_gtg_flicc(f, mat = "length"); r <- spr_gtg_flicc(f)
  data.frame(year = y, spr = a$spr, spr_len = a$spr_gtg, spr_rel = r$spr_gtg, adj = r$adj,
             lwr = r$lwr, upr = r$upr)
})
res <- do.call(rbind, res)
res$LBSPR <- lbspr_est$SPR
print(res, digits = 4)
cat("max |spr_len - spr| =", max(abs(res$spr_len - res$spr)), "\n")
cat("mean |spr - LBSPR| =", mean(abs(res$spr - res$LBSPR)),
    " mean |spr_rel - LBSPR| =", mean(abs(res$spr_rel - res$LBSPR)), "\n")
stopifnot(max(abs(res$spr_len - res$spr)) < 1e-6)
