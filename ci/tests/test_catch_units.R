# catch_by_gear in numbers (default) vs weight
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
lfd200 <- lfdess(lfd, ess.g = 200)
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)

# catch share of each gear (numbers and weight) at the fitted state, last year
shares <- function(fit) {
  eng <- FLicc:::pr_engine_flicc(fit); st <- FLicc:::gear_state_flicc(fit)
  S <- lapply(st$gears, function(g) sel_curve_flicc(st$Lmid, st$type[[g]], st$pars[[g]]))
  Fl <- Reduce(`+`, Map(function(f, s) f * s, st$F, S))
  d <- eng(Fl, detail = TRUE)
  w <- as.numeric(FLCore::yearMeans(fit$report$w_at_len[, tail(dimnames(fit$report$N)$year, 1)]))
  Cn <- d$Cw / w
  out <- sapply(seq_along(st$gears), function(i) {
    part <- ifelse(Fl > 0, st$F[[i]] * S[[i]] / Fl, 0)
    c(numbers = sum(Cn * part) / sum(Cn), weight = sum(d$Cw * part) / sum(d$Cw))
  })
  colnames(out) <- st$gears; round(out, 3)
}
res <- list()
for (om in c("mn", "dm")) for (cu in c("numbers", "weight")) {
  nm <- paste(om, cu, sep = "_")
  f <- fiticc(lfd200, stklen, sel_fun = sf, catch_by_gear = cbg,
              settings = list(CVL = 0.1, obs_model = om, catch_units = cu))
  res[[nm]] <- f
  cat("\n==", nm, " conv", f$opt$convergence, "maxgrad", signif(f$max_gradient, 2),
      " SPR last", round(sprcur_flicc(f), 3), " F by gear", round(unlist(FLicc:::gear_state_flicc(f)$F), 4), "\n")
  print(shares(f))
}
cat("\nSPR by year\n")
print(round(sapply(res, function(f) c(f$report$spr)), 3))
# default unchanged
stopifnot(isTRUE(all.equal(res$dm_numbers$opt$objective,
  fiticc(lfd200, stklen, sel_fun = sf, catch_by_gear = cbg, settings = list(CVL = 0.1, obs_model = "dm"))$opt$objective)))
# with weight, the predicted weight shares match the input closely
sw <- shares(res$dm_weight)
stopifnot(abs(sw["weight", "Trawl"] - 0.7) < 0.05)
od <- if (dir.exists("out/tests")) "out/tests" else "."
ggplot2::ggsave(file.path(od, "catch_units_spr.png"), plot_spr(res, ci = TRUE), width = 8, height = 4, dpi = 110)
