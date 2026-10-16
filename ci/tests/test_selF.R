# selF surfaces: consistency with the fit and timing
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
lfd200 <- lfdess(lfd, ess.g = 200)
for (pm in c("gtg", "gamma")) {
  fit <- fiticc(lfd200, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
                settings = list(CVL = 0.1, obs_model = "dm", pop_model = pm))
  st <- FLicc:::gear_state_flicc(fit)
  for (g in st$gears) {
    s1 <- sel_curve_flicc(st$Lmid, st$type[[g]], st$pars[[g]]); s2 <- as.numeric(fit$report$sel_gear[[g]])
    cat(pm, g, "max |sel diff|", signif(max(abs(s1 - s2)), 3), "\n"); stopifnot(max(abs(s1 - s2)) < 1e-3)
  }
  t <- system.time(sf <- selF_surface_flicc(fit))[["elapsed"]]
  cat(pm, "surface seconds:", t, "\n"); print(sf$current)
  spr_rep <- tail(c(fit$report$spr), 1)
  ypr_ref <- ypr_flicc(fit, F = as.numeric(tail(c(fit$report$Fap), 1)))
  cat(pm, "SPR engine", sf$current$SPR[1], "report", spr_rep, " YPR engine", sf$current$YPR[1], "ypr_flicc", ypr_ref, "\n")
  stopifnot(abs(sf$current$SPR[1] - spr_rep) < 1e-4, abs(sf$current$SPR[2] - spr_rep) < 1e-4,
            abs(sf$current$YPR[1] / ypr_ref - 1) < 1e-3)
  print(summary(sf$surface[, c("loc", "F", "SPR", "YPRrel")]))
  p <- plot_selF_flicc(sf)
  od <- if (dir.exists("out/tests")) "out/tests" else "."
  ggplot2::ggsave(file.path(od, paste0("selF_", pm, ".png")), p, width = 9, height = 8, dpi = 110)
  ggplot2::ggsave(file.path(od, paste0("selF_", pm, "_gillnet.png")), plot_selF_flicc(sf, gear = "Gillnet"),
                  width = 9, height = 4.4, dpi = 110)
}
