# plot_growth(), growth_phase_flicc(), plot_growth_lfd()
library(FLicc)
od <- if (dir.exists("out/tests")) "out/tests" else "."
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)[, c("len", "year", "data", "qname")]),
                  unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")

# 1. plot_growth before and after a fit
p1 <- plot_growth(lhpar, CVL = 0.1, lfd = lfd)
p2 <- plot_growth(lhpar, CVL = 0.1, lfd = lfd, composition = TRUE)
fit <- fiticc(lfd, stklen, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3))
p3 <- plot_growth(fit)
ggplot2::ggsave(file.path(od, "growth.png"), p1, width = 8, height = 5, dpi = 90)

# 2. recovery of a known phase from simulated mixtures (free weights)
set.seed(3)
sim_phase <- 0.35; k <- 0.25; linf <- 50; cv <- 0.08
L <- 0:60
sim <- FLQuants(G = FLQuant(0, dimnames = list(len = L, year = 2001:2006)))
for (y in 2001:2006) {
  mu <- linf * (1 - exp(-k * (0:12 + sim_phase))); sdv <- sqrt((cv * mu)^2 + 1)
  w <- exp(-0.4 * (0:12)) * exp(rnorm(13, 0, 0.5))
  pr <- sapply(seq_along(mu), function(j) pnorm(L + 1, mu[j], sdv[j]) - pnorm(L, mu[j], sdv[j])) %*% w
  sim$G[, as.character(y)] <- rmultinom(1, 2000, pr)
}
lh_sim <- FLPar(linf = linf, k = k, M = 0.3, L50 = 30)
gp_sim <- growth_phase_flicc(lh_sim, lfd = sim, CVL = cv / 0.99, sd0 = 1)
print(gp_sim$best); print(gp_sim$t0)
stopifnot(abs(gp_sim$best$phase - sim_phase) <= 0.06)
gp_k <- growth_phase_flicc(lh_sim, lfd = sim, CVL = cv / 0.99, sd0 = 1, k = c(0.15, 0.2, 0.25, 0.3, 0.35),
                           phase = seq(0, 0.95, 0.05))
print(gp_k$best)
stopifnot(gp_k$best$k == 0.25)

# 3. alfonsino example
t_g <- system.time(gp <- growth_phase_flicc(fit))
print(t_g)
print(gp$best); print(gp$t0); cat("contrast", gp$contrast, "\n")
print(round(gp$profile$dLL, 1))
ggplot2::ggsave(file.path(od, "growth_phase.png"), plot_growth_phase(gp), width = 6, height = 4, dpi = 90)
ggplot2::ggsave(file.path(od, "growth_lfd.png"), plot_growth_lfd(gp), width = 9, height = 7, dpi = 90)
gpk <- growth_phase_flicc(fit, k = c(0.06, 0.07, 0.08, 0.09, 0.1, 0.12), phase = seq(0, 0.9, 0.1))
print(gpk$best)
ggplot2::ggsave(file.path(od, "growth_phase_k.png"), plot_growth_phase(gpk), width = 6, height = 4, dpi = 90)
