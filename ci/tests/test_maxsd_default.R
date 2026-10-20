# Default GTG grid (maxsd = 3, ngtg = 19, FLicc >= 1.0.11) vs LBSPR's (2, 13): results and speed
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)[, c("len", "year", "data", "qname")]),
                  unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)

f_def <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg)
stopifnot(f_def$tmb_data$ngtg == 19, length(f_def$tmb_data$gtg_z) == 19, max(f_def$tmb_data$gtg_z) == 3)

time_fit <- function(set, se = TRUE, n = 3) {
  t <- replicate(n, system.time(f <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg,
                                             settings = set, se = se))[["elapsed"]])
  list(fit = f, time = median(t))
}
res <- list(
  sd2_13 = time_fit(list(maxsd = 2, ngtg = 13)),
  sd3_19 = time_fit(list(maxsd = 3, ngtg = 19)),
  sd2_13_noSE = time_fit(list(maxsd = 2, ngtg = 13), se = FALSE),
  sd3_19_noSE = time_fit(list(maxsd = 3, ngtg = 19), se = FALSE))
tab <- data.frame(run = names(res),
                  time_s = sapply(res, `[[`, "time"),
                  SPR = sapply(res, function(r) sprcur_flicc(r$fit)),
                  nll = sapply(res, function(r) r$fit$opt$objective),
                  conv = sapply(res, function(r) paste(flicc_convergence(r$fit)$passed, collapse = ",")))
print(tab, row.names = FALSE)
cat("time ratio 3/19 vs 2/13: with SE", round(tab$time_s[2] / tab$time_s[1], 2),
    " without SE", round(tab$time_s[4] / tab$time_s[3], 2), "\n")
stopifnot(abs(tab$SPR[1] - tab$SPR[2]) < 0.05)
