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
  t <- numeric(n)
  for (i in seq_len(n))
    t[i] <- system.time(f <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg,
                                    settings = set, se = se))[["elapsed"]]
  list(fit = f, time = median(t))
}
res <- list(
  sd2_13 = time_fit(list(maxsd = 2, ngtg = 13)),
  sd3_19 = time_fit(list(maxsd = 3, ngtg = 19)),
  sd3_13 = time_fit(list(maxsd = 3, ngtg = 13)),
  sd3_15 = time_fit(list(maxsd = 3, ngtg = 15)),
  sd2_13_noSE = time_fit(list(maxsd = 2, ngtg = 13), se = FALSE),
  sd3_19_noSE = time_fit(list(maxsd = 3, ngtg = 19), se = FALSE))
tab <- data.frame(run = names(res),
                  time_s = sapply(res, `[[`, "time"),
                  SPR = sapply(res, function(r) sprcur_flicc(r$fit)),
                  nll = sapply(res, function(r) r$fit$opt$objective),
                  conv = sapply(res, function(r) paste(flicc_convergence(r$fit)$passed, collapse = ",")))
print(tab, row.names = FALSE)
tab$time_rel <- round(tab$time_s / tab$time_s[1], 2)
print(tab[, c("run", "time_rel", "SPR", "nll")], row.names = FALSE)
# predicted length composition of the gillnet in the tail: smoothness with fewer groups
tail_pred <- try(sapply(res[c("sd2_13", "sd3_19", "sd3_15", "sd3_13")], function(r) {
  q <- r$fit$report$predlen_gear[[2]]
  p <- apply(q@.Data, 1, mean); L <- as.numeric(dimnames(q)$len)
  round(1000 * p[L >= 45 & L <= 62] / sum(p), 3)            # per mille, gillnet, 45-62 cm
}))
print(tail_pred)
stopifnot(abs(tab$SPR[1] - tab$SPR[2]) < 0.05)
