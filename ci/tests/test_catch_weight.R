# catch_by_gear is matched in weight: predicted shares (numbers and weight) and SPR
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
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
for (case in list(list(nm = "raw mn", lfd = lfd, set = list(CVL = 0.1)),
                  list(nm = "ess200 dm", lfd = lfdess(lfd, ess.g = 200), set = list(CVL = 0.1, obs_model = "dm")))) {
  f <- fiticc(case$lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = case$set)
  cs <- f$obj$report(f$obj$env$last.par.best)$catch_share
  cat("\n==", case$nm, " conv", f$opt$convergence, "maxgrad", signif(f$max_gradient, 2),
      " SPR last", round(sprcur_flicc(f), 3), " F by gear (last yr)", round(unlist(FLicc:::gear_state_flicc(f)$F), 4), "\n")
  cat("reported catch_share (weight), mean over years:", round(colMeans(cs), 3), "\n")
  print(shares(f))
  print(round(c(f$report$spr), 3))
  stopifnot(max(abs(sweep(cs, 2, cbg))) < 0.05)
}
