# Newton polishing and DM prior SD: old (newton = 0, dm_prior_sd = 2) vs new (newton = 3, dm_prior_sd = 1)
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
se_last <- function(f) { ci <- spr_ci_flicc(f); if (is.null(ci)) NA else tail(ci$se_log, 1) }
ess_values <- c(50, 100, 200, 500, 1000, NA)          # NA = raw counts
variants <- list(old = list(newton = 0, sd = 2), new = list(newton = 3, sd = 1))
res <- list()
for (dt in c("saturating", "linear", "fixed")) for (ess in ess_values) for (v in names(variants)) {
  d <- if (is.na(ess)) lfd else lfdess(lfd, ess.g = ess)
  vv <- variants[[v]]
  t0 <- system.time(f <- tryCatch(
    fiticc(d, stklen, sel_fun = sf, catch_by_gear = cbg, newton = vv$newton, verbose = FALSE,
           settings = list(CVL = 0.1, obs_model = "dm", dm_theta = dt, dm_prior_sd = vv$sd)),
    error = function(e) e))[["elapsed"]]
  if (inherits(f, "error")) { cat(dt, ess, v, "ERROR", conditionMessage(f), "\n"); next }
  ne <- neff_flicc(f)
  ld <- f$obj$env$parList(par = f$obj$env$last.par.best)$log_dm
  dg <- f$diagnostics
  res[[length(res) + 1]] <- data.frame(
    dm = dt, ESS = ifelse(is.na(ess), "raw", ess), variant = v,
    conv = f$opt$convergence, nlminb = if (is.null(f$opt$convergence_nlminb)) f$opt$convergence else f$opt$convergence_nlminb,
    pdHess = isTRUE(f$rep$pdHess), tier = dg$conv_tier, newton = dg$newton_steps,
    restarts = f$n_restarts, maxgrad = signif(f$max_gradient, 2), decr = signif(dg$decrement, 2),
    SPR_last = round(sprcur_flicc(f), 4), se_log = round(se_last(f), 3),
    par_T = signif(exp(ld[1]), 3), par_G = signif(exp(ld[2]), 3),
    neff_T = round(mean(ne$neff[ne$gear == "Trawl"]), 1),
    neff_G = round(mean(ne$neff[ne$gear == "Gillnet"]), 1),
    at_bound = if (is.null(dg$dm)) NA else paste(dg$dm$at_bound, collapse = "/"),
    sec = round(t0, 2))
}
tab <- do.call(rbind, res)
options(width = 250)
print(tab, row.names = FALSE)
saveRDS(tab, if (dir.exists("out/tests")) "out/tests/newton_dm.rds" else "newton_dm.rds")

# speed without standard errors (as in an MSE), multinomial and DM, 3 repeats each
sp <- list()
for (om in c("mn", "dm")) for (ess in c(200, NA)) for (v in names(variants)) {
  d <- if (is.na(ess)) lfd else lfdess(lfd, ess.g = ess)
  vv <- variants[[v]]
  st <- list(CVL = 0.1, obs_model = om, dm_prior_sd = vv$sd)
  if (om == "dm") st$dm_theta <- "saturating"
  tt <- replicate(3, system.time(f <<- fiticc(d, stklen, sel_fun = sf, catch_by_gear = cbg, se = FALSE,
                                             newton = vv$newton, verbose = FALSE, settings = st))[["elapsed"]])
  sp[[length(sp) + 1]] <- data.frame(obs = om, ESS = ifelse(is.na(ess), "raw", ess), variant = v,
                                     sec = round(median(tt), 3), restarts = f$n_restarts,
                                     maxgrad = signif(f$max_gradient, 2), tier = f$diagnostics$conv_tier,
                                     hess_pd = f$diagnostics$hess_pd, conv = f$opt$convergence)
}
print(do.call(rbind, sp), row.names = FALSE)

# checks
new <- tab[tab$variant == "new", ]; old <- tab[tab$variant == "old", ]
key <- function(x) paste(x$dm, x$ESS)
m <- match(key(new), key(old))
cat("max |SPR new - old|, raw counts:",
    max(abs(new$SPR_last[new$ESS == "raw"] - old$SPR_last[m][new$ESS == "raw"])), "\n")
# Newton steps never increase the objective: compare objectives of a fit with and without
f0 <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, newton = 0, verbose = FALSE,
             settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
f3 <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, verbose = FALSE,
             settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
stopifnot(f3$opt$objective <= f0$opt$objective + 1e-8)
print(flicc_convergence(f3))
# messages: printed with verbose = TRUE, silent with verbose = FALSE
d200 <- lfdess(lfd, ess.g = 200)
msg <- character()
withCallingHandlers(fiticc(d200, stklen, sel_fun = sf, catch_by_gear = cbg,
                           settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating")),
                    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("messages (verbose = TRUE):\n"); print(msg)
msg2 <- character()
withCallingHandlers(fiticc(d200, stklen, sel_fun = sf, catch_by_gear = cbg, verbose = FALSE,
                           settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating")),
                    message = function(m) { msg2 <<- c(msg2, conditionMessage(m)); invokeRestart("muffleMessage") })
stopifnot(length(msg2) == 0, any(grepl("upper bound", msg)))
