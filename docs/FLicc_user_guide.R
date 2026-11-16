# FLicc user guide: R code
# All code of the user guide (docs/FLicc_user_guide.html), chapter by chapter, extracted with
# knitr::purl() from inst/guide/FLicc_user_guide.Rmd. Run from the top; later chapters use objects
# made in earlier ones. Lines starting with "## " are code not run in the guide (eval = FALSE).
# Regenerate: setwd("inst/guide"); have_lbb <- TRUE
#             knitr::purl("FLicc_user_guide.Rmd", output = "../../docs/FLicc_user_guide.R", documentation = 1)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 7, fig.height = 4.5, message = FALSE, warning = FALSE)


## ----child = "chapters/01_introduction.Rmd"-----------------------------------




## ----child = "chapters/02_installation.Rmd"-----------------------------------

## ----install, eval = FALSE----------------------------------------------------
## install.packages("FLicc", repos = c("https://henning-winker.r-universe.dev",
##                                     "https://flr.r-universe.dev",
##                                     "https://cloud.r-project.org"))


## ----install-github, eval = FALSE---------------------------------------------
## install.packages("TMB", type = "source")
## install.packages(c("FLCore", "ggplotFL"),
##                  repos = c("https://flr.r-universe.dev", "https://cloud.r-project.org"))
## remotes::install_github("henning-winker/FLicc")


## ----load---------------------------------------------------------------------
library(FLicc)
packageVersion("FLicc")



## ----child = "chapters/03_quick_start.Rmd"------------------------------------

## ----qs-data------------------------------------------------------------------
library(FLicc)
data("alfonsino")
lfd_df <- as.data.frame(lfd_alfonsino)[, c("len", "year", "data", "qname")]
head(lfd_df)
unique(lfd_df$qname)       # gears
range(lfd_df$len)          # length range (cm)
range(lfd_df$year)         # years


## ----qs-own-data, eval = FALSE------------------------------------------------
## lfd_df <- read.csv("my_length_data.csv")              # columns: length, year, count, fleet
## names(lfd_df) <- c("len", "year", "data", "qname")


## ----qs-flquants--------------------------------------------------------------
lfd <- FLQuantLen(lfd_long_to_wide(lfd_df), unit = "cm", midL = FALSE)
names(lfd)


## ----qs-plot-lfd, fig.height = 5----------------------------------------------
plot_lfd(lfd)


## ----qs-lhpar-----------------------------------------------------------------
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859,
               a = 0.004721956 / 1000, b = 3.146168)


## ----qs-stocklen, fig.height = 3.5--------------------------------------------
stklen <- stocklen(lfd, lhpar, m_model = "constant")
plot_m(stklen)
plot_mat(stklen)


## ----qs-fit-------------------------------------------------------------------
fit1 <- fiticc(lfd, stklen,
               sel_fun = c("dsnormal", "logistic"),
               catch_by_gear = c(0.95, 0.05))


## ----qs-convergence-----------------------------------------------------------
flicc_convergence(fit1)


## ----qs-plot-len, fig.height = 5----------------------------------------------
plot_len(fit1)


## ----qs-plot-lengear, fig.height = 5------------------------------------------
plot_len(fit1, by_gear = TRUE,year=2020:2024)


## ----qs-plot-sel, fig.height = 3.5--------------------------------------------
plot_sel(fit1)


## ----qs-plot-spr, fig.height = 3.5--------------------------------------------
plot_spr(fit1, ci = TRUE)
sprcur_flicc(fit1)



## ----child = "chapters/04_workflow.Rmd"---------------------------------------

## ----reference-fit------------------------------------------------------------
sel_fun <- c("dsnormal", "logistic")
cbg <- c(0.95, 0.05)
fit <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
              settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
sprcur_flicc(fit)



## ----child = "chapters/05_length_data.Rmd"------------------------------------

## ----data---------------------------------------------------------------------
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
names(lfd)
range(as.numeric(dimnames(lfd[[1]])$len)); range(as.numeric(dimnames(lfd[[1]])$year))


## ----check-n, fig.height = 5--------------------------------------------------
plot_lfdn(lfd)


## ----check-lfd, fig.height = 5------------------------------------------------
plot_lfd(lfd, type = "relmax")


## ----check-tail---------------------------------------------------------------
linf_in <- 55.7
sapply(lfd, function(x) {
  L <- as.numeric(dimnames(x)$len); n <- apply(x@.Data, 1, sum)
  c(max_length = max(L[n > 0]), pct_ge_0.9Linf = round(100 * sum(n[L >= 0.9 * linf_in]) / sum(n), 1),
    pct_ge_Linf = round(100 * sum(n[L >= linf_in]) / sum(n), 2))
})


## ----check-zlcc, fig.height = 5-----------------------------------------------
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859,
               a = 0.004721956 / 1000, b = 3.146168)
plot_zlcc(lfd$Trawl, lhpar)
plot_zlcc(lfd$Gillnet, lhpar)


## ----missing------------------------------------------------------------------
lfd_miss <- lfd
lfd_miss[["Gillnet"]][, c("2016", "2017","2022")] <- 0
lfd_miss[["Trawl"]][, c("2020","2022")] <- 0
fit_miss <- fiticc(lfd_miss, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                   settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
flicc_convergence(fit_miss)


## ----missing-plot, fig.height = 4---------------------------------------------
plot_lfdn(lfd_miss)
plot_spr(list(complete = fit, gaps = fit_miss), ci = TRUE)



## ----child = "chapters/06_life_history.Rmd"-----------------------------------

## ----stocklen-----------------------------------------------------------------
stklen


## ----growth, fig.height = 4.5-------------------------------------------------
plot_growth(lhpar, CVL = 0.1, lfd = lfd, composition = TRUE)


## ----growth-phase, fig.height = 3.5-------------------------------------------
gph <- growth_phase_flicc(fit, spawn = 0, sample = 0.5)   # spawning at the start, sampling mid-year
gph$best
gph$t0
plot_growth_phase(gph)


## ----growth-lfd, fig.height = 6-----------------------------------------------
plot_growth_lfd(gph)
plot_growth_lfd(gph, ages = c(2, 12))     # selected ages only


## ----growth-sim, fig.height = 3.5---------------------------------------------
set.seed(3)
L <- 0:60
lfd_sim <- FLQuants(Sim = FLQuant(0, dimnames = list(len = L, year = 2001:2006)))
for (y in 2001:2006) {
  mu <- 50 * (1 - exp(-0.25 * (0:12 + 0.35)))                   # mean length at age
  s  <- sqrt((0.08 * mu)^2 + 1)
  w  <- exp(-0.4 * (0:12)) * exp(rnorm(13, 0, 0.5))             # year classes
  p  <- sapply(seq_along(mu), function(j) pnorm(L + 1, mu[j], s[j]) - pnorm(L, mu[j], s[j])) %*% w
  lfd_sim$Sim[, as.character(y)] <- rmultinom(1, 2000, p)
}
gph_sim <- growth_phase_flicc(FLPar(linf = 50, k = 0.25, M = 0.3, L50 = 30), lfd = lfd_sim,
                              CVL = 0.08, sd0 = 1)
gph_sim$best
plot_growth_lfd(gph_sim)



## ----child = "chapters/07_selectivity.Rmd"------------------------------------

## ----sel----------------------------------------------------------------------
plot_sel(fit)
selpars_flicc(fit)


## ----sel-fixed, eval = FALSE--------------------------------------------------
## fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
##        settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"),
##        sel_fixed = list(Gillnet = list(SL50 = 30, SL95 = 34)))   # logistic
## # dsnormal: list(mode = 24, lsd = 3, rsd = 1000)  (rsd = 1000: flat right limb)
## # normal:   list(mode = 24, sd = 4)
## # binormal:   list(mode1 = 15, sd1 = 2, mode2 = 30, sd2 = 5, h = 0.5)
## # bidsnormal: list(mode1 = 15, sd1 = 2, mode2 = 30, lsd2 = 4, rsd2 = 1000, h = 0.5)


## ----selfun, fig.height = 4.5-------------------------------------------------
gn <- lfd[["Gillnet"]]
Lmid_g <- as.numeric(dimnames(gn)$len) + 0.5
Lmean_g <- round(sum(Lmid_g * apply(gn, 1, sum)) / sum(gn), 1)
Lmean_g

ref_L <- c(L50mat = 31.2, Lmean = Lmean_g, Linf = 55.7)
plot_selfun("dsnormal", list(mode = 34, lsd = 4, rsd = c(1000, round(c(0.85, 0.7, 0.55, 0.4) * 55.7, 1))),
            lfd = gn, ref = ref_L)
plot_selfun("logistic", list(SL50 = c(28, 31, 34), SL95 = 37), lfd = gn, ref = ref_L[c(1,3)])


## ----combine-lfd, fig.height = 3.5--------------------------------------------
yr_last  <- tail(as.numeric(dimnames(lfd[[1]])$year), 1)
lfd_fleet <- combine_lfd(window(lfd, start = yr_last), catch_by_gear = c(0.7,0.3) ,
                         lhpar = lhpar, ess.n = 500)
plot_lfd(window(lfd, start = yr_last))
plot_lfd(lfd_fleet)


## ----selfun-twomodes, fig.height = 4.5----------------------------------------
plot_selfun("binormal", list(mode1 = 22, sd1 = 3.3, mode2 = 31, sd2 = c(3, 8), h = 0.8),
            lfd = lfd_fleet, ref = ref_L)
plot_selfun("bidsnormal", list(mode1 = 22, sd1 = 3.3, mode2 = 31, lsd2 = 3,
                               rsd2 = c(1000, 10, 4), h = 0.7),
            lfd = lfd_fleet, ref = ref_L)


## ----dome-sens, fig.height = 4.5----------------------------------------------
ds <- dome_sens_flicc(fit, gear = "Gillnet")
ds$table
plot_dome_sens(ds)
plot_dome_sens(ds, "spr", ci = TRUE)


## ----dome-sens-cut------------------------------------------------------------
fit_cut_sel <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                     settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", tail_cut = 47))
ds_cut <- dome_sens_flicc(fit_cut_sel, gear = "Gillnet", se = FALSE)
data.frame(scenario = ds$table$scenario, SPR_no_cut = ds$table$SPR_last,
           SPR_cut47 = ds_cut$table$SPR_last)


## ----dome-ensemble, fig.height = 4--------------------------------------------
ens <- ensemble_flicc(ds$fits[-1])     # the five dome levels, equal weights
ens$members
tail(ens$table[, c("year", "spr", "lwr", "upr", "sd_within", "sd_between")], 3)
plot_ensemble_flicc(ens)


## ----compare-sel, fig.height = 4.5--------------------------------------------
gill_fun <- c(logistic = "logistic", dsnormal_rsd1000 = "dsnormal",
              normal = "normal", dsnormal = "dsnormal")
gill_fix <- list(NULL, list(rsd = 1000), NULL, NULL)
fits_sel <- Map(function(fun, fix) {
  fiticc(lfd, stklen, sel_fun = c("dsnormal", fun), catch_by_gear = cbg,
         sel_fixed = if (is.null(fix)) NULL else list(Gillnet = fix),
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
}, gill_fun, gill_fix)
plot_compare_sel(fits_sel)
sel_tab <- data.frame(fit = names(fits_sel),
                      AIC = sapply(fits_sel, function(f) round(AIC(LLflicc(f)), 1)),
                      SPR_last = sapply(fits_sel, function(f) round(sprcur_flicc(f), 3)),
                      row.names = NULL)
sel_tab$dAIC <- sel_tab$AIC - min(sel_tab$AIC)
sel_tab



## ----child = "chapters/08_fitting.Rmd"----------------------------------------

## ----fit----------------------------------------------------------------------
system.time(fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                   settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating")))
sigmaF_flicc(fit)
neff_flicc(fit)[1:4, ]


## ----check-conv---------------------------------------------------------------
flicc_convergence(fit)


## ----check-len, fig.height = 6------------------------------------------------
plot_len(fit, by_gear = TRUE, year = 2020:2024)
plot_len(fit)


## ----resid, fig.height = 4.5--------------------------------------------------
res_osa <- resid_flicc(fit, type = "osa")
c(mean = mean(res_osa$resid, na.rm = TRUE), sd = sd(res_osa$resid, na.rm = TRUE))
plot_resid_qq(res_osa)


## ----resid-bubble, fig.height = 5---------------------------------------------
plot_resid_bubble(res_osa, fit = fit, cohorts = TRUE)


## ----resid-len-pooled, fig.height = 3.5---------------------------------------
plot_resid_len(res_osa, pool_years = TRUE)


## ----resid-len-years, fig.height = 6------------------------------------------
plot_resid_len(res_osa, year = 2022:2024)
plot_resid_len(fit, type = "pearson", year = 2022:2024)


## ----lmean-resid, fig.height = 3.5--------------------------------------------
lmr <- lmean_resid_flicc(fit)
lmr$runs
lmr$plot


## ----ppc----------------------------------------------------------------------
ppc <- ppc_flicc(fit, nsim = 300)
ppc


## ----ppc-lfd, fig.height = 3.5------------------------------------------------
plot_ppc_lfd(ppc, pool_years = TRUE)


## ----ppc-numbers, include = FALSE---------------------------------------------
pb <- setNames(round(100 * ppc$coverage_bins$share_bins_outside_95),
               as.character(ppc$coverage_bins$gear))
yr4 <- tail(seq_along(dimnames(lfd[[1]])$year), 4)
dy <- attr(plot_ppc_lfd(ppc, year = yr4), "data")
red <- aggregate(outside ~ gear + year, data = dy, FUN = sum)
worst <- red[which.max(red$outside), ]
rl <- range(dy$len[dy$outside & dy$gear == worst$gear & dy$year == worst$year])


## ----ppc-lfd-years, fig.height = 7--------------------------------------------
plot_ppc_lfd(ppc, year = yr4)


## ----ppc-stats, fig.height = 7------------------------------------------------
ppc$plot


## ----resid-numbers, include = FALSE-------------------------------------------
osa_sd <- round(sd(res_osa$resid, na.rm = TRUE), 2)
sdz <- setNames(lmr$runs$sd_z, as.character(lmr$runs$gear))
out_lmean <- with(ppc$coverage[ppc$coverage$stat == "Lmean", ],
                  setNames(round(100 * share_outside_95), as.character(gear)))


## ----check-retro, fig.height = 4----------------------------------------------
yrs <- as.numeric(dimnames(lfd[[1]])$year)
peels <- 1:3

retro <- Map(function(p) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"),
         years = head(yrs, -p), se = FALSE)
}, peels)
names(retro) <- paste0("peel", peels)
plot_spr(c(list(full = fit), retro))

spr_full <- c(fit$report$spr)
spr_peel <- sapply(retro, sprcur_flicc)               # terminal SPR of each peel
spr_same <- spr_full[length(spr_full) - peels]        # same years in the full fit
mohn_rho <- mean((spr_peel - spr_same) / spr_same)
round(mohn_rho, 3)



## ----child = "chapters/09_reference_points.Rmd"-------------------------------

## ----refpts-------------------------------------------------------------------
fspr_flicc(fit)                 # stored Fspr40
fspr_flicc(fit, spr = 30)       # computed
fcur_flicc(fit); sprcur_flicc(fit)


## ----eqcurves, fig.height = 5-------------------------------------------------
eqstk <- eqstklen(fit, s = 0.75)
eqstk@refpts
plot_eqcurves(eqstk)


## ----prbrp--------------------------------------------------------------------
prbrp <- prbrp_flicc(fit)
fmsy_flicc(prbrp, s = 0.75)
fspr_flicc(fit, spr = 35)


## ----steepness, fig.height = 5------------------------------------------------
round(spr2brel(c(0.2, 0.3, 0.4), s = 0.75, sprref = 40), 3)
plot_LBAdvice(FLStocks(SPR = flicc2FLStockR(fit, rel = TRUE),
                       relB_s0.75 = flicc2FLStockR(fit, rel = TRUE, s = 0.75)), panel = 1)


## ----status-plots, fig.height = 4---------------------------------------------
plot_spr(fit, ci = TRUE)
plot_lbfao(fit)


## ----advice, fig.height = 5---------------------------------------------------
stk <- flicc2FLStockR(fit)
stk@refpts
plot_LBAdvice(stk)


## ----advice-rel, fig.height = 4-----------------------------------------------
stk_rel <- flicc2FLStockR(fit, rel = TRUE)
plot_LBAdvice(stk_rel, panels = 1) + ylim(0, 1.5)


## ----contour, fig.height = 8.5, fig.width = 9---------------------------------
cf <- contour_flicc(fit)
cf$current
plot_contour_flicc(cf)


## ----contour-share, fig.height = 4.4, fig.width = 9---------------------------
cs <- contour_flicc(fit, view = "share", gear = "Gillnet")
plot_contour_flicc(cs)


## ----contour-s, fig.height = 8.5, fig.width = 9-------------------------------
plot_contour_flicc(cf, s = 0.75)


## ----contour-numbers, include = FALSE-----------------------------------------
eng <- FLicc:::pr_engine_flicc(fit); 
st <- FLicc:::gear_state_flicc(fit)
sel_g <- function(g, loc = NULL) {
  p <- st$pars[[g]]
  if (!is.null(loc)) p <- FLicc:::shift_sel_flicc(st$type[[g]], p, loc)
  sel_curve_flicc(st$Lmid, st$type[[g]], p)
}
pr_gn <- function(loc, Fg) eng(st$F[["Trawl"]] * sel_g("Trawl") + Fg * sel_g("Gillnet", loc))
gn_loc <- st$pars[["Gillnet"]]$SL50
gn_F <- st$F[["Gillnet"]]
gn_now <- pr_gn(gn_loc, gn_F)
gn_lmat <- pr_gn(cf$Lmat, gn_F)
f40 <- tryCatch(uniroot(function(f) pr_gn(gn_loc, f)[["SPR"]] - 0.4, c(0, gn_F))$root,
                error = function(e) NA)
gn_cut <- if (is.na(f40)) "more than 100" else round(100 * (1 - f40 / gn_F))
gn_ypr_change <- round(100 * (gn_lmat[["YPR"]] / gn_now[["YPR"]] - 1))



## ----child = "chapters/10_single_year.Rmd"------------------------------------

## ----one-year-fits------------------------------------------------------------
yr <- 2024
set_ref <- list(CVL = 0.1, obs_model = "mn", FM_max = 4)

fit1 <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
               settings = set_ref, years = yr)

rsd <- c(dome_0.85 = 47, dome_0.7 = 39, dome_0.55 = 31)       # x Linf = 55.7
dome_runs <- lapply(rsd, function(r)
  fiticc(lfd, stklen, sel_fun = c("dsnormal", "dsnormal"), catch_by_gear = cbg,
         settings = set_ref, years = yr, sel_fixed = list(Gillnet = list(rsd = r))))

cvl <- c(CVL_0.08 = 0.08, CVL_0.12 = 0.12)
cvl_runs <- lapply(cvl, function(cv)
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = modifyList(set_ref, list(CVL = cv)), years = yr))

linf <- c(Linf_50 = 50.1, Linf_61 = 61.3)
linf_runs <- lapply(linf, function(li) {
  lh <- lhpar; lh["linf"] <- li
  fiticc(lfd, stocklen(lfd, lh, m_model = "constant"), sel_fun = sel_fun,
         catch_by_gear = cbg, settings = set_ref, years = yr)
})

runs1 <- c(list(reference = fit1), dome_runs, cvl_runs, linf_runs)
tab1 <- fits_summary(runs1)
tab1


## ----resid-singleyr, fig.height = 3.5-----------------------------------------
plot_resid_len(fit1, type = "pearson")
plot_resid_len(fit1, type = "osa")


## ----one-year-ensemble, fig.height = 4----------------------------------------
ok1 <- runs1[tab1$conv == "ok"]

ens_eq  <- ensemble_flicc(ok1, weights = "equal")
ens_aic <- ensemble_flicc(ok1, weights = "aic")
ens_eq$members
ens_aic$members
rbind(equal = ens_eq$table, aic = ens_aic$table)[, c("year", "spr", "lwr", "upr", "sd_within", "sd_between")]
plot_ensemble_flicc(ens_eq)
plot_ensemble_flicc(ens_aic)


## ----one-fleet-fits-----------------------------------------------------------

fit2gear <-fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
              settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))

lfd_fleet <- combine_lfd(lfd, catch_by_gear = c(0.95,0.05) ,
                         lhpar = lhpar, ess.n = 500)


stk_fleet <- stocklen(lfd_fleet, lhpar, m_model = "constant")




fleet_fun <- c(logistic = "logistic", binormal = "binormal",
               bidsnormal = "bidsnormal", bidsnormal_flat = "bidsnormal")
fleet_fix <- list(NULL, NULL, NULL, list(fleet = list(rsd2 = 1000)))

set_fleet <- modifyList(set_ref, list(dm_theta = "fixed"))

fleet_runs <- Map(function(fun, fix) {
  fiticc(lfd_fleet, stk_fleet, sel_fun = fun, catch_by_gear = 1,
         settings = set_fleet, sel_fixed = fix)
}, fleet_fun, fleet_fix)

tab_fleet <- fits_summary(c(list(ref2gear=fit2gear),fleet_runs))
tab_fleet
plot_compare_sel(fleet_runs)


## ----one-fleet-len, fig.height = 3.5------------------------------------------
ok_fleet   <- tab_fleet$conv == "ok"
best_fleet <- if (any(ok_fleet)) tab_fleet$run[ok_fleet][which.min(tab_fleet$AIC[ok_fleet])] else
  tab_fleet$run[which.min(tab_fleet$AIC)]
best_fleet
plot_len(fleet_runs[[best_fleet]])


## ----one-fleet-pars-----------------------------------------------------------
selpars_flicc(fleet_runs$bidsnormal)
selpars_flicc(fleet_runs$bidsnormal_flat)



## ----child = "chapters/11_indicators.Rmd"-------------------------------------

## ----lbi-empirical, fig.height = 8--------------------------------------------
lbi <- lbi_flicc(lfd, lhpar, pool = TRUE, catch_by_gear = c(Trawl = 0.95, Gillnet = 0.05))
lbi$values[lbi$values$year == 2024, c("gear", "Lc", "Lmean", "LFM", "L95", "Lmax5", "Pmat", "Pmega", "Z")]
lbi$values[lbi$values$year == 2024, c("gear", "Pobj", "cp_selectivity", "cp_test", "cp_status")]
plot_lbi(lbi)


## ----lbi-trend, fig.height = 6------------------------------------------------
plot_lbi(lbi, type = "trend")


## ----lbi-relative, fig.height = 6---------------------------------------------
lbi_ref <- lbi_flicc(lfd, lhpar, ref_years = 2013:2015)
plot_lbi(lbi_ref, type = "relative")


## ----lbispr-------------------------------------------------------------------
lbi <- LBIspr(fit, spr = 40, thresh = 0.75)
attr(lbi, "Lref")
plot_LBIspr(fit, thresh = 0.75)
plot_LBIp(fit, year = 2024)


## ----lbimean------------------------------------------------------------------
LBImean(fit, ref = "FM")



## ----child = "chapters/12_sample_size.Rmd"------------------------------------

## ----ess-demo-----------------------------------------------------------------
sel_fun <- c("dsnormal", "logistic")
cbg <- c(0.95, 0.05)
stklen <- stocklen(lfd, lhpar, m_model = "constant")

ess_values <- c(50, 100, 200, 500, 1000)

# multinomial: one fit per ESS
fits_mn <- Map(function(ess) {
  fiticc(lfdess(lfd, ess.g = ess), stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "mn"))
}, ess_values)
names(fits_mn) <- paste0("mn_", ess_values)

# Dirichlet-multinomial with fixed precision: one fit per ESS
fits_dm <- Map(function(ess) {
  fiticc(lfdess(lfd, ess.g = ess), stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "fixed"))
}, ess_values)
names(fits_dm) <- paste0("dm_", ess_values)

# raw counts (no rescaling)
fit_raw_mn <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                     settings = list(CVL = 0.1, obs_model = "mn"))
fit_raw_dm <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                     settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "fixed"))
raw_n <- round(mean(quantSums(lfd[["Trawl"]])))   # fish per year, trawl

fits_mn <- c(fits_mn, raw = list(fit_raw_mn))
fits_dm <- c(fits_dm, raw = list(fit_raw_dm))

# standard error of log SPR in the last year
se_logspr <- function(fit) tail(spr_ci_flicc(fit)$se_log, 1)
ess_table <- data.frame(
  ESS        = c(ess_values, paste0("raw (~", raw_n, ")")),
  SPR_mn     = round(sapply(fits_mn, sprcur_flicc), 3),
  SPR_dm     = round(sapply(fits_dm, sprcur_flicc), 3),
  se_mn      = round(sapply(fits_mn, se_logspr), 3),
  se_dm      = round(sapply(fits_dm, se_logspr), 3),
  maxgrad_dm = signif(sapply(fits_dm, function(fit) fit$max_gradient), 2),
  row.names  = NULL)
ess_table


## ----ess-plot, fig.height = 4-------------------------------------------------
plot_spr(fits_dm[c("dm_50", "dm_200", "dm_1000", "raw")], ci = TRUE)


## ----dm-sat-------------------------------------------------------------------
dm_inputs <- list(ESS_500 = lfdess(lfd, ess.g = 500), ESS_1000 = lfdess(lfd, ess.g = 1000), raw = lfd)
fits_dms <- Map(function(d) {
  fiticc(d, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
}, dm_inputs)

neff_gear <- function(fit, g) { ne <- neff_flicc(fit); round(mean(ne$neff[ne$gear == g])) }
dm_table <- data.frame(
  input        = names(fits_dms),
  SPR_last     = round(sapply(fits_dms, sprcur_flicc), 3),
  se_log       = round(sapply(fits_dms, se_logspr), 3),
  neff_Trawl   = sapply(fits_dms, neff_gear, g = "Trawl"),
  neff_Gillnet = sapply(fits_dms, neff_gear, g = "Gillnet"),
  row.names = NULL)
dm_table
plot_spr(c(fixed_ESS_200 = list(fits_dm$dm_200), fits_dms), ci = TRUE)



## ----child = "chapters/13_fishing_mortality.Rmd"------------------------------

## ----f-options----------------------------------------------------------------
f_settings <- list(
  fixed_0.3 = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", sigmaF = 0.3),                       # default
  prior     = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", prior_sigmaF = c(log(0.5), 0.3, 1)), # penalised
  F_re      = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", F_re = TRUE),                        # random effects
  no_RW     = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", prior_sigmaF = FALSE))               # no random walk

fits_F <- Map(function(set) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg, settings = set)
}, f_settings)

data.frame(option   = names(fits_F),
           sigmaF   = round(sapply(fits_F, sigmaF_flicc), 3),
           SPR_last = round(sapply(fits_F, sprcur_flicc), 3),
           se_log   = round(sapply(fits_F, se_logspr), 3),
           row.names = NULL)


## ----f-options-plot, fig.height = 4-------------------------------------------
plot_spr(fits_F, ci = TRUE)

F_by_option <- FLQuants(fixed_0.3 = fits_F$fixed_0.3$report$Fap,
                        prior     = fits_F$prior$report$Fap,
                        F_re      = fits_F$F_re$report$Fap,
                        no_RW     = fits_F$no_RW$report$Fap)
ggplot(F_by_option, aes(year, data, colour = qname)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  labs(x = "Year", y = "Apical F", colour = NULL) + theme_bw()


## ----sigmaF-range, fig.height = 4---------------------------------------------
sigmaF_values <- c(0.05, 0.1, 0.3, 0.6, 1)

fits_sigmaF <- Map(function(sF) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", sigmaF = sF), se = FALSE)
}, sigmaF_values)
names(fits_sigmaF) <- paste0("sigmaF_", sigmaF_values)

plot_spr(fits_sigmaF)
data.frame(sigmaF   = sigmaF_values,
           SPR_last = round(sapply(fits_sigmaF, sprcur_flicc), 3),
           neg_logL = round(sapply(fits_sigmaF, function(fit) fit$opt$objective), 1))


## ----f-re-time----------------------------------------------------------------
system.time(fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg, settings = f_settings$fixed_0.3))
system.time(fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg, settings = f_settings$F_re))


## ----by-year, fig.height = 4--------------------------------------------------
fit_by_year <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                      settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"), by_year = TRUE)

plot_spr(list(joint = fit, by_year = fit_by_year))
plot_sel(flicc_stklen(fit_by_year))



## ----child = "chapters/14_growth_variability.Rmd"-----------------------------

## ----cvl-sens, fig.height = 4-------------------------------------------------
cvl_values <- c(0.07, 0.10, 0.13)

fits_cvl <- Map(function(cvl) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = cvl, obs_model = "dm", dm_theta = "saturating"), se = FALSE)
}, cvl_values)
names(fits_cvl) <- paste0("CVL_", cvl_values)

fits_cvl_cut <- Map(function(cvl) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = cvl, obs_model = "dm", dm_theta = "saturating", tail_cut = 47), se = FALSE)
}, cvl_values)
names(fits_cvl_cut) <- paste0("CVL_", cvl_values, "_cut47")

cvl_table <- data.frame(CVL       = cvl_values,
                        SPR_no_cut = round(sapply(fits_cvl, sprcur_flicc), 3),
                        SPR_cut47  = round(sapply(fits_cvl_cut, sprcur_flicc), 3),
                        row.names = NULL)
cvl_table
plot_spr(c(fits_cvl, fits_cvl_cut))


## ----maxsd-cv-----------------------------------------------------------------
z <- seq(-2, 2, length.out = 13)          # maxsd = 2, ngtg = 13 (LBSPR)
w <- dnorm(z) / sum(dnorm(z))
sqrt(sum(w * z^2))                        # realised CV factor

z <- seq(-3, 3, length.out = 13)          # maxsd = 3, ngtg = 13 (default)
w <- dnorm(z) / sum(dnorm(z))
sqrt(sum(w * z^2))

55.7 * (1 + c(2, 3) * 0.1)                # largest group Linf for maxsd 2 and 3


## ----maxsd-sens---------------------------------------------------------------
maxsd_values <- c(2, 3, 3)
ngtg_values  <- c(13, 13, 19)

fits_maxsd <- Map(function(maxsd, ngtg) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", maxsd = maxsd, ngtg = ngtg), se = FALSE)
}, maxsd_values, ngtg_values)

fits_maxsd_cut <- Map(function(maxsd, ngtg) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
         settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", maxsd = maxsd, ngtg = ngtg,
                         tail_cut = 47), se = FALSE)
}, maxsd_values, ngtg_values)

data.frame(maxsd      = maxsd_values,
           ngtg       = ngtg_values,
           SPR_no_cut = round(sapply(fits_maxsd, sprcur_flicc), 3),
           SPR_cut47  = round(sapply(fits_maxsd_cut, sprcur_flicc), 3))


## ----priors, fig.height = 4---------------------------------------------------
fit_pri <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                  settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating",
                                  linf.sd = 0.05, Mk.sd = 0.1, CVL.sd = 0.25))

est <- fit_pri$par[fit_pri$par$par %in% c("log_Linf", "log_Mk", "log_Galpha"), ]
est
exp(est$mpd[est$par == "log_Linf"])            # Linf (input 55.7)
exp(est$mpd[est$par == "log_Mk"])              # M/k  (input 0.162 / 0.08 = 2.03)
1 / sqrt(exp(est$mpd[est$par == "log_Galpha"])) # CVL  (input 0.1)

plot_spr(list(fixed = fit, priors = fit_pri), ci = TRUE)


## ----priors-numbers, include = FALSE------------------------------------------
pp <- function(nm, col = "mpd") est[[col]][est$par == nm][1]



## ----child = "chapters/14b_sex_structure.Rmd"---------------------------------

## ----sex-data-----------------------------------------------------------------
data(dimorph)
lhpar_dimorph
lapply(lfd_dimorph, function(q) dimnames(q)$unit)


## ----sex-long, eval = FALSE---------------------------------------------------
## lfd <- FLQuantLen(lfd_long_to_wide(df_long_with_sex), unit = "cm")


## ----sex-lfd, fig.height = 6--------------------------------------------------
plot_lfd(FLQuants(lapply(lfd_dimorph, function(q) q[, ac(2021:2025)])), type = "prop")


## ----sex-lhpar----------------------------------------------------------------
lh_F <- FLPar(lhpar_dimorph@.Data[, "F", 1], params = dimnames(lhpar_dimorph)$params)
lhpar_sex(lh_F, linf = 45, k = 0.25, L50 = 28)


## ----sex-stock----------------------------------------------------------------
stk_sex <- stocklen_sex(lfd_dimorph, lhpar_dimorph, m_model = "constant")
plot_growth(stk_sex, sex = "M")


## ----sex-fits-----------------------------------------------------------------
sf2 <- c("dsnormal", "logistic")
set_s <- list(obs_model = "mn", surveys = "Survey")

trawl_U <- FLQuants(Trawl = lfd_dimorph$Trawl[, , "U"])
fit_1sex <- fiticc(trawl_U, stocklen(trawl_U, lh_F, m_model = "constant"), sel_fun = "dsnormal",
                   catch_by_gear = 1, settings = list(obs_model = "mn"))
fit_2sex_U <- fiticc(FLQuants(Trawl = lfd_dimorph$Trawl), stocklen_sex(trawl_U, lhpar_dimorph, m_model = "constant"),
                     sel_fun = "dsnormal", catch_by_gear = 1, settings = list(obs_model = "mn"))
fit_sex <- fiticc(lfd_dimorph, stk_sex, sel_fun = sf2, catch_by_gear = 1, settings = set_s)
fit_sex_est <- fiticc(lfd_dimorph, stocklen_sex(lfd_dimorph, lh_F, male = list(linf = 52, k = 0.25, L50 = 28),
                                              m_model = "constant"),
                      sel_fun = sf2, catch_by_gear = 1, settings = modifyList(set_s, list(linf_m.sd = 0.3)))

runs <- list(one_sex = fit_1sex, two_sex_trawl = fit_2sex_U, two_sex_survey = fit_sex,
             male_Linf_est = fit_sex_est)
tab <- data.frame(year = truth_dimorph$year, true = round(truth_dimorph$spr, 3),
                  sapply(runs, function(f) round(as.numeric(f$report$spr), 3)))
tail(tab, 5)
fits_summary(runs)[, c("run", "conv", "SPR", "lwr", "upr", "Linf_M", "Mk_M")]


## ----sex-ratio, fig.height = 3.5----------------------------------------------
plot_sexratio(fit_sex)


## ----sex-len, fig.height = 6--------------------------------------------------
plot_len(fit_sex, year = 2021:2025)
plot_resid_bubble(fit_sex, type = "osa")


## ----sex-eq, fig.height = 5---------------------------------------------------
eq_sex <- eqstklen(fit_sex, s = 0.75)
attr(eq_sex, "refpts")
plot_eqcurves(eq_sex)


## ----sex-lbi, fig.height = 5--------------------------------------------------
plot_lbi(lbi_flicc(fit_sex))



## ----child = "chapters/14c_dimorph_example.Rmd"-------------------------------

## ----dim-data-----------------------------------------------------------------
data(dimorph)
lhpar_dimorph


## ----dim-stock, fig.height = 3.2----------------------------------------------
stk_dim <- stocklen_sex(lfd_dimorph, lhpar_dimorph, m_model = "constant")
plot_mat(stk_dim)
plot_lw(stk_dim)
plot_m(stk_dim)


## ----dim-growth, fig.height = 4-----------------------------------------------
plot_growth(stk_dim, sex = "F") + ggplot2::ggtitle("Females")
plot_growth(stk_dim, sex = "M") + ggplot2::ggtitle("Males")


## ----dim-lfd, fig.height = 6--------------------------------------------------
plot_lfd(FLQuants(lapply(lfd_dimorph, function(q) q[, ac(2019:2025)])), type = "prop")


## ----dim-fit------------------------------------------------------------------
fit_dim <- fiticc(lfd_dimorph, stk_dim, sel_fun = c("dsnormal", "logistic"), catch_by_gear = 1,
                  settings = list(obs_model = "mn", surveys = "Survey"))
fit_dim$opt$convergence


## ----dim-len, fig.height = 7--------------------------------------------------
plot_len(fit_dim, year = 2019:2025)


## ----dim-sexratio, fig.height = 3.5-------------------------------------------
plot_sexratio(fit_dim)


## ----dim-resid, fig.height = 5------------------------------------------------
plot_resid_bubble(fit_dim, type = "osa")
plot_resid_qq(fit_dim)
plot_resid_len(fit_dim, type = "osa", pool_years = TRUE)


## ----dim-ppc, fig.height = 6--------------------------------------------------
ppc_dim <- ppc_flicc(fit_dim, nsim = 300)
plot_ppc_lfd(ppc_dim, year = 2023:2025)


## ----dim-n-sex, fig.height = 3.5----------------------------------------------
nsx <- as.data.frame(fit_dim$report$N_sex[, "2025"])
ggplot2::ggplot(nsx, ggplot2::aes(as.numeric(as.character(len)), data, colour = unit)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::scale_colour_manual(values = c(F = "#D7191C", M = "#2C7BB6"), name = "Sex") +
  ggplot2::labs(x = "Length", y = "Numbers per recruit, 2025") + ggplot2::theme_bw()


## ----dim-phase, fig.height = 3.5----------------------------------------------
plot_growth_phase(growth_phase_flicc(fit_dim, sex = "F", gear = "Survey"))
plot_growth_phase(growth_phase_flicc(fit_dim, sex = "M", gear = "Survey"))


## ----dim-status, fig.height = 4-----------------------------------------------
plot_spr(fit_dim, ci = TRUE, ci_type = "ribbon")
eq_dim <- eqstklen(fit_dim, s = 0.75)
attr(eq_dim, "refpts")
plot_eqcurves(eq_dim)


## ----dim-lbi, fig.height = 5--------------------------------------------------
plot_lbi(lbi_flicc(fit_dim, sex = "each"))


## ----dim-pool-data------------------------------------------------------------
pool <- function(q) {
  u <- FLCore::unitSums(q)
  dimnames(u)$unit <- "U"
  u
}
lfd_pool <- FLQuants(lapply(lfd_dimorph, pool))
lh_F <- FLPar(lhpar_dimorph@.Data[, "F", 1], params = dimnames(lhpar_dimorph)$params)
lfd_pool1 <- FLQuants(lapply(lfd_pool, function(q) { dimnames(q)$unit <- "unique"; q }))
fit_1sex <- fiticc(lfd_pool1, stocklen(lfd_pool1, lh_F, m_model = "constant"),
                   sel_fun = c("dsnormal", "logistic"), catch_by_gear = 1,
                   settings = list(obs_model = "mn", surveys = "Survey"))


## ----dim-1sex-len, fig.height = 6---------------------------------------------
plot_len(fit_1sex, by_gear = TRUE, year = 2019:2025)


## ----dim-2sex-pool------------------------------------------------------------
fit_2pool <- fiticc(lfd_pool, stocklen_sex(lfd_pool, lhpar_dimorph, m_model = "constant"),
                    sel_fun = c("dsnormal", "logistic"), catch_by_gear = 1,
                    settings = list(obs_model = "mn", surveys = "Survey"))


## ----dim-2pool-len, fig.height = 6--------------------------------------------
plot_len(fit_2pool, year = 2019:2025)


## ----dim-compare, fig.height = 4.5--------------------------------------------
runs <- list("Ignoring sex (female LH, pooled)" = fit_1sex,
             "Sex-structured, pooled data" = fit_2pool,
             "Sex-structured, sexed survey" = fit_dim)
plot_spr(runs, ci = TRUE, ci_type = "ribbon") +
  ggplot2::geom_line(data = data.frame(year = truth_dimorph$year, spr = truth_dimorph$spr),
                     ggplot2::aes(x = year, y = spr), inherit.aes = FALSE, linetype = 2, linewidth = 0.8)


## ----dim-compare-tab----------------------------------------------------------
err <- sapply(runs, function(f) c(f$report$spr) / truth_dimorph$spr - 1)
round(rbind("2011-2015" = colMeans(err[1:5, ]), "2016-2020" = colMeans(err[6:10, ]),
            "2021-2025" = colMeans(err[11:15, ])), 2)


## ----dim-compare-sel, fig.height = 4------------------------------------------
L <- seq(1, 79, by = 0.5)
true_sel <- rbind(data.frame(L = L, sel = sel_curve_flicc(L, "dsnormal", list(mode = 34, lsd = 6, rsd = 20)),
                             gear = "Trawl"),
                  data.frame(L = L, sel = sel_curve_flicc(L, "logistic", list(SL50 = 12, SL95 = 16)),
                             gear = "Survey"))
plot_compare_sel(runs, colours = c("#D55E00", "#E69F00", "#0072B2")) +
  ggplot2::geom_line(data = true_sel, ggplot2::aes(L, sel), inherit.aes = FALSE,
                     linetype = 2, linewidth = 0.8)



## ----child = "chapters/15_tail.Rmd"-------------------------------------------

## ----tailcut, fig.height = 4--------------------------------------------------
fit_cut <- fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                  settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating", tail_cut = 47))  # ~0.85 x Linf
fit_cut$tmb_data$LLB[fit_cut$tmb_data$cut_bin + 1]
plot_spr(list(no_cut = fit, cut_47 = fit_cut), ci = TRUE)



## ----child = "chapters/16_model_equations.Rmd"--------------------------------

## ----model-compare, fig.height = 4--------------------------------------------
model_settings <- list(
  gtg_mn   = list(CVL = 0.1, pop_model = "gtg",   obs_model = "mn"),
  gtg_dm   = list(CVL = 0.1, pop_model = "gtg",   obs_model = "dm", dm_theta = "saturating"),
  gtg_nb   = list(CVL = 0.1, pop_model = "gtg",   obs_model = "nb"),
  gamma_nb = list(CVL = 0.1, pop_model = "gamma", obs_model = "nb"))   # as fishblicc

fits_model <- Map(function(set) {
  fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg, settings = set)
}, model_settings)

data.frame(model    = names(fits_model),
           SPR_last = round(sapply(fits_model, sprcur_flicc), 3),
           se_log   = round(sapply(fits_model, se_logspr), 3),
           row.names = NULL)
plot_spr(fits_model, ci = TRUE)


## ----model-time---------------------------------------------------------------
system.time(fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                   settings = model_settings$gtg_nb))
system.time(fiticc(lfd, stklen, sel_fun = sel_fun, catch_by_gear = cbg,
                   settings = model_settings$gamma_nb))



## ----child = "chapters/17_reproduce.Rmd"--------------------------------------

## ----lbspr-data---------------------------------------------------------------
data("lbspr_test")   # lfd_lbspr, lhpar_lbspr, lbspr_est
lfd_lbspr
lhpar_lbspr          # k only scales F and M; F/M and SPR depend on M/k = 1.5
stk_lbspr <- stocklen(lfd_lbspr, lhpar_lbspr, m_model = "constant")


## ----lbspr-fits---------------------------------------------------------------
set_lbspr <- list(CVL = 0.1, ngtg = 13, maxsd = 2, obs_model = "mn")
yrs_lbspr <- dimnames(stk_lbspr)$year
fits_lbspr_y <- Map(function(y) {
  lfd_y <- FLQuants(lapply(lfd_lbspr, function(x) x[, y]))
  fiticc(lfd_y, stocklen(lfd_y, lhpar_lbspr, m_model = "constant"), sel_fun = "logistic",
         catch_by_gear = 1, settings = set_lbspr, verbose = FALSE)
}, yrs_lbspr)
fit_lbspr_j <- fiticc(lfd_lbspr, stk_lbspr, sel_fun = "logistic", catch_by_gear = 1,
                      settings = set_lbspr)
fit_lbspr_d <- fiticc(lfd_lbspr, stk_lbspr, sel_fun = "logistic", catch_by_gear = 1,
                      settings = list(CVL = 0.1))


## ----lbspr-compare------------------------------------------------------------
lbspr_out <- lbspr_est[, c("year", "SL50", "SL95", "FM", "SPR")]
lbspr_smooth <- setNames(lbspr_est[, c("year", "SL50_smooth", "SL95_smooth", "FM_smooth", "SPR_smooth")],
                         c("year", "SL50", "SL95", "FM", "SPR"))
est_flicc <- function(f) {
  sm <- exp(f$obj$env$parList(par = f$obj$env$last.par.best)$Sm) * 100   # relative to Linf
  data.frame(year = as.numeric(f$tmb_data$year_names), SL50 = sm[1], SL95 = sm[1] + sm[2],
             FM = as.numeric(f$report$Fap) / 0.3, SPR = as.numeric(f$report$spr))
}
lbspr_runs <- list("LBSPR, by year" = lbspr_out,
                   "LBSPR, smoothed" = lbspr_smooth,
                   "FLicc mn, by year" = do.call(rbind, lapply(fits_lbspr_y, est_flicc)),
                   "FLicc mn, joint" = est_flicc(fit_lbspr_j),
                   "FLicc default" = est_flicc(fit_lbspr_d))
lbspr_tab <- do.call(rbind, Map(function(x, nm) cbind(run = nm, round(x, 3)), lbspr_runs, names(lbspr_runs)))
rownames(lbspr_tab) <- NULL
lbspr_tab


## ----lbspr-plot, fig.height = 3.5---------------------------------------------
lbspr_long <- do.call(rbind, lapply(c("SL50", "FM", "SPR"), function(v)
  data.frame(lbspr_tab[, c("run", "year")], what = v, value = lbspr_tab[[v]])))
lbspr_long$run <- factor(lbspr_long$run, levels = names(lbspr_runs))
ggplot(lbspr_long, aes(year, value, colour = run)) +
  geom_line() + geom_point() + facet_wrap(~what, scales = "free_y") +
  scale_colour_manual(values = c("black", "grey55", "#E69F00", "#0072B2", "#009E73")) +
  labs(x = "Year", y = NULL, colour = NULL) + theme_bw() + theme(legend.position = "bottom")
plot_spr(fit_lbspr_d, ci = TRUE)


## ----lbspr-sprgtg-------------------------------------------------------------
spr_y <- function(f) {
  data.frame(year = as.numeric(f$tmb_data$year_names),
             FLicc = spr_gtg_flicc(f, mat = "length")$spr,
             mat_rel = spr_gtg_flicc(f)$spr_gtg,
             mat_rel_lag = spr_gtg_flicc(f, sel_lag = 1)$spr_gtg)
}
spr_cmp <- do.call(rbind, lapply(fits_lbspr_y, spr_y))
cbind(round(spr_cmp, 3), LBSPR = round(lbspr_out$SPR, 3))
spr_gtg_flicc(fit_lbspr_d, sel_lag = 1)   # joint default fit, LBSPR-type SPR with interval


## ----lbb-data-----------------------------------------------------------------
url <- "https://raw.githubusercontent.com/SISTA16/LBB/master/Example_Dat.csv"
lbb_dat <- tryCatch(read.csv(url), error = function(e) NULL)
have_lbb <- !is.null(lbb_dat)
if (have_lbb) {
  d <- subset(lbb_dat, Stock == "Ille_coi_AD" & !is.na(CatchNo))
  d <- aggregate(CatchNo ~ Year + Length, data = d, FUN = sum)
  d$data <- d$CatchNo / ave(d$CatchNo, d$Year, FUN = sum) * 1000
  lfd_illex <- FLQuantLen(lfd_long_to_wide(data.frame(len = d$Length / 10, year = d$Year,
                                                      data = d$data, qname = "Illex")),
                          unit = "cm", midL = FALSE)
}
yrs_illex <- c(1996, 1999, 2001, 2002, 2009, 2014, 2016)


## ----lbb-fits, eval = have_lbb------------------------------------------------
fit_lbb0 <- lbb_flicc(lfd_illex, Linf = 23, Lm50 = 16, years = yrs_illex,
                      pile = 0, ess = 20000)
fit_lbb1 <- lbb_flicc(lfd_illex, Linf = 23, Lm50 = 16, years = yrs_illex,
                      pile = 1, ess = 20000)
fit_lbbj <- lbb_flicc(lfd_illex, Linf = 23, Lm50 = 16, years = yrs_illex,
                      by_year = FALSE)
fit_lbb0


## ----lbb-compare, eval = have_lbb, fig.height = 4.5---------------------------
lbb_out <- data.frame(year = yrs_illex,
                      ZK  = c(2.32, 3.12, 2.57, 3.72, 3.15, 5.24, 6.85),
                      FM  = c(0.265, 0.878, 0.732, 1.467, 0.910, 2.713, 3.570),
                      BB0 = c(0.699, 0.267, 0.293, 0.138, 0.294, 0.065, 0.044))
runs <- list("LBB (JAGS)" = lbb_out, "FLicc, Pile = 0" = fit_lbb0$table,
             "FLicc, Pile = 1" = fit_lbb1$table, "FLicc, Pile = 1, joint" = fit_lbbj$table)
cmp <- do.call(rbind, Map(function(x, nm) data.frame(run = nm, year = x$year, ZK = x$ZK,
                                                     FM = x$FM, BB0 = x$BB0), runs, names(runs)))
cmp$run <- factor(cmp$run, levels = names(runs))
long <- rbind(data.frame(cmp[, 1:2], what = "Z/K", value = cmp$ZK),
              data.frame(cmp[, 1:2], what = "F/M", value = cmp$FM),
              data.frame(cmp[, 1:2], what = "B/B0", value = cmp$BB0))
ggplot(long, aes(year, value, colour = run)) +
  geom_line() + geom_point() +
  facet_wrap(~what, scales = "free_y") +
  scale_colour_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  labs(x = "Year", y = NULL, colour = NULL) +
  theme_bw() + theme(legend.position = "bottom")
dif <- function(x) round(mean(x$ZK - lbb_out$ZK), 2)


## ----lbb-epr-spr, eval = have_lbb, fig.height = 3.5---------------------------
yr_flq <- function(tab, col) FLQuant(matrix(tab[[col]], nrow = 1), dimnames = list(year = tab$year))
bfm <- attr(fit_lbb1$table, "BFM1B0")
plot_epr(FLQuants("LBB (JAGS)" = yr_flq(lbb_out, "BB0"),
                  "FLicc, Pile = 0" = yr_flq(fit_lbb0$table, "EPR"),
                  "FLicc, Pile = 1" = yr_flq(fit_lbb1$table, "EPR"),
                  "FLicc, Pile = 1, joint" = yr_flq(fit_lbbj$table, "EPR")),
         ref_levels = c(bfm, bfm / 2), title = "B/B0 (LBB) and EPR (FLicc)")
plot_spr(FLQuants("FLicc, Pile = 0" = yr_flq(fit_lbb0$table, "SPR"),
                  "FLicc, Pile = 1" = yr_flq(fit_lbb1$table, "SPR"),
                  "FLicc, Pile = 1, joint" = yr_flq(fit_lbbj$table, "SPR")))


## ----lbb-plot, eval = have_lbb, fig.height = 7--------------------------------
plot_lbb(fit_lbb1)


## ----fishblicc-data-----------------------------------------------------------
data("fishblicc_example")
lhpar_fishblicc <- rbind(lhpar_fishblicc,
                         FLPar(M = lhpar_fishblicc["k"] * lhpar_fishblicc["Mk"]))
stklen_fishblicc <- stocklen(lfd_fishblicc, lhpar_fishblicc,
                             m_model = "inverse", reflen = 21.9)


## ----fishblicc-fit------------------------------------------------------------
fit_fishblicc <- fiticc(lfd_fishblicc, stklen_fishblicc,
                        sel_fun = c("dsnormal", "dsnormal", "dsnormal"),
                        catch_by_gear = c(0.1802070, 0.2101353, 0.6096577),
                        settings = list(pop_model = "gamma", obs_model = "nb",
                                        CVL = 0.14, CVL.sd = 0.1,
                                        linf.sd = 2 / 40, Mk.sd = 0.1,
                                        catch_units = "numbers"))
flicc_convergence(fit_fishblicc)
sprcur_flicc(fit_fishblicc)


## ----fishblicc-plots, fig.height = 4.5----------------------------------------
plot_sel(fit_fishblicc)
plot_len(fit_fishblicc, by_gear = TRUE)


## ----fishblicc-eq, fig.height = 5---------------------------------------------
eqstk_fishblicc <- eqstklen(fit_fishblicc, s = 0.7)
eqstk_fishblicc@refpts
plot_eqcurves(eqstk_fishblicc)



## ----child = "chapters/18_mse.Rmd"--------------------------------------------

## ----mse, eval = FALSE--------------------------------------------------------
## est <- mseCtrl(method = flicc.sa, args = list(
##   lhpar = lhpar, sel_fun = sel_fun, ess.g = NULL,   # raw counts
##   settings = list(obs_model = "dm", dm_theta = "saturating", CVL = 0.1,
##                   sigmaF = 0.3,      # fixed (or the F_re estimate from the real data)
##                   tail_cut = 47,     # absolute length, the same in every iteration
##                   spr_ref = 40),
##   spr = 40))



## ----child = "chapters/19_summary.Rmd"----------------------------------------




## ----child = "chapters/20_references.Rmd"-------------------------------------



