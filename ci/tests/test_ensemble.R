# ensemble_flicc: mixture of fits
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stklen <- stocklen(lfd, lhpar, m_model = "constant")
sf <- c("dsnormal", "logistic"); cbg <- c(0.7, 0.3)
set0 <- list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating")
fit <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = set0)
ds <- dome_sens_flicc(fit, gear = "Gillnet")
print(ds$table)
ens <- ensemble_flicc(ds$fits[-1])
print(ens$members); print(tail(ens$table, 3))
ens_aic <- ensemble_flicc(ds$fits[-1], weights = "aic", s = 0.75)
print(ens_aic$members); print(tail(ens_aic$table[, c("year", "spr", "lwr", "upr", "brel", "brel_lwr", "brel_upr")], 3))
# sanity 1: identical members reproduce the single-fit interval
e1 <- ensemble_flicc(list(a = fit, b = fit))
ci <- spr_ci_flicc(fit)
stopifnot(max(abs(e1$table$lwr - ci$lwr)) < 1e-6, max(abs(e1$table$upr - ci$upr)) < 1e-6,
          max(e1$table$sd_between) < 1e-10)
# sanity 2: ensemble interval at least as wide as the narrowest member and contains the median
stopifnot(all(ens$table$lwr <= ens$table$spr), all(ens$table$spr <= ens$table$upr))
# sanity 3: AIC weights refused for fits with different data (tail cut)
fcut <- fiticc(lfd, stklen, sel_fun = sf, catch_by_gear = cbg, settings = c(set0, list(tail_cut = 47)))
r <- tryCatch(ensemble_flicc(list(a = fit, b = fcut), weights = "aic"), error = function(e) "refused")
stopifnot(identical(r, "refused"))
print(ensemble_flicc(list(no_cut = fit, cut47 = fcut))$members)
od <- if (dir.exists("out/tests")) "out/tests" else "."
ggplot2::ggsave(file.path(od, "ensemble_dome.png"), plot_ensemble_flicc(ens), width = 8, height = 4, dpi = 100)
ggplot2::ggsave(file.path(od, "ensemble_dome_brel.png"), plot_ensemble_flicc(ens_aic, "brel"), width = 8, height = 4, dpi = 100)

# few years: error-bar layout
yrs <- dimnames(lfd[[1]])$year
y1 <- yrs[which.max(sapply(yrs, function(y) min(sapply(lfd, function(x) sum(x[, y])))))]
lfd1 <- FLQuants(lapply(lfd, function(x) x[, y1]))
st1 <- stocklen(lfd1, lhpar, m_model = "constant")
f1 <- fiticc(lfd1, st1, sel_fun = sf, catch_by_gear = cbg, settings = set0)
ds1 <- dome_sens_flicc(f1, gear = "Gillnet", rsd = c(1000, 40))
e1 <- ensemble_flicc(ds1$fits)
print(e1$table)
ggplot2::ggsave(file.path(od, "ensemble_1yr.png"), plot_ensemble_flicc(e1), width = 6, height = 4, dpi = 100)
ggplot2::ggsave(file.path(od, "ensemble_1yr_brel.png"),
                plot_ensemble_flicc(ensemble_flicc(ds1$fits, s = 0.75), "brel"), width = 6, height = 4, dpi = 100)

# plot_spr with ci: error bars at the same dodge positions as their points
fl <- list(zz = ds$fits[[2]], aa = ds$fits[[3]], mm = ds$fits[[1]])   # list order != alphabetical
b <- ggplot2::ggplot_build(plot_spr(fl, ci = TRUE))
lay <- vapply(b$plot$layers, function(l) class(l$geom)[1], character(1))
eb <- b$data[[which(lay == "GeomErrorbar")]]; pt <- b$data[[which(lay == "GeomPoint")]]
xe <- tapply(eb$x, eb$colour, mean); xp <- tapply(pt$x, pt$colour, mean)
print(rbind(errorbar = xe, point = xp[names(xe)]))
stopifnot(max(abs(xe - xp[names(xe)])) < 1e-8)
