# lbi_flicc(): empirical length-based indicators
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)[, c("len", "year", "data", "qname")]),
                  unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)

x <- lbi_flicc(lfd, lhpar, catch_by_gear = c(Trawl = 0.95, Gillnet = 0.05), ref_years = 2013:2015)
v <- x$values
stopifnot(inherits(x, "flicc_lbi"), all(levels(v$gear) == c("Trawl", "Gillnet", "pooled")),
          all(v$Pmat >= 0 & v$Pmat <= 1), all(v$Pmega <= v$Popt + v$Pmat + 1e-12),
          all(v$L25 <= v$Lbar), all(v$Lbar <= v$L95), all(v$L95 <= v$Lmax5 + 1),
          all(is.finite(v$Lmax5)), all(v$Lmean >= v$Lbar - 1e-9), all(is.finite(v$Z)),
          all(abs(v$LFM - (0.75 * v$Lc + 0.25 * 55.7)) < 1e-9))
print(v[v$year == 2024, ])

# a known composition: uniform 10-50 cm, 1-cm classes, 100 fish per class
L <- 0:60; n <- ifelse(L >= 10 & L < 50, 100, 0)
u <- FLQuants(U = FLQuant(n, dimnames = list(len = L, year = 2000)))
xu <- lbi_flicc(u, FLPar(linf = 55, L50 = 30))$values
stopifnot(abs(xu$Lbar - 30) < 1e-9, abs(xu$L25 - 20) < 1e-9, abs(xu$L95 - 48) < 1e-9,
          abs(xu$Lmax5 - 49) < 1e-9, abs(xu$Pmat - 0.5) < 1e-9, xu$Lc == 10.5)

# from a fit; relative indicators equal 1 on average in the reference period
fit <- fiticc(lfd, stocklen(lfd, lhpar, m_model = "constant"), sel_fun = c("dsnormal", "logistic"),
              catch_by_gear = c(0.95, 0.05))
xf <- lbi_flicc(fit)
stopifnot(nrow(xf$values) == nrow(v[v$gear != "pooled", ]))
rr <- x$relative[x$relative$year %in% 2013:2015, ]
stopifnot(all(abs(tapply(rr$Lmean, rr$gear, mean) - 1) < 1e-9, na.rm = TRUE))
for (ty in c("traffic", "trend", "relative")) ggplot2::ggplot_build(plot_lbi(x, ty))
# Cope and Punt decision tree
cp <- FLicc:::cope_punt_flicc(Pmat = c(0.1, 0.5, 0.3, 0.97, 1, 1), Popt = c(0, 0.2, 0.2, 0.3, 0.5, 1),
                              Pmega = c(0, 0.1, 0.1, 0.2, 0.5, 0), Lmat_Lopt = 0.7)
print(cp)
stopifnot(cp$cp_status == c("SB < target", "SB >= target", "SB < target", "SB >= target",
                            "SB >= target", "no information"))
print(table(x$values$gear, x$values$cp_status))
# three-colour status consistent with reference and limit
rr2 <- x$ratios[is.finite(x$ratios$ref), ]
stopifnot(all(rr2$status[rr2$value >= rr2$ref] == "green"),
          all(rr2$status[rr2$value < rr2$lim] == "red"),
          all(rr2$status[rr2$value >= rr2$lim & rr2$value < rr2$ref] == "orange"))
xl <- lbi_flicc(lfd, lhpar, limits = c("Lmean/LFM" = 0.95))
stopifnot(all(xl$ratios$lim[xl$ratios$indicator == "Lmean/LFM"] == 0.95))
print(table(x$ratios$indicator, x$ratios$status))
# pooling: explicit pool, numbers, FLQuants by year, from a fit
xp1 <- lbi_flicc(lfd, lhpar, pool = TRUE, catch_by_gear = c(0.95, 0.05))
stopifnot(isTRUE(all.equal(xp1$values[xp1$values$gear == "pooled", "Lmean"],
                           x$values[x$values$gear == "pooled", "Lmean"])))
xp2 <- lbi_flicc(lfd, lhpar, catch_by_gear = c(Trawl = 0.95, Gillnet = 0.05), catch_units = "numbers")
yrs <- dimnames(lfd[[1]])$year
cq <- FLQuants(Trawl = FLQuant(0.95, dimnames = list(year = yrs)), Gillnet = FLQuant(0.05, dimnames = list(year = yrs)))
xp3 <- lbi_flicc(lfd, lhpar, catch_by_gear = cq)
stopifnot(isTRUE(all.equal(xp3$values$Lmean, x$values$Lmean)))
xp4 <- lbi_flicc(fit, pool = TRUE)                     # fit's catch shares
stopifnot("pooled" %in% xp4$values$gear)
xp5 <- suppressMessages(lbi_flicc(lfd, lhpar, pool = TRUE))   # equal shares
stopifnot("pooled" %in% xp5$values$gear, !"pooled" %in% lbi_flicc(lfd, lhpar)$values$gear)
cat("lbi tests passed\n")
