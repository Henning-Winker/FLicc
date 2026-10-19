# LBB as a special case of FLicc
library(FLicc)
options(width = 200)
od <- if (dir.exists("out/tests")) "out/tests" else "."

# 1. Port of LBB's BH(): values for the check case in LBB_33a.R comments
#    (Linf 120, Lc/Linf 0.2917, alpha 60, M/K 1.5385, F/K 0.7692: Y/R' 0.0332, B/B0 0.467 in LBB)
for (w in c(0.5, 1, 2.5)) print(c(width = w, lbb_bh_flicc(seq(22.5, 120, by = w), 120, 1.5385, 0.7692, FALSE, 0.2917, 60)))

# 2. Alfonsino gillnet (logistic), LBB defaults and user Linf/M/K
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
gn <- lfd["Gillnet"]
pr <- lbb_priors_flicc(gn); str(pr[c("Linf", "Linf.sd", "ZK", "MK", "Lc", "alpha", "Lstart", "Lmax", "Lmax.med")])
t1 <- system.time(lb <- lbb_flicc(gn, Lm50 = 31.2))
print(lb); cat("sec", t1[["elapsed"]], "\n")
print(lb$table[, c("year", "Linf", "Lc", "alpha", "MK", "FK", "FM", "FM_lcl", "FM_ucl", "BB0", "BB0_lcl", "BB0_ucl", "YR", "SPR", "Lmean_Lopt", "Lc_Lc_opt")])
lbu <- lbb_flicc(gn, Linf = 55.7, MK = 0.162 / 0.08, Lm50 = 31.2)
print(lbu)
lbj <- lbb_flicc(gn, Linf = 55.7, Lm50 = 31.2, by_year = FALSE)
print(lbj)
stopifnot(inherits(lb, "flicc_lbb"), nrow(lb$table) == 12, all(is.finite(lb$table$BB0)),
          all(lb$table$BB0 > 0), nrow(lbj$table) == 12)
ggplot2::ggsave(file.path(od, "lbb_alfonsino.png"), plot_lbb(lb), width = 10, height = 7, dpi = 100)
ggplot2::ggsave(file.path(od, "lbb_alfonsino_joint.png"), plot_lbb(lbj), width = 10, height = 7, dpi = 100)

# 3. LBB example data (Illex coindetii, Adriatic), settings of Example_ID.csv:
#    Years.user 1996,1999,2001,2002,2009,2014,2016; Linf.user 23 cm; Lm50 16 cm; lengths in mm
url <- "https://raw.githubusercontent.com/SISTA16/LBB/master/Example_Dat.csv"
dat <- tryCatch(utils::read.csv(url), error = function(e) NULL)
if (!is.null(dat)) {
  d <- dat[dat$Stock == "Ille_coi_AD" & !is.na(dat$CatchNo), ]
  df <- aggregate(data ~ len + year + qname,
                  data.frame(len = d$Length / 10, year = d$Year, data = d$CatchNo, qname = "Illex"), sum)
  il <- FLQuantLen(lfd_long_to_wide(df), unit = "cm", midL = FALSE)
  li <- lbb_flicc(il, Linf = 23, Lm50 = 16, years = c(1996, 1999, 2001, 2002, 2009, 2014, 2016))
  print(li)
  print(li$table[, c("year", "Linf", "Lc", "MK", "FM", "ZK", "BB0", "YR", "Lmean_Lopt")])
  li0 <- lbb_flicc(il, Lm50 = 16, years = c(1996, 1999, 2001, 2002, 2009, 2014, 2016))
  print(li0)
  ggplot2::ggsave(file.path(od, "lbb_illex.png"), plot_lbb(li), width = 10, height = 7, dpi = 100)
} else cat("LBB example data not reachable; skipped.\n")

# 4. EPR from the FLicc engine vs LBB's B/B0 (closed form), and plot_epr
print(lb$table[, c("year", "BB0", "EPR", "SPR")])
stopifnot(all(is.finite(lb$table$EPR)), cor(lb$table$BB0, lb$table$EPR) > 0.9)
data("alfonsino")
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
stk <- stocklen(lfd, lhpar, m_model = "constant")
fit2 <- fiticc(lfd, stk, sel_fun = c("dsnormal", "logistic"), catch_by_gear = c(0.7, 0.3),
               settings = list(CVL = 0.1, obs_model = "dm", dm_theta = "saturating"))
print(epr_ci_flicc(fit2)); print(spr_ci_flicc(fit2)[1:3, ])
ggplot2::ggsave(file.path(od, "epr_alfonsino.png"), plot_epr(fit2, ci = TRUE), width = 7, height = 4, dpi = 100)
ggplot2::ggsave(file.path(od, "epr_lbb.png"),
                plot_epr(lb$fits[["2024"]], ref_levels = attr(lb$table, "BFM1B0") * c(1, 0.5)), width = 7, height = 4, dpi = 100)

# 5. LBB-like weighting: ESS 20000 with fixed DM (LBB Dirichlet concentration 10^4)
lbe <- lbb_flicc(gn, Lm50 = 31.2, ess = 20000)
print(lbe)

# 6. LBB summary for any FLicc fit (two gears, dome trawl)
tb2 <- lbbref_flicc(fit2)
print(tb2[, c("year", "Linf", "Lc", "alpha", "MK", "FK", "FM", "Lopt", "Lc_opt", "Lmean", "BB0", "BB0_lcl", "BB0_ucl", "SPR")])
cat("BFM1B0", attr(tb2, "BFM1B0"), "\n")
stopifnot(nrow(tb2) == 12, all(is.finite(tb2$BB0)), all(is.finite(tb2$Lc)))
ggplot2::ggsave(file.path(od, "lbb_anyfit.png"), plot_lbb(fit2), width = 10, height = 7, dpi = 100)
