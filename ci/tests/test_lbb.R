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
