# Build data/lbspr_test.rda: LBSPR's multi-year example in FLicc format (written to the CI output)
library(FLicc)
od <- if (dir.exists("out/tests")) "out/tests" else "."
lbspr_csv <- "
Length,2010,2011,2012,2013,2014,2015
2.5,0,0,0,0,0,0
7.5,0,0,0,0,0,0
12.5,0,0,0,0,0,0
17.5,0,6,0,0,0,0
22.5,0,0,0,0,0,0
27.5,6,6,6,0,0,6
32.5,0,6,6,18,12,18
37.5,36,42,30,60,30,36
42.5,78,96,66,42,66,78
47.5,126,102,126,102,114,150
52.5,156,174,222,192,174,180
57.5,252,246,216,228,288,198
62.5,192,192,144,174,234,168
67.5,126,162,150,180,156,168
72.5,180,138,210,114,108,132
77.5,138,84,114,132,78,90
82.5,48,72,66,78,108,96
87.5,36,48,30,48,24,48
92.5,42,24,24,42,24,36
97.5,6,24,12,12,0,24
102.5,12,12,12,12,18,0
107.5,0,0,0,0,0,6
112.5,0,0,0,0,0,0
117.5,0,0,0,0,0,0
122.5,0,0,0,0,0,0
127.5,0,0,0,0,0,0"
lb <- read.csv(text = lbspr_csv, check.names = FALSE)
# length classes: lower bounds (the file gives mid-points of 5-mm classes)
lfd_lbspr <- FLQuants(LBSPR = FLQuant(as.matrix(lb[, -1]),
                      dimnames = list(len = lb$Length - 2.5, year = names(lb)[-1]), units = "mm"))
# LBSPR vignette life history; k only scales F and M (F/M and SPR depend on M/k = 1.5)
lhpar_lbspr <- FLPar(linf = 100, k = 0.2, M = 0.3, L50 = 66, L95 = 70, a = 1e-5, b = 3)
# LBSPR estimates for these data (LBSPRfit()): yearly point estimates
# (@SL50, @SL95, @FM, @SPR) and the Kalman-smoothed estimates (@Ests)
lbspr_est <- data.frame(year = 2010:2015,
                        SL50 = c(49.10, 48.92, 48.47, 49.43, 50.05, 46.35),
                        SL95 = c(62.90, 64.58, 61.44, 64.73, 63.37, 59.79),
                        FM   = c(0.56, 0.55, 0.56, 0.49, 0.76, 0.37),
                        SPR  = c(0.4372388, 0.4390794, 0.4322175, 0.4830121, 0.3507105, 0.5404657),
                        SL50_smooth = c(48.69, 48.70, 48.69, 48.70, 48.63, 48.42),
                        SL95_smooth = c(62.82, 62.88, 62.76, 62.78, 62.60, 62.34),
                        FM_smooth   = c(0.55, 0.55, 0.55, 0.55, 0.55, 0.54),
                        SPR_smooth  = c(0.44, 0.44, 0.44, 0.45, 0.45, 0.45))
save(lfd_lbspr, lhpar_lbspr, lbspr_est, file = file.path(od, "lbspr_test.rda"), compress = "xz")
# check
rm(lfd_lbspr, lhpar_lbspr, lbspr_est)
load(file.path(od, "lbspr_test.rda"))
print(lfd_lbspr); print(lhpar_lbspr); print(lbspr_est)
stk <- stocklen(lfd_lbspr, lhpar_lbspr, m_model = "constant")
f <- fiticc(lfd_lbspr, stk, sel_fun = "logistic", catch_by_gear = 1, settings = list(CVL = 0.1, obs_model = "mn"))
print(sprcur_flicc(f))
