# make_hake.R (FLicc data(hake))
# FLicc inputs for southern hake (Merluccius merluccius, hke.27.8c9a) from the SS3 assessment
# (SS3 3.30.22, two sexes, 4 seasons, 1960-2024), read with r4ss and ss3diags2.
#
# Fleets in SS3: 4 fisheries (trawlers with retained and discard partitions; volpal, artisanal,
# cdTrw whole catch), 3 surveys with length compositions (SpSurv and PtSurv sexed, CdSurv unsexed),
# and 2 CPUE indices (SpCPUE_trawlers, SpCPUE_volpal) whose selectivity MIRRORS trawlers and volpal
# (checked: identical size selectivity). The CPUE fleets have no length compositions; their
# selectivity is that of the mirrored fisheries, so they are not separate FLicc gears.
#
# Output (data/hake.rda):
#   lfd_hake    FLQuants, 1994-2024: fisheries unit "U" (annual catch at length: quarterly compositions
#               raised by SS3 catch numbers by season and partition), SpSurv and PtSurv units "F","M"
#               (joint sexed samples; below 21 cm the sexes are combined and stored as F, as in SS3),
#               CdSurv unit "U"
#   catch_hake  FLQuants, dead catch weight (t) of the four fisheries by year
#   lhpar_hake  FLPar params x unit (F, M): growth, maturity, weight, M (inverse M at length with
#               reflen = female L50) from SS3; male growth and M from the SS3 exponential offsets
#   ss3_hake    list: SS3 SPR, SSB/SSB0 and F by year; selectivity at length (females, last year);
#               sex_minlen of the sexed surveys; CPUE fleets and the fishery they mirror; notes
# Run from the FLicc root with ss3file set to the r4ss output of the assessment.
suppressMessages({library(FLCore); library(data.table); library(ss3diags2)})
ss3file <- if (exists("ss3file")) ss3file else "hke.27.8c9a.rds"

ss <- readRDS(ss3file)    # r4ss SS_output() of the SS3 assessment (not distributed with FLicc)
fisheries <- c("trawlers", "volpal", "artisanal", "cdTrw")
surveys_sex <- c("SpSurv", "PtSurv"); surveys_u <- "CdSurv"
yrs <- 1994:2024

# --- mirrored selectivity (CPUE fleets) ---
sz <- as.data.table(ss$sizeselex)[Factor == "Lsel" & Yr == ss$endyr & Sex == 1]
selm <- as.matrix(sz[, -(1:5)]); rownames(selm) <- sz$Fleet
fl <- ss$FleetNames
mirror <- do.call(rbind, lapply(which(ss$fleet_type == 3 & !(fl %in% c(surveys_sex, surveys_u))), function(f) {
  m <- which(sapply(seq_len(nrow(selm)), function(g) g != f && isTRUE(all.equal(selm[ac(f), ], selm[ac(g), ]))))
  data.frame(fleet = fl[f], mirrors = if (length(m)) fl[m[1]] else NA)
}))

# --- length compositions (proportions x adjusted sample size) ---
ld <- as.data.table(ss$lendbase)
ld <- ld[, .(Yr, Seas, Fleet, Sex, Part, len = Bin, obs = Obs, nsamp = Nsamp_adj)]
ld[, name := fl[Fleet]]
lens <- sort(unique(ld$len))

# fishery catch numbers by season and partition (retained, discarded)
ct <- as.data.table(ss$catch)[Fleet_Name %in% fisheries,
        .(name = Fleet_Name, Yr, Seas, ret = ret_num, dis = pmax(sel_num - ret_num, 0), dead_bio)]

# annual catch at length of a fishery: each sampled season x partition composition raised by its
# catch numbers; a year is used only if every season-partition with catch was sampled
annual_fishery <- function(f) {
  d <- ld[name == f & Yr %in% yrs]
  d[, prop := obs / sum(obs), by = .(Yr, Seas, Part)]
  cf <- melt(ct[name == f], id.vars = c("name", "Yr", "Seas"), measure.vars = c("ret", "dis"),
             variable.name = "pt", value.name = "num")
  cf[, Part := ifelse(pt == "dis", 1L, ifelse(f == "trawlers", 2L, 0L))]
  cf <- cf[num > 0 & !(f != "trawlers" & pt == "dis")]
  smp <- unique(d[, .(Yr, Seas, Part, nsamp)])
  cov <- merge(cf, smp, by = c("Yr", "Seas", "Part"), all.x = TRUE)
  okyrs <- cov[, .(ok = all(!is.na(nsamp))), by = Yr][ok == TRUE, Yr]
  d <- merge(d[Yr %in% okyrs], cf[, .(Yr, Seas, Part, num)], by = c("Yr", "Seas", "Part"))
  a <- d[, .(n = sum(prop * num)), by = .(Yr, len)]
  ns <- smp[Yr %in% okyrs, .(ess = sum(nsamp)), by = Yr]
  a <- merge(a, ns, by = "Yr")[, n := n / sum(n) * ess, by = Yr]
  q <- FLQuant(0, dimnames = list(len = lens, year = yrs, unit = "U"), units = "1")
  for (y in okyrs) q[match(a[Yr == y, len], lens), ac(y), "U"] <- a[Yr == y, n]
  q
}

survey_q <- function(f, sexed) {
  d <- ld[name == f & Yr %in% yrs]
  units <- if (sexed) c("F", "M") else "U"
  q <- FLQuant(0, dimnames = list(len = lens, year = yrs, unit = units), units = "1")
  for (i in seq_len(nrow(d))) {
    u <- if (sexed) c("F", "M")[d$Sex[i]] else "U"
    q[match(d$len[i], lens), ac(d$Yr[i]), u] <- q[match(d$len[i], lens), ac(d$Yr[i]), u] + d$obs[i] * d$nsamp[i]
  }
  q
}

lfd <- FLQuants(c(setNames(lapply(fisheries, annual_fishery), fisheries),
                  setNames(lapply(surveys_sex, survey_q, sexed = TRUE), surveys_sex),
                  setNames(lapply(surveys_u, survey_q, sexed = FALSE), surveys_u)))
pool <- function(q) { u <- unitSums(q); dimnames(u)$unit <- "U"; u }
lfd_pooled <- FLQuants(lapply(lfd, pool))

catch_by_gear <- FLQuants(setNames(lapply(fisheries, function(f) {
  v <- ct[name == f, .(c = sum(dead_bio)), by = Yr][match(yrs, Yr), c]; v[is.na(v)] <- 0
  FLQuant(v, dimnames = list(age = "all", year = yrs), units = "t")
}), fisheries))

# --- life history by sex ---
p <- as.data.table(ss$parameters)[, .(Label, Value)]
pv <- function(l) p[Label == l, Value]
offset <- function(fem, off) fem * exp(off)          # SS3 offset approach 2 (exponential offsets)
LinfF <- pv("L_at_Amax_Fem_GP_1"); kF <- pv("VonBert_K_Fem_GP_1"); L1F <- pv("L_at_Amin_Fem_GP_1")
A1 <- ss$Growth_Parameters$A1[1]
LinfM <- offset(LinfF, pv("L_at_Amax_Mal_GP_1")); kM <- offset(kF, pv("VonBert_K_Mal_GP_1"))
L1M <- offset(L1F, pv("L_at_Amin_Mal_GP_1"))
t0 <- function(Linf, k, L1) A1 + log(1 - L1 / Linf) / k
L50 <- pv("Mat50%_Fem_GP_1"); slope <- pv("Mat_slope_Fem_GP_1")
L95 <- L50 + log(19) / abs(slope)
# M at age (SS3 breakpoints, offsets for males) -> M at length; FLicc inverse M: M(L) = M * reflen / L
eg <- as.data.table(ss$endgrowth)[!duplicated(paste(Sex, Age_Beg))]
m_at <- function(s, L) approx(eg[Sex == s, Len_Mid], eg[Sex == s, M], L, rule = 2)$y
reflen <- L50
lh <- function(linf, k, t0, M) c(linf = linf, k = k, t0 = t0, M = M, L50 = L50, L95 = L95,
                                 a = pv("Wtlen_1_Fem_GP_1"), b = pv("Wtlen_2_Fem_GP_1"), reflen = reflen)
lhF <- lh(LinfF, kF, t0(LinfF, kF, L1F), m_at(1, reflen))
lhM <- lh(LinfM, kM, t0(LinfM, kM, L1M), m_at(2, reflen))
lhpar <- FLPar(array(c(lhF, lhM), dim = c(length(lhF), 2, 1),
                     dimnames = list(params = names(lhF), unit = c("F", "M"), iter = 1)))

# --- SS3 benchmark ---
dq <- as.data.table(ss$derived_quants)
ssb0 <- dq[Label == "SSB_unfished", Value]
ssb <- dq[grepl("^SSB_[0-9]{4}$", Label), .(year = as.integer(sub("SSB_", "", Label)), ssb = Value)]
sp <- as.data.table(ss$sprseries)[Era == "TIME"]
fstd <- grep("F_std", names(sp), value = TRUE)[1]
ss3 <- merge(data.frame(year = sp$Yr, spr = sp$SPR, F_std = sp[[fstd]]),
             data.frame(year = ssb$year, ssb_ssb0 = ssb$ssb / ssb0), by = "year")
ss3$spr_msy <- dq[Label == "SPR_MSY", Value]

notes <- c(
  "Source: SS3 3.30.22 southern hake (hke.27.8c9a), r4ss output; two sexes, 4 seasons.",
  "Fishery compositions (quarterly; trawlers retained + discards) raised to annual catch at length with SS3 catch numbers by season and partition; years used only when all season-partitions with catch were sampled; sample size = sum of SS3 adjusted sample sizes.",
  "SpSurv and PtSurv are sexed (joint F+M samples, SS3 gender 3); CdSurv unsexed; all surveys in season 4.",
  "SpCPUE_trawlers and SpCPUE_volpal mirror the selectivity of trawlers and volpal and have no length data: not separate FLicc gears.",
  "Male growth from SS3 exponential offsets: Linf 110 F, 73.7 M; k 0.14 both. ss3_lhpar() does not apply the offsets (male = female), so they are applied here.",
  "M declines with age in SS3 (breakpoints, male offsets); FLicc inverse M at length with reflen = female L50 matches it within about 0.03 over the selected sizes.",
  "Maturity and weight from SS3 female parameters (male maturity is not used: SPR is female)."
)

# SS3 selectivity at length in the last year (females), by fleet with length data
szl <- as.data.table(ss$sizeselex)[Factor == "Lsel" & Yr == ss$endyr & Sex == 1 & Fleet <= 7]
selm <- melt(szl[, -(c(1, 3:5))], id.vars = "Fleet", variable.name = "len", value.name = "sel")
selm[, len := as.numeric(as.character(len))][, gear := ss$FleetNames[Fleet]]
# sexes combined below this length in the sexed surveys (smallest length with male counts)
sex_minlen <- ld[name %in% surveys_sex & Sex == 2 & obs > 0, .(m = min(len)), by = name]
sex_minlen <- setNames(sex_minlen$m, sex_minlen$name)[surveys_sex]

lfd_hake <- lfd
catch_hake <- catch_by_gear
lhpar_hake <- lhpar
ss3_hake <- list(spr = ss3, sel = as.data.frame(selm[, .(gear, len, sel)]), sex_minlen = sex_minlen,
                 mirror = mirror, spr_msy = ss3$spr_msy[1], notes = notes)
save(lfd_hake, catch_hake, lhpar_hake, ss3_hake, file = "data/hake.rda", compress = "xz")
