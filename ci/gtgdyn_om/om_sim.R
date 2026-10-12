# OM-like truth simulator for gtg.dyn tests (base R; mimics MSE_OM_Alfonsino_tune_v2)
#
# Annual age-structured model (ages 0-20+) for splendid alfonsino with
#  - growth: non-persistent (one group, length-at-age CV 0.1, as the OM) or
#    persistent growth-type groups (Linf_g = Linf (1 + 0.1 z), as FLicc's GTG)
#  - recruitment: annual pulse at age 0 (as the OM) or spread within the year
#    (4 sub-cohorts; approximately continuous, as FLicc)
#  - Beverton-Holt recruitment (h = 0.75), optional AR1 log deviations
#  - two gears: Trawl dome (dnormal 23.35, 3.19, 11.22; f_mult 0.9),
#    Gillnet logistic (L50 29.7, L95 34.1; f_mult 0.1); F partitioned on the
#    Fbar scale (ages 3-19) as in build_gear() / f_age_gear()
#  - F history shaped like FLRef::fudc(fhi = 3, flo = 1.5, breaks = c(0.6, 0.8)),
#    lognormal year noise 0.15 and an iteration multiplier (lognormal 0.2), or flat
#  - two-stage length sampling at t = 0.5 (age-cell sample of n_age, then a
#    length sample of n_len) by gear, as lfd.sim()
# Truth: SPR at the year's F (same cell structure, per recruit).

om_lh <- list(linf = 55.662, k = 0.08, t0 = -1, M = 0.162, a = 3.24e-05, b = 2.932,
              tm = 6, h = 0.75)

om_structure <- function(persist = FALSE, continuous = FALSE, lh = om_lh) {
  ages <- 0:20
  lb <- 5:67                                   # lower bounds, last bin plus group
  if (persist) {
    z <- seq(-3, 3, length.out = 19); wz <- dnorm(z); wz <- wz / sum(wz)
    linf_g <- lh$linf * (1 + 0.1 * z); cvw <- 0.02
  } else {
    wz <- 1; linf_g <- lh$linf; cvw <- 0.1
  }
  nsub <- if (continuous) 4 else 1
  delta <- if (nsub == 1) 0 else 0.5 - (seq_len(nsub) - 0.5) / nsub   # age shift at sampling
  cells <- expand.grid(g = seq_along(linf_g), c = seq_len(nsub), a = ages)
  cells$w <- wz[cells$g] / nsub                               # recruitment share
  agef <- pmax(cells$a + 0.5 + delta[cells$c], 0.05)          # age at mid-year sampling
  cells$L <- linf_g[cells$g] * (1 - exp(-lh$k * (agef - lh$t0)))
  cells$wt <- lh$a * cells$L^lh$b
  cells$mat <- ifelse(cells$a == 0, 0, 1 / (1 + exp(-log(19) * (cells$a - lh$tm) / (0.1 * lh$tm))))
  # length distribution of each cell over the bins (below the first bin dropped)
  up <- c(lb[-1], Inf)
  P <- t(sapply(seq_len(nrow(cells)), function(i) {
    sd <- cvw * cells$L[i]
    pnorm(up, cells$L[i], sd) - pnorm(lb, cells$L[i], sd)
  }))
  lmid <- lb + 0.5
  sel_len <- list(
    Trawl   = { d <- lmid - 23.35; exp(-d^2 / (2 * ifelse(d < 0, 3.19, 11.22)^2)) },
    Gillnet = 1 / (1 + exp(-log(19) * (lmid - 29.7) / (34.1 - 29.7))))
  sel_len <- lapply(sel_len, function(s) s / max(s))
  f_mult <- c(Trawl = 0.9, Gillnet = 0.1)
  # selectivity of each cell by gear and conditional length distribution
  s_cell <- sapply(sel_len, function(s) as.numeric(P %*% s))
  cond <- lapply(sel_len, function(s) { x <- sweep(P, 2, s, "*"); x / pmax(rowSums(x), 1e-300) })
  # Fbar-scale partition: gear selectivity-at-age (population mean over cells)
  # scaled to mean 1 over ages 3-19
  fb <- cells$a %in% 3:19
  scale_g <- sapply(names(sel_len), function(g) {
    sa <- tapply(cells$w * s_cell[, g], cells$a, sum) / tapply(cells$w, cells$a, sum)
    mean(sa[as.character(3:19)])
  })
  Fsel <- sweep(s_cell, 2, f_mult / scale_g, "*")   # F_cell,gear = Fbar * Fsel[, gear]
  list(cells = cells, lb = lb, P = P, cond = cond, Fsel = Fsel, ages = ages, lh = lh,
       gears = names(sel_len), nsub = nsub, ng = length(linf_g))
}

# equilibrium numbers per unit recruitment by cell at Fbar
om_eq <- function(st, Fbar) {
  cl <- st$cells; Z <- st$lh$M + Fbar * rowSums(st$Fsel)
  N <- numeric(nrow(cl))
  key <- paste(cl$g, cl$c)
  for (k in unique(key)) {
    i <- which(key == k); i <- i[order(cl$a[i])]
    s <- exp(-c(0, cumsum(Z[i][-length(i)])))
    s[length(i)] <- s[length(i)] / (1 - exp(-Z[i][length(i)]))
    N[i] <- cl$w[i] * s
  }
  N
}
om_sbpr <- function(st, Fbar) sum(om_eq(st, Fbar) * st$cells$mat * st$cells$wt)
om_spr  <- function(st, Fbar) om_sbpr(st, Fbar) / om_sbpr(st, 0)

om_fmsy <- function(st) {
  h <- st$lh$h; sb0 <- om_sbpr(st, 0)
  ypr <- function(F) { N <- om_eq(st, F); Fc <- F * rowSums(st$Fsel); Z <- st$lh$M + Fc
    sum(N * Fc / Z * (1 - exp(-Z)) * st$cells$wt) }
  req <- function(F) { sb <- om_sbpr(st, F); max((4 * h * sb - sb0 * (1 - h)) / ((5 * h - 1) * sb), 0) }
  # grid first: the objective is flat (0) beyond the crash F, which can trap
  # a bounded line search
  obj <- function(F) ypr(F) * req(F)
  fg <- seq(0.005, 1, by = 0.005); i <- which.max(sapply(fg, obj))
  optimize(function(F) -obj(F), c(fg[max(i - 1, 1)], fg[min(i + 1, length(fg))]))$minimum
}

om_fhist <- function(ny, Fmsy, shape = "fudc", fhi = 3, flo = 1.5, f0 = 0.2,
                     breaks = c(0.6, 0.8), sdF = 0.15, sdIter = 0.2) {
  if (shape == "flat") {
    m <- rep(flo, ny)
  } else {
    i <- seq_len(ny); b1 <- round(breaks[1] * ny); b2 <- round(breaks[2] * ny)
    m <- approx(c(1, b1, b2, ny), c(f0, fhi, flo, flo), xout = i)$y
  }
  m * Fmsy * exp(rnorm(ny, 0, sdF)) * exp(rnorm(1, 0, sdIter))
}

om_simulate <- function(st, years = 1971:2024, data_years = 2012:2024, sigR = 0, rho = 0.3,
                        fshape = "fudc", n_age = c(Trawl = 50, Gillnet = 30),
                        n_len = c(Trawl = 150, Gillnet = 100)) {
  ny <- length(years); cl <- st$cells; lh <- st$lh
  Fmsy <- om_fmsy(st)
  Fbar <- om_fhist(ny, Fmsy, shape = fshape)
  eps <- numeric(ny)
  for (y in seq_len(ny)) eps[y] <- (if (y > 1) rho * eps[y - 1] else 0) + sqrt(1 - rho^2) * rnorm(1, 0, sigR)
  dev <- exp(eps - sigR^2 / 2)
  sb0 <- om_sbpr(st, 0); h <- lh$h
  bh <- function(ssb) 4 * h * ssb / (sb0 * (1 - h) + (5 * h - 1) * ssb)
  # initial equilibrium at the first year's F
  sb1 <- om_sbpr(st, Fbar[1]); R1 <- max((4 * h * sb1 - sb0 * (1 - h)) / ((5 * h - 1) * sb1), 0.05)
  N <- R1 * om_eq(st, Fbar[1])
  key <- paste(cl$g, cl$c); amax <- max(st$ages)
  catch_cell <- array(0, c(nrow(cl), length(st$gears), ny))
  ssb <- numeric(ny)
  for (y in seq_len(ny)) {
    ssb[y] <- sum(N * cl$mat * cl$wt)
    N[cl$a == 0] <- bh(ssb[y]) * dev[y] * cl$w[cl$a == 0]
    Fc <- Fbar[y] * st$Fsel; Z <- lh$M + rowSums(Fc)
    catch_cell[, , y] <- N * Fc / Z * (1 - exp(-Z))
    Ns <- N * exp(-Z); Nn <- numeric(length(N))
    for (k in unique(key)) {
      i <- which(key == k); i <- i[order(cl$a[i])]; n <- length(i)
      Nn[i[2:n]] <- Ns[i[1:(n - 1)]]; Nn[i[n]] <- Nn[i[n]] + Ns[i[n]]
    }
    N <- Nn
  }
  # two-stage length samples in the data years
  kd <- match(data_years, years)
  lfd <- setNames(lapply(seq_along(st$gears), function(g) {
    sapply(kd, function(y) {
      p <- catch_cell[, g, y]; p <- p / sum(p)
      na <- as.numeric(rmultinom(1, n_age[st$gears[g]], p)) / n_age[st$gears[g]]
      pl <- as.numeric(na %*% st$cond[[g]])
      as.numeric(rmultinom(1, n_len[st$gears[g]], pl / sum(pl)))
    })
  }), st$gears)
  cw <- sapply(seq_along(st$gears), function(g) colSums(catch_cell[, g, ] * cl$wt))
  colnames(cw) <- st$gears
  list(years = years, data_years = data_years, Fbar = Fbar, Fmsy = Fmsy, ssb = ssb, dev = dev,
       spr = sapply(Fbar, function(F) om_spr(st, F)), lfd = lfd, catch_w = cw, lb = st$lb)
}
