# OM simulator for FLicc-SPM tests: om_simulate() of om_sim.R with exploitable
# biomass by gear (start of year, for CPUE), true SB/SB0, B/Bmsy and F/Fmsy.
# Needs om_sim.R sourced first.

om_simulate_spm <- function(st, years = 1971:2024, data_years = 2012:2024, sigR = 0, rho = 0.3,
                            fshape = "fudc", n_age = c(Trawl = 50, Gillnet = 30),
                            n_len = c(Trawl = 150, Gillnet = 100),
                            index_gear = "Gillnet", index_years = 1990:2024, index_cv = 0.2, q = 1e-3) {
  ny <- length(years); cl <- st$cells; lh <- st$lh
  Fmsy <- om_fmsy(st)
  Fbar <- om_fhist(ny, Fmsy, shape = fshape)
  eps <- numeric(ny)
  for (y in seq_len(ny)) eps[y] <- (if (y > 1) rho * eps[y - 1] else 0) + sqrt(1 - rho^2) * rnorm(1, 0, sigR)
  dev <- exp(eps - sigR^2 / 2)
  sb0 <- om_sbpr(st, 0); h <- lh$h
  bh <- function(ssb) 4 * h * ssb / (sb0 * (1 - h) + (5 * h - 1) * ssb)
  sbmsy <- { sb <- om_sbpr(st, Fmsy); max((4 * h * sb - sb0 * (1 - h)) / ((5 * h - 1) * sb), 0) * sb }
  s_cell <- sweep(st$Fsel, 2, st$Fap, "/")            # selectivity of each cell by gear (max 1 at length)
  sj <- rowSums(st$Fsel); sj <- sj / max(sj)          # joint fishery selectivity (apical 1)
  eqH <- function(F) {                                # equilibrium harvest rates at F (start / mean of year)
    N <- om_eq(st, F); Fc <- F * rowSums(st$Fsel); Z <- lh$M + Fc
    C <- sum(N * Fc / Z * (1 - exp(-Z)) * cl$wt)
    c(start = C / sum(N * cl$wt * sj), mean = C / sum(N * cl$wt * sj * (1 - exp(-Z)) / Z))
  }
  Hmsy <- eqH(Fmsy)
  sb1 <- om_sbpr(st, Fbar[1]); R1 <- max((4 * h * sb1 - sb0 * (1 - h)) / ((5 * h - 1) * sb1), 0.05)
  N <- R1 * om_eq(st, Fbar[1])
  key <- paste(cl$g, cl$c)
  catch_cell <- array(0, c(nrow(cl), length(st$gears), ny))
  ssb <- numeric(ny); eb <- matrix(0, ny, length(st$gears), dimnames = list(years, st$gears))
  ebj <- ebjm <- numeric(ny)
  for (y in seq_len(ny)) {
    ssb[y] <- sum(N * cl$mat * cl$wt)
    N[cl$a == 0] <- bh(ssb[y]) * dev[y] * cl$w[cl$a == 0]
    eb[y, ] <- colSums(N * cl$wt * s_cell)
    Fc <- Fbar[y] * st$Fsel; Z <- lh$M + rowSums(Fc)
    ebj[y] <- sum(N * cl$wt * sj); ebjm[y] <- sum(N * cl$wt * sj * (1 - exp(-Z)) / Z)
    catch_cell[, , y] <- N * Fc / Z * (1 - exp(-Z))
    Ns <- N * exp(-Z); Nn <- numeric(length(N))
    for (k in unique(key)) {
      i <- which(key == k); i <- i[order(cl$a[i])]; n <- length(i)
      Nn[i[2:n]] <- Ns[i[1:(n - 1)]]; Nn[i[n]] <- Nn[i[n]] + Ns[i[n]]
    }
    N <- Nn
  }
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
  ki <- match(index_years, years)
  index <- setNames(q * eb[ki, index_gear] * exp(rnorm(length(ki), -index_cv^2 / 2, index_cv)), index_years)
  list(years = years, data_years = data_years, Fbar = Fbar, Fmsy = Fmsy, ssb = ssb, dev = dev,
       sb0 = sb0, sbmsy = sbmsy, depletion = ssb / sb0, bbmsy = ssb / sbmsy, ffmsy = Fbar / Fmsy,
       spr = sapply(Fbar, function(F) om_spr(st, F)), lfd = lfd, catch_w = cw, lb = st$lb,
       eb = eb, index = index, catch = setNames(rowSums(cw), years),
       Hmsy = Hmsy, H_start = rowSums(cw) / ebj, H_mean = rowSums(cw) / ebjm,
       hhmsy_start = rowSums(cw) / ebj / Hmsy["start"], hhmsy_mean = rowSums(cw) / ebjm / Hmsy["mean"])
}
