# combine_lfd(): pooling gears raised by catch
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)[, c("len", "year", "data", "qname")]),
                  unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
print(names(lfd))
cf <- combine_lfd(lfd, c(0.7, 0.3), lhpar, ess.n = 500)
stopifnot(inherits(cf, "FLQuants"), names(cf) == "fleet",
          all(abs(apply(cf[[1]]@.Data, 2, sum) - 500) < 1e-6))
# all catch on one gear -> its own composition
cf1 <- combine_lfd(lfd, c(1, 0), lhpar, ess.n = 500)
p1 <- lfdess(lfd, 500)[[1]]
stopifnot(max(abs(cf1[[1]]@.Data - p1@.Data)) < 1e-8)
# named and FLQuants catch give the same as the vector
cf2 <- combine_lfd(lfd, setNames(c(0.7, 0.3), names(lfd)), lhpar)
cq <- FLQuants(lapply(c(0.7, 0.3), function(v) FLQuant(v, dimnames = list(year = dimnames(lfd[[1]])$year))))
names(cq) <- names(lfd)
cf3 <- combine_lfd(lfd, cq, lhpar)
stopifnot(max(abs(cf2[[1]]@.Data - cf[[1]]@.Data)) < 1e-8, max(abs(cf3[[1]]@.Data - cf[[1]]@.Data)) < 1e-8)
# one-fleet fit with two-mode selectivity, last year
yr <- tail(dimnames(lfd[[1]])$year, 1)
lf <- window(cf, start = as.numeric(yr))
stk <- stocklen(lf, lhpar, m_model = "constant")
for (sf in c("logistic", "binormal", "bidsnormal")) {
  f <- fiticc(lf, stk, sel_fun = sf, catch_by_gear = 1)
  cat(sf, ": SPR", round(sprcur_flicc(f), 3), " conv", f$diagnostics$conv_tier, "\n")
}
