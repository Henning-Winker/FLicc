# binormal selectivity: two gears pooled into one bimodal composition
library(FLicc)
data("alfonsino")
lfd <- FLQuantLen(lfd_long_to_wide(as.data.frame(lfd_alfonsino)), unit = "cm", midL = FALSE)
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859, a = 0.004721956 / 1000, b = 3.146168)
lfdc <- FLQuants(Pooled = lfd[[1]] + lfd[[2]])
stk <- stocklen(lfdc, lhpar, m_model = "constant")
od <- if (dir.exists("out/tests")) "out/tests" else "."
set <- list(CVL = 0.1)
f_ds <- fiticc(lfdc, stk, sel_fun = "dsnormal", catch_by_gear = 1, settings = set)
f_bn <- fiticc(lfdc, stk, sel_fun = "binormal", catch_by_gear = 1, settings = set)
print(f_bn$diagnostics)
sp <- selpars_flicc(f_bn); print(sp)
tab <- data.frame(model = c("dsnormal", "binormal"),
                  conv = c(f_ds$diagnostics$conv_tier, f_bn$diagnostics$conv_tier),
                  AIC = c(AIC(LLflicc(f_ds)), AIC(LLflicc(f_bn))),
                  SPR_last = c(tail(c(f_ds$report$spr), 1), tail(c(f_bn$report$spr), 1)))
print(tab)
p <- as.list(as.numeric(sp$Pooled)); names(p) <- dimnames(sp$Pooled)$params
stopifnot(p$mode2 > p$mode1, all(is.finite(unlist(p))))
# R curve matches the template's selectivity
L <- as.numeric(f_bn$tmb_data$Lmid)
s_r <- sel_curve_flicc(L, "binormal", p)
s_t <- as.numeric(f_bn$obj$report(f_bn$obj$env$last.par.best)$Sel[, 1])
cat("max |R - TMB| selectivity:", max(abs(s_r - s_t)), "\n")
stopifnot(max(abs(s_r - s_t)) < 1e-6, abs(max(s_t) - 1) < 1e-8)
# sel_fixed with binormal names
f_fx <- fiticc(lfdc, stk, sel_fun = "binormal", catch_by_gear = 1, settings = set,
               sel_fixed = list(Pooled = list(h = 1, sd2 = p$sd2)))
print(selpars_flicc(f_fx)); stopifnot(abs(as.numeric(selpars_flicc(f_fx)$Pooled["h"]) - 1) < 1e-8)
ggplot2::ggsave(file.path(od, "binormal_len.png"),
                plot_len(f_bn), width = 8, height = 4, dpi = 100)
ggplot2::ggsave(file.path(od, "binormal_sel.png"),
                plot_compare_sel(list(dsnormal = f_ds, binormal = f_bn)), width = 6, height = 4, dpi = 100)
ggplot2::ggsave(file.path(od, "binormal_resid.png"),
                plot_resid_len(f_ds), width = 6, height = 4, dpi = 100)

# bidsnormal: second mode double-sided; rsd2 = 1000 gives an asymptotic second mode
f_bd <- fiticc(lfdc, stk, sel_fun = "bidsnormal", catch_by_gear = 1, settings = set,
               sel_fixed = list(Pooled = list(rsd2 = 1000)))
print(f_bd$diagnostics)
spd <- selpars_flicc(f_bd); print(spd)
pd <- as.list(as.numeric(spd$Pooled)); names(pd) <- dimnames(spd$Pooled)$params
stopifnot(abs(pd$rsd2 - 1000) < 1e-6, pd$mode2 > pd$mode1)
s_r <- sel_curve_flicc(L, "bidsnormal", pd)
s_t <- as.numeric(f_bd$obj$report(f_bd$obj$env$last.par.best)$Sel[, 1])
cat("bidsnormal max |R - TMB|:", max(abs(s_r - s_t)), "\n")
stopifnot(max(abs(s_r - s_t)) < 1e-6)
tab2 <- rbind(tab, data.frame(model = "bidsnormal (rsd2 = 1000)", conv = f_bd$diagnostics$conv_tier,
                              AIC = AIC(LLflicc(f_bd)), SPR_last = tail(c(f_bd$report$spr), 1)))
print(tab2)
ggplot2::ggsave(file.path(od, "binormal_sel.png"),
                plot_compare_sel(list(dsnormal = f_ds, binormal = f_bn, bidsnormal = f_bd)),
                width = 6, height = 4, dpi = 100)
