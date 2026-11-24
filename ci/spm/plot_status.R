suppressPackageStartupMessages({ library(FLCore); library(FLicc); library(ggplot2) })
sc <- data.frame(seed = c(2001, 11, 12, 11, 12, 11, 12), shape = c("fudc", "rebuild", "rebuild", "decline", "decline", "flat", "flat"))
lab_sh <- c(fudc = "Up-down", rebuild = "Falling F", decline = "Rising F", flat = "Flat F")
cfgs <- c(spr = "SPR link", both = "SPR + F links")
qs <- c("SB/SB0", "B/Bmsy", "F/Fmsy")
L <- list(); P <- list()
for (i in seq_len(nrow(sc))) {
  scen <- sprintf("%s (seed %d)", lab_sh[sc$shape[i]], sc$seed[i])
  for (cf in names(cfgs)) {
    x <- readRDS(sprintf("ci/spm/fits/case_%s_%s_%d.rds", cf, sc$shape[i], sc$seed[i]))
    sim <- x$sim; f <- x$f; s <- spm_flicc(f)
    if (cf == "spr") {
      L[[length(L) + 1]] <- data.frame(scen = scen, year = sim$years, quantity = rep(qs, each = length(sim$years)),
        est = c(sim$depletion, sim$bbmsy, sim$ffmsy), lci = NA, uci = NA, series = "True (OM)")
    }
    ts <- s$ts[s$ts$quantity %in% qs, ]
    L[[length(L) + 1]] <- data.frame(scen = scen, year = ts$year, quantity = ts$quantity, est = ts$est,
      lci = ts$lci, uci = ts$uci, series = cfgs[cf])
    if (cf == "both") {
      d <- s$ts[s$ts$quantity == "D_spr", ]
      P[[length(P) + 1]] <- data.frame(scen = scen, year = d$year, quantity = "SB/SB0", est = d$est,
                                       series = "D_spr from FLicc SPR")
    }
    attr(L, scen) <- TRUE
  }
}
D <- do.call(rbind, L); Pt <- do.call(rbind, P)
lev <- unique(D$scen); D$scen <- factor(D$scen, lev); Pt$scen <- factor(Pt$scen, lev)
D$quantity <- factor(D$quantity, qs); Pt$quantity <- factor(Pt$quantity, qs)
D$series <- factor(D$series, c("True (OM)", cfgs))
ref <- data.frame(quantity = factor(c("B/Bmsy", "F/Fmsy"), qs), y = 1)
cols <- c("True (OM)" = "#1a1a19", "SPR link" = "#2a78d6", "SPR + F links" = "#eb6834", "D_spr from FLicc SPR" = "#1baf7a")
g <- ggplot(D, aes(year, est, colour = series)) +
  annotate("rect", xmin = 2011.5, xmax = 2024.5, ymin = -Inf, ymax = Inf, fill = "#e8e7e2", alpha = 0.6) +
  geom_hline(data = ref, aes(yintercept = y), colour = "#8a8980", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(data = subset(D, series == "SPR + F links"), aes(ymin = lci, ymax = uci, fill = series),
              colour = NA, alpha = 0.18, show.legend = FALSE) +
  geom_line(aes(linewidth = series)) +
  geom_point(data = Pt, aes(year, est, colour = series), shape = 21, fill = "white", size = 1.8, stroke = 0.8) +
  scale_colour_manual(values = cols, breaks = names(cols), name = NULL) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_linewidth_manual(values = c("True (OM)" = 0.9, "SPR link" = 0.7, "SPR + F links" = 0.7), guide = "none") +
  facet_wrap(~ scen + quantity, ncol = 3, scales = "free_y", labeller = function(l) list(paste0(l$scen, " | ", l$quantity))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = NULL, title = "FLicc-SPM fits vs the alfonsino OM (dev-spm, 18 Nov)",
       subtitle = "Shaded: length-data years 2012-2024 (CPUE from 1990, catch from 1971). Ribbon: 95% CI of SPR + F links.") +
  guides(colour = guide_legend(override.aes = list(linewidth = c(0.9, 0.7, 0.7, 0), shape = c(NA, NA, NA, 21)))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ecebe6"),
        strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold"), panel.spacing = unit(0.6, "lines"))
ggsave("flicc_spm_status.png", g, width = 11, height = 15, dpi = 130, bg = "white")
cat("ok\n")
