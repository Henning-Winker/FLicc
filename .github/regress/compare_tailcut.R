b <- readRDS("base.rds"); n <- readRDS("new.rds"); fail <- FALSE
cat("\n=== no cut: main vs dev-tailcut ===\n")
for (k in names(b)) {
  d <- c(obj0 = abs(b[[k]]$obj0 - n[[k]]$obj0), obj = abs(b[[k]]$obj - n[[k]]$obj),
         spr = max(abs(b[[k]]$spr - n[[k]]$spr)), fspr = abs(b[[k]]$fspr - n[[k]]$fspr))
  ok <- d["obj0"] < 1e-8 && d["obj"] < 1e-6 && d["spr"] < 1e-5
  cat(sprintf("%-9s %s  %s\n", k, paste(names(d), signif(d, 3), sep = "=", collapse = " "), if (ok) "OK" else "DIFFERS"))
  if (!ok) fail <- TRUE
}
if (fail) quit(status = 1)
