# Write summary.md from status.txt and the comparison table.
args <- commandArgs(trailingOnly = TRUE); out <- args[1]; branch <- args[2]; sha <- args[3]
st <- if (file.exists(file.path(out, "status.txt"))) readLines(file.path(out, "status.txt")) else character()
kv <- do.call(rbind, strsplit(st, "=(?=[^=]*$)", perl = TRUE))
lab <- c("0" = "passed", "1" = "FAILED", "2" = "changed (see table)")
s <- c(sprintf("# FLicc CI: %s @ %s", branch, substr(sha, 1, 7)), "",
       sprintf("Run: %s UTC", format(Sys.time(), tz = "UTC")), "",
       "| step | result |", "|---|---|",
       if (length(st)) sprintf("| %s | %s |", kv[, 1], lab[kv[, 2]]) else "| (none) | |", "")
if (file.exists(file.path(out, "compare_main.md")))
  s <- c(s, "## Regression against main (alfonsino, default settings)", "",
         readLines(file.path(out, "compare_main.md")), "")
s <- c(s, "Logs: build.log, tests/*.log, compare_main.log, render_*.log; rendered files in html/.")
writeLines(s, file.path(out, "summary.md")); cat(s, sep = "\n")
