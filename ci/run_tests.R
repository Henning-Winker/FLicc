# Run each test script in its own R process with the branch build.
# Usage: Rscript ci/run_tests.R <repo> <lib> <outdir>
# Scripts: test_dev.R (repo root) and ci/tests/*.R, in that order.
args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]; lib <- normalizePath(args[2]); out <- args[3]
scripts <- c(file.path(repo, "test_dev.R"),
             sort(list.files(file.path(repo, "ci", "tests"), "\\.R$", full.names = TRUE)))
scripts <- scripts[file.exists(scripts)]
dir.create(file.path(out, "tests"), showWarnings = FALSE, recursive = TRUE)
for (s in scripts) {
  log <- file.path(out, "tests", paste0(tools::file_path_sans_ext(basename(s)), ".log"))
  expr <- sprintf(".libPaths(c('%s', .libPaths())); source('%s', echo = TRUE, max.deparse.length = 300)",
                  lib, normalizePath(s))
  t <- system.time(st <- system2("Rscript", c("-e", shQuote(expr)), stdout = log, stderr = log))
  cat(sprintf("test:%s=%d\n", basename(s), as.integer(st != 0)),
      file = file.path(out, "status.txt"), append = TRUE)
  cat(sprintf("%s: %s (%.0f s)\n", basename(s), if (st == 0) "passed" else "FAILED", t[["elapsed"]]))
}
