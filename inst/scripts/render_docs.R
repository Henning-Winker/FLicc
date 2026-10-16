# Render the FLicc guideline documents to docs/ (GitHub Pages) -----------------
# Sources live in inst/rmd/, rendered HTML goes to docs/ (served by GitHub
# Pages from the /docs folder). Run from the package folder, after installing
# FLicc with the current R/ (the Rmd calls the exported functions):
#   source("inst/scripts/render_docs.R")
#
# Slow MSE runs (mse-02 onwards) are cached as .rds in rdata/ (not pushed;
# listed in .gitignore and .Rbuildignore). Delete a cached file to re-run it.

dir.create("rdata", showWarnings = FALSE)
options(FLicc.rdata   = normalizePath("rdata", mustWork = TRUE),
        FLicc.run_mse = TRUE)

docs <- c(
  mse_modules = "inst/rmd/mse-01-flicc-modules.Rmd"
  # mse_alfonsino  = "inst/rmd/mse-02-alfonsino-mse.Rmd",  # OM + sampling.lfd.oem + flicc.sa
  # mse_robustness = "inst/rmd/mse-03-robustness.Rmd"      # misspecified M, Linf, CVL, selectivity
)

for (d in docs) {
  rmarkdown::render(d, output_dir = "docs", envir = new.env(),
                    output_options = list(self_contained = TRUE))
}

# GitHub Pages: serve files as-is (no Jekyll processing)
if (!file.exists("docs/.nojekyll")) file.create("docs/.nojekyll")
