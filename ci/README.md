# FLicc CI

Every push to a branch named `dev-*` or `ci-*` runs `.github/workflows/flicc-ci.yml`
(it can also be started by hand from the Actions tab):

1. builds the branch and `main`;
2. runs `test_dev.R` and every script in `ci/tests/` (each in its own R session);
3. compares default fits on the alfonsino example with `main`
   (`ci/compare_main.R`): identical at the same parameters, or flagged as changed;
4. renders `.Rmd` files changed on the branch. The user guide is split into
   chapters (`inst/guide/chapters/`): a push renders only the chapters it changed,
   with `inst/guide/preview.Rmd` (a few minutes). The full guide is rendered only
   when the commit message contains `[render]` (e.g. before merging into main)
   or when the workflow is started by hand with `render` ticked;
5. publishes `summary.md`, logs and rendered html to the branch `ci-out/<branch>`.

Read the results without opening GitHub:

```
git fetch origin ci-out/<branch>
git show origin/ci-out/<branch>:summary.md
```

Add a test by putting an R script in `ci/tests/`; it fails if the script errors
(use `stopifnot()` for checks). R packages are cached between runs.

## Editing the user guide

The guide `inst/guide/FLicc_user_guide.Rmd` only holds the header and includes
the chapters from `inst/guide/chapters/`; edit the chapters. To check one or
two chapters locally:

```r
rmarkdown::render("inst/guide/preview.Rmd",
                  params = list(chapters = "07_selectivity, 10_single_year"))
```

The preview first runs the quick start and the reference fit (chapters 03 and
04) silently. A chapter that uses objects from another chapter names it in its
first line, e.g. `<!-- needs: 12_sample_size -->`. Before merging, push with
`[render]` in the commit message to render the whole guide.
