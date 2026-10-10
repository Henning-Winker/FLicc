# FLicc CI

Every push to a branch named `dev-*` or `ci-*` runs `.github/workflows/flicc-ci.yml`
(it can also be started by hand from the Actions tab):

1. builds the branch and `main`;
2. runs `test_dev.R` and every script in `ci/tests/` (each in its own R session);
3. compares default fits on the alfonsino example with `main`
   (`ci/compare_main.R`): identical at the same parameters, or flagged as changed;
4. renders `.Rmd` files changed on the branch;
5. publishes `summary.md`, logs and rendered html to the branch `ci-out/<branch>`.

Read the results without opening GitHub:

```
git fetch origin ci-out/<branch>
git show origin/ci-out/<branch>:summary.md
```

Add a test by putting an R script in `ci/tests/`; it fails if the script errors
(use `stopifnot()` for checks). R packages are cached between runs.
