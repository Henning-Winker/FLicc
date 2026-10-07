# Contributing to FLicc

Thank you for your interest in FLicc. Contributions of all kinds are welcome:
bug reports, documentation fixes, new features, tests and validation datasets.

## Before you start

- **Bugs:** open an issue with a minimal reproducible example (the
  [reprex](https://reprex.tidyverse.org) package helps) and the output of
  `sessionInfo()`. If possible, use one of the example datasets shipped with
  the package (e.g. `data("alfonsino")`).
- **New features or larger changes:** open an issue first to discuss the idea,
  especially if it changes model assumptions, defaults or outputs.
- **Questions about the science:** issues are fine; please cite the relevant
  literature where possible.

## Workflow

1. Fork the repository and create a branch from `main`
   (e.g. `fix-sel-param-check`, `feature-new-selectivity`).
2. Make your changes, with tests where possible.
3. Run `devtools::document()` and `devtools::check()` locally. The check
   should produce no errors or warnings. If you change `src/FLicc.cpp`, make
   sure the TMB model still compiles and that `test.R` runs end to end.
4. Open a pull request describing what changed and why.

All changes reach `main` only through reviewed pull requests, including the
maintainer's own.

## Code conventions

- Keep the style consistent with the surrounding code.
- Document exported functions with roxygen2, including `@param`, `@return`,
  `@examples` and, where relevant, `@references`.
- State units and conventions explicitly (e.g. length in cm, mortality as
  instantaneous annual rate) in documentation and argument checks.
- Keep FLR compatibility: functions should accept and return FLR classes
  (`FLQuant`, `FLPar`, `FLStockLen`, ...) where that is natural.

## Tests and validation

- New functions should come with tests; bug fixes should come with a test
  that fails without the fix.
- Where a function implements a published method (e.g. LBSPR, fishblicc),
  compare its results with the original implementation or a worked example
  from the source publication.
- Use `set.seed()` for anything stochastic so results are reproducible.

## Use of generative AI

AI tools are allowed, but please:

- say in your pull request whether and how you used them;
- make sure you understand, and have tested, every line you submit;
- do not paste code from other packages without checking that its licence is
  compatible with FLicc's.

You are responsible for your contribution regardless of the tools used.

## Licence

FLicc is released under the European Union Public Licence (EUPL), as stated
in the `DESCRIPTION` file. By contributing, you agree that your contributions
are released under the same licence.
