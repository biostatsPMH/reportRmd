# reportRmd 0.1.4

## New features

- new `predsum()` summarises binary prediction model performance (sensitivity, specificity, PPV, NPV, accuracy, AUC, AUC-PR, F1, Brier score) with Clopper-Pearson or bootstrap CIs; `output` presets "epi", "data-science", "all" or a custom set of statistics
- support for Firth-penalised logistic regression (`logistf`) models in `rm_mvsum()`, with penalised likelihood ratio global p-values and profile penalised likelihood CIs
- `rm_uvsum()` and `rm_mvsum()` now support clustered `coxph` models, reporting robust Wald tests
- `rm_survsum()` and `rm_cifsum()` now accept factor or character status variables; `rm_cifsum()` works when there are no censored observations
- `rm_covsum()` and `rm_compactsum()` Monte Carlo p-values and `predsum()` bootstrap CIs are now reproducible (`B` and `seed` arguments; RNG state is restored)
- VIF calculation now uses `type = "predictor"` in the presence of interactions

## Minor improvements and fixes

- `rm_mvsum(include_unadjusted = TRUE)`: fixed failures with `coxph` models, non-95% `CIwidth`, formulas stored in variables, dot formulas and transformed terms; continuous covariates and continuous:continuous interactions are no longer dropped from the unadjusted column; `digits`, `CIwidth`, `whichp` and `showEvent` are now passed to the univariate models
- `rm_uvsum()`: excluded covariates (dates, the response, single-value variables) are no longer fit; `reflevel` now works with ordered factors; `removeInf` is now applied
- `rm_covsum()` and `rm_compactsum()` now use expected rather than observed cell counts to choose the exact test, so the two functions agree; the grouping variable is dropped from covariates; non-numeric p-values (e.g. "excl") are no longer adjusted
- `rm_covsum()`: fixed error with date variables on Windows
- `ggkmcif2()`: fixed curves being cut off when `times` was shorter than the data range; accepts US spelling of arguments; works with ggplot2 >= 3.4.0 (now the stated minimum)
- `psthr()` displays "Not Estimable" only when all values are NA
- fixed display of special characters in column headers
- fixed unescaped variable names in `getVarLevels()`, Global p-values being adjusted and double-formatted in bolded output, and a tidyselect deprecation warning for character column names
- consistent tidy evaluation (`enquo`) for variable selection
- vignette: corrected `include_unadjusted` section; new Prediction Metrics and Reproducibility sections
- expanded test suite

# reportRmd 0.1.1

## New features

- new more compact and flexible reporting function (rm_compactsum) that allows more control over each variable summarised and displays text summarising which tests and effect sizes were computed
- new function to apply variable labels to ggplots (replace_plot_labels)
- more comprehensive model summary function to allow plot extensions
- incorporated tidyselect into main summary functions
- support for tidycmprsk models
- ability to output to a scrolling table in HTML formats


## Minor improvements and fixes

- fixes to special characters in HTML output
- removed ability to force Wald confidence intervals, confidence intervals now computed by the updated confint in base R
- rm_mvsum and rm_uvsum have been updated in the backend to enable easier extendibility to different model types
- bug fix in nestTable so variables order properly with repeat level names
- numerous bug fixes in formatting of tables, especially with duplicated variable and level names
- added automated testing of the proportional hazards assumption when reporting coxph models
- bug fix for computing global p-values in models with offset terms
- documentation updates


# reportRmd 0.1.0

## New features

- cumulative incidence can now be summarised by event type and time
- revised survival curves are now properly aligned and return ggplots
- variable labels are automatically displayed in summary tables
- forest plots now support poisson and negative binomial models
- options to globally set the number of digits for a report

## Minor improvements and fixes

- updated excelCol to handle column ranges (ie A-Z, instead of A,B,C,...)
- updated mvsum to accept negative binomial models from MASS
- updated uvsum to fit negative binomial models
- differentiated univariate and multivariable estimates in forest plots with shapes
- added functionality to export tables in html and latex formats
- added functionality to summarise date variables
- previously the digit argument did not work for variables with more than one class, this has been fixed
- long data calls (as may be created with tidy statements) are now supported in summary functions
- table bolding for interaction terms is now consistent
- Total N,  Number of Events now properly reported when interactions are present
- nestTable now works properly with two-column tables
- covTitle argument restored in rm_covsum
  

# reportRmd 0.0.2

- Added functionality to include effect sizes in summary tables
- forestplot2 now works with geeglm objects

Bug fixes:
  - fixed bug where some p-values didn't output after update to R 4.2.2
  - special characters now better shown in tables
  - fixed bug so that NA are always shown as blank
  - fixed bug in nestTable so that variable indenting is consistent

# reportRmd 0.0.1

* 5 October 2022 first CRAN release

