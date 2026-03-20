# reportRmd 0.1.3

## New features

- `rm_mvsum()` and `forestplotMV()` now support `include_unadjusted` parameter to display univariate (unadjusted) estimates alongside multivariable (adjusted) estimates in a single table or plot
- `forestplotMV()` data parameter is now optional - data will be automatically extracted from the model object if not provided
- new `extract_package_details()` function to extract R packages used by a script or markdown file, along with citations.
- new `scrolling_table()` creates scrollable HTML tables for long outputs or keeps regular tables for non-HTML formats

## Minor improvements and fixes

- improved plotting of survival curves
- Cox proportional hazards clustering: fixed bug when clustering is used in coxph models
- Date handling in Word tables fixed bug where tables with date columns containing NA values caused errors in Word output
- rm_compactsum improvements: fixed output formatting issues, fixed bug with difftime variable handling, improved bolding for corrected p-values, better handling of negative numbers
- rm_uvsum fixes: fixed handling of integer predictors, various bug fixes for edge cases
- more robust handling of tidy evaluation
- HTML output improvements: improved display for special characters in HTML rendered tables

# reportRmd 0.1.1

## New features

- new more compact and flexible reporting function (rm_compactsum) that allows more control over each variable summarised and displays text summarising which tests and effect sizes were computed
- new function to apply variable labels to ggplots (replace_plot_labels)
- more comprehensive model summary function to allow plot extensions
- incorporated tidyselect into main summary functions
- support for tidycmprsk models

## Minor improvements and fixes

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

