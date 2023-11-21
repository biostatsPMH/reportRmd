# reportRmd 0.1.0

## New features

- cumulative incidence can now be summarised by event type and time
- revised survival curves are now properly aligned and return ggplots
- variable labels are automatically displayed in summary tables
- forest plots now suppose poisson and negative binomial models
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

