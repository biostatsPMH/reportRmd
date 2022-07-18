
<!-- This file is used to create README.md
Note that the README.md document may need updating to change
'\<0.001' to '<0.001'. 
-->

# reportRmd

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/reportRmd)](https://CRAN.R-project.org/package=reportRmd)
<!-- badges: end -->

The goal of reportRmd is to automate the reporting of clinical data in
Rmarkdown environments. Functions include table one-style summary
statistics, compilation of multiple univariate models, tidy output of
multivariable models and side by side comparisons of univariate and
multivariable models. Plotting functions include customisable survival
curves, forest plots, and automated bivariate plots.

## Installation

You can install the development version of reportRmd from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("biostatsPMH/reportRmd")
```

## Examples

### Summary statistics by Sex

``` r
library(reportRmd)
data("pembrolizumab")
rm_covsum(data=pembrolizumab, maincov = 'sex',
covs=c('age','pdl1','change_ctdna_group'),
show.tests=TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Full Sample (n=94)
</th>
<th style="text-align:right;">
Female (n=58)
</th>
<th style="text-align:right;">
Male (n=36)
</th>
<th style="text-align:right;">
p-value
</th>
<th style="text-align:right;">
StatTest
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">age</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
Wilcoxon Rank Sum
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Mean (sd)
</td>
<td style="text-align:right;">
57.9 (12.8)
</td>
<td style="text-align:right;">
56.9 (12.6)
</td>
<td style="text-align:right;">
59.3 (13.1)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Median (Min,Max)
</td>
<td style="text-align:right;">
59.1 (21.1,81.8)
</td>
<td style="text-align:right;">
56.6 (34.1,78.2)
</td>
<td style="text-align:right;">
61.2 (21.1,81.8)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">pdl1</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
Wilcoxon Rank Sum
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Mean (sd)
</td>
<td style="text-align:right;">
13.9 (29.2)
</td>
<td style="text-align:right;">
15.0 (30.5)
</td>
<td style="text-align:right;">
12.1 (27.3)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Median (Min,Max)
</td>
<td style="text-align:right;">
0 (0,100)
</td>
<td style="text-align:right;">
0.5 (0.0,100.0)
</td>
<td style="text-align:right;">
0 (0,100)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Missing
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">change ctdna group</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
Fisher Exact
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Decrease from baseline
</td>
<td style="text-align:right;">
33 (45)
</td>
<td style="text-align:right;">
19 (48)
</td>
<td style="text-align:right;">
14 (42)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Increase from baseline
</td>
<td style="text-align:right;">
40 (55)
</td>
<td style="text-align:right;">
21 (52)
</td>
<td style="text-align:right;">
19 (58)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Missing
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
</tbody>
</table>

### Multiple Univariate Regression Analyses

``` r
rm_uvsum(data=pembrolizumab, response='orr',
covs=c('age','pdl1','change_ctdna_group'))
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
OR(95%CI)
</th>
<th style="text-align:right;">
p-value
</th>
<th style="text-align:right;">
N
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">age</span>
</td>
<td style="text-align:right;">
0.96 (0.91,1.00)
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
94
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">pdl1</span>
</td>
<td style="text-align:right;">
0.97 (0.95,0.98)
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">\<0.001</span>
</td>
<td style="text-align:right;">
93
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">change ctdna group</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">0.002</span>
</td>
<td style="text-align:right;">
73
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Decrease from baseline
</td>
<td style="text-align:right;">
Reference
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
33
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Increase from baseline
</td>
<td style="text-align:right;">
28.74 (5.20,540.18)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
40
</td>
</tr>
</tbody>
</table>

### Tidy multivariable analysis

``` r
glm_fit <- glm(orr~change_ctdna_group+pdl1+cohort,
               family='binomial',
               data = pembrolizumab)
rm_mvsum(glm_fit)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
OR(95%CI)
</th>
<th style="text-align:right;">
p-value
</th>
<th style="text-align:right;">
Global p-value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">change ctdna group</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">0.009</span>
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Decrease from baseline
</td>
<td style="text-align:right;">
Reference
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Increase from baseline
</td>
<td style="text-align:right;">
19.99 (2.08,191.60)
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">pdl1</span>
</td>
<td style="text-align:right;">
0.97 (0.95,1.00)
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">cohort</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">0.004</span>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
A
</td>
<td style="text-align:right;">
Reference
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
B
</td>
<td style="text-align:right;">
2.6e+07 (0e+00,Inf)
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
C
</td>
<td style="text-align:right;">
4.2e+07 (0e+00,Inf)
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
D
</td>
<td style="text-align:right;">
0.07 (4.2e-03,1.09)
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
E
</td>
<td style="text-align:right;">
0.44 (0.04,5.10)
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
</td>
</tr>
</tbody>
</table>

### Combining univariate and multivariable models

``` r
uvsumTable <- rm_uvsum(data=pembrolizumab, response='orr',
covs=c('age','sex','pdl1','change_ctdna_group'),tableOnly = TRUE)
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...

glm_fit <- glm(orr~change_ctdna_group+pdl1,
               family='binomial',
               data = pembrolizumab)
mvsumTable <- rm_mvsum(glm_fit,tableOnly = TRUE)

rm_uv_mv(uvsumTable,mvsumTable)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Unadjusted OR(95%CI)
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
p
</th>
<th style="text-align:right;">
Adjusted OR(95%CI)
</th>
<th style="text-align:right;">
p (adj)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">age</span>
</td>
<td style="text-align:right;">
0.96 (0.91,1.00)
</td>
<td style="text-align:right;">
94
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">sex</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
94
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Female
</td>
<td style="text-align:right;">
Reference
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Male
</td>
<td style="text-align:right;">
0.41 (0.13,1.22)
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">pdl1</span>
</td>
<td style="text-align:right;">
0.97 (0.95,0.98)
</td>
<td style="text-align:right;">
93
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">\<0.001</span>
</td>
<td style="text-align:right;">
0.98 (0.96,1.00)
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">0.02</span>
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style=" font-weight: bold;    ">change ctdna group</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
73
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">0.002</span>
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
<span style=" font-weight: bold;    ">0.004</span>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Decrease from baseline
</td>
<td style="text-align:right;">
Reference
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
Reference
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Increase from baseline
</td>
<td style="text-align:right;">
28.74 (5.20,540.18)
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
24.71 (2.87,212.70)
</td>
<td style="text-align:right;">
</td>
</tr>
</tbody>
</table>

### Survival Probabilities at different times by sex

``` r
rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
 strata="sex",survtimes=seq(12,72,12),survtimeunit='months')
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Time (months)
</th>
<th style="text-align:right;">
At Risk
</th>
<th style="text-align:right;">
Events
</th>
<th style="text-align:right;">
Censored
</th>
<th style="text-align:right;">
Survival Rate (95% CI)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Overall
</td>
<td style="text-align:right;">
94
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
12
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.53 (0.44,0.64)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
24
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.33 (0.24,0.44)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
36
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.28 (0.20,0.40)
</td>
</tr>
<tr>
<td style="text-align:left;">
Female
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
12
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.55 (0.44,0.69)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
24
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.34 (0.24,0.50)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
36
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.29 (0.18,0.45)
</td>
</tr>
<tr>
<td style="text-align:left;">
Male
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
12
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.50 (0.36,0.69)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
24
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.31 (0.18,0.52)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
36
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.27 (0.15,0.48)
</td>
</tr>
</tbody>
</table>

### Logrank test on KM Survival Curves by sex

``` r
rm_survdiff(data=pembrolizumab,time='os_time',status='os_status', 
            covs='sex',digits=1)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
group
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Observed
</th>
<th style="text-align:right;">
Expected
</th>
<th style="text-align:right;">
Median (95%CI)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Overall
</td>
<td style="text-align:right;">
94
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
14.0 (9.0,20.1)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Female
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
41.8
</td>
<td style="text-align:right;">
14.3 (9.7,23.8)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Male
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
22.2
</td>
<td style="text-align:right;">
11.2 (6.1,25.3)
</td>
</tr>
<tr>
<td style="text-align:left;">
Log Rank Test
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
ChiSq = 0.5 on 1 df
</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
p-value = 0.54
</td>
</tr>
</tbody>
</table>

### Logrank test on Survival Curves by sex stratifying by cohort

``` r
rm_survdiff(data=pembrolizumab,time='os_time',status='os_status', 
            covs='sex',strata='cohort',digits=1)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
group
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Observed
</th>
<th style="text-align:right;">
Expected
</th>
<th style="text-align:right;">
Median (95%CI)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Overall
</td>
<td style="text-align:right;">
94
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
14.0 (9.0,20.1)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Female
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
43.0
</td>
<td style="text-align:right;">
14.3 (9.7,23.8)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
Male
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
21.0
</td>
<td style="text-align:right;">
11.2 (6.1,25.3)
</td>
</tr>
<tr>
<td style="text-align:left;">
Log Rank Test
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
ChiSq = 1.9 on 1 df
</td>
</tr>
<tr>
<td style="text-align:left;">
stratified by cohort
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
p-value = 0.83
</td>
</tr>
</tbody>
</table>

### Plotting survival curves

``` r
ggkmcif(response = c('os_time','os_status'),
cov='cohort',
data=pembrolizumab)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

### Plotting odds ratios

``` r
require(ggplot2)
#> Loading required package: ggplot2
forestplot2(glm_fit)
#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Plotting bivariate relationships

These plots are designed for quick inspection of many variables, not for
publication.

``` r
require(ggplot2)
plotuv(data=pembrolizumab, response='orr',
covs=c('age','cohort','pdl1','change_ctdna_group'))
#> Boxplots not shown for categories with fewer than 20 observations.
#> Boxplots not shown for categories with fewer than 20 observations.
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />
