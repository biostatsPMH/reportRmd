# ggkmcif3 Testing Scripts
# Comprehensive test suite for all parameters in ggkmcif3Parameters

# Load required libraries
library(survival)
library(ggplot2)
library(dplyr)
library(cowplot)

# Create sample datasets for testing ----

# Sample KM data (time-to-event with censoring)
set.seed(123)
n <- 200
km_data <- data.frame(
  time = rexp(n, rate = 0.1) |> round(2),
  status = rbinom(n, 1, 0.7),  # 70% event rate
  treatment = factor(sample(c("Control", "Treatment"), n, replace = TRUE)),
  age_group = factor(sample(c("Young", "Old"), n, replace = TRUE)),
  age_cont = rnorm(n, 60, 15) |> round(1),
  sex = factor(sample(c("Male", "Female"), n, replace = TRUE))
)

forestplotUV(response="orr", covs=c("change_ctdna_group", "sex", "age", "l_size"),
              data=pembrolizumab, family='binomial')
km_data2 <- data.frame(
  time = rexp(n, rate = 0.1) |> round(2),
  status = rbinom(n, 1, 0.7),  # 70% event rate

  # 2 levels
  treatment = factor(sample(c("Control", "Treatment"), n, replace = TRUE)),

  # 3 levels
  stage = factor(sample(c("Early", "Intermediate", "Advanced"), n, replace = TRUE)),

  # 4 levels
  region = factor(sample(c("North", "South", "East", "West"), n, replace = TRUE)),

  # 5 levels
  risk_score = factor(sample(c("Very Low", "Low", "Moderate", "High", "Very High"),
                             n, replace = TRUE,
                             prob = c(0.15, 0.25, 0.3, 0.2, 0.1))),  # Weighted sampling

  # Keep your continuous variables
  age_cont = rnorm(n, 60, 15) |> round(1),
  sex = factor(sample(c("Male", "Female"), n, replace = TRUE))
)
# Sample CIF data (competing risks)
cif_data <- data.frame(
  time = rexp(n, rate = 0.1) |> round(2),
  status = sample(0:2, n, replace = TRUE, prob = c(0.3, 0.4, 0.3)),  # 0=censor, 1=event1, 2=event2
  treatment = factor(sample(c("Control", "Treatment"), n, replace = TRUE)),
  age_group = factor(sample(c("Young", "Old"), n, replace = TRUE)),
  risk_score = rnorm(n, 0, 1) |> round(2)
)

# Quick Test
data("pembrolizumab")
data=pembrolizumab;cov="cohort";response=c("os_time","os_status")
ggkmcif3(response,cov,data,censor.size=2.5,censor.stroke=7)
ggkmcif3(response,cov,data,times=seq(0,42,6))
ggkmcif3(response,cov,data,times=seq(0,42,6),xlim=c(0,30))

ggkmcif3(response=c("time","status"),cov="sex",km_data)
ggkmcif2(response=c("time","status"),cov="age_group",cif_data)

p <- ggkmcif3(response=c("time","status"),cov="risk_score",km_data2,returns = T)
get_theme_base_size <- function() {
  current_theme <- theme_get()

  # The base_size is stored in the text element
  if (!is.null(current_theme$text) && !is.null(current_theme$text$size)) {
    return(current_theme$text$size)
  }
  # Fallback to ggplot2 default
  return(11)
}
line_size_in_inches <- get_theme_base_size()/72
dev_height_inches <- dev.size("in")[2]
p_risk_lines <- length(unique(p$table$data$strata))+1.5
p2_risk_height <- line_size_in_inches*p_risk_lines
p1_height <-  dev_height_inches - p2_risk_height
rel_height = c(1,(dev_height_inches-p1_height)/p1_height)
cowplot::plot_grid(p$plot,p$table,rel_heights = rel_height,ncol=1)


df <- ggkmcif2(response=c("time","status"),cov="age_group",cif_data,censor.size=2.5,censor.stroke=7)

# Find which file contains the broken system2 call ----
fs::dir_ls(getwd(), recurse = TRUE, glob = "*.R") |>
  purrr::keep(\(f) any(grepl("quarto.*-V|TMPDIR", readLines(f, warn = FALSE)))) |>
  print()
