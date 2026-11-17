# Packages, parameters, plot theme and shared helpers. Sourced first by every
# page, by run_all.R and by the tests. Paths are relative to the project root.

library(tidyverse)
library(scales)
library(patchwork)
library(gt)

params <- yaml::read_yaml("params.yaml")
source("R/functions.R", local = environment())

set.seed(params$seed$value)

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = "#1a1a1a"),
      plot.subtitle = element_text(color = "#555555", size = 10.5),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
)

pct1 <- function(x) sprintf("%.1f%%", 100 * x)
pp1 <- function(x) sprintf("%.1f pp", 100 * x)
