# Load parameters, functions and the closed-form model (fast; no simulations)
# into the test environment. Scripts use paths relative to the project root.
project_root <- normalizePath(testthat::test_path("..", ".."))
withr::with_dir(project_root, {
  source("R/00_setup.R", local = TRUE)
  source("R/01_model.R", local = TRUE)
})
