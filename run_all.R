# Rebuild everything from the project root: Rscript run_all.R
#   1. Check that installed packages match renv.lock.
#   2. Run the unit tests; stop on any failure.
#   3. Re-execute and render every page. Pages are rendered one at a time
#      because a project-wide `quarto render` skips pages whose .qmd file is
#      unchanged (freeze: auto), even when R/ or params.yaml changed.

status <- renv::status()
if (!isTRUE(status$synchronized)) {
  stop("Packages differ from renv.lock. Run renv::restore() first.", call. = FALSE)
}

testthat::test_dir("tests/testthat", stop_on_failure = TRUE)

pages <- yaml::read_yaml("_quarto.yml")$project$render
for (page in pages) {
  message("Rendering ", page)
  if (system2("quarto", c("render", page)) != 0) {
    stop("quarto render failed for ", page, call. = FALSE)
  }
}

message("Done: site written to docs/")
