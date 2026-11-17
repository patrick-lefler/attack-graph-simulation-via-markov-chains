test_that("every transient row is a probability distribution", {
  for (P in list(P_baseline_raw, P_hardened_raw)) {
    expect_true(all(P >= 0 & P <= 1))
    expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-12)
  }
})

test_that("absorbing rows are the identity", {
  for (P in list(P_baseline_raw, P_hardened_raw)) {
    expect_equal(unname(P[absorbing_names, ]), cbind(matrix(0, r_dim, t_dim), diag(r_dim)))
  }
})

test_that("matrices keep the documented structure", {
  # Only S5 and S6 reach the catastrophic outcomes, and the S4 -> S3 back-edge is the only regression
  for (P in list(P_baseline_raw, P_hardened_raw)) {
    expect_true(all(P[1:4, t_dim + 2:3] == 0))
    lower <- P[1:t_dim, 1:t_dim][lower.tri(diag(t_dim))]
    expect_equal(sum(lower > 0), 1)
    expect_gt(P["S4: Lateral Movement", "S3: Privilege Escalation"], 0)
  }
})

test_that("matrix_from_params rejects malformed rows", {
  bad <- params$transition_matrices$baseline
  bad$S1 <- bad$S1[-1]
  expect_error(matrix_from_params(bad, transient_names, absorbing_names))
})

test_that("scalar parameters carry a type label and plausible values", {
  scalars <- list(
    params$seed, params$step_hours, params$sensitivity$shift_delta, params$sensitivity$prop_rel,
    params$dirichlet$alpha_conc, params$dirichlet$alpha_grid, params$dirichlet$n_draws,
    params$dirichlet$n_sweep, params$dirichlet$n_rank, params$dirichlet$shared_sd_log,
    params$simulation$n_sim
  )
  for (p in scalars) {
    expect_true(p$type %in% c("assumption", "design"))
    expect_true(all(unlist(p$value) > 0))
  }
  expect_true(params$transition_matrices$type == "assumption")
  expect_true(params$dirichlet$alpha_conc$value %in% unlist(params$dirichlet$alpha_grid$value))
})
