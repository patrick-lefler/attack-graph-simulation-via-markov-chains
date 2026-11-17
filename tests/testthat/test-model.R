test_that("absorption probabilities sum to 1 from every stage", {
  for (sol in list(base_solution, hard_solution)) {
    expect_equal(unname(rowSums(sol$B)), rep(1, t_dim), tolerance = 1e-12)
  }
})

test_that("B is unchanged by removing self-loops (jump-chain invariance)", {
  for (P in list(P_baseline_raw, P_hardened_raw)) {
    J <- jump_chain(P)
    B_direct <- solve_absorbing_markov(P[1:t_dim, 1:t_dim], P[1:t_dim, -(1:t_dim)], step_hours)$B
    B_jump <- solve_absorbing_markov(J[1:t_dim, 1:t_dim], J[1:t_dim, -(1:t_dim)], step_hours)$B
    expect_equal(B_jump, B_direct, tolerance = 1e-10)
  }
})

test_that("the outcome distribution by step matches the closed form", {
  pmf <- outcome_pmf(Q_base, R_base, start = 1)
  expect_equal(sum(pmf), 1, tolerance = 1e-10)
  expect_equal(unname(colSums(pmf)), unname(base_solution$B[1, ]), tolerance = 1e-10)
  expect_equal(sum(seq_len(nrow(pmf)) * rowSums(pmf)), unname(base_solution$steps[1]), tolerance = 1e-8)
})

test_that("no catastrophic outcome is possible before four steps from S1", {
  pmf <- outcome_pmf(Q_base, R_base, start = 1)
  expect_true(all(pmf[1:3, 2:3] == 0))
  expect_gt(sum(pmf[4, 2:3]), 0)
})

test_that("the improved posture raises containment from every stage", {
  expect_true(all(hard_solution$B[, 1] > base_solution$B[, 1]))
})

test_that("headline values quoted in the text are reproduced", {
  expect_equal(unname(round(100 * res$cont_b[1], 1)), 72.5)
  expect_equal(unname(round(100 * res$cont_h[1], 1)), 97.3)
  expect_equal(unname(round(100 * res$cat_b[1], 1)), 27.5)
  expect_equal(unname(round(100 * res$cat_b[3], 1)), 52.9)
  expect_equal(unname(round(res$t_cat_b[1], 1)), 25.0)
  expect_equal(floor_cat_h, 16)
  expect_equal(unname(cross_stage), 3L)
})
